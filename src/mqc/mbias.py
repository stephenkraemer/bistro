import os
import pickle
import itertools
import warnings

from abc import ABCMeta, abstractmethod
from itertools import product
from math import floor, ceil
from pathlib import Path
from typing import Dict

import numpy as np
import pandas as pd

import matplotlib

matplotlib.use('Agg')  # import before pyplot import!
import matplotlib.pyplot as plt
import seaborn as sns
import plotnine as gg

import mqc.flag_and_index_values as mfl
from mqc.pileup.bsseq_pileup_read import BSSeqPileupRead
from mqc.pileup.pileup import MotifPileup
from mqc.utils import convert_array_to_df
from mqc.visitors import Counter

b_inds = mfl.bsseq_strand_indices
b_na_ind = mfl.bsseq_strand_na_index
m_flags = mfl.methylation_status_flags

idxs = pd.IndexSlice


class MbiasCounter(Counter):
    """Count stratified M-bias stats

    *Implementation notes:*
    The Fragment length dimension includes the length 0, so that it can be
    indexed by 1-based values. The read position indexes on the other hand
    are zero-based, for better interaction with the C/cython parts of the
    program.
    """

    def __init__(self, config: Dict):

        self.save_stem = config["paths"]["mbias_counts"]
        self.max_read_length = config["data_properties"]["max_read_length_bp"]
        sts = config["stats"]
        self.max_flen = sts["max_flen"]
        self.max_single_flen = sts["max_flen_with_single_flen_resolution"]
        self.flen_bin_size = sts["flen_bin_size"]
        self.max_phred = sts["max_phred"]
        self.phred_bin_size = sts["phred_bin_size"]
        # TODO: get from IndexFile
        self.seq_context_size = config["stats"]["seq_context_size"]

        dim_names = ["seq_context", "bs_strand", "flen",
                     "phred", "pos", "meth_status"]

        flen_breaks = (
            list(range(0, self.max_single_flen + 1)) +
            list(range(self.max_single_flen + 1,
                       self.max_flen + 2,
                       self.flen_bin_size)))
        if flen_breaks[-1] != self.max_flen + 1:
            flen_breaks.append(self.max_flen + 1)

        flen_intervals = (pd.IntervalIndex
                          .from_breaks(flen_breaks, closed="left")
                          .values.tolist())
        self.last_flen_bin_idx = len(flen_intervals) - 1

        phred_breaks = list(range(0, self.max_phred + 2, self.phred_bin_size))
        if phred_breaks[-1] != self.max_phred + 1:
            phred_breaks.append(self.max_phred + 1)
        phred_intervals = (pd.IntervalIndex
                           .from_breaks(phred_breaks, closed="left")
                           .values.tolist())
        self.last_phred_bin_idx = len(phred_intervals) - 1

        self.seq_ctx_idx_dict, self.binned_motif_to_index_dict = get_sequence_context_to_array_index_table(
            self.seq_context_size)

        # Note: 1-based position labels in dataframe, 0-based position indices
        # in array
        ordered_seq_contexts = [seq_context
                                for seq_context, idx in sorted(
                self.binned_motif_to_index_dict.items(),
                key=lambda x: x[1])]

        def get_categorical(ordered_str_labels):
            return pd.Categorical(ordered_str_labels,
                                  categories=ordered_str_labels,
                                  ordered=True)

        dim_levels = [get_categorical(ordered_seq_contexts),
                      get_categorical(['c_bc', 'c_bc_rv', 'w_bc', 'w_bc_rv']),
                      flen_intervals,
                      phred_intervals,
                      range(1, self.max_read_length + 1),
                      get_categorical(['n_meth', 'n_unmeth'])]

        array_shape = [len(ordered_seq_contexts),
                       4,  # BSSeq-strands
                       len(flen_intervals),
                       len(phred_intervals),
                       self.max_read_length,
                       2]  # meth status
        counter_array = np.zeros(array_shape, dtype='u8')

        super().__init__(dim_names=dim_names,
                         dim_levels=dim_levels,
                         counter_array=counter_array,
                         save_stem=self.save_stem)

    def process(self, motif_pileup: MotifPileup):
        """Extract M-bias stats from MotifPileup

        Reads are discarded if

        - they have a qc_fail_flag
        - their bsseq strand could not be determined
        - they have methylation calling status: NA, SNP or Ref

        Stratified by:

        - motif
        - BSSeq-strand
        - fragment length
        - position in read
        - methylation status
        """

        # Can't handle Ns and too short motifs
        try:
            seq_ctx_idx = self.seq_ctx_idx_dict[
                motif_pileup.idx_pos.seq_context]
        except KeyError:
            return  # can't count this MotifPileup

        curr_read: BSSeqPileupRead
        for curr_read in motif_pileup.reads:
            # TODO: currently this sorts out any qc_fail, including phred
            # score fails, phred score fails should be kept here
            if (curr_read.qc_fail_flag
                or curr_read.bsseq_strand_ind == b_na_ind):
                continue

            meth_status_flag = curr_read.meth_status_flag
            if meth_status_flag == m_flags.is_methylated:
                meth_status_index = 0
            elif meth_status_flag == m_flags.is_unmethylated:
                meth_status_index = 1
            else:  # SNP, Ref, NA
                continue

            tlen = abs(curr_read.alignment.template_length)
            if tlen < self.max_single_flen:
                tlen_idx = tlen
            elif tlen > self.max_flen:
                tlen_idx = self.last_flen_bin_idx
            else:
                tlen_idx = (self.max_single_flen
                            + ceil((tlen - self.max_single_flen)
                                   / self.flen_bin_size))

            phred = curr_read.baseq_at_pos
            if phred > self.max_phred:
                phred_idx = self.last_phred_bin_idx
            else:
                phred_idx = floor(phred / self.phred_bin_size)

            event_class = (seq_ctx_idx,
                           curr_read.bsseq_strand_ind,
                           tlen_idx,
                           phred_idx,
                           curr_read.pos_in_read,
                           meth_status_index)

            self.counter_array[event_class] += 1


def get_sequence_context_to_array_index_table(motif_size: int):
    """ Return dicts used in mapping sequence contexts to counter array indices

    Parameters
    ----------
    motif_size: int
        size of motifs (e.g. 3 or 5 bp)

    Returns
    -------
    Dict[str, str]
        mapping of full motif [ATCG] to array integer index. Several motifs
        may map to the same index
    Dict[str, int]
        mapping of binned motif to array integer index. Every mapping is
        unique.
    """

    if motif_size % 2 != 1:
        raise ValueError("Motif size must be an uneven number")

    all_bases = ['C', 'G', 'T', 'A']
    three_letter_bases = ['C', 'G', 'W']

    n_bp_per_side = (motif_size - 1) // 2
    binned_bases_set = ([three_letter_bases] * n_bp_per_side
                        + [['C']] + [three_letter_bases] * n_bp_per_side)

    # note that indicies are given in alphabetical sorting order
    all_binned_motifs = sorted([''.join(motif)
                                for motif in product(*binned_bases_set)])

    binned_motif_to_idx_mapping = {motif: i
                                   for i, motif in
                                   enumerate(all_binned_motifs)}

    l2 = [all_bases] * n_bp_per_side + [['C']] + [all_bases] * n_bp_per_side
    all_5bp_motifs = [''.join(motif) for motif in product(*l2)]

    _5bp_to_three_letter_motif_index_mapping = {
        motif: binned_motif_to_idx_mapping[
            motif.translate(str.maketrans('CGTA', 'CGWW'))]
        for motif in all_5bp_motifs}

    return _5bp_to_three_letter_motif_index_mapping, binned_motif_to_idx_mapping


def map_seq_ctx_to_motif(seq_ctx, use_classical=True):
    """Map sequence context strings containing [ACGTW] to motifs

    Motifs may be classical: [CG, CHG, CHH] or extended (composed of C,G,W)
    # TODO: Ns? at the end of chroms?
    """

    middle_idx = len(seq_ctx) // 2
    if seq_ctx[middle_idx:(middle_idx + 2)] == 'CG':
        return 'CG'

    if use_classical:
        base_mapping = str.maketrans('ACTW', 'HHHH')  # G unchanged
    else:
        base_mapping = str.maketrans('AT', 'WW')  # CGW unchanged

    seq_suffix = seq_ctx[(middle_idx + 1):(middle_idx + 3)]
    motif_suffix = seq_suffix.translate(base_mapping)

    motif = 'C' + motif_suffix
    return motif


class CuttingSites(metaclass=ABCMeta):
    # TODO: Make sure that plateau intervals are 1-based labels, not zero-based indices
    @abstractmethod
    def get_array(self):
        pass

    @abstractmethod
    def get_df(self):
        pass


class FixedRelativeCuttingSites(CuttingSites):
    """Calculate and provide M-bias cutting sites in dataframe and array format"""

    def __init__(self, config):

        self._minimal_cutting_sites_array = np.array([])
        self._minimal_cutting_sites_df = pd.DataFrame()

        self.max_read_length_bp = config["data_properties"][
            "max_read_length_bp"]
        self.max_flen = config["trimming"]["max_flen_considered_for_trimming"]
        self.n_bp_to_discard_at_frag_ends = config["trimming"][
            "relative_to_fragment_ends_dict"]

        self.dim_levels = [b_inds._fields,
                           range(0, self.max_flen + 1),
                           ['left_cut_end', 'right_cut_end']]
        self.dim_names = ['bs_strand', 'flen', 'cut_end']
        self.value_column_name = 'cut_pos'

    def get_array(self):
        """
        Given the number of bp to be removed from the 5' and 3' fragment ends,
        define cutting sites expressed as zero-based positions in the read.

        The cutting sites are given as plateau start and end, i.e. as the first
        and last acceptable position in the read.

        Cutting sites are given for every (bsseq_strand, flen) stratum with
        flen \in [0, max_flen]

        The (input) 'relative cutting sites' may look like this:
        {
            'w_bc': [0, 9]
            'c_bc': [0, 9]
            'w_bc_rv': [9, 0]
            'c_bc_rv': [9, 0]
        }

        The, the result should be a multidimensional array, which in table format would look like this
        bsseq_strand      where       flen         pos_in_read
        'w_bc'            'start'     1            10
                                      2            10
                                      ...
                                      max_flen
                          end'        1            90
                                      ...
        'w_bc_rv'         ...         ...
        """

        if not self._minimal_cutting_sites_array.size:
            # Flen is 1-based
            min_trimmsite_arr = np.zeros([
                4, (self.max_flen + 1), 2], dtype=np.uint32)

            for bsseq_strand_name, bsseq_strand_index in b_inds._asdict().items():

                # Set start position
                # Start position is the zero-based plateau start in the read
                # the relative_cutting_site_dict gives the number of bases
                # to be removed from the start of the read
                min_trimmsite_arr[bsseq_strand_index, :, 0] = (
                    self.n_bp_to_discard_at_frag_ends[bsseq_strand_name][0])

                # Set last included position (zero-based)
                for curr_flen in range(0, self.max_flen):
                    max_allowed_pos_in_fragment = (
                        curr_flen -
                        self.n_bp_to_discard_at_frag_ends[bsseq_strand_name][
                            1] - 1)
                    max_allow_pos_in_read = (max_allowed_pos_in_fragment
                                             if max_allowed_pos_in_fragment <= (
                        self.max_read_length_bp - 1)
                                             else (
                        self.max_read_length_bp - 1))
                    min_trimmsite_arr[
                        bsseq_strand_index, curr_flen, 1] = max_allow_pos_in_read

            # cache
            self._minimal_cutting_sites_array = min_trimmsite_arr
            return min_trimmsite_arr

        else:
            # return cached value
            return self._minimal_cutting_sites_array

    def get_df(self):
        if self._minimal_cutting_sites_df.empty:
            self._minimal_cutting_sites_df = convert_array_to_df(
                arr=self.get_array(),
                dim_levels=self.dim_levels,
                dim_names=self.dim_names,
                value_column_name=self.value_column_name
            )
        return self._minimal_cutting_sites_df


class AdjustedCuttingSites(CuttingSites):
    def __init__(self, mbias_df, config):
        """Computes adjusted cutting sites during init
        ~ 1 min
        """
        self.max_flen = config['trimming']['max_flen_considered_for_trimming']
        self._df = self._compute_df(mbias_df, config)
        self._arr = np.array([])

    def get_df(self) -> pd.DataFrame:
        return self._df

    def get_array(self):
        if not self._arr.size:
            self._arr = convert_cutting_sites_df_to_array(self._df,
                                                          self.max_flen)
        return self._arr

    @staticmethod
    def _compute_df(mbias_df, config):
        print("Calculating cutting sites")
        df = (mbias_df
              .loc['CG', :]
              .groupby(axis='index',
                       level=['bs_strand', 'flen'])
              .apply(fit_normalvariate_plateau, config))
        df.columns.name = 'cut_end'
        df = df.stack().to_frame('cut_pos')
        df = df.sort_index()
        return df


def convert_cutting_sites_df_to_array(df, max_flen):
    # May be shared between AdjustedCuttingSites and
    # UserCuttingSites in the future
    df = df.reset_index()
    level_name_to_index_mappings = {
        'bs_strand': b_inds._asdict(),
        'cut_end': {'left_cut_end': 0, 'right_cut_end': 1},
    }
    df_int_idx = df.replace(level_name_to_index_mappings)
    arr = np.zeros([4, max_flen + 1, 2])
    arr[df_int_idx['bs_strand'].tolist(),
        df_int_idx['flen'].tolist(),
        df_int_idx['cut_end'].tolist()] = df_int_idx['cut_pos']
    return arr


def fit_normalvariate_plateau(group_df: pd.DataFrame, config) -> pd.Series:
    """Find the longest possible plateau of good quality

    Good quality: std along the plateau below a threshold, ends of the plateau
    not different from other read positions (not implemented yet)

    Algorithm:
    1. Start with longest possible plateau length (full read)
    2. Check if the plateau has good quality
    3. If yes: use the plateau and break
    4. Decrease the plateau length by one
    5. Check all possibilities of finding a plateau of the given length within the read
    6. If there are one or more way of fitting the plateau with good quality: break and return the best solution
    7. If not: go to step 4
    """

    min_perc = config["trimming"]["min_plateau_perc"]
    max_std = config["trimming"]["max_std_within_plateau"]
    min_flen = config['trimming']['min_flen_considered_for_trimming']

    beta_values = group_df['beta_value']

    effective_read_length = len(beta_values)
    if effective_read_length < min_flen:
        return pd.Series([0, 0],
                         index=['left_cut_end', 'right_cut_end'])

    if beta_values.isnull().all():
        return pd.Series([0, 0],
                         index=['left_cut_end', 'right_cut_end'])

    min_plateau_length = int(effective_read_length * min_perc)

    for plateau_length in range(effective_read_length, min_plateau_length - 1,
                                -1):
        max_start_pos = effective_read_length - plateau_length
        std_to_beat = max_std
        best_start = None
        for start_pos in range(0, max_start_pos):
            end_pos = start_pos + plateau_length
            curr_beta_values = beta_values.iloc[start_pos:end_pos]
            curr_std = curr_beta_values.std()
            if curr_std < std_to_beat:
                plateau_height = curr_beta_values.mean()
                left_end_bad = (abs(curr_beta_values[
                                    0:4] - plateau_height) > 2 * curr_std).any()
                right_end_bad = (abs(curr_beta_values[
                                     -4:] - plateau_height) > 2 * curr_std).any()
                if not (left_end_bad or right_end_bad):
                    std_to_beat = curr_std
                    best_start = start_pos
                    best_end = end_pos
        if best_start is not None:
            break
    else:
        best_start = 0
        best_end = 0

    return pd.Series([best_start, best_end],
                     index=['left_cut_end', 'right_cut_end'])


def cutting_sites_plot(cutting_sites_df, config):
    g = sns.FacetGrid(data=cutting_sites_df.reset_index(),
                      col='bs_strand', col_wrap=2,
                      col_order=['w_bc', 'w_bc_rv', 'c_bc', 'c_bc_rv'],
                      hue='cut_end')
    g.set(xticks=range(0, 500, 50))
    g = g.map(plt.plot, 'flen', 'cut_pos', alpha=0.7)
    g.savefig(config['paths']['adj_cutting_sites_plot'])


def compute_mbias_stats_df(mbias_counter_fp_str):
    """Compute DataFrame of Mbias-Stats

    Parameters
    ----------
    mbias_counter_fp_str: str
        Path to pickle of MbiasCounter

    Returns
    -------
    pd.DataFrame:
        Index: ['motif', 'seq_context', 'bs_strand', 'flen', 'phred', 'pos']
        Columns: ['n_meth', 'n_unmeth', 'beta_value']

    takes approx. 7 min, mainly due to index sorting, which may be unnecessary.
    Sorting to be future-proof at the moment.
    """

    print("Creating mbias stats dataframe")

    # only necessary while interval levels are coded in MbiasCounter.__init__
    mbias_counter = pd.read_pickle(mbias_counter_fp_str)
    mbias_counter.dim_levels[-3] = [x.right for x in
                                    mbias_counter.dim_levels[-3]]
    mbias_counter.dim_levels[-4] = [x.right for x in
                                    mbias_counter.dim_levels[-4]]

    print("Reading data in, removing impossible strata (pos > flen)")
    # convert MbiasCounter.counter_array to dataframe
    # remove positions > flen
    # takes approx. 2 min
    mbias_stats_df = (mbias_counter
                      .get_dataframe()
                      # TEST
                      # for interactive testing
                      # .loc[['WWCGW', 'CCCGC', 'GGCGG', 'CGCGC', 'GCCGG'], :]
                      # \TEST
                      .groupby(['flen', 'pos'])
                      .filter(lambda group_df: (group_df.name[0] - 1
                                                >= group_df.name[1]))
                      )

    print("Adding beta values")
    # Add beta value, approx. 1 min
    mbias_stats_df_with_meth = (mbias_stats_df
                                .loc[:, 'counts']
                                .unstack('meth_status'))
    # columns is categorical index,
    # replace with object index for easier extension
    mbias_stats_df_with_meth.columns = ['n_meth', 'n_unmeth']
    mbias_stats_df_with_meth['beta_value'] = compute_beta_values(
        mbias_stats_df_with_meth)

    print("Adding motif labels to index")
    # prepend motif level to index (~30s)
    mbias_stats_df_with_motif_idx = prepend_motif_level_to_index(
        mbias_stats_df_with_meth)

    print("Sorting")
    # sort - takes approx. 3 min with IntegerIndices for flen and phred
    mbias_stats_df_with_motif_idx.sort_index(inplace=True)

    return mbias_stats_df_with_motif_idx


def prepend_motif_level_to_index(mbias_stats_df):
    # first, create motif column. This way is much faster than using
    # apply(map_seq_ctx_to_motif) on seq_contexts (15s)

    # create motif column, assuming all rows have CHH seq_contexts
    mbias_stats_df["motif"] = pd.Categorical(
        ['CHH'] * len(mbias_stats_df),
        categories=['CG', 'CHG', 'CHH'],
        ordered=True)

    # find the seq_context labels which are CG and CHG
    seq_contexts = (mbias_stats_df.index.levels[0].categories.tolist())
    motifs = [map_seq_ctx_to_motif(x) for x in seq_contexts]
    motif_is_cg = [x == 'CG' for x in motifs]
    motif_is_chg = [x == 'CHG' for x in motifs]
    cg_seq_contexts = itertools.compress(seq_contexts, motif_is_cg)
    chg_seq_contexts = itertools.compress(seq_contexts, motif_is_chg)

    # replace CHH motif label with correct labels at CG / CHG seq_contexts
    mbias_stats_df.loc[cg_seq_contexts, "motif"] = "CG"
    mbias_stats_df.loc[chg_seq_contexts, "motif"] = "CHG"

    # Then set as index and prepend
    # This takes 15s with interval flen and phred indices
    # *with IntervalIndices for flen and phred, it takes 6-7 min*
    index_cols = ['motif', 'seq_context', 'bs_strand', 'flen', 'phred', 'pos']
    return (mbias_stats_df
            .set_index(["motif"], append=True)
            .reorder_levels(index_cols, axis=0))


def compute_classic_mbias_stats_df(mbias_stats_df):
    """Compute Mbias-Stats df with only 'standard' levels

    Takes ~30s
    """

    print("Creating classic mbias stats dataframe")
    return (mbias_stats_df
            .groupby(level=['motif', 'bs_strand', 'flen', 'pos'])
            .sum()
            .groupby(['flen', 'pos'])
            .filter(lambda group_df: (group_df.name[0] - 1
                                      >= group_df.name[1]))
            .assign(beta_value=compute_beta_values)
            )


def mask_mbias_stats_df(mbias_stats_df, cutting_sites_df):
    """Cover pos in trimming zones (depends on bs_strand, flen...) with NA

    takes approx. 2:15 min
    """

    print("Masking dataframe")

    def mask_trimming_zones(group_df, cutting_sites_df):
        bs_strand, flen, pos = group_df.name
        left, right = cutting_sites_df.loc[
            (bs_strand, flen, ['left_cut_end', 'right_cut_end']), 'cut_pos']
        if left <= pos <= right:
            return True
        else:
            return False

    return (mbias_stats_df
            .groupby(['bs_strand', 'flen', 'pos'])
            .filter(mask_trimming_zones, dropna=False,
                    cutting_sites_df=cutting_sites_df)
            )


def pos_vs_beta_plots(mbias_stats_dfs_dict, config, aes_mappings):
    trunk_path = config['paths']['mbias_plots_trunk']
    for (curr_name, curr_df), curr_aes_mapping in product(
            mbias_stats_dfs_dict.items(), aes_mappings):

        print(f"Creating plot for dataframe {curr_name}, with aesthetics: "
              f"{curr_aes_mapping}")

        groupby_vars = (
            [val for val in curr_aes_mapping.values() if val is not None]
            + ['pos'])

        if 'flen' in groupby_vars:
            flen_sel = config['plots']['mbias_flens_to_display']
            curr_df = curr_df.loc[idxs[:, :, :, flen_sel], :]

        plot_df = (curr_df
                   .groupby(level=groupby_vars)
                   .sum()
                   .assign(beta_value=compute_beta_values)
                   .reset_index()
                   )

        p = (sns.FacetGrid(plot_df, **curr_aes_mapping)
             .map(plt.plot, 'pos', 'beta_value', alpha=0.3)
             .set(ylim = (.3, 1))
             .add_legend()
             )
        strat_name = '_'.join([f"{aes}-{var}"
                               for aes, var in curr_aes_mapping.items()])
        p.savefig(f"{trunk_path}_{curr_name}_{strat_name}.png")


def freq_over_pos_plot_per_motif(mbias_stats_dfs_dict, config):
    trunk_path = config['paths']['cg_occurence_plot_trunk']

    # this must always contain the meth_status, otherwise it will
    # affect the sum aggregation below
    aes_mappings = [
        {'row': 'bs_strand', 'col': 'meth_status', 'hue': None},
        {'row': 'bs_strand', 'col': 'meth_status', 'hue': 'flen'},
        {'row': 'bs_strand', 'col': 'meth_status', 'hue': 'phred'},
    ]

    for (curr_name, curr_df), curr_aes_mapping in product(
            mbias_stats_dfs_dict.items(), aes_mappings):

        print(f"Creating plot for dataframe {curr_name}, with aesthetics: "
              f"{curr_aes_mapping}")

        curr_df.columns.name = 'meth_status'

        groupby_vars = (
            ['motif']
            + [val for val in curr_aes_mapping.values()
               # meth status level does not exist yet, still in columns
               if val is not None and val != 'meth_status']
            + ['pos'])

        print('Filtering for displayed levels')
        if 'flen' in groupby_vars:
            flens = config['plots']['mbias_flens_to_display']
            curr_df = curr_df.loc[idxs[:, :, :, flens], :]
        # if 'phred' in groupby_vars:
        #     phreds = config['plots']['mbias_phreds_to_display']
        #     curr_df = curr_df.loc[idxs[:, :, :, :, phreds], :]

        print('Aggregating unused levels')
        curr_df = curr_df.groupby(level=groupby_vars).sum()

        print('Creating long format DataFrame')
        curr_df = (curr_df[['n_meth', 'n_unmeth']]
                   .stack(dropna=False)
                   .to_frame('counts'))

        for curr_motif, curr_motif_df in curr_df.groupby(level='motif'):

            plot_df = curr_motif_df.reset_index()

            print('Plotting')
            p = sns.FacetGrid(plot_df, sharey=False, margin_titles=True,
                              legend_out=True, **curr_aes_mapping)
            p.map(plt.plot, 'pos', 'counts')
            p.add_legend()
            p.fig.tight_layout()

            strat_name = '_'.join(
                [f"{aes}-{var}" for aes, var in curr_aes_mapping.items()])
            p.savefig(
                f"{trunk_path}_{curr_name}_{strat_name}_{curr_motif}.png")


def create_freq_with_agg_pos_plots(mbias_stats_df, config):
    trunk_path = config['paths']['freq_with_agg_pos_plot']

    mbias_stats_df.columns.name = "meth_status"
    plot_df = (mbias_stats_df
               .loc[:, ["n_meth", "n_unmeth"]]
               .groupby(['bs_strand', 'phred'])
               .sum()
               .stack()
               .to_frame("counts")
               .reset_index())

    g = sns.factorplot("phred", "counts",
                       row="bs_strand", col="meth_status", data=plot_df)
    out_fp = trunk_path + ".png"
    g.savefig(out_fp)


def seq_context_beta_value_plots(mbias_stats_df, config):
    trunk_path = config['paths']['freq_with_agg_pos_plot']

    plot_df = (mbias_stats_df
               .loc['CG', :]
               .groupby(['seq_context', 'bs_strand'])
               .sum()
               .assign(beta_value=compute_beta_values)
               .reset_index()
               )

    plot_df['seq_context'].cat.remove_unused_categories(inplace=True)


    g = (gg.ggplot(data=plot_df, mapping=gg.aes(
        x = "seq_context",
        y = "beta_value",
        color="bs_strand",
        group="bs_strand"))
         + gg.geom_point()
         + gg.geom_line()
         + gg.theme(axis_text_x = gg.element_text(angle=90, size=8))
         )
    out_fp = trunk_path + "seq-context-beta_value-plot.png"
    g.save(out_fp, height=12, width= 28, units='cm')

    # g = sns.factorplot('seq_context', 'beta_value', hue='bs_strand',
    #                    data=plot_df)
    # g.savefig(out_fp)

def seq_context_beta_value_plots_by_phred(mbias_stats_df, config):
    trunk_path = config['paths']['freq_with_agg_pos_plot']

    df1 = (mbias_stats_df
           .loc[idxs['CG', :, ['c_bc', 'c_bc_rv']], :]
           .groupby(['seq_context', 'bs_strand', 'phred'])
           .sum()
           .assign(beta_value=compute_beta_values)
           )

    df2 = (df1
           .drop(['n_meth', 'n_unmeth'], axis=1)
           .sort_index()
           .groupby(['seq_context', 'phred'])
           .agg(lambda ser: ser.iloc[0] - ser.iloc[1])
           .reset_index()
           )

    df2['seq_context'].cat.remove_unused_categories(inplace=True)

    g = (gg.ggplot(data=df2, mapping=gg.aes(
        x = "seq_context",
        y = "beta_value",
        color="phred",
        group="phred"))
         + gg.geom_point()
         + gg.geom_line()
         + gg.theme(axis_text_x = gg.element_text(angle=90, size=8))
         )
    out_fp = trunk_path + "seq-context-beta_value-plot-with-phred-groups.png"
    g.save(out_fp, height=12, width= 28, units='cm')

    def get_filtering_result(ser):
        cum_ser = ser.cumsum()
        return cum_ser[-1] - cum_ser
    df1 = (mbias_stats_df
           .loc[idxs['CG', :, ['c_bc', 'c_bc_rv']], ['n_meth', 'n_unmeth']]
           .groupby(['seq_context', 'bs_strand', 'phred'])
           .sum()
           .groupby(['seq_context', 'bs_strand'])
           .transform(get_filtering_result)
           .assign(beta_value=compute_beta_values)
           )

    df2 = (df1
           .drop(['n_meth', 'n_unmeth'], axis=1)
           .sort_index()
           .groupby(['seq_context', 'phred'])
           .agg(lambda ser: ser.iloc[0] - ser.iloc[1])
           .reset_index()
           )

    df2['seq_context'].cat.remove_unused_categories(inplace=True)

    g = (gg.ggplot(data=df2, mapping=gg.aes(
        x = "seq_context",
        y = "beta_value",
        color="phred",
        group="phred"))
         + gg.geom_point()
         + gg.geom_line()
         + gg.theme(axis_text_x = gg.element_text(angle=90, size=8))
         )

    out_fp = trunk_path + "seq-context-beta_value-plot-with-phred-filtering.png"
    g.save(out_fp, height=12, width= 28, units='cm')

    def get_filtering_result(ser):
        cum_ser = ser.cumsum()
        return cum_ser[-1] - cum_ser
    df1 = (mbias_stats_df
           .loc[idxs['CG', :, ['c_bc', 'c_bc_rv']], ['n_meth', 'n_unmeth']]
           .groupby(['seq_context', 'phred'])
           .sum()
           .groupby(['seq_context'])
           .transform(get_filtering_result)
           .assign(beta_value=compute_beta_values)
           .reset_index()
           )

    df1['seq_context'].cat.remove_unused_categories(inplace=True)

    g = (gg.ggplot(data=df1, mapping=gg.aes(
        x = "seq_context",
        y = "beta_value",
        color="phred",
        group="phred"))
         + gg.geom_point()
         + gg.geom_line()
         + gg.theme(axis_text_x = gg.element_text(angle=90, size=8))
         )

    out_fp = trunk_path + "seq-context-beta_value-plot-for-diff-ptresh.png"
    g.save(out_fp, height=12, width= 28, units='cm')


def freq_over_seq_context(mbias_stats_df, config):

    mbias_stats_df.columns.name = "meth_status"
    df1 = (mbias_stats_df
           .drop(["beta_value"], axis=1)
           .groupby(['bs_strand'])
           .sum()
           .stack()
           .to_frame("counts")
           )
    (df1
     .sort_index()
     .groupby("bs_strand")
     .agg(lambda ser: ser.iloc[0] / ser.iloc[1])
     )

    df1 = (mbias_stats_df
           .loc[idxs['CG', :, ['c_bc', 'c_bc_rv']], ['n_meth', 'n_unmeth']]
           .groupby(['bs_strand', 'seq_context'])
           .sum()
           .reset_index()
           )
    df1['seq_context'].cat.remove_unused_categories(inplace=True)


def create_mbias_stats_plots(mbias_stats_dfs_dict, config):
    print("Generating plots")
    aes_mappings = [
        # {'row': 'motif', 'col': 'bs_strand', 'hue': None},
        # {'row': 'motif', 'col': 'bs_strand', 'hue': 'flen'},
        # {'row': 'motif', 'col': 'bs_strand', 'hue': 'phred'},
        {'row': 'phred', 'col': 'bs_strand', 'hue': 'seq_context'},
    ]
    pos_vs_beta_plots(mbias_stats_dfs_dict, config, aes_mappings)  # ~ 1 min
    freq_over_pos_plot_per_motif(mbias_stats_dfs_dict, config)
    create_phred_filtering_mbias_plots(
        mbias_stats_dfs_dict['full'], config)  # ~ 1 min


def create_phred_filtering_mbias_plots(mbias_stats_df, config):
    phred_filter_dfs = compute_phred_filtering_dfs(
        mbias_stats_df)

    print('Plot M-bias plots for different phred filtering thresholds')
    plot_data = phred_filter_dfs['counts_by_phred_filtering'].reset_index()
    g = (sns.FacetGrid(plot_data, row='motif', col='bs_strand', hue='phred')
         .map(plt.plot, 'pos', 'beta_value')
         .add_legend())
    g.savefig(config['paths'][
                  'mbias_plots_trunk'] + '_phred_filtering_effect_by_pos.png')

    print('Plot per strand and global meth. levels')
    plot_data = phred_filter_dfs[
        'counts_by_phred_filtering_agg_pos'].reset_index()
    g = sns.factorplot(x='bs_strand', y='beta_value', hue='phred',
                       data=plot_data)
    g.savefig(config['paths']['mbias_plots_trunk']
              + '_phred_filtering_effect_by_strand.png')

    plot_data = phred_filter_dfs[
        'counts_by_phred_filtering_global'].reset_index()
    g = sns.factorplot(x='phred', y='beta_value', data=plot_data)
    g.savefig(
        config['paths']['mbias_plots_trunk'] + '_phred_filtering_global.png')


def compute_phred_filtering_dfs(mbias_stats_df):
    res = {}

    print('Calculating cumulative phred counts')
    res['cum_counts_by_phred'] = (
        mbias_stats_df
            .groupby(['motif', 'bs_strand', 'phred', 'pos'])
            .sum()
        # group_keys=False prevents index doubling by expanding().sum()
            .groupby(['motif', 'bs_strand', 'pos'], group_keys=False)
            .expanding(1)
            .sum()
            .sort_index()
    )

    print('Calculating counts for simulated phred filtering thresholds')
    res['counts_by_phred_filtering'] = (
        res['cum_counts_by_phred']
            .groupby(['motif', 'bs_strand', 'pos'])
            .transform(lambda ser: -1 * ser + ser.iloc[-1])
            .assign(beta_value=compute_beta_values)
            .sort_index()
    )

    print('Calculationg per strand statistics (aggregate pos)')
    res['counts_by_phred_filtering_agg_pos'] = (
        res['counts_by_phred_filtering']
            .groupby(['motif', 'bs_strand', 'phred'])
            .sum()
            .assign(beta_value=compute_beta_values))

    print('Calculating global statistics')
    res['counts_by_phred_filtering_global'] = (
        res['counts_by_phred_filtering_agg_pos']
            .groupby('phred')
            .sum()
            .assign(beta_value=compute_beta_values))

    return res


def create_tile_plots(mbias_stats_df, config):
    motif = 'CG'
    displayed_flens = config['plots']['mbias_flens_to_display']
    displayed_phreds = config['plots']['mbias_phreds_to_display']
    aes_mapping = {
        'x': 'phred',
        'y': 'seq_context',
        'fill': 'beta_value',
        # 'label': 'cell_label',
        'col': 'bs_strand',
        # 'row': 'flen',
    }
    # facetting_tuple = (aes_mapping.pop('row'), aes_mapping.pop('col'))
    facet_var = aes_mapping.pop('col')

    motif_df = (mbias_stats_df
                .loc[idxs[motif, :, :, displayed_flens, displayed_phreds], :]
                #         .groupby(
                # ['motif', 'seq_context', 'bs_strand', 'flen', 'phred'])
                .groupby(
        ['motif', 'seq_context', 'bs_strand', 'phred'])
                .sum()
                .assign(beta_value=compute_beta_values)
                .reset_index()
                )

    motif_df['seq_context'].cat.remove_unused_categories(inplace=True)
    motif_df['n_total'] = motif_df['n_meth'] + motif_df['n_unmeth']
    motif_df['beta_value_label'] = motif_df['beta_value'].apply(
        lambda fl: f"{fl:.2f}")
    motif_df['cell_label'] = motif_df[['beta_value', 'n_total']].apply(
        lambda ser: f"{ser.beta_value:.2f}\n({(ser.n_total/1000):.1f})",
        axis=1)

    g = (gg.ggplot(data=motif_df, mapping=gg.aes(**aes_mapping))
         # + gg.facet_grid(facetting_tuple)
         + gg.facet_wrap(facet_var)
         + gg.geom_tile()
         # + gg.geom_text(size=9, )
         )

    trunk_path = config['paths']['mbias_stats_heatmap_trunk']
    aes_mapping_str = '_'.join(
        [f"{key}-{value}" for key, value in aes_mapping.items()])
    target_path = f"{trunk_path}_{aes_mapping_str}.png"
    g.save(target_path, height=15, width=15, units='in')


def seq_context_meth_levels(mbias_stats_df, config):
    aes_mapping = {
        'x': 'seq_context',
        'y': 'beta_value',
        'fill': 'bs_strand'
    }
    plot_df = (mbias_stats_df
               .groupby(['seq_context', 'bs_strand'])
               .sum()
               .assign(beta_value=compute_beta_values)
               .reset_index()
               )
    g = (gg.ggplot(data=plot_df, mapping=gg.aes(**aes_mapping))
         + gg.geom_point())

    trunk_path = config['paths']['mbias_stats_heatmap_trunk']  # wrong!
    target_path = f"{trunk_path}_test.png"
    g.save(target_path)


"""
For testing analyze_mbias_counts:

from mqc.config import assemble_config_vars
from mqc.utils import get_resource_abspath
cli_params = {'motifs_str': 'CG',
              'sample_name': 'sample1',
              'sample_meta': None,
              'output_dir': (
                  "/icgc/dkfzlsdf/analysis/B080/kraemers/projects/mbias"
                  "/results_per_sample/prostate/wgbs/"
                  "BPH203_st-normal-luminal_se-x_me-wgbs_sw-34038/bseqtools"),
              'no_cache': False}
default_config_file = get_resource_abspath('config.default.toml')
config = assemble_config_vars(cli_params,
                              default_config_file_path=default_config_file)
config['plots']['mbias_flens_to_display'] = [60, 75, 90, 100, 125, 150, 171, 191, 251, 311, 411]
# config['paths']['mbias_counts'] = (
#    "/icgc/dkfzlsdf/analysis/B080/kraemers/projects/mbias"
#    "/results_per_sample/prostate/wgbs/"
#    "BPH171_st-normal-basal_se-x_me-wgbs_sw-34038/bseqtools/"
#    "qc_stats/BPH171_st-normal-basal_se-x_me-wgbs_sw-34038_mbias-counts_CG")
config['paths']['mbias_counts'] = (
    "/icgc/dkfzlsdf/analysis/B080/kraemers/projects/mbias/results_per_sample/"
    "prostate/wgbs/BPH203_st-normal-luminal_se-x_me-wgbs_sw-34038/bseqtools/"
    "qc_stats/BPH203_st-normal-luminal_se-x_me-wgbs_sw-34038_mbias-counts_CG")
"""


def analyze_mbias_counts(config):
    """ Run standard analysis on Mbias-Stats

    Parameters
    ----------
    config: dict
        config['run']['no_cache'] -> don't read previously computed
         mbias stats df from disc
    """

    fps = config['paths']

    # compute and save mbias stats should return the dict with all mbias stats dfs
    if config['run']['no_cache']:
        print('Not looking for cached results, directly computing stats')
        mbias_stats_df, masked_mbias_stats_df = compute_and_save_mbias_stats(config)
    else:
        print('Looking for M-bias stats data on disc')
        try:
            mbias_stats_df = pd.read_pickle(fps['mbias_stats_p'])
            masked_mbias_stats_df = pd.read_pickle(fps['mbias_stats_masked_p'])
            print('Data found and loaded')
        except FileNotFoundError:
            print('No data found, computing them now')
            mbias_stats_df, masked_mbias_stats_df = compute_and_save_mbias_stats(config)

    stat_dfs_dict = {'full': mbias_stats_df,
                     'trimmed': masked_mbias_stats_df}

    print("Adding phred filtering info")
    stat_dfs_dict['full_flen-agg_cg_phred-threshold'] = convert_phred_bins_to_thresholds(
        mbias_stats_df)

    stat_dfs_dict['trimmed_flen-agg_cg_phred-threshold'] = (
        convert_phred_bins_to_thresholds(masked_mbias_stats_df))

    print('Creating plots')

    create_mbias_stats_plots(stat_dfs_dict, config)


def convert_phred_bins_to_thresholds(mbias_stats_df):
    """Takes full Mbias stats df and returns CG only, flen agg df
    with phred bins converted to phred thresholds

    The operation would take prohibitively long on the entire dataframe
    (> 500 min ?)

    Important
    ---------
    This function assumes that mbias_stats_df is sorted
    """


    # *This function assumes that mbias_stats_df is sorted*

    non_flen_levels = [n for n in mbias_stats_df.index.names if n != "flen"]
    non_flen_phred_levels = [n for n in non_flen_levels if n != "phred"]

    mbias_stats_df_cg = mbias_stats_df.loc[idxs['CG':'CG'], :]
    mbias_stats_df_cg_flen_agg = (mbias_stats_df_cg
                                  .drop(['beta_value'], axis=1)
                                  .groupby(non_flen_levels)
                                  .sum())


    def compute_phred_threshold_counts_for_group(ser):
        """Will work on n_meth and n_unmeth"""
        cum_ser = ser.cumsum()
        return cum_ser[-1] - cum_ser
    res = (mbias_stats_df_cg_flen_agg
           .groupby(non_flen_phred_levels)
           .transform(compute_phred_threshold_counts_for_group)
           .assign(beta_value=compute_beta_values)
           )

    return res


def compute_and_save_mbias_stats(config):
    """Compute different M-bias stats

    Computes and saves
    - Extended Mbias-Stats DataFrame
    - Minimal Mbias-Stats DataFrame
    - bs_strand and flen adjusted cutting sites
    - various QC plots

    Returns
    -------

    """
    # TODO: switch to trunk paths in config file
    # TODO: switch to integer labels in MbiasCounterfps = config['paths']
    fps = config['paths']
    os.makedirs(fps['qc_stats_dir'], exist_ok=True, mode=0o770)

    mbias_stats_df = compute_mbias_stats_df(fps['mbias_counts'] + '.p')

    save_df_to_trunk_path(mbias_stats_df, fps['mbias_stats_p'].replace('.p', ''))

    classic_mbias_stats_df = compute_classic_mbias_stats_df(mbias_stats_df)
    save_df_to_trunk_path(classic_mbias_stats_df, fps['mbias_stats_classic_p'].replace('.p', ''))

    adjusted_cutting_sites = AdjustedCuttingSites(
        classic_mbias_stats_df, config)

    with open(fps['adjusted_cutting_sites_obj_p'], 'wb') as fobj:
        pickle.dump(adjusted_cutting_sites, fobj)
    adjusted_cutting_sites_df = adjusted_cutting_sites.get_df()
    save_df_to_trunk_path(adjusted_cutting_sites_df, fps["adjusted_cutting_sites_df_p"].replace('.p', ''))
    # cutting_sites_plot(adjusted_cutting_sites_df, config)

    masked_mbias_stats_df = mask_mbias_stats_df(
        mbias_stats_df, adjusted_cutting_sites_df)
    save_df_to_trunk_path(masked_mbias_stats_df, fps['mbias_stats_masked_p'].replace('.p', ''))

    return mbias_stats_df, masked_mbias_stats_df


def compute_beta_values(df):
    return df['n_meth'] / (df['n_meth'] + df['n_unmeth'])


def save_df_to_trunk_path(df, trunk_path):
    print("Saving pickle")
    df.to_pickle(trunk_path + '.p')
    print("Saving feather")
    df.reset_index().to_feather(trunk_path + '.feather')
