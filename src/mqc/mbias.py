import os
import pickle
import warnings

from abc import ABCMeta, abstractmethod
from itertools import product
from math import floor, ceil
from typing import Dict

import numpy as np
import pandas as pd

import matplotlib

matplotlib.use('Agg')  # import before pyplot import!
import matplotlib.pyplot as plt
import seaborn as sns

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
    index_cols = ['motif', 'seq_context', 'bs_strand', 'flen', 'phred', 'pos',
                  'meth_status']
    mbias_stats_df = (pd.read_pickle(mbias_counter_fp_str)
                      .get_dataframe()
                      # for interactive testing
                      # .query("seq_context in "
                      #        "['WWCGW', 'CCCGC', 'GGCGG', 'CGCGC', 'GCCGG']")
                      .groupby(['flen', 'pos'])
                      .filter(lambda group_df: (group_df.name[0].right - 1
                                                >= group_df.name[1]))
                      .assign(motif=lambda df:
    df['seq_context'].apply(map_seq_ctx_to_motif))
                      .set_index(index_cols)
                      # do not use .loc[:, 'counts'] to convert to series
                      # then you can't access group keys via group_df.name
                      .loc[:, 'counts']
                      .unstack('meth_status')
                      )
    # columns is categorical index, can't be extended
    # mbias_stats_df.columns = mbias_stats_df.columns.tolist()
    mbias_stats_df.columns = ['n_meth', 'n_unmeth']
    mbias_stats_df['beta_value'] = compute_beta_values(mbias_stats_df)
    return mbias_stats_df


def compute_classic_mbias_stats_df(mbias_stats_df):
    return (mbias_stats_df
        .groupby(level=['motif', 'bs_strand', 'flen', 'pos'])
        .sum()
        .groupby(['flen', 'pos'])
        .filter(lambda group_df: (group_df.name[0].right - 1
                                  >= group_df.name[1]))
        .assign(beta_value=(
        lambda df: df['n_meth'] / (df['n_unmeth'] + df['n_meth'])))
    )


def mask_mbias_stats_df(df: pd.DataFrame, cutting_sites_df: pd.DataFrame):
    def mask_positions(group_df):
        bs_strand, flen = group_df.name
        left, right = cutting_sites_df.loc[
            (bs_strand, flen, ['left_cut_end', 'right_cut_end']), 'cut_pos']
        group_df.loc[idxs[:, bs_strand, flen, 1:(left - 1)], :] = np.nan
        group_df.loc[idxs[:, bs_strand, flen, (right + 1):], :] = np.nan
        return group_df

    return (df.groupby(level=['bs_strand', 'flen'])
            .apply(mask_positions))


def mask_mbias_stats_df(mbias_stats_df, cutting_sites_df):
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
             .map(plt.plot, 'pos', 'beta_value')
             .add_legend()
             )
        strat_name = '_'.join([f"{aes}-{var}"
                               for aes, var in curr_aes_mapping.items()])
        p.savefig(f"{trunk_path}_{curr_name}_{strat_name}.png")


def freq_plot_per_motif(mbias_stats_dfs_dict, config):
    trunk_path = config['paths']['cg_occurence_plot_trunk']

    aes_mappings = [
        {'row': 'bs_strand', 'col': 'meth_status', 'hue': None},
        {'row': 'bs_strand', 'col': 'meth_status', 'hue': 'flen'},
        # {'row': 'bs_strand', 'col': 'meth_status', 'hue': 'Phred'},
    ]

    for (curr_name, curr_df), curr_aes_mapping in product(
            mbias_stats_dfs_dict.items(), aes_mappings):

        curr_df = (curr_df[['n_meth', 'n_unmeth']]
                   .stack(dropna=False)
                   .to_frame('counts'))

        for curr_motif, curr_motif_df in curr_df.groupby(level='motif'):

            groupby_vars = (['motif']
                            + [val for val in curr_aes_mapping.values() if
                               val is not None]
                            + ['pos'])
            agg_df = curr_motif_df.groupby(level=groupby_vars).sum()

            plot_df = agg_df.reset_index()

            if 'flen' in groupby_vars:
                row_is_displayed = plot_df['flen'].isin(
                    config['plots']['mbias_flens_to_display'])
                plot_df = plot_df.loc[row_is_displayed, :]
            if 'phred' in groupby_vars:
                row_is_displayed = plot_df['phred'].isin(
                    config['plots']['mbias_phreds_to_display'])
                plot_df = plot_df.loc[row_is_displayed, :]

            p = sns.FacetGrid(plot_df, sharey=False, margin_titles=True,
                              legend_out=True, **curr_aes_mapping)
            p.map(plt.plot, 'pos', 'counts')
            p.add_legend()
            p.fig.tight_layout()

            strat_name = '_'.join(
                [f"{aes}-{var}" for aes, var in curr_aes_mapping.items()])
            p.savefig(
                f"{trunk_path}_{curr_name}_{strat_name}_{curr_motif}.png")


def create_mbias_stats_plots(mbias_stats_dfs_dict, config):
    aes_mappings = [
        {'row': 'motif', 'col': 'bs_strand', 'hue': None},
        {'row': 'motif', 'col': 'bs_strand', 'hue': 'flen'},
        {'row': 'motif', 'col': 'bs_strand', 'hue': 'phred'},
    ]
    pos_vs_beta_plots(mbias_stats_dfs_dict, config, aes_mappings)
    # freq_plot_per_motif(mbias_stats_dfs_dict, config)
    create_phred_filtering_mbias_plots(
        mbias_stats_dfs_dict['full'], config)


def create_phred_filtering_mbias_plots(mbias_stats_df, config):
    phred_filter_dfs = compute_phred_filtering_dfs(
        mbias_stats_df)

    plot_data = phred_filter_dfs['counts_by_phred_filtering'].reset_index()
    g = (sns.FacetGrid(plot_data, row='motif', col='bs_strand', hue='phred')
         .map(plt.plot, 'pos', 'beta_value')
         .add_legend())
    g.savefig(config['paths'][
                  'mbias_plots_trunk'] + '_phred_filtering_effect_by_pos.png')

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

    res['cum_counts_by_phred'] = (
        mbias_stats_df
            .groupby(['motif', 'bs_strand', 'phred', 'pos'])
            .sum()
            .groupby(['motif', 'bs_strand', 'pos'], group_keys=False)
            .expanding(1)
            .sum()
            .sort_index()
    )

    res['counts_by_phred_filtering'] = (
        res['cum_counts_by_phred']
            .groupby(['motif', 'bs_strand', 'pos'])
            .transform(lambda ser: -1 * ser + ser.iloc[-1])
            .assign(beta_value=compute_beta_values)
            .sort_index()
    )

    res['counts_by_phred_filtering_agg_pos'] = (
        res['counts_by_phred_filtering']
            .groupby(['motif', 'bs_strand', 'phred'])
            .sum()
            .assign(beta_value=compute_beta_values))

    res['counts_by_phred_filtering_global'] = (
        res['counts_by_phred_filtering_agg_pos']
            .groupby('phred')
            .sum()
            .assign(beta_value=compute_beta_values))

    return res


def create_tile_plots(mbias_stats_df, config):
    motif = 'CG'
    # displayed_flens = [80, 100, 120, 150, 200, 300, 400]
    displayed_flens = [120, 200]
    displayed_phreds = [10, 20, 30]
    aes_mapping = {
        'x': 'phred',
        'y': 'seq_context',
        'col': 'bs_strand',
        'row': 'flen',
        'fill': 'beta_value',
        'label': 'cell_label',
    }
    facetting_tuple = (aes_mapping.pop('row'), aes_mapping.pop('col'))

    trunk_path = config['paths']['mbias_stats_heatmap_trunk']
    aes_mapping_str = '_'.join(
        [f"{key}-{value}" for key, value in aes_mapping.items()])
    target_path = f"{trunk_path}_{aes_mapping_str}.png"

    motif_df = (mbias_stats_df
                .loc[idxs[motif, :, :, displayed_flens, displayed_phreds], :]
                .groupby(
        ['motif', 'seq_context', 'bs_strand', 'flen', 'phred'])
                .sum()
                .assign(beta_value=lambda df:
    df['n_meth'] / (df['n_meth'] + df['n_unmeth']))
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
         + gg.facet_grid(facetting_tuple)
         + gg.geom_tile()
         + gg.geom_text(size=9, )
         )
    g.save(target_path, height=15, width=15, units='in')


def analyze_mbias_counts(config):
    os.makedirs(config['paths']['qc_stats_dir'],
                exist_ok=True, mode=0o770)
    mbias_stats_df = compute_mbias_stats_df(
        config['paths']['mbias_counts'] + '.p')

    mbias_stats_df.to_pickle(config['paths']['mbias_stats_p'])
    mbias_stats_df.to_csv(config['paths']['mbias_stats_tsv'],
                          header=True, index=True, sep='\t')

    classic_mbias_stats_df = compute_classic_mbias_stats_df(mbias_stats_df)
    classic_mbias_stats_df.to_pickle(config['paths']['mbias_stats_classic_p'])

    adjusted_cutting_sites = AdjustedCuttingSites(
        classic_mbias_stats_df, config)
    del classic_mbias_stats_df

    with open(config['paths']['adjusted_cutting_sites_obj_p'], 'wb') as fobj:
        pickle.dump(adjusted_cutting_sites, fobj)
    adjusted_cutting_sites_df = adjusted_cutting_sites.get_df()
    del adjusted_cutting_sites

    adjusted_cutting_sites_df.to_pickle(config["paths"]
                                        ["adjusted_cutting_sites_df_p"])
    adjusted_cutting_sites_df.to_csv(
        config["paths"]['adjusted_cutting_sites_df_tsv'],
        header=True, index=True, sep="\t")
    # cutting_sites_plot(adjusted_cutting_sites_df, config)

    masked_mbias_stats_df = mask_mbias_stats_df(
        mbias_stats_df, adjusted_cutting_sites_df)

    del adjusted_cutting_sites_df

    masked_mbias_stats_df.to_pickle(config['paths']['mbias_stats_masked_p'])
    masked_mbias_stats_df.to_csv(config['paths']['mbias_stats_masked_tsv'],
                                 sep='\t', header=True, index=True)

    mbias_stats_dfs_dict = {'full': mbias_stats_df,
                            'trimmed': masked_mbias_stats_df}

    create_mbias_stats_plots(mbias_stats_dfs_dict, config)


def compute_beta_values(df):
    return df['n_meth'] / (df['n_meth'] + df['n_unmeth'])
