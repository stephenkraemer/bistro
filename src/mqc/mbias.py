#-
import importlib
import itertools
import json
import os
import os.path as op
import pickle
import shelve
from pathlib import Path

import toolz as tz
from abc import ABCMeta, abstractmethod
from collections import namedtuple
from copy import deepcopy
from itertools import product
from math import floor, ceil
from typing import Dict, List, Union, Tuple, Optional, Sequence, Set
from dataclasses import dataclass

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use('Agg')  # import before pyplot import!
import matplotlib.pyplot as plt
import importlib
import seaborn as sns
import plotnine as gg
from plotnine import *
from sklearn import linear_model

import mqc.flag_and_index_values as mfl
# noinspection PyUnresolvedReferences
from mqc.pileup.bsseq_pileup_read import BSSeqPileupRead
from mqc.pileup.pileup import MotifPileup
from mqc.utils import convert_array_to_df, hash_dict, get_resource_abspath
from mqc.visitors import Counter
import mqc.filepaths
from mqc.utils import (update_nested_dict, NamedIndexSlice,
                       assert_match_between_variables_and_index_levels,
                       subset_dict)
import more_itertools

import altair as alt

b_inds = mfl.bsseq_strand_indices
b_na_ind = mfl.bsseq_strand_na_index
m_flags = mfl.methylation_status_flags

idxs = pd.IndexSlice
nidxs = NamedIndexSlice

# required if plots are to be generated with ggplot2
# import rpy2
# import rpy2.robjects as ro
# # automatically convert arrays and dfs
# from rpy2.robjects.numpy2ri import numpy2ri
# rpy2.robjects.numpy2ri.activate()
# from rpy2.robjects import pandas2ri
# pandas2ri.activate()
# import rpy2.robjects.lib.ggplot2 as gg
#-

MBIAS_STATS_DIMENSION_PLOT_LABEL_MAPPING = dict(
    pos='Position',
    flen='Fragment length',
    beta_value='Beta value',
    counts='Frequency',
    bs_strand='BS-strand',
    seq_context='Sequence context',
    motif='Motif',
    phred='Phred score',
    phred_threshold='Phred >',
)


class MbiasCounter(Counter):
    """Counter for multidimensional M-bias stats

    Implementation notes
    --------------------

    This counter is pretty close to the maximal object size for serialization
    with older pickle protocols. If you deviate too much from the datasets
    this was tested on (unusually high read lengths or fragment lengths etc.)
    it may become too big. You would then need to switch to a suitable
    serialization protocol. This will likely also be the case if you want
    to add another dimension, or if you want to add 'snp' and 'reference'
    levels to the methylation status dimension.

    *counter_array indexing*
    - the fragment length dimension includes the length 0, so that length N has index N.
    - the read position indexes are zero-based, for better interaction with the C/cython parts of the program.

    *counter dataframe index levels*
    - seq context: 3 letter motifs, e.g. CWCGW (size is parameter) [categorical]
    - bs strands are labelled in lowercase, e.g. c_bc [categorical]
    - single fragment lengths are labelled with the flen as int
    - binned fragment lengths are labelled with the rightmost flen in the bin
    - phred scores are always binned, and are also labelled with the rightmost phred score in the bin
    - pos labels: 1-based (note that array indexing is 0-based)
    - meth. status labels: n_meth, n_unmeth [categorical]
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

        self.seq_ctx_idx_dict, self.binned_motif_to_index_dict = \
            get_sequence_context_to_array_index_table(self.seq_context_size)

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

        flen_bin_labels = (
            list(range(0, self.max_single_flen + 1))
            + list(range(
            self.max_single_flen + self.flen_bin_size, self.max_flen + 1, self.flen_bin_size)))
        if flen_bin_labels[-1] < self.max_flen:
            flen_bin_labels.append(self.max_flen)

        phred_bin_labels = list(
            range(self.phred_bin_size - 1, self.max_phred + 1, self.phred_bin_size))
        if phred_bin_labels[-1] < self.max_phred:
            phred_bin_labels.append(self.max_phred)

        dim_levels = [get_categorical(ordered_seq_contexts),
                      get_categorical(['c_bc', 'c_bc_rv', 'w_bc', 'w_bc_rv']),
                      flen_bin_labels,
                      phred_bin_labels,
                      range(1, self.max_read_length + 1),
                      get_categorical(['n_meth', 'n_unmeth'])]

        array_shape = [len(ordered_seq_contexts),
                       4,  # BSSeq-strands
                       len(flen_bin_labels),
                       len(phred_bin_labels),
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
            if tlen <= self.max_single_flen:
                tlen_idx = tlen
            elif tlen > self.max_flen:
                tlen_idx = -1
            else:
                tlen_idx = (self.max_single_flen
                            + ceil((tlen - self.max_single_flen)
                                   / self.flen_bin_size))

            phred = curr_read.baseq_at_pos
            if phred > self.max_phred:
                phred_idx = -1
            else:
                phred_idx = floor(phred / self.phred_bin_size)

            self.counter_array[seq_ctx_idx,
                               curr_read.bsseq_strand_ind,
                               tlen_idx,
                               phred_idx,
                               curr_read.pos_in_read,
                               meth_status_index] += 1


def get_sequence_context_to_array_index_table(motif_size: int):
    """ Return dicts mapping sequence contexts to counter array indices

    Parameters
    ----------
    motif_size: int
        size of motifs (e.g. 3 or 5 bp)

    Returns
    -------
    Dict[str, str]
        mapping of four letter motif to array integer index of the corresponding
        three letter motif. I.e. several motifs may map to the same index.
    Dict[str, int]
        mapping of three-letter motif [CGW] to array integer index.
        Every mapping is unique.
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

    all_bases_set = [all_bases] * n_bp_per_side + [['C']] + [all_bases] * n_bp_per_side
    all_5bp_motifs = [''.join(motif) for motif in product(*all_bases_set)]

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
        run times: [1.5 min, 9 min]
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
              # .apply(fit_normalvariate_plateau, config)
              .apply(fit_percentiles, config)
              )
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


def fit_percentiles(group_df: pd.DataFrame, config) -> pd.Series:

    min_perc = 0.5  # min_plateau_length = effective_read_length * min_perc
    percentiles = (0.02, 0.98)
    min_percentile_delta = 0.1
    min_flen = 45
    max_read_length = 101

    max_slope = min_percentile_delta/max_read_length

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
        percentile_delta_to_beat = min_percentile_delta
        best_start = None
        for start_pos in range(0, max_start_pos):
            end_pos = start_pos + plateau_length
            curr_beta_values = beta_values.iloc[start_pos:end_pos]
            low_percentile, high_percentile = curr_beta_values.quantile(percentiles)
            curr_percentile_delta = high_percentile - low_percentile
            if curr_percentile_delta < percentile_delta_to_beat:
                # plateau_height = curr_beta_values.mean()
                # left_end_ok = (curr_beta_values[0:4] > low_percentile).all()
                # right_end_ok = (curr_beta_values[-4:] < high_percentile).all()
                curr_beta_values_arr = curr_beta_values.values
                plateau_end_deltas = (curr_beta_values_arr[0:4, np.newaxis] -
                                      curr_beta_values_arr[np.newaxis, -4:])
                both_ends_ok = (plateau_end_deltas < min_percentile_delta).all()
                assert isinstance(both_ends_ok, np.bool_), type(both_ends_ok)
                if both_ends_ok:
                    regr = linear_model.LinearRegression()
                    X = np.arange(0, plateau_length)[:, np.newaxis]
                    Y = curr_beta_values_arr[:, np.newaxis]
                    regr.fit(X, Y)
                    if regr.coef_[0, 0] <= max_slope:
                        percentile_delta_to_beat = curr_percentile_delta
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

    print("Reading data in, removing impossible strata (pos > flen)")
    # convert MbiasCounter.counter_array to dataframe
    # remove positions > flen
    # takes approx. 2.5 min
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
    # Add beta value, approx. 1.5 min
    mbias_stats_df_with_meth = (mbias_stats_df
                                .loc[:, 'counts']
                                .unstack('meth_status'))
    # columns is categorical index,
    # replace with object index for easier extension
    mbias_stats_df_with_meth.columns = ['n_meth', 'n_unmeth']
    mbias_stats_df_with_meth['beta_value'] = compute_beta_values(
        mbias_stats_df_with_meth)

    print("Adding motif labels to index")
    # prepend motif level to index (~40s)
    mbias_stats_df_with_motif_idx = prepend_motif_level_to_index(
        mbias_stats_df_with_meth)

    print("Sorting")
    # sort - takes approx. 3-4.5 min with IntegerIndices for flen and phred
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

    takes between 2:15 min and 4:30 min
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


# def get_plotnine_poster_theme():
#     t = (theme_classic() +
#          theme(text = element_text(family = "DejaVu Sans", size = 50),
#                title = element_text(size = 40),
#                axis_title = element_text(size = 26),
#                axis_text = element_text(size = 22),
#                # plot_margin = 0.5,
#                strip_background = element_rect(color="white"),
#                strip_text = element_text(size=26),
#                axis_line = element_line(color="black"),
#                legend_title = element_text(size=26),
#                legend_text = element_text(size=22)
#                ))
#     return t

def get_plotnine_poster_theme():
    t = (theme_classic() +
         theme(text = element_text(family = "DejaVu Sans", size = 8),
               title = element_text(size = 40),
               axis_title = element_text(size = 24),
               axis_text = element_text(size = 20),
               # plot_margin = 0.5,
               strip_background = element_rect(color="white"),
               strip_text = element_text(size=24),
               axis_line = element_line(color="black"),
               legend_title = element_text(size=20),
               legend_text = element_text(size=20)
               ))
    return t

PLOTNINE_THEMES = dict(
    poster=(theme_classic() +
            theme(text = element_text(family = "DejaVu Sans", size = 8),
                  title = element_text(size = 40),
                  axis_title = element_text(size = 24),
                  axis_text = element_text(size = 20),
                  # plot_margin = 0.5,
                  strip_background = element_rect(color="white"),
                  strip_text = element_text(size=24),
                  axis_line = element_line(color="black"),
                  legend_title = element_text(size=20),
                  legend_text = element_text(size=20)
                  )),
)


def plot_mappings_for_dfs(mbias_stats_dfs_dict, df_names,
                          config, aes_mappings,
                          plot_vars={}, motifs=["CG"]):
    """

    Parameters
    ----------
    mbias_stats_dfs_dict: Dict
        contains 'mbias_stats_df' and 'masked_mbias_stats_df'
    df_names: List[str]
        mbias_stats_df_dict keys pointing to DFs to be used for plotting
    config: Dict
        Uses config['paths']['mbias_plots_trunk']
    aes_mappings: List[Dict]
        List of aesthetic mappings for seaborn FacetGrid,
        - required are "row", "col", "hue".
        - To leave out an aesthetic, assign None to it
        - You can also set the x and y aesthetics,
          which default to "pos" and "beta_value".
    alpha: int
        passed to ggplot
    motifs: Union[List[str], None]
        Restrict the mbias stats df to motifs, if not None

    Returns
    -------
    Plots as png, svg, pdf

    Notes
    -----

    - The y-axis limits are indicated in the filename, e.g. ylim-0-1
    - alpha for colors is hardcoded to 0.5

    ToDos
    -----
    - discard strata where coverage is too low
    - allowing creating one plot per motif using a "split" aesthetic or something similar
      (Then I could in principle also split by bs_strand or other aesthetics
    - allow restriction to a given motif
    """

    # Set defaults for x and y aesthetics
    default_mapping = {"x": "pos",
                       "y": "beta_value"}

    for (curr_name, curr_df), curr_aes_mapping in product(
            mbias_stats_dfs_dict.items(), aes_mappings):

        if curr_name not in df_names:
            continue

        # Motifs is List[str], so keeps index dimensions
        if motifs:
            curr_df = curr_df.loc[motifs, :]

        complete_aes = default_mapping.copy()
        complete_aes.update(curr_aes_mapping)

        mbias_plots_for_df_aes_combination(curr_df, curr_name, complete_aes,
                                           config, plot_vars)


def mbias_plots_for_df_aes_combination(mbias_stats_df, df_name,
                                       complete_aes, config, plot_vars,
                                       use_r_ggplot=True):

    print(f"Working on {df_name} with aes {list(complete_aes.items())}")

    flen_sel = config['plots']['mbias_flens_to_display']
    trunk_path = config['paths']['mbias_plots_trunk']

    plot_df_var_names = {"motif": "Motif", "seq_context": "Sequence context",
                         "mate": "Mate", "bs_strand": "BS-strand",
                         "flen": "Frag.\nlength", "phred": "Phred \u2265",
                         "pos": "Position", "beta_value": "Beta value",
                         "n_meth": "#methylated", "n_unmeth": "#unmethylated"}

    plot_df = calculate_plot_df(mbias_stats_df, complete_aes, flen_sel,
                                plot_df_var_names=plot_df_var_names)

    plot_vars = add_plot_aes_mappings(plot_vars, complete_aes, plot_df_var_names)

    plot_vars = add_derived_plot_vars(plot_vars, plot_df)

    for ylim_tuple in plot_vars["ylim_tuples"]:

        strat_name = '_'.join([f"{aes}-{var}"
                               for aes, var in complete_aes.items()])
        file_path_pattern = (f"{trunk_path}_{df_name}_{strat_name}"
                             f"_ylim-{ylim_tuple[0]}-{ylim_tuple[1]}"
                             ".{file_format}")  # note: file format not expanded here!

        plot_vars["ylim_tuple"] = ylim_tuple

        if use_r_ggplot:
            # plot_with_ggplot(plot_df, file_path_pattern, **plot_vars)
            raise NotImplemented
        else:
            plot_with_plotnine(plot_df, file_path_pattern, **plot_vars)


def add_plot_aes_mappings(plot_vars, complete_aes, plot_df_var_names):

    plot_vars["row_aes"] = complete_aes["row"]
    plot_vars["col_aes"] = complete_aes["col"]

    plot_aes_full = {k: v
                     for k, v in complete_aes.items()
                     if v is not None}

    plot_aes_no_grid_vars = {k: v for k, v in plot_aes_full.items()
                             if k not in ["row", "col"]}

    if "color" in plot_aes_no_grid_vars:
        plot_aes_no_grid_vars["group"] = plot_aes_no_grid_vars["color"]

    plot_vars["plot_aes"] = plot_aes_no_grid_vars

    plot_vars["xlab"] = plot_df_var_names[plot_aes_full["x"]]
    plot_vars["ylab"] = plot_df_var_names[plot_aes_full["y"]]
    try:
        plot_vars["legend_title"] = plot_df_var_names[plot_aes_full["color"]]
        # if plot_aes_full["color"] == "flen":
        #     plot_vars["legend_title"] = "Frag.\nlength"
    except KeyError:
        plot_vars["legend_title"] = "No color defined"

    return plot_vars


def add_derived_plot_vars(plot_vars, plot_df):

    pv = plot_vars
    row_aes = pv["row_aes"]
    col_aes = pv["col_aes"]
    # Generate plots in different formats and with different y-axis limits
    pv["need_row_facet"] = (bool(row_aes)
                            and (plot_df[row_aes].unique().shape[0] > 1))
    pv["need_col_facet"] = (bool(col_aes)
                            and (plot_df[col_aes].unique().shape[0] > 1))


    pv["plot_height"] = pv["panel_height"]
    pv["plot_width"] = pv["panel_width"]
    if row_aes:
        pv["plot_height"] = pv["panel_height"] * len(plot_df[row_aes].unique())
    if col_aes:
        pv["plot_width"] = pv["panel_width"] * len(plot_df[col_aes].unique())


    x_aes = pv["plot_aes"]["x"]
    pv["x_is_categorical"] = (plot_df[x_aes].dtype.name == "category")
    pv["min_x"] = plot_df[x_aes].min()
    pv["max_x"] = plot_df[x_aes].max()

    return pv


def plot_with_plotnine(plot_df, file_path_pattern,
                       plot_aes, col_aes, row_aes,
                       need_row_facet, need_col_facet,
                       x_is_categorical,
                       plot_height, plot_width,
                       ylim_tuple=[(0,1)],
                       y_breaks=list(np.arange(0, 1.01, 0.1)),
                       min_x=None, max_x=None,
                       x_breaks=None,
                       alpha=0.7,
                       rotate_x_labels=False,
                       theme=None,
                       ):


    if theme == "poster":
        ggtheme = get_plotnine_poster_theme()
    else:
        ggtheme = theme_classic()

    gp = (gg.ggplot(gg.aes(**plot_aes), plot_df)
          + gg.geom_line(alpha=alpha)
          + gg.scale_y_continuous(limits=ylim_tuple,
                                  breaks=y_breaks)
          + ggtheme
          )
    if x_is_categorical:
        gp += gg.scale_x_discrete()
    else:
        gp += gg.scale_x_continuous(limits=(min_x, max_x),
                                    breaks=x_breaks)

    if need_row_facet and need_col_facet:
        gp += gg.facet_grid((row_aes, col_aes))
    elif need_row_facet:
        gp += gg.facet_wrap(row_aes, ncol=1)
    elif need_col_facet:
        gp += gg.facet_wrap(col_aes, nrow=1)

    if x_is_categorical:
        gp += gg.geom_point(alpha=alpha)

    if rotate_x_labels:
        gp += gg.theme(axis_text_x=gg.element_text(angle=90,
                                                   hjust=0.5))

    for file_format in ["png", "svg", "pdf"]:
        gp.save(file_path_pattern.format(file_format=file_format),
                height=plot_height, width=plot_width, units="cm")


# def plot_with_ggplot(plot_df, file_path_pattern,
#                      plot_aes, col_aes, row_aes,
#                      need_row_facet, need_col_facet,
#                      x_is_categorical,
#                      plot_height, plot_width,
#                      xlab, ylab, legend_title=None,
#                      ylim_tuple=(0,1),
#                      y_breaks=list(np.arange(0, 1.01, 0.1)),
#                      min_x=None, max_x=None,
#                      x_breaks=None,
#                      alpha=0.7,
#                      rotate_x_labels=False,
#                      theme=None, **xargs):
#
#     # **xargs catches unused plot vars, which were used to compute other plot vars,
#     # e.g. panel_height. Would be cleaner to remove these redundant vars
#     # need to delete xargs now, can't pass it to R with other locals because
#     # xargs also contains ylim_tuples, see above
#     # it is a dict
#     del xargs
#
#     ylim_tuple = ro.vectors.FloatVector(ylim_tuple)
#     y_breaks = ro.vectors.FloatVector(y_breaks)
#     x_breaks = ro.vectors.FloatVector(x_breaks) if x_breaks else ro.NULL
#     plot_aes = {k: v.replace(" ", ".") for k,v in plot_aes.items()}
#     plot_aes = {k: v.replace("#", "") for k,v in plot_aes.items()}
#     plot_aes = ro.vectors.ListVector(plot_aes)
#
#     plot_vars = locals().copy()
#
#     # r_compatible_dtypes = plot_df.dtypes.replace({"uint64": "int64"}).to_dict()
#     # plot_df = plot_df.astype(r_compatible_dtypes)
#     # del r_compatible_dtypes
#
#     tmppath = f"/home/kraemers/temp/{op.basename(file_path_pattern)}.feather"
#     plot_df = plot_df.rename(columns=lambda s: s.replace(" ", "."))
#     plot_df = plot_df.rename(columns=lambda s: s.replace("#", ""))
#     plot_df.to_feather(tmppath)
#     plot_vars["plot_df_path"] = tmppath
#
#     for k,v in plot_vars.items():
#         if v is None:
#             v = ro.NULL
#         if k == "plot_df":
#             continue
#         ro.globalenv[k] = v
#
#
#     # ro.r("""\
#     # library(ggplot2)
#     # library(feather)
#     # library(glue)
#     # plot_df = read_feather(plot_df_path)
#     # print(head(plot_df))
#     # plot_theme = theme_classic()
#     # g = ggplot(plot_df, aes_string("Position", "Beta.value")) +
#     #      geom_line(alpha=0.7) +
#     #      scale_y_continuous(limits = c(0,1), breaks = c(0,1)) +
#     #      plot_theme
#     # ggsave("/home/kraemers/temp/plots/test.png")
#     # """)
#
#     ro.r("""\
#     library(ggplot2)
#     library(feather)
#     library(glue)
#
#     plot_df = read_feather(plot_df_path)
#
#     if ("flen" %in% colnames(plot_df)) {
#         unique_flens = unique(plot_df["flen"])
#         print(unique_flens)
#         n_unique_flens = nrow(unique_flens)
#         gg_color_hue <- function(n) {
#           hues = seq(15, 375, length = n + 1)
#           hcl(h = hues, l = 65, c = 100)[1:n]
#         }
#         cols = gg_color_hue(n_unique_flens)
#         names(cols) = as.character(unique_flens[["flen"]])
#         print(cols)
#         color_scale = scale_color_manual(
#           values = cols, guide = guide_legend(title = legend_title))
#     } else {
#         color_scale = scale_color_discrete(
#         guide = guide_legend(title = legend_title))
#     }
#
#
#
#     if (is.null(theme)) {
#       theme = "poster"
#     }
#
#     if (theme == "poster") {
#
#         plot_theme = theme_classic(base_family = "sans") +
#              theme(text = element_text(size = 40),
#                    title = element_text(size = 40),
#                    axis.title = element_text(size = 28),
#                    axis.text = element_text(color="black", size = 20),
#                    strip.background = element_rect(color="white"),
#                    strip.text = element_text(size=28),
#                    axis.line = element_line(color="black"),
#                    legend.title = element_text(size=24),
#                    legend.text = element_text(size=18),
#                    panel.spacing = unit(1, "lines")
#                    )
#
#     } else if (theme == "poster_small") {
#
#         plot_theme = theme_classic(base_family = "sans") +
#              theme(text = element_text(size = 40),
#                    title = element_text(size = 40),
#                    axis.title = element_text(size = 28),
#                    axis.text = element_text(color="black", size = 10),
#                    strip.background = element_rect(color="white"),
#                    strip.text = element_text(size=28),
#                    axis.line = element_line(color="black"),
#                    legend.title = element_text(size=24),
#                    legend.text = element_text(size=18),
#                    panel.spacing = unit(1, "lines")
#                    )
#     }
#
#     mapping = do.call(aes_string, plot_aes)
#
#     print(xlab)
#     print(ylab)
#
#     g = ggplot(plot_df, mapping) +
#       geom_line(alpha=0.7) +
#       scale_y_continuous(limits = ylim_tuple, breaks = y_breaks)  +
#       labs(x = xlab, y = ylab) +
#       color_scale +
#       # coord_flip() +
#       plot_theme
#     if (x_is_categorical) {
#       g = g + scale_x_discrete() +
#              geom_point(alpha=0.7)
#     } else {
#       g = g + scale_x_continuous(breaks = x_breaks)
#     }
#     print("Adding facets in next step")
#     if (need_row_facet && need_col_facet) {
#         g = g + facet_grid(paste(row_aes, "~", col_aes))
#     } else if (need_row_facet) {
#         g = g + facet_wrap(row_aes, ncol=1)
#     } else if (need_col_facet) {
#         g = g + facet_wrap(col_aes, nrow=1)
#     }
#     if (rotate_x_labels) {
#       g = g + theme(axis.text.x=element_text(angle=90, hjust=0.5, vjust=0.5))
#     }
#
#     print("Saving")
#     for (file_format in c("png", "svg", "pdf")) {
#       ggsave(filename = glue(file_path_pattern), plot = g,
#       height = plot_height, width = plot_width, units="cm")
#     }
#     file_format = "rds"
#     saveRDS(g, glue(file_path_pattern))
#     """)
#
#     # #
#     # # print(plot_aes)
#     # # print(class(plot_aes))
#     #
#     # # print("mapping")
#     # # print(mapping)
#     #
#     #
#     #
#     # print("base plot done")
#     #
#     #
#     # #
#     # #
#     # #
#     # #
#     # #
#     #
#     # a = 3
#     # print("Done")
#     # """)
#     #
#



def freq_over_pos_plot_per_motif(mbias_stats_dfs_dict, config):
    """

    Parameters
    ----------
    mbias_stats_dfs_dict
    config

    Returns
    -------

    ToDos
    -----

    Normalize frequencies

    """
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


# TODO: remove
# def create_mbias_stats_plots(mbias_stats_dfs_dict, config):
#     """
#
#
#     """
#
#
#     print("Generating pos vs. beta plots")
#     pos_vs_beta_aes_mappings = [
#         # {'row': 'motif', 'col': 'bs_strand', 'color': None},
#         # {'row': 'motif', 'col': 'bs_strand', 'color': 'flen'},
#         # {'row': 'motif', 'col': 'bs_strand', 'color': 'phred'},
#         # {'row': 'motif', 'col': 'bs_strand', 'color': 'seq_context'},
#         {'row': 'motif', 'col': 'mate', 'color': None},
#         {'row': 'motif', 'col': 'mate', 'color': 'flen'},
#         # {'row': 'motif', 'col': 'mate', 'color': 'phred'},
#         # {'row': 'motif', 'col': 'mate', 'color': 'seq_context'},
#     ]
#     plot_mappings_for_dfs(mbias_stats_dfs_dict,
#                           df_names=[
#                               'full',
#                               'trimmed'
#                                     ],
#                           aes_mappings=pos_vs_beta_aes_mappings,
#                           motifs=["CG"],
#                           config=config,
#                           plot_vars=dict(
#                               rotate_x_labels=True,
#                               panel_height=11, panel_width=9,
#                               theme="poster",
#                               ylim_tuples=[(0, 1), (0.3, 0.9)],
#                               y_breaks=list(np.arange(0,1.01, 0.2)),
#                               x_breaks=list(range(0, 151, 30)),
#                           )
#                           )  # ~ 4 min
#
#     # phred filtering plots
#     phred_filtering_aes_mappings = [
#         # {"row": "motif", "col": "bs_strand", "color": "phred"},
#         {"row": "motif", "col": "mate", "color": "phred"},
#     ]
#     plot_mappings_for_dfs(mbias_stats_dfs_dict,
#                           df_names=[
#                           # "trimmed_flen-agg_cg_phred-threshold",
#                           "full_flen-agg_cg_phred-threshold"
#                       ],
#                           aes_mappings=phred_filtering_aes_mappings,
#                           motifs=["CG"],
#                           config=config,
#                           plot_vars=dict(
#                               rotate_x_labels=True,
#                               panel_height=11, panel_width=9,
#                               theme="poster",
#                               ylim_tuples=[(0,1), (0.2, 0.9)],
#                               y_breaks = list(np.arange(0,1.01, 0.2)),
#                               x_breaks = list(range(0, 151, 30)),
#                           ))  # ~ 1 min
#
#     # 5pb motif plots
#     # needs motif selection
#     aes_mappings = [
#         # {"x": "seq_context", "y": "beta_value",
#         #  "row": None, "col": None, "color": "phred"},
#         {"x": "seq_context", "y": "beta_value",
#          "row": None, "col": "mate", "color": "phred"},
#         # {"x": "seq_context", "y": "beta_value",
#         #  "row": None, "col": "bs_strand", "color": "phred"},
#     ]
#     plot_mappings_for_dfs(mbias_stats_dfs_dict,
#                           df_names=[
#                               # "trimmed_flen-agg_cg_phred-threshold",
#                               "full_flen-agg_cg_phred-threshold"
#                           ],
#                           aes_mappings=aes_mappings,
#                           motifs=["CG"],
#                           config=config,
#                           plot_vars=dict(
#                               rotate_x_labels=False,
#                               panel_height=11, panel_width=9,
#                               theme="poster_small",
#                               ylim_tuples=[(0,1), (0.3, 0.9)],
#                               y_breaks = list(np.arange(0,1.01, 0.1)),
#                               x_breaks = None,
#                           ))  # ~ 1 min
#
#     # create_phred_filtering_mbias_plots(
#     #     mbias_stats_dfs_dict['full'], config)  # ~ 1 min
#
#     # freq_over_pos_plot_per_motif(mbias_stats_dfs_dict, config)


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
For testing compute_mbias_stats:

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


def compute_mbias_stats(config):
    """ Run standard analysis on Mbias-Stats

    Parameters
    ----------
    config: dict
        config['run']['no_cache'] -> don't read previously computed
         mbias stats df from disc
    """

    # TODO: get phred threshold dfs from cache
    fps = config['paths']

    os.makedirs(fps['qc_stats_dir'], exist_ok=True, mode=0o770)

    print("Computing M-bias stats from M-bias counter array")
    mbias_stats_df = compute_mbias_stats_df(fps['mbias_counts'] + '.p')
    mbias_stats_df = add_mate_info(mbias_stats_df)
    mbias_stats_df = mbias_stats_df.sort_index()

    # TODO: Important: fix hardcoding
    mbias_stats_df = mbias_stats_df.loc[idxs[:, :, :, :, :, :, 1:150], :]

    # discard phreds which are not present in bins
    print('Discarding unused phred scores')
    n_total = mbias_stats_df["n_meth"] + mbias_stats_df["n_unmeth"]
    phred_group_sizes = n_total.groupby("phred").sum()
    phred_bin_has_counts = (phred_group_sizes > 0)
    existing_phred_scores = phred_group_sizes.index.values[phred_bin_has_counts]
    mbias_stats_df = mbias_stats_df.loc[idxs[:, :, :, :, :, existing_phred_scores], :]
    mbias_stats_df.index = mbias_stats_df.index.remove_unused_levels()

    # # interactive testing
    # mbias_stats_df = pd.read_pickle(fps["mbias_stats_p"])
    # flen_sel = config['plots']['mbias_flens_to_display']
    # # careful with slice, depends on whether mate is already present
    # mbias_stats_df = mbias_stats_df.loc[idxs[["CG"], :, :, :, flen_sel], :].copy()

    save_df_to_trunk_path(mbias_stats_df, fps["mbias_stats_trunk"])

    # Could be parallelized
    compute_derived_mbias_stats(mbias_stats_df, config)

    print("DONE")



def add_mate_info(mbias_stats_df):
    # Alternatively merge mate1 and mate2 calls
    # index_level_order = list(mbias_stats_df.index.names)
    # res = mbias_stats_df.reset_index("bs_strand")
    # res["bs_strand"] = res["bs_strand"].replace({"c_bc": "Read 1", "c_bc_rv": "Read 2",
    #                                              "w_bc": "Read 1", "w_bc_rv": "Read 2"})
    # res["dummy"] = 1
    # res = res.set_index(["bs_strand", "dummy"], append=True)
    # res = res.reorder_levels(index_level_order + ["dummy"])
    # res = res.groupby(level=index_level_order).sum().assign(beta_value = compute_beta_values)

    print("Adding mate info")
    bs_strand_to_read_mapping = {"c_bc": "Mate 1", "c_bc_rv": "Mate 2",
                                 "w_bc": "Mate 1", "w_bc_rv": "Mate 2"}

    # ~ 1.5 min for CG only
    mbias_stats_df["mate"] = (mbias_stats_df
                              .index
                              .get_level_values("bs_strand")
                              .to_series()
                              .replace(bs_strand_to_read_mapping)
                              .values
                              )

    # quick
    mbias_stats_df = (mbias_stats_df
                      .set_index("mate", append=True)
                      .reorder_levels(["motif", "seq_context",
                                       "mate", "bs_strand",
                                       "flen", "phred", "pos"])
                      )

    return(mbias_stats_df)


def compute_derived_mbias_stats(mbias_stats_df, config):
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

    print("Computing cutting sites")
    classic_mbias_stats_df = compute_classic_mbias_stats_df(mbias_stats_df)

    adjusted_cutting_sites = AdjustedCuttingSites(
        classic_mbias_stats_df, config)
    adjusted_cutting_sites_df = adjusted_cutting_sites.get_df()

    print("Computing masked M-bias stats DF")
    masked_mbias_stats_df = mask_mbias_stats_df(
        mbias_stats_df, adjusted_cutting_sites_df)

    print("Adding phred filtering info")
    phred_threshold_df = convert_phred_bins_to_thresholds(mbias_stats_df)

    phred_threshold_df_trimmed = convert_phred_bins_to_thresholds(
        masked_mbias_stats_df)

    print("Saving results")
    fps = config['paths']
    os.makedirs(fps['qc_stats_dir'], exist_ok=True, mode=0o770)

    # TODO: clean up
    fps["mbias_stats_classic_trunk"] = fps['mbias_stats_classic_p'].replace('.p', '')
    fps["mbias_stats_masked_trunk"] = fps['mbias_stats_masked_p'].replace('.p', '')
    fps["adjusted_cutting_sites_df_trunk"] = fps["adjusted_cutting_sites_df_p"].replace('.p', '')

    with open(fps['adjusted_cutting_sites_obj_p'], 'wb') as fobj:
        pickle.dump(adjusted_cutting_sites, fobj)
    save_df_to_trunk_path(adjusted_cutting_sites_df,
                          fps["adjusted_cutting_sites_df_trunk"])
    save_df_to_trunk_path(classic_mbias_stats_df,
                          fps["mbias_stats_classic_trunk"])
    save_df_to_trunk_path(masked_mbias_stats_df,
                          fps["mbias_stats_masked_trunk"])
    save_df_to_trunk_path(phred_threshold_df,
                          fps['mbias_stats_phred_threshold_trunk'])
    save_df_to_trunk_path(phred_threshold_df_trimmed,
                          fps['mbias_stats_masked_phred_threshold_trunk'])


def mbias_stat_plots(
        output_dir: str, sample_meta: dict,
        dataset_name_to_fp: dict,
        compact_mbias_plot_config_dict_fp: Optional[str] = None):
    """ M-bias plot generation manager

    Creates all M-bias plot specified through the config dictionary.
    The config file refers to shorthand dataset names. These
    names are mapped to filepaths in the dataset_name_to_fp dict.

    This function is accessible through the mbias_plots tool. Here, the
    config dictionary is defined via the path to a JSON config file. The
    datasets are given as comma-separated key=value pairs.

    Args:
        output_dir:
            All output filepaths are relative to the output_dir. It must
            be possible to make the output_dir the working directory.
        compact_mbias_plot_config_dict_fp:
            Path to the python module containing the M-bias plot config
            dict, plus the name of the desired config dict, append
            with '::'. See mqc/resources/default_mbias_plot_config.py
            for an example. The default config file also contains
            more documentation about the structure of the config file.
        dataset_name_to_fp:
            dataset name to path mapping

    Notes:
        Roadmap:
            - add sample metadata to output dfs? Currently the sample
              metadata are unused


    """
    # TODO: check that dataset filepaths are absolute? Or check that files exist right in the begining? document that relative to output_dir?
    # TODO: assert that all dataset files are there
    # TODO: check json config for 'None' strings - likely null was meant, raise
    # All filepaths are relative paths, need to go to correct output dir
    try:
        os.chdir(output_dir)
    except:
        raise OSError('Cannot enter the specified working directory')

    # Load compact mbias plot config dict from python module path
    # given as '/path/to/module::dict_name'
    if compact_mbias_plot_config_dict_fp is None:
        compact_mbias_plot_config_dict_fp = (
                get_resource_abspath('default_mbias_plot_config.py')
                + '::default_config')
    mbias_plot_config_module_path, config_dict_name = (
        compact_mbias_plot_config_dict_fp.split('::')
    )
    spec = importlib.util.spec_from_file_location(
        'cf', mbias_plot_config_module_path)
    cf = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(cf)
    compact_mbias_plot_config_dict = getattr(cf, config_dict_name)

    mbias_plot_configs = get_plot_configs(compact_mbias_plot_config_dict,
                                          dataset_name_to_fp=dataset_name_to_fp)

    aggregated_mbias_stats = AggregatedMbiasStats()
    create_aggregated_tables(mbias_plot_configs=mbias_plot_configs,
                             aggregated_mbias_stats=aggregated_mbias_stats,
                             dataset_filepath_mapping=dataset_name_to_fp)

    create_mbias_stats_plots(mbias_plot_configs=mbias_plot_configs,
                             aggregated_mbias_stats=aggregated_mbias_stats)
    # (Path(output_dir) / 'test.png').touch()


# def map_dataset_specs_to_filepaths(mbias_plot_config) -> None:
#     for curr_dict in mbias_plot_config.values():
#         # curr_dict is either plot param defaults or plot group dict
#         try:
#             curr_dict['datasets'] = [getattr(mqc.filepaths, name)
#                                      for name in curr_dict['datasets']]
#         except KeyError:
#             continue


class MbiasPlotAxisDefinition:
    def __init__(self, share: bool = True,
                 breaks: Union[int, List[float]] = 5,
                 limits: Optional[dict] = None,
                 rotate_labels: bool = True):
        self.breaks = breaks
        if limits is None:
            self.limits = {}
        self.limits = {'default': 'auto'}
        self.share = share
        self.rotate_labels = rotate_labels

    def __eq__(self, other):
        if isinstance(other, MbiasPlotAxisDefinition):
            return self.__dict__ == other.__dict__
        return False


class MbiasPlotParams:
    def __init__(self,
                 y_axis: MbiasPlotAxisDefinition,
                 x_axis: MbiasPlotAxisDefinition,
                 panel_height_cm: int = 6, panel_width_cm: int = 6,
                 theme: str = 'paper',
                 plot: Optional[List[str]] = None,
                 ):
        self.theme = theme
        self.panel_width_cm = panel_width_cm
        self.panel_height_cm = panel_height_cm
        self.y_axis = y_axis
        self.x_axis = x_axis
        if plot is None:
            self.plot = ['line']
        else:
            self.plot = plot

    def __eq__(self, other):
        if isinstance(other, type(self)):
            return self.__dict__ == other.__dict__
        return False


class MbiasPlotMapping:
    """Contains encoding and facetting information"""
    def __init__(self, x: str, y: str,
                 color: Optional[str] = None,
                 detail: Optional[str] = None,
                 column: Optional[str] = None,
                 row: Optional[str] = None,
                 # Implemented later
                 # column: Optional[Union[str, List, Tuple]] = None,
                 # row: Optional[Union[str, List, Tuple]] = None,
                 wrap: Optional[int] = None):
        """

        Facetting (row and col) may be done on combinations of variables. The
        facets are sorted based on the order of the facetting variables in the
        col and row values. Row and col may be specified as str, list or
        tuple, but will always be converted to tuple during init

        Args:
            x
            y
            color
            column
            row
            wrap: wrap determines the number of columns/rows if only row or
                column facetting is specified
        """
        self.x = x
        self.y = y
        self.color = color
        self.detail = detail

        # Multiple fields for facetting will be implemented later
        # ------------------------------
        # if column is None:
        #     self.column = column
        # elif isinstance(column, str):
        #     self.column = (column,)
        # elif isinstance(column, (tuple, list)):
        #     self.column = tuple(column)
        # else:
        #     raise TypeError('Unexpected type for col')
        #
        # if row is None:
        #     self.row = row
        # elif isinstance(row, str):
        #     self.row = (row, )
        # elif isinstance(row, (tuple, list)):
        #     self.row = tuple(row)
        # else:
        #     raise TypeError('Unexpected type for row')

        self.column = column
        self.row = row

        self.wrap = wrap

    # noinspection PyIncorrectDocstring
    def get_plot_aes_dict(self, include_facetting: bool =False,
                          x_name='x', y_name='y', color_name='color',
                          column_name='column', row_name='row',
                          detail_name='detail') -> dict:
        """Get mapping of column names to encoding channels (aes. mappings)

        Args:
            include_facetting: Include variables for row facetting
                and column facetting. Will not include wrap parameter
            *_name: name of the encoding channel in the plotting tool for which
                this mapping is produced

        Returns:
            dict with channel name -> column name mapping

        Notes:
            Vega encoding channels are equivalent to plotnine aesthetic
            mappings as well as all the other similar concepts in related
            tools. To adapt to different tools, the channel names can be set.

        """
        base_dict = {x_name: self.x,
                     y_name: self.y,
                     color_name: self.color,
                     detail_name: self.detail}

        if include_facetting:
            base_dict.update({row_name: self.row,
                              column_name: self.column})

        base_dict = {k: v for k, v in base_dict.items()
                     if v is not None}

        return base_dict

    def get_all_agg_variables_unordered(self) -> set:
        """Get dictionary with all variables to be used in agg. groupby

        Notes:
            - Returns set (which has no order) to facilitate comparison
            of variable lists, also in tests. Otherwise, tests would
            have to know the hardcoded order.
            - the 'dataset' variable is not considered, because
            aggregation currently happens per dataset
        """
        return {x for x in more_itertools.collapse(
            [self.x, self.column, self.row, self.color, self.detail])
                if x is not None and x != 'dataset'}

    # # TODO: adjust for tuples in col and row
    # def get_facetting_cmd(self):
    #     """Get facetting command for plotnine"""
    #     if self.wrap and self.col and self.row:
    #         raise ValueError("Can't set wrap, row and col all at once")
    #     elif self.wrap and not (self.col or self.row):
    #         raise ValueError("Wrap is set, but not a row or col variable")
    #     elif self.col and self.row:
    #         return gg.facet_grid(self.row, self.col)
    #     # self.wrap may be None
    #     elif self.col:
    #         return gg.facet_wrap(self.col, nrow=self.wrap)
    #     elif self.row:
    #         return gg.facet_wrap(self.row, ncol=self.wrap)

    def __eq__(self, other):
        if isinstance(other, MbiasPlotMapping):
            return self.__dict__ == other.__dict__
        return False


class MbiasPlotConfig:
    def __init__(self,
                 datasets: Dict[str, str],
                 aes_mapping: MbiasPlotMapping,
                 plot_params: MbiasPlotParams,
                 pre_agg_filters: Union[dict, None] = None,
                 post_agg_filters: Union[dict, None] = None):

        # Sanity checks
        if pre_agg_filters:
            self._check_filter_dicts(pre_agg_filters)
        if post_agg_filters:
            self._check_filter_dicts(post_agg_filters)
        assert isinstance(datasets, dict)

        self.datasets = datasets
        self.aes = aes_mapping
        self.plot_params = plot_params
        self.pre_agg_filters = pre_agg_filters
        self.post_agg_filters = post_agg_filters

    def __eq__(self, other):
        if isinstance(other, MbiasPlotConfig):
            return self.__dict__ == other.__dict__
        return False

    def __hash__(self):
        # TODO: test and improve
        return sum([hash_dict(x) for x in [
            tuple(self.datasets.values()),
            self.aes.__dict__,
            self.plot_params.__dict__,
            self.pre_agg_filters,
            self.post_agg_filters
        ]])

    def get_str_repr_for_filename_construction(self):
        return str(self.__hash__())


    def __str__(self):
        return str(self.__dict__)

    @staticmethod
    def _check_filter_dicts(filter_dict):
        """Make sure that filter dicts only contain lists of scalars

        Filtering is not meant to remove index levels. The filtering values
        are used for indexing. Because scalar values would remove their corresponding
        index levels, they are not allowed.
        """
        for filter_value in filter_dict.values():
            if not isinstance(filter_value, list):
                raise TypeError(f'Filter values must be list, but got this'
                                f' value of type {type(filter_value)}:'
                                f' {filter_value}')


def get_plot_configs(mbias_plot_config: dict,
                     dataset_name_to_fp: dict) -> List[MbiasPlotConfig]:
    """ Get MbiasPlotConfig objects for plotting function

    Args:
        mbias_plot_config: User-specified dictionary detailing the required
            M-bias plots in compact form
        dataset_name_to_fp: Mapping of the shorthand dataset names
            from the compact plot config to the actual filepaths


    Returns:
        List with one MbiasPlotConfig object per plotting task. Note that one
        task may act on several separate datasets at the same time.

    Note:
        Because MbiasPlotConfig objects may use several datasets at the
        same time, this list is not suited for generating the required
        aggregated variants of the used datasets. The function ... can
        extract aggregation task configs from the returned
        MbiasPlotConfig objects provided by this function
    """

    mbias_plot_config = deepcopy(mbias_plot_config)
    try:
        defaults = mbias_plot_config.pop('defaults')
    except KeyError:
        defaults = {}

    mbias_plot_configs = []
    for plot_group_name, plot_group_dict in mbias_plot_config.items():
        plot_group_dict = update_nested_dict(
            base_dict=defaults, custom_dict=plot_group_dict,
            allow_new_keys=True)

        try:
            for curr_datasets, curr_aes_mapping in product(
                    plot_group_dict['datasets'], plot_group_dict['aes_mappings']):
                assert isinstance(curr_datasets, (str, list)), \
                    f'Curr_dataset spec {curr_datasets} is not str, list'
                mbias_plot_configs.append(MbiasPlotConfig(
                    datasets=subset_dict(curr_datasets, dataset_name_to_fp),
                    aes_mapping=MbiasPlotMapping(
                        **curr_aes_mapping
                    ),
                    pre_agg_filters = plot_group_dict.get('pre_agg_filters'),
                    post_agg_filters = plot_group_dict.get('post_agg_filters'),
                    plot_params=MbiasPlotParams(**plot_group_dict['plot_params'])
                ))
        except (KeyError, TypeError) as e:
            raise ValueError(f'Incomplete configuration for {plot_group_name}. '
                             f'Original message: {e}')

    return mbias_plot_configs


@dataclass()
class MbiasStatAggregationParams:
    dataset_fp: str
    variables: Set[str]
    pre_agg_filters: dict

    def __hash__(self):
        return hash(json.dumps(self.pre_agg_filters, sort_keys=True)
                    + json.dumps(sorted(list(self.variables)))
                    + self.dataset_fp)

    # def __eq__(self, other):
    #     if (isinstance(other, type(self))
    #         and other.dataset_fp == self.dataset_fp
    #         and other.variables == self.variables
    #         and other.pre_agg_filters.equals(self.pre_agg_filters)):
    #         return True
    #     return False
    #

class AggregatedMbiasStats:
    """ Precompute and retrieve aggregated M-bias stats

    Aggregation calculations are only performed if the output file is
    not yet present or older than the input file. Currently, this class
    is written under the assumption that aggregated stats are
    precomputed in a first step, and may be retrieved several times in
    subsequent steps
    """

    def __init__(self):
        self.output_dir = mqc.filepaths.aggregated_mbias_stats_dir

    def _create_agg_dataset_fp(self, params: MbiasStatAggregationParams) -> Path:
        basename = (params.dataset_fp.replace('/', '-') + '_'
                    + '-'.join(params.variables) + '_'
                    + json.dumps(params.pre_agg_filters).replace(' ', '-')
                    + '.p')
        return self.output_dir / basename

    def get_aggregated_stats(self, params: MbiasStatAggregationParams) -> pd.DataFrame:
        """Retrieve precomputed M-bias stats"""
        fp = self._create_agg_dataset_fp(params)
        return pd.read_pickle(fp)

    def precompute_aggregated_stats(
            self, params: MbiasStatAggregationParams, mbias_stats_df) -> None:
        fp = self._create_agg_dataset_fp(params)
        if (not fp.exists()
                or fp.stat().st_mtime < Path(params.dataset_fp).stat().st_mtime):
            aggregated_df = aggregate_table(
                mbias_stats_df,
                unordered_variables=params.variables, pre_agg_filters=params.pre_agg_filters)

            fp.parent.mkdir(parents=True, exist_ok=True)

            aggregated_df.to_pickle(fp)


def create_aggregated_tables(
        mbias_plot_configs: List[MbiasPlotConfig],
        aggregated_mbias_stats: AggregatedMbiasStats,
        dataset_filepath_mapping: dict) -> None:
    """Compute aggregation tasks derived from MbiasPlotConfigs

    Aggregation tasks per dataset are determined from the
    MbiasPlotConfigs. Note that one plot config may contain several
    datasets, but aggregation results are computed and cached per
    (single) dataset. Aggregation tasks are computed by passing the task
    definitions to the AggregatedMbiasStats.precompute_aggregated_stats
    method. This method will take care of storing the aggregated data
    and will recognize if the data are already cached in a sufficienlty
    recent version.

    Args:
        mbias_plot_configs
        aggregated_mbias_stats
        dataset_filepath_mapping: Datasets are referenced by a shorthand
            name in the MbiasPlotConfigs. This mapping allows finding the
            corresponding dataset. Lookup for already existing results is
            always done based on the absolute dataset filepath, not based on
            the shorthand name

    Note:
        Aggregation tasks are performed explicitely (as opposed to
        implicitely within the plotting tasks) because this saves
        significant time during plot optimization. Perhaps even more
        importantly this allows quick combination of aggregated
        dataframes across many samples for cohort-wide plots
    """
    # TODO: make sure that absolute dataset path is used for hashing/storage
    # or change docstring

    # Collect individual (per dataset) aggregation tasks
    # M-bias plot configs may include more than one dataset, but aggregation
    # is performed per dataset, so we need to further split the MbiasPlotConfigs
    single_dataset_plot_configs = []
    for curr_mbias_plot_config in mbias_plot_configs:
        for curr_dataset_fp in curr_mbias_plot_config.datasets.values():
            single_dataset_plot_configs.append(
                MbiasStatAggregationParams(
                    dataset_fp=curr_dataset_fp,
                    variables=curr_mbias_plot_config.aes.get_all_agg_variables_unordered(),
                    pre_agg_filters=curr_mbias_plot_config.pre_agg_filters
                ))

    # avoid duplicates to avoid interfering writes when parallel processing
    unique_single_dataset_plot_configs = set(single_dataset_plot_configs)

    # Group by dataset - allows either loading datasets sequentially to
    # save memory, or to start one process per dataset
    plots_by_dataset = tz.groupby(lambda x: x.dataset_fp,
                                  unique_single_dataset_plot_configs)
    for dataset_fp, agg_configs in plots_by_dataset.items():
        agg_configs: List[MbiasStatAggregationParams]
        curr_mbias_stats_df = pd.read_pickle(dataset_fp)
        for curr_agg_config in agg_configs:
            aggregated_mbias_stats.precompute_aggregated_stats(
                params=curr_agg_config,
                mbias_stats_df=curr_mbias_stats_df)


def aggregate_table(mbias_stats_df: pd.DataFrame,
                    unordered_variables: Set[str],
                    pre_agg_filters: Optional[dict] = None) -> pd.DataFrame:
    """ Prepare aggregated M-bias stats variant for the corresponding plot

    Will first apply filtering as per the pre_agg_filter dict. This dict
    may only contain list-based specification of level values which are
    to be retained for a given index level. Therefore, index levels are
    not removed by filtering.

    In a second step, groupby the plot variables (maintaining the
    original order in the index) and do a sum aggregation. Beta values
    are recomputed.

    Args:
        mbias_stats_df: Arbitrary index levels. Must contain columns
            n_meth and n_unmeth. May contain beta_value column.
        unordered_variables: All variables which will be used in the
            corresponding plot
        pre_agg_filters: dict mapping index levels to lists of values
            which should exclusively be retained

    Returns: The filtered and aggregated DataFrame.
    """

    # perform sanity checks, this is the first point where we know which index
    # levels we are actually dealing with
    assert_match_between_variables_and_index_levels(
        variables=unordered_variables,
        index_level_names=mbias_stats_df.index.names
    )
    ordered_groupby_vars = [x for x in mbias_stats_df.index.names
                            if x in unordered_variables]

    if pre_agg_filters:  # may be None
        assert_match_between_variables_and_index_levels(
            variables=list(pre_agg_filters.keys()),
            index_level_names=mbias_stats_df.index.names
        )

        mbias_stats_df = mbias_stats_df.loc[nidxs(**pre_agg_filters), :]

    agg_df = (mbias_stats_df
              .groupby(ordered_groupby_vars)
              .sum()
              .assign(beta_value=compute_beta_values)
              )
    return agg_df


def calculate_plot_df(mbias_stats_df, complete_aes, flens_to_display,
                      plot_df_var_names):

    # We need to aggregate variables which are not part of the plot
    # ourselves, seaborn can't do that
    # For this purpose, find all variables we use for aesthetics,
    # then use these to do a groupby.sum
    groupby_vars = (
        [val for key, val in complete_aes.items()
         if (val is not None) and key != "y"])

    # We can't show all flens. Selection of flens is config param
    if 'flen' in groupby_vars:
        # TODO: make robust against changes in index levels
        idx = idxs[:, :, :, :, flens_to_display]
        mbias_stats_df = mbias_stats_df.loc[idx, :]

    # Aggregate over variables which are not aesthetics in the plot
    # Recompute beta values
    plot_df = (mbias_stats_df
               .groupby(groupby_vars)
               .sum()
               .assign(beta_value=compute_beta_values)
               .reset_index()
               )

    # Remove unused levels in index and categoricals
    # TODO: make more general
    if "seq_context" in plot_df:
        plot_df["seq_context"].cat.remove_unused_categories(inplace=True)
    if "phred" in plot_df:
        plot_df["phred"] = pd.Categorical(plot_df["phred"])
    if "flen" in plot_df:
        plot_df["flen"] = pd.Categorical(plot_df["flen"])

    # drop nas so that lines are drawn "across" NAs
    # https://stackoverflow.com/questions/9617629/connecting-across-missing-values-with-geom-line
    plot_df = plot_df.loc[~plot_df[complete_aes["y"]].isnull(), :]

    # plot_df = plot_df.rename(columns=plot_df_var_names)

    # need to discard index because it may have "holes"
    # then I can't serialize to feather, and not use R ggplot
    plot_df = plot_df.reset_index(drop=True)

    return plot_df


# TODO: remove duplicate function
def create_mbias_stats_plots(mbias_plot_configs: List[MbiasPlotConfig],
                             aggregated_mbias_stats: AggregatedMbiasStats) -> None:
    for plot_config in mbias_plot_configs:
        per_dataset_dfs = []
        for curr_dataset_fp in plot_config.datasets.values():
            per_dataset_dfs.append(aggregated_mbias_stats.get_aggregated_stats(
                params = MbiasStatAggregationParams(
                    dataset_fp=curr_dataset_fp,
                    variables=plot_config.aes.get_all_agg_variables_unordered(),
                    pre_agg_filters=plot_config.pre_agg_filters
                )))
        full_plot_df = pd.concat(
            per_dataset_dfs, axis=0,
            keys=plot_config.datasets.keys(),
            names=['dataset'] + per_dataset_dfs[0].index.names)
        chart = create_single_mbias_stat_plot(plot_config, full_plot_df)
        target_fp = (mqc.filepaths.mbias_plots_trunk.with_name(
            mqc.filepaths.mbias_plots_trunk.name
            + plot_config.get_str_repr_for_filename_construction()
            + '.html'))
        print(target_fp)
        target_fp.parent.mkdir(parents=True, exist_ok=True)
        chart.save(str(target_fp))


def create_single_mbias_stat_plot(plot_config: MbiasPlotConfig,
                                  df: pd.DataFrame) -> alt.Chart:

    chart = (alt.Chart(df.reset_index())
             .mark_line()
             .encode(**plot_config.aes.get_plot_aes_dict(include_facetting=True))
             ).interactive()
    return chart



# def mbias_stat_plots(config):
#     """Create M-bias plots
#
#     Creates all standard M-bias plots as documented in the bseqtools QC section.
#
#     This function can be called from the command line using bseqtools mbias_plots
#     """
#
#     fps = config['paths']
#
#     mbias_stats_dfs_dict = {}
#
#     print("Reading full M-bias stats DF")
#     mbias_stats_dfs_dict["full"] = pd.read_pickle(fps["mbias_stats_trunk"] + ".p")
#     print("Reading trimmed M-bias stats DF")
#     mbias_stats_dfs_dict["trimmed"] = pd.read_pickle(fps["mbias_stats_masked_trunk"] + ".p")
#     print("Reading phred threshold stats DF")
#     mbias_stats_dfs_dict['full_flen-agg_cg_phred-threshold'] = pd.read_pickle(fps['mbias_stats_phred_threshold_trunk'] + '.p')
#     mbias_stats_dfs_dict['trimmed_flen-agg_cg_phred-threshold'] = pd.read_pickle(fps['mbias_stats_masked_phred_threshold_trunk'] + '.p')
#
#     print('Creating plots')
#     create_mbias_stats_plots(mbias_stats_dfs_dict, config)


def convert_phred_bins_to_thresholds(mbias_stats_df):
    """Takes full Mbias stats df and returns CG only, flen agg df
    with phred bins converted to phred thresholds

    The operation would take prohibitively long on the entire dataframe
    (> 500 min ?) in the current implementation

    Takes approx. 2 min

    Important
    ---------
    This function assumes that mbias_stats_df is sorted
    """


    # *This function assumes that mbias_stats_df is sorted*

    non_flen_levels = [n for n in mbias_stats_df.index.names if n != "flen"]
    non_flen_phred_levels = [n for n in non_flen_levels if n != "phred"]

    print("Aggregating counts with same flen")
    mbias_stats_df_cg = mbias_stats_df.loc[idxs['CG':'CG'], :]
    # ~ 15s
    mbias_stats_df_cg_flen_agg = (mbias_stats_df_cg
                                  .drop(['beta_value'], axis=1)
                                  .groupby(non_flen_levels)
                                  .sum())


    print("Computing phred threshold scores")
    # ~ 1.5 min
    def compute_phred_threshold_counts_for_group(ser):
        """Will work on n_meth and n_unmeth"""
        cum_ser = ser.cumsum()
        return cum_ser[-1] - cum_ser
    res = (mbias_stats_df_cg_flen_agg
           .groupby(non_flen_phred_levels)
           .transform(compute_phred_threshold_counts_for_group)
           .assign(beta_value=compute_beta_values)
           )

    # Discard the highest phred bin - it has no events left after filtering
    phred_idx = res.index.get_level_values("phred").unique()[:-2]
    res = res.query("phred in @phred_idx")

    return res


def compute_beta_values(df):
    return df['n_meth'] / (df['n_meth'] + df['n_unmeth'])


def save_df_to_trunk_path(df, trunk_path):
    print(f"Saving {op.basename(trunk_path)}")
    print("Saving pickle")
    df.to_pickle(trunk_path + '.p')
    print("Saving feather")
    df.reset_index().to_feather(trunk_path + '.feather')
