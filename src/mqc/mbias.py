import importlib
import itertools
from itertools import product, count
import json
import os
import os.path as op
from copy import deepcopy
from dataclasses import dataclass
# from dpcontracts import invariant
from math import floor, ceil
import more_itertools
from pathlib import Path
import toolz as tz
from typing import Dict, List, Union, Optional, Set, Any, Tuple, Mapping

import matplotlib

matplotlib.use('Agg')  # import before pyplot import!
import altair as alt

import numpy as np
import pandas as pd
import scipy.signal
import scipy.stats
from sklearn import linear_model
from figure_report import Report
idxs = pd.IndexSlice

from mqc.flag_and_index_values import bsseq_strand_indices as bstrand_idxs
# noinspection PyUnresolvedReferences
from mqc.pileup.bsseq_pileup_read import BSSeqPileupRead
from mqc.pileup.pileup import MotifPileup
from mqc.utils import hash_dict, get_resource_abspath
from mqc.visitors import Counter
import mqc.filepaths
from mqc.utils import (update_nested_dict, NamedIndexSlice,
                       assert_match_between_variables_and_index_levels,
                       subset_dict)
nidxs = NamedIndexSlice

from mqc.flag_and_index_values import (
    bsseq_strand_na_index as b_na_ind,
    methylation_status_flags as m_flags
)


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

ConfigDict = Dict[str, Any]


class MbiasCounter(Counter):
    """Counter for multidimensional M-bias stats

    Implementation notes:

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
        - the last fragment length bin is guaranteed to end with the max_flen
        - the last bin may therefore be smaller than the specified bin size
        - phred scores are always binned, and are also labelled with the rightmost phred score in the bin
        - pos labels: 1-based (note that array indexing is 0-based)
        - meth. status labels: n_meth, n_unmeth [categorical]
     """

    def __init__(self, config: ConfigDict) -> None:

        self.save_stem = config["paths"]["mbias_counts"]
        self.max_read_length = config["run"]["read_length"]
        sts = config["stats"]
        self.max_flen = sts["max_flen"]
        self.max_single_flen = sts["max_flen_with_single_flen_resolution"]
        self.flen_bin_size = sts["flen_bin_size"]
        self.max_phred = sts["max_phred"]
        self.phred_bin_size = sts["phred_bin_size"]
        # TODO-api_change: get from IndexFile
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

        def get_categorical(ordered_str_labels: List[str]):
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

    def process(self, motif_pileup: MotifPileup) -> None:
        """Extract M-bias stats from MotifPileup

        The entire MotifPileup is skipped if it has an undefined
        sequence context (e.g. containing N)

        Reads are discarded if

        - they have a qc_fail_flag
        - their bsseq strand could not be determined
        - they have methylation calling status NA, SNP or Ref

        Reads are not discarded if they are phred score fails. The phred
        score information is recorded as part of the M-bias stats. In
        the standard MbiasDeterminationRun, the PhredFilter is therefore
        not applied at all. You could however apply it if you use this
        Counter as part of a larger workflow, because the
        phred_fail_flag is not checked here.

        M-bias stat dimensions recorded:
        - motif
        - sequence context
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

        # noinspection PyUnusedLocal
        curr_read: BSSeqPileupRead
        for curr_read in motif_pileup.reads:
            # note: this does explicitely not discard phred_score_fail
            # we stratify the M-bias stats by phred score and do
            # pseudo-filtering as appropriate for the different M-bias plots
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


def get_sequence_context_to_array_index_table(motif_size: int) \
        -> Tuple[Dict[str, int], Dict[str, int]]:
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


def map_seq_ctx_to_motif(seq_ctx: str, use_classical: bool = True) -> str:
    """Map sequence context strings containing [ACGTW] to motifs

    Motifs may be classical: [CG, CHG, CHH] or extended (composed of C,G,W)
    # TODO-minor: Ns? at the end of chroms?
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


def fit_normalvariate_plateau(group_df: pd.DataFrame, config: ConfigDict) -> pd.Series:
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

    # Pycharm false positive due to for...else flow
    # noinspection PyUnboundLocalVariable
    return pd.Series([best_start, best_end],
                     index=['left_cut_end', 'right_cut_end'])


# noinspection PyUnreachableCode
def fit_percentiles(group_df: pd.DataFrame) -> pd.Series:

    raise NotImplementedError
    # remove hardcoding before using this
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

    # Pycharm false positive due to for...else flow
    # noinspection PyUnboundLocalVariable
    return pd.Series([best_start, best_end],
                     index=['left_cut_end', 'right_cut_end'])


def compute_mbias_stats_df(mbias_counter_fp_str: str) -> pd.DataFrame:
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


def prepend_motif_level_to_index(mbias_stats_df: pd.DataFrame) -> pd.DataFrame:
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


def compute_classic_mbias_stats_df(mbias_stats_df: pd.DataFrame) -> pd.DataFrame:
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


def mask_mbias_stats_df(mbias_stats_df: pd.DataFrame, cutting_sites_df: pd.DataFrame) -> pd.DataFrame:
    """Cover pos in trimming zones (depends on bs_strand, flen...) with NA"""

    print("Masking dataframe")

    def mask_trimming_zones(group_df, cutting_sites_df):
        bs_strand, flen, pos = group_df.name
        left, right = cutting_sites_df.loc[(bs_strand, flen), ['start', 'end']]
        # left, right are indicated as slice(left, right) for 0-based
        # position coordinates
        left += 1
        return left <= pos <= right

    return (mbias_stats_df
            .groupby(['bs_strand', 'flen', 'pos'])
            .filter(mask_trimming_zones, dropna=False,
                    cutting_sites_df=cutting_sites_df)
            )




def compute_mbias_stats(config: ConfigDict) -> None:
    """ Run standard analysis on M-bias stats"""

    fps = config['paths']

    os.makedirs(fps['qc_stats_dir'], exist_ok=True, mode=0o770)

    if (Path(fps['mbias_stats_trunk'] + '.p').exists()
        and config['run']['use_cached_mbias_stats']):
        print('Reading mbias stats from previously computed pickle')
        mbias_stats_df = pd.read_pickle(fps['mbias_stats_trunk'] + '.p')
        print(mbias_stats_df.head())
    else:
        print("Computing M-bias stats from M-bias counter array")
        mbias_stats_df = compute_mbias_stats_df(fps['mbias_counts'] + '.p')
        mbias_stats_df = add_mate_info(mbias_stats_df)
        mbias_stats_df = mbias_stats_df.sort_index()

        print('Discarding unused phred scores')
        n_total = mbias_stats_df["n_meth"] + mbias_stats_df["n_unmeth"]
        phred_group_sizes = n_total.groupby("phred").sum()
        phred_bin_has_counts = (phred_group_sizes > 0)
        existing_phred_scores = phred_group_sizes.index.values[phred_bin_has_counts]
        mbias_stats_df = mbias_stats_df.loc[idxs[:, :, :, :, :, existing_phred_scores], :]
        mbias_stats_df.index = mbias_stats_df.index.remove_unused_levels()

        save_df_to_trunk_path(mbias_stats_df, fps["mbias_stats_trunk"])

    print('Computing derived stats')
    compute_derived_mbias_stats(mbias_stats_df, config)

    print("DONE")



def add_mate_info(mbias_stats_df: pd.DataFrame) -> pd.DataFrame:
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

    return mbias_stats_df


def compute_derived_mbias_stats(mbias_stats_df: pd.DataFrame, config: ConfigDict) -> None:
    """Compute various objects derived from the M-bias stats dataframe

    All objects are stored on disc, downstream evaluation is done in
    separate steps, e.g. by using the mbias_plots command.

    Notes:
        - the plateau detection is performed only based on CG context
          calls, also if information for CHG and CHH context is present.
          The cutting sites determined in this way will then be applied
          across all sequence contexts
    """

    plateau_detection_params = config['run']['plateau_detection']
    # Will be used when more algorithms are implemented
    unused_plateau_detection_algorithm = plateau_detection_params.pop('algorithm')

    print("Computing cutting sites")
    classic_mbias_stats_df = compute_classic_mbias_stats_df(mbias_stats_df)
    classic_mbias_stats_df_with_n_total = classic_mbias_stats_df.assign(
        n_total = lambda df: df['n_meth'] + df['n_unmeth']
    )
    cutting_sites = CuttingSites.from_mbias_stats(
            mbias_stats_df=classic_mbias_stats_df_with_n_total.loc['CG', :],
            **plateau_detection_params)

    print("Computing masked M-bias stats DF")
    masked_mbias_stats_df = mask_mbias_stats_df(
        mbias_stats_df, cutting_sites.df)

    print("Adding phred filtering info")
    phred_threshold_df = convert_phred_bins_to_thresholds(mbias_stats_df)

    phred_threshold_df_trimmed = convert_phred_bins_to_thresholds(
        masked_mbias_stats_df)

    print("Saving results")
    fps = config['paths']
    os.makedirs(fps['qc_stats_dir'], exist_ok=True, mode=0o770)

    # TODO-refactor: clean up
    fps["mbias_stats_classic_trunk"] = fps['mbias_stats_classic_p'].replace('.p', '')
    fps["mbias_stats_masked_trunk"] = fps['mbias_stats_masked_p'].replace('.p', '')
    fps["adjusted_cutting_sites_df_trunk"] = fps["adjusted_cutting_sites_df_p"].replace('.p', '')

    save_df_to_trunk_path(cutting_sites.df,
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
        output_dir: str,
        dataset_name_to_fp: dict,
        compact_mbias_plot_config_dict_fp: Optional[str] = None) -> None:
    """ Analysis workflow creating the M-bias plots and report

    Creates all M-bias plots specified through the compact
    dict representation of the MbiasAnalysisConfig. The config
    file refers to shorthand dataset names. These names
    are mapped to filepaths in the dataset_name_to_fp dict.

    This function is accessible through the cli mbias_plots tool.
    The datasets are given as comma-separated key=value pairs.

    Args:
        output_dir:
            All output filepaths are relative to the output_dir. It must
            be possible to make the output_dir the working directory.
        compact_mbias_plot_config_dict_fp:
            Path to the python module containing the M-bias plot config
            dict, plus the name of the desired config dict, appended
            with '::'. See mqc/resources/default_mbias_plot_config.py
            for an example. The default config file also contains
            more documentation about the structure of the config file.
        dataset_name_to_fp:
            dataset name to path mapping

    Notes:
        Roadmap:
            - this function will be refactored into a general PlotAnalysis class
            - add sample metadata to output dfs? Currently the sample
              metadata are unused
    """

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
        compact_mbias_plot_config_dict_fp.split('::')  # type: ignore
    )

    # Pycharm cant deal with importlib.util attribute
    # noinspection PyUnresolvedReferences
    spec = importlib.util.spec_from_file_location(
        'cf', mbias_plot_config_module_path)
    # Pycharm cant deal with importlib.util attribute
    # noinspection PyUnresolvedReferences
    cf = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(cf)  # type: ignore
    compact_mbias_plot_config_dict = getattr(cf, config_dict_name)

    mbias_plot_configs = get_plot_configs(compact_mbias_plot_config_dict,
                                          dataset_name_to_fp=dataset_name_to_fp)

    aggregated_mbias_stats = AggregatedMbiasStats()
    create_aggregated_tables(mbias_plot_configs=mbias_plot_configs,
                             aggregated_mbias_stats=aggregated_mbias_stats)

    create_mbias_stats_plots(mbias_plot_configs=mbias_plot_configs,
                             aggregated_mbias_stats=aggregated_mbias_stats)

    mbias_report_config = MbiasAnalysisConfig.from_compact_config_dict(
        compact_mbias_plot_config_dict,
        dataset_name_to_fp
    ).get_report_config()

    Report({'M-bias': mbias_report_config}).generate(
        mqc.filepaths.qc_report_dir
    )

    # (Path(output_dir) / 'test.png').touch()




class MbiasPlotAxisDefinition:
    def __init__(self, share: bool = True,
                 breaks: Union[int, List[float]] = 5,
                 limits: Optional[dict] = None,
                 rotate_labels: bool = True) -> None:
        self.breaks = breaks
        if limits is None:
            self.limits: Dict = {}
        self.limits = {'default': 'auto'}
        self.share = share
        self.rotate_labels = rotate_labels

    def __eq__(self, other: Any) -> bool:
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
                 ) -> None:
        self.theme = theme
        self.panel_width_cm = panel_width_cm
        self.panel_height_cm = panel_height_cm
        self.y_axis = y_axis
        self.x_axis = x_axis
        if plot is None:
            self.plot = ['line']
        else:
            self.plot = plot

    def __eq__(self, other: Any) -> bool:
        if isinstance(other, type(self)):
            return self.__dict__ == other.__dict__
        return False

    def __str__(self):
        return str(self.__dict__)

    def __repr__(self):
        return str(self.__dict__)


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
                 wrap: Optional[int] = None) -> None:
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
    def get_plot_aes_dict(self, include_facetting: bool = False,
                          x_name: str = 'x', y_name: str = 'y', color_name: str = 'color',
                          column_name: str = 'column', row_name: str = 'row',
                          detail_name: str = 'detail') -> Dict:
        """Get mapping of column names to encoding channels (aes. mappings)

        Args:
            include_facetting: Include variables for row facetting
                and column facetting. Will not include wrap parameter

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

    def get_facetting_vars(self, column_name: str = 'column', row_name:str = 'row'):
        return {k: v for k, v in {column_name: self.column,
                                  row_name: self.row}.items()
                if v is not None}

    def get_all_agg_variables_unordered(self) -> Set[str]:
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
                if x is not None and x not in ['dataset', 'statistic']}

    def __eq__(self, other: Any) -> bool:
        if isinstance(other, MbiasPlotMapping):
            return self.__dict__ == other.__dict__
        return False

    def __str__(self):
        return str(self.__dict__)

    def __repr__(self):
        return str(self.__dict__)

class MbiasPlotConfig:
    def __init__(self,
                 datasets: Dict[str, str],
                 aes_mapping: MbiasPlotMapping,
                 plot_params: MbiasPlotParams,
                 pre_agg_filters: Union[dict, None] = None,
                 post_agg_filters: Union[dict, None] = None) -> None:

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

    def __eq__(self, other: Any) -> bool:
        if isinstance(other, MbiasPlotConfig):
            return self.__dict__ == other.__dict__
        return False

    def __hash__(self):
        return sum([hash_dict(x) for x in [
            tuple(self.datasets.values()),
            self.aes.__dict__,
            self.plot_params.__dict__,
            self.pre_agg_filters,
            self.post_agg_filters
        ]])

    def get_str_repr_for_filename_construction(self) -> str:
        return str(self.__hash__())


    def __str__(self):
        return str(self.__dict__)

    @staticmethod
    def _check_filter_dicts(filter_dict: Dict) -> None:
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


def get_plot_configs(mbias_plot_config: Dict,
                     dataset_name_to_fp: Dict[str, str]) -> List[MbiasPlotConfig]:
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
    pre_agg_filters: Optional[Dict[str, List]]

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

    def __init__(self) -> None:
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
            self, params: MbiasStatAggregationParams, mbias_stats_df: pd.DataFrame) -> None:
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
        aggregated_mbias_stats: AggregatedMbiasStats) -> None:
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

    Note:
        Aggregation tasks are performed explicitely (as opposed to
        implicitely within the plotting tasks) because this saves
        significant time during plot optimization. Perhaps even more
        importantly this allows quick combination of aggregated
        dataframes across many samples for cohort-wide plots
    """
    # TODO-minor: make sure that absolute dataset path is used for hashing/storage or change docstring

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
    # noinspection PyUnusedLocal
    agg_configs: List[MbiasStatAggregationParams]
    for dataset_fp, agg_configs in plots_by_dataset.items():
        curr_mbias_stats_df = pd.read_pickle(dataset_fp)
        for curr_agg_config in agg_configs:
            aggregated_mbias_stats.precompute_aggregated_stats(
                params=curr_agg_config,
                mbias_stats_df=curr_mbias_stats_df)


def aggregate_table(
        mbias_stats_df: pd.DataFrame,
        unordered_variables: Set[str],
        pre_agg_filters: Optional[Dict[str, List]] = None) -> pd.DataFrame:
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
        pre_agg_filters: dict mapping index levels to *lists* of values
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


def calculate_plot_df(mbias_stats_df: pd.DataFrame, complete_aes: Dict, flens_to_display: List[int]) -> pd.DataFrame:

    # We need to aggregate variables which are not part of the plot
    # ourselves, seaborn can't do that
    # For this purpose, find all variables we use for aesthetics,
    # then use these to do a groupby.sum
    groupby_vars = (
        [val for key, val in complete_aes.items()
         if (val is not None) and key != "y"])

    # We can't show all flens. Selection of flens is config param
    if 'flen' in groupby_vars:
        # TODO-minor: make robust against changes in index levels
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
    # TODO-minor: make more general
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

class MbiasAnalysisConfig:
    """Configuration of the analysis workflow creating M-bias plots"""
    def __init__(self, grouped_mbias_plot_configs: Dict,
                 dataset_name_to_fp: Mapping[str, str]) -> None:
        """MbiasAnalysiConfig default constructor

        In most cases, one of the convenience
        constructors is the better choice.


        Args:
            grouped_mbias_plot_configs:
                Mapping of section names to List[MbiasPlotConfig].
                The section names are used during plot generation.
            dataset_name_to_fp:
                Mapping of dataset shorthand names to dataset filepaths

        Notes:
            - Roadmap
                - allow multi-level (nested) sections. This will require switching
                  to recursive processing function in the clients of this class
        """
        self.grouped_mbias_plot_configs = grouped_mbias_plot_configs
        self.dataset_name_to_fp = dataset_name_to_fp

    @staticmethod
    def from_compact_config_dict(mbias_plot_config: Dict, dataset_name_to_fp: Mapping[str, str]):
        """Constructor based on compact config dict for the analysis

        Args:
            mbias_plot_config: User-specified dictionary detailing the required
                M-bias plots in compact form. See default_mbias_plot_config.py
                for example
            dataset_name_to_fp: Mapping of the shorthand dataset names
                from the compact plot config to the actual filepaths

        Returns:
            MbiasAnalysisConfig instance
        """


        mbias_plot_config = deepcopy(mbias_plot_config)
        try:
            defaults = mbias_plot_config.pop('defaults')
        except KeyError:
            defaults = {}

        grouped_mbias_plot_configs: Dict[str, Any] = {}
        for plot_group_name, plot_group_dict in mbias_plot_config.items():

            grouped_mbias_plot_configs[plot_group_name] = []
            plot_group_dict = update_nested_dict(
                base_dict=defaults, custom_dict=plot_group_dict,
                allow_new_keys=True)

            try:
                for curr_datasets, curr_aes_mapping in product(
                        plot_group_dict['datasets'], plot_group_dict['aes_mappings']):
                    assert isinstance(curr_datasets, (str, list)), \
                        f'Curr_dataset spec {curr_datasets} is not str, list'
                    grouped_mbias_plot_configs[plot_group_name].append(MbiasPlotConfig(
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

        return MbiasAnalysisConfig(grouped_mbias_plot_configs=grouped_mbias_plot_configs,
                                   dataset_name_to_fp=dataset_name_to_fp)

    def get_report_config(self) -> Dict[str, Dict]:
        """Return config dict in the format used by figure_report.Report

        Returns:
            A mapping of sections to their contained figures, with
            titles. The M-bias plot config objects are turned
            into figure definition dicts, ie are represented by
            the filepath where the plot can be found, as well as
            a title derived from the MbiasPlotConfig information.
        """

        report_dict: Dict = {}

        def get_figure_config(mbias_plot_config: MbiasPlotConfig):
            title = 'Datasets: ' + ', '.join(mbias_plot_config.datasets) + '\n\n'
            title += json.dumps(mbias_plot_config
                                    .aes
                                    .get_plot_aes_dict(include_facetting=True),
                                sort_keys=True) + '\n'
            # When the MbiasPlotStrip class is created, this can be replaced
            # by calling its filepath query method
            path = (mqc.filepaths.mbias_plots_trunk.name
                    + mbias_plot_config.get_str_repr_for_filename_construction()
                    + '.json')
            figure_config = {'path': path, 'title': title}
            return figure_config

        for group_name, mbias_plot_configs_list in self.grouped_mbias_plot_configs.items():
            report_dict[group_name] = {'figures': []}
            for mbias_plot_config in mbias_plot_configs_list:
                report_dict[group_name]['figures'].append(
                    get_figure_config(mbias_plot_config)
                )

        return report_dict


def create_mbias_stats_plots(mbias_plot_configs: List[MbiasPlotConfig],
                             aggregated_mbias_stats: AggregatedMbiasStats) -> None:
    for curr_plot_config in mbias_plot_configs:
        per_dataset_dfs = []
        for curr_dataset_fp in curr_plot_config.datasets.values():
            curr_df = aggregated_mbias_stats.get_aggregated_stats(
                params = MbiasStatAggregationParams(
                    dataset_fp=curr_dataset_fp,
                    variables=curr_plot_config.aes.get_all_agg_variables_unordered(),
                    pre_agg_filters=curr_plot_config.pre_agg_filters
                ))
            if curr_plot_config.post_agg_filters:
                relevant_post_agg_filters = {
                    k: v for k, v in curr_plot_config.post_agg_filters.items()
                    if k in curr_df.index.names}
                if relevant_post_agg_filters:
                    # Raise if elements which should pass the filter
                    # are missing (This will not be necessary in
                    # the future, as this error will be raised
                    # automatically in one of the next pandas versions
                    for index_col_name, elem_list in relevant_post_agg_filters.items():
                        unknown_elems = (
                                set(elem_list)
                                - set(curr_df.index.get_level_values(index_col_name).unique()))
                        if unknown_elems:
                            raise ValueError(
                                'Post-agg filter contains unknown elements.\n'
                                f'Filter: {relevant_post_agg_filters}'
                                f'Data head: {curr_df.head()}')
                    curr_df = curr_df.loc[nidxs(**relevant_post_agg_filters), :]
            per_dataset_dfs.append(curr_df)
        full_plot_df = pd.concat(
            per_dataset_dfs, axis=0,
            keys=curr_plot_config.datasets.keys(),
            names=['dataset'] + per_dataset_dfs[0].index.names)
        chart = create_single_mbias_stat_plot(curr_plot_config, full_plot_df)
        target_fp = (mqc.filepaths.mbias_plots_trunk.with_name(
            mqc.filepaths.mbias_plots_trunk.name
            + curr_plot_config.get_str_repr_for_filename_construction()
            + '.json'))
        print(target_fp)
        target_fp.parent.mkdir(parents=True, exist_ok=True)
        try:
            chart.save(str(target_fp))
        except ValueError as e:
            print('Error working on: ', curr_plot_config)
            raise e



def create_single_mbias_stat_plot(plot_config: MbiasPlotConfig,
                                  df: pd.DataFrame) -> alt.Chart:

    assert set(df.columns.values) == {'n_meth', 'n_unmeth', 'beta_value'}
    if plot_config.aes.y == 'value':
        df.columns.name = 'statistic'
        df = df.stack().to_frame('value')

    x_zoom =alt.selection_interval(bind='scales', encodings=['x'],
                               zoom="wheel![event.shiftKey]")
    y_zoom =alt.selection_interval(bind='scales', encodings=['y'])

    chart = (alt.Chart(df.reset_index())
             .mark_line()
             .encode(**plot_config.aes.get_plot_aes_dict(include_facetting=False))
             .properties(selection=x_zoom + y_zoom)
             )

    facet_dict = plot_config.aes.get_facetting_vars()
    if facet_dict:
        chart = (chart
                 .facet(**facet_dict)
                 .resolve_scale(x='shared', y='independent')
                 )
    return chart



def convert_phred_bins_to_thresholds(mbias_stats_df: pd.DataFrame) -> pd.DataFrame:
    """Takes full Mbias stats df and returns CG only, flen agg df
    with phred bins converted to phred thresholds

    The operation would take prohibitively long on the entire dataframe
    (> 500 min ?) in the current implementation

    Takes approx. 2 min

    **Important:** This function assumes that mbias_stats_df is sorted
    """


    # *This function assumes that mbias_stats_df is sorted*

    non_flen_levels = [n for n in mbias_stats_df.index.names if n != "flen"]
    non_flen_phred_levels = [n for n in non_flen_levels if n != "phred"]

    print("Aggregating counts with same flen")
    mbias_stats_df_cg = mbias_stats_df.loc[idxs['CG':'CG'], :] # type: ignore
    # ~ 15s
    mbias_stats_df_cg_flen_agg = (mbias_stats_df_cg
                                  .drop(['beta_value'], axis=1)
                                  .groupby(non_flen_levels)
                                  .sum())


    print("Computing phred threshold scores")
    # ~ 1.5 min
    def compute_phred_threshold_counts_for_group(ser: pd.Series) -> pd.Series:
        """Will work on n_meth and n_unmeth"""
        cum_ser = ser.cumsum()
        return cum_ser[-1] - cum_ser
    res = (mbias_stats_df_cg_flen_agg
           .groupby(non_flen_phred_levels)
           .transform(compute_phred_threshold_counts_for_group)
           .assign(beta_value=compute_beta_values)
           )

    # Discard the highest phred bin - it has no events left after filtering
    # noinspection PyUnusedLocal
    # pylint: disable=unused-variable
    phred_idx = res.index.get_level_values("phred").unique()[:-2]  # pylint: disable=unused-variable
    res = res.query("phred in @phred_idx")

    return res


def compute_beta_values(df: pd.DataFrame) -> pd.Series:
    return df['n_meth'] / (df['n_meth'] + df['n_unmeth'])


def save_df_to_trunk_path(df: pd.DataFrame, trunk_path: str) -> None:
    print(f"Saving {op.basename(trunk_path)}")
    print("Saving pickle")
    df.to_pickle(trunk_path + '.p')
    print("Saving feather")
    df.reset_index().to_feather(trunk_path + '.feather')


def cutting_sites_df_has_correct_format(cutting_sites_df: pd.DataFrame) -> bool:
    # False positive for all() call
    # noinspection PyUnresolvedReferences
    return (not cutting_sites_df.empty
            and list(cutting_sites_df.columns.values) == ['start', 'end']
            and (cutting_sites_df.dtypes == np.int64).all()
            and ['bs_strand', 'flen'] == list(cutting_sites_df.index.names)
            and cutting_sites_df.columns.name == 'cutting_site'
            and cutting_sites_df.index.get_level_values('bs_strand').dtype.name == 'category'
            )
# Adding more index levels is in general not a problem
# Some methods would have to be adapted. With no guarantee of completeness
# - as_array()
# @invariant('Meets cutting site df format',
#            lambda inst: cutting_sites_df_has_correct_format(inst.df))
class CuttingSites:
    """CuttingSites for M-bias trimming

    Start and end indicate the python slice containing the plateau
    positions. I.e. all zero-based positions in [start, end[ are in
    the plateau.

    See the invariant checks for further format specs.
    """

    def __init__(self, cutting_sites_df: pd.DataFrame, max_read_length: int) -> None:
        self.df = cutting_sites_df
        self.max_read_length = max_read_length

    @staticmethod
    def from_mbias_stats(mbias_stats_df: pd.DataFrame, **kwargs) -> 'CuttingSites':
        """Detect CuttingSites in M-bias stats

        Args:
            **kwargs: passed on to the plateau detection algorithm.
                The choice of algorithm is not yet implemented. At the
                moment, this is directly passed on to the default choice
                (:class:`.BinomPvalueBasedCuttingSiteDetermination`)


        """

        try:
            cutting_sites_df = BinomPvalueBasedCuttingSiteDetermination(
                mbias_stats_df=mbias_stats_df, **kwargs).compute_cutting_site_df()
        except TypeError:
            print(f'The plateau detection parameters {kwargs} do not match the'
                  f'signature of the selected plateau detection algorithm')
            raise


        # this line is duplicated in BinomPvalueBasedCuttingSiteDetermination
        # should be moved to - yet to be done - MbiasStats class
        max_read_length = mbias_stats_df.index.get_level_values('pos').max()
        return CuttingSites(
            cutting_sites_df=cutting_sites_df,
            max_read_length=max_read_length)

    @staticmethod
    def from_rel_to_frag_end_cutting_sites(
            cut_site_spec: Mapping[str, Tuple[int, int]], max_read_length: int) \
            -> 'CuttingSites':
        """Construct from #bp to be removed from the fragment ends

        The maximum fragment length requiring a specific definition
        of the plateau positions is automatically determined
        based on the max_read_length and the cut_site_spec.
        (All fragment lengths larger than the maximum fragment length
        defined in a CuttingSites object are automatically
        treated like the maximum defined fragment length.)

        Examples:
            >>> cut_site_spec = dict(
            >>>     w_bc = (0, 9),
            >>>     c_bc = (0, 9),
            >>>     w_bc_rv = (9, 0),
            >>>     c_bc_rv = (9, 0),
            >>> )
            >>> cutting_sites = CuttingSites.from_rel_to_frag_end_cutting_sites(
            >>>     cut_site_spec=cut_site_spec, max_read_length=125)
        """

        max_n_bp_removed_from_end = max([elem[1] for elem in cut_site_spec.values()])
        max_flen = max_read_length + max_n_bp_removed_from_end

        midx = pd.MultiIndex.from_product([
            pd.Categorical(bstrand_idxs._fields, ordered=True),
            range(0, max_flen + 1)], names=['bs_strand', 'flen'])
        end_col = np.tile(list(range(0, max_flen + 1)), 4).astype('i8')
        cutting_sites_df = pd.DataFrame(
            {'start': -1, 'end': end_col}, index=midx)
        cutting_sites_df.columns.name = 'cutting_site'

        for strand in bstrand_idxs._fields:
            cutting_sites_df.loc[strand, 'start'] = cut_site_spec[strand][0]
            cutting_sites_df.loc[[strand], 'end'] = (cutting_sites_df.loc[[strand], 'end'] -
                                                    cut_site_spec[strand][1])
        bad_fragment_pos_are_outside_of_read_length = cutting_sites_df['end'] >= max_read_length
        cutting_sites_df.loc[bad_fragment_pos_are_outside_of_read_length, 'end'] = \
            max_read_length
        cutting_sites_df.loc[cutting_sites_df['start'] >= cutting_sites_df['end'], :] = 0
        # implies:
        # cutting_sites_df.loc[cutting_sites_df['end'] < 0, 'end'] = 0
        # cutting_sites_df.loc[cutting_sites_df['end'] == 0, 'start'] = 0

        return CuttingSites(
            cutting_sites_df=cutting_sites_df,
            max_read_length=max_read_length
        )

    def as_array(self) -> np.ndarray:
        full_midx = pd.MultiIndex.from_product([
            pd.Categorical(bstrand_idxs._fields, ordered=True),
            range(self.df.index.get_level_values('flen').min(),
                  self.df.index.get_level_values('flen').max() + 1)
        ], names=['bs_strand', 'flen'])
        df = self.df.reindex(index=full_midx).fillna(method='bfill')

        df = df.astype('i8', errors='raise', copy=True)
        df = df.stack().to_frame('cut_pos').reset_index()

        df['bs_strand'].cat.categories = range(4)
        df['cutting_site'] = df['cutting_site'].replace({'start': 0, 'end': 1})

        max_flen = df['flen'].max()
        arr = np.zeros((4, (max_flen + 1), 2), dtype='i8')
        arr[df['bs_strand'], df['flen'], df['cutting_site']] = df['cut_pos']

        return arr

    def plot(self, path: str) -> None:
        """Create CuttingSites line plot

        Plot start and end values for all fragment lengths, with
        BS-Seq strands facetted across the columns.

        The CuttingSites may optionally be stratified by motif, i.e.
        a motif index may be present or absent, and it may have arbitrary
        levels. If a motif index is present, the motifs will be displayed
        via row facetting.
        """
        # TODO-minor: adjust the start and end values to account for
        # their definition as slice boundaries
        if 'motif' in self.df.index.names:
            row_facetting_dict = {'row': 'motif'}
        else:
            row_facetting_dict = {}
        plot_df = self.df.copy()
        plot_df.columns.name = 'cut_end'
        plot_df = plot_df.stack().to_frame('cut_pos').reset_index()
        (alt.Chart(plot_df)
         .mark_line()
         .encode(x='flen:Q', y='cut_pos:Q', color='cut_end:N')
         .facet(column='bs_strand:N',
                **row_facetting_dict
                )
         ).interactive().save(path, webdriver='firefox')


class BinomPvalueBasedCuttingSiteDetermination:
    """Detect unlikely deviations from estimated global methylation level

    Often, even if M-bias is present, it mainly affects certain
    BS-Seq strands, and certain fragment lengths. This algorithm
    first estimates the unbiased global methylation level based on the
    set of fragment lengths and BS-Seq strands defined as 'unbiased'
    by **plateau_flen** and *plateau_bs_strands**. Then it detects
    significant deviations from this plateau and translates them
    into BS-Seq strand and fragment length-specific cutting sites.

    Args:
        min_plateau_length: Minimum number of bp unaffected of M-bias
            required to include a stratum into the methylation calls.
        allow_slope: Whether the M-bias curve is allowed to have a slope.
        max_slope: Maximum allowed slope of
            the M-bias curve, if slope is allowed.
        plateau_flen: Fragment lengths >= plateau_flen are
            used to estimate the correct global methylation level
        plateau_bs_strands: one or more of [w_bc c_bc
            w_bc_rv c_bc_rv] Only these BS-Seq strands are used
            to estimate the correct global methylation level
        always_accept_distance_from_plateau: beta values within
           estimated global methylation level +- always_accept_distance_from_plateau
           are never rejected
        plateau_detection_step_size: increase to improve the
            estimate of the true global methylation level and slope of
            the M-bias curve. Computation time increases with O(n^3).
    """

    def __init__(self, mbias_stats_df: pd.DataFrame, allow_slope: bool = True,
                 min_plateau_length: int = 30, max_slope: float = 0.0006,
                 plateau_flen: int = 210,
                 plateau_bs_strands: Tuple[str, ...] = ('c_bc', 'w_bc'),
                 always_accept_distance_from_plateau: float = 0.02,
                 plateau_detection_step_size: int = 1000) -> None:
        self.plateau_flen = plateau_flen
        self.max_slope = max_slope
        # this line is duplicated in CuttingSites.from_mbias_stats
        # should be moved to - yet to be done - MbiasStats class
        self.max_read_length = mbias_stats_df.index.get_level_values('pos').max()
        self.plateau_bs_strands = plateau_bs_strands
        self.min_plateau_length = min_plateau_length
        self.allow_slope = allow_slope
        self.mbias_stats_df = mbias_stats_df
        self.always_accept_distance_from_plateau = always_accept_distance_from_plateau
        self.plateau_detection_step_size = plateau_detection_step_size


    def compute_cutting_site_df(self) -> pd.DataFrame:

        fn_mbias_stats_df = self.mbias_stats_df.copy()
        fn_mbias_stats_df = fn_mbias_stats_df.fillna(0)
        n_meth_wide = fn_mbias_stats_df['n_meth'].unstack(level='pos', fill_value=0)
        n_total_wide = fn_mbias_stats_df['n_total'].unstack(level='pos', fill_value=0)
        beta_values_wide = fn_mbias_stats_df['beta_value'].unstack('pos', fill_value=0)


        intercept_, slope_ = self._estimate_plateau_height(
                fn_mbias_stats_df, allow_slope=self.allow_slope)
        plateau_heights_ = (intercept_ + np.arange(n_meth_wide.shape[1]) * slope_)
        strip_low_bound = plateau_heights_ - self.always_accept_distance_from_plateau
        strip_high_bound = plateau_heights_ + self.always_accept_distance_from_plateau

        # Currently impossible because we only search up until the max. slope
        # but search may go beyond max slope in the future, so let's be defensive here
        if slope_ > self.max_slope:
            print('The estimated M-bias slope of this sample exceeds the threshold set '
                  'by you. This sample can\'t be processed further with the current '
                  'algorithm and parameters.')
            raise ValueError

        p_values_wide_low = pd.DataFrame(self._get_p_values(k=n_meth_wide,
                                                            n=n_total_wide,
                                                            p=strip_low_bound),
                                         index=n_meth_wide.index,
                                         columns=n_meth_wide.columns)
        p_values_wide_high = pd.DataFrame(self._get_p_values(k=n_meth_wide,
                                                             n=n_total_wide,
                                                             p=strip_high_bound),
                                          index=n_meth_wide.index,
                                          columns=n_meth_wide.columns)
        p_values_wide_integrated = self.integrate_p_values(
            beta_values=beta_values_wide,
            high_p_values=p_values_wide_high,
            low_p_values=p_values_wide_low,
            strip_low_bound=strip_low_bound,
            strip_high_bound=strip_high_bound,
        )

        cutting_sites_list = []
        try:
            flen_idx: Optional[int] = list(fn_mbias_stats_df.index.names).index('flen')
        except ValueError:
            flen_idx = None
        for i in range(p_values_wide_integrated.shape[0]):
            if flen_idx:
                flen = p_values_wide_integrated.index[i][flen_idx]
            else:
                flen = None
            cutting_sites_list.append(self._find_breakpoints2(
                # neighbors=predecessor_p_values_integrated.iloc[i, :],
                neighbors=p_values_wide_integrated.iloc[i, :],
                row=p_values_wide_integrated.iloc[i, :], flen=flen))
        df = pd.DataFrame(cutting_sites_list, index=n_meth_wide.index)

        group_levels = list(df.index.names)
        group_levels.remove('flen')
        # bug in groupby in pandas 0.23, need to take elaborate construct
        min_flen = df.index.get_level_values('flen')[0]
        def fn(df: pd.DataFrame) -> pd.DataFrame:
            window_size = 21
            return (df
                    .rolling(window=window_size, center=True, min_periods=1)
                    .apply(self._smooth_cutting_sites, raw=False,
                           kwargs=dict(window_size=window_size,
                                       min_flen=min_flen)))
        df = (df.groupby(group_levels, group_keys=False)
              .apply(fn))

        df.loc[df['end'] - df['start'] < self.min_plateau_length, :] = 0
        df = df.astype('i8')
        df.columns.name = 'cutting_site'

        return df


    @staticmethod
    def integrate_p_values(beta_values: pd.DataFrame, high_p_values: pd.DataFrame,
                           low_p_values: pd.DataFrame, strip_low_bound: pd.DataFrame,
                           strip_high_bound: pd.DataFrame) -> pd.DataFrame:
        res = high_p_values.copy()
        max_mask = low_p_values > high_p_values
        res[max_mask] = low_p_values[max_mask]
        in_window_mask = beta_values.apply(
            lambda ser: ser.between(strip_low_bound, strip_high_bound),
            axis=1
        )
        res[in_window_mask] = 1
        return res

    @staticmethod
    def _get_p_values(k: np.ndarray, n: np.ndarray, p: np.ndarray) -> np.ndarray:
        binom_prob = scipy.stats.binom.cdf(k=k, n=n, p=p)
        binom_prob[binom_prob > 0.5] = 1 - binom_prob[binom_prob > 0.5]
        return binom_prob * 2

    def _estimate_plateau_height(self, mbias_stats_df: pd.DataFrame, allow_slope: bool) \
            -> Tuple[float, float]:

        mbias_curve_for_plateau_detection = (mbias_stats_df
                                             .loc[nidxs(bs_strand=self.plateau_bs_strands), :]
                                             .query('flen >= @self.plateau_flen')
                                             .groupby(level='pos').sum())
        mbias_curve_for_plateau_detection = (
            mbias_curve_for_plateau_detection.assign(
                beta_value=compute_beta_values))

        if allow_slope:
            intercept_, slope_ = self._estimate_plateau_with_slope(
                n_meth=mbias_curve_for_plateau_detection['n_meth'],
                n_total=mbias_curve_for_plateau_detection['n_total'],
                step_size=self.plateau_detection_step_size)
        else:
            intercept_, slope_ = self._estimate_plateau(
                n_meth=mbias_curve_for_plateau_detection['n_meth'],
                coverage_arr=mbias_curve_for_plateau_detection['n_total'])

        return intercept_, slope_

    def _smooth_cutting_sites(self, ser: pd.Series, window_size: int, min_flen: int) -> float:

        # rolling windows don't prepend/append nan at the edges of the df
        # have to prepend ourselves, makes following logic easier
        if ser.shape[0] < window_size:
            if ser.index.get_level_values('flen')[0] == min_flen:
                arr = np.concatenate([np.tile(np.nan, window_size - ser.shape[0]),
                                      ser.values])
            else:
                arr = np.concatenate([ser.values,
                                      np.tile(np.nan, window_size - ser.shape[0])])
        else:
            arr = ser.values

        window_half_size = np.int(window_size/2 - 1/2)
        left_repr_value, unused_left_repr_value_counts = self.get_representative_value(
            arr[0:window_half_size])
        right_repr_value, unused_right_repr_value_counts = self.get_representative_value(
            arr[-window_half_size:])
        if left_repr_value == right_repr_value:
            return left_repr_value
        return arr[window_half_size]

    @staticmethod
    def get_representative_value(arr: np.ndarray) -> Tuple[float, int]:
        arr = arr[~np.isnan(arr)].astype(int)
        repr_value = np.nan
        repr_value_counts = np.nan
        if arr.shape[0] >= 3:
            counts = np.bincount(arr)
            max_counts = np.max(counts)
            n_max = np.sum(counts == max_counts)
            if n_max == 1 and max_counts >= 3:
                repr_value = np.argmax(counts)
                repr_value_counts = max_counts
        return repr_value, repr_value_counts

    # noinspection PyUnusedLocal
    # pylint: disable=unused-argument
    @staticmethod
    def _find_breakpoints2(neighbors: pd.Series, row: pd.Series,
                           p_value_threshold: float = 10**-6,
                           n_consecutive_plateau_points: int = 5,
                           flen: Optional[int] = None) -> pd.Series:


        def go_from_seed_to_plateau(seed_idx: int, direction: str,
                                    p_value_threshold: float = 0.1) -> int:
            if direction == 'forward':
                step = 1
            else:
                step = -1
            n_good_points = 0
            for i in count(seed_idx + step, step):
                try:
                    if neighbors.iloc[i] > p_value_threshold:
                        n_good_points += 1
                        if n_good_points > n_consecutive_plateau_points:
                            if direction == 'forward':
                                return i - n_consecutive_plateau_points
                            # direction == reverse
                            return neighbors.shape[0] + (i + n_consecutive_plateau_points) + 1
                    else:
                        n_good_points = 0
                except IndexError:
                    return -1

            # TODO-minor: improve by using more explicit loop
            # this return is just to make linters happy, it cannot be reached
            return -1

        if neighbors.iloc[0] < p_value_threshold:
            left_bkp = go_from_seed_to_plateau(0, 'forward')
            if left_bkp == -1:
                return pd.Series(dict(start=0, end=0))
        else:
            left_bkp = 0

        if flen and flen >= neighbors.shape[0]:
            last_pos_idx = -1
        elif flen and flen < neighbors.shape[0]:
            last_pos_idx = -(neighbors.shape[0] - flen + 1)
        else:
            last_pos_idx = -1

        # TODO-plateau_detection: coverage of zero -> p value of 1
        # TODO-plateau_detection: skip nan inwards
        if neighbors.iloc[last_pos_idx] < p_value_threshold:
            right_bkp = go_from_seed_to_plateau(last_pos_idx, 'backward')
            if right_bkp == -1:
                return pd.Series(dict(start=0, end=0))
        else:
            right_bkp = neighbors.shape[0]

        return pd.Series((dict(start=left_bkp, end=right_bkp)))

    @staticmethod
    def _estimate_plateau(n_meth: np.ndarray, coverage_arr: np.ndarray) -> Tuple[float, float]:
        plateau_heights = np.linspace(0, 1, 1000)
        p_value_mat = scipy.stats.binom.cdf(k=n_meth[np.newaxis, :],
                                            n=coverage_arr[np.newaxis, :],
                                            p=plateau_heights[:, np.newaxis])
        p_value_mat[p_value_mat > 0.5] = 1 - p_value_mat[p_value_mat > 0.5]
        # TODO-plateau_detection: only take elements times 2 where it makes sense
        p_value_mat *= 2
        log_p_value_mat = np.log10(p_value_mat + 10**-30)
        # TODO-plateau_detection: discarding rows if more than half of the row is bad may be too harsh
        row_is_ok = np.sum(log_p_value_mat < -4, axis=1) < n_meth.shape[0] / 2
        log_p_value_mat = log_p_value_mat[row_is_ok, :]
        remaining_indices = np.arange(0, 1000)[row_is_ok]
        log_p_value_mat[log_p_value_mat < -8] = 0

        try:
            plateau_height_ = remaining_indices[np.argmax(np.sum(log_p_value_mat, axis=1))] / 1000
        except ValueError:
            # if log_p_value_mat is empty due to filtering above
            print('Cannot identify a plateau with binomp algorithm. Change parameters'
                  'or use another plateau detection algorithm. One place to start would'
                  'be allowing a slope.')
            raise

        return plateau_height_, 0

    def _estimate_plateau_with_slope(self, n_meth: np.ndarray, n_total: np.ndarray,
                                     step_size: int = 1000, inclusion_log_p_threshold = -8)\
            -> Tuple[float, float]:

        # TODO-plateau_detection: this will raise warning on zero divide and then
        # fail at np.median. problematic when in read and when there are zero
        # coverage positions at the end of the read
        beta_values = n_meth / n_total
        initial_intercept_ = np.median(n_meth / n_total)
        intercept_search_window = np.linspace(initial_intercept_ - 0.1, initial_intercept_ + 0.1, step_size)
        slope_search_window = np.linspace(0, self.max_slope, step_size)
        pos = np.arange(n_meth.shape[0])
        added_by_slope = slope_search_window[:, np.newaxis] * pos[np.newaxis, :]
        model_beta_value_arr = intercept_search_window[:, np.newaxis, np.newaxis] + added_by_slope[np.newaxis, :, :]
        model_beta_value_arr[model_beta_value_arr > 1] = 1
        model_beta_value_arr[model_beta_value_arr < 0] = 0

        p_value_mat = scipy.stats.binom.cdf(k=n_meth[np.newaxis, np.newaxis, :],
                                            n=n_total[np.newaxis, np.newaxis, :],
                                            p=model_beta_value_arr)
        p_value_mat[p_value_mat > 0.5] = 1 - p_value_mat[p_value_mat > 0.5]
        p_value_mat *= 2
        log_p_value_mat = np.log10(p_value_mat + 10**-30)

        beta_value_arr = n_meth[np.newaxis, np.newaxis, :] / n_total[np.newaxis, np.newaxis, :]
        is_within_strip = np.abs(beta_value_arr - model_beta_value_arr) < self.always_accept_distance_from_plateau
        log_p_value_mat[is_within_strip] = 0
        log_p_value_mat[log_p_value_mat < inclusion_log_p_threshold] = np.nan

        # TODO: don't use point in strip, use point above inclusion p threshold
        try:
            n_ok_pos_arr = np.nansum(log_p_value_mat > inclusion_log_p_threshold, axis=2)
            if np.max(n_ok_pos_arr) < n_meth.shape[0] / 2:
                raise  ValueError
            ph_idx, slope_idx = np.unravel_index(np.argmax(n_ok_pos_arr), n_ok_pos_arr.shape)
            temp_intercept_ = intercept_search_window[ph_idx]
            slope_ = slope_search_window[slope_idx]
            model_beta_values = temp_intercept_ + np.arange(0, n_meth.shape[0]) * slope_
            beta_is_ok = np.abs(beta_values - model_beta_values) < self.always_accept_distance_from_plateau
            final_intercept_ = np.median(beta_values[beta_is_ok])
        except ValueError:
            # argmax of empty sequence -> no plateau found
            print('Cannot identify a plateau with binomp algorithm. Change parameters'
                  'or use another plateau detection algorithm.')
            raise

        print(final_intercept_, slope_)
        return final_intercept_, slope_

        # row_is_ok = np.sum(log_p_value_mat < inclusion_log_p_threshold, axis=2) < n_meth.shape[0] / 2
        # log_p_value_mat[log_p_value_mat < inclusion_log_p_threshold] = np.nan
        # # get probability of observations by adding up the log probabilities
        # sums = np.nansum(log_p_value_mat, axis=2)
        # idx = np.nonzero(row_is_ok)
        # ok_events = sums[row_is_ok]
        # #
        # try:
        #     plateau_height = intercept_search_window[idx[0][np.argmax(ok_events)]]
        #     slope = slope_search_window[idx[1][np.argmax(ok_events)]]
        # except ValueError:
        #     # argmax of empty sequence -> no plateau found
        #     print('Cannot identify a plateau with binomp algorithm. Change parameters'
        #           'or use another plateau detection algorithm.')
        #     raise
        # return plateau_height, slope

        # remaining_indices = np.arange(0, 1000)[row_is_ok]
        # log_p_value_mat[log_p_value_mat < -8] = 0
        # plateau_height_ = remaining_indices[np.argmax(np.sum(log_p_value_mat, axis=1))] / 1000

    @staticmethod
    def _convert_mbias_curve_to_pvalues(n_meth: np.ndarray, coverage_arr: np.ndarray,
                                        plateau_heights_: Union[int, np.ndarray]) -> np.ndarray:

        p_values = scipy.stats.binom.cdf(n=coverage_arr,
                                         p=plateau_heights_, k=n_meth)
        p_values[p_values > 0.5] = 1 - p_values[p_values > 0.5]
        p_values *= 2
        return p_values
