"""Test M-bias module"""
import json
import os
import pickle
import shutil
import tempfile
from collections import namedtuple, defaultdict, Sequence
from io import StringIO
from pathlib import Path
from textwrap import dedent
from unittest.mock import MagicMock, patch, call
from typing import List, Tuple, Union, Callable, DefaultDict, Dict, Any, cast
import pandas as pd
import toml
import toolz as tz


idxs = pd.IndexSlice
import numpy as np
import pytest
import subprocess

from itertools import product

from mqc.config import assemble_config_vars
from mqc.pileup.pileup import MotifPileup
from mqc.mbias import (
    MbiasCounter, get_sequence_context_to_array_index_table,
    fit_normalvariate_plateau,
    fit_percentiles,
    mask_mbias_stats_df, map_seq_ctx_to_motif,
    compute_beta_values,
    get_plot_configs, MbiasPlotConfig,
    aggregate_table, create_aggregated_tables,
    MbiasPlotParams, MbiasPlotMapping, MbiasAnalysisConfig,
    AggregatedMbiasStats,
    MbiasStatAggregationParams, CuttingSites,
    BinomPvalueBasedCuttingSiteDetermination,
    cutting_sites_df_has_correct_format,
)
from mqc.utils import get_resource_abspath, NamedIndexSlice, subset_dict, TmpChdir
import mqc.filepaths
import mqc.flag_and_index_values as mfl

# import matplotlib
# matplotlib.use('Agg')  # import before pyplot import!
# from matplotlib.axes import Axes
# import matplotlib.pyplot as plt
# import seaborn as sns



# TODO: I currently test that any non-zero qc_fail_flag leads to discard from M-bias
# stats counting. When I update the behavior so that phred score fails are kept
# in the stats, the tests here also need to be updated accordingly
b_inds = mfl.bsseq_strand_indices
b_na_ind = mfl.bsseq_strand_na_index
m_flags = mfl.methylation_status_flags
qc_flags = mfl.qc_fail_flags

CONFIG: DefaultDict[str, Dict[str, Any]] = defaultdict(dict)
CONFIG['paths']['mbias_counts'] = None
CONFIG['data_properties']['max_read_length_bp'] = 101
CONFIG['stats']['max_flen'] = 500
CONFIG['stats']['max_flen_with_single_flen_resolution'] = 150
CONFIG['stats']['flen_bin_size'] = 10
CONFIG['stats']['max_phred'] = 40
CONFIG['stats']['phred_bin_size'] = 5
CONFIG["stats"]["seq_context_size"] = 5
SAMPLE_NAME = 'hsc_1'

AlignmentStub = namedtuple('AlignmentStub',
                           'template_length')

IndexPositionStub = namedtuple('IndexPositionStub',
                               'seq_context')

SEQ_CONTEXT_TO_IDX_MAPPING = {
    'CCCCC': 0,
    'GGCGG': 1,
    'AACAA': 2,
    'TTCTT': 2,
}
BINNED_SEQ_CONTEXT_TO_IDX_MAPPING = {
    'CCCCC': 0,
    'GGCGG': 1,
    'WWCWW': 2,
}


class BSSeqPileupReadStub:
    """Stub for BSSeqPileupRead with AlignmentStub"""
    def __init__(self, strand_idx, flen, pos, mflag, qc_fail_flag, phred,
                 expected_flen_bin, expected_phred_bin, phred_fail_flag):
        self.qc_fail_flag = qc_fail_flag
        self.phred_fail_flag = phred_fail_flag
        self.bsseq_strand_ind = strand_idx
        self.meth_status_flag = mflag
        self.alignment = AlignmentStub(template_length=flen)
        self.pos_in_read = pos
        self.baseq_at_pos = phred
        self.expected_phred_bin = expected_phred_bin
        self.expected_flen_bin = expected_flen_bin

    def get_idx_tuple(self):
        """Get tuple for the M-bias counts array element of this read"""
        if self.meth_status_flag == m_flags.is_methylated:
            meth_status_index = 0
        elif self.meth_status_flag == m_flags.is_unmethylated:
            meth_status_index = 1
        else:
            raise ValueError("This is a failed read, "
                             "can't give event_class indices")

        return (self.bsseq_strand_ind,
                self.expected_flen_bin,
                self.expected_phred_bin,
                self.pos_in_read,
                meth_status_index)


bstub = BSSeqPileupReadStub


def bstub_with(**kwargs):
    event_class = dict(strand_idx=b_inds.c_bc,
                       flen=180,
                       expected_flen_bin=153,
                       pos=10,
                       mflag=m_flags.is_methylated,
                       qc_fail_flag=0,
                       phred_fail_flag=0,
                       phred=30,
                       expected_phred_bin=6)
    event_class.update(kwargs)
    return bstub(**event_class)


class MotifPileupStub:
    def __init__(self, seq_context, reads):
        self.idx_pos = IndexPositionStub(seq_context)
        self.reads = reads


class TestSeqContextToBinMapping:
    def test_even_number_motif_sizes_are_not_allowed(self):
        with pytest.raises(ValueError):
            get_sequence_context_to_array_index_table(4)

    @pytest.mark.parametrize("N", [3, 5])
    def test_seq_ctxs_are_mapped_to_idxs_in_alpha_order_of_binned_motifs(
            self, N):
        """Sequence contexts are mapped to bin numbers, and bin numbers
        are given according to the alphabetical order of the binned
        representations of the sequence contexts"""

        n_possible_cgw_motifs = 3 ** (N - 1)
        seq_ctx_to_idx, _ = get_sequence_context_to_array_index_table(
            N)
        if N == 3:
            motif_start = "CCC"
            motif_end = "TCT"
        else:  # N == 5
            motif_start = "CCCCC"  # binned: "CCCCC"
            motif_end = "TTCTT"  # binned: "WWCWW"
        assert seq_ctx_to_idx[motif_start] == 0
        assert seq_ctx_to_idx[motif_end] == n_possible_cgw_motifs - 1

    @pytest.mark.parametrize("N", [3, 5])
    def test_seq_context_does_not_map_illegal_motifs(self, N):
        seq_ctx_to_idx, _ = (
            get_sequence_context_to_array_index_table(N))
        if N == 3:
            bad_motif = "NCT"
        else:  # N == 5
            bad_motif = "ANCTA"
        with pytest.raises(KeyError):
            assert seq_ctx_to_idx[bad_motif]

    def test_equivalent_motifs_are_mapped_to_the_same_bin(self):
        """Because motifs are binned, several full sequence contexts may
        correspond to the same binned motif, and thus map to the same
        array index"""
        seq_ctx_to_idx, _ = (
            get_sequence_context_to_array_index_table(5))
        assert seq_ctx_to_idx["AACAA"] == seq_ctx_to_idx["TTCAA"]
        assert seq_ctx_to_idx["AACAA"] == seq_ctx_to_idx["ATCAT"]


@pytest.mark.parametrize("seq_ctx,  use_classical, exp_motif",
                         [
                             ("ACA", True, "CH"),
                             ("ACA", False, "CW"),
                             ("TCG", True, "CG"),
                             ("TCG", False, "CG"),
                             ("TCC", True, "CH"),
                             ("TCC", False, "CC"),
                             ("TCW", False, "CW"),
                             ("GGCGG", True, "CG"),
                             ("GGCGG", False, "CG"),
                             ("AACAA", True, "CHH"),
                             ("AACAA", False, "CWW"),
                             ("AACAG", True, "CHG"),
                             ("AACAG", False, "CWG"),
                             ("AACCG", True, "CHG"),
                             ("AACCG", False, "CCG"),
                             ("AACWW", True, "CHH"),
                             ("AACWW", False, "CWW"),
                             ("AACWG", False, "CWG"),
                             ("AACWG", True, "CHG"),
                             ("TACCW", True, "CHH"),
                             ("TACCW", False, "CCW"),
                             ("CGGCGGA", True, "CG"),
                             ("AAACAGA", True, "CHG"),
                             ("AAACAGA", False, "CWG"),
                         ])
def test_seq_ctx_to_motif_mapping(seq_ctx, use_classical, exp_motif):
    mapped_motif = map_seq_ctx_to_motif(seq_ctx, use_classical=use_classical)
    assert mapped_motif == exp_motif


def test_mbias_counter_get_dataframe(mocker):
    """This  is a test against the combination of MbiasCounter.__init__
    and get_dataframe, since they are closely coupled"""

    config = defaultdict(dict)
    config["data_properties"]["max_read_length_bp"] = 10
    config["paths"]["mbias_counts"] = None
    config["stats"]["max_flen"] = 10
    config["stats"]["max_flen_with_single_flen_resolution"] = 5
    config["stats"]["flen_bin_size"] = 2
    config["stats"]["max_phred"] = 10
    config["stats"]["phred_bin_size"] = 5
    config["stats"]["seq_context_size"] = 5

    map_fn = mocker.patch('mqc.mbias.'
                          'get_sequence_context_to_array_index_table')
    map_fn.return_value = (SEQ_CONTEXT_TO_IDX_MAPPING,
                           BINNED_SEQ_CONTEXT_TO_IDX_MAPPING)

    flen_levels = [0, 1, 2, 3, 4, 5, 7, 9, 10]
    phred_levels = [4, 9, 10]
    midx = pd.MultiIndex.from_product([
        pd.Categorical(['CCCCC', 'GGCGG', 'WWCWW'], ordered=True),
        pd.Categorical(['c_bc', 'c_bc_rv', 'w_bc', 'w_bc_rv'], ordered=True),
        flen_levels,
        phred_levels,
        list(range(1, 11)),
        pd.Categorical(["n_meth", "n_unmeth"], ordered=True)
    ], names=["seq_context", "bs_strand", "flen", "phred", "pos",
              "meth_status"])
    exp_df = pd.DataFrame(index=midx)
    # counts must be uint64. But if we set this now, we have to add
    # every count below as uint64 to avoid type casting
    # better to do manual type casting to uint64 once at the end
    exp_df["counts"] = 0

    mbias_counter = MbiasCounter(config)

    # random middle values
    exp_df.loc[idxs["WWCWW", "w_bc", 7, 4, 5, "n_meth":"n_meth"], "counts"] = 1
    mbias_counter.counter_array[2, 2, 6, 0, 4, 0] = 1
    exp_df.loc[idxs["CCCCC", "c_bc", 9, 9, 3, "n_unmeth":"n_unmeth"],
               "counts"] = 2
    mbias_counter.counter_array[0, 0, 7, 1, 2, 1] = 2
    # max values
    exp_df.loc[idxs["CCCCC", "c_bc_rv", 10, 10, 10, "n_unmeth":"n_unmeth"],
               "counts"] = 1
    mbias_counter.counter_array[0, 1, -1, -1, -1, 1] = 1
    # min values
    exp_df.loc[idxs["GGCGG", "w_bc_rv", 0, 4, 1, "n_meth":"n_meth"],
               "counts"] = 1
    mbias_counter.counter_array[1, 3, 0, 0, 0, 0] = 1

    computed_df = mbias_counter.get_dataframe()
    # noinspection PyTypeChecker
    pd.testing.assert_frame_equal(exp_df.astype('u8'), computed_df)


class TestMbiasCounterMotifPileupProcessing:
    def test_phred_flen_seqcontext_in_allowed_range_are_binned(
            self, mocker) -> None:
        reads = [bstub(strand_idx=b_inds.w_bc,
                       flen=100,
                       expected_flen_bin=100,
                       pos=10,
                       mflag=m_flags.is_methylated,
                       qc_fail_flag=0,
                       phred=20,
                       expected_phred_bin=4,
                       phred_fail_flag=0),
                 bstub(strand_idx=b_inds.c_bc,
                       flen=180,
                       expected_flen_bin=153,
                       pos=10,
                       mflag=m_flags.is_methylated,
                       qc_fail_flag=0,
                       phred=30,
                       expected_phred_bin=6,
                       phred_fail_flag=0), ]

        map_fn = mocker.patch('mqc.mbias.'
                              'get_sequence_context_to_array_index_table')
        map_fn.return_value = (
            SEQ_CONTEXT_TO_IDX_MAPPING, BINNED_SEQ_CONTEXT_TO_IDX_MAPPING)

        mbias_counter = MbiasCounter(CONFIG)
        mbias_counter.seq_ctx_idx_dict = SEQ_CONTEXT_TO_IDX_MAPPING

        seq_context1 = 'CCCCC'
        seq_context1_idx = SEQ_CONTEXT_TO_IDX_MAPPING[seq_context1]
        seq_context2 = 'GGCGG'
        seq_context2_idx = SEQ_CONTEXT_TO_IDX_MAPPING[seq_context2]

        motif_pileup1 = MotifPileupStub(seq_context=seq_context1,
                                        reads=reads)
        motif_pileup2 = MotifPileupStub(seq_context=seq_context2,
                                        reads=reads)
        mbias_counter.process(cast(MotifPileup, motif_pileup1))
        mbias_counter.process(cast(MotifPileup, motif_pileup2))

        read_strata_for_advanced_indexing = [
            read.get_idx_tuple() for read in reads]

        full = [(curr_seq_ctx,) + read_tuple
                for curr_seq_ctx in [seq_context1_idx, seq_context2_idx]
                for read_tuple in read_strata_for_advanced_indexing]

        full2 = list(zip(*full))

        target_subarray = mbias_counter.counter_array[full2]
        assert (target_subarray == 1).all()
        mbias_counter.counter_array[full2] = 0
        assert (mbias_counter.counter_array == 0).all()

    def test_unusable_reads_are_discarded(self, mocker) -> None:
        """Test that unreliable reads are discarded

        Reads are discarded if
            - they have a qc_fail_flag
            - they have no methylation calling status: NA, SNP or Ref

        Reads are not discared if they have a phred_fail_flag, because
        the phred score is recorded as dimension of the M-bias stats
        """

        map_fn = mocker.patch('mqc.mbias.'
                              'get_sequence_context_to_array_index_table')
        map_fn.return_value = (
            SEQ_CONTEXT_TO_IDX_MAPPING, BINNED_SEQ_CONTEXT_TO_IDX_MAPPING)

        reads = [bstub_with(),
                 bstub_with(phred_fail_flag=1),
                 bstub_with(mflag=m_flags.is_ref),
                 bstub_with(mflag=m_flags.is_snp),
                 bstub_with(mflag=m_flags.is_na),
                 bstub_with(qc_fail_flag=qc_flags.overlap_fail),
                 bstub_with(qc_fail_flag=qc_flags.sam_flag_fail),
                 bstub_with(qc_fail_flag=qc_flags.mapq_fail),
                 ]

        motif_pileup = MotifPileupStub(seq_context='CCCCC',
                                       reads=reads)

        mbias_counter = MbiasCounter(CONFIG)

        mbias_counter.process(cast(MotifPileup, motif_pileup))

        base_event_class = ((SEQ_CONTEXT_TO_IDX_MAPPING['CCCCC'],)
                            + bstub_with().get_idx_tuple())
        assert mbias_counter.counter_array[base_event_class] == 2
        mbias_counter.counter_array[base_event_class] = 0
        assert (mbias_counter.counter_array == 0).all()

    def test_flen_and_phred_are_pooled_in_highest_bin_if_above_max(
            self, mocker) -> None:
        map_fn = mocker.patch('mqc.mbias.'
                              'get_sequence_context_to_array_index_table')
        map_fn.return_value = (
            SEQ_CONTEXT_TO_IDX_MAPPING, BINNED_SEQ_CONTEXT_TO_IDX_MAPPING)

        reads = [bstub_with(),
                 bstub_with(flen=1000),
                 bstub_with(phred=70)]

        motif_pileup = MotifPileupStub(seq_context='GGCGG',
                                       reads=reads)

        mbias_counter = MbiasCounter(CONFIG)
        mbias_counter.process(cast(MotifPileup, motif_pileup))

        base_event_class = ((SEQ_CONTEXT_TO_IDX_MAPPING['GGCGG'],)
                            + bstub_with().get_idx_tuple())

        flen_dim = mbias_counter.dim_names.index('flen')
        phred_dim = mbias_counter.dim_names.index('phred')
        max_flen_bin = mbias_counter.counter_array.shape[flen_dim] - 1
        max_phred_bin = mbias_counter.counter_array.shape[phred_dim] - 1
        pooled_flen_event = (base_event_class[0:flen_dim]
                             + (max_flen_bin,) + base_event_class[
                                                 flen_dim + 1:])
        pooled_phred_event = (base_event_class[0:phred_dim]
                              + (max_phred_bin,) + base_event_class[
                                                   phred_dim + 1:])
        exp_counts_idx = list(zip(base_event_class, pooled_flen_event,
                                  pooled_phred_event))
        assert (mbias_counter.counter_array[exp_counts_idx] == 1).all()
        mbias_counter.counter_array[exp_counts_idx] = 0
        assert (mbias_counter.counter_array == 0).all()

    def test_undefined_seq_contexts_are_discarded(
            self, mocker):
        """Undefined seq_contexts are given if the seq_context
        is unknown, i.e. not contained in MbiasCounter.seq_ctx_idx_dict

        This may happen because the seq_context
           - is shorter than the expected size (at ends of chromosomes)
           - contains Ns
           - in some cases, not all possible seq_context of a given length
             may be intended to be counted, and therefore will not be
             present in seq_ctx_idx_dict defined by the user
        """

        map_fn = mocker.patch('mqc.mbias.'
                              'get_sequence_context_to_array_index_table')
        # we are only interested in these three motifs
        map_fn.return_value = (SEQ_CONTEXT_TO_IDX_MAPPING,
                               BINNED_SEQ_CONTEXT_TO_IDX_MAPPING)

        reads = [bstub_with()]

        motif_pileups = [
            MotifPileupStub('GGCGG', reads),  # ok
            MotifPileupStub('NGCGG', reads),  # N
            MotifPileupStub('GCGG', reads),  # too short
            MotifPileupStub('TTCAA', reads),  # this context not of interest
        ]

        mbias_counter = MbiasCounter(CONFIG)
        for mp in motif_pileups:
            mbias_counter.process(cast(MotifPileup, mp))

        assert mbias_counter.counter_array[
                   (SEQ_CONTEXT_TO_IDX_MAPPING['GGCGG'],)
                   + bstub_with().get_idx_tuple()] == 1
        assert (mbias_counter.counter_array > 0).sum() == 1

    def test_motif_pileups_are_added_incrementally(
            self, mocker) -> None:
        map_fn = mocker.patch('mqc.mbias.'
                              'get_sequence_context_to_array_index_table')
        # we are only interested in these three motifs
        map_fn.return_value = (
            SEQ_CONTEXT_TO_IDX_MAPPING, BINNED_SEQ_CONTEXT_TO_IDX_MAPPING)

        reads = [bstub_with(flen=400,
                            expected_flen_bin=175,
                            phred=38,
                            expected_phred_bin=7),
                 bstub_with(flen=300,
                            expected_flen_bin=165,
                            phred=8,
                            expected_phred_bin=1)]

        event_classes_adv_idx = list(zip(*[
            (SEQ_CONTEXT_TO_IDX_MAPPING['AACAA'],) + curr_read.get_idx_tuple()
            for curr_read in reads]))

        mbias_counter = MbiasCounter(CONFIG)
        mbias_counter.counter_array[event_classes_adv_idx] += 1

        mbias_counter.process(cast(MotifPileup, MotifPileupStub('AACAA', reads)))
        mbias_counter.process(cast(MotifPileup, MotifPileupStub('AACAA', reads)))

        assert (mbias_counter.counter_array[event_classes_adv_idx] == 3).all()
        mbias_counter.counter_array[event_classes_adv_idx] = 0
        assert (mbias_counter.counter_array == 0).all()


@pytest.fixture()
def config_fit_normal_variate_plateau():
    config = defaultdict(defaultdict)
    config["trimming"]["min_plateau_perc"] = 0.8
    config["trimming"]["max_std_within_plateau"] = 0.1
    config["trimming"]["min_flen_considered_for_trimming"] = 50
    return config


#-
# Test M-bias plateau fitting
# ===========================

# prepare test shapes. Each shape is associated with 'correct' cutting sites
# either these are integers, e.g. if cuts at position 8 and 90 are expected:
# [8, 90]
# or these are ranges, e.g. if the left cut can be anywhere between position 5 and 15
# [(5,15), 90]

def add_noise_to_mbias_shape(arr, sd=0.002):
    np.random.seed(123)
    return arr + (np.random.randn(*arr.shape) * sd)

def sigmoid_curve(steepness, min_height, delta, longer_plateau_length, reverse=False):
    sigmoid_area = np.round((101 - longer_plateau_length) * 2)
    res = np.concatenate([
        1 / (1 + np.exp(-steepness*np.linspace(-5, 5, sigmoid_area))) * delta + min_height,
        np.ones(101 - sigmoid_area) * min_height + delta])
    if reverse:
        res = res[::-1]
    return res


def michaelis_menten_curve(vmax, km, steepness, flat_plateau_start=None, reverse=False):
    res = vmax * np.arange(101) * steepness / (km + np.arange(101) * steepness)
    if flat_plateau_start:
        res[flat_plateau_start:] = res[flat_plateau_start]
    if reverse:
        res = res[::-1]
    return res

#-
test_mbias_shapes = {
    'perfect_shape': (
        np.ones(101) * 0.75,
        [0, 101]),
    'left_gap_repair_nucleotides': (
        np.concatenate([np.zeros(9), np.ones(92) * 0.75]),
        [9, 101]),
    'right_gap_repair_nucleotides': (
        np.concatenate([np.ones(92) * 0.75, np.zeros(9)]),
        [0, 92]),
    'gaps_at_both_sites': (
        np.concatenate([np.zeros(5), np.ones(87) * 0.75, np.zeros(9)]),
        [5, 92]),

    'sigmoid_curve_1': [
        sigmoid_curve(steepness=3, min_height=0.65,
                      delta=0.2, longer_plateau_length=60),
        [(42, 46), 101]],
    'sigmoid_curve_2': [
        sigmoid_curve(steepness=3, min_height=0.65, delta=0.1,
                      longer_plateau_length=70, reverse=True),
        [0, (65, 75)]
    ],
    'slight_slope': [
        np.linspace(0.78, 0.8, 101),
        [(0, 5), (95, 101)],
    ],
    'medium_slope': [
        np.linspace(0.75, 0.8, 101),
        [(0, 20), (80, 101)],
    ],
    'bad_slope': [
        np.linspace(0.70, 0.8, 101),
        [0, 0]
    ],
    'very_bad_slope': [
        np.linspace(0.60, 0.8, 101),
        [0, 0]
    ],
    'intermittent_drops': [
        np.repeat([.7, 0, .7, 0, .7, 0, .7],
                  [20, 5, 30, 10, 20, 10, 5]),
        [0, 0]
    ],
    'michaelis-menten-curve1': [
        michaelis_menten_curve(vmax=0.75, km=20, steepness=10),
        [(10, 40), 101]
    ],
    'michaelis-menten-curve2': [
        michaelis_menten_curve(vmax=0.7, km=20, steepness=20),
        [(0, 20), 101]
    ],
    'michaelis-menten-curve3': [
        michaelis_menten_curve(vmax=0.7, km=20, steepness=2),
        [0, 0]
    ],
    'michaelis-menten-curve4': [
        michaelis_menten_curve(vmax=0.7, km=20, steepness=3,
                               flat_plateau_start=40, reverse=True),
        [0, (60,90)]
    ],
}
test_mbias_shapes_with_low_coverage_slope = ['michaelis-menten-curve1']
#-

# def visualize_mbias_test_shapes() -> None:
#     """Plot M-bias test function using interactive backend"""
#
#     axes: List[Axes]
#     fig, axes = plt.subplots(int(np.ceil(len(test_mbias_shapes) / 2)), 2,
#                              constrained_layout=True)
#     sd = 0.05
#     for i, (shape_name, (arr, exp_cutting_sites)) in enumerate(test_mbias_shapes.items(), 1):
#         axes_coord = (int(np.ceil(i / 2)) - 1, 1 - i % 2)
#         axes[axes_coord].plot(add_noise_to_mbias_shape(arr, sd))
#         axes[axes_coord].set_ylim([0,1])
#         axes[axes_coord].set_title(shape_name)
#     fig.show()



noise_levels = [0.005, 0.01]

# Each plateau fitter function is associated with a config dict and with
# tests which are expected to fail
plateau_fitter_functions = [fit_percentiles]

plateau_fitter_func_configs: Dict[Callable, Dict] = {
    fit_percentiles: {},
}

# Every plateau fitter function can have a list of shape names for which it is
# expected to fail
xfails_by_plateau_fitter_func_mapping: Dict[Callable, List[str]] = {
    fit_percentiles: [],
}

plateau_fitting_test_parametrizations = []
for plateau_fitter_func, test_mbias_shape_item, sd in product(
    plateau_fitter_functions, test_mbias_shapes.items(),
    noise_levels
):

    curr_marks: List[Any] = []  # List of pytest marks
    if (test_mbias_shape_item[0]
            in xfails_by_plateau_fitter_func_mapping[plateau_fitter_func]):
        curr_marks = [pytest.mark.xfail]

    plateau_fitting_test_parametrizations.append(
        pytest.param(plateau_fitter_func, *test_mbias_shape_item[1], sd,
                     plateau_fitter_func_configs[plateau_fitter_func],
                     marks=curr_marks, id=f'{test_mbias_shape_item[0]}_sd-{sd}')
    )


@pytest.mark.xfail(skip=True)
@pytest.mark.parametrize('plateau_caller_fn,beta_value_arr,'
                         'expected_cutting_sites_tuple,sd,config',
                         plateau_fitting_test_parametrizations)
def test_plateau_fitter(plateau_caller_fn, beta_value_arr,
                        expected_cutting_sites_tuple, sd, config):
    df_stub = pd.DataFrame(
        {'beta_value': add_noise_to_mbias_shape(beta_value_arr, sd=sd)})
    res_ser = plateau_caller_fn(df_stub, config)
    cut_end_ok_ser = pd.Series()
    for i, cut_end in enumerate(['left_cut_end', 'right_cut_end']):
        elem = expected_cutting_sites_tuple[i]
        if isinstance(elem, Sequence):
            assert len(elem) == 2, \
                    'Unexpected input for expected cutting sites'
            cut_end_ok_ser[cut_end] = elem[0] <= res_ser[cut_end] <= elem[1]
        else:
            assert isinstance(elem, int), \
                    'Unexpected input for expected cutting sites'
            cut_end_ok_ser[cut_end] = res_ser[cut_end] == elem
    message = f'Computed:\n{res_ser}\nExpected:{expected_cutting_sites_tuple}'
    assert cut_end_ok_ser.all(), message


class TestFitNormalVariatePlateau:
    def test_left_gap_nucleotides_are_recognized(self,
                                                 config_fit_normal_variate_plateau):
        df_stub = pd.DataFrame({'beta_value': np.ones(101) * 0.75})
        df_stub.loc[0:8, 'beta_value'] = 0
        res_ser = fit_normalvariate_plateau(df_stub,
                                            config_fit_normal_variate_plateau)
        assert res_ser.left_cut_end == 9
        assert res_ser.right_cut_end == 100

    def test_right_gap_nucleotides_are_recognized(self,
                                                  config_fit_normal_variate_plateau):
        df_stub = pd.DataFrame({'beta_value': np.ones(101) * 0.75})
        df_stub.loc[92:100, 'beta_value'] = 0
        res_ser = fit_normalvariate_plateau(df_stub,
                                            config_fit_normal_variate_plateau)
        assert res_ser.left_cut_end == 0
        assert res_ser.right_cut_end == 92

    def test_hill_curves_are_recognized(self,
                                        config_fit_normal_variate_plateau):
        df_stub = pd.DataFrame({'beta_value': np.ones(101) * 0.75})
        df_stub.loc[0:19, 'beta_value'] = np.linspace(0, 0.75, 20)
        res_ser = fit_normalvariate_plateau(df_stub,
                                            config_fit_normal_variate_plateau)
        print(res_ser)
        assert 10 < res_ser.left_cut_end < 20
        assert res_ser.right_cut_end == 100

    def test_plateaus_on_less_than_T_percent_of_the_read_length_are_discarded(
            self, config_fit_normal_variate_plateau):
        df_stub = pd.DataFrame({'beta_value': np.ones(80) * 0.75})
        df_stub.loc[0:19, 'beta_value'] = np.linspace(0, 0.3, 20)
        res_ser = fit_normalvariate_plateau(df_stub,
                                            config_fit_normal_variate_plateau)
        print(res_ser)
        assert res_ser.left_cut_end == 0.0
        assert res_ser.right_cut_end == 0.0

# TODO: adapt to extended Mbias stats counter
def test_compute_mbias_stats_df_converts_mbias_counts_to_mbias_stats_df():
    """
    from mqc.mbias import compute_mbias_stats_df
    mbias_counts_df_stub = pd.DataFrame([
        ['CG', 'w_bc', 1, 1, 'n_meth', 10],
        ['CG', 'w_bc', 1, 1, 'n_unmeth', 10],
        ['CG', 'w_bc', 1, 2, 'n_meth', 10],
        ['CG', 'w_bc', 1, 2, 'n_unmeth', 10],
        ['CG', 'w_bc', 1, 3, 'n_meth', 10],
        ['CG', 'w_bc', 1, 3, 'n_unmeth', 10],
        ['CG', 'w_bc', 2, 1, 'n_meth', 10],
        ['CG', 'w_bc', 2, 1, 'n_unmeth', 10],
        ['CG', 'w_bc', 2, 2, 'n_meth', 10],
        ['CG', 'w_bc', 2, 2, 'n_unmeth', 10],
        ['CG', 'w_bc', 2, 3, 'n_meth', 10],
        ['CG', 'w_bc', 2, 3, 'n_unmeth', 10],
        ['CG', 'w_bc', 3, 1, 'n_meth', 10],
        ['CG', 'w_bc', 3, 1, 'n_unmeth', 10],
        ['CG', 'w_bc', 3, 2, 'n_meth', 10],
        ['CG', 'w_bc', 3, 2, 'n_unmeth', 10],
        ['CG', 'w_bc', 3, 3, 'n_meth', 10],
        ['CG', 'w_bc', 3, 3, 'n_unmeth', 10],
    ], columns=['motif', 'bs_strand', 'flen', 'pos', 'meth_status',
                'counts']).set_index(
        ['motif', 'bs_strand', 'flen', 'pos', 'meth_status'])

    exp_mbias_stats_df = (pd.DataFrame([
        ['CG', 'w_bc', 1, 1, 0.5, 10, 10],
        ['CG', 'w_bc', 2, 1, 0.5, 10, 10],
        ['CG', 'w_bc', 2, 2, 0.5, 10, 10],
        ['CG', 'w_bc', 3, 1, 0.5, 10, 10],
        ['CG', 'w_bc', 3, 2, 0.5, 10, 10],
        ['CG', 'w_bc', 3, 3, 0.5, 10, 10], ],
        columns=['motif', 'bs_strand', 'flen', 'pos', 'beta_value', 'n_meth',
                 'n_unmeth'])
                          .set_index(['motif', 'bs_strand', 'flen', 'pos']))
    exp_mbias_stats_df

    comp_mbias_stats_df = compute_mbias_stats_df(mbias_counts_df_stub)

    assert exp_mbias_stats_df.equals(comp_mbias_stats_df)
    """
    pass




def test_mask_mbias_stats_df_sets_positions_in_trimming_zone_to_nan():

    mbias_stats_df_stub = (pd.DataFrame([['CG', 'w_bc', 1, 1, 0.5, 10, 10],
                                         ['CG', 'w_bc', 2, 1, 0.5, 10, 10],
                                         ['CG', 'w_bc', 2, 2, 0.5, 10, 10],
                                         ['CG', 'w_bc', 3, 1, 0.5, 10, 10],
                                         ['CG', 'w_bc', 3, 2, 0.5, 10, 10],
                                         ['CG', 'w_bc', 3, 3, 0.5, 10, 10], ],
                                        columns=['motif', 'bs_strand', 'flen', 'pos',
                                                 'beta_value', 'n_meth', 'n_unmeth'])
                           .set_index(['motif', 'bs_strand', 'flen', 'pos']))

    cutting_sites_df_stub = (pd.DataFrame([['w_bc', 1, 0, 1],
                                           ['w_bc', 2, 1, 2],
                                           ['w_bc', 3, 0, 2]],
                                          columns=['bs_strand', 'flen', 'start', 'end'])
                             .set_index(['bs_strand', 'flen'])
                             )
    cutting_sites_df_stub.columns.name = 'cutting_site'

    exp_masked_df = (pd.DataFrame([['CG', 'w_bc', 1, 1, 0.5, 10, 10],
                                   ['CG', 'w_bc', 2, 1, np.nan, np.nan, np.nan],
                                   ['CG', 'w_bc', 2, 2, 0.5, 10, 10],
                                   ['CG', 'w_bc', 3, 1, 0.5, 10, 10],
                                   ['CG', 'w_bc', 3, 2, 0.5, 10, 10],
                                   ['CG', 'w_bc', 3, 3, np.nan, np.nan, np.nan]],
                                  columns=['motif', 'bs_strand', 'flen', 'pos',
                                           'beta_value', 'n_meth', 'n_unmeth'])
                     .set_index(['motif', 'bs_strand', 'flen', 'pos']))

    computed_masked_df = mask_mbias_stats_df(mbias_stats_df_stub,
                                             cutting_sites_df_stub)

    assert exp_masked_df.equals(computed_masked_df)


@pytest.fixture()
def max_flen_cutting_df_to_array_conversion():
    return 100





# -----------------------------------------------------------------------------
# Acceptance tests
# -----------------------------------------------------------------------------

@pytest.fixture(scope='module')
def user_config_file():
    tmpdir = Path(tempfile.mkdtemp())
    user_config_file_path = tmpdir / "user_config_file.toml"
    user_config_file_path.write_text(dedent(f"""\
            [paths]
                mbias_counts = "{tmpdir}/mbias-counter"
            [data_properties]
                max_read_length_bp = 101
            [trimming]
                max_flen_considered_for_trimming = 60
                min_plateau_perc = 0.8
                max_std_within_plateau = 0.1
                min_flen_considered_for_trimming = 60
            [stats]
                max_flen = 300
                max_flen_with_single_flen_resolution = 110
                flen_bin_size = 30
                max_phred = 40
                phred_bin_size = 20
                seq_context_size = 5
                
            [plots]
                mbias_flens_to_display = [100, 150, 190]
            """))
    yield str(user_config_file_path)
    shutil.rmtree(tmpdir)


@pytest.fixture(scope='module')
def test_mbias_counter(user_config_file):
    with open(user_config_file, 'rt') as fin:
        config = toml.load(fin)

    with patch('mqc.mbias.' 'get_sequence_context_to_array_index_table') as map_fn:
        map_fn.return_value = (SEQ_CONTEXT_TO_IDX_MAPPING,
                               BINNED_SEQ_CONTEXT_TO_IDX_MAPPING)

        mbias_counter = MbiasCounter(config)
        mbias_counter.counter_array[:, :, :, :, :, 0] = 700
        mbias_counter.counter_array[:, :, :, :, :, 1] = 300
        mbias_counter.counter_array[:, :, :, :, 0:10, 0] = 0
        mbias_counter.counter_array[:, :, :, :, 0:10, 1] = 1000

        with open(config['paths']['mbias_counts'] + '.p', 'wb') as fout:
            pickle.dump(mbias_counter, fout)

        return mbias_counter
# to test the fixture:
# def test_mbias_counts(test_mbias_counter: MbiasCounter):
#     print(test_mbias_counter.get_dataframe()
#           .loc[idxs['GGCGG', 'c_bc',
#                16, 0, :, 'n_meth'], :])


# test_mbias_counter fixture must be called to place the counter
# at the appropriate file path
# noinspection PyUnusedLocal
@pytest.fixture(scope='module',
                params=['CG']
                # params=['CG', 'CG-CHG-CHH']
                )
def run_evaluate_mbias_then_return_config(
        request, user_config_file, test_mbias_counter):  # pylint: disable=unused-argument
    output_dir = tempfile.mkdtemp()

    subprocess.check_call(['mqc', 'evaluate_mbias',
                           '--config_file', user_config_file,
                           '--motifs', request.param,
                           '--sample_name', SAMPLE_NAME,
                           '--output_dir', output_dir])

    default_config_file = get_resource_abspath('config.default.toml')
    cli_params = {'motifs_str': request.param,
                  'sample_name': SAMPLE_NAME,
                  'sample_meta': None,
                  'output_dir': output_dir}
    config = assemble_config_vars(cli_params,
                                  default_config_file_path=default_config_file,
                                  user_config_file_path=user_config_file)

    yield config
    shutil.rmtree(output_dir)


@pytest.mark.acceptance_test
class TestMbiasEvaluate:
    def test_runs_through(self, run_evaluate_mbias_then_return_config):
        config = run_evaluate_mbias_then_return_config
        mbias_stats_df = pd.read_pickle(config['paths']['mbias_stats_p'])
        assert isinstance(mbias_stats_df, pd.DataFrame)

# # redundant with CLI-based acceptance test, but great for interactive debugging
# # in IDE
# @pytest.mark.acceptance_test
# class TestAnalyseMbiasCounts:
#     def test_runs_through(self, user_config_file):
#         tmpdir = tempfile.TemporaryDirectory()
#
#         config = assemble_config_vars(
#             command_line_args_dict={'sample_name': SAMPLE_NAME,
#                                     'sample_meta': 'population=hsc',
#                                     'output_dir': tmpdir.name,
#                                     'motifs_str': 'CG',
#                                     'no_cache': False},
#             default_config_file_path=get_resource_abspath(
#                 'config.default.toml'),
#             user_config_file_path=str(user_config_file)
#         )
#
#         compute_mbias_stats(config)

# M-bias plots
# ==============================================================================


#-
class MbiasStatsDfStub:

    def __init__(self):
        self.df = pd.DataFrame()

    @staticmethod
    def from_subsequent_strata_with_noise_overlay(
            idx_levels: List[List[Union[str, int]]],
            idx_level_names: List[str],
            beta_value_int_or_arr: Union[int, np.ndarray],
            coverage_int_or_arr: Union[int, np.ndarray],
            n_pos=100,
            seed=123,
            do_flen_adjustment=False):

        self = MbiasStatsDfStub()

        if isinstance(beta_value_int_or_arr, np.ndarray):
            n_pos = beta_value_int_or_arr.shape[0]
            beta_value_arr = beta_value_int_or_arr
        else:
            beta_value_arr = np.tile(beta_value_int_or_arr, n_pos)

        if not isinstance(coverage_int_or_arr, np.ndarray):
            coverage_arr = np.tile(coverage_int_or_arr, n_pos)
        else:
            coverage_arr = coverage_int_or_arr

        number_of_strata = np.prod([len(sublist) for sublist in idx_levels])

        coverage_column = np.tile(coverage_arr, number_of_strata)
        np.random.seed(seed)
        n_meth_column = np.random.binomial(
            n=coverage_column,
            p=np.tile(beta_value_arr, number_of_strata)
        )

        idx_levels_with_pos = idx_levels + [list(range(n_pos))]
        idx_level_names_with_pos = idx_level_names + ['pos']
        midx = pd.MultiIndex.from_product(idx_levels_with_pos,
                                          names=idx_level_names_with_pos)

        self.df = pd.DataFrame({'n_meth': n_meth_column,
                                'n_total': coverage_column,
                                'n_unmeth': coverage_column - n_meth_column},
                               index=midx)
        self.df['beta_value'] = self.df['n_meth'] / self.df['n_total']

        def adjust_for_flen(group_df, n_pos):
            flen = group_df.index.get_level_values('flen')[0]

            if flen >= n_pos:
                return group_df

            n_undef_pos = n_pos - flen
            middle_pos = n_pos // 2
            lower_bound = middle_pos - np.int(np.ceil(n_undef_pos / 2))
            upper_bound = middle_pos + np.int(np.floor(n_undef_pos / 2))
            high_bound = (n_pos - n_undef_pos)
            group_df.iloc[lower_bound:high_bound, :] = group_df.iloc[upper_bound:, :].values
            group_df.iloc[high_bound:, :] = np.nan

            return group_df

        if do_flen_adjustment:
            self.df = (self.df
                       .groupby(level=idx_level_names)
                       .apply(adjust_for_flen, n_pos)
             )

        return self


    @staticmethod
    def create_test_mbias_stats_df(
            strata_curve_tuples:
            List[Tuple[Union[Tuple, Callable],
                       Union[np.ndarray, int], Union[np.ndarray, int]]]) \
            -> pd.DataFrame:
        """ Create small mbias stats dataframe for testing

        Args:
            strata_curve_tuples: specifies index slices (tuple
            or callable) where beta values (np.ndarray or int)
            should be entered with coverage (np.ndarray or int)

        Returns:
            mbias stats dataframe with mbias stats entered into the
            appropriate slices. Undefined strata have n_meth = n_unmeth =
            10 and thus a beta value of 0.5
        """

        def get_categorical(ordered_str_labels):
            return pd.Categorical(ordered_str_labels,
                                  categories=ordered_str_labels,
                                  ordered=True)

        ordered_seq_contexts = ['CCCCC', 'GGCGG']
        flen_bin_labels = [50, 100, 150, 200, 250]
        phred_bin_labels = [20, 30, 40]
        max_read_length = 100

        dim_levels = [get_categorical(['CG', 'CHH']),
                      get_categorical(ordered_seq_contexts),
                      get_categorical(['c_bc', 'c_bc_rv', 'w_bc', 'w_bc_rv']),
                      flen_bin_labels,
                      phred_bin_labels,
                      range(1, max_read_length + 1),
                      get_categorical(['n_meth', 'n_unmeth'])]

        midx = pd.MultiIndex.from_product(dim_levels,
                                          names=["motif", "seq_context", "bs_strand", "flen",
                                                 "phred", "pos", "meth_status"])
        mbias_stats_df = pd.DataFrame(index=midx)

        mbias_stats_df['counts'] = 10
        mbias_stats_df = mbias_stats_df.unstack(level=-1)
        mbias_stats_df.columns = ['n_meth', 'n_unmeth']
        mbias_stats_df['beta_value'] = 0.5

        for idx, beta, coverage in strata_curve_tuples:
            n_meth = np.round((beta * coverage)).astype(int)  # type: ignore
            n_unmeth = coverage - n_meth
            if not isinstance(n_meth, np.ndarray):
                n_meth = np.tile(n_meth, max_read_length)
                n_unmeth = np.tile(n_unmeth, max_read_length)

            rep_times_float = mbias_stats_df.loc[idx, :].shape[0] / max_read_length
            assert rep_times_float.is_integer()
            mbias_stats_df.loc[idx, 'n_meth'] = np.tile(n_meth, int(rep_times_float))
            mbias_stats_df.loc[idx, 'n_unmeth'] = np.tile(n_unmeth, int(rep_times_float))
            mbias_stats_df = mbias_stats_df.assign(beta_value = compute_beta_values)

        return mbias_stats_df
#-


@pytest.fixture()
def mbias_plot_config_fix():
    return {
        "defaults": {
            'plot_params': {
                "panel_height_cm": 6,
                "panel_width_cm": 6,
                "theme": "paper",
                'y_axis': {'breaks': 5,
                           'limits': {
                               'beta_value': [(0, 1), (0.5, 1), 'auto'],
                               'n_meth': 'auto',
                               'n_unmeth': 'auto', },
                           'share': True,
                           'rotate_labels': False,
                           },
                'x_axis': {'breaks': [0, 30, 60, 90, 120, 150],
                           'limits': {
                               'phred': 45,
                           },
                           'share': True,
                           'rotate_labels': True},
                'plot': ['line'],
            },
            "post_agg_filters": {
                'flen': [60, 75, 90, 100, 120, 150, 190, 290, 490],
                'phred': [0, 5, 10, 15, 20, 25, 30, 35, 40],
            },
        },
        "group1": {
            "datasets": ["full", ["full", "trimmed"]],
            # necessary because currently only trimming for CG
            "pre_agg_filters": {
                'motif': ['CG']
            },
            'post_agg_filters': {
                'flen': [100, 150],
            },
            "aes_mappings": [
                {'x': 'pos', 'y': 'beta_value',
                 "row": "dataset", "column": ("mate", "bs_strand"), "color": "None"},
            ],
        },
        'group2': {
            "datasets": ["trimmed"],
            "pre_agg_filters": {
                'motif': ['CG', 'CHH'],
                'seq_context': ['CCCCC']
            },
            "plot_params": {
                'y_axis': {'limits': {'beta_value': [(0.2, 0.9)]}}
            },
            "aes_mappings": [
                {'x': 'pos', 'y': 'value',
                 "row": 'dataset', 'color': 'flen', "column": ("mate", "statistic")},
            ],
        },
    }

@pytest.fixture()
def expected_hierarchical_mbias_plot_config_dict(mbias_plot_config_fix):

    dataset_name_to_fp = {'trimmed': '/path/to/trimmed',
                          'full': '/path/to/full'}

    group1_config = dict(
        pre_agg_filters=mbias_plot_config_fix['group1']['pre_agg_filters'],
        post_agg_filters=tz.merge(mbias_plot_config_fix['defaults']['post_agg_filters'],
                                    mbias_plot_config_fix['group1']['post_agg_filters']),
        plot_params=MbiasPlotParams(
            **mbias_plot_config_fix['defaults']['plot_params']),
        aes_mapping=MbiasPlotMapping(
            **mbias_plot_config_fix['group1']['aes_mappings'][0]),
    )

    expected_config_dict = {
        'group1': [
            MbiasPlotConfig(
                datasets=subset_dict('full', dataset_name_to_fp),
                **group1_config
            ),
            MbiasPlotConfig(
                datasets=subset_dict(['full', 'trimmed'], dataset_name_to_fp),
                **group1_config,
            ),
        ],
        'group2': [
            MbiasPlotConfig(
                datasets=subset_dict('trimmed', dataset_name_to_fp),
                plot_params=MbiasPlotParams(
                    **tz.assoc_in(
                        mbias_plot_config_fix['defaults']['plot_params'],
                        ['y_axis', 'limits', 'beta_value'],
                        [(0.2, 0.9)])
                ),
                aes_mapping=MbiasPlotMapping(
                    **mbias_plot_config_fix['group2']['aes_mappings'][0]),
                post_agg_filters = mbias_plot_config_fix['defaults']['post_agg_filters'],
                pre_agg_filters = mbias_plot_config_fix['group2']['pre_agg_filters'],
            ),
        ]
    }

    return expected_config_dict


class TestMbiasAnalysisConfig:
    def test_construct_from_compact_dict_representation(
            self, mbias_plot_config_fix,
            expected_hierarchical_mbias_plot_config_dict):

        dataset_name_to_fp = {'trimmed': '/path/to/trimmed',
                              'full': '/path/to/full'}

        computed_hierarchical_config = (
            MbiasAnalysisConfig
                .from_compact_config_dict(mbias_plot_config_fix,
                                          dataset_name_to_fp)
                .grouped_mbias_plot_configs)

        assert (computed_hierarchical_config
                == expected_hierarchical_mbias_plot_config_dict)

    def test_get_report_config(
            self,
            mbias_plot_config_fix: Dict[str, Any],
            expected_hierarchical_mbias_plot_config_dict) -> None:

        def get_figure_config(mbias_plot_config: MbiasPlotConfig):
            title = 'Datasets: ' + ', '.join(mbias_plot_config.datasets) + '\n\n'
            title += json.dumps(mbias_plot_config.aes.get_plot_aes_dict(include_facetting=True),
                                sort_keys=True) + '\n'
            figure_config = {'path': (
                    mqc.filepaths.mbias_plots_trunk.name +
                    mbias_plot_config.get_str_repr_for_filename_construction()
                    + '.json'), 'title': title}
            return figure_config

        d = expected_hierarchical_mbias_plot_config_dict
        expected_plot_config = {
            'group1': {'figures': [
                get_figure_config(d['group1'][0]),
                get_figure_config(d['group1'][1]),
            ]},
            'group2': {'figures': [
                get_figure_config(d['group2'][0]),
            ]}
        }

        dataset_name_to_fp = {'trimmed': '/path/to/trimmed',
                              'full': '/path/to/full'}

        computed_report_config = (
            MbiasAnalysisConfig
                .from_compact_config_dict(mbias_plot_config_fix,
                                          dataset_name_to_fp)
                .get_report_config())

        assert expected_plot_config == computed_report_config

# TODO: force absolute dataset filepaths
class TestGetMbiasPlotConfigs:
    def test_produces_list_of_mbiasplotconfigs_from_defaults_and_group_definitions(
            self, mbias_plot_config_fix):

        dataset_name_to_fp = {'trimmed': '/path/to/trimmed',
                              'full': '/path/to/full'}

        computed_table_of_required_plots = get_plot_configs(mbias_plot_config_fix,
                                                            dataset_name_to_fp)

        group1_config = dict(
            pre_agg_filters = mbias_plot_config_fix['group1']['pre_agg_filters'],
            post_agg_filters = tz.merge(mbias_plot_config_fix['defaults']['post_agg_filters'],
                                        mbias_plot_config_fix['group1']['post_agg_filters']),
            plot_params=MbiasPlotParams(
                **mbias_plot_config_fix['defaults']['plot_params']),
            aes_mapping=MbiasPlotMapping(
                **mbias_plot_config_fix['group1']['aes_mappings'][0]),
        )

        expected_table_of_required_plots = [
            MbiasPlotConfig(
                datasets=subset_dict('full', dataset_name_to_fp),
                **group1_config
            ),
            MbiasPlotConfig(
                datasets=subset_dict(['full', 'trimmed'], dataset_name_to_fp),
                **group1_config,
            ),
            MbiasPlotConfig(
                datasets=subset_dict('trimmed', dataset_name_to_fp),
                plot_params=MbiasPlotParams(
                    **tz.assoc_in(
                        mbias_plot_config_fix['defaults']['plot_params'],
                        ['y_axis', 'limits', 'beta_value'],
                        [(0.2, 0.9)])
                    ),
                aes_mapping=MbiasPlotMapping(
                    **mbias_plot_config_fix['group2']['aes_mappings'][0]),
                post_agg_filters = mbias_plot_config_fix['defaults']['post_agg_filters'],
                pre_agg_filters = mbias_plot_config_fix['group2']['pre_agg_filters'],
            ),
            ]

        message = (str(expected_table_of_required_plots[0])
                   + '\n' + str(computed_table_of_required_plots[0]))

        assert (expected_table_of_required_plots
                == computed_table_of_required_plots), message


    @pytest.mark.parametrize('missing_key',
                             [('group1', 'datasets'),
                              ('group1', 'aes_mappings', 0, 'x'),
                              ('group2', 'aes_mappings')])
    def test_raises_when_plot_group_config_is_incomplete(self, missing_key,
                                                         mbias_plot_config_fix):
        dict_to_modify = tz.get_in(missing_key[:-1],
                                   mbias_plot_config_fix, no_default=True)
        dict_to_modify.pop(missing_key[-1])

        with pytest.raises(ValueError):
            get_plot_configs(mbias_plot_config_fix,
                             dataset_name_to_fp={})



# def test_map_dataset_specs_to_filepaths():
#     config = {'default_plot_params': {'a': 1,
#                                       'b': 2 },
#               'plot_group1': {'datasets': ['mbias_stats', 'mbias_stats_masked'],
#                               'other_stuff': [1,2,3] },
#               'plot_group2': {'datasets': ['mbias_stats_classic'],
#                               'other_stuff': [1, 2, 3]},
#               }
#     map_dataset_specs_to_filepaths(config)
#
#     assert config['plot_group1']['datasets'] == [
#         mqc.filepaths.mbias_stats,
#         mqc.filepaths.mbias_stats_masked
#     ]
#
#     assert config['plot_group2']['datasets'] == [
#         mqc.filepaths.mbias_stats_classic
#     ]


def test_create_aggregated_tables(mocker, tmpdir,
                                  mbias_plot_config_fix: Dict[str, Dict]) -> None:

    dataset1 = pd.DataFrame({'a': [1, 2, 3]})
    dataset1_fp = str(tmpdir / 'dataset1.p')
    dataset1.to_pickle(dataset1_fp)
    dataset2 = pd.DataFrame({'b': [2, 2, 3]})
    dataset2_fp = str(tmpdir / 'dataset2.p')
    dataset2.to_pickle(dataset2_fp)
    dataset_name_to_fp = dict(dataset1=dataset1_fp,
                              dataset2=dataset2_fp)

    aes_dict = dict(x='x_col', y='y_col', color='color_col', row='row_col')

    # To determine the correct aggregation tasks, MbiasPlotConfig functionality,
    # in particular determining the aggregation variables, is necessary -
    # so don't use a Mock here
    plot_configs = [
        MbiasPlotConfig(datasets=subset_dict('dataset1', dataset_name_to_fp),
                        aes_mapping=MbiasPlotMapping(**aes_dict, column='col1'),  # type: ignore
                        plot_params=mbias_plot_config_fix['defaults']['plot_params'],
                        pre_agg_filters={'motifs': ['CG']}
                        ),
        MbiasPlotConfig(datasets=subset_dict('dataset1', dataset_name_to_fp),
                        aes_mapping=MbiasPlotMapping(**aes_dict, column='col2'),   # type: ignore
                        plot_params=mbias_plot_config_fix['defaults']['plot_params'],
                        pre_agg_filters={'motifs': ['CHG']}
                        ),
        MbiasPlotConfig(datasets=subset_dict('dataset2', dataset_name_to_fp),
                        aes_mapping=MbiasPlotMapping(**aes_dict, column='col3'),   # type: ignore
                        plot_params=mbias_plot_config_fix['defaults']['plot_params'],
                        pre_agg_filters={'motifs': ['CG']}
                        ),
        MbiasPlotConfig(datasets=subset_dict('dataset2', dataset_name_to_fp),
                        aes_mapping=MbiasPlotMapping(**aes_dict, column='col4'),   # type: ignore
                        plot_params=mbias_plot_config_fix['defaults']['plot_params'],
                        pre_agg_filters={'motifs': ['CHG']}
                        ),
    ]

    mocker.patch(
        'mqc.mbias.AggregatedMbiasStats.precompute_aggregated_stats')
    aggregated_mbias_stats = AggregatedMbiasStats()

    create_aggregated_tables(mbias_plot_configs=plot_configs,
                             aggregated_mbias_stats=aggregated_mbias_stats)

    shared_variables = {'x_col', 'color_col', 'row_col'}
    expected_calls = [
        dict(params=MbiasStatAggregationParams(
            dataset_fp=dataset1_fp,
            variables=shared_variables | {'col1'},
            pre_agg_filters={'motifs': ['CG']}),
            mbias_stats_df=dataset1),
        dict(params=MbiasStatAggregationParams(
            dataset_fp=tmpdir / 'dataset1.p',
            variables=shared_variables | {'col2'},
            pre_agg_filters={'motifs': ['CHG']}),
            mbias_stats_df=dataset1),
        dict(params=MbiasStatAggregationParams(
            dataset_fp=tmpdir / 'dataset2.p',
            variables=shared_variables | {'col3'},
            pre_agg_filters={'motifs': ['CG']}),
            mbias_stats_df=dataset2),
        dict(params=MbiasStatAggregationParams(
            dataset_fp=tmpdir / 'dataset2.p',
            variables=shared_variables | {'col4'},
            pre_agg_filters={'motifs': ['CHG']}),
            mbias_stats_df=dataset2),
    ]

    # precompute..stats is patched
    # noinspection PyUnresolvedReferences
    call_args_list = aggregated_mbias_stats.precompute_aggregated_stats.call_args_list  # type: ignore

    assert len(call_args_list) == 4

    for unused_args, kwargs in call_args_list:
        for curr_expected_call in expected_calls:
            if (curr_expected_call['mbias_stats_df'].equals(kwargs['mbias_stats_df'])
                and curr_expected_call['params'] == kwargs['params']):
                break
        else:
            # False positive
            # noinspection PyUnboundLocalVariable
            raise ValueError('Did not find call: ',
                             curr_expected_call, call_args_list)


@pytest.fixture()
def aggregate_table_mbias_stats_df_fix():
    return pd.read_csv(StringIO(dedent("""\
            motif    bs_strand    flen    pos    n_meth  n_unmeth
            CG             c_bc         100     1      1       1
            CG             c_bc         100     2      2       2
            CG             c_bc         200     1      3       3
            CG             c_bc         200     2      4       4
            CG             w_bc         100     1      5       5
            CG             w_bc         100     2      6       6
            CG             w_bc         200     1      7       7
            CG             w_bc         200     2      8       8
            CH             c_bc         100     1      1       1
            CH             c_bc         100     2      2       2
            CH             c_bc         200     1      3       3
            CH             c_bc         200     2      4       4
            CH             w_bc         100     1      5       5
            CH             w_bc         100     2      6       6
            CH             w_bc         200     1      7       7
            CH             w_bc         200     2      8       8
            """)),
                       sep=r'\s+', header=0,
                       index_col=['motif', 'bs_strand', 'flen', 'pos'])


class TestAggregateTable:

    # noinspection PyShadowingNames
    def test_applies_filters_and_sum_aggregation(
            self, aggregate_table_mbias_stats_df_fix):

        agg_data = aggregate_table(
            mbias_stats_df=aggregate_table_mbias_stats_df_fix,
            unordered_variables={'pos', 'bs_strand', 'motif'},
            pre_agg_filters={'motif': ['CG']})

        # Aggregation always triggers a recomputation of beta values, therefore
        # we don't need beta values in the input data, and do still expect them
        # in the aggregation result
        expected_agg_data = pd.read_csv(StringIO(dedent(
            '''\
            motif    bs_strand    pos    n_meth n_unmeth   beta_value
            CG        c_bc         1      4      4          0.5
            CG        c_bc         2      6     6         0.5
            CG        w_bc         1      12     12         0.5
            CG        w_bc         2      14     14         0.5
            ''')), sep=r'\s+', header=0, index_col=['motif', 'bs_strand', 'pos'])

        pd.testing.assert_frame_equal(agg_data, expected_agg_data)

    def test_raises_when_filters_or_variables_do_not_match_index_level_names(
            self, aggregate_table_mbias_stats_df_fix):

        with pytest.raises(ValueError):
            _ = aggregate_table(
                    mbias_stats_df=aggregate_table_mbias_stats_df_fix,
                    unordered_variables={'pos', 'bs_strand', 'motif'},
                    pre_agg_filters={'MOTIF': ['CG']})

        with pytest.raises(ValueError):
            _ = aggregate_table(
                    mbias_stats_df=aggregate_table_mbias_stats_df_fix,
                    unordered_variables={'POS', 'bs_strand', 'motif'},
                    pre_agg_filters={'motif': ['CG']})


class TestMbiasPlotMapping:
    @pytest.mark.xfail(skip=True)
    def test_get_facetting_cmd(self):
        raise NotImplementedError

    @pytest.mark.parametrize(
        'kwargs,correct_set',
        [
            (dict(x='pos', y='beta_value', column=('mate', 'bs_strand')),
             {'pos', 'mate', 'bs_strand'}),
            (dict(x='pos', y='value', column=('statistic', 'bs_strand')),
             {'pos', 'bs_strand'}),
            (dict(x='pos', y='beta_value'),
             {'pos'}),
            # don't include dataset variable!
            (dict(x='pos', y='beta_value', row='dataset'),
             {'pos'}),
        ])
    def test_get_all_variables_unordered(self, kwargs, correct_set):
        computed_set = (MbiasPlotMapping(**kwargs)
                        .get_all_agg_variables_unordered())
        assert computed_set == correct_set
    #
    # def test1(self):
    #     raise NotImplemented




class TestAggregatedMbiasStats:

    def test_fp_creation(self) -> None:
        aggregated_mbias_stats = AggregatedMbiasStats()
        params = MbiasStatAggregationParams(
            dataset_fp='/path/to/dataset.p',
            variables={'varA', 'varB'},
            pre_agg_filters={'flen': [1, 2, 3]}
        )
        fp = aggregated_mbias_stats._create_agg_dataset_fp(params)
        assert fp == (mqc.filepaths.aggregated_mbias_stats_dir
        / (f'{params.dataset_fp.replace("/", "-")}'
           f'_{"-".join(params.variables)}'
           f'_{json.dumps(params.pre_agg_filters).replace(" ", "-")}.p'))

    def test_stats_are_only_computed_if_there_is_no_recent_copy_present(
            self, mocker, tmpdir) -> None:

        with TmpChdir(tmpdir):
            aggregated_mbias_stats = AggregatedMbiasStats()
            mbias_stats_df = pd.DataFrame({'a': [1,2,3]})

            aggregate_mbias_stats_mock: MagicMock = mocker.patch(
                'mqc.mbias.aggregate_table', return_value=mbias_stats_df)
            params = MbiasStatAggregationParams(
                dataset_fp='path/to/dataset.p',
                variables={'varA', 'varB'},
                pre_agg_filters={'flen': [1, 2, 3]}
            )
            Path(params.dataset_fp).parent.mkdir(parents=True, exist_ok=False)
            Path(params.dataset_fp).touch()
            for _ in range(2):
                aggregated_mbias_stats.precompute_aggregated_stats(
                    params=params,
                    mbias_stats_df=mbias_stats_df
                )
            assert aggregate_mbias_stats_mock.called_once_with(
                call(mbias_stats_df=mbias_stats_df,
                     variables=params.variables,
                     pre_agg_filters=params.pre_agg_filters)
            )



@pytest.mark.xfail(skip=True)
class TestMbiasPlotConfig:
    def test_raises_if_filters_are_not_defined_correctly(self):
        raise NotImplementedError
    def test_raises_if_datasets_are_not_defined_correctly(self):
        raise NotImplementedError
    def test_hash(self):
        raise NotImplementedError

@pytest.fixture()
def test_mbias_plot_config_file() -> str:
    """Return path to test M-bias plot config file

    Note that this does not append the specification of the config dict
    to be used (e.g. '*::default_config)
    """
    return str(Path(__file__).parent / 'test_files' / 'test_mbias_plot_config.py')


@pytest.fixture()
def full_mbias_stats_fp():
    tmpdir = tempfile.mkdtemp()
    fp = os.path.join(tmpdir, 'full_mbias_stats_df.p')
    df = MbiasStatsDfStub.create_test_mbias_stats_df(
        [
            (NamedIndexSlice(motif='CHH'), 0.05, 20),
            (NamedIndexSlice(motif='CG'), 0.8, 20),
            (NamedIndexSlice(motif='CG', flen=slice(1, 110)),
             np.concatenate([np.tile(0.01, 9), np.tile(0.8, 91)]), 30),
        ]
    )
    df.to_pickle(fp)
    yield fp
    shutil.rmtree(tmpdir)

@pytest.fixture()
def trimmed_mbias_stats_fp(tmpdir):
    fp = tmpdir / 'full_mbias_stats_df.p'
    df = MbiasStatsDfStub.create_test_mbias_stats_df(
        [
            (NamedIndexSlice(motif='CHH'), 0.05, 20),
            (NamedIndexSlice(motif='CG'), 0.8, 20),
            (NamedIndexSlice(motif='CG', flen=slice(1, 110)),
             np.concatenate([np.tile(np.nan, 9), np.tile(0.8, 91)]), 30),
        ]
    )
    df.to_pickle(fp)
    yield fp


# Once it is possible to give a dict name to --mbias_plot_config,
# I can add a default dict suitable for testing to the default
# mbias_plot_config module. Until then, it is not trivial
# to test running mbias_plots without the mbias_plot_config option
# because the current default config file has post_agg_filters with too many
# levels, which will raise a ValueError
@pytest.fixture(params=[test_mbias_plot_config_file() + '::default_config'])  # type: ignore
def run_mbias_stats_plots(
        request, full_mbias_stats_fp: str, trimmed_mbias_stats_fp: str) -> None:


    output_dir = tempfile.mkdtemp()

    command_list = ['mqc', 'mbias_plots',
                    '--output_dir', output_dir,
                    '--datasets',
                    f'full={full_mbias_stats_fp},trimmed={trimmed_mbias_stats_fp}',
                    ]

    if request.param is not None:
        command_list += ['--mbias_plot_config', request.param]

    subprocess.run(command_list, check=True)

    yield output_dir

    shutil.rmtree(output_dir)

@pytest.mark.xfail(skip=True)
def test_create_mbias_stats_plots():
    raise NotImplementedError

# noinspection PyMethodMayBeStatic
@pytest.mark.acceptance_test
class TestMbiasPlots:
    def test_runs_through(self, run_mbias_stats_plots,
                          make_interactive):
        """Trigger fixture to run mbias_plots, then check results

        The run_mbias_stats_plots fixture runs mbias_plots and returns
        the base output dir used for the test run.

        In interactive mode, the output_dir is opened in firefox, to
        allow viewing the figures
        """
        if make_interactive:
            subprocess.run(['firefox', run_mbias_stats_plots])
            while True:
                answer=input('Everything ok?: y/n\n')
                if answer in ['y', 'n']:
                    if answer == 'n':
                        raise AssertionError
                    break
            assert True
        else:
            with TmpChdir(run_mbias_stats_plots):
                plots = list(mqc.filepaths.mbias_plots_trunk.parent.glob('*.json'))
            assert len(plots) >= 1

    # TODO: improve; currently this fails because of several reasons, ...
    # at least:
    # - the default config file expects phred-threshold datasets
    # - the post_agg filters have levels which don't exist in the test data
    def test_fails_with_inappropriate_config_file(self, tmpdir,
                                                  full_mbias_stats_fp,
                                                  trimmed_mbias_stats_fp):

        # Leave --mbias_plot_config to default. The default config file
        # has too many elements in the post_agg_filters
        # This must raise a ValueError
        command_list = ['mqc', 'mbias_plots',
                        '--output_dir', tmpdir,
                        '--datasets',
                        f'full={full_mbias_stats_fp},trimmed={trimmed_mbias_stats_fp}',
                        ]
        with pytest.raises(subprocess.CalledProcessError):
            subprocess.run(command_list, check=True)


class TestCuttingSitesPlots:
    def test_cutting_sites_line_plot(self, tmpdir):

        cutting_sites_df = pd.DataFrame(
            {'start': 9, 'end': 92},
            index=pd.MultiIndex.from_product(
                [pd.Categorical('c_bc c_bc_rv'.split(), ordered=True), list(range(50, 150))],
                names='bs_strand flen'.split())
        )
        cutting_sites_df.loc[('c_bc', slice(70, 90)), 'start'] = 15
        cutting_sites_df.columns.name = 'cutting_site'
        cutting_sites = CuttingSites(
            cutting_sites_df=cutting_sites_df,
            max_read_length=101,
        )

        path = str(tmpdir.join('test.html'))
        cutting_sites.plot(path)
        # subprocess.run(['firefox', path])

        path = str(tmpdir.join('test.svg'))
        cutting_sites.plot(path)
        # subprocess.run(['firefox', path])

class TestCuttingSites:
    @pytest.mark.xfail(skip=True)
    def test_invariant_detects_wrong_cutting_sites_df_format_during_init(self):
        with pytest.raises(AssertionError):
            _ = CuttingSites(
                pd.DataFrame({'start': [1, 2, 3]}),
                max_read_length=100,
            )

        with pytest.raises(AssertionError):
            _ = CuttingSites(
                cutting_sites_df=pd.DataFrame(
                    index=pd.MultiIndex.from_product(
                        ['c_bc w_bc'.split(),
                         range(100, 150)], names='bs_strand flen'.split())
                ).assign(
                    start=9,
                    end=92.0
                ),
                max_read_length=100
            )
    @pytest.mark.xfail(skip=True)
    def test_fails_if_mbias_stats_df_does_not_match_max_read_length(self):
        raise NotImplementedError

    def test_from_relative_to_frag_end_cutting_sites(self) -> None:

        cut_site_spec = dict(
            w_bc = (0, 1),
            c_bc = (0, 2),
            w_bc_rv = (5, 0),
            c_bc_rv = (1, 0),
        )
        cutting_sites = CuttingSites.from_rel_to_frag_end_cutting_sites(
            cut_site_spec=cut_site_spec, max_read_length=2)

        # max flen is the first flen where the read cannot be affected any
        # more in any strand, based on cut_site_spec and max_read_length
        expected_df = pd.DataFrame([
            ['c_bc', 0, 0, 0],
            ['c_bc', 1, 0, 0],
            ['c_bc', 2, 0, 0],
            ['c_bc', 3, 0, 1],
            ['c_bc', 4, 0, 2],
            ['c_bc_rv', 0, 0, 0],
            ['c_bc_rv', 1, 0, 0],
            ['c_bc_rv', 2, 1, 2],
            ['c_bc_rv', 3, 1, 2],
            ['c_bc_rv', 4, 1, 2],
            ['w_bc', 0, 0, 0],
            ['w_bc', 1, 0, 0],
            ['w_bc', 2, 0, 1],
            ['w_bc', 3, 0, 2],
            ['w_bc', 4, 0, 2],
            ['w_bc_rv', 0, 0, 0],
            ['w_bc_rv', 1, 0, 0],
            ['w_bc_rv', 2, 0, 0],
            ['w_bc_rv', 3, 0, 0],
            ['w_bc_rv', 4, 0, 0],
        ], columns=['bs_strand', 'flen', 'start', 'end'])
        expected_df['bs_strand'] = pd.Categorical(expected_df['bs_strand'], ordered=True)
        expected_df = expected_df.set_index(['bs_strand', 'flen'])
        expected_df.columns.name = 'cutting_site'

        pd.testing.assert_frame_equal(cutting_sites.df, expected_df)


    def test_as_array(self):

        cutting_sites_df = pd.DataFrame(index=pd.MultiIndex.from_product(
            [pd.Categorical('c_bc c_bc_rv w_bc w_bc_rv'.split()),
             list(range(20, 31)) + [35, 40, 60] ],
            names='bs_strand flen'.split()))
        cutting_sites_df = cutting_sites_df.assign(start=9, end=92)
        cutting_sites_df.loc[idxs[:, 25:59], 'start'] = 5
        cutting_sites_df.loc[idxs[:, 25:59], 'end'] = 95
        cutting_sites_df.loc[idxs[:, 60], 'start'] = 0
        cutting_sites_df.loc[idxs[:, 60], 'end'] = 101

        cutting_sites_df.columns.name = 'cutting_site'
        cutting_sites = CuttingSites(
            cutting_sites_df=cutting_sites_df, max_read_length=100)

        expected_arr = np.zeros((4, 61, 2), dtype='i8')
        expected_arr[0:4, 20:25, 0] = 9
        expected_arr[0:4, 20:25, 1] = 92
        expected_arr[0:4, 20:25, 0] = 9
        expected_arr[0:4, 20:25, 1] = 92
        expected_arr[0:4, 25:41, 0] = 5
        expected_arr[0:4, 25:41, 1] = 95
        expected_arr[0:4, 25:41, 0] = 5
        expected_arr[0:4, 25:41, 1] = 95
        expected_arr[0:4, 41:61, 0] = 0
        expected_arr[0:4, 41:61, 1] = 101
        expected_arr[0:4, 41:61, 0] = 0
        expected_arr[0:4, 41:61, 1] = 101

        computed_arr = cutting_sites.as_array()

        assert (computed_arr == expected_arr).all()



class TestBinomPvalueBasedCuttingSiteDetermination:
    @pytest.mark.parametrize('coverage_int_or_arr', [
        500,
        2000
    ])
    @pytest.mark.parametrize('seed', list(range(1)))
    @pytest.mark.parametrize('shape_name', [
        'perfect_shape',
        'left_gap_repair_nucleotides',
        'right_gap_repair_nucleotides',
        'gaps_at_both_sites',
        # 'michaelis-menten-curve1',
        'michaelis-menten-curve4',
        # 'slight_slope',
        # 'medium_slope',
        # 'bad_slope',
        # 'very_bad_slope',
    ])
    @pytest.mark.parametrize('plateau_height_factor', [1])
    @pytest.mark.parametrize('allow_slope', [True])
    def test_identify_plateaus_in_various_full_read_length_shapes(
            self, coverage_int_or_arr, seed, shape_name,
            plateau_height_factor, allow_slope, make_interactive):

        flens = list(range(120, 220))
        self._cutting_sites_detection_test(shape_name, coverage_int_or_arr,
                                           flens, plateau_height_factor, seed,
                                           allow_slope, make_interactive,
                                           adjust_cutting_site_calls_from_small_flens=False)


    @pytest.mark.parametrize('coverage_int_or_arr', [
        500,
        2000
    ])
    @pytest.mark.parametrize('seed', list(range(1)))
    @pytest.mark.parametrize('shape_name', [
        'perfect_shape',
        'left_gap_repair_nucleotides',
        'right_gap_repair_nucleotides',
        'gaps_at_both_sites',
    ])
    @pytest.mark.parametrize('plateau_height_factor', [1])
    @pytest.mark.parametrize('allow_slope', [True])
    def test_finds_gap_repair_nucleotides_in_flen_smaller_than_read_length(
            self, coverage_int_or_arr, seed, shape_name,
            plateau_height_factor, make_interactive, allow_slope):

        flens = list(range(50, 190))
        self._cutting_sites_detection_test(shape_name, coverage_int_or_arr,
                                           flens, plateau_height_factor, seed,
                                           allow_slope, make_interactive,
                                           adjust_cutting_site_calls_from_small_flens=True)

    def _cutting_sites_detection_test(self, shape_name, coverage_int_or_arr,
                                      flens, plateau_height_factor, seed,
                                      allow_slope, make_interactive,
                                      adjust_cutting_site_calls_from_small_flens):
        beta_value_arr, (correct_start, correct_end) = (
            test_mbias_shapes[shape_name])
        cutting_sites_df, mbias_stats_stub = self._create_mbias_stub_based_cutting_sites(
            allow_slope, beta_value_arr, coverage_int_or_arr, flens,
            plateau_height_factor, seed, adjust_cutting_site_calls_from_small_flens)

        pd.testing.assert_index_equal(cutting_sites_df.index,
                                      mbias_stats_stub.df.unstack('pos').index)

        # low_coverage_strata_i_idx = list(range(10)) + list(range(-10, 0))
        # self._remove_low_coverage_strata(cutting_sites,
        #                                  low_coverage_strata_i_idx)

        if make_interactive:
            raise NotImplementedError
            # plot_stem = f'{shape_name}_cov-{coverage_int_or_arr}_seed-{seed}_fact-{plateau_height_factor}_high-flen'
            # self._plot_cutting_sites(mbias_stats_df=mbias_stats_stub,
            #                          cutting_sites_df=cutting_sites_df,
            #                          stem=plot_stem)

        if adjust_cutting_site_calls_from_small_flens:
            adjustment_summand = beta_value_arr.shape[0] - np.array(flens)
            adjustment_summand[adjustment_summand < 0] = 0
            adjustment_summand = pd.Series(adjustment_summand,
                                           index=flens)
            cutting_sites_df['end'] = cutting_sites_df['end'].add(adjustment_summand,
                                                                  axis=0,
                                                                  level='flen')
            cutting_sites_df.loc[cutting_sites_df['end'] > beta_value_arr.shape[0], 'end'] = beta_value_arr.shape[0]

        is_low_coverage_situation = plateau_height_factor < 0.8 or coverage_int_or_arr <= 500
        self._assert_cutting_site_correctness(cutting_sites_df, correct_start,
                                              correct_end,
                                              is_low_coverage_situation,
                                              shape_name)

    def _create_mbias_stub_based_cutting_sites(self, allow_slope,
                                               beta_value_arr,
                                               coverage_int_or_arr, flens,
                                               plateau_height_factor, seed,
                                               adjust_for_small_flens):
        mbias_stats_stub = self._get_mbias_stats_stub(beta_value_arr,
                                                      coverage_int_or_arr,
                                                      flens,
                                                      plateau_height_factor,
                                                      seed,
                                                      adjust_small_flens=adjust_for_small_flens)
        cutting_sites = BinomPvalueBasedCuttingSiteDetermination(
            mbias_stats_df=mbias_stats_stub.df,
            max_read_length=100,
            allow_slope=allow_slope, plateau_flen=170).compute_cutting_site_df()

        return cutting_sites, mbias_stats_stub

    @staticmethod
    def _assert_cutting_site_correctness(cutting_sites_df, correct_start,
                                         correct_end, is_low_coverage_situation,
                                         shape_name) -> None:
        assert cutting_sites_df_has_correct_format(cutting_sites_df)
        if (is_low_coverage_situation
                and shape_name in test_mbias_shapes_with_low_coverage_slope
                and (cutting_sites_df == 0).all().all()):
            return

        if (correct_start, correct_end) == (0, 0):
            assert (all(cutting_sites_df['start'] == 0) and all(
                cutting_sites_df['end']) == 0), f'{cutting_sites_df}\n0,0'

        message = f'{cutting_sites_df["start"]}\n{correct_start}'
        if isinstance(correct_start, int):
            if not is_low_coverage_situation:
                assert all(
                    cutting_sites_df['start'] == correct_start), str(
                    cutting_sites_df['start'].loc[
                    cutting_sites_df['start'] != correct_start, :])
            else:
                assert all(cutting_sites_df[
                               'start'] >= correct_start), message
        else:
            if not is_low_coverage_situation:
                is_correct = cutting_sites_df['start'].between(correct_start[0], correct_start[1])
                assert all(is_correct), cutting_sites_df.loc[~is_correct, :]
            else:
                # Pycharm does not understand that we are comparing Series
                # noinspection PyTypeChecker
                assert all(cutting_sites_df['start'] >= correct_start[
                    0]), message
        message = f'{cutting_sites_df["end"]}\n{correct_end}'
        if isinstance(correct_end, int):
            if not is_low_coverage_situation:
                pd.testing.assert_series_equal(cutting_sites_df['end'],
                                               pd.Series(correct_end,
                                                         index=cutting_sites_df['end'].index),
                                               check_dtype=False, check_names=False)
            else:
                assert all(
                    cutting_sites_df['end'] <= correct_end), message
        else:
            if not is_low_coverage_situation:
                is_correct = cutting_sites_df['end'].between(correct_end[0], correct_end[1])
                assert all(is_correct), cutting_sites_df.loc[~is_correct, :]
            else:
                # Pycharm does not understand that we are comparing Series
                # noinspection PyTypeChecker
                assert all(
                    cutting_sites_df['end'] <= correct_end[1]), message

    @staticmethod
    def _get_mbias_stats_stub(beta_value_arr, coverage_int_or_arr, flens,
                              plateau_height_factor, seed, adjust_small_flens):
        idx_levels = [pd.Categorical(['c_bc', 'w_bc']), flens]
        idx_level_names = ['bs_strand', 'flen']
        mbias_stats_df = (
        MbiasStatsDfStub.from_subsequent_strata_with_noise_overlay(
            idx_levels=idx_levels, idx_level_names=idx_level_names,
            beta_value_int_or_arr=beta_value_arr * plateau_height_factor,
            coverage_int_or_arr=coverage_int_or_arr, seed=seed,
            do_flen_adjustment=adjust_small_flens))
        return mbias_stats_df



    # @staticmethod
    # def _plot_cutting_sites(mbias_stats_df, cutting_sites_df, stem):
    #     index_names = list(mbias_stats_df.df.index.names)
    #     index_names.remove('pos')
    #     target_fp = f'/home/stephen/temp/test_plots/{stem}.html'
    #     Path(target_fp).parent.mkdir(parents=True, exist_ok=True)
    #     plot_df = mbias_stats_df.df.reset_index()
    #     plot_df['group_label'] = plot_df[index_names].astype(str).apply(lambda ser: ser.str.cat(sep='_'), axis=1)
    #
    #     mbias_lines_layer = (alt.Chart(plot_df)
    #                          .mark_line(opacity=0.05, color='black')
    #                          .encode(x='pos', y='beta_value', detail='group_label:N')
    #                          )
    #
    #     cutting_sites_layer = (alt.Chart(cutting_sites_df.stack()
    #                                      .to_frame('cutting_site')
    #                                      .reset_index())
    #                            .mark_rule(opacity=0.05)
    #                            .encode(x='cutting_site')
    #                            )
    #
    #     full_beta_value_chart = (mbias_lines_layer + cutting_sites_layer).interactive()
    #     # full_chart = cutting_sites_layer
    #
    #     p_value_plot_df = (cutting_sites_df
    #                        .predecessor_p_values
    #                        .stack()
    #                        .to_frame('p_value')
    #                        .reset_index())
    #     p_value_plot_df['group_label'] = (p_value_plot_df[index_names]
    #                                       .astype(str)
    #                                       .apply(lambda ser: ser.str.cat(sep='_'), axis=1))
    #     # p_value_plot_df['p_value'] = np.log10(p_value_plot_df['p_value'] + 10**-16)
    #     p_value_plot_df['p_value'] = (p_value_plot_df['p_value'] + 10**-16)
    #
    #     p_value_chart = (alt.Chart(p_value_plot_df)
    #         .mark_line(opacity=0.05)
    #         .encode(
    #         # y='p_value',
    #         alt.Y('p_value',
    #               scale=alt.Scale(type='log', base=10),
    #               axis=alt.Axis(orient='right')
    #               ),
    #         x='pos',
    #         detail='group_label'
    #     )
    #     )
    #
    #     full_p_value_chart = (p_value_chart + cutting_sites_layer).interactive()
    #
    #     full_chart = alt.vconcat(full_beta_value_chart, full_p_value_chart)
    #
    #
    #     (full_chart
    #      .configure(background='white')
    #      # .interactive()
    #      # .save(target_fp, webdriver='firefox')
    #      .save(target_fp)
    #      )
