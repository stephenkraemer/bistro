import os
import shutil
import tempfile
from collections import namedtuple, defaultdict
from copy import deepcopy
from itertools import product
from pathlib import Path
from textwrap import dedent
from unittest.mock import MagicMock
import os.path as op
import pytoml
from itertools import product

import pickle
import pandas as pd
import time

idxs = pd.IndexSlice
import numpy as np
from numpy.testing import assert_array_equal
import pytest
import subprocess

from mqc.config import assemble_config_vars
# from mqc.mbias import MbiasCounter
from mqc.mbias import MbiasCounter, create_mbias_stats_plots
from mqc.mbias import \
    get_sequence_context_to_array_index_table, map_seq_ctx_to_motif
from mqc.mbias import FixedRelativeCuttingSites
from mqc.mbias import convert_cutting_sites_df_to_array
from mqc.mbias import fit_normalvariate_plateau
from mqc.mbias import AdjustedCuttingSites
from mqc.mbias import mask_mbias_stats_df
from mqc.mbias import map_seq_ctx_to_motif
from mqc.mbias import convert_cutting_sites_df_to_array
from mqc.mbias import analyze_mbias_counts
from mqc.utils import get_resource_abspath

import matplotlib

matplotlib.use('Agg')  # import before pyplot import!
import matplotlib.pyplot as plt
import seaborn as sns

import mqc.flag_and_index_values as mfl
import plotnine as gg

# TODO: I currently test that any non-zero qc_fail_flag leads to discard from M-bias stats counting. When I update the behavior so that phred score fails are kept in the stats, the tests here also need to be updated accordingly
b_inds = mfl.bsseq_strand_indices
b_na_ind = mfl.bsseq_strand_na_index
m_flags = mfl.methylation_status_flags
qc_flags = mfl.qc_fail_flags

CONFIG = defaultdict(dict)
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
    def __init__(self, strand_idx, flen, pos, mflag,
                 fail, phred, expected_flen_bin, expected_phred_bin):
        self.qc_fail_flag = fail
        self.bsseq_strand_ind = strand_idx
        self.meth_status_flag = mflag
        self.alignment = AlignmentStub(template_length=flen)
        self.pos_in_read = pos
        self.baseq_at_pos = phred
        self.expected_phred_bin = expected_phred_bin
        self.expected_flen_bin = expected_flen_bin

    def get_idx_tuple(self):
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
                       fail=0,
                       phred=30,
                       expected_phred_bin=6)
    event_class.update(kwargs)
    return bstub(**event_class)


class MotifPileupStub:
    def __init__(self, seq_context, reads):
        self.idx_pos = IndexPositionStub(seq_context)
        self.reads = reads


class TestSeqContextToBinMapping():
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
        seq_ctx_to_idx, binned_motif_to_idx = get_sequence_context_to_array_index_table(
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
        seq_ctx_to_idx, binned_motif_to_idx = (
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
        seq_ctx_to_idx, binned_motif_to_idx = (
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

    flen_levels = pd.IntervalIndex.from_breaks(
        [0, 1, 2, 3, 4, 5, 6, 8, 10, 11],
        closed="left")
    phred_levels = pd.IntervalIndex.from_breaks(
        [0, 5, 10, 11], closed='left')
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
    exp_df["counts"] = np.uint64(0)

    mbias_counter = MbiasCounter(config)

    # random middle values
    exp_df.loc[idxs["WWCWW", "w_bc", 7, 2, 5, "n_meth":"n_meth"], "counts"] = 1
    mbias_counter.counter_array[2, 2, 6, 0, 4, 0] = 1
    exp_df.loc[idxs["CCCCC", "c_bc", 9, 6, 3, "n_unmeth":"n_unmeth"],
               "counts"] = 2
    mbias_counter.counter_array[0, 0, 7, 1, 2, 1] = 2
    # max values
    exp_df.loc[idxs["CCCCC", "c_bc_rv", 10, 10, 10, "n_unmeth":"n_unmeth"],
               "counts"] = 1
    mbias_counter.counter_array[0, 1, -1, -1, -1, 1] = 1
    # min values
    exp_df.loc[idxs["GGCGG", "w_bc_rv", 0, 0, 1, "n_meth":"n_meth"],
               "counts"] = 1
    mbias_counter.counter_array[1, 3, 0, 0, 0, 0] = 1

    computed_df = mbias_counter.get_dataframe()
    pd.testing.assert_frame_equal(exp_df, computed_df)


class TestMbiasCounterMotifPileupProcessing:
    def test_phred_flen_seqcontext_in_allowed_range_are_binned(
            self, mocker):
        reads = [bstub(strand_idx=b_inds.w_bc,
                       flen=100,
                       expected_flen_bin=100,
                       pos=10,
                       mflag=m_flags.is_methylated,
                       fail=0,
                       phred=20,
                       expected_phred_bin=4),
                 bstub(strand_idx=b_inds.c_bc,
                       flen=180,
                       expected_flen_bin=153,
                       pos=10,
                       mflag=m_flags.is_methylated,
                       fail=0,
                       phred=30,
                       expected_phred_bin=6), ]

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
        mbias_counter.process(motif_pileup1)
        mbias_counter.process(motif_pileup2)

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

    def test_unusable_reads_are_discarded(self, mocker):
        """ Reads are discarded if

            - they have a qc_fail_flag
            - they have no methylation calling status: NA, SNP or Ref
        """

        map_fn = mocker.patch('mqc.mbias.'
                              'get_sequence_context_to_array_index_table')
        map_fn.return_value = (
            SEQ_CONTEXT_TO_IDX_MAPPING, BINNED_SEQ_CONTEXT_TO_IDX_MAPPING)

        reads = [bstub_with(),
                 bstub_with(mflag=m_flags.is_ref),
                 bstub_with(mflag=m_flags.is_snp),
                 bstub_with(mflag=m_flags.is_na),
                 bstub_with(fail=qc_flags.overlap_fail),
                 bstub_with(fail=qc_flags.sam_flag_fail),
                 bstub_with(fail=qc_flags.phred_score_fail),
                 bstub_with(fail=qc_flags.mapq_fail),
                 ]

        motif_pileup = MotifPileupStub(seq_context='CCCCC',
                                       reads=reads)

        mbias_counter = MbiasCounter(CONFIG)

        mbias_counter.process(motif_pileup)

        base_event_class = ((SEQ_CONTEXT_TO_IDX_MAPPING['CCCCC'],)
                            + bstub_with().get_idx_tuple())
        assert mbias_counter.counter_array[base_event_class] == 1
        mbias_counter.counter_array[base_event_class] = 0
        assert (mbias_counter.counter_array == 0).all()

    def test_flen_and_phred_are_pooled_in_highest_bin_if_above_max(
            self, mocker):
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
        mbias_counter.process(motif_pileup)

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
            mbias_counter.process(mp)

        assert mbias_counter.counter_array[
                   (SEQ_CONTEXT_TO_IDX_MAPPING['GGCGG'],)
                   + bstub_with().get_idx_tuple()] == 1
        assert (mbias_counter.counter_array > 0).sum() == 1

    def test_motif_pileups_are_added_incrementally(
            self, mocker):
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

        motif_pileup = MotifPileupStub('AACAA', reads)

        event_classes_adv_idx = list(zip(*[
            (SEQ_CONTEXT_TO_IDX_MAPPING['AACAA'],) + curr_read.get_idx_tuple()
            for curr_read in reads]))

        mbias_counter = MbiasCounter(CONFIG)
        mbias_counter.counter_array[event_classes_adv_idx] += 1

        mbias_counter.process(MotifPileupStub('AACAA', reads))
        mbias_counter.process(MotifPileupStub('AACAA', reads))

        assert (mbias_counter.counter_array[event_classes_adv_idx] == 3).all()
        mbias_counter.counter_array[event_classes_adv_idx] = 0
        assert (mbias_counter.counter_array == 0).all()


def test_fixed_cutting_sites_are_computed_correctly():
    config = defaultdict(defaultdict)
    config['trimming']['relative_to_fragment_ends_dict'] = dict(
        w_bc=[10, 10],
        c_bc=[0, 10],
        w_bc_rv=[10, 0],
        c_bc_rv=[10, 10],
    )
    config['data_properties']['max_read_length_bp'] = 101
    config['trimming']['max_flen_considered_for_trimming'] = 500
    fixed_cutting_sites = FixedRelativeCuttingSites(config)
    cutting_arr = fixed_cutting_sites.get_array()
    cutting_df = fixed_cutting_sites.get_df()

    expected_result = (pd.DataFrame([['w_bc', 70, 'left_cut_end', 10],
                                     ['c_bc_rv', 90, 'right_cut_end', 79],
                                     ['w_bc', 110, 'right_cut_end', 99],
                                     ['w_bc', 111, 'right_cut_end', 100],
                                     ['w_bc', 115, 'right_cut_end', 100],
                                     ['c_bc_rv', 300, 'left_cut_end', 10]],
                                    columns=['bs_strand', 'flen', 'cut_end',
                                             'cut_pos'])
                       .set_index(['bs_strand', 'flen', 'cut_end'])
                       .astype(np.uint32)
                       )

    computed_result = cutting_df.loc[expected_result.index, :]

    pd.testing.assert_frame_equal(expected_result, computed_result)


def test_cutting_site_df_to_cutting_site_array_conversion():
    df = pd.DataFrame([
        ['w_bc', 100, 'left_cut_end', 10],
        ['c_bc', 110, 'right_cut_end', 20],
    ], columns=['bs_strand', 'flen', 'cut_end', 'cut_pos']).set_index(
        ['bs_strand', 'flen', 'cut_end'])
    max_flen = 120
    arr = convert_cutting_sites_df_to_array(df, max_flen)
    assert arr[b_inds.w_bc, 100, 0] == 10
    assert arr[b_inds.c_bc, 110, 1] == 20


@pytest.fixture()
def config_fit_normal_variate_plateau():
    config = defaultdict(defaultdict)
    config["trimming"]["min_plateau_perc"] = 0.8
    config["trimming"]["max_std_within_plateau"] = 0.1
    config["trimming"]["min_flen_considered_for_trimming"] = 50
    return config


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


class TestAdjustedCuttingSites:
    def test_compute_df_provides_df_in_standard_cutting_sites_df_format(self,
                                                                        mocker):
        mbias_df_stub = pd.DataFrame([
            ['CG', 'w_bc', 100, 1, 0.8],
            ['CG', 'w_bc', 100, 2, 0.8],
            ['CG', 'c_bc', 101, 1, 0.7],
            ['CG', 'c_bc', 101, 2, 0.7],
            ['CG', 'w_bc', 101, 1, 0.8],
            ['CG', 'w_bc', 101, 2, 0.8],
        ], columns=['motif', 'bs_strand', 'flen', 'pos',
                    'beta_value']).set_index(
            ['motif', 'bs_strand', 'flen', 'pos'])

        fit_mock = MagicMock(side_effect=[
            pd.Series([10, 11], index=['left_cut_end', 'right_cut_end']),
            pd.Series([10, 11], index=['left_cut_end', 'right_cut_end']),
            pd.Series([10, 11], index=['left_cut_end', 'right_cut_end']),
        ])
        mocker.patch('mqc.mbias.fit_normalvariate_plateau', fit_mock)

        computed_cutting_sites_ser = AdjustedCuttingSites._compute_df(
            mbias_df_stub, config={})

        # Note that this is sorted
        expected_cutting_sites_ser = pd.DataFrame([
            ['c_bc', 101, 'left_cut_end', 10],
            ['c_bc', 101, 'right_cut_end', 11],
            ['w_bc', 100, 'left_cut_end', 10],
            ['w_bc', 100, 'right_cut_end', 11],
            ['w_bc', 101, 'left_cut_end', 10],
            ['w_bc', 101, 'right_cut_end', 11],
        ], columns=['bs_strand', 'flen', 'cut_end', 'cut_pos']).set_index(
            ['bs_strand', 'flen', 'cut_end'])

        assert computed_cutting_sites_ser.equals(expected_cutting_sites_ser)


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
    mbias_stats_df_stub = (pd.DataFrame([
        ['CG', 'w_bc', 1, 1, 0.5, 10, 10],
        ['CG', 'w_bc', 2, 1, 0.5, 10, 10],
        ['CG', 'w_bc', 2, 2, 0.5, 10, 10],
        ['CG', 'w_bc', 3, 1, 0.5, 10, 10],
        ['CG', 'w_bc', 3, 2, 0.5, 10, 10],
        ['CG', 'w_bc', 3, 3, 0.5, 10, 10], ],
        columns=['motif', 'bs_strand', 'flen', 'pos', 'beta_value', 'n_meth',
                 'n_unmeth'])
                           .set_index(['motif', 'bs_strand', 'flen', 'pos']))

    cutting_sites_df_stub = (pd.DataFrame([
        ['w_bc', 1, 'left_cut_end', 1],
        ['w_bc', 1, 'right_cut_end', 1],
        ['w_bc', 2, 'left_cut_end', 2],
        ['w_bc', 2, 'right_cut_end', 2],
        ['w_bc', 3, 'left_cut_end', 2],
        ['w_bc', 3, 'right_cut_end', 3],
    ], columns=['bs_strand', 'flen', 'cut_end', 'cut_pos'])
                             .set_index(['bs_strand', 'flen', 'cut_end']))

    exp_masked_df = (pd.DataFrame([
        ['CG', 'w_bc', 1, 1, 0.5, 10, 10],
        ['CG', 'w_bc', 2, 1, np.nan, np.nan, np.nan],
        ['CG', 'w_bc', 2, 2, 0.5, 10, 10],
        ['CG', 'w_bc', 3, 1, np.nan, np.nan, np.nan],
        ['CG', 'w_bc', 3, 2, 0.5, 10, 10],
        ['CG', 'w_bc', 3, 3, 0.5, 10, 10], ],
        columns=['motif', 'bs_strand', 'flen', 'pos', 'beta_value', 'n_meth',
                 'n_unmeth'])
                     .set_index(['motif', 'bs_strand', 'flen', 'pos']))

    computed_masked_df = mask_mbias_stats_df(mbias_stats_df_stub,
                                             cutting_sites_df_stub)

    assert exp_masked_df.equals(computed_masked_df)


@pytest.fixture()
def max_flen_cutting_df_to_array_conversion():
    return 100


@pytest.fixture()
def cutting_sites_array_computed_from_df(
        max_flen_cutting_df_to_array_conversion):
    cutting_sites_df_stub = (pd.DataFrame([
        ['w_bc', 1, 'left_cut_end', 1],
        ['w_bc', 1, 'right_cut_end', 1],
        ['w_bc', 2, 'left_cut_end', 2],
        ['w_bc', 2, 'right_cut_end', 2],
        ['w_bc', 3, 'left_cut_end', 2],
        ['w_bc', 3, 'right_cut_end', 3],
    ], columns=['bs_strand', 'flen', 'cut_end', 'cut_pos'])
                             .set_index(['bs_strand', 'flen', 'cut_end']))

    cutting_sites_array = convert_cutting_sites_df_to_array(
        cutting_sites_df_stub, max_flen_cutting_df_to_array_conversion)

    return cutting_sites_array


class TestConvertCuttingSitesDfToArray:
    def test_cutting_sites_array_covers_all_possible_strata(
            self, cutting_sites_array_computed_from_df,
            max_flen_cutting_df_to_array_conversion
    ):
        assert cutting_sites_array_computed_from_df.shape == (
            4, max_flen_cutting_df_to_array_conversion + 1, 2)

    @pytest.mark.parametrize('stratum,exp_value',
                             (((2, 1, 0), 1,),
                              ((2, 1, 1), 1,),
                              ((2, 2, 0), 2,),
                              ((2, 2, 1), 2),
                              ((2, 3, 0), 2),
                              ((2, 3, 1), 3)))
    def test_cutting_sites_array_is_filled_based_on_df_where_data_are_available(
            self, stratum, exp_value, cutting_sites_array_computed_from_df):
        assert cutting_sites_array_computed_from_df[stratum] == exp_value

    def test_cutting_sites_array_is_0_0_where_no_data_are_available(
            self, cutting_sites_array_computed_from_df):
        idx = list(zip((2, 1, 0),
                       (2, 1, 1),
                       (2, 2, 0),
                       (2, 2, 1),
                       (2, 3, 0),
                       (2, 3, 1)))
        cutting_sites_array_computed_from_df[idx] = 0
        assert (cutting_sites_array_computed_from_df == 0).all()


# -----------------------------------------------------------------------------
# Acceptance tests
# -----------------------------------------------------------------------------

@pytest.fixture(scope='module')
def user_config_file():
    test_file = (Path(__file__).parent
                 / "test_files/test-evaluate-stats_mbias-counter")
    tmpdir = Path(tempfile.mkdtemp())
    user_config_file_path = tmpdir / "user_config_file.toml"
    user_config_file_path.write_text(dedent(f"""\
            [paths]
                mbias_counts = "{test_file}"
            [data_properties]
                max_read_length_bp = 101
            [trimming]
                max_flen_considered_for_trimming = 30
                min_plateau_perc = 0.8
                max_std_within_plateau = 0.1
                min_flen_considered_for_trimming = 10
            [stats]
                max_flen = 30
                max_flen_with_single_flen_resolution = 5
                flen_bin_size = 10
                max_phred = 40
                phred_bin_size = 5
                seq_context_size = 5
                
            [plots]
                mbias_flens_to_display = [100, 150, 190]
            """))
    yield str(user_config_file_path)
    shutil.rmtree(tmpdir)


@pytest.fixture(scope='module',
                params=['CG', 'CG-CHG-CHH'])
def run_evaluate_mbias_then_return_config(request, user_config_file):
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

# redundant with CLI-based acceptance test, but great for interactive debugging
# in IDE
@pytest.mark.acceptance_test
class TestAnalyseMbiasCounts:
    def test_runs_through(self, user_config_file):
        tmpdir = tempfile.TemporaryDirectory()

        config = assemble_config_vars(
            command_line_args_dict={'sample_name': SAMPLE_NAME,
                                    'sample_meta': 'population=hsc',
                                    'output_dir': tmpdir.name,
                                    'motifs_str': 'CG',
                                    'no_cache': False},
            default_config_file_path=get_resource_abspath(
                'config.default.toml'),
            user_config_file_path=str(user_config_file)
        )

        analyze_mbias_counts(config)

