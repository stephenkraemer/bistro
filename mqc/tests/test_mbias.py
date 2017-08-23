import shutil
import tempfile
from collections import namedtuple, defaultdict
from copy import deepcopy
from itertools import product
from textwrap import dedent
from unittest.mock import MagicMock
import os.path as op
import pytoml


import pandas as pd
import numpy as np
import pytest
import subprocess

from mqc.config import assemble_config_vars
from mqc.mbias import MbiasCounter
from mqc.mbias import FixedRelativeCuttingSites
from mqc.mbias import convert_cutting_sites_df_to_array
from mqc.mbias import fit_normalvariate_plateau
from mqc.mbias import AdjustedCuttingSites
from mqc.mbias import mask_mbias_stats_df
from mqc.mbias import convert_cutting_sites_df_to_array
from mqc.mbias import analyze_mbias_counts

import mqc.flag_and_index_values as mfl

#TODO: I currently test that any non-zero qc_fail_flag leads to discard from M-bias stats counting. When I update the behavior so that phred score fails are kept in the stats, the tests here also need to be updated accordingly
b_inds = mfl.bsseq_strand_indices
b_na_ind = mfl.bsseq_strand_na_index
m_flags = mfl.methylation_status_flags

MAX_FLEN = 500
CONFIG = defaultdict(dict)
CONFIG['data_properties']['max_read_length_bp'] = 101
CONFIG['trimming']['max_flen_considered_for_trimming'] = MAX_FLEN
CONFIG['run']['motifs'] = ['cg', 'chg']

AlignmentStub = namedtuple('AlignmentStub',
                           'template_length')

IndexPositionStub = namedtuple('IndexPositionStub',
                               'motif')

class BSSeqPileupReadStub:
    def __init__(self, strand_idx, tlen, pos, mflag, fail):
        self.qc_fail_flag = fail
        self.bsseq_strand_ind = strand_idx
        self.meth_status_flag = mflag
        self.alignment = AlignmentStub(template_length = tlen)
        self.pos_in_read = pos
    # TODO: change to method
    def __iter__(self):
        if self.meth_status_flag == m_flags.is_methylated:
            meth_status_index = 0
        elif self.meth_status_flag == m_flags.is_unmethylated:
            meth_status_index = 1
        else:
            raise ValueError("This is a failed read, "
                             "can't give event_class indices")

        return iter((self.bsseq_strand_ind,
                     self.alignment.template_length,
                     self.pos_in_read,
                     meth_status_index))

bstub = BSSeqPileupReadStub

class MotifPileupStub:
    def __init__(self, motif, reads):
        self.idx_pos = IndexPositionStub(motif=motif)
        self.reads = reads

class TestMbiasCounterMotifPileupProcessing:

    def test_updates_counter_array(self):
        reads = [
            bstub(strand_idx=b_inds.w_bc, tlen=100, pos=50,
                  mflag=m_flags.is_methylated, fail=0),
            bstub(strand_idx=b_inds.c_bc_rv, tlen=120, pos=50,
                  mflag=m_flags.is_unmethylated, fail=0)
        ]

        mbias_counter = MbiasCounter(CONFIG)

        # add reads to CG and CHG motif stratum
        for motif_ind, curr_read in product([0,1], reads):
                mbias_counter.counter_array[(motif_ind,) + tuple(curr_read)] = 2

        motif_pileup_cg = MotifPileupStub(motif = 'cg',
                                          reads=reads)
        motif_pileup_chg = MotifPileupStub(motif = 'chg',
                                           reads=reads)
        mbias_counter.process(motif_pileup_cg)
        mbias_counter.process(motif_pileup_chg)

        # assert counts for CG
        assert mbias_counter.counter_array[(0,) + tuple(reads[0])] == 3
        assert mbias_counter.counter_array[(0,) + tuple(reads[1])] == 3

        # assert counts for CHG
        assert mbias_counter.counter_array[(1,) + tuple(reads[0])] == 3
        assert mbias_counter.counter_array[(1,) + tuple(reads[1])] == 3

    def test_discards_reads_w_bad_qcflag_or_na_bsstrand_or_mcall_fail_ref_snp(self):
        """ Reads are discarded if

            - they have a qc_fail_flag
            - their bsseq strand could not be determined
            - they have methylation calling status: NA, SNP or Ref
        """
        reads = [
            bstub(strand_idx=b_inds.w_bc, tlen=100, pos=50, mflag=m_flags.is_methylated, fail=1),  # qc fail
            bstub(strand_idx=b_na_ind,    tlen=100, pos=50, mflag=m_flags.is_methylated, fail=0), # bsseq strand not identifiable
            bstub(strand_idx=b_inds.w_bc, tlen=100, pos=50, mflag=m_flags.is_na, fail=0), # na meth
            bstub(strand_idx=b_inds.w_bc, tlen=100, pos=50, mflag=m_flags.is_snp, fail=0),  # snp meth
            bstub(strand_idx=b_inds.w_bc, tlen=100, pos=50, mflag=m_flags.is_ref, fail=0),  # ref meth
            bstub(strand_idx=b_inds.c_bc, tlen=200, pos=90, mflag=m_flags.is_unmethylated, fail=0),  # this one should count
        ]
        motif_pileup = MotifPileupStub(motif = 'cg', reads=reads)
        mbias_counter = MbiasCounter(CONFIG)
        mbias_counter.process(motif_pileup)
        assert mbias_counter.counter_array[(0,) + tuple(reads[-1])] == 1
        mbias_counter.counter_array[(0,) + tuple(reads[-1])] = 0
        assert (mbias_counter.counter_array == 0).all()

    def test_threshold_exceeding_flens_are_added_to_max_flen_bin(self):
        too_long_flen_read_properties = dict(
            strand_idx=b_inds.w_bc_rv, tlen=MAX_FLEN + 200, pos=20,
            mflag=m_flags.is_methylated, fail=0)
        reads = [bstub(**too_long_flen_read_properties)]
        motif_pileup = MotifPileupStub(motif = 'chg', reads=reads)
        mbias_counter = MbiasCounter(CONFIG)
        mbias_counter.process(motif_pileup)

        capped_flen_read_properties = deepcopy(too_long_flen_read_properties)
        capped_flen_read_properties['tlen'] = MAX_FLEN
        capped_read_idx = tuple(bstub(**capped_flen_read_properties))
        assert mbias_counter.counter_array[(1,) + capped_read_idx] == 1

def test_fixed_cutting_sites_are_computed_correctly():
    config = defaultdict(defaultdict)
    config['trimming']['relative_to_fragment_ends_dict'] = dict(
        w_bc = [10, 10],
        c_bc = [0, 10],
        w_bc_rv = [10, 0],
        c_bc_rv = [10, 10],
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
                                   columns=['bs_strand', 'flen', 'cut_end', 'cut_pos'])
                       .set_index(['bs_strand', 'flen', 'cut_end'])
                       )

    computed_result = cutting_df.loc[expected_result.index, :]

    msg = f"""
    Expected result
    ---------------
    {expected_result}
    
    Computed result
    ---------------
    {computed_result}
    
    """

    assert computed_result.equals(expected_result), msg

def test_cutting_site_df_to_cutting_site_array_conversion():
    df = pd.DataFrame([
        ['w_bc', 100, 'left_cut_end',  10],
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
    def test_left_gap_nucleotides_are_recognized(self, config_fit_normal_variate_plateau):
        df_stub = pd.DataFrame({'beta_value': np.ones(101) * 0.75})
        df_stub.loc[0:8, 'beta_value'] = 0
        res_ser = fit_normalvariate_plateau(df_stub, config_fit_normal_variate_plateau)
        assert res_ser.left_cut_end == 9
        assert res_ser.right_cut_end == 100

    def test_right_gap_nucleotides_are_recognized(self, config_fit_normal_variate_plateau):
        df_stub = pd.DataFrame({'beta_value': np.ones(101) * 0.75})
        df_stub.loc[92:100, 'beta_value'] = 0
        res_ser = fit_normalvariate_plateau(df_stub, config_fit_normal_variate_plateau)
        assert res_ser.left_cut_end == 0
        assert res_ser.right_cut_end == 92

    def test_hill_curves_are_recognized(self, config_fit_normal_variate_plateau):
        df_stub = pd.DataFrame({'beta_value': np.ones(101) * 0.75})
        df_stub.loc[0:19, 'beta_value'] = np.linspace(0, 0.75, 20)
        res_ser = fit_normalvariate_plateau(df_stub, config_fit_normal_variate_plateau)
        print(res_ser)
        assert 10 < res_ser.left_cut_end < 20
        assert res_ser.right_cut_end == 100

    def test_plateaus_on_less_than_T_percent_of_the_read_length_are_discarded(self, config_fit_normal_variate_plateau):
        df_stub = pd.DataFrame({'beta_value': np.ones(80) * 0.75})
        df_stub.loc[0:19, 'beta_value'] = np.linspace(0, 0.3, 20)
        res_ser = fit_normalvariate_plateau(df_stub, config_fit_normal_variate_plateau)
        print(res_ser)
        assert res_ser.left_cut_end == 0.0
        assert res_ser.right_cut_end == 0.0

class TestAdjustedCuttingSites:
    def test_compute_df_provides_df_in_standard_cutting_sites_df_format(self, mocker):

        mbias_df_stub = pd.DataFrame([
            ['CG', 'w_bc', 100, 1, 0.8],
            ['CG', 'w_bc', 100, 2, 0.8],
            ['CG', 'c_bc', 101, 1, 0.7],
            ['CG', 'c_bc', 101, 2, 0.7],
            ['CG', 'w_bc', 101, 1, 0.8],
            ['CG', 'w_bc', 101, 2, 0.8],
        ], columns=['motif', 'bs_strand', 'flen', 'pos', 'beta_value']).set_index(
            ['motif', 'bs_strand', 'flen', 'pos'])

        fit_mock = MagicMock(side_effect = [
            pd.Series([10,11], index=['left_cut_end', 'right_cut_end']),
            pd.Series([10,11], index=['left_cut_end', 'right_cut_end']),
            pd.Series([10,11], index=['left_cut_end', 'right_cut_end']),
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

from mqc.mbias import compute_mbias_stats_df
def test_compute_mbias_stats_df_converts_mbias_counts_to_mbias_stats_df():
    mbias_counts_df_stub = pd.DataFrame([
        ['CG', 'w_bc', 1, 1, 'n_meth',   10],
        ['CG', 'w_bc', 1, 1, 'n_unmeth', 10],
        ['CG', 'w_bc', 1, 2, 'n_meth',   10],
        ['CG', 'w_bc', 1, 2, 'n_unmeth', 10],
        ['CG', 'w_bc', 1, 3, 'n_meth',   10],
        ['CG', 'w_bc', 1, 3, 'n_unmeth', 10],
        ['CG', 'w_bc', 2, 1, 'n_meth',   10],
        ['CG', 'w_bc', 2, 1, 'n_unmeth', 10],
        ['CG', 'w_bc', 2, 2, 'n_meth',   10],
        ['CG', 'w_bc', 2, 2, 'n_unmeth', 10],
        ['CG', 'w_bc', 2, 3, 'n_meth',   10],
        ['CG', 'w_bc', 2, 3, 'n_unmeth', 10],
        ['CG', 'w_bc', 3, 1, 'n_meth',   10],
        ['CG', 'w_bc', 3, 1, 'n_unmeth', 10],
        ['CG', 'w_bc', 3, 2, 'n_meth',   10],
        ['CG', 'w_bc', 3, 2, 'n_unmeth', 10],
        ['CG', 'w_bc', 3, 3, 'n_meth',   10],
        ['CG', 'w_bc', 3, 3, 'n_unmeth', 10],
    ], columns=['motif', 'bs_strand', 'flen', 'pos', 'meth_status', 'counts']).set_index(
        ['motif', 'bs_strand', 'flen', 'pos', 'meth_status'])


    exp_mbias_stats_df = (pd.DataFrame([
        ['CG', 'w_bc', 1, 1, 0.5, 10, 10],
        ['CG', 'w_bc', 2, 1, 0.5, 10, 10],
        ['CG', 'w_bc', 2, 2, 0.5, 10, 10],
        ['CG', 'w_bc', 3, 1, 0.5, 10, 10],
        ['CG', 'w_bc', 3, 2, 0.5, 10, 10],
        ['CG', 'w_bc', 3, 3, 0.5, 10, 10], ],
        columns=['motif', 'bs_strand', 'flen', 'pos', 'beta_value', 'n_meth', 'n_unmeth'])
        .set_index(['motif', 'bs_strand', 'flen', 'pos']))
    exp_mbias_stats_df

    comp_mbias_stats_df = compute_mbias_stats_df(mbias_counts_df_stub)

    assert exp_mbias_stats_df.equals(comp_mbias_stats_df)


def run_evaluate_mbias(motifs_str, output_dir):
    user_config_file_fp = op.join(output_dir, 'user_config.toml')

    user_config = {'paths': {'mbias_counts_p': f"/home/stephen2/projects/mqc/mqc/mqc/tests/test_files/hsc_1_mbias-counts_{motifs_str}.p", }}

    with open(user_config_file_fp, 'wt') as fout:
        pytoml.dump(fout, user_config)

    subprocess.check_call(['mqc', 'evaluate_mbias',
                           '--config_file', user_config_file_fp,
                           '--motifs', motifs_str,
                           '--sample_name', 'hsc_1',
                           '--output_dir', output_dir])

def get_evaluate_mbias_config(motif_str, output_dir):
    test_dir = op.abspath(op.dirname(__file__))
    default_config_file = op.join(test_dir, '../config.default.toml')
    user_config_file = op.join(output_dir, 'user_config.toml')
    cli_params = {'motifs_str': motif_str,
                  'sample_name': 'hsc_1',
                  'sample_meta': None,
                  'output_dir': output_dir}
    config = assemble_config_vars(cli_params,
                                  default_config_file_path=default_config_file,
                                  user_config_file_path=user_config_file)
    return config

@pytest.fixture(scope='module')
def evaluate_mbias_run_all_motifs():
    """
    mbias_counts_df = pd.read_pickle(config['paths']['mbias_counts_p'])
    mbias_counts_df = mbias_counts_df.rename(columns={'Counts': 'counts'})
    mbias_counts_df.index.names = 'motif bs_strand flen pos meth_status'.split()
    mbias_counts_df = mbias_counts_df.reset_index()
    mbias_counts_df['bs_strand'] = mbias_counts_df['bs_strand'].str.lower()
    mbias_counts_df = mbias_counts_df.replace({'meth_status': {'m': 'n_meth', 'u': 'n_unmeth'}})
    mbias_counts_df = mbias_counts_df.set_index(['motif', 'bs_strand', 'flen', 'pos', 'meth_status'])
    # mbias_counts_df = mbias_counts_df.query("flen in [60, 80, 100, 120, 200]")
    mbias_counts_df.head()
    mbias_counts_df.to_pickle(config['paths']['mbias_counts_p'])
    mbias_counts_df.reset_index().to_csv(config['paths']['mbias_counts_tsv'],
                                         sep='\t', header=True, index=False)
    mbias_counts_df_all_motifs = pd.read_pickle("/home/stephen2/projects/mqc/mqc/mqc/tests/test_files/hsc_1_mbias-counts_CG-CHG-CHH.p")
    flens = list(range(50,120,10)) + list(range(120, 500, 40))
    idx = pd.IndexSlice
    mbias_counts_df_all_motifs = mbias_counts_df_all_motifs.loc[idx[:, :, flens], :]
    mbias_counts_df_all_motifs.to_pickle("/home/stephen2/projects/mqc/mqc/mqc/tests/test_files/hsc_1_mbias-counts_CG-CHG-CHH.p")
    mbias_counts_cg_only = mbias_counts_df_all_motifs.loc[['CG'], :]
    mbias_counts_cg_only.head()
    mbias_counts_cg_only.to_pickle("/home/stephen2/projects/mqc/mqc/mqc/tests/test_files/hsc_1_mbias-counts_CG.p")
    """

    output_dir = tempfile.mkdtemp()
    run_evaluate_mbias('CG-CHG-CHH', output_dir)
    config = get_evaluate_mbias_config('CG-CHG-CHH', output_dir)
    yield config
    shutil.rmtree(output_dir)

@pytest.fixture(scope='module')
def run_evaluate_mbias_for_cg_only():
    output_dir = tempfile.mkdtemp()
    run_evaluate_mbias('CG', output_dir)
    config = get_evaluate_mbias_config('CG', output_dir)
    yield config
    shutil.rmtree(output_dir)

@pytest.mark.acceptance_test
class TestMbiasEvaluateForAllMotifs:
    def test_converts_mbias_counts_to_mbias_stats_df(self, evaluate_mbias_run_all_motifs):
        config = evaluate_mbias_run_all_motifs
        mbias_stats_df = pd.read_pickle(config['paths']['mbias_stats_p'])
        assert mbias_stats_df.loc[('CG', 'w_bc', 100, 1), 'beta_value'] == 1872 / (1872 + 354)
        # TODO: add more integration tests!
        #       - masked dataframe
        #       - correct cutting sites df (which also is interpreted as correct AdjustedCuttingSites object)

@pytest.mark.acceptance_test
class TestMbiasEvaluateForCGonly:
    def test_converts_mbias_counts_to_mbias_stats_df(self, run_evaluate_mbias_for_cg_only):
        config = run_evaluate_mbias_for_cg_only
        mbias_stats_df = pd.read_pickle(config['paths']['mbias_stats_p'])
        assert mbias_stats_df.loc[('CG', 'w_bc', 100, 1), 'beta_value'] == 1872 / (1872 + 354)


def test_mask_mbias_stats_df_sets_positions_in_trimming_zone_to_nan():

    mbias_stats_df_stub = (pd.DataFrame([
        ['CG', 'w_bc', 1, 1, 0.5, 10, 10],
        ['CG', 'w_bc', 2, 1, 0.5, 10, 10],
        ['CG', 'w_bc', 2, 2, 0.5, 10, 10],
        ['CG', 'w_bc', 3, 1, 0.5, 10, 10],
        ['CG', 'w_bc', 3, 2, 0.5, 10, 10],
        ['CG', 'w_bc', 3, 3, 0.5, 10, 10], ],
        columns=['motif', 'bs_strand', 'flen', 'pos', 'beta_value', 'n_meth', 'n_unmeth'])
                          .set_index(['motif', 'bs_strand', 'flen', 'pos']))

    cutting_sites_df_stub = (pd.DataFrame([
        ['w_bc', 1, 'left_cut_end',  1],
        ['w_bc', 1, 'right_cut_end', 1],
        ['w_bc', 2, 'left_cut_end',  2],
        ['w_bc', 2, 'right_cut_end', 2],
        ['w_bc', 3, 'left_cut_end',  2],
        ['w_bc', 3, 'right_cut_end', 3],
    ], columns = ['bs_strand', 'flen', 'cut_end', 'cut_pos'])
    .set_index(['bs_strand', 'flen', 'cut_end']))

    exp_masked_df = (pd.DataFrame([
        ['CG', 'w_bc', 1, 1, 0.5, 10, 10],
        ['CG', 'w_bc', 2, 1, np.nan, np.nan, np.nan],
        ['CG', 'w_bc', 2, 2, 0.5, 10, 10],
        ['CG', 'w_bc', 3, 1, np.nan, np.nan, np.nan],
        ['CG', 'w_bc', 3, 2, 0.5, 10, 10],
        ['CG', 'w_bc', 3, 3, 0.5, 10, 10], ],
        columns=['motif', 'bs_strand', 'flen', 'pos', 'beta_value', 'n_meth', 'n_unmeth'])
                           .set_index(['motif', 'bs_strand', 'flen', 'pos']))

    computed_masked_df = mask_mbias_stats_df(mbias_stats_df_stub, cutting_sites_df_stub)

    assert exp_masked_df.equals(computed_masked_df)

@pytest.fixture()
def max_flen_cutting_df_to_array_conversion():
    return 100

@pytest.fixture()
def cutting_sites_array_computed_from_df(max_flen_cutting_df_to_array_conversion):

    cutting_sites_df_stub = (pd.DataFrame([
        ['w_bc', 1, 'left_cut_end',  1],
        ['w_bc', 1, 'right_cut_end', 1],
        ['w_bc', 2, 'left_cut_end',  2],
        ['w_bc', 2, 'right_cut_end', 2],
        ['w_bc', 3, 'left_cut_end',  2],
        ['w_bc', 3, 'right_cut_end', 3],
    ], columns = ['bs_strand', 'flen', 'cut_end', 'cut_pos'])
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
                               ((2, 1, 1), 1, ),
                               ((2, 2, 0), 2, ),
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
