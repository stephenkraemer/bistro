import os.path as op
import subprocess
import tempfile
from collections import defaultdict
from unittest.mock import MagicMock

import pandas as pd
import pytest
from mqc.coverage import CoverageCounter

CONFIG = defaultdict(dict)
CONFIG['run']['motifs'] = ['cg', 'chg']
CONFIG['coverage_analysis']['max_per_motif_cov'] = 20
CONFIG['paths']['cov_counts'] = None

TESTS_DIR = op.dirname(__file__)
TEST_BAM = op.join(TESTS_DIR, 'test_files', 'test_mcall.sorted.bam')
SAMPLE_NAME = 'hsc_rep1'
TEST_FILES_DIR = op.join(TESTS_DIR, 'test_files')
QC_STATS_DIR = 'qc_stats'

MOTIF_SETS = ['CG', 'CG-CHG-CHH']
TEST_CHROMS = ['1', '2']

EXPECTED_COVERAGE = {'CG': pd.DataFrame({'motif': ['CG', 'CG', 'CG'],
                                         'coverage': [1, 2, 8],
                                         'counts': [0, 2, 1]}),
                     'CG-CHG-CHH': pd.DataFrame({'motif': ['CG', 'CG', 'CHG', 'CHG', 'CHH'],
                                                 'coverage': [2, 8, 4, 2, 1],
                                                 'counts': [2, 1, 1, 0, 1]})}

class TestCoverageCounterProcessing:
    def test_updates_cov_counter_array(self):

        idx_motif_tuples = enumerate(CONFIG['run']['motifs'])
        motif_idx_dict = {motif: i for i, motif in idx_motif_tuples}
        coverage_counter = CoverageCounter(CONFIG)

        for n_meth, n_total, motif in ([3, 5, 'cg'],
                                        [4, 5, 'cg'],
                                        [3, 10, 'chg'],
                                        [20, 40, 'chg']):
            motif_pileup = MagicMock(n_meth=n_meth, n_total=n_total)
            motif_pileup.idx_pos.motif = motif
            coverage_counter.process(motif_pileup)

        assert coverage_counter.counter_array[motif_idx_dict['cg'], 5] == 2
        assert coverage_counter.counter_array[motif_idx_dict['cg'], 10] == 0
        assert coverage_counter.counter_array[motif_idx_dict['chg'], 10] == 1
        assert coverage_counter.counter_array[motif_idx_dict['chg'], 20] == 1

    def test_cov_counter_pools_overshoot(self):

        idx_motif_tuples = enumerate(CONFIG['run']['motifs'])
        motif_idx_dict = {motif: i for i, motif in idx_motif_tuples}
        coverage_counter = CoverageCounter(CONFIG)

        motif_pileup = MagicMock(n_meth=20, n_total=40)
        motif_pileup.idx_pos.motif = 'chg'
        coverage_counter.process(motif_pileup)

        assert coverage_counter.counter_array[motif_idx_dict['chg'], 20] == 1
        with pytest.raises(IndexError):
            assert coverage_counter.counter_array[motif_idx_dict['chg'], 40]


    @pytest.mark.acceptance_test
    @pytest.mark.parametrize('motifs', MOTIF_SETS)
    def test_cov_counter_file_output(self, motifs):
        with tempfile.TemporaryDirectory() as tmpdir:
            idx_files = [op.join(TEST_FILES_DIR, 'test_mcall_' + motifs + '_' + chrom + '.bed.gz') for chrom in TEST_CHROMS]
            output_dir = op.join(tmpdir, SAMPLE_NAME)

            subprocess.run(['mqc', 'call', '--bam', TEST_BAM, idx_files[0], idx_files[1], '--output_dir', output_dir, '--sample_name', SAMPLE_NAME])

            counter_df_tsv = pd.read_csv(op.join(output_dir, QC_STATS_DIR, f"{SAMPLE_NAME}_coverage-counts_{motifs}.tsv"), sep='\t').set_index(['motif', 'coverage'])
            counter_df_p = pd.read_pickle(op.join(output_dir, QC_STATS_DIR, f"{SAMPLE_NAME}_coverage-counts_{motifs}.p"))

            expected_df = EXPECTED_COVERAGE[motifs].set_index(['motif', 'coverage'])

            for df_type, computed_df in [['tsv', counter_df_tsv], ['pickled', counter_df_p]]:

                msg = f"Computed ({df_type} dataframe): {computed_df} \nExpected coverage: {expected_df}"
                assert computed_df.loc[expected_df.index].equals(expected_df), msg


