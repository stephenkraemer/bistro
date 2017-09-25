from mqc.beta_value import StratifiedBetaCounter
from unittest.mock import MagicMock
import numpy as np
import pytest
import tempfile
import os.path as op
import pandas as pd
from collections import defaultdict
import subprocess
import pytoml

TESTS_DIR = op.dirname(__file__)
TEST_BAM = op.join(TESTS_DIR, 'test_files', 'test_mcall.sorted.bam')
SAMPLE_NAME = 'hsc_rep1'
TEST_FILES_DIR = op.join(TESTS_DIR, 'test_files')
QC_STATS_DIR = 'qc_stats'
MOTIF_SETS = ['CG', 'CG-CHG-CHH']
TEST_CHROMS = ['1', '2']

CONFIG = defaultdict(dict)
CONFIG['paths']['stratified_beta_counts'] = ''

class TestStratifiedBetaCounter:
    def test_stratified_beta_counter_array(self):

        strat_beta_counter = StratifiedBetaCounter(CONFIG)

        for beta_arr in ([0, 0, 0, 0, 0, 0, 1],
                         [0, 0, 0, 0, 0, 0, 0.5],
                         [0, 0, 0, 0, 0, 0.5, 0.5],
                         [0, 0, 0, 0.2, 0, 0, 0]):
            motif_pileup = MagicMock(strat_beta_arr=np.array(beta_arr))
            strat_beta_counter.process(motif_pileup)

        assert strat_beta_counter.counter_array[6, 1000] == 1
        assert strat_beta_counter.counter_array[6, 500] == 2
        assert strat_beta_counter.counter_array[3, 200] == 1
        assert strat_beta_counter.counter_array[5, 1000] == 0

    def test_handles_inf_nan_correctly(self):

        strat_beta_counter = StratifiedBetaCounter(CONFIG)

        motif_pileup = MagicMock(strat_beta_arr=np.array([np.inf, 0, 0, 0, 0, 0, np.nan]))
        strat_beta_counter.process(motif_pileup)

        assert strat_beta_counter.counter_array[0].sum() == 0
        assert strat_beta_counter.counter_array[6].sum() == 0


    @pytest.mark.acceptance_test
    @pytest.mark.parametrize('motifs', MOTIF_SETS)
    def test_strat_beta_counter_file_output(self, motifs):
        with tempfile.TemporaryDirectory() as tmpdir:
            idx_files = [op.join(TEST_FILES_DIR, 'test_mcall_' + motifs + '_' + chrom + '.bed.gz') for chrom in TEST_CHROMS]
            output_dir = op.join(tmpdir, SAMPLE_NAME)

            test_config_path = op.join(tmpdir, 'test_config.toml')
            computed_df_basename = f"{SAMPLE_NAME}_stratified-beta-counts_{motifs}.p"
            user_config = {'paths': {'stratified_beta_counts': f"{QC_STATS_DIR}/{computed_df_basename}", }}

            with open(test_config_path, 'w') as fout:
                pytoml.dump(fout, user_config)

            subprocess.run(['mqc', 'call', '--bam', TEST_BAM, idx_files[0], idx_files[1], '--config_file', test_config_path, '--output_dir', output_dir, '--sample_name', SAMPLE_NAME, '--strat_beta_dist'])

            counter_df_tsv = pd.read_csv(op.join(output_dir, QC_STATS_DIR, f'{computed_df_basename}.tsv'), sep='\t').set_index(['strat_idx', 'beta_value'])
            counter_df_p = pd.read_pickle(op.join(output_dir, QC_STATS_DIR, f'{computed_df_basename}.p'))

            expected_df = pd.read_pickle(op.join(TEST_FILES_DIR, f'test_strat_beta_counts_{motifs}.p'))

            for df_type, computed_df in [['tsv', counter_df_tsv], ['pickled', counter_df_p]]:

                msg = f"Computed ({df_type} dataframe): {computed_df} \nExpected coverage: {expected_df}"
                assert computed_df.loc[expected_df.index].equals(expected_df), msg
