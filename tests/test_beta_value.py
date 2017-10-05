from mqc.beta_value import StratifiedBetaCounter, BedFile
from unittest.mock import MagicMock, patch, mock_open
import numpy as np
import pytest
import tempfile
import os.path as op
import pandas as pd
from collections import defaultdict
import subprocess
import mqc.flag_and_index_values as mfl
from pandas.util.testing import assert_frame_equal

import pytoml

TESTS_DIR = op.dirname(__file__)
TEST_BAM = op.join(TESTS_DIR, 'test_files', 'test_mcall.sorted.bam')
SAMPLE_NAME = 'hsc_rep1'
TEST_FILES_DIR = op.join(TESTS_DIR, 'test_files')
QC_STATS_DIR = 'qc_stats'
MOTIF_SETS = ['CG', 'CG-CHG-CHH']
TEST_CHROMS = ['1', '2']

roi_strats = mfl.region_of_interest_indices

CONFIG = defaultdict(dict)
CONFIG['paths']['stratified_beta_counts'] = ''
CONFIG['plots']['stratified_beta_hist_bins'] = 1000
CONFIG['run']['motifs'] = ['CG', 'CHG', 'CHH']
CONFIG['run']['roi_index_files'] = []

class TestRoiBedReader:
    def test_chrom_filtering(self):
        chrom='1'
        with patch('mqc.beta_value.open', mock_open(read_data='1\t1\t4\n1\t2\t5\n2\t3\t6', )):
            roi_iterator = BedFile(bedPath='bedfile.bed', chrom=chrom, idx=1)
            filtered_tuples = [x for x in roi_iterator.intervals]
            assert filtered_tuples == [(1, 4), (2, 5)]


class TestStratifiedBetaCounter:
    def test_stratified_beta_counter_array(self):

        with patch('mqc.beta_value.open', mock_open(read_data='\n')):

            strat_beta_counter = StratifiedBetaCounter(CONFIG, '1')
            motif_idx = strat_beta_counter.motif_idx_dict['CG']

            for beta_arr in ([0, 0, 0, 0, 0, 0, 1],
                             [0, 0, 0, 0, 0, 0, 0.5],
                             [0, 0, 0, 0, 0, 0.5, 0.5],
                             [0, 0, 0, 0.2, 0, 0, 0]):
                motif_pileup = MagicMock(strat_beta_arr=np.array(beta_arr))
                motif_pileup.idx_pos.motif = 'CG'
                strat_beta_counter.process(motif_pileup)

            assert strat_beta_counter.counter_array[motif_idx, roi_strats.whole_genome, 6, 1000] == 1
            assert strat_beta_counter.counter_array[motif_idx, roi_strats.whole_genome, 6, 500 ] == 2
            assert strat_beta_counter.counter_array[motif_idx, roi_strats.whole_genome, 3, 200 ] == 1
            assert strat_beta_counter.counter_array[motif_idx, roi_strats.whole_genome, 5, 1000] == 0

    def test_correct_roi_stratifying(self):

        #TODO: create additional test to make sure bedfiles with both chr prefix and without it are supported
        with patch('mqc.beta_value.open', mock_open(read_data='1\t10\t20\nchr1\t30\t40\n')):
            config_with_roi = CONFIG.copy()
            config_with_roi['run'] = {'roi_index_files': 'bedfile.bed', 'motifs': CONFIG['run']['motifs']}

            strat_beta_counter = StratifiedBetaCounter(config_with_roi, '1')
            motif_idx = strat_beta_counter.motif_idx_dict['CG']
            roi_idx = strat_beta_counter.roi_list[0].idx

            for start, end, beta in [(10, 15, 1), (15, 25, 1), (20, 25, 0), (45, 50, 1)]:
                motif_pileup = MagicMock(strat_beta_arr=np.array((0, 0, 0, 0, 0, 0, beta)))
                motif_pileup.idx_pos.start = start
                motif_pileup.idx_pos.end = end
                motif_pileup.idx_pos.motif = 'CG'
                strat_beta_counter.process(motif_pileup)

            assert strat_beta_counter.counter_array[motif_idx, roi_idx, 6, 1000] == 2
            assert strat_beta_counter.counter_array[motif_idx, roi_strats.whole_genome, 6, 1000] == 3
            assert strat_beta_counter.counter_array[motif_idx, roi_strats.whole_genome, 6, 0] == 1
            assert strat_beta_counter.counter_array[motif_idx, roi_idx, 6, 0] == 0

    def test_correct_roi_stratifying_with_multiple_files(self):

        beddata_1 = "chr1 5 25\nchr1 45 50\n".replace(' ', '\t')
        beddata_2 = "chr1 15 35\nchr1 45 60\n".replace(' ', '\t')

        m_open = bedhandle_1 = mock_open(read_data=beddata_1)
        bedhandle_2 = mock_open(read_data=beddata_2)
        m_open.side_effect = [bedhandle_1.return_value, bedhandle_2.return_value]

        with patch('mqc.beta_value.open', m_open):
            config_with_roi = CONFIG.copy()
            config_with_roi['run'] = {'roi_index_files': 'bed1.bed,bed2.bed', 'motifs': CONFIG['run']['motifs']}

            strat_beta_counter = StratifiedBetaCounter(config_with_roi, '1')
            motif_idx = strat_beta_counter.motif_idx_dict['CG']
            roi_indices = {roi.name: roi.idx for roi in strat_beta_counter.roi_list}
            roi_indices['whole_genome'] = mfl.region_of_interest_indices.whole_genome

            pileup_data = [(10, 0), (20, 0.1), (30, 0.2), (30, 0.3), (40, 0.4), (55, 0.5), (65, 0.6)]

            for pos, beta in pileup_data:
                motif_pileup = MagicMock(strat_beta_arr=np.array((0, 0, 0, 0, 0, 0, beta)))
                motif_pileup.idx_pos.start = pos
                motif_pileup.idx_pos.end = pos+1
                motif_pileup.idx_pos.motif = 'CG'
                strat_beta_counter.process(motif_pileup)

            expected_slices = {'whole_genome': [1, 1, 1, 1, 1, 1, 1],
                               'bed1': [1, 1, 0, 0, 0, 0, 0],
                               'bed2': [0, 1, 1, 1, 0, 1, 0]}

            for name, expected_array_slice in expected_slices.items():
                computed_array_slice = strat_beta_counter.counter_array[motif_idx,
                                                                        roi_indices[name],
                                                                        6,
                                                                        [int(x[1]*strat_beta_counter.bins) for x in pileup_data]]
                np.testing.assert_array_equal(computed_array_slice, np.array(expected_array_slice))



    def test_correct_motif_stratifying(self):

        with patch('mqc.beta_value.open', mock_open(read_data='\n')):
            strat_beta_counter = StratifiedBetaCounter(CONFIG, '1')
            strat_beta_counter.roi_list = []
            motif_idx = strat_beta_counter.motif_idx_dict['CHG']
            wrong_motif_idx = strat_beta_counter.motif_idx_dict['CHH']

            motif_pileup = MagicMock(strat_beta_arr=np.array((0, 0, 0, 0, 0, 0, 1)))
            motif_pileup.idx_pos.motif = 'CHG'
            strat_beta_counter.process(motif_pileup)

            assert strat_beta_counter.counter_array[motif_idx, roi_strats.whole_genome, 6, 1000] == 1
            assert strat_beta_counter.counter_array[wrong_motif_idx, roi_strats.whole_genome, 6, 1000] == 0


    def test_handles_inf_nan_correctly(self):

        strat_beta_counter = StratifiedBetaCounter(CONFIG, '1')
        strat_beta_counter.roi_list = []
        motif_idx = strat_beta_counter.motif_idx_dict['CG']

        motif_pileup = MagicMock(strat_beta_arr=np.array([np.inf, 0, 0, 0, 0, 0, np.nan]))
        motif_pileup.idx_pos.motif = 'CG'
        strat_beta_counter.process(motif_pileup)

        assert strat_beta_counter.counter_array[motif_idx, roi_strats.whole_genome, 0].sum() == 0
        assert strat_beta_counter.counter_array[motif_idx, roi_strats.whole_genome, 6].sum() == 0


    def test_writes_metadata_to_counter_dataframe(self):

        config_with_meta = CONFIG.copy()
        config_with_meta['sample']['metadata_key1'] = 'value1'
        config_with_meta['sample']['metadata_key2'] = 'value2'

        strat_beta_counter = StratifiedBetaCounter(config_with_meta)
        df = strat_beta_counter.get_dataframe()

        assert df.index.names == ['metadata_key1', 'metadata_key2', 'motif', 'roi_type', 'bs_strand', 'beta_value']
        df = df.reset_index()

        assert (df['metadata_key1'] == 'value1').all()
        assert (df['metadata_key2'] == 'value2').all()

    # TODO: Create additional acceptance test with roi stratifying

    @pytest.mark.acceptance_test
    @pytest.mark.parametrize('motifs', MOTIF_SETS)
    def test_strat_beta_counter_file_output(self, motifs):
        with tempfile.TemporaryDirectory() as tmpdir:
            idx_files = [op.join(TEST_FILES_DIR, 'test_mcall_' + motifs + '_' + chrom + '.bed.gz') for chrom in TEST_CHROMS]
            output_dir = '/home/mattausc/Documents/sandbox05'#tmpdir

            test_config_path = op.join(tmpdir, 'test_config.toml')
            computed_df_basename = f"{SAMPLE_NAME}_stratified-beta-counts_{motifs}"
            user_config = {'paths': {'stratified_beta_counts': f"{QC_STATS_DIR}/{computed_df_basename}", }}

            with open(test_config_path, 'w') as fout:
                pytoml.dump(fout, user_config)

            subprocess.run(['mqc', 'call',
                            '--bam', TEST_BAM,
                            idx_files[0],
                            idx_files[1],
                            '--config_file', test_config_path,
                            '--output_dir', output_dir,
                            '--sample_name', SAMPLE_NAME,
                            '--strat_beta_dist'])

            counter_df_tsv = pd.read_csv(op.join(output_dir, QC_STATS_DIR, f'{computed_df_basename}.tsv'), sep='\t')\
                .set_index(['name', 'motif', 'roi_type', 'bs_strand', 'beta_value'])
            counter_df_p = pd.read_pickle(op.join(output_dir, QC_STATS_DIR, f'{computed_df_basename}.p'))

            expected_df = pd.read_pickle(op.join(TEST_FILES_DIR, f'test_strat_beta_counts_{motifs}.p'))

            for computed_df in [counter_df_tsv, counter_df_p]:

                assert_frame_equal(computed_df, expected_df, check_names=False)

