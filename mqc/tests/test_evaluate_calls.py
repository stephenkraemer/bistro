from textwrap import dedent

import pytest
import os.path as op
import tempfile
import subprocess

TESTS_DIR = op.dirname(__file__)
TEST_FILES_DIR = op.join(TESTS_DIR, 'test_files')
TEST_BAM = op.join(TEST_FILES_DIR, 'test_mcall.sorted.bam')

SAMPLE_NAME = 'hsc_1'
MOTIF_SETS = ['CG', 'CG-CHG-CHH']
TEST_CHROMS = ['1', '2']


@pytest.mark.parametrize('motifs', MOTIF_SETS)
def test_evaluate_calls_cmd(motifs):
    with tempfile.TemporaryDirectory() as out_dir:

        test_config_file = op.join(out_dir, 'test_conf.toml')
        with open(test_config_file, mode='w') as config_fobj:
            config_fobj.write(dedent(f"""\
            [paths]
            qc_stats_dir = 'qc_stats'
            cov_counts_p = '{TEST_FILES_DIR}/test_coverage-counts_{motifs}.p'
            """))

        ret_code = subprocess.run(['mqc',
                                   'evaluate_calls',
                                   '--config_file', test_config_file,
                                   '--output_dir', out_dir,
                                   '--sample_name', SAMPLE_NAME,
                                   '--motifs', motifs.replace('-', ',')
                                   ], check=False).returncode
        assert ret_code == 0

