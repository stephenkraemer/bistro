import gzip
import shutil
import tempfile
from textwrap import dedent

import pytest
import os.path as op
import subprocess

from mqc.config import assemble_config_vars

""" All motifs test
"""
TESTS_DIR = op.dirname(__file__)
TEST_BAM = op.join(TESTS_DIR, 'test_files', 'test_mcall.sorted.bam'),
SAMPLE_NAME = 'hsc_rep1'
SAMPLE_META = 'population=hsc,rep=1,model=blk6'
DEFAULT_CONFIG_FILE = op.join(TESTS_DIR, '../config.default.toml')
TEST_FILES_DIR = op.join(TESTS_DIR, 'test_files')
INDEX_FILES = {
    'CG': [op.join(TEST_FILES_DIR, 'test_mcall_CG_1.bed.gz'),
           op.join(TEST_FILES_DIR, 'test_mcall_CG_2.bed.gz')],
    'CG-CHG-CHH': [op.join(TEST_FILES_DIR, 'test_mcall_CG-CHG-CHH_1.bed.gz'),
                   op.join(TEST_FILES_DIR, 'test_mcall_CG-CHG-CHH_2.bed.gz')]
}
USER_FLEN_MAX=300

@pytest.fixture(scope='module')
def expected_results_dict():
    return {
        'CG': {
            1: dedent("""\
                #chrom	start	end	motif	score	strand	beta	n_meth	u_unmeth
                1	11298399	11298400	CG	.	+	0.875	7	1
                1	11298400	11298401	CG	.	-	1	3	0
                1	11299330	11299331	CG	.	+	1	4	0
                1	11299331	11299332	CG	.	-	0.6666666666	6	3"""),
            2: dedent("""\
                #chrom	start	end	motif	score	strand	beta	n_meth	u_unmeth
                2	9042611	9042612	CG	.	+	0.8	8	2
                2	9042612	9042613	CG	.	-	0.875	7	1""")
        }
    }


@pytest.fixture(scope='module')
def output_dir():
    tmpdir = tempfile.mkdtemp()
    yield tmpdir
    # TODO: is this safe against exceptions?
    shutil.rmtree(tmpdir)

# TODO: do with and without user config file
@pytest.fixture(scope='module')
def user_config_file():
    tmpdir = tempfile.mkdtemp()
    user_config_file_path = op.join(tmpdir, 'user_config.toml')
    with open(user_config_file_path, 'wt') as fobj:
        fobj.write(dedent(f"""\
            [trimming]
            max_flen_considered_for_trimming = {USER_FLEN_MAX}
            
            """))
    yield user_config_file_path
    # TODO: does this also remove after an exception occured?
    shutil.rmtree(tmpdir)

@pytest.fixture(scope='module')
def cg_only_mcall_run(request, user_config_file, output_dir):
    index_files = INDEX_FILES['CG']
    # Note that this output dir does not exist before execution of
    # mqc call. This should be fine and a common case!
    command_list = ['mqc', 'call',
                    '--bam', TEST_BAM,
                    '--config_file', user_config_file,
                    '--output_dir', output_dir,
                    '--sample_name', 'hsc_1',
                    '--sample_meta', 'population=HSC,rep=rep1',
                    '--cores', '2'] + index_files,
    subprocess.check_call(command_list)
    return 0

@pytest.fixture(scope='module')
def default_paths_cg_only_run(output_dir):
    """Reconstructs config file used in mqc stats run to allow usage of config['paths']"""
    config = assemble_config_vars(
        command_line_args_dict=dict(output_dir=output_dir,
                                    sample_name = SAMPLE_NAME,
                                    sample_meta=SAMPLE_META,
                                    # motifs_str is computed and supplied to
                                    # assemble_config_vars within collect_stats
                                    motifs_str='CG'),
        default_config_file_path=DEFAULT_CONFIG_FILE,
        user_config_file_path='')
    return config['paths']

def read_calls(fp):
    with gzip.open(fp, 'rt') as fobj:
        res = fobj.read()
    return res

@pytest.mark.acceptance_test
@pytest.mark.parametrize("motif_str,chrom",
                         ('CG', '1'),
                         ('CG', '2'),
                         ('CG-CHG-CHH', '1'),
                         ('CG-CHG-CHH', '2'),
                         )
def test_makes_correct_calls_in_CG_only_run(
        motif_str, chrom, cg_only_mcall_run,
        default_paths_cg_only_run, expected_results_dict):

    mcall_fp = (default_paths_cg_only_run['mcall_base_name']
                + f"{motif_str}_{chrom}.bed.gz")
    computed_text = read_calls(mcall_fp)
    expected_text = expected_results_dict[motif_str][chrom]

    msg = dedent(f"""\
    Computed
    --------
    {computed_text}
    
    Expected
    --------
    {expected_text}
    """)

    assert expected_text == computed_text, msg





