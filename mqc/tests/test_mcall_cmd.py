import gzip
import shutil
import tempfile
from textwrap import dedent

import pytest
import os
import os.path as op
import subprocess

from mqc.config import assemble_config_vars

""" All motifs test
"""
# __file__ = '/home/stephen2/projects/mqc/mqc/mqc/tests/blabla.py'
TESTS_DIR = op.dirname(__file__)
TEST_BAM = op.join(TESTS_DIR, 'test_files', 'test_mcall.sorted.bam')
SAMPLE_NAME = 'hsc_rep1'
SAMPLE_META = 'population=hsc,rep=1,model=blk6'
DEFAULT_CONFIG_FILE = op.join(TESTS_DIR, '../config.default.toml')
TEST_FILES_DIR = op.join(TESTS_DIR, 'test_files')
# TODO: add ref genome to index files to match standard index file path template
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
            '1': dedent(f"""\
                #chrom	start	end	motif	score	strand	beta_value	n_meth	n_unmeth
                1	11298399	11298400	CG	.	+	{7/8:.6f}	7	1
                1	11298400	11298401	CG	.	-	{1:.6f}	3	0
                1	11298418	11298419	CG	.	+	{1:.6f}	9	0
                1	11298419	11298420	CG	.	-	{1:.6f}	4	0
                """),
                # 1	11298883	11298884	CG	.	+	{1:.6f}	12	0
                # 1	11298884	11298885	CG	.	-	{1:.6f}	14	0
                # 1	11299330	11299331	CG	.	+	{1:.6f}	4	0
                # 1	11299331	11299332	CG	.	-	{6/9:.6f}	6	3
            '2': dedent(f"""\
                #chrom	start	end	motif	score	strand	beta_value	n_meth	n_unmeth
                2	9042611	9042612	CG	.	+	{8/10:.6f}	8	2
                2	9042612	9042613	CG	.	-	{7/8:.6f}	7	1
                2	9042613	9042614	CG	.	+	{7/10:.6f}	7	3
                2	9042614	9042615	CG	.	-	{7/8:.6f}	7	1
                """)
        },
        'CG-CHG-CHH': {
            '1': {'CG': dedent(f"""\
                            #chrom	start	end	motif	score	strand	beta_value	n_meth	n_unmeth
                            """),
                  'CHG': dedent(f"""\
                            #chrom	start	end	motif	score	strand	beta_value	n_meth	n_unmeth
                            1	3258371	3258372	CHG	.	+	{0:.6f}	0	6
                            1	3258373	3258374	CHG	.	-	{0:.6f}	0	4
                            """),
                  'CHH': dedent(f"""\
                            #chrom	start	end	motif	score	strand	beta_value	n_meth	n_unmeth
                            1	3258374	3258375	CHH	.	-	{0:.6f}	0	4
                            1	3258376	3258377	CHH	.	+	{0:.6f}	0	7
                            """),
                  },
            '2': {'CG': dedent(f"""\
                            #chrom	start	end	motif	score	strand	beta_value	n_meth	n_unmeth
                            2	3281429	3281430	CG	.	+	{1:.6f}	2	0
                            2	3281430	3281431	CG	.	-	{1:.6f}	6	0
                            """),
                  'CHG': dedent(f"""\
                            2	3281431	3281432	CHG	.	-	{0:.6f}	0	6
                            """),
                  'CHH': dedent(f"""\
                            2	3281432	3281433	CHH	.	-	{0:.6f}	0	6
                            2	3281434	3281435	CHH	.	+	{0:.6f}	0	1
                            """),
                  }
        },

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
def base_command_list(user_config_file):
    # Note that this output dir does not exist before execution of
    # mqc call. This should be fine and a common case!
    return ['mqc', 'call',
            '--bam', TEST_BAM,
            '--config_file', user_config_file,
            '--sample_name', SAMPLE_NAME,
            '--sample_meta', SAMPLE_META,
            '--cores', '2']

@pytest.fixture(scope='module')
def cg_only_mcall_run(output_dir, base_command_list):
    cg_only_output_dir = op.join(output_dir, 'cg_only')
    subprocess.check_call(base_command_list
                          + ['--output_dir', cg_only_output_dir]
                          + INDEX_FILES['CG'])
    return cg_only_output_dir

@pytest.fixture(scope='module')
def all_motifs_mcall_run(output_dir, base_command_list):
    all_motifs_output_dir = op.join(output_dir, 'all_motifs')
    subprocess.check_call(base_command_list
                          + ['--output_dir', all_motifs_output_dir]
                          + INDEX_FILES['CG-CHG-CHH'])
    return all_motifs_output_dir

@pytest.fixture(scope='module')
def default_paths_cg_only_run(cg_only_mcall_run):
    """Reconstructs config file used in mqc stats run to allow usage of config['paths']"""
    config = assemble_config_vars(
        command_line_args_dict=dict(output_dir=cg_only_mcall_run,
                                    sample_name = SAMPLE_NAME,
                                    sample_meta=SAMPLE_META,
                                    # motifs_str is computed and supplied to
                                    # assemble_config_vars within collect_stats
                                    motifs_str='CG'),
        default_config_file_path=DEFAULT_CONFIG_FILE,
        command='call',
        user_config_file_path='')
    return config['paths']

@pytest.fixture(scope='module')
def default_paths_all_motifs_run(all_motifs_mcall_run):
    """Reconstructs config file used in mqc stats run to allow usage of config['paths']"""
    config = assemble_config_vars(
        command_line_args_dict=dict(output_dir=all_motifs_mcall_run,
                                    sample_name = SAMPLE_NAME,
                                    sample_meta=SAMPLE_META,
                                    # motifs_str is computed and supplied to
                                    # assemble_config_vars within collect_stats
                                    motifs_str='CG-CHG-CHH'),
        default_config_file_path=DEFAULT_CONFIG_FILE,
        command='call',
        user_config_file_path='')
    return config['paths']

def read_calls(fp, end_line, start_line=0):
    with gzip.open(fp, 'rt') as fobj:
        lines = fobj.readlines()
    text = ''.join(lines[start_line:end_line])
    return text

@pytest.mark.acceptance_test
@pytest.mark.parametrize("motif_str,chrom",
                         (('CG', '1'),
                          ('CG', '2')))
def test_makes_correct_calls_for_cg_only_run(
        motif_str, chrom, cg_only_mcall_run,
        default_paths_cg_only_run, expected_results_dict):

    mcall_fp = (default_paths_cg_only_run['call']['meth_calls_basepath']
                + f"_{SAMPLE_NAME}_{motif_str}_{chrom}.bed.gz")
    computed_text = read_calls(mcall_fp, end_line=5)
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

@pytest.mark.acceptance_test
@pytest.mark.parametrize("motif_str,chrom,start_line,end_line",
                         (('CG', '1', 0, 5),
                         ('CHG', '1', 0, 3),
                         ('CHH', '1', 0, 3),
                          ('CG', '2', 0, 3),
                          ('CHG', '2', 3, 4),
                          ('CHH', '2', 24, 26)
                          ))
def test_makes_correct_calls_for_different_motifs(
        motif_str, chrom, start_line, end_line, default_paths_all_motifs_run,
        expected_results_dict):

    mcall_fp = (default_paths_all_motifs_run['call']['meth_calls_basepath']
                + f"_{SAMPLE_NAME}_{motif_str}_{chrom}.bed.gz")
    computed_text = read_calls(mcall_fp, start_line=start_line, end_line=end_line)
    expected_text = expected_results_dict['CG-CHG-CHH'][chrom][motif_str]

    msg = dedent(f"""\
    Computed
    --------
    {computed_text}
    
    Expected
    --------
    {expected_text}
    """)

    assert expected_text == computed_text, msg
