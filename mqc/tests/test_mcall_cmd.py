import gzip
import os.path as op
from itertools import product

import pytest
import shutil
import subprocess
import tempfile
from textwrap import dedent

from mqc.config import assemble_config_vars

""" All motifs test
"""
# __file__ = '/home/stephen2/projects/mqc/mqc/mqc/tests/blabla.py'
TESTS_DIR = op.dirname(__file__)
TEST_BAM = op.join(TESTS_DIR, 'test_files', 'test_mcall.sorted.bam')
DEFAULT_CONFIG_FILE = op.join(TESTS_DIR, '../config.default.toml')
SAMPLE_NAME = 'hsc_rep1'
SAMPLE_META = 'population=hsc,rep=1,model=blk6'
TEST_FILES_DIR = op.join(TESTS_DIR, 'test_files')
MIN_MAPQ = 20
MIN_PHRED = 20

@pytest.fixture(scope='module')
def index_file_paths_dict():
    # TODO: add ref genome to index files to match standard index file path template
    return {
        'CG': [op.join(TEST_FILES_DIR, 'test_mcall_CG_1.bed.gz'),
               op.join(TEST_FILES_DIR, 'test_mcall_CG_2.bed.gz')],
        'CG-CHG-CHH': [op.join(TEST_FILES_DIR, 'test_mcall_CG-CHG-CHH_1.bed.gz'),
                       op.join(TEST_FILES_DIR, 'test_mcall_CG-CHG-CHH_2.bed.gz')]
    }

@pytest.fixture(scope='module')
def expected_results_dict():
    return {'1': {'CG': dedent(f"""\
                            #chrom	start	end	motif	score	strand	beta_value	n_meth	n_unmeth
                            1	11298399	11298400	CG	.	+	{7/8:.6f}	7	1
                            1	11298400	11298401	CG	.	-	{1:.6f}	3	0
                            1	11299330	11299331	CG	.	+	{1:.6f}	4	0
                            1	11299331	11299332	CG	.	-	{6/9:.6f}	6	3
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
                            2	9042611	9042612	CG	.	+	{8/10:.6f}	8	2
                            2	9042612	9042613	CG	.	-	{7/8:.6f}	7	1
                            2	9042613	9042614	CG	.	+	{7/10:.6f}	7	3
                            2	9042614	9042615	CG	.	-	{7/8:.6f}	7	1
                            """),
                  'CHG': dedent(f"""\
                            #chrom	start	end	motif	score	strand	beta_value	n_meth	n_unmeth
                            2	3281431	3281432	CHG	.	-	{0:.6f}	0	6
                            """),
                  'CHH': dedent(f"""\
                            #chrom	start	end	motif	score	strand	beta_value	n_meth	n_unmeth
                            2	3281432	3281433	CHH	.	-	{0:.6f}	0	6
                            2	3281434	3281435	CHH	.	+	{0:.6f}	0	1
                            """),
                  }
            }


@pytest.fixture(scope='module')
def config_file_path():
    config_file_dir = tempfile.mkdtemp()
    user_config_file_path = op.join(config_file_dir, 'user_config.toml')
    with open(user_config_file_path, 'wt') as fobj:
        fobj.write(dedent(f"""\
            [basic_quality_filtering]
              min_mapq = {MIN_MAPQ}
              min_phred_score = {MIN_PHRED}
            """))
    yield user_config_file_path
    shutil.rmtree(config_file_dir)

@pytest.fixture(scope='module')
def base_command_list(config_file_path):
    return ['mqc', 'call',
            '--bam', TEST_BAM,
            '--config_file', config_file_path,
            '--sample_name', SAMPLE_NAME,
            '--sample_meta', SAMPLE_META,
            '--cores', '2']

def get_text_of_call_files(output_dir, motifs_str, user_config_file):
    config = assemble_config_vars(
        command_line_args_dict=dict(output_dir=output_dir,
                                    sample_name = SAMPLE_NAME,
                                    sample_meta=SAMPLE_META,
                                    motifs_str=motifs_str),
        default_config_file_path=DEFAULT_CONFIG_FILE,
        command='call',
        user_config_file_path=user_config_file)

    def read_calls(chrom, motif):
        mcall_fp = (config['paths']['call']['meth_calls_basepath']
                    + f"_{SAMPLE_NAME}_{motif}_{chrom}.bed.gz")
        with gzip.open(mcall_fp, 'rt') as fobj:
            file_content = fobj.read()
        return file_content

    motifs = motifs_str.split('-')
    computed_calls = {(chrom, curr_motif): read_calls(chrom, curr_motif)
                      for chrom, curr_motif in product(['1', '2'], motifs)}

    return computed_calls

@pytest.fixture(scope='module')
def cg_run_file_contents(base_command_list, index_file_paths_dict, config_file_path):

    # using a non-existent output dir is intentional
    # mqc should create the output dir recursively
    tmpdir = tempfile.mkdtemp()
    cg_only_output_dir = op.join(tmpdir, 'hsc_1')
    subprocess.check_call(base_command_list
                          + ['--output_dir', cg_only_output_dir]
                          + index_file_paths_dict['CG'])
    computed_calls = get_text_of_call_files(cg_only_output_dir,
                                            motifs_str='CG',
                                            user_config_file=config_file_path)
    yield computed_calls
    shutil.rmtree(tmpdir)

@pytest.fixture(scope='module')
def all_motifs_run_file_contents(base_command_list, index_file_paths_dict, config_file_path):

    # using a non-existent output dir is intentional
    # mqc should create the output dir recursively
    tmpdir = tempfile.mkdtemp()
    all_motifs_output_dir = op.join(tmpdir, 'all-motifs')
    subprocess.check_call(base_command_list
                          + ['--output_dir', all_motifs_output_dir]
                          + index_file_paths_dict['CG-CHG-CHH'])
    computed_calls = get_text_of_call_files(all_motifs_output_dir,
                                            motifs_str='CG-CHG-CHH',
                                            user_config_file=config_file_path)
    yield computed_calls
    shutil.rmtree(tmpdir)

@pytest.mark.acceptance_test
@pytest.mark.parametrize('motif_str,chrom',
                         (('CG', '1'),
                          ('CG', '2')))
def test_makes_correct_calls_for_cg_only_run(
        motif_str, chrom, cg_run_file_contents, expected_results_dict):

    computed_calls = cg_run_file_contents[(chrom, motif_str)]
    expected_calls = expected_results_dict[chrom][motif_str]

    msg = dedent(f"""\
    Computed
    --------
    {computed_calls}
    
    Expected
    --------
    {expected_calls}
    """)

    assert computed_calls == expected_calls, msg

@pytest.mark.acceptance_test
@pytest.mark.parametrize('motif_str,chrom',
                         (('CG',  '1'),
                          ('CG',  '2'),
                          ('CHG', '1'),
                          ('CHG', '2'),
                          ('CHH', '1'),
                          ('CHH', '2'),)
                         )
def test_makes_correct_calls_for_all_motifs_run(
        motif_str, chrom, all_motifs_run_file_contents, expected_results_dict):

    computed_calls = all_motifs_run_file_contents[(chrom, motif_str)]
    expected_calls = expected_results_dict[chrom][motif_str]

    msg = dedent(f"""\
    Computed
    --------
    {computed_calls}
    
    Expected
    --------
    {expected_calls}
    """)

    assert computed_calls == expected_calls, msg
