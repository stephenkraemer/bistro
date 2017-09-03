import gzip
import os.path as op
from itertools import product

import pickle
import pytest
import shutil
import subprocess
import tempfile
from textwrap import dedent

import pytoml

from mqc.config import assemble_config_vars
from mqc.mbias import FixedRelativeCuttingSites
from mqc.utils import get_resource_abspath

""" All motifs test
"""
# __file__ = '/home/stephen2/projects/mqc/mqc/mqc/tests/blabla.py'
TESTS_DIR = op.dirname(__file__)
TEST_BAM = op.join(TESTS_DIR, 'test_files', 'test_mcall.sorted.bam')
DEFAULT_CONFIG_FILE = get_resource_abspath('config.default.toml')
SAMPLE_NAME = 'hsc_rep1'
SAMPLE_META = 'population=hsc,rep=1,model=blk6'
TEST_FILES_DIR = op.join(TESTS_DIR, 'test_files')
MIN_MAPQ = 41
MIN_PHRED = 35

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
                            1	11298399	11298400	CG	.	+	{1:.6f}	5	0
                            1	11298400	11298401	CG	.	-	nan	0	0
                            1	11299330	11299331	CG	.	+	{1:.6f}	3	0
                            1	11299331	11299332	CG	.	-	{1:.6f}	4	0
                            """),
                  'CHG': dedent(f"""\
                            #chrom	start	end	motif	score	strand	beta_value	n_meth	n_unmeth
                            1	3258371	3258372	CHG	.	+	{0:.6f}	0	3
                            1	3258373	3258374	CHG	.	-	{0:.6f}	0	3
                            """),
                  'CHH': dedent(f"""\
                            #chrom	start	end	motif	score	strand	beta_value	n_meth	n_unmeth
                            1	3258374	3258375	CHH	.	-	{0:.6f}	0	3
                            1	3258376	3258377	CHH	.	+	{0:.6f}	0	3
                            """),
                  },
            '2': {'CG': dedent(f"""\
                            #chrom	start	end	motif	score	strand	beta_value	n_meth	n_unmeth
                            2	9042611	9042612	CG	.	+	{1:.6f}	4	0
                            2	9042612	9042613	CG	.	-	{1:.6f}	1	0
                            2	9042613	9042614	CG	.	+	{1:.6f}	5	0
                            2	9042614	9042615	CG	.	-	{1:.6f}	1	0
                            """),
                  'CHG': dedent(f"""\
                            #chrom	start	end	motif	score	strand	beta_value	n_meth	n_unmeth
                            2	3281431	3281432	CHG	.	-	{0:.6f}	0	3
                            """),
                  'CHH': dedent(f"""\
                            #chrom	start	end	motif	score	strand	beta_value	n_meth	n_unmeth
                            2	3281432	3281433	CHH	.	-	{0:.6f}	0	3
                            2	3281434	3281435	CHH	.	+	{0:.6f}	0	1
                            """),
                  }
            }

@pytest.fixture(scope="module")
def cutting_sites_obj_fp():
    config_stub = pytoml.loads(dedent("""
            [trimming]
              max_flen_considered_for_trimming = 500
              [trimming.relative_to_fragment_ends_dict]
                w_bc    = [0, 20]
                c_bc    = [0, 20]
                w_bc_rv = [15, 0]
                c_bc_rv = [15, 0]
                
            [data_properties]
              max_read_length_bp = 101
    """))
    fixed_cutting_sites = FixedRelativeCuttingSites(config_stub)
    tmpdir = tempfile.mkdtemp()
    fixed_cutting_sites_fp = op.join(tmpdir, 'cutting_sites.p')
    with open(fixed_cutting_sites_fp, 'wb') as fout:
        pickle.dump(fixed_cutting_sites, fout)
    yield fixed_cutting_sites_fp
    shutil.rmtree(tmpdir)

@pytest.fixture(scope='module',
                params=['correct_fixed_cut_sites', 'all_zero_fixed_cut_sites'])
def config_file_path(request, cutting_sites_obj_fp):

    config_file_dir = tempfile.mkdtemp()
    user_config_file_path = op.join(config_file_dir, f"{request.param}.toml")

    config_file_text = dedent(f"""
            [basic_quality_filtering]
              min_mapq = {MIN_MAPQ}
              min_phred_score = {MIN_PHRED}
            
            [data_properties]
              max_read_length_bp = 101
            """)

    if request.param == 'correct_fixed_cut_sites':
        config_file_text += dedent("""
            [trimming]
              [trimming.relative_to_fragment_ends_dict]
                w_bc    = [0, 20]
                c_bc    = [0, 20]
                w_bc_rv = [15, 0]
                c_bc_rv = [15, 0]
                """)
    else:
        config_file_text += dedent(f"""
            [trimming]
              [trimming.relative_to_fragment_ends_dict]
                w_bc    = [0, 0]
                c_bc    = [0, 0]
                w_bc_rv = [0, 0]
                c_bc_rv = [0, 0]
            
            [paths]
              adjusted_cutting_sites_obj_p  = "{cutting_sites_obj_fp}"
              
            """)

    with open(user_config_file_path, "wt") as fobj:
        fobj.write(config_file_text)

    yield user_config_file_path

    shutil.rmtree(config_file_dir)

@pytest.fixture(scope="module")
def base_command_list(config_file_path):

    base_command_list = ['mqc', 'call',
                         '--bam', TEST_BAM,
                         '--config_file', config_file_path,
                         '--sample_name', SAMPLE_NAME,
                         '--sample_meta', SAMPLE_META,
                         '--cores', '2']

    if 'all_zero_fixed_cut_sites' in config_file_path:
        base_command_list += ['--use_mbias_fit']

    return base_command_list


def get_text_of_call_files(output_dir, motifs_str, user_config_file):
    config = assemble_config_vars(
        command_line_args_dict=dict(output_dir=output_dir,
                                    sample_name = SAMPLE_NAME,
                                    sample_meta=SAMPLE_META,
                                    motifs_str=motifs_str),
        default_config_file_path=DEFAULT_CONFIG_FILE,
        user_config_file_path=user_config_file)

    def read_calls(chrom, motif):
        mcall_fp = (config['paths']['meth_calls_basepath']
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
