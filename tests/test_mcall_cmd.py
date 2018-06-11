import gzip
import os.path as op
import re
from itertools import product

import pickle
from pathlib import Path

import pytest
import shutil
import subprocess
import tempfile
from textwrap import dedent

import toml

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
                            #chrom	start	end	motif	score	strand	beta_value	n_meth	n_total
                            1	11298399	11298400	CG	.	+	{1:.6f}	5	5
                            1	11298400	11298401	CG	.	-	nan	0	0
                            1	11299330	11299331	CG	.	+	{1:.6f}	3	3
                            1	11299331	11299332	CG	.	-	{1:.6f}	4	4
                            """),
                  'CHG': dedent(f"""\
                            #chrom	start	end	motif	score	strand	beta_value	n_meth	n_total
                            1	3258371	3258372	CHG	.	+	{0:.6f}	0	3
                            1	3258373	3258374	CHG	.	-	{0:.6f}	0	3
                            """),
                  'CHH': dedent(f"""\
                            #chrom	start	end	motif	score	strand	beta_value	n_meth	n_total
                            1	3258374	3258375	CHH	.	-	{0:.6f}	0	3
                            1	3258376	3258377	CHH	.	+	{0:.6f}	0	3
                            """),
                  },
            '2': {'CG': dedent(f"""\
                            #chrom	start	end	motif	score	strand	beta_value	n_meth	n_total
                            2	9042611	9042612	CG	.	+	{1:.6f}	4	4
                            2	9042612	9042613	CG	.	-	{1:.6f}	1	1
                            2	9042613	9042614	CG	.	+	{1:.6f}	5	5
                            2	9042614	9042615	CG	.	-	{1:.6f}	1	1
                            """),
                  'CHG': dedent(f"""\
                            #chrom	start	end	motif	score	strand	beta_value	n_meth	n_total
                            2	3281431	3281432	CHG	.	-	{0:.6f}	0	3
                            """),
                  'CHH': dedent(f"""\
                            #chrom	start	end	motif	score	strand	beta_value	n_meth	n_total
                            2	3281432	3281433	CHH	.	-	{0:.6f}	0	3
                            2	3281434	3281435	CHH	.	+	{0:.6f}	0	1
                            """),
                  }
            }


def space_to_tab(s):
    return re.sub(' +', '\t', s)

@pytest.fixture(scope="module")
def cutting_sites_obj_fp():
    config_stub = toml.loads(dedent("""
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
                                    sample_name=SAMPLE_NAME,
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


@pytest.mark.parametrize('output_formats', 'bismark bismark,bed bed'.split())
@pytest.mark.parametrize('motifs_str', ['CG', 'CG-CHG-CHH'])
def test_call_tool(tmpdir, config_file_path, index_file_paths_dict,
                   motifs_str, output_formats, make_interactive):

    tmpdir = str(tmpdir)

    command_args = ['mqc', 'call',
                    '--bam', TEST_BAM,
                    '--config_file', config_file_path,
                    '--sample_name', SAMPLE_NAME,
                    '--sample_meta', SAMPLE_META,
                    '--output_formats', output_formats,
                    '--output_dir', tmpdir,
                    '--cores', '2']
    # add index files as positional arguments
    command_args += index_file_paths_dict[motifs_str]

    if Path(config_file_path).stem == 'all_zero_fixed_cut_sites':
        command_args += ['--use_mbias_fit']

    subprocess.run(command_args, check=True)

    if not make_interactive:
        config = assemble_config_vars(
            command_line_args_dict=dict(
                output_dir=tmpdir,
                sample_name=SAMPLE_NAME,
                sample_meta=SAMPLE_META,
                # required because we currently need to expand the entire paths
                # dict, which has {motifs_str} fields
                motifs_str=motifs_str),
            default_config_file_path=DEFAULT_CONFIG_FILE,
            user_config_file_path=config_file_path)
        output_path_templates = dict(
            bed=(config['paths']['meth_calls_basepath'] +
                 f'_{SAMPLE_NAME}_{{motif}}_{{chrom}}.bed.gz'),
            bismark=(config['paths']['bismark_calls_by_chrom_motif']
                     .replace('[', '{').replace(']', '}'))
        )

        all_motifs = 'CG CHG CHH'.split()
        all_output_formats = 'bed bismark'.split()
        all_chroms = '1 2'.split()
        all_possible_file_spec_tuples = set(
            product(all_motifs, all_output_formats, all_chroms))

        motifs = motifs_str.split('-')
        output_formats = output_formats.split(',')
        chroms = ['1', '2']
        spec_tuples_for_expected_output = set(product(
            motifs, output_formats, chroms))

        for motif, output_format, chrom in spec_tuples_for_expected_output:
            assert Path(output_path_templates[output_format].format(motif=motif, chrom=chrom)).exists()

        for motif, output_format, chrom in all_possible_file_spec_tuples - spec_tuples_for_expected_output:
            assert not Path(output_path_templates[output_format].format(motif=motif, chrom=chrom)).exists()

        computed_file_contents = {}
        expected_file_contents = {}
        for motif, output_format, chrom in spec_tuples_for_expected_output:
            computed_file_contents[motif, output_format, chrom] = (
                gzip.open(output_path_templates[output_format].format(motif=motif, chrom=chrom), 'rt').read())
            expected_file_contents[motif, output_format, chrom] = (
                EXPECTED_RESULTS_DICT2[motif, output_format, chrom]
            )
        assert computed_file_contents == expected_file_contents

    else:
        subprocess.run(['firefox', tmpdir])
        ans = input('Everything ok? y/n')
        if ans == 'y':
            assert True
        else:
            assert False

# ============================================================================ #
# Expected methylation calling outputs                                         #
# ============================================================================ #
# Contains correct calls for all combinations of
# - output formats: bed, stratified_bed, bismark  [more to be added, e.g. vcf]
# - motifs: CG, CHG, CHH
# - test data chromosomes: 1, 2

EXPECTED_RESULTS_DICT2 = {

    ('CG', 'bed', '1'): dedent(f"""\
        #chrom	start	end	motif	score	strand	beta_value	n_meth	n_total
        1	11298399	11298400	CG	.	+	{1:.6f}	5	5
        1	11298400	11298401	CG	.	-	nan	0	0
        1	11299330	11299331	CG	.	+	{1:.6f}	3	3
        1	11299331	11299332	CG	.	-	{1:.6f}	4	4
        """),

    ('CG', 'bismark', '1'): space_to_tab(dedent(f"""\
        HWI-ST1153:88:D1E30ACXX:4:1211:7814:62078       +       1       11298399        Z
        HWI-ST1153:88:D1E30ACXX:1:2301:18367:58510      +       1       11298399        Z
        HWI-ST1153:88:D1E30ACXX:4:1102:21273:32039      +       1       11298399        Z
        HWI-ST1153:88:D1E30ACXX:4:1101:5201:8139        +       1       11298399        Z
        HWI-ST1153:88:D1E30ACXX:1:2308:15221:47689      +       1       11298399        Z
        HWI-ST1153:88:D1E30ACXX:4:1216:4865:1375        +       1       11299330        Z
        HWI-ST1153:88:D1E30ACXX:1:1207:1233:41678       +       1       11299330        Z
        HWI-ST1153:88:D1E30ACXX:1:1203:16922:31249      +       1       11299330        Z
        HWI-ST1153:88:D1E30ACXX:4:2315:9833:30242       +       1       11299331        Z
        HWI-ST1153:88:D1E30ACXX:1:2316:10863:87110      +       1       11299331        Z
        HWI-ST1153:88:D1E30ACXX:1:1311:7616:6898        +       1       11299331        Z
        HWI-ST1153:88:D1E30ACXX:1:1313:8267:73811       +       1       11299331        Z
        """)),

    ('CHG', 'bed', '1'): dedent(f"""\
                            #chrom	start	end	motif	score	strand	beta_value	n_meth	n_total
                            1	3258371	3258372	CHG	.	+	{0:.6f}	0	3
                            1	3258373	3258374	CHG	.	-	{0:.6f}	0	3
                            """),

    ('CHG', 'bismark', '1'): space_to_tab(dedent(f"""\
        HWI-ST1153:88:D1E30ACXX:1:2315:8385:43128       -       1       3258371 x
        HWI-ST1153:88:D1E30ACXX:4:2302:5844:61601       -       1       3258371 x
        HWI-ST1153:88:D1E30ACXX:4:1107:1986:53479       -       1       3258371 x
        HWI-ST1153:88:D1E30ACXX:1:2106:16572:67528      -       1       3258373 x
        HWI-ST1153:88:D1E30ACXX:1:2316:19549:96321      -       1       3258373 x
        HWI-ST1153:88:D1E30ACXX:4:1307:17993:25897      -       1       3258373 x
        """)),

    ('CHH', 'bed', '1'): dedent(f"""\
        #chrom	start	end	motif	score	strand	beta_value	n_meth	n_total
        1	3258374	3258375	CHH	.	-	{0:.6f}	0	3
        1	3258376	3258377	CHH	.	+	{0:.6f}	0	3
        """),

    ('CHH', 'bismark', '1'): space_to_tab(dedent(f"""\
        HWI-ST1153:88:D1E30ACXX:1:2106:16572:67528      -       1       3258374 h
        HWI-ST1153:88:D1E30ACXX:1:2316:19549:96321      -       1       3258374 h
        HWI-ST1153:88:D1E30ACXX:4:1307:17993:25897      -       1       3258374 h
        HWI-ST1153:88:D1E30ACXX:1:2315:8385:43128       -       1       3258376 h
        HWI-ST1153:88:D1E30ACXX:4:2302:5844:61601       -       1       3258376 h
        HWI-ST1153:88:D1E30ACXX:4:1107:1986:53479       -       1       3258376 h
        """)),

    ('CG', 'bed', '2'): dedent(f"""\
         #chrom	start	end	motif	score	strand	beta_value	n_meth	n_total
         2	9042611	9042612	CG	.	+	{1:.6f}	4	4
         2	9042612	9042613	CG	.	-	{1:.6f}	1	1
         2	9042613	9042614	CG	.	+	{1:.6f}	5	5
         2	9042614	9042615	CG	.	-	{1:.6f}	1	1
         """),

    ('CG', 'bismark', '2'): space_to_tab(dedent(f"""\
        HWI-ST1153:88:D1E30ACXX:4:1105:17809:24267      +       2       9042611 Z
        HWI-ST1153:88:D1E30ACXX:4:2306:3628:53316       +       2       9042611 Z
        HWI-ST1153:88:D1E30ACXX:1:1314:16939:69498      +       2       9042611 Z
        HWI-ST1153:88:D1E30ACXX:1:1201:10436:81016      +       2       9042611 Z
        HWI-ST1153:88:D1E30ACXX:4:2307:6311:10125       +       2       9042612 Z
        HWI-ST1153:88:D1E30ACXX:4:1105:17809:24267      +       2       9042613 Z
        HWI-ST1153:88:D1E30ACXX:1:2311:15387:62265      +       2       9042613 Z
        HWI-ST1153:88:D1E30ACXX:4:2306:3628:53316       +       2       9042613 Z
        HWI-ST1153:88:D1E30ACXX:1:1314:16939:69498      +       2       9042613 Z
        HWI-ST1153:88:D1E30ACXX:1:1303:2047:52877       +       2       9042613 Z
        HWI-ST1153:88:D1E30ACXX:4:2307:6311:10125       +       2       9042614 Z
        """)),

    ('CHG', 'bed', '2'): dedent(f"""\
        #chrom	start	end	motif	score	strand	beta_value	n_meth	n_total
        2	3281431	3281432	CHG	.	-	{0:.6f}	0	3
        """),

    ('CHG', 'bismark', '2'): space_to_tab(dedent(f"""\
        HWI-ST1153:88:D1E30ACXX:1:2115:1487:37860       -       2       3281431 x
        HWI-ST1153:88:D1E30ACXX:4:1302:7825:95345       -       2       3281431 x
        HWI-ST1153:88:D1E30ACXX:4:1315:8590:5265        -       2       3281431 x
        """)),


    ('CHH', 'bed', '2'): dedent(f"""\
        #chrom	start	end	motif	score	strand	beta_value	n_meth	n_total
        2	3281432	3281433	CHH	.	-	{0:.6f}	0	3
        2	3281434	3281435	CHH	.	+	{0:.6f}	0	1
        """),

    ('CHH', 'bismark', '2'): space_to_tab(dedent(f"""\
        HWI-ST1153:88:D1E30ACXX:1:2115:1487:37860       -       2       3281432 h
        HWI-ST1153:88:D1E30ACXX:4:1302:7825:95345       -       2       3281432 h
        HWI-ST1153:88:D1E30ACXX:4:1315:8590:5265        -       2       3281432 h
        HWI-ST1153:88:D1E30ACXX:1:1114:17740:51566      -       2       3281434 h
        """)),
}














































