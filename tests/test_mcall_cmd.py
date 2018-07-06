"""Test the call tool"""

import gzip
import json
import os.path as op
import re
import shutil
import subprocess
import tempfile
from itertools import product
from pathlib import Path
from textwrap import dedent
from typing import Dict, List, Any

import pytest

from mqc.config import assemble_config_vars
from mqc.mbias import CuttingSitesReplacementClass
from mqc.utils import get_resource_abspath

TESTS_DIR = op.dirname(__file__)
TEST_BAM = op.join(TESTS_DIR, 'test_files', 'test_mcall.sorted.bam')
DEFAULT_CONFIG_FILE = get_resource_abspath('config.default.toml')
SAMPLE_NAME = 'hsc_rep1'
SAMPLE_META = 'population=hsc,rep=1,model=blk6'
TEST_FILES_DIR = op.join(TESTS_DIR, 'test_files')
MIN_MAPQ = 41
MIN_PHRED = 35


@pytest.fixture(scope='module')
def index_file_paths_dict_fix() -> Dict[str, List[str]]:
    """Provide index files by motif"""
    return {
        'CG': [op.join(TEST_FILES_DIR, 'test_mcall_CG_1.bed.gz'),
               op.join(TEST_FILES_DIR, 'test_mcall_CG_2.bed.gz')],
        'CG-CHG-CHH': [op.join(TEST_FILES_DIR, 'test_mcall_CG-CHG-CHH_1.bed.gz'),
                       op.join(TEST_FILES_DIR, 'test_mcall_CG-CHG-CHH_2.bed.gz')]
    }


def space_to_tab(s: str) -> str:
    return re.sub(' +', '\t', s)


@pytest.fixture(scope="module")
def cutting_sites_obj_fp_fix():
    """Provide FixedRelativeCuttingSites as fixture"""
    cutting_sites = CuttingSitesReplacementClass.from_rel_to_frag_end_cutting_sites(
        cut_site_spec=dict(w_bc=[0, 20], c_bc=[0, 20],
                           w_bc_rv=[15, 0], c_bc_rv=[15, 0]),
        max_read_length=101)
    tmpdir = tempfile.mkdtemp()
    cutting_sites_df_fp = op.join(tmpdir, 'cutting_sites_df.p')
    cutting_sites.df.to_pickle(cutting_sites_df_fp)
    yield cutting_sites_df_fp
    shutil.rmtree(tmpdir)


# noinspection PyIncorrectDocstring,PyShadowingNames
@pytest.fixture(scope='module',
                params=['correct_fixed_cut_sites', 'all_zero_fixed_cut_sites'])
def config_file_path(request):
    """Provide custom config files for testing fixed and adjusted cutting sites

    Provides two variants of a user-supplied config
    file for the acceptance tests of the call tool.

    The config file 'correct_fixed_cut_sites.toml' provides correct
    definitions of the fixed cutting sites that were used for
    generation of the expected results. The other config file
    variant, 'all_zero_fixed_cut_sites', deliberately provides
    wrong cutting sites (no cutting at all). This variant is
    intended to be used together with a cutting sites object,
    as would be created from the M-bias stats analysis. If the
    provided cutting sites object (the config file contains the
    correct path) is not used appropriately, there is no correct
    fallback option in the config file and the test will fail.

    Returns:
        Paths to two config files:
            '{tmpdir}/correct_fixed_cut_sites.toml',
            '{tmpdir}/all_zero_fixed_cut_sites.toml'.

    """

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
              
            """)

    with open(user_config_file_path, "wt") as fobj:
        fobj.write(config_file_text)

    yield user_config_file_path

    shutil.rmtree(config_file_dir)


# noinspection PyShadowingNames
@pytest.mark.acceptance_test
@pytest.mark.parametrize('output_formats_csv', 'bismark bismark,bed bed stratified_bed bismark,stratified_bed'.split())
@pytest.mark.parametrize('motifs_str', ['CG', 'CG-CHG-CHH'])
def test_call_tool(tmpdir: Any, config_file_path: str,
                   index_file_paths_dict_fix: Dict[str, List[str]],
                   motifs_str: str, output_formats_csv: str,
                   make_interactive: bool, cutting_sites_obj_fp_fix: str) -> None:
    """Test call tool across various combinations of options

    Args:
        config_file_path: fixture, providing paths to custom config
            files. The config files are used to support testing the
            trimming functionality. See fixture doc for details.
        index_file_paths_dict_fix: CG or CG-CHG-CHH (combined) index
            files for chromosomes 1 and 2.
    """

    # When new output formats are added:
    # - add to output_formats parametrization (above)
    # - add to all_output_formats variable (below)
    # - add to output_path_templates (below)

    tmpdir = str(tmpdir)

    pos_to_remove_from_frag_end = dict(
        w_bc=[0, 20],
        c_bc=[0, 20],
        w_bc_rv=[15, 0],
        c_bc_rv=[15, 0],
    )

    command_args = ['mqc', 'call',
                    '--bam', TEST_BAM,
                    '--config_file', config_file_path,
                    '--sample_name', SAMPLE_NAME,
                    '--sample_meta', SAMPLE_META,
                    '--output_formats', output_formats_csv,
                    '--output_dir', tmpdir,
                    '--trimming', f"frag_end::{json.dumps(pos_to_remove_from_frag_end)}",
                    '--max_read_length', '101',
                    '--cores', '2']
    # add index files as positional arguments
    command_args += index_file_paths_dict_fix[motifs_str]

    # If the config file with deliberately wrong config of the fixed
    # cutting sites is used, we must use the cutting sites object
    # defined in this config file
    if Path(config_file_path).stem == 'all_zero_fixed_cut_sites':
        command_args += ['--trimming', f'cutting_sites::{cutting_sites_obj_fp_fix}']

    subprocess.run(command_args, check=True)

    # interactive testing allows inspecting results in the browser
    # which may be helpful during development
    # for automated testing, don't set the make_interactive flag of
    # pytest!
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

        # The mcall files will be present per motif and chrom
        # we always process chromosomes 1 and 2 on the test BAM file
        output_path_templates = dict(
            bed=(config['paths']['meth_calls_basepath'] +
                 f'_{SAMPLE_NAME}_{{motif}}_{{chrom}}.bed.gz'),
            stratified_bed=(config['paths']['strat_bed_calls_by_chrom_motif']
                         .replace('[', '{').replace(']', '}')),
            bismark=(config['paths']['bismark_calls_by_chrom_motif']
                     .replace('[', '{').replace(']', '}'))
        )

        # Check that all expected files are there, but no unexpected
        # files were produced
        # --------------------------------------------------------------
        # To find out which mcall files should not be there:
        # 1. Compute all files which may possibly be created
        # 2. Compute files which should be created
        # 3. Find the complement 1 - 2
        all_motifs = 'CG CHG CHH'.split()
        all_output_formats = 'bed bismark stratified_bed'.split()
        all_chroms = '1 2'.split()
        all_possible_file_spec_tuples = set(
            product(all_motifs, all_output_formats, all_chroms))

        motifs = motifs_str.split('-')
        output_formats = output_formats_csv.split(',')
        chroms = ['1', '2']
        spec_tuples_for_expected_output = set(product(
            motifs, output_formats, chroms))

        for motif, output_format, chrom in spec_tuples_for_expected_output:
            assert (Path(output_path_templates[output_format]
                         .format(motif=motif, chrom=chrom)).exists())

        for motif, output_format, chrom in (
                all_possible_file_spec_tuples - spec_tuples_for_expected_output):
            assert not (Path(output_path_templates[output_format]
                             .format(motif=motif, chrom=chrom)).exists())

        # Check that file contents are correct
        # --------------------------------------------------------------
        # All output files are collected in a flat dict (subset of the
        # elements in EXPECTED_RESULTS_DICT2) and then compared at once
        # This way, the diff provided by pytest shows all failures and
        # does not abort after the first faulty file. Achieving this
        # through a further parametrization of this test function
        # would be quite complicated - this is the better trade off.

        computed_file_contents = {}
        expected_file_contents = {}
        for motif, output_format, chrom in spec_tuples_for_expected_output:
            computed_file_contents[motif, output_format, chrom] = (
                gzip.open(output_path_templates[output_format]
                          .format(motif=motif, chrom=chrom), 'rt').read())
            expected_file_contents[motif, output_format, chrom] = (
                EXPECTED_RESULTS_DICT2[motif, output_format, chrom]
            )
        assert computed_file_contents == expected_file_contents, \
            computed_file_contents

    else:
        # interactive test to be used during development only
        # run pytest with --make_interactive flag
        subprocess.run(['firefox', tmpdir])
        ans = input('Everything ok? y/n')
        if ans == 'y':
            assert True
        else:
            assert False

# ============================================================================ #
# Expected methylation calling outputs                                         #
# ============================================================================ #
#
# Contains correct calls for all combinations of
# - output formats: bed, stratified_bed, bismark [more to come, e.g. vcf]
# - motifs: CG, CHG, CHH
# - test data chromosomes: 1, 2
#
# The correct file contents were created manually, based on the
# cutting sites defined in the config file. The same cutting sites
# are defined in both the fixed cutting sites section of the 'correct
# config file' and in the CuttingSites object provided through
# the 'faulty config file'. See config file fixture for details.


EXPECTED_RESULTS_DICT2 = {

    ('CG', 'bed', '1'): dedent(f"""\
        #chrom	start	end	motif	score	strand	beta_value	n_meth	n_total
        1	11298399	11298400	CG	.	+	{1:.8f}	5	5
        1	11298400	11298401	CG	.	-	nan	0	0
        1	11299330	11299331	CG	.	+	{1:.8f}	3	3
        1	11299331	11299332	CG	.	-	{1:.8f}	4	4
        """),

    # verified manually
    ('CG', 'stratified_bed', '1'): space_to_tab(dedent(f"""\
        #chrom	start	end	motif	score	strand	beta_value	n_meth	n_total	c_bc_beta_value	c_bc_n_meth	c_bc_n_total	c_bc_rv_beta_value	c_bc_rv_n_meth	c_bc_rv_n_total	w_bc_beta_value	w_bc_n_meth	w_bc_n_total	w_bc_rv_beta_value	w_bc_rv_n_meth	w_bc_rv_n_total	mate1_beta_value	mate1_n_meth	mate1_n_total	mate2_beta_value	mate2_n_meth	mate2_n_total
        1	11298399	11298400	CG	.	+	1.00000000	5	5	nan	0	0	nan	0	0	1.00000000	3	3	1.00000000	2	2	1.00000000	3	3	1.00000000	2	2
        1	11298400	11298401	CG	.	-	nan	0	0	nan	0	0	nan	0	0	nan	0	0	nan	0	0	nan	0	0	nan	0	0
        1	11299330	11299331	CG	.	+	1.00000000	3	3	nan	0	0	nan	0	0	1.00000000	2	2	1.00000000	2	2	1.00000000	2	2	1.00000000	2	2
        1	11299331	11299332	CG	.	-	1.00000000	4	4	1.00000000	4	4	nan	0	0	nan	0	0	nan	0	0	1.00000000	4	4	nan	0	0
        """)),

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
                            1	3258371	3258372	CHG	.	+	{0:.8f}	0	3
                            1	3258373	3258374	CHG	.	-	{0:.8f}	0	3
                            """),

    ('CHG', 'stratified_bed', '1'): space_to_tab(dedent(f"""\
        #chrom	start	end	motif	score	strand	beta_value	n_meth	n_total	c_bc_beta_value	c_bc_n_meth	c_bc_n_total	c_bc_rv_beta_value	c_bc_rv_n_meth	c_bc_rv_n_total	w_bc_beta_value	w_bc_n_meth	w_bc_n_total	w_bc_rv_beta_value	w_bc_rv_n_meth	w_bc_rv_n_total	mate1_beta_value	mate1_n_meth	mate1_n_total	mate2_beta_value	mate2_n_meth	mate2_n_total
        1	3258371	3258372	CHG	.	+	0.00000000	0	3	nan	0	0	nan	0	0	0.00000000	0	1	0.00000000	0	2	0.00000000	0	1	0.00000000	0	2
        1	3258373	3258374	CHG	.	-	0.00000000	0	3	0.00000000	0	2	0.00000000	0	1	nan	0	0	nan	0	0	0.00000000	0	2	0.00000000	0	1
        """)),

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
        1	3258374	3258375	CHH	.	-	{0:.8f}	0	3
        1	3258376	3258377	CHH	.	+	{0:.8f}	0	3
        """),

    ('CHH', 'stratified_bed', '1'): space_to_tab(dedent("""\
        #chrom	start	end	motif	score	strand	beta_value	n_meth	n_total	c_bc_beta_value	c_bc_n_meth	c_bc_n_total	c_bc_rv_beta_value	c_bc_rv_n_meth	c_bc_rv_n_total	w_bc_beta_value	w_bc_n_meth	w_bc_n_total	w_bc_rv_beta_value	w_bc_rv_n_meth	w_bc_rv_n_total	mate1_beta_value	mate1_n_meth	mate1_n_total	mate2_beta_value	mate2_n_meth	mate2_n_total
        1	3258374	3258375	CHH	.	-	0.00000000	0	3	0.00000000	0	2	0.00000000	0	1	nan	0	0	nan	0	0	0.00000000	0	2	0.00000000	0	1
        1	3258376	3258377	CHH	.	+	0.00000000	0	3	nan	0	0	nan	0	0	0.00000000	0	1	0.00000000	0	3	0.00000000	0	1	0.00000000	0	3
        """)),


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
        2	9042611	9042612	CG	.	+	{1:.8f}	4	4
        2	9042612	9042613	CG	.	-	{1:.8f}	1	1
        2	9042613	9042614	CG	.	+	{1:.8f}	5	5
        2	9042614	9042615	CG	.	-	{1:.8f}	1	1
        """),

    # verified manually
    ('CG', 'stratified_bed', '2'): space_to_tab(dedent(f"""\
        #chrom	start	end	motif	score	strand	beta_value	n_meth	n_total	c_bc_beta_value	c_bc_n_meth	c_bc_n_total	c_bc_rv_beta_value	c_bc_rv_n_meth	c_bc_rv_n_total	w_bc_beta_value	w_bc_n_meth	w_bc_n_total	w_bc_rv_beta_value	w_bc_rv_n_meth	w_bc_rv_n_total	mate1_beta_value	mate1_n_meth	mate1_n_total	mate2_beta_value	mate2_n_meth	mate2_n_total
        2	9042611	9042612	CG	.	+	1.00000000	4	4	nan	0	0	nan	0	0	1.00000000	1	1	1.00000000	4	4	1.00000000	1	1	1.00000000	4	4
        2	9042612	9042613	CG	.	-	1.00000000	1	1	1.00000000	1	1	nan	0	0	nan	0	0	nan	0	0	1.00000000	1	1	nan	0	0
        2	9042613	9042614	CG	.	+	1.00000000	5	5	nan	0	0	nan	0	0	1.00000000	1	1	1.00000000	5	5	1.00000000	1	1	1.00000000	5	5
        2	9042614	9042615	CG	.	-	1.00000000	1	1	1.00000000	1	1	nan	0	0	nan	0	0	nan	0	0	1.00000000	1	1	nan	0	0
        """)),

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
        2	3281431	3281432	CHG	.	-	{0:.8f}	0	3
        """),

    ('CHG', 'stratified_bed', '2'): space_to_tab(dedent(f"""\
        #chrom	start	end	motif	score	strand	beta_value	n_meth	n_total	c_bc_beta_value	c_bc_n_meth	c_bc_n_total	c_bc_rv_beta_value	c_bc_rv_n_meth	c_bc_rv_n_total	w_bc_beta_value	w_bc_n_meth	w_bc_n_total	w_bc_rv_beta_value	w_bc_rv_n_meth	w_bc_rv_n_total	mate1_beta_value	mate1_n_meth	mate1_n_total	mate2_beta_value	mate2_n_meth	mate2_n_total
        2	3281431	3281432	CHG	.	-	0.00000000	0	3	0.00000000	0	1	0.00000000	0	2	nan	0	0	nan	0	0	0.00000000	0	1	0.00000000	0	2
        """)),

    ('CHG', 'bismark', '2'): space_to_tab(dedent(f"""\
        HWI-ST1153:88:D1E30ACXX:1:2115:1487:37860       -       2       3281431 x
        HWI-ST1153:88:D1E30ACXX:4:1302:7825:95345       -       2       3281431 x
        HWI-ST1153:88:D1E30ACXX:4:1315:8590:5265        -       2       3281431 x
        """)),


    ('CHH', 'bed', '2'): dedent(f"""\
        #chrom	start	end	motif	score	strand	beta_value	n_meth	n_total
        2	3281432	3281433	CHH	.	-	{0:.8f}	0	3
        2	3281434	3281435	CHH	.	+	{0:.8f}	0	1
        """),

    ('CHH', 'stratified_bed', '2'): space_to_tab(dedent(f"""\
        #chrom	start	end	motif	score	strand	beta_value	n_meth	n_total	c_bc_beta_value	c_bc_n_meth	c_bc_n_total	c_bc_rv_beta_value	c_bc_rv_n_meth	c_bc_rv_n_total	w_bc_beta_value	w_bc_n_meth	w_bc_n_total	w_bc_rv_beta_value	w_bc_rv_n_meth	w_bc_rv_n_total	mate1_beta_value	mate1_n_meth	mate1_n_total	mate2_beta_value	mate2_n_meth	mate2_n_total
        2	3281432	3281433	CHH	.	-	0.00000000	0	3	0.00000000	0	1	0.00000000	0	2	nan	0	0	nan	0	0	0.00000000	0	1	0.00000000	0	2
        2	3281434	3281435	CHH	.	+	0.00000000	0	1	nan	0	0	nan	0	0	0.00000000	0	1	nan	0	0	0.00000000	0	1	nan	0	0
        """)),

    ('CHH', 'bismark', '2'): space_to_tab(dedent(f"""\
        HWI-ST1153:88:D1E30ACXX:1:2115:1487:37860       -       2       3281432 h
        HWI-ST1153:88:D1E30ACXX:4:1302:7825:95345       -       2       3281432 h
        HWI-ST1153:88:D1E30ACXX:4:1315:8590:5265        -       2       3281432 h
        HWI-ST1153:88:D1E30ACXX:1:1114:17740:51566      -       2       3281434 h
        """)),
}
