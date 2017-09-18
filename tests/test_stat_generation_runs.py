import os.path as op
import pytest
import pytoml
import shutil
import subprocess
import tempfile
import warnings

from textwrap import dedent

import pandas as pd

import mqc.flag_and_index_values as mfl
from mqc.config import assemble_config_vars
from mqc.utils import get_resource_abspath

b_inds = mfl.bsseq_strand_indices

TESTS_DIR = op.dirname(__file__)
SAMPLE_NAME = 'hsc_rep1'
SAMPLE_META = 'population=hsc,rep=1,model=blk6'
DEFAULT_CONFIG_FILE = get_resource_abspath('config.default.toml')
columns = ['motif', 'bs_strand', 'flen', 'pos', 'meth_status', 'counts']
index_cols = columns[0:-1]
USER_FLEN_MAX = 200


@pytest.fixture(scope='module',
                params=[True, False],
                ids=['with_user_config', 'wo_user_config'])
def user_config_file(request):
    if request.param:
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
    else:
        yield ''


@pytest.fixture(scope='module')
def computed_mbias_counts_df(user_config_file):
    """Runs mqc stats and provides path to results dir"""

    tmpdir = tempfile.mkdtemp()
    command_args = [
        'mqc', 'stats',
        '--bam', op.join(TESTS_DIR, 'test_files', 'test_mbias-stats.bam'),
        '--output_dir', tmpdir,
        '--sample_name', SAMPLE_NAME,
        '--sample_meta', SAMPLE_META,
        '--cores', '1',
        '--motifs', 'CG,CHG',
        op.join(TESTS_DIR, 'test_files', 'test_index_mbias_chr11.bed.gz'),
        op.join(TESTS_DIR, 'test_files', 'test_index_mbias_chr11.bed.gz'),
    ]

    if user_config_file:
        command_args += ['--config_file', user_config_file]

    subprocess.check_call(command_args)

    config = assemble_config_vars(
        command_line_args_dict=dict(output_dir=tmpdir,
                                    sample_name=SAMPLE_NAME,
                                    sample_meta=SAMPLE_META,
                                    motifs_str='CG-CHG'),
        default_config_file_path=DEFAULT_CONFIG_FILE,
        user_config_file_path=user_config_file)

    mbias_stats_p_path = config["paths"]["mbias_counts"] + ".p"
    df = pd.read_pickle(mbias_stats_p_path)
    df.columns.name = ('with_user_config'
                       if user_config_file
                       else 'no_user_config')
    yield df
    shutil.rmtree(tmpdir)


expected_counts_df = (pd.DataFrame(
    [('CHG', 'w_bc', 116, 24, 'n_unmeth', 2),
     ('CHG', 'w_bc_rv', 116, 93, 'n_unmeth', 2),
     ('CG', 'c_bc_rv', 146, 81, 'n_meth', 2),
     ('CG', 'w_bc', 116, 25, 'n_meth', 2),
     ('CG', 'w_bc_rv', 116, 92, 'n_meth', 2),
     ('CG', 'c_bc', 146, 66, 'n_meth', 2),
     ('CG', 'w_bc', 116, 64, 'n_unmeth', 2),
     ('CG', 'w_bc_rv', 116, 53, 'n_meth', 2),
     ('CG', 'c_bc', 146, 27, 'n_meth', 2), ],
    columns=columns)
                      .set_index(columns[:-1])
                      .sort_index())


def test_stats_run_through(computed_mbias_counts_df):
    print(computed_mbias_counts_df.head())
    # assert type(computed_mbias_counts_df) == pd.DataFrame
    assert False



# @pytest.mark.acceptance_test
# @pytest.mark.parametrize("stratum_idx,values",
#                          expected_counts_df.iterrows())
# def test_stats_generation_on_small_bam_with_undefined_properties(stratum_idx,
#                                                                  values,
#                                                                  mbias_stats_df):
#     assert mbias_stats_df.loc[stratum_idx, 'counts'] == values['counts'], repr(
#         stratum_idx)
#
#
# @pytest.mark.acceptance_test
# def test_no_counts_in_strata_not_hit_by_bam(mbias_stats_df):
#     mbias_stats_df.loc[expected_counts_df.index, 'counts'] = 0
#     false_counts = str(mbias_stats_df.loc[mbias_stats_df['counts'] > 0, :])
#     message = "The following events are unexpected counts:\n" + false_counts
#     assert (mbias_stats_df['counts'] == 0).all(), message
#
#
# @pytest.mark.acceptance_test
# def test_provides_strat_mbias_counts_as_pickle_and_tsv(
#         mbias_stats_df, filled_mbias_stats_results_dir, default_paths):
#     # TODO: FutureWarning: elementwise comparison failed; returning scalar instead, but in the future will perform elementwise comparison mask |= (ar1 == a)
#     mbias_stats_tsv = default_paths['mbias_counts'] + '.tsv'
#     mbias_stats_df_from_tsv = pd.read_csv(mbias_stats_tsv,
#                                           sep='\t', header=0,
#                                           dtype=None,
#                                           index_col=index_cols,
#                                           nrows=None, usecols=None,
#                                           comment=None)
#     msg = f"""
#
#     Pickled dataframe
#     -----------------
#     {mbias_stats_df.head()}
#
#     Read from TSV
#     -------------
#     {mbias_stats_df_from_tsv.head()}
#
#     """
#
#     # TODO: does not catch warning
#     with warnings.catch_warnings():
#         assert mbias_stats_df.equals(mbias_stats_df_from_tsv), msg
#
#
# @pytest.mark.acceptance_test
# def test_user_config_is_used(mbias_stats_df: pd.DataFrame):
#     computed_max_flen = mbias_stats_df.index.get_level_values('flen').max()
#
#     if mbias_stats_df.columns.name == 'with_user_config':
#         expected_max_flen = USER_FLEN_MAX
#     else:  # no user config
#         with open(DEFAULT_CONFIG_FILE) as f:
#             default_config_file_dict = pytoml.load(f)
#         expected_max_flen = default_config_file_dict['trimming'][
#             'max_flen_considered_for_trimming']
#
#     assert expected_max_flen == computed_max_flen
