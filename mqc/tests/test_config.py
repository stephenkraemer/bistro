import os
import pytest
import shutil
import tempfile
import textwrap

from mqc.config import assemble_config_vars

POPULATION = "HSC"
OUTPUT_DIR = "/abs/path/to/out_dir"
SAMPLE_NAME = "hsc_rep1"


@pytest.fixture(scope='module')
def command_line_args_dict():
    """This dict is provided by click parsing the command line args"""
    return {'sample_name': SAMPLE_NAME,
            'sample_meta': f"population={POPULATION},replicate=1",
            'output_dir': OUTPUT_DIR,
            'run_parameter1': "value1"}


@pytest.fixture(scope='module')
def user_config_file():
    tmp_dir_path = tempfile.mkdtemp()
    tmp_file_path = os.path.join(tmp_dir_path, 'user_config_file.toml')
    with open(tmp_file_path, 'wt') as fobj:
        fobj.write(textwrap.dedent("""\
            [config_section1]
            value = "custom_value"
                [config_section1.nested_section]
                value = "custom_value"
                [config_section1.nested_section.nested_section]
                value = "custom_value"
            
            [paths]
            dir1 = "/path/to/dir1"
            dir2 = "{dir1}/dir2"
            dir1b = "path/to/dir1b"
            dir2b = "{dir1b}/dir2b"
            results_dir = "relpath/to/results_dir"
            coverage_dir = "{results_dir}/coverage_dir"
            coverage_counts_p = "{coverage_dir}/{name}_{population}.cov.counts.p"
            """))
    yield tmp_file_path
    shutil.rmtree(tmp_dir_path)


@pytest.fixture(scope='module')
def default_config_file():
    tmp_dir_path = tempfile.mkdtemp()
    tmp_file_path = os.path.join(tmp_dir_path, 'user_config_file.toml')
    with open(tmp_file_path, 'wt') as fobj:
        fobj.write(textwrap.dedent("""\
            [config_section1]
            value = "default_value"
            value2 = "default_value2"
                [config_section1.nested_section]
                value = "default_value"
                value2 = "default_value2"
                    [config_section1.nested_section.nested_section]
                    value = "default_value"
                    value2 = "default_value2"
            
            [config_section2]
            value = "default_value"
            
            [paths]
            dir1 = "default"
            dir2 = "default"
            dir1b = "default"
            dir2b = "default"
            results_dir = "default"
            coverage_dir = "default"
            coverage_counts_p = "default"
            """))
    yield tmp_file_path
    shutil.rmtree(tmp_dir_path)


@pytest.fixture()
def config(default_config_file,
           user_config_file,
           command_line_args_dict):
    config = assemble_config_vars(command_line_args_dict,
                                  default_config_file_path=default_config_file,
                                  user_config_file_path=user_config_file, )
    return config


class TestAssembleConfigVarsFn:
    def test_default_to_base_config_file(self, config):
        assert config['config_section1']['value2'] == "default_value2"
        assert (config['config_section1']
                ['nested_section']['nested_section'][
                    'value2'] == "default_value2")
        assert config['config_section2']['value'] == "default_value"

    def test_overwrite_base_config_file_defaults_with_user_config_file(self, config):
        # Check for user config file precedence in top level variables
        assert config['config_section1']['value'] == "custom_value"
        # Check for user config file precedence in nested sections
        assert (config['config_section1']
                ['nested_section']['nested_section'][
                    'value'] == "custom_value")

    def test_do_not_allow_keys_in_user_config_file_which_are_not_in_base_config_file(
            self, default_config_file, user_config_file, command_line_args_dict):
        illegal_config_text = textwrap.dedent("""\
        illegal_key = "illegal_value"
        
        """)
        tmp_dir_path, illegal_user_config_file = add_to_config_file(
            user_config_file, illegal_config_text
        )
        with pytest.raises(KeyError) as excinfo:
            assemble_config_vars(command_line_args_dict,
                                 default_config_file_path=default_config_file,
                                 user_config_file_path=illegal_user_config_file)
        assert ("Config file contains unknown element 'illegal_key'"
                in str(excinfo.value))
        shutil.rmtree(tmp_dir_path)

    def store_sample_metadata_passed_through_cli_in_own_config_section(self, config):
        assert config['sample']['population'] == POPULATION
        assert config['sample']['name'] == SAMPLE_NAME

    def test_store_run_params_passed_through_cli_in_own_config_section(
            self, config, command_line_args_dict):
        assert config['run']['output_dir'] == OUTPUT_DIR
        assert config['run']['run_parameter1'] == "value1"

    def test_interpret_relpaths_as_rel_to_output_dir(
            self, config):
        """Only for config variables in 'paths' section"""
        # absolute path is left as is
        assert config['paths']['dir2'] == "/path/to/dir1/dir2"
        # relpath is relative to output dir
        assert config['paths']['dir2b'] == f"{OUTPUT_DIR}/path/to/dir1b/dir2b"

    def test_expand_fields_referring_to_sample_metadata_and_other_paths(
            self, config):
        """Only for config variables in 'paths' section

        This test demonstrates the key features of path expansion:

        - relpaths are considered relative to the output dir and are expanded to abspaths
        - braced fields can access:

            - other variables defined in the path section
            - sample metadata

        - field lookup may be recursive

        """
        expected_path_under_coverage_counts_p_key = (
            '{output_dir}/relpath/to/results_dir/coverage_dir/'
            '{name}_{population}.cov.counts.p'.format(
                output_dir=OUTPUT_DIR,
                population=POPULATION,
                name=SAMPLE_NAME))
        assert (config['paths']['coverage_counts_p'] ==
                expected_path_under_coverage_counts_p_key)

    def test_raises_when_base_config_file_contains_sample_or_run_section(
        self, default_config_file, user_config_file, command_line_args_dict):
        illegal_section = textwrap.dedent("""\
            [run]
            run_param1 = "abc"
            
            """)
        tmp_dir_path, illegal_default_config_file = add_to_config_file(
            default_config_file, illegal_section)
        with pytest.raises(KeyError) as excinfo:
            assemble_config_vars(
                command_line_args_dict,
                default_config_file_path=illegal_default_config_file,
                user_config_file_path=user_config_file)
        expected_err_msg = ("Config file may not contain reserved "
                            "sections 'run' and 'sample'")
        assert (expected_err_msg in str(excinfo.value))
        shutil.rmtree(tmp_dir_path)

    def test_raises_when_default_and_custom_types_do_not_match(
            self, default_config_file, user_config_file, command_line_args_dict):

        default_section_to_append = textwrap.dedent("""\
        
        [appended_section]
        param1 = "value1"
        
        """)

        custom_value_to_prepend = textwrap.dedent("""\
        
        appended_section = "value1"
        
        """)

        tmp_dir_path1, new_default_config_file = add_to_config_file(
            default_config_file, default_section_to_append)
        tmp_dir_path2, new_user_config_file = add_to_config_file(
            user_config_file, custom_value_to_prepend, prepend=True)

        with pytest.raises(TypeError):
            assemble_config_vars(
                command_line_args_dict,
                default_config_file_path=new_default_config_file,
                user_config_file_path=new_user_config_file)

        shutil.rmtree(tmp_dir_path1)
        shutil.rmtree(tmp_dir_path2)

def add_to_config_file(config_file_path, text_to_add, prepend=False):

    with open(config_file_path) as fobj:
        orig_config_file_text = fobj.read()

    tmp_dir_path = tempfile.mkdtemp()
    new_config_file = os.path.join(tmp_dir_path,
                                   'new_config_file.toml')

    with open(new_config_file, 'wt') as fobj:
        if prepend:
            fobj.write(text_to_add)
        fobj.write(orig_config_file_text + '\n')
        if not prepend:
            fobj.write(text_to_add)

    return tmp_dir_path, new_config_file
