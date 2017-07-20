"""Config file and command line parsing

WGBS data processing tasks are usually run over a cohort of samples.
Thereby, there a configuration values which

1. Change from sample to sample within the same project
2. change from project to project, but remain the same between samples of the same project
3. Can often be left at defaults across projects, but sometimes have to be adapted

This is handled as follows

1. Sample-specific config variables are passed through the CLI
2. Project-specific config variables are defined in a TOML config file provided by the user
3. Nothing is hardcoded, all parameters of the program from file paths to plotting parameters are accessible through a single global config file providing sensible defaults whereever possible

This function parses and integrates all these configuration levels and provides
them in one dict for use in the program. The dict always
contains the sections 'sample' and 'run' with sample metadata and run
parameters passed through the CLI, and additionally all sections specified
in the default config file.

Sample metadata and run parameter specification
-----------------------------------------------

This function expects a CLI parameter named sample_name. In addition,
the run parameter 'output_dir' should be specified. Optionally, a string
specifying further sample metadata can be given as CL arg, it will be parsed
assuming this format: key=value[, key=value...]. All metadata are made
available in the section config[ 'sample'], e.g.

| config['sample']['name'] = 'sample_name'
| config[ 'sample']['metadata_key'] = 'metadata_value'

All other parameters passed through the CLI are collected in config['run']

Path specification and expansion
--------------------------------

All relative paths are assumed to be relative to output dir, which is
assumed to differ between samples and must be given through the CLI.
Therefore, all relative paths will be prepended with the path to
output_dir. The path to output_dir must be absolute.

Moreover, this function expands fields ({}) in path config variables.
Fields may point to either other paths or sample metadata (e.g. {name}).
Currently, circular references are not catched, so be careful. Note that
it is not possible to protect fields using double braces, as one would
do in python format strings. It is also possible to have braces in the
output_dir variable.

"""

import copy
import os
import pytoml
import re
from collections import OrderedDict


def assemble_config_vars(command_line_args_dict,
                         default_config_file_path,
                         user_config_file_path=None):
    """ Construct dict of all configuration parameters

    Parameters
    ----------
    command_line_args_dict: dict
     command line params
    default_config_file_path: str
    user_config_file_path : str

    Returns
    -------
    dict
        dict integrating all configuration variables
    """

    # Pop sample metadata from command line args dict
    # The remaining command line args are run parameters
    command_line_args_dict = copy.deepcopy(command_line_args_dict)
    sample_name = command_line_args_dict.pop('sample_name')
    sample_meta_str = command_line_args_dict.pop('sample_meta')

    # Initialize config dict with sample metadata and run params sections
    sample_info_dict = get_sample_info_dict(sample_name, sample_meta_str)
    run_params_dict = command_line_args_dict
    config = {
        'sample': sample_info_dict,
        'run': run_params_dict
    }

    # Add config variables from the default config file
    with open(default_config_file_path) as f:
        default_config_file_dict = pytoml.load(f)

    # 'run' and 'sample' section names are reserved for CLI args
    has_illegal_section = ('run' in default_config_file_dict
                           or 'sample' in default_config_file_dict)
    if has_illegal_section:
        raise KeyError("Config file may not contain reserved sections"
                       " 'run' and 'sample'")

    # Note that this update is not nested, because there can be no overlap
    # between the sections in the config file and the sections already present
    config.update(default_config_file_dict)

    # Add config variables from user config file if present
    # update_nested_dict will raise a KeyError if it encounters a key
    # not defined in the base dict
    if user_config_file_path:
        with open(user_config_file_path) as f:
            user_config_file_dict = pytoml.load(f)
        update_nested_dict(base_dict=config,
                           custom_dict=user_config_file_dict)

    # Copy output_dir from run parameters to the paths dict
    # Required for expansion of paths with expand_path
    output_dir = config['run']['output_dir']
    if not os.path.isabs(output_dir):
        raise ValueError("The output_dir path must be absolute.")
    config['paths']['output_dir'] = output_dir

    # Expand braced fields in paths
    path_names = list(config['paths'].keys())
    for curr_path_name in path_names:
        config['paths'][curr_path_name] = expand_path(config, curr_path_name)

    # Prepend output dir to paths if they are relative
    for curr_path_name in path_names:
        if not os.path.isabs(config['paths'][curr_path_name]):
            config['paths'][curr_path_name] = os.path.join(
                config['paths']['output_dir'], config['paths'][curr_path_name])

    return config


def expand_path(config, path_name):
    field_pattern = r'{(.*?)}'
    path_to_expand = config['paths'][path_name]
    field_names_unexpanded_path = re.findall(field_pattern,
                                             path_to_expand)
    for curr_field_name in field_names_unexpanded_path:
        if curr_field_name in config['sample']:
            curr_field_value = config['sample'][curr_field_name]
        elif curr_field_name in config['paths']:
            curr_field_value = config['paths'][curr_field_name]
            if re.search(field_pattern, curr_field_value):
                curr_field_value = expand_path(config, curr_field_name)
                config['paths'][curr_field_name] = curr_field_value
        else:
            raise ValueError(f"Can't expand path {path_to_expand}")
        path_to_expand = path_to_expand.replace(
            '{' + curr_field_name + '}',
            curr_field_value)

    return path_to_expand


def get_sample_info_dict(sample_name, sample_meta_str=None):
    sample_metadata = OrderedDict(name=sample_name)
    if sample_meta_str:
        add_infos = OrderedDict((x.split('=')[0], x.split('=')[1])
                                for x in sample_meta_str.split(','))
        sample_metadata.update(add_infos)
    return sample_metadata


def update_nested_dict(base_dict, custom_dict):
    for key, custom_value in custom_dict.items():

        try:
            default_value = base_dict[key]
        except KeyError:
            raise KeyError(f"Config file contains unknown element '{key}'")

        if type(custom_value) != type(default_value):
            raise TypeError(f"Can't overwrite type {type(default_value)}"
                             f" with type {type(custom_value)} for"
                             f"variable '{key}'")

        if isinstance(custom_value, dict):
            recursively_updated_dict = update_nested_dict(
                default_value, custom_value)
            base_dict[key] = recursively_updated_dict
        else:  # key points to 'scalar' value, not mapping
            base_dict[key] = custom_dict[key]

    return base_dict
