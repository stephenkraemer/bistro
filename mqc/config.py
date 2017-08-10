"""Config file and command line parsing"""

import copy
import os
import pytoml
import re
from collections import OrderedDict


def assemble_config_vars(command_line_args_dict: dict,
                         default_config_file_path: str,
                         command: str,
                         user_config_file_path: str = '', ):
    """Assemble configuration variables in one dict and expand paths

    Parameters
    ----------
    command_line_args_dict:
        Command line arguments, provided as dict by command line parser
    default_config_file_path: str
        Required, typically every mqc based tool will have its own default
        config file
    user_config_file_path : str
        Optional path to user config file

    Returns
    -------
    dict:
        dict integrating all configuration variables


    **Implementation notes**

    *Sample metadata and run parameter specification*

    This function expects two CLI parameters: sample_name and output_dir.
    The path to output_dir must be absolute. Optionally, a string specifying
    further sample metadata can be given as CL arg, it will be parsed
    assuming this format: key=value[, key=value...]. All sample metadata are
    made available in the section config[ 'sample'], e.g.

    | config['sample']['name'] = 'sample_name'
    | config[ 'sample']['metadata_key'] = 'metadata_value'

    All command line args except the sample_name and sample_metadata are
    collected in config['run']

    *Path expansion*

    All file paths must be defined in the section 'paths'.

    All relative paths are assumed to be relative to output dir, and expanded
    to absolute paths accordingly.

    Fields ({}) in path config variables may point to either other path
    variables or sample metadata (e.g. {name}). Fields are expanded
    recursively. Currently, circular references are not catched,
    so be careful. Note that it is not possible to protect fields using
    double braces, as one would do in python format strings.
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
    config['paths'][command]['output_dir'] = output_dir

    # Expand braced fields in paths
    path_names = list(config['paths'][command].keys())
    for curr_path_name in path_names:
        config['paths'][command][curr_path_name] = _expand_path(config, curr_path_name, command)

    # Prepend output dir to paths if they are relative
    for curr_path_name in path_names:
        if not os.path.isabs(config['paths'][command][curr_path_name]):
            config['paths'][command][curr_path_name] = os.path.join(
                config['paths'][command]['output_dir'], config['paths'][command][curr_path_name])

    return config


def _expand_path(config, path_name, command):
    """Expand single path

    Path expansion is explained in the docstring for assemble_config_vars
    """
    field_pattern = r'{(.*?)}'
    path_to_expand = config['paths'][command][path_name]
    field_names_unexpanded_path = re.findall(field_pattern,
                                             path_to_expand)
    for curr_field_name in field_names_unexpanded_path:
        if curr_field_name in config['sample']:
            curr_field_value = config['sample'][curr_field_name]
        elif curr_field_name in config['run']:
            curr_field_value = config['run'][curr_field_name]
        elif curr_field_name in config['paths'][command]:
            curr_field_value = config['paths'][command][curr_field_name]
            if re.search(field_pattern, curr_field_value):
                curr_field_value = _expand_path(config, curr_field_name, command)
                config['paths'][command][curr_field_name] = curr_field_value
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
    """Update nested dict with rather strict rules

    Raises
    ------
    TypeError
        if type of the update value does not match the original value
    KeyError
        if custom_dict contains a key not present in base_dict (also
        within nested dicts)
    """
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
