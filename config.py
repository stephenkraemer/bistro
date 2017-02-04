from collections import OrderedDict
import pytoml
import os


def get_config_dict(args_dict, config_file_path=None):

    sample_name = args_dict.pop('sample_name')
    sample_meta_str = args_dict.pop('sample_meta')
    sample_info_dict = get_sample_info_dict(sample_name, sample_meta_str)

    run_params_dict = args_dict

    config = {
    'sample': sample_info_dict,
    'run': run_params_dict
    }

    package_dir = os.path.dirname(__file__)
    default_config_file_path = os.path.join(package_dir, 'config.default.toml')
    with open(default_config_file_path) as f:
        default_config_file_dict = pytoml.load(f)
    config.update(default_config_file_dict)

    if config_file_path:
        with open(config_file_path) as f:
            config_file_dict = pytoml.load(f)

        for k in config.keys():
            config[k].update(config_file_dict.get(k, {}))

    paths_dict = config.pop('paths')
    paths_dict['output_dir'] = config['run']['output_dir']

    def get_abs_path_obj(rel_path_template_str):
        rel_path_str = rel_path_template_str.format(**paths_dict, **config['sample'])
        return os.path.join(paths_dict['output_dir'], rel_path_str)

    """
    Note that this creates a new path dict, where all paths are absolute,
    inclusively the coverage_dir, mbias_dir etc. paths for the subdirectories.
    However, because we iterate over the elements in the old paths_dict and then
    generate a new expanded_paths_dict, this expansion of the subdirectory paths
    is unproblematic: during expansion, they are always retrieved from the old,
    unmodified paths_dict
    """
    expanded_paths_dict = { name: get_abs_path_obj(rel_path_template_str)
        for name, rel_path_template_str in paths_dict.items() }

    config['paths'] = expanded_paths_dict

    return config


def get_sample_info_dict(sample_name, sample_meta_str=None):
    sample_metadata = OrderedDict(name= sample_name)
    if sample_meta_str:
        add_infos = OrderedDict((x.split('=')[0], x.split('=')[1])
                                           for x in sample_meta_str.split(','))
        sample_metadata.update(add_infos)
    return sample_metadata
