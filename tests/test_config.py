import os
import mqc
import pytest


def test_get_config(test_output_dir):

    test_config_file = os.path.join(os.path.dirname(__file__),
                                    'custom.config.toml')

    config = mqc.config.get_config_dict(
            args_dict={'sample_name': 'hsc_rep1',
                       'sample_meta': 'population=hsc,replicate=1',
                       'output_dir': test_output_dir},
            config_file_path=test_config_file)

    # Test that sample metadata are saved
    assert config['sample']['replicate'] == '1'

    # Test that defaults are used
    assert config['trimming']['max_flen_considered_for_trimming'] == 500

    # Test that defaults are overwritten by the custom config file
    assert config['conversion_error_calling']['max_number_of_unconv_cytosines'] == 10

    # Test that paths are constructed correctly

    ## Default path
    assert (config['paths']['coverage_counts_tsv']
            == '{output_dir}/coverage_dir/hsc_rep1.cov.counts.tsv'.format(
        output_dir=config['paths']['output_dir']
    ))

    # Custom path
    assert (config['paths']['coverage_counts_p'] ==
            '{coverage_dir}/{population}.cov.counts.p'.format(
                output_dir=config['paths']['output_dir'],
                coverage_dir=config['paths']['coverage_dir'],
                population=config['sample']['population']
            ))

def test_get_sample_info_dict():
    d1 = mqc.config.get_sample_info_dict('hsc_rep1', 'pop=a,rep=1')
    assert d1['pop'] == 'a'
    assert d1['name'] == 'hsc_rep1'

    d2 = mqc.config.get_sample_info_dict('b_cells_rep1')
    assert  d2['name'] == 'b_cells_rep1'
    with pytest.raises(KeyError):
        d2['pop']
