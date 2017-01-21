import click
import pytoml

def get_config_dict(config_dict_path):
    with open(config_dict_path) as f:
        config_dict = pytoml.load(f)
    return config_dict
#
#
# def set_config_file(ctx, config_option_name, config_option_value):
#     with open(config_option_value) as f:
#         config_dict = pytoml.load(f)
#     print(config_dict)
#     ctx.default_map = config_dict['command_line_options']


@click.group()
def cli():
    pass


@cli.command()
@click.help_option()
@click.option('--index')
@click.option('--bam', required=True)
@click.option('--output_dir', required=True)
@click.option('--config_file', required=True)
@click.option('--sample_name')
# @click.option('--config', required=True, is_eager=True, callback=set_config_file)
def qc_run(bam, index, config_file, output_dir, sample_name):
    from mqc.qc_run import qc_run
    config = get_config_dict(config_file)
    qc_run(bam_path=bam,
           index_file_path=index,
           config=config,
           output_dir_abspath=output_dir,
           sample_name=sample_name)

if __name__ == '__main__':
    cli()
