import click
import multiprocessing as mp
import pytoml


# def get_config_dict(config_dict_path):
#     with open(config_dict_path) as f:
#         config_dict = pytoml.load(f)
#     return config_dict
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
@click.option('--output_file', required=True)
@click.option('--config')
@click.option('--pos_per_process', default=1000)
# @click.option('--config', required=True, is_eager=True, callback=set_config_file)
@click.option('--n_cores', type=int, default=1)
def call(bam, index, config, n_cores, output_file, pos_per_process):
    lock = mp.Lock()
    # config = get_config_dict(config)
    from mqc import call_methylation_parallelized
    call_methylation_parallelized(
        bam_file_path=bam,
        index_file_path=index,
        n_cores=n_cores,
        output_file=output_file,
        pos_per_process=pos_per_process,
        lock=lock, )


if __name__ == '__main__':
    cli()
