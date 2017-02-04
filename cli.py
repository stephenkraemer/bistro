import mqc
import click


@click.group()
def cli():
    pass

@cli.command()
@click.help_option()
@click.option('--index', required=True,
              type=click.Path(exists=True, readable=True))
@click.option('--bam', required=True,
              type=click.Path(exists=True, readable=True))
@click.option('--output_dir', required=True,
              type=click.Path(writable=True))
@click.option('--config_file', required=True,
              type=click.Path(exists=True, readable=True))
@click.option('--sample_name', required=True)
@click.option('--sample_meta')
@click.pass_context
def qc_run(ctx, bam, index, config_file, output_dir, sample_name, sample_meta):

    config = mqc.config.get_config_dict(args_dict=ctx.params,
                                        config_file_path=config_file)

    mqc.qc_run.qc_run(bam_path=bam,
                      index_file_path=index,
                      config=config,
                      meth_metrics_dir_abs=output_dir,
                      sample_name=sample_name)


if __name__ == '__main__':
    cli()


# TODO-learn: --test, default=True seems to have created string option when used with --test False?
# @click.option('--config', required=True, is_eager=True, callback=set_config_file)
#
#
# def set_config_file(ctx, config_option_name, config_option_value):
#     with open(config_option_value) as f:
#         config_dict = pytoml.load(f)
#     print(config_dict)
#     ctx.default_map = config_dict['command_line_options']
