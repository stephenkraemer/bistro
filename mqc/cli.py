import click
import os.path as op

from collections import OrderedDict

import copy


from mqc.index import start_parallel_index_generation
from mqc.config import assemble_config_vars
from mqc.mcall_run import collect_stats


@click.group()
def mqc():
    pass

input_click_path = click.Path(exists=True, readable=True,
                              dir_okay=False, resolve_path=True)

# mqc stats
# =========
@mqc.command()
@click.option('--bam', required=True,
              type=input_click_path)
@click.argument('index_files', nargs=-1,
                type=click.Path(resolve_path=True, exists=True, dir_okay=False))
# TODO: which checks are necessary
@click.option('--output_dir', required=True,
              type=click.Path(exists=False, file_okay=False,
                              writable=True, resolve_path=True))
@click.option('--config_file', type=input_click_path, help='[optional]')
@click.option('--sample_name', required=True)
@click.option('--sample_meta',
              help="Pass additional metadata as"
                   " 'key=value,key2=value2' [optional]")
@click.option('--cores', default=1)
@click.option('--motifs', required=True,
              help=("Comma separated list of the motifs in the index files."
                    "For example: CG,CHG"))
@click.pass_context
def stats(ctx, bam, index_files,
          output_dir, config_file,
          sample_name, sample_meta, cores, motifs):
    """Gather Mbias stats"""

    package_top_level_dir = op.abspath(op.dirname(__file__))
    default_config_file = op.join(package_top_level_dir, 'config.default.toml')

    user_config_file = config_file if config_file else ''

    cli_params = copy.deepcopy(ctx.params)
    cli_params['motifs'] = cli_params['motifs'].upper().split(',')


    config = assemble_config_vars(cli_params,
                                  default_config_file_path=default_config_file,
                                  user_config_file_path=user_config_file)

    collect_stats(config=config)

    print(f"Stats collection for {sample_name} successful."
          f"See {output_dir} for results.")


# mqc make_index
# ==============
@mqc.command()
@click.help_option()
@click.option('--genome_fasta', required=True)
@click.option('--output_dir', required=True)
@click.option('--cores', default = 1)
@click.option('--cg', is_flag=True,
              help='Include CG positions in index file')
@click.option('--chg', is_flag=True)
@click.option('--chh', is_flag=True)
@click.option('--cwg', is_flag=True,
              help="Include CHG positions, but distinguish between CCG"
                   " and CWG. Cannot be used together with --chg")

@click.option('--triplet_seq', is_flag=True, help='Add info about base tripletts')
@click.option('--seq_context', default=0, type=int,
              help="Add info about sequence context, with INTEGER bases before and "
                   "after the cytosine. Set INTEGER to 0 to "
                   "disable sequence context [default].")
def make_index(genome_fasta, output_dir, cores,
               cg, chg, chh, cwg, triplet_seq, seq_context):

    if chg and cwg:
        raise ValueError("--chg and --cwg are mutually exclusive!")

    motifs = []
    if cg: motifs.append('CG')
    if chg: motifs.append('CHG')
    if chh: motifs.append('CHH')
    if cwg: motifs.append('CWG')

    if not motifs:
        raise ValueError('You have to select at least one motif.')

    annotations = OrderedDict((('triplet_seq', triplet_seq),
                              ('seq_context', seq_context)))

    start_parallel_index_generation(genome_fasta=genome_fasta,
                                    index_output_dir=output_dir,
                                    motifs=motifs,
                                    annotations=annotations,
                                    cores=cores)

