import copy
import os.path as op
from collections import OrderedDict

import click
from mqc.config import assemble_config_vars
from mqc.index import start_parallel_index_generation
from mqc.mcall_run import collect_stats, run_mcalling
from mqc.utils import get_resource_abspath


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

    default_config_file = get_resource_abspath('config.default.toml')

    user_config_file = config_file if config_file else ''

    cli_params = copy.deepcopy(ctx.params)
    cli_params['motifs'] = cli_params['motifs'].upper().split(',')
    cli_params['motifs_str'] = '-'.join(cli_params['motifs'])


    config = assemble_config_vars(cli_params,
                                  default_config_file_path=default_config_file,
                                  user_config_file_path=user_config_file)

    collect_stats(config=config)

    print(f"Stats collection for {sample_name} successful."
          f"See {output_dir} for results.")

#==============================================================================
#                             mqc evaluate_mbias
#==============================================================================
from mqc.mbias import analyze_mbias_counts
@mqc.command()
@click.option('--config_file', type=input_click_path, help='[optional]')
@click.option('--motifs', help='e.g. CG or CG,CHG,CHH')
@click.option('--output_dir', required=True,
              type=click.Path(exists=False, file_okay=False,
                              writable=True, resolve_path=True))
@click.option('--sample_name', required=True)
@click.option('--sample_meta',
              help="Pass additional metadata as"
                   " 'key=value,key2=value2' [optional]")
@click.pass_context
def evaluate_mbias(ctx, config_file, motifs, output_dir, sample_name, sample_meta):

    default_config_file = get_resource_abspath('config.default.toml')
    user_config_file = config_file if config_file else ''
    cli_params = copy.deepcopy(ctx.params)
    cli_params['motifs'] =  motifs.split(',')
    cli_params['motifs_str'] = '-'.join(cli_params['motifs'])
    config = assemble_config_vars(cli_params,
                                  default_config_file_path=default_config_file,
                                  user_config_file_path=user_config_file)
    analyze_mbias_counts(config)


#==============================================================================
#                             mqc call
#==============================================================================
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
@click.option('--use_mbias_fit', is_flag=True)
@click.option('--strat_beta_dist', is_flag=True)

@click.pass_context
def call(ctx, bam, index_files,
          output_dir, config_file,
          sample_name, sample_meta, cores, use_mbias_fit, strat_beta_dist):
    """Methylation calling"""

    default_config_file = get_resource_abspath('config.default.toml')

    user_config_file = config_file if config_file else ''

    cli_params = copy.deepcopy(ctx.params)
    #TODO: avoid hard coding!
    # output_file_template = f"{index_output_dir}/{genome_name}_{motifs_str}_{{chrom}}.bed.gz"
    motifs_str = index_files[0].split('_')[-2]
    cli_params['motifs_str'] = motifs_str
    cli_params['motifs'] = motifs_str.split('-')

    config = assemble_config_vars(cli_params,
                                  default_config_file_path=default_config_file,
                                  user_config_file_path=user_config_file)

    run_mcalling(config=config)


# mqc make_index
# =============================================================================
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

#==============================================================================
#                             mqc evaluate_calls
#==============================================================================
from mqc.coverage import analyze_coverage
from mqc.beta_value import analyze_stratified_beta_values
@mqc.command()
@click.option('--config_file', type=input_click_path, help='[optional]')
@click.option('--motifs', help='e.g. CG or CG,CHG,CHH')
@click.option('--output_dir', required=True,
              type=click.Path(exists=False, file_okay=False,
                              writable=True, resolve_path=True))
@click.option('--sample_name', required=True)
@click.option('--sample_meta',
              help="Pass additional metadata as"
                   " 'key=value,key2=value2' [optional]")
@click.option('--strat_beta_dist', is_flag=True)
@click.pass_context
def evaluate_calls(ctx, config_file, motifs, output_dir, sample_name, sample_meta, strat_beta_dist):

    default_config_file = get_resource_abspath('config.default.toml')
    user_config_file = config_file if config_file else ''

    cli_params = copy.deepcopy(ctx.params)
    cli_params['motifs'] = motifs.split(',')
    cli_params['motifs_str'] = '-'.join(cli_params['motifs'])
    config = assemble_config_vars(cli_params,
                                  default_config_file_path=default_config_file,
                                  user_config_file_path=user_config_file)
    analyze_coverage(config)
    if(cli_params['strat_beta_dist']): analyze_stratified_beta_values(config)
