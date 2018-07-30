"""Command line interface for tools

Registered as console script entry point (mqc)

Based on click package
"""

# pylint: disable=unused-argument
import argparse
import copy
import json
import sys
from collections import OrderedDict
from textwrap import dedent
from typing import Dict

import click

from mqc.config import assemble_config_vars
from mqc.index import start_parallel_index_generation
from mqc.mcall_run import collect_mbias_stats, run_mcalling
from mqc.merge_strands import BedMethCallsMerger
from mqc.utils import get_resource_abspath
from mqc.mbias import compute_mbias_stats, mbias_stat_plots
from mqc.coverage import analyze_coverage


@click.group()
def mqc():
    pass


input_click_path = click.Path(exists=True, readable=True,
                              dir_okay=False, resolve_path=True)


# mqc stats
# =========
# noinspection PyUnusedLocal
@mqc.command()
@click.option('--bam', required=True,
              type=input_click_path)
@click.argument('index_files', nargs=-1,
                required=True,
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
@click.option('--read_length', type=int, required=True)
@click.pass_context
def stats(ctx, bam, index_files,
          output_dir, config_file,
          sample_name, sample_meta, cores, motifs, read_length) -> None:
    """Gather Mbias stats"""

    default_config_file = get_resource_abspath('config.default.toml')

    user_config_file = config_file if config_file else ''

    cli_params = copy.deepcopy(ctx.params)
    cli_params['motifs'] = cli_params['motifs'].upper().split(',')
    cli_params['motifs_str'] = '-'.join(cli_params['motifs'])

    config = assemble_config_vars(cli_params,
                                  default_config_file_path=default_config_file,
                                  user_config_file_path=user_config_file)

    collect_mbias_stats(config=config)

    print(f"Stats collection for {sample_name} successful."
          f"See {output_dir} for results.")


# ==============================================================================
#                             mqc evaluate_mbias
# ==============================================================================


# noinspection PyUnusedLocal
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
@click.option('--use_cached_mbias_stats', is_flag=True)
@click.option('--plateau_detection', required=True,
              help=('Plateau detection parameters. Expects JSON dict containing the name'
                    'of the algorithm as well as parameters for the algorithm. See documentation'
                    'for details'))
@click.pass_context
def evaluate_mbias(ctx, config_file, motifs, output_dir,
                   sample_name, sample_meta, use_cached_mbias_stats, plateau_detection) -> None:
    """Process M-bias stats"""
    default_config_file = get_resource_abspath('config.default.toml')
    user_config_file = config_file if config_file else ''

    cli_params = copy.deepcopy(ctx.params)

    cli_params['motifs'] = motifs.split(',')
    cli_params['motifs_str'] = '-'.join(cli_params['motifs'])

    cli_params['plateau_detection'] = json.loads(cli_params['plateau_detection'])
    try:
        assert isinstance(cli_params['plateau_detection']['algorithm'], str)
    except (TypeError, KeyError, AssertionError) as e:
        raise ValueError('Can not parse the plateau_detection JSON dict correctly') from e

    config = assemble_config_vars(cli_params,
                                  default_config_file_path=default_config_file,
                                  user_config_file_path=user_config_file)

    compute_mbias_stats(config)


# ==============================================================================
#                             mqc mbias_plots
# ==============================================================================


@mqc.command()
@click.option('--output_dir', required=True,
              type=click.Path(exists=False, file_okay=False,
                              writable=True, resolve_path=True))
@click.option('--mbias_plot_config')
@click.option('--datasets', required=True)
def mbias_plots(output_dir, mbias_plot_config, datasets) -> None:
    """Plot processed M-bias stats"""
    datasets_dict = dict_from_kwarg_cli_option(datasets)
    mbias_stat_plots(output_dir=output_dir,
                     compact_mbias_plot_config_dict_fp=mbias_plot_config,
                     dataset_name_to_fp=datasets_dict)


# ==============================================================================
#                             mqc call
# ==============================================================================
# noinspection PyUnusedLocal
@mqc.command()
@click.option('--bam', required=True, type=input_click_path)
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
@click.option('--max_read_length', required=True, type=int)
@click.option('--trimming', help=(
        'Desired trimming in the form command::param. Possible commands:'
        'frag_end, read_end, cutting_sites'
        'frag_end and read_end are parametrized through a json object'
        'cutting_sites expects the path to a dataframe, either created '
        'with evaluate_stats or user-specified. See manual for details.')
              )
@click.option('--output_formats', default='bed', help=(
        'any csv combination of [bed bismark stratified_bed].'
        ' The options bed and stratified_bed are mutually exclusive.'
        ' Further output formats will soon be added, in particular VCF'
        ' for methylation calling with SNP calls')
              )
@click.pass_context
def call(ctx, bam, index_files,
         output_dir, config_file,
         sample_name: str, sample_meta: str, cores: str, trimming: str,
         max_read_length: int, output_formats: str) -> None:
    """Methylation calling tool

    Command line tool for generating QC-filtered methylation calls
    in various output formats, using common filtering strategies.

    This tool wraps the base QcAndMethCallingRun.
    """

    # To add new output formats
    # 1. Add to help message of --output_formats option
    # 2. Add to allowed_output_formats_set
    # 3. Modify QcAndMethCallingRun._get_visitors to include appropriate
    #    Visitors when required

    default_config_file = get_resource_abspath('config.default.toml')

    user_config_file = config_file if config_file else ''

    try:
        trimm_command, trimm_param = trimming.split('::')  # pylint: disable=unused-variable
    except ValueError as e:
        print('Could not find parameters for trimming command. '
              'Did you separate with a "::"?', file=sys.stderr)
        raise e
    if trimm_command not in 'frag_end read_end cutting_sites'.split():
        raise ValueError('The trimming command is not known.')

    cli_params = copy.deepcopy(ctx.params)
    # TODO: avoid hard coding!
    # output_file_template = f"{index_output_dir}/{genome_name}_{motifs_str}_{{chrom}}.bed.gz"
    motifs_str = index_files[0].split('_')[-2]
    cli_params['motifs_str'] = motifs_str
    cli_params['motifs'] = motifs_str.split('-')

    allowed_output_formats_set = set('bed bismark stratified_bed'.split())
    output_formats_list = output_formats.split(',')
    assert set(output_formats_list) <= allowed_output_formats_set, \
        f'Aborting: unknown output formats in {output_formats_list}'
    cli_params['output_formats'] = output_formats_list

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
@click.option('--cores', default=1)
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
               cg, chg, chh, cwg, triplet_seq, seq_context) -> None:
    """Create augmented index files

    Can add triplet seq and more general sequence context information.
    """
    if chg and cwg:
        raise ValueError("--chg and --cwg are mutually exclusive!")

    motifs = []
    if cg:
        motifs.append('CG')
    if chg:
        motifs.append('CHG')
    if chh:
        motifs.append('CHH')
    if cwg:
        motifs.append('CWG')

    if not motifs:
        raise ValueError('You have to select at least one motif.')

    annotations = OrderedDict((('triplet_seq', triplet_seq),
                               ('seq_context', seq_context)))

    start_parallel_index_generation(genome_fasta=genome_fasta,
                                    index_output_dir=output_dir,
                                    motifs=motifs,
                                    annotations=annotations,
                                    cores=cores)


# ==============================================================================
#                             mqc evaluate_calls
# ==============================================================================


# noinspection PyUnusedLocal
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
def evaluate_calls(ctx, config_file, motifs, output_dir, sample_name, sample_meta) -> None:
    """Process statistics gathered during methylation calling"""
    default_config_file = get_resource_abspath('config.default.toml')
    user_config_file = config_file if config_file else ''

    cli_params = copy.deepcopy(ctx.params)
    cli_params['motifs'] = motifs.split(',')
    cli_params['motifs_str'] = '-'.join(cli_params['motifs'])
    config = assemble_config_vars(cli_params,
                                  default_config_file_path=default_config_file,
                                  user_config_file_path=user_config_file)
    analyze_coverage(config)


# ==============================================================================
# Util functions
# ==============================================================================
def dict_from_kwarg_cli_option(value: str) -> Dict[str, str]:
    return {s.split('=')[0]: s.split('=')[1] for s in value.split(',')}


# new argparse-based interface
# ============================

tools_short_help = {
    'meth_calls create': 'Methylation calling workhorse, various output formats and processing options',
    'meth_calls merge': 'Merge strand-resolved calls',
    'meth_calls qc': 'Create QC plots from (optionally stratified) methylation calls for one or more samples',
    'index create': 'Create index, optionally including sequence context and other annotations',
    'index merge': 'Merge strand-resolved index positions'
}

def bistro_parser_help_message(unused_args):
    print(dedent(f'''\
    
    Welcome to bistro - methylation calls a la carte
    
    usage: bistro <subcommand> <args>
    
    [ Methylation calling ]
        meth_calls create        {tools_short_help["meth_calls create"]}
        meth_calls qc            {tools_short_help["meth_calls qc"]}
    
    [ Index ]
        index create             {tools_short_help["index create"]}
        index merge              {tools_short_help["index merge"]}
    
    [ General help ]
        --help, -h               Print this help menu
        --version
    '''))
    sys.exit(1)

def meth_calls_parser_help_message(unused_args):
    print(dedent(f'''\
    usage: bistro meth_calls <command> <args>
    
    Creation and processing of methylation calls.
    
    [ commands ]
        meth_calls create        {tools_short_help["meth_calls create"]}
        meth_calls merge         {tools_short_help["meth_calls merge"]}
        meth_calls qc            {tools_short_help["meth_calls qc"]}
    
    '''))
    sys.exit(1)

def get_parser() -> argparse.ArgumentParser:
    """Create ArgumentParser for bistro command line tool

    This function can be used to create automatic documentation
    with the sphinx argparse extension.
    """

    bistro = argparse.ArgumentParser('bistro', add_help=False)
    bistro.add_argument('-h', '--help', action='store_true')
    bistro.set_defaults(func=bistro_parser_help_message)
    bistro_subparsers = bistro.add_subparsers(title='commands')

    add_meth_calls_tool(bistro_subparsers)

    return bistro

def add_meth_calls_tool(bistro_subparsers: argparse._SubParsersAction) -> None:
    meth_calls = bistro_subparsers.add_parser('meth_calls', add_help=False)
    meth_calls.add_argument('-h', '--help', action='store_true')
    meth_calls.set_defaults(func=meth_calls_parser_help_message)
    meth_calls_subparser = meth_calls.add_subparsers(title='commands')

    meth_calls_merge = meth_calls_subparser.add_parser(
            'merge',
            usage='bistro meth_calls merge [<options>] <strand_resolved> <merged>',
            description=dedent(f'''\
            {tools_short_help["meth_calls merge"]}
            
            Merge strand-resolved calls for symmetric motifs. Currently, only merging
            of CG motifs is implemented. Merging of CHG motifs will be added in the future.
            If you need CHG motif merging, please get in contact via github.
            
            Filepaths must be given as absolute paths.
            '''),
            help=tools_short_help['meth_calls merge'],
            formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    meth_calls_merge.set_defaults(func=meth_calls_merge_runner)
    meth_calls_merge_required = meth_calls_merge.add_argument_group('Required')
    meth_calls_merge_required.add_argument(
            '--strand_resolved',
            help='Path to strand-resolved BED or stratified-BED file',
            required=True)
    meth_calls_merge_required.add_argument(
            '--merged',
            help='Path for the created merged methylation calls',
            required=True)
    meth_calls_merge_optional = meth_calls_merge.add_argument_group('Optional')
    meth_calls_merge_optional.add_argument('--verbose', '-v',
                                           action='store_true', default=False)


def meth_calls_merge_runner(args):
    merger = BedMethCallsMerger(strand_resolved_meth_calls=args.strand_resolved,
                                merged_meth_calls=args.merged, verbose=args.verbose)
    merger.run()

def bistro():
    args = get_parser().parse_args()
    args.func(args)
    print('Done - thanks for using bistro.')
