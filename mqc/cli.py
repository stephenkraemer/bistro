from collections import OrderedDict

from os.path import isabs

from mqc.index import parallel_index_generation
import click


@click.group()
def mqc():
    pass

@mqc.command()
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
def qc_run(ctx, bam, index, config_file, output_dir, sample_name):

    config = mqc.config.get_config_dict(args_dict=ctx.params,
                                        config_file_path=config_file)

    mqc.qc_run.qc_run(bam_path=bam,
                      index_file_path=index,
                      config=config,
                      meth_metrics_dir_abs=output_dir,
                      sample_name=sample_name)

@mqc.command()
@click.help_option()
@click.option('--fasta_path_template', required=True,
              help="Template path for fasta files, indicating the"
                   " chromosome field, e.g. 'mm10_chr{chr}.fa.gz'")
@click.option('--output_path_template', required=True,
              help="Template path for index files to be generated, "
                   "with field for chromosome indication, e.g."
                   "'mm10_chr{chr}_CG-CHH-CHG.bed.gz."
                   " Suffix must be specified in the template,"
                   " note that file will be gzipped.")
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
def make_index(fasta_path_template, output_path_template, cores,
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

    parallel_index_generation(fasta_path_template=fasta_path_template,
                              output_path_template=output_path_template,
                              motifs=motifs,
                              annotations=annotations,
                              cores=cores)

