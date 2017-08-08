""" Demo workflow using mqc

How do I run this?
-----------------
snakemake \
--snakefile /home/kraemers/projects/mqc/mqc/tests/demo_snakefile.py \
-n
# --cluster "qsub -l walltime={params.walltime},mem={params.mem},nodes=1:ppn={params.cores} -N {params.name}" \
# --jobs 20 \
"""
# shell.executable("/bin/bash")
# shell.prefix("source ~/.bashrc; echo 'This job is running in bash'")

import os.path as op

sandbox_dir = "/icgc/dkfzlsdf/analysis/B080/kraemers/projects/mcall_qc/sandbox"
genome_dir = f"{sandbox_dir}/genomes/GRCm38mm10_PhiX_Lambda"
grcm38_fa =  f"{genome_dir}/GRCm38mm10_PhiX_Lambda.fa"
user_config_file = "/icgc/dkfzlsdf/analysis/B080/kraemers/projects/mcall_qc/sandbox/user_config.toml"
hsc_rep1_dir = f"{sandbox_dir}/hsc_rep1"
autosomes = [str(i) for i in range(1,20)]
other_chroms = ['X', 'Y', 'MT', 'phix', 'L']
bam_template = "/icgc/dkfzlsdf/analysis/hs_ontogeny/results/wgbs/results_per_pid/hsc_rep1/alignment/blood_hsc_rep1_merged.mdup.bam"

rule all:
    input:
        index_files = expand("{genome_dir}/GRCm38mm10_PhiX_Lambda_{motif_str}_{chrom}.bed.gz",
                             genome_dir=genome_dir,
                             # Indices for other contigs will also be created
                             # not referenced here, because unused
                             chrom=autosomes + other_chroms,
                             motif_str=['CG', 'CG-CHG-CHH']),

        mbias_stats = expand("{hsc_rep1_dir}/qc_stats/hsc_1_mbias-stats_{motif_str}.{ext}",
                             hsc_rep1_dir=hsc_rep1_dir,
                             motif_str=['CG', 'CG-CHG-CHH'], ext=['tsv', 'p']),

        cg_mcalls = expand("{hsc_rep1_dir}/mcalls/hsc_1_{motif_str}_{chrom}.bed.gz)",
                           hsc_rep1_dir=hsc_rep1_dir,
                           motif_str = 'CG',
                           chrom=autosomes),

        # all_motif_mcalls = expand("{hsc_rep1_dir}/mcalls/hsc_1_{motif_str}_{chrom}.bed.gz)",
        #                           hsc_rep1_dir=hsc_rep1_dir,
        #                           motif_str=['CG', 'CHG', 'CHH'],
        #                           chrom=autosomes),

rule make_all_motifs_index:
    input: grcm38_fa
    params:
        output_dir = genome_dir,
        walltime = '08:00:00',
        mem = '8g',
        cores = '12',
        name = 'make_all_motifs_index',
    output:
        index_files = expand("{genome_dir}/GRCm38mm10_PhiX_Lambda_CG-CHG-CHH_{chrom}.bed.gz",
                             genome_dir=genome_dir,
                             chrom=autosomes + other_chroms)
    benchmark: f"{genome_dir}/index-benchmark_all-motifs.txt"
    shell:
        """
        mqc make_index --genome_fasta {input} \
        --output_dir {params.output_dir} \
        --cores 12 \
        --cg --chg --chh
        """

rule make_cg_index:
    input: grcm38_fa
    output:
        index_files = expand("{genome_dir}/GRCm38mm10_PhiX_Lambda_CG_{chrom}.bed.gz",
                             genome_dir=genome_dir,
                             chrom=autosomes + other_chroms)
    params:
        output_dir = genome_dir,
        walltime = '08:00:00',
        mem = '8g',
        cores = '12',
        name = 'make_cg_motifs_index',
    benchmark: f"{genome_dir}/index-benchmark_cg.txt"
    shell:
        """
        mqc make_index --genome_fasta {input} \
        --output_dir {params.output_dir} \
        --cores 12 \
        --cg
        """

#TODO: log message for mbias stats generation
localrules: get_stats
rule get_stats:
    input:
        bam=bam_template,
        index_files = expand("{genome_dir}/GRCm38mm10_PhiX_Lambda_{{motif_str}}_{chrom}.bed.gz",
                             genome_dir=genome_dir,
                             chrom=autosomes),
    params:
        config_file = user_config_file,
        output_dir = hsc_rep1_dir,
        motif_csv = lambda wildcards: ','.join(wildcards.motif_str.split('-')),
        walltime = '48:00:00',
        mem = '12g',
        cores = '10',
        name = 'get_stats_{motif_str}',
    output:
        mbias_stats = expand("{hsc_rep1_dir}/qc_stats/hsc_1_mbias-stats_{{motif_str}}.{ext}",
                              hsc_rep1_dir = hsc_rep1_dir,
                              ext=['tsv', 'p']),
    shell:
        """
        mqc stats \
            --bam {input.bam} \
            --config_file {params.config_file} \
            --output_dir {params.output_dir} \
            --sample_name hsc_1 \
            --sample_meta population=HSC,rep=rep1 \
            --cores 10 \
            --motifs {params.motif_csv} \
            {input.index_files}
        """


common_meth_call_params = dict(
    config_file = user_config_file,
    output_dir = hsc_rep1_dir,
    mem = '12g',
    cores = '10',
)

mcall_command = (
    'mqc call'
    ' --bam {input.bam}'
    ' --config_file {params.config_file}'
    ' --output_dir {params.output_dir}'
    ' --sample_name hsc_1'
    ' --sample_meta population=HSC,rep=rep1'
    ' --cores 10'
    ' {input.index_files}'
)

rule call_meth_for_CG:
    input:
        bam=bam_template,
        #TODO: can expand deal with scalar strings?
        index_files = expand("{genome_dir}/GRCm38mm10_PhiX_Lambda_{motif_str}_{chrom}.bed.gz",
                             genome_dir=genome_dir,
                             motif_str='CG',
                             chrom=autosomes),
    params:
        **common_meth_call_params,
        walltime = '60:00:00',
        name = 'mcall_CG',
    output:
        expand("{hsc_rep1_dir}/mcalls/hsc_1_{motif_str}_{chrom}.bed.gz)",
               hsc_rep1_dir=hsc_rep1_dir,
               motif_str = 'CG',
               chrom=autosomes),
    shell: mcall_command


## Calling rules for all motifs and CG are ambiguous,
## can't both be active at the same time
## could probably be handled with if blocks, but
## good enough for now
# rule call_meth_for_all_motifs:
#     input:
#         bam=bam_template,
#         index_files = expand("{genome_dir}/GRCm38mm10_PhiX_Lambda_{motif_str}_{chrom}.bed.gz",
#                              genome_dir=genome_dir,
#                              motif_str='CG-CHG-CHH',
#                              chrom=autosomes),
#     params:
#         **common_meth_call_params,
#         walltime = '60:00:00',
#         name = 'mcall_CG-CHG-CHH',
#     output:
#         expand("{hsc_rep1_dir}/mcalls/hsc_1_{motif_str}_{chrom}.bed.gz)",
#                 hsc_rep1_dir=hsc_rep1_dir,
#                 motif_str=['CG', 'CHG', 'CHH'],
#                 chrom=autosomes),
#     shell: mcall_command

## Obsolete, because we are now using the reference genome in ngs share
# which is also used by the pipeline
# localrules: download_fa
# rule download_ucsc_mm10:
#     # TODO: the glob sorts chr10 before chr1, is that a problem?
#     # only worked on tbi-worker, not in interactive session - why?
#     # No mm10.fa.gz with all chromosomes in one file available
#     # for human, just download the entire genome fasta at once:
#     # http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
#     output:
#         f"{genome_dir}/mm10.fa.gz"
#     params:
#         genome_dir=genome_dir
#     shell:
#         """
#         cd {params.genome_dir}
#         rsync -avzP rsync://hgdownload.cse.ucsc.edu/goldenPath/mm10/chromosomes/ .
#         chrom_files=(*.fa.gz)
#         cat ${{chrom_files[@]}} > {output}
#         rm ${{chrom_files[@]}}
#         echo 'Done: Genoma fasta in {output}'
#         """
