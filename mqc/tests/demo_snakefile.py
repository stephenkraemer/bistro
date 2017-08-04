"""
Install on cluster
Download chromosome fastas (or later: the whole genome)
Create index files
Run mqc stats on one chromosome
Create M-bias plots
Inspect M-bias plots
"""
import os.path as op

sandbox_dir = "/icgc/dkfzlsdf/analysis/B080/kraemers/projects/mcall_qc/sandbox"
mm10_genome_dir = f"{sandbox_dir}/genomes/mm10"
user_config_file = "/icgc/dkfzlsdf/analysis/B080/kraemers/projects/mcall_qc/sandbox/user_config.toml"
hsc_rep1_dir = f"{sandbox_dir}/hsc_rep1"


chroms = [19]
rule all:
    input:
        mm10_genome = f"{mm10_genome_dir}/mm10.fa.gz",
        # all_motifs_index_files_done_flag = f"{mm10_genome_dir}/.all_motifs.index_done.mm10.fa.gz",
        cg_index_files_done_flag = f"{mm10_genome_dir}/.cg.index_done.mm10.fa.gz",
        mbias_stats_tsv = f"{hsc_rep1_dir}/qc_stats/mbias_stats.tsv",
        mbias_stats_p = f"{hsc_rep1_dir}/qc_stats/mbias_stats.p",

rule download_fa:
    # TODO: the glob sorts chr10 before chr1, is that a problem?
    # only worked on tbi-worker, not in interactive session - why?
    # No mm10.fa.gz with all chromosomes in one file available
    # for human, just download the entire genome fasta at once:
    # http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
    output:
        f"{mm10_genome_dir}/mm10.fa.gz"
    params:
        mm10_genome_dir=mm10_genome_dir
    shell:
        """
        cd {params.mm10_genome_dir}
        rsync -avzP rsync://hgdownload.cse.ucsc.edu/goldenPath/mm10/chromosomes/ .
        chrom_files=(*.fa.gz)
        cat ${{chrom_files[@]}} > {output}
        rm ${{chrom_files[@]}}
        echo 'Done: Genoma fasta in {output}'
        """

rule make_all_motifs_index:
    input: f"{mm10_genome_dir}/mm10.fa.gz"
    params:
        output_dir = mm10_genome_dir
    output: touch(f"{mm10_genome_dir}/.all_motifs.index_done.mm10.fa.gz")
    benchmark: f"{mm10_genome_dir}/all_motifs_index_generation_benchmark.txt"
    shell:
        """
        mqc make_index --genome_fasta {input} \
        --output_dir {params.output_dir} \
        --cores 8 \
        --cg --chg --chh
        """

rule make_cg_index:
    input: f"{mm10_genome_dir}/mm10.fa.gz"
    params:
        output_dir = mm10_genome_dir
    output: touch(f"{mm10_genome_dir}/.cg.index_done.mm10.fa.gz")
    benchmark: f"{mm10_genome_dir}/cg_index_generation_benchmark.txt"
    shell:
        """
        mqc make_index --genome_fasta {input} \
        --output_dir {params.output_dir} \
        --cores 8 \
        --cg
        """

# TODO: make sure that output dir is created, also if not run with snakemake (which will create it automatically)
rule get_stats:
    input:
        bam="/icgc/dkfzlsdf/analysis/hs_ontogeny/results/wgbs/results_per_pid/hsc_rep1/alignment/blood_hsc_rep1_merged.mdup.bam",
        index_files = ["/icgc/dkfzlsdf/analysis/B080/kraemers/projects/mcall_qc/sandbox/genomes/mm10/mm10_CG_chr18.bed.gz",
                       "/icgc/dkfzlsdf/analysis/B080/kraemers/projects/mcall_qc/sandbox/genomes/mm10/mm10_CG_chr19.bed.gz"],
    params:
        config_file = user_config_file,
        output_dir = hsc_rep1_dir,
    output:
        mbias_stats_tsv = f"{hsc_rep1_dir}/qc_stats/mbias_stats.tsv",
        mbias_stats_p = f"{hsc_rep1_dir}/qc_stats/mbias_stats.p",
    shell:
        """
        mqc stats \
        --bam {input.bam} \
        --config_file {params.config_file} \
        --output_dir {params.output_dir} \
        --sample_name hsc_1 \
        --sample_meta population=HSC,rep=rep1 \
        --cores 2 \
        --motifs CG \
        {input.index_files}
        """
