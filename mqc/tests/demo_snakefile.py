""" Demo workflow using mqc

How do I run this?
-----------------

snakemake \
--snakefile /home/kraemers/projects/mqc/mqc/tests/demo_snakefile.py \
--config pids=mpp1_rep1,mpp1_rep2 \
--cluster "qsub -S /bin/bash -l walltime={params.walltime},mem={params.mem},nodes=1:ppn={params.cores} -N {params.name}" \
--jobs 40 \
--jobscript /home/kraemers/projects/mqc/mqc/tests/jobscript.sh \
--dryrun

"""
# TODO: remove these lines and test again
shell.executable("/bin/bash")
shell.prefix("module load python/3.6.0; source /home/kraemers/programs/python_virtualenvs/mqc_test/bin/activate; echo loaded mqc_test; ")

import os.path as op

sandbox_dir = "/icgc/dkfzlsdf/analysis/B080/kraemers/projects/mcall_qc/sandbox"
output_rpp_dir = f"{sandbox_dir}/results_per_pid"
index_dir = f"{sandbox_dir}/genomes/GRCm38mm10_PhiX_Lambda"
user_config_file = f"{sandbox_dir}/user_config.toml"
alignment_rpp_dir = "/icgc/dkfzlsdf/analysis/hs_ontogeny/results/wgbs/results_per_pid"

# autosomes = [str(i) for i in range(1,20)]
# other_chroms = ['X', 'Y', 'MT', 'phix', 'L']

autosomes = [str(i) for i in range(19,20)]
other_chroms = ['L']

config['pids'] = config['pids'].split(',')
motifs_str_list = ['CG', 'CG-CHG-CHH']

rule all:
    input:
        # mqc make_index
        index_files = expand("{index_dir}/GRCm38mm10_PhiX_Lambda_{motifs_str}_{chrom}.bed.gz",
                             index_dir=index_dir,
                             chrom=autosomes + other_chroms,
                             motifs_str=motifs_str_list),

        # mqc stats
        mbias_counts = expand("{output_rpp_dir}/{pid}/meth/qc_stats/{pid}_mbias-counts_{motifs_str}.{ext}",
                              output_rpp_dir=output_rpp_dir,
                              pid=config['pids'],
                              motifs_str=motifs_str_list,
                              ext=['p', 'tsv']),

        # mqc evaluate_mbias
        mbias_stats = expand(["{output_rpp_dir}/{pid}/meth/qc_stats/{pid}_mbias-stats_{motifs_str}.p",
                              "{output_rpp_dir}/{pid}/meth/qc_stats/{pid}_mbias-stats_masked_{motifs_str}.p",
                              "{output_rpp_dir}/{pid}/meth/qc_stats/{pid}_adjusted_cutting_sites_obj_{motifs_str}.p",
                              "{output_rpp_dir}/{pid}/meth/qc_stats/{pid}_adjusted_cutting_sites_df_{motifs_str}.p",
                              "{output_rpp_dir}/{pid}/meth/qc_stats/{pid}_adjusted-cutting-sites_barplot_{motifs_str}.png",
                              "{output_rpp_dir}/{pid}/meth/qc_stats/.done_{motifs_str}_{pid}_mbias-line-plot",
                              "{output_rpp_dir}/{pid}/meth/qc_stats/.done_{motifs_str}_{pid}_freq-line-plot",],
                             output_rpp_dir=output_rpp_dir,
                             pid=config['pids'],
                             motifs_str=motifs_str_list),

        # mqc call
        meth_calls = expand("{output_rpp_dir}/{pid}/meth/meth_calls/mcalls_{pid}_{single_motif}_{chrom}.bed.gz",
                            output_rpp_dir=output_rpp_dir,
                            pid=config['pids'],
                            single_motif=['CG', 'CHG', 'CHH'],
                            chrom=autosomes + other_chroms),

        meth_calls_cg_only = expand("{output_rpp_dir}/{pid}_cg-only/meth/meth_calls/mcalls_{pid}_cg-only_CG_{chrom}.bed.gz",
                                    output_rpp_dir=output_rpp_dir,
                                    pid=config['pids'],
                                    chrom=autosomes + other_chroms),

rule make_index:
    input: f"{index_dir}/GRCm38mm10_PhiX_Lambda.fa",
    params:
        output_dir = index_dir,
        walltime = '08:00:00',
        mem = '8g',
        cores = '12',
        name = 'make_index_{motifs_str}',
        motifs_flags = lambda wildcards: '--cg' if wildcards.motifs_str == 'CG' else '--cg --chg --chh',
    output:
        index_files = expand("{index_dir}/GRCm38mm10_PhiX_Lambda_{{motifs_str}}_{chrom}.bed.gz",
                             index_dir=index_dir,
                             chrom=autosomes + other_chroms),
    shell:
        """
        mqc make_index --genome_fasta {input} \
        --output_dir {params.output_dir} \
        --cores {params.cores} \
        {params.motifs_flags}
        """


rule get_stats:
    input:
        bam=f"{alignment_rpp_dir}/{{pid}}/alignment/blood_{{pid}}_merged.mdup.bam",
        index_files = expand("{index_dir}/GRCm38mm10_PhiX_Lambda_{{motifs_str}}_{chrom}.bed.gz",
                             index_dir=index_dir,
                             chrom=autosomes),
    params:
        config_file = user_config_file,
        output_dir = f"{output_rpp_dir}/{{pid}}/meth/",
        motif_csv = lambda wildcards: ','.join(wildcards.motifs_str.split('-')),
        walltime = lambda wildcards: '01:30:00' if wildcards.motifs_str == 'CG' else '08:00:00',
        mem = '2g',
        cores = '12',
        name = 'get_stats_{pid}_{motifs_str}',
        sample_meta = lambda wildcards: f"population={wildcards.pid.split('_')[0]},rep={wildcards.pid.split('_')[-1]}",
    output:
        mbias_counts_p = f"{output_rpp_dir}/{{pid}}/meth/qc_stats/{{pid}}_mbias-counts_{{motifs_str}}.p",
        mbias_counts_tsv = f"{output_rpp_dir}/{{pid}}/meth/qc_stats/{{pid}}_mbias-counts_{{motifs_str}}.tsv",
    shell:
        """
        mqc stats \
            --bam {input.bam} \
            --config_file {params.config_file} \
            --output_dir {params.output_dir} \
            --sample_name {wildcards.pid} \
            --sample_meta {params.sample_meta} \
            --cores {params.cores} \
            --motifs {params.motif_csv} \
            {input.index_files}
        """



rule evaluate_mbias:
    input:
        mbias_counts_p = f"{output_rpp_dir}/{{pid}}/meth/qc_stats/{{pid}}_mbias-counts_{{motifs_str}}.p",
    output:
        mbias_stats_p                = f"{output_rpp_dir}/{{pid}}/meth/qc_stats/{{pid}}_mbias-stats_{{motifs_str}}.p",
        mbias_stats_masked_p         = f"{output_rpp_dir}/{{pid}}/meth/qc_stats/{{pid}}_mbias-stats_masked_{{motifs_str}}.p",
        adjusted_cutting_sites_obj_p = f"{output_rpp_dir}/{{pid}}/meth/qc_stats/{{pid}}_adjusted_cutting_sites_obj_{{motifs_str}}.p",
        adjusted_cutting_sites_df_p  = f"{output_rpp_dir}/{{pid}}/meth/qc_stats/{{pid}}_adjusted_cutting_sites_df_{{motifs_str}}.p",
        adj_cutting_sites_plot_done  = f"{output_rpp_dir}/{{pid}}/meth/qc_stats/{{pid}}_adjusted-cutting-sites_barplot_{{motifs_str}}.png",
        mbias_line_plots_done        = touch(f"{output_rpp_dir}/{{pid}}/meth/qc_stats/.done_{{motifs_str}}_{{pid}}_mbias-line-plot"),
        freq_line_plot_done          = touch(f"{output_rpp_dir}/{{pid}}/meth/qc_stats/.done_{{motifs_str}}_{{pid}}_freq-line-plot"),
    params:
        config_file = user_config_file,
        output_dir = f"{output_rpp_dir}/{{pid}}/meth/",
        motif_csv = lambda wildcards: ','.join(wildcards.motifs_str.split('-')),
        walltime = '00:30:00',
        mem = '6g',
        cores = '2',
        name = 'evalute_mbias_{pid}_{motifs_str}',
        sample_meta = lambda wildcards: f"population={wildcards.pid.split('_')[0]},rep={wildcards.pid.split('_')[-1]}",
    shell:
        """
        mqc evaluate_mbias \
        --config_file {params.config_file} \
        --motifs {params.motif_csv} \
        --sample_name {wildcards.pid} \
        --sample_meta {params.sample_meta} \
        --output_dir {params.output_dir}
        """



mcall_command = (
    'mqc call'
    ' --bam {input.bam}'
    ' --config_file {params.config_file}'
    ' --output_dir {params.output_dir}'
    ' --sample_name {params.pid}'
    ' --sample_meta {params.sample_meta}'
    ' --cores {params.cores}'
    ' {input.index_files}'
)

rule call_meth_for_CG:
    input:
        bam=f"{alignment_rpp_dir}/{{pid}}/alignment/blood_{{pid}}_merged.mdup.bam",
        index_files = expand("{index_dir}/GRCm38mm10_PhiX_Lambda_CG_{chrom}.bed.gz",
                             index_dir=index_dir,
                             chrom=autosomes + other_chroms)
    params:
        config_file = user_config_file,
        mem = '8g',
        cores = '12',
        walltime = '04:00:00',
        output_dir = f"{output_rpp_dir}/{{pid}}_cg-only/meth/",
        name = 'mcall_CG_{pid}_cg-only',
        sample_meta = lambda wildcards: f"population={wildcards.pid.split('_')[0]},rep={wildcards.pid.split('_')[-1]}",
        pid = "{pid}_cg-only",
    output:
        expand("{output_rpp_dir}/{{pid}}_cg-only/meth/meth_calls/mcalls_{{pid}}_cg-only_CG_{chrom}.bed.gz",
               output_rpp_dir=output_rpp_dir,
               chrom=autosomes + other_chroms),
    shell: mcall_command


rule call_meth_for_all_motifs:
    input:
        bam=f"{alignment_rpp_dir}/{{pid}}/alignment/blood_{{pid}}_merged.mdup.bam",
        index_files = expand("{index_dir}/GRCm38mm10_PhiX_Lambda_CG-CHG-CHH_{chrom}.bed.gz",
                             index_dir=index_dir,
                             chrom=autosomes + other_chroms)
    params:
        config_file = user_config_file,
        mem = '8g',
        cores = '12',
        output_dir = f"{output_rpp_dir}/{{pid}}/meth/",
        walltime = '42:00:00',
        name = 'mcall_CG-CHG-CHH_{pid}',
        sample_meta = lambda wildcards: f"population={wildcards.pid.split('_')[0]},rep={wildcards.pid.split('_')[-1]}",
        pid = "{pid}",
    output:
        expand("{output_rpp_dir}/{{pid}}/meth/meth_calls/mcalls_{{pid}}_{single_motif}_{chrom}.bed.gz",
               output_rpp_dir=output_rpp_dir,
               chrom=autosomes + other_chroms,
               single_motif=['CG', 'CHG', 'CHH']),
    shell: mcall_command

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
