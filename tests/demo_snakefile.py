""" Demo workflow using mqc

How do I run this?
-----------------

snakemake \
--snakefile /home/kraemers/projects/mqc/tests/demo_snakefile.py \
--config pids=mpp1_rep1,mpp1_rep2 motifs=CG \
--cluster "qsub -S /bin/bash -l walltime={params.walltime},mem={params.mem},nodes=1:ppn={params.cores} -N {params.name}" \
--jobs 40 \
--jobscript /home/kraemers/projects/mqc/tests/jobscript.sh \
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

autosomes = ['19']  #[str(i) for i in range(1,20)]
other_chroms = []  #['X', 'Y', 'MT', 'phix', 'L']

config['pids'] = config['pids'].split(',')
motifs_str = config['motifs']
single_motifs = motifs_str.split('-')
motif_csv = ','.join(motifs_str.split('-'))

rule all:
    input:
        # mqc make_index
        index_files = expand("{index_dir}/GRCm38mm10_PhiX_Lambda_{motifs_str}_{chrom}.bed.gz",
                             index_dir=index_dir,
                             chrom=autosomes + other_chroms,
                             motifs_str=motifs_str),

        # mqc stats
        mbias_counts = expand("{output_rpp_dir}/{pid}/meth/qc_stats/{pid}_mbias-counts_{motifs_str}.{ext}",
                              output_rpp_dir=output_rpp_dir,
                              pid=config['pids'],
                              motifs_str=motifs_str,
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
                             motifs_str=motifs_str),

        # mqc evaluate_calls
        coverage_hist = expand("{output_rpp_dir}/{pid}/meth/qc_stats/{pid}_coverage-hist_{motifs_str}.png",
                            output_rpp_dir=output_rpp_dir,
                            pid=config['pids'],
                            motifs_str=motifs_str),

        # mqc call
        meth_calls = expand("{output_rpp_dir}/{pid}/meth/meth_calls/mcalls_{pid}_{single_motif}_{chrom}.bed.gz",
                            output_rpp_dir=output_rpp_dir,
                            pid=config['pids'],
                            single_motif=single_motifs,
                            chrom=autosomes + other_chroms),


rule make_index:
    input: f"{index_dir}/GRCm38mm10_PhiX_Lambda.fa",
    params:
        output_dir = index_dir,
        walltime = '08:00:00',
        mem = '8g',
        cores = '12',
        name = f'make_index_{motifs_str}',
        motifs_flags = '--cg' if motifs_str == 'CG' else '--cg --chg --chh',
    output:
        index_files = expand("{index_dir}/GRCm38mm10_PhiX_Lambda_{motifs_str}_{chrom}.bed.gz",
                             index_dir=index_dir,
                             chrom=autosomes + other_chroms,
                             motifs_str=motifs_str),
    shell:
        """
        mqc make_index --genome_fasta {input} \
        --output_dir {params.output_dir} \
        --cores {params.cores} \
        --seq_context 2 \
        {params.motifs_flags}
        """


rule get_stats:
    input:
        bam=f"{alignment_rpp_dir}/{{pid}}/alignment/blood_{{pid}}_merged.mdup.bam",
        index_files = expand("{index_dir}/GRCm38mm10_PhiX_Lambda_{motifs_str}_{chrom}.bed.gz",
                             index_dir=index_dir,
                             motifs_str=motifs_str,
                             chrom=autosomes),
    params:
        config_file = user_config_file,
        output_dir = f"{output_rpp_dir}/{{pid}}/meth/",
        motif_csv = motif_csv,
        walltime = '01:30:00' if motifs_str == 'CG' else '08:00:00',
        mem = '2g',
        cores = '12',
        name = f'get_stats_{{pid}}_{motifs_str}',
        sample_meta = lambda wildcards: f"population={wildcards.pid.split('_')[0]},rep={wildcards.pid.split('_')[-1]}",
    output:
        mbias_counts_p = f"{output_rpp_dir}/{{pid}}/meth/qc_stats/{{pid}}_mbias-counts_{motifs_str}.p",
        mbias_counts_tsv = f"{output_rpp_dir}/{{pid}}/meth/qc_stats/{{pid}}_mbias-counts_{motifs_str}.tsv",
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
        mbias_counts_p = f"{output_rpp_dir}/{{pid}}/meth/qc_stats/{{pid}}_mbias-counts_{motifs_str}.p",
    output:
        mbias_stats_p                = f"{output_rpp_dir}/{{pid}}/meth/qc_stats/{{pid}}_mbias-stats_{motifs_str}.p",
        mbias_stats_masked_p         = f"{output_rpp_dir}/{{pid}}/meth/qc_stats/{{pid}}_mbias-stats_masked_{motifs_str}.p",
        adjusted_cutting_sites_obj_p = f"{output_rpp_dir}/{{pid}}/meth/qc_stats/{{pid}}_adjusted_cutting_sites_obj_{motifs_str}.p",
        adjusted_cutting_sites_df_p  = f"{output_rpp_dir}/{{pid}}/meth/qc_stats/{{pid}}_adjusted_cutting_sites_df_{motifs_str}.p",
        adj_cutting_sites_plot_done  = f"{output_rpp_dir}/{{pid}}/meth/qc_stats/{{pid}}_adjusted-cutting-sites_barplot_{motifs_str}.png",
        mbias_line_plots_done        = touch(f"{output_rpp_dir}/{{pid}}/meth/qc_stats/.done_{motifs_str}_{{pid}}_mbias-line-plot"),
        freq_line_plot_done          = touch(f"{output_rpp_dir}/{{pid}}/meth/qc_stats/.done_{motifs_str}_{{pid}}_freq-line-plot"),
    params:
        config_file = user_config_file,
        output_dir = f"{output_rpp_dir}/{{pid}}/meth/",
        motif_csv = motif_csv,
        walltime = '00:30:00',
        mem = '6g',
        cores = '2',
        name = f'evalute_mbias_{{pid}}_{motifs_str}',
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
    ' --sample_name {wildcards.pid}'
    ' --sample_meta {params.sample_meta}'
    ' --use_mbias_fit'
    ' --cores {params.cores}'
    ' {input.index_files}'
)

rule call:
    input:
        bam=f"{alignment_rpp_dir}/{{pid}}/alignment/blood_{{pid}}_merged.mdup.bam",
        index_files = expand("{index_dir}/GRCm38mm10_PhiX_Lambda_{motifs_str}_{chrom}.bed.gz",
                             index_dir=index_dir,
                             motifs_str=motifs_str,
                             chrom=autosomes + other_chroms),
        adj_cut_sites_obj = f"{output_rpp_dir}/{{pid}}/meth/qc_stats/{{pid}}_adjusted_cutting_sites_obj_{motifs_str}.p",
    params:
        config_file = user_config_file,
        mem = '8g',
        cores = '12',
        output_dir = f"{output_rpp_dir}/{{pid}}/meth/",
        walltime = '42:00:00',
        name = f'mcall_{motifs_str}_{{pid}}',
        sample_meta = lambda wildcards: f"population={wildcards.pid.split('_')[0]},rep={wildcards.pid.split('_')[-1]}",
    output:
        mcall_files = expand("{output_rpp_dir}/{{pid}}/meth/meth_calls/mcalls_{{pid}}_{single_motif}_{chrom}.bed.gz",
                            output_rpp_dir=output_rpp_dir,
                            chrom=autosomes + other_chroms,
                            single_motif=single_motifs),
        coverage_counts = expand("{output_rpp_dir}/{{pid}}/meth/qc_stats/{{pid}}_coverage-counts_{motifs_str}.{ext}",
                            output_rpp_dir=output_rpp_dir,
                            motifs_str=motifs_str,
                            ext=['p', 'tsv']),
    shell: mcall_command


rule evaluate_calls:
    input:
        coverage_counts = f"{output_rpp_dir}/{{pid}}/meth/qc_stats/{{pid}}_coverage-counts_{motifs_str}.p",
    output:
        coverage_hist  = f"{output_rpp_dir}/{{pid}}/meth/qc_stats/{{pid}}_coverage-hist_{motifs_str}.png",
    params:
        config_file = user_config_file,
        output_dir = f"{output_rpp_dir}/{{pid}}/meth/",
        motif_csv = motif_csv,
        walltime = '00:30:00',
        mem = '6g',
        cores = '2',
        name = f'evaluate_calls_{{pid}}_{motifs_str}',
        sample_meta = lambda wildcards: f"population={wildcards.pid.split('_')[0]},rep={wildcards.pid.split('_')[-1]}",
    shell:
        """
        mqc evaluate_calls \
        --config_file {params.config_file} \
        --motifs {params.motif_csv} \
        --sample_name {wildcards.pid} \
        --output_dir {params.output_dir}
        """
