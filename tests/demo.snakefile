# TODO: remove these lines and test again
shell.executable("/bin/bash")
shell.prefix("module load python/3.6.0; source /home/kraemers/programs/python_virtualenvs/mqc_test/bin/activate; echo loaded mqc_test; ")

""" Demo workflow using mqc

How do I run this?
-----------------
/home/kraemers/projects/mqc/tests/demo_snakefile.py

snakemake \
--snakefile /home/stephen2/projects/mqc/tests/demo_snakefile.py \
--config pids=mpp1_rep1,mpp1_rep2 motifs=CG \
--cluster "qsub -S /bin/bash -l walltime={params.walltime},mem={params.mem},nodes=1:ppn={params.cores} -N {params.name}" \
--jobs 40 \
--jobscript /home/kraemers/projects/mqc/tests/jobscript.sh \
--dryrun

"""

import os.path as op
from pathlib import Path
import sys

script_dir = op.dirname(__file__)
sys.path.append(script_dir)
from snakefile_config import *
sys.path.remove(script_dir)

# assemble config vars
config['pids'] = config['pids'].split(',')
all_chroms = autosomes + other_chroms
reference_genome_name = Path(reference_genome).name.replace('.fa', '').replace('.gz', '')
output_dir_by_pid = f"{output_rpp_dir}/{{pid}}/meth/"
config['config_file'] = mqc_config_file
config['motifs_csv'] = motifs_csv_str

# fill out patterns with config vars, leaving only wildcards

## make_index
index_file_pattern_by_chrom = f"{index_dir}/{reference_genome_name}_{motifs_msv_str}_{chr_prefix}{{chrom}}.bed.gz"
all_index_files = expand(index_file_pattern_by_chrom, chrom = all_chroms),

## stats
mbias_counter_pattern_by_pid = f"{output_rpp_dir}/{{pid}}/meth/qc_stats/{{pid}}_mbias-counts_{motifs_msv_str}.p",
all_stats_files = expand(mbias_counter_pattern_by_pid, pid=config['pids']),

## evaluate_mbias
adj_cut_sites_obj_by_pid = f"{output_rpp_dir}/{{pid}}/meth/qc_stats/{{pid}}_adjusted_cutting_sites_obj_{motifs_msv_str}.p"
mbias_stats_patterns_by_pid = [
    adj_cut_sites_obj_by_pid,
    f"{output_rpp_dir}/{{pid}}/meth/qc_stats/{{pid}}_mbias-stats_{motifs_msv_str}.p",
    f"{output_rpp_dir}/{{pid}}/meth/qc_stats/{{pid}}_mbias-stats_masked_{motifs_msv_str}.p",
    f"{output_rpp_dir}/{{pid}}/meth/qc_stats/{{pid}}_adjusted_cutting_sites_df_{motifs_msv_str}.p",
    f"{output_rpp_dir}/{{pid}}/meth/qc_stats/{{pid}}_adjusted-cutting-sites_barplot_{motifs_msv_str}.png",
    f"{output_rpp_dir}/{{pid}}/meth/qc_stats/.done_{motifs_msv_str}_{{pid}}_mbias-line-plot",
    f"{output_rpp_dir}/{{pid}}/meth/qc_stats/.done_{motifs_msv_str}_{{pid}}_freq-line-plot",
    ]
all_evaluate_stats_files = expand(mbias_stats_patterns_by_pid, pid=config['pids']),

## call
call_file_patterns_by_pid = expand("{output_rpp_dir}/{{pid}}/meth/meth_calls/"
                                   "mcalls_{{pid}}_{single_motif}_{chrom}.bed.gz",
                                   output_rpp_dir=output_rpp_dir,
                                   single_motif=single_motifs,
                                   chrom=all_chroms)
coverage_count_patterns_by_pid = [f"{output_rpp_dir}/{{pid}}/meth/qc_stats/{{pid}}_coverage-counts_{motifs_msv_str}.tsv",
                                  f"{output_rpp_dir}/{{pid}}/meth/qc_stats/{{pid}}_coverage-counts_{motifs_msv_str}.p"]
all_mcall_files = expand(call_file_patterns_by_pid, pid=config['pids']),
all_coverage_counter_files = expand(coverage_count_patterns_by_pid, pid=config['pids']),

## evaluate_calls
evaluate_calls_patterns_by_pid = [
    f"{output_rpp_dir}/{{pid}}/meth/qc_stats/{{pid}}_coverage-hist_{motifs_msv_str}.png",
]
all_evaluate_calls_files = expand(evaluate_calls_patterns_by_pid, pid=config['pids']),

rule all:
    input:
        all_index_files,
        all_stats_files,
        # all_evaluate_stats_files,
        # all_mcall_files,
        # all_coverage_counter_files,
        # all_evaluate_calls_files,


rule make_index:
    input: reference_genome
    params:
        output_dir = index_dir,
        walltime = '08:00:00',
        mem = '8g',
        cores = '12',
        name = f'make_index_{motifs_csv_str}',
        motifs_flags = '--cg' if motifs_csv_str == 'CG' else '--cg --chg --chh',
    output:
        index_files = all_index_files,
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
        bam=bam_pattern_by_pid,
        index_files = expand(index_file_pattern_by_chrom, chrom=autosomes),
    params:
        output_dir = output_dir_by_pid,
        walltime = '01:30:00' if motifs_csv_str == 'CG' else '08:00:00',
        mem = '12g',
        cores = '12',
        name = f'get_stats_{{pid}}_{motifs_msv_str}',
        sample_meta = lambda wildcards: f"population={wildcards.pid.split('_')[0]},rep={wildcards.pid.split('_')[-1]}",
    output:
        mbias_counter = mbias_counter_pattern_by_pid,
    shell:
        """
        mqc stats \
            --bam {input.bam} \
            --config_file {config[config_file]} \
            --output_dir {params.output_dir} \
            --sample_name {wildcards.pid} \
            --sample_meta {params.sample_meta} \
            --cores {params.cores} \
            --motifs {config[motifs_csv]} \
            {input.index_files}
        """

rule evaluate_mbias:
    input:
        mbias_counter = mbias_counter_pattern_by_pid
    output:
        mbias_stats_patterns_by_pid
    params:
        output_dir = output_dir_by_pid,
        walltime = '00:30:00',
        mem = '6g',
        cores = '2',
        name = f'evalute_mbias_{{pid}}_{motifs_msv_str}',
        sample_meta = lambda wildcards: f"population={wildcards.pid.split('_')[0]},rep={wildcards.pid.split('_')[-1]}",
    shell:
        """
        mqc evaluate_mbias \
        --config_file {config[config_file]} \
        --motifs {config[motifs_csv]} \
        --sample_name {wildcards.pid} \
        --sample_meta {params.sample_meta} \
        --output_dir {params.output_dir}
        """


mcall_command = (
    'mqc call'
    ' --bam {input.bam}'
    ' --config_file {config[config_file]}'
    ' --output_dir {params.output_dir}'
    ' --sample_name {wildcards.pid}'
    ' --sample_meta {params.sample_meta}'
    ' --cores {params.cores}'
    ' {input.index_files}'
)
# ' --use_mbias_fit'

rule call:
    input:
        bam=bam_pattern_by_pid,
        index_files = all_index_files,
        # adj_cut_sites_obj = adj_cut_sites_obj_by_pid,
    params:
        mem = '8g',
        cores = '12',
        output_dir = output_dir_by_pid,
        walltime = '42:00:00',
        name = f'mcall_{motifs_msv_str}_{{pid}}',
        sample_meta = lambda wildcards: f"population={wildcards.pid.split('_')[0]},rep={wildcards.pid.split('_')[-1]}",
    output:
        mcall_files = call_file_patterns_by_pid,
        coverage_files = coverage_count_patterns_by_pid,
    shell: mcall_command


rule evaluate_calls:
    input:
        coverage_counts = coverage_count_patterns_by_pid,
    output:
        evaluate_calls_patterns_by_pid,
    params:
        output_dir = output_dir_by_pid,
        walltime = '00:30:00',
        mem = '6g',
        cores = '2',
        name = f'evaluate_calls_{{pid}}_{motifs_csv_str}',
        sample_meta = lambda wildcards: f"population={wildcards.pid.split('_')[0]},rep={wildcards.pid.split('_')[-1]}",
    shell:
        """
        mqc evaluate_calls \
        --config_file {config[config_file]} \
        --motifs {config[motifs_csv]} \
        --sample_name {wildcards.pid} \
        --output_dir {params.output_dir}
        """
