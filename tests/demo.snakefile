""" Demo workflow using mqc

How do I run this?
-----------------
/home/kraemers/projects/mqc/tests/demo_snakefile.py

cd /icgc/dkfzlsdf/project/mouse_hematopoiesis/sequencing/whole_genome_bisulfite_tagmentation_sequencing/view-by-pid
pids_csv=$(echo mpp?_? hsc_? | tr ' ' ,)
echo $pids_csv
cd

pids_csv='monos_3,pdc_3,mep_3'

snakemake \
--snakefile ~/projects/mqc/tests/demo.snakefile \
--configfile /home/kraemers/projects/mqc/tests/snakefile_config.json \
--config pids=$pids_csv motifs=CG \
--jobs 1000 \
--jobscript /home/kraemers/projects/mqc/tests/jobscript_lsf.sh \
--cluster "bsub -R rusage[mem={params.mem}] -M {params.mem} -n {params.cores} -J {params.name} -W {params.walltime} -o /home/kraemers/temp/logs/" \
--dryrun

--keep-going \
--forcerun get_stats \

--forcerun get_stats \
--cluster "qsub -S /bin/bash -l walltime={params.walltime},mem={params.mem}g,nodes=1:ppn={params.cores} -N {params.name}" \

"""

import os.path as op
import json
import pickle
from collections import defaultdict
from pathlib import Path
import sys
import re

# with open('/home/kraemers/projects/mqc/tests/snakefile_config.json') as fin:
#     config = json.load(fin)

output_rpp_dir = config['output_rpp_dir']
index_dir = config['index_dir']
reference_genome = config['reference_genome']
mqc_config_file = config['mqc_config_file']
in_rpp_dir = config['in_rpp_dir']
bam_pattern_by_pid = config['bam_pattern_by_pid']
autosomes = config['autosomes']
other_chroms = config['other_chroms']
single_motifs = config['single_motifs']
motifs_msv_str = config['motifs_msv_str']
motifs_csv_str = config['motifs_csv_str']
chr_prefix = config['chr_prefix']

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
all_stats_files = expand(mbias_counter_pattern_by_pid, pid=config['pids'])

## evaluate_mbias
adj_cut_sites_obj_by_pid = f"{output_rpp_dir}/{{pid}}/meth/qc_stats/{{pid}}_adjusted_cutting_sites_obj_{motifs_msv_str}.p"
adj_cut_sites_df_by_pid = f"{output_rpp_dir}/{{pid}}/meth/qc_stats/{{pid}}_adjusted_cutting_sites_df_{motifs_msv_str}.p"
full_mbias_stats_pattern_by_pid = f"{output_rpp_dir}/{{pid}}/meth/qc_stats/{{pid}}_mbias-stats_{motifs_msv_str}.p"
trimmed_mbias_stats_pattern_by_pid = f"{output_rpp_dir}/{{pid}}/meth/qc_stats/{{pid}}_mbias-stats_masked_{motifs_msv_str}.p"
full_phredfiltered_mbias_stats_by_pid = f"{output_rpp_dir}/{{pid}}/meth/qc_stats/{{pid}}_mbias-stats_phred-threshold_{motifs_msv_str}.p.p"
trimmed_phredfiltered_mbias_stats_by_pid = f"{output_rpp_dir}/{{pid}}/meth/qc_stats/{{pid}}_mbias-stats_masked_phred-threshold_{motifs_msv_str}.p.p"
all_evaluate_stats_files = expand(trimmed_mbias_stats_pattern_by_pid, pid=config['pids'])

## plot mbias
test_mbias_plot_config_json = op.expanduser('~/projects/mqc/src/mqc/resources/mbias_plots_config.json')
mbias_plot_done_by_pid = f"{output_rpp_dir}/{{pid}}/meth/qc_stats/mbias_stats/plots.done"
## call
mcall_bed_patterns_by_pid = expand("{output_rpp_dir}/{{pid}}/meth/meth_calls/"
                                   "mcalls_{{pid}}_{single_motif}_{chrom}.bed.gz",
                                   output_rpp_dir=output_rpp_dir,
                                   single_motif=single_motifs,
                                   chrom=all_chroms)

mcall_bismark_patterns_by_pid = expand("{output_rpp_dir}/{{pid}}/meth/meth_calls/"
                                       "{{pid}}_{single_motif}_chr-{chrom}.bismark.gz",
                                       output_rpp_dir=output_rpp_dir,
                                       single_motif=single_motifs,
                                       chrom=all_chroms)

mcall_stratifiedbed_patterns_by_pid = expand("{output_rpp_dir}/{{pid}}/meth/meth_calls/"
                                             "{{pid}}_{single_motif}_chr-{chrom}_stratified.bed.gz",
                                             output_rpp_dir=output_rpp_dir,
                                             single_motif=single_motifs,
                                             chrom=all_chroms)


coverage_count_patterns_by_pid = [f"{output_rpp_dir}/{{pid}}/meth/qc_stats/{{pid}}_coverage-counts_{motifs_msv_str}.tsv",
                                  f"{output_rpp_dir}/{{pid}}/meth/qc_stats/{{pid}}_coverage-counts_{motifs_msv_str}.p"]
all_mcall_files = expand(mcall_bed_patterns_by_pid + mcall_bismark_patterns_by_pid,
                         pid=config['pids']),
all_coverage_counter_files = expand(coverage_count_patterns_by_pid, pid=config['pids']),

## evaluate_calls
evaluate_calls_patterns_by_pid = [
    f"{output_rpp_dir}/{{pid}}/meth/qc_stats/{{pid}}_coverage-hist_{motifs_msv_str}.png",
]
all_evaluate_calls_files = expand(evaluate_calls_patterns_by_pid, pid=config['pids']),

def get_sample_metadata(wildcards) -> str:
    """Extract sample name and replicate ID from sample ID

    This function assumes the standard mouse hematopoiesis project
    nomenclature

    Returns:
        String with comma-separated metadata fields
        sample=<sample_name>,rep=<repID>
    """
    # example: hsc_mu-d3a-ki_1
    sample_name_pattern = r'^(.*)_(\d+)$'
    try:
        population, rep = re.search(sample_name_pattern, wildcards.pid).groups()
    except KeyError:
        raise ValueError(f"Can't extract sample metadata from pid {wildcards.pid}")
    return (f"population={population},"
            f"rep={rep}")

with open('/home/kraemers/projects/mqc/tests/test_files/hs_data_read_length_mapping.p', 'rb') as fin:
  pid_to_read_length_mapping = pickle.load(fin)  # Dict[str, int], e.g. {'t-cells_1': 101}
# originally created with:
# /home/stephen/projects/mouse_hematopoiesis/my_docs/labbook/attachments/sample_id_to_read_length_mapping.ipynb
# Add samples which came later
new_read_lengths = {
    'hsc-aged_10': 125,
    'hsc-aged_11': 125,
    'megas_4': 125,
    'megas_5': 125,
    'megas_6': 125,
    'hsc_tr-pbs_5': 125,
    'lsk_mu-d3a-r882h_tr-aza_4': 125,
    'lsk_mu-d3a-r882h_tr-aza_5': 125,
    'lsk_mu-d3a-r882h_tr-aza_7': 125,
    'lsk_mu-d3a-r882h_tr-none_5': 125,
    'lsk_mu-d3a-r882h_tr-none_6': 125,
    'lsk_mu-d3a-r882h_tr-none_7': 125,
    'lsk_mu-d3a-r882h_tr-none_8': 125,
    'lsk_mu-d3a-wt_tr-aza_4': 125,
    'lsk_mu-d3a-wt_tr-aza_5': 125,
    'lsk_mu-d3a-wt_tr-aza_6': 125,
    'lsk_mu-d3a-wt_tr-aza_7': 125,
    'lsk_mu-d3a-wt_tr-none_4': 125,
    'lsk_mu-d3a-wt_tr-none_5': 125,
    'lsk_mu-d3a-wt_tr-none_6': 125,
    'lsk_mu-d3a-wt_tr-none_7': 125,
    '5N_pbat': 125,
    '5N_pbat_reprocessed': 125,
    '5N_tagmentation': 125,
    '5N_xten': 150,
    '5T_pbat': 125,
    '5T_pbat_reprocessed': 125,
    '5T_tagmentation': 125,
    '5T_xten': 150,
    '6N_pbat': 125,
    '6N_pbat_reprocessed': 125,
    '6N_tagmentation': 125,
    '6N_xten': 150,
    '6T_pbat': 125,
    '6T_pbat_reprocessed': 125,
    '6T_tagmentation': 125,
    '6T_xten': 150,
}
pid_to_read_length_mapping.update(new_read_lengths)

rule all:
    input:
        all_index_files,
        all_stats_files,
        all_evaluate_stats_files,
        expand(mbias_plot_done_by_pid, pid=config['pids']),
        all_mcall_files,
        # all_coverage_counter_files,
        # all_evaluate_calls_files,


rule make_index:
    input: reference_genome
    params:
        output_dir = index_dir,
        walltime = '08:00',
        mem = 18000,
        cores = '25',
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
        walltime = '01:30' if motifs_csv_str == 'CG' else '08:00',
        mem = 12000,
        cores = '8',
        name = f'get_stats_{{pid}}_{motifs_msv_str}',
        sample_meta = get_sample_metadata,
        read_length = lambda wildcards: pid_to_read_length_mapping[wildcards.pid],
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
            --read_length {params.read_length} \
            {input.index_files}
        """


rule evaluate_mbias:
    input:
        mbias_counter = mbias_counter_pattern_by_pid
    output:
        # full_mbias_stats_pattern_by_pid,
        trimmed_mbias_stats_pattern_by_pid,
        full_phredfiltered_mbias_stats_by_pid,
        trimmed_phredfiltered_mbias_stats_by_pid,
        adj_cut_sites_df_by_pid,
    params:
        output_dir = output_dir_by_pid,
        walltime = '01:00',
        mem = 40000,
        cores = '8',
        name = f'evalute_mbias_{{pid}}_{motifs_msv_str}',
        sample_meta = get_sample_metadata,
        plateau_detection_params = json.dumps(
            {"algorithm": "binomp",
             "allow_slope": True,
             # "min_plateau_length": 30,
             # "max_slope": 0.001,
             # "plateau_flen": 130,
             # "plateau_bs_strands": ["w_bc", "c_bc"],
             # "always_accept_distance_from_plateau": 0.02,
            })
    shell:
        """
        mqc evaluate_mbias \
        --config_file {config[config_file]} \
        --motifs {config[motifs_csv]} \
        --sample_name {wildcards.pid} \
        --sample_meta {params.sample_meta} \
        --output_dir {params.output_dir} \
        --plateau_detection '{params.plateau_detection_params}' \
        --use_cached_mbias_stats
        """

def get_datasets_str(wildcards, input):
    full_mbias_stats = full_mbias_stats_pattern_by_pid.format(pid=wildcards.pid)
    res =(f'full={full_mbias_stats},trimmed={input.trimmed_mbias_stats},'
          f'full_phred-threshold={input.full_phredfiltered_mbias_stats},'
          f'trimmed_phred-threshold={input.trimmed_phredfiltered_mbias_stats}')
    return res

# Use default mbias plot config
rule plot_mbias:
    input:
        # full_mbias_stats=full_mbias_stats_pattern_by_pid,
        trimmed_mbias_stats=trimmed_mbias_stats_pattern_by_pid,
        full_phredfiltered_mbias_stats=full_phredfiltered_mbias_stats_by_pid,
        trimmed_phredfiltered_mbias_stats=trimmed_phredfiltered_mbias_stats_by_pid,
        mbias_plot_config=test_mbias_plot_config_json,
    output:
        touch(mbias_plot_done_by_pid)
    params:
        output_dir = output_dir_by_pid,
        walltime = '00:25',
        mem = 40000,
        cores = '1',
        name = f'plot_mbias_{{pid}}_{motifs_msv_str}',
        sample_meta = get_sample_metadata,
        datasets = get_datasets_str
    shell:
        """
        mqc mbias_plots \
        --output_dir {params.output_dir} \
        --datasets {params.datasets} \
        """


mcall_command = (
    'mqc call'
    ' --bam {input.bam}'
    ' --config_file {config[config_file]}'
    ' --output_dir {params.output_dir}'
    ' --sample_name {wildcards.pid}'
    ' --sample_meta {params.sample_meta}'
    ' --cores {params.cores}'
    ' --output_formats bed,bismark,stratified_bed'
    ' --trimming cutting_sites::{input.cutting_sites_df}'
    ' --max_read_length {params.read_length}'
    ' {input.index_files}'
)
rule call:
    input:
        bam=bam_pattern_by_pid,
        index_files = all_index_files,
        cutting_sites_df = adj_cut_sites_df_by_pid,
    params:
        mem = 26000,
        cores = '12',
        output_dir = output_dir_by_pid,
        walltime = '8:00',
        name = f'mcall_{motifs_msv_str}_{{pid}}',
        sample_meta = get_sample_metadata,
        read_length = lambda wildcards: pid_to_read_length_mapping[wildcards.pid],
    output:
        mcall_bed_patterns_by_pid,
        mcall_bismark_patterns_by_pid,
        mcall_stratifiedbed_patterns_by_pid,
        # coverage_files = coverage_count_patterns_by_pid,
    shell: mcall_command


rule evaluate_calls:
    input:
        coverage_counts = coverage_count_patterns_by_pid,
    output:
        evaluate_calls_patterns_by_pid,
    params:
        output_dir = output_dir_by_pid,
        walltime = '00:30:00',
        mem = 6000,
        cores = '2',
        name = f'evaluate_calls_{{pid}}_{motifs_csv_str}',
        sample_meta = get_sample_metadata
    shell:
        """
        mqc evaluate_calls \
        --config_file {config[config_file]} \
        --motifs {config[motifs_csv]} \
        --sample_name {wildcards.pid} \
        --output_dir {params.output_dir}
        """
