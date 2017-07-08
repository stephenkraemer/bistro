#!/usr/bin/env bash

idx_file=/icgc/dkfzlsdf/analysis/B080/kraemers/projects/mcall_qc/mcall_qc_test_data/indices/chr11.10_000.cg.bed.gz
#idx_file=/icgc/dkfzlsdf/analysis/B080/kraemers/projects/mcall_qc/mcall_qc_test_data/indices/chr11.cg.bed.gz
bam_file=/icgc/dkfzlsdf/analysis/hs_ontogeny/results/wgbs/results_per_pid/hsc_rep1/alignment/blood_hsc_rep1_merged.mdup.bam
output_dir_job_results=/home/kraemers/temp/qc_run_interactive_test_results/
sample_name='hsc_rep1'
config_file=/home/kraemers/projects/mqc/mqc/config.default.toml

module load python/3.6.0
source ~/programs/python_virtualenvs/python_3.6.0/bin/activate
export PYTHONPATH=${PYTHONPATH}:/home/kraemers/projects/mqc/

cd /home/kraemers/projects/mqc/mqc/

#python3 cli.py qc_run \
python3 -m pdb cli.py qc_run \
--index ${idx_file} \
--bam ${bam_file} \
--output_dir ${output_dir_job_results} \
--sample_name ${sample_name} \
--config_file ${config_file}
