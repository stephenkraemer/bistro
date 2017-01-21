""" Benchmark the methylation calling and qc_run functionality of mqc

ToDo
- command line interface to select tests (methylation calling, qc_run or both)
- add methylation calling tests
  - may be easier if the index files are presented in a table format?
- wait for jobs to finish and merge results?
"""
import textwrap
import time
import shutil
import os
import subprocess

benchmark_dir = '/home/kraemers/temp/mqc_benchmarks/'
os.makedirs(benchmark_dir, exist_ok=True)

timestamp = time.strftime('%m_%d_%H-%M-%S')
run_folder = os.path.join(benchmark_dir, timestamp)
os.makedirs(run_folder)
os.makedirs(os.path.join(run_folder, 'qc_run'))
os.makedirs(os.path.join(run_folder, 'mcall'))
code_snapshot_dir = os.path.join(run_folder, 'mqc')
shutil.copytree('/home/kraemers/projects/mqc/mqc/', dst=code_snapshot_dir)


config_file = '/home/kraemers/projects/mqc/mqc/config.default.toml'

index_files = [
    ('/icgc/dkfzlsdf/analysis/B080/kraemers/projects/mcall_qc/mcall_qc_test_data/indices/chr11.10_000.cg.bed.gz',  10_000, '01:00:00'),
    ('/icgc/dkfzlsdf/analysis/B080/kraemers/projects/mcall_qc/mcall_qc_test_data/indices/chr11.100_000.cg.bed.gz', 100_000, '01:00:00'),
    ('/icgc/dkfzlsdf/analysis/B080/kraemers/projects/mcall_qc/mcall_qc_test_data/indices/chr11.cg.bed.gz',         772_956, '02:00:00'),
]

bam_files = [
    ('/icgc/dkfzlsdf/analysis/hs_ontogeny/results/wgbs/results_per_pid/hsc_rep1/alignment/blood_hsc_rep1_merged.mdup.bam', 'hsc_rep1'),
]

for idx_file, n_pos_in_index, walltime in index_files:
    for bam_file, sample_name in bam_files:

        job_name = sample_name + '_' + str(n_pos_in_index)
        output_dir_job_results = os.path.join(run_folder, 'qc_run', job_name)
        os.makedirs(output_dir_job_results, exist_ok=True)

        script_path = os.path.join(output_dir_job_results, 'cluster_commands.sh')
        with open(script_path, 'wt') as f:
            bash_code = textwrap.dedent(f'''\
            module load python/3.6.0
            source ~/programs/python_virtualenvs/python_3.6.0/bin/activate
            export PYTHONPATH=$PYTHONPATH:{run_folder}
            cd {code_snapshot_dir}
            START=$(date +%s.%N)
            python3 cli.py qc_run \\
            --index {idx_file} \\
            --bam {bam_file} \\
            --output_dir {output_dir_job_results} \\
            --sample_name {sample_name} \\
            --config_file {config_file}
            END=$(date +%s.%N)
            diff=$(echo "$END - $START" | bc)
            echo "run time for {sample_name} with {n_pos_in_index} index positions:\\n$diff"''')
            f.write(bash_code)

        out_log_file = os.path.join(run_folder, 'qc_run', job_name + '.o.log')
        err_log_file = os.path.join(run_folder, 'qc_run', job_name + '.e.log')
        command_list = ['qsub', '-N', job_name ,'-o', out_log_file,
                        '-e', err_log_file, '-l', f'walltime={walltime},nodes=1:ppn=1',
                        script_path]
        subprocess.check_call(command_list)
