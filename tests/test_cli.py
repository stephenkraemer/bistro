import os
import mqc
import textwrap
import subprocess

from mqc.tests.conftest import bam_path, index_file_path, sample_name
test_output_dir = mqc.tests.conftest.test_output_dir()

def test_cli_local():
    script_path = os.path.join(test_output_dir, 'cli_script.sh')
    with open(script_path, 'wt') as f:
        bash_code = textwrap.dedent(f'''\
                source ~/programs/python_virtualenvs/python3.6/bin/activate
                export PYTHONPATH=$PYTHONPATH:/home/kraemers/project/mqc
                cd /home/kraemers/projects/mqc/mqc
                python3 cli.py qc_run \\
                --index {index_file_path} \\
                --bam {bam_path} \\
                --output_dir {test_output_dir} \\
                --sample_name {sample_name} \\
                --sample_meta population=hsc,replicate=1 \\
                --config_file /home/kraemers/projects/mqc/mqc/config.default.toml''')
        f.write(bash_code)

    exit_code = subprocess.call(['bash', script_path])
    assert exit_code == 0
