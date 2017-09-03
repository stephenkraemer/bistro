#!/bin/bash
# properties = {properties}
module load python/3.6.0 && echo 'loaded module python'
source /home/kraemers/programs/python_virtualenvs/mqc_test/bin/activate && echo 'sourced virtualenv'
{exec_job}
