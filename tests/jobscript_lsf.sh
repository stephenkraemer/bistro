#!/bin/bash
# properties = {properties}
module load python/3.6.1
module load R/3.4.0
module load gcc/5.4.0
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/tbi/software/x86_64/gcc/gcc-5.4.0/el7/lib64
source /odcf/cluster/virtualenvs/kraemers/bseqtools2/bin/activate
{exec_job}
