# load module
# load virtualenv
# export PYTHONPATH=$PYTHONPATH:/home/kraemers/projects/mqc
# cd /home/kraemers/projects/mqc/mqc

cython bsseq_pileup_read.pyx

gcc -shared -pthread -fPIC -fwrapv -O2 -Wall -fno-strict-aliasing \
-I /home/kraemers/programs/python_virtualenvs/python_3.6.0/include/python3.6m \
-I /home/kraemers/programs/python_virtualenvs/python_3.6.0/lib/python3.6/site-packages/pysam/ \
-I /home/kraemers/programs/python_virtualenvs/python_3.6.0/lib/python3.6/site-packages/pysam/include/htslib \
-o bsseq_pileup_read.so bsseq_pileup_read.c

echo $?
