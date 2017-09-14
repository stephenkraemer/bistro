Develop install
---------------

Install cython and pysam before attempting install of mqc, e.g.

conda install pysam cython

Install pytoml 0.1.11, version is important (latest version has bug), e.g.

pip install pytoml==0.1.11

If you get mysterious errors, try removing
mqc/src/mqc/pileup/bsseq_pileup_read.c
and
mqc/src/mqc/pileup/bsseq_pileup_read.*.so
