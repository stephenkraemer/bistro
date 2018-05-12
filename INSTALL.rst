User install
------------
If cython and pysam are installed into the environment, the program
will attempt to compile the extension modules. If they are not installed,
the precompiled extension modules will be used. In this case, the required
pysam version will be automatically installed with the package.

Try to install tested cython and pysam version.
If you use untested versions, a warning will be displayed,
but compilation will still be attempted. The currently recommended cython and pysam
versions are in requirements.txt

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
