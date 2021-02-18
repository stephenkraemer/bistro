# Introduction

Bistro is a methylation caller specialized in the detection and removal of various biases affecting the many specialized bisulfite sequencing protocols. In particular, bistro has several features for the detailed characterization and removal of non-standard M-bias in protocols such as tagmentation-based whole-genome bisulfite sequencing (T-WGBS) or PBAT-based protocols.

# Package maturity

This package is unreleased, unpublished pre-alpha software, but it is relatively mature and stable. We use it often in in-house projects, but cannot yet provide support for external users. We did not yet have a chance to test the package outside of the DKFZ HPC cluster. Also, given the pre-alpha status, we do change the API from time to time, without regard for health and safety of external users.

*Note: I am in the process of renaming the package. I started with the github repo :) The python package is currently still named mqc - apologies for the confusion.*


# Installation

First, install all dependencies. Please use the frozen dependencies detailed in requirements.yml. This is the only currently supported and tested environment. Use conda to install all dependencies.

```
conda env create -f ~/projects/mqc/requirements.yml
conda activate bistro-0.2.0
```

The package has one dependency which is not yet released as pypi or conda package: figure_report. Install it with

```
pip install git+https://github.com/stephenkraemer/figure_report
```

Finally, install bistro:

```
pip install git+https://github.com/stephenkraemer/bistro.git
```

## Notes for experienced users

This is a cython-package. If cython and pysam are installed into the environment, the program
will attempt to compile the extension modules. If they are not installed,
the precompiled extension modules will be used. In this case, the required
pysam version will be automatically installed with the package.

Try to install tested cython and pysam version. If you use untested versions, a warning will be displayed,
but compilation will still be attempted. The currently recommended cython and pysam versions are in requirements.txt and requirements.yml

*Develop install*

General
- Install cython and pysam before attempting install of mqc

If you get mysterious errors, try removing
mqc/src/mqc/pileup/bsseq_pileup_read.c
and
mqc/src/mqc/pileup/bsseq_pileup_read.*.so

to test the code, you need to install the develop dependencies (see setup.py)


# Supported operating systems

We only support Linux at the moment. The package may work on MacOS, but this is untested.

# Usage

Usage examples are provided [here](./doc/usage.ipynb)

*Note: I am in the process of renaming the package. I started with the github repo :) The python package is currently still named mqc - apologies for the confusion.*

We are working on providing an updated sphinx-based documentation, bear with us. Some rst files are already in doc, but they are not yet complete.
