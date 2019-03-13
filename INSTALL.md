# Current situation

It is not possible to install from the custom conda channel, or from the repo using the setup.py dependency list. The problem is that bistro has not been updated to deal with depecrations introduced by updates to several of the dependencies. The setup.py and meta.yaml are not specific enough at the moment, and don't fix all critical dependencies to the allowed maximum version numbers. 

This also means that meta.yml and setup.py are not up to date with respect to the dependencies they define. Currently i only work with freezed dependencies in requirements.txt and requirements.yml, these should be up to date.

Known problems with newer versions of the dependencies:
click >= 7.0 has a problem with underscores in command names (it expects the command to contain dashes instead of underscores)
pytest >= 4.0 does not allow calling fixtures
figure_report usage in mqc is broken by the changes in the working tree, but the latest commit works

once these two problems are fixed, another problem occurs in evaluate_mbias. the command will fail and complain that itertools.compress (or so) has no len(). i assume this means that another dependency tries to call len on the generator, and did not do that in older versions. Maybe pandas?
NOTE: there is a deprecation warning in potentially relevant parts of mbias.py, when I run bistro from the current freezed environment. maybe that will point me to the problem that occurs when I try to update the bistro dependencies? 

## Install with pip and virtualenv

```
module load python/3.6.1
mkdir /odcf/cluster/virtualenvs/kraemers/bistro_bak2
virtualenv /odcf/cluster/virtualenvs/kraemers/bistro_bak2
source /odcf/cluster/virtualenvs/kraemers/bistro_bak2/bin/activate
pip install -r ~/projects/mqc/requirements.txt
pip install git+file:///home/kraemers/projects/figure_report@fb64c91
rm ~/projects/mqc/src/mqc/pileup/*.{c,so}
pip install -e ~/projects/mqc
# Latest tested commit:
# see below for test commands
```

## Install with conda

```
conda env create -f ~/projects/mqc/requirements.yml
conda activate bistro-0.2.0
pip install git+file:///home/kraemers/projects/figure_report@fb64c91
rm ~/projects/mqc/src/mqc/pileup/*.{c,so}
pip install -e ~/projects/mqc
# Latest tested commit:
# run tests (see below)
```

# Notes

This package has a dependency which is not on pip or conda: figure_report. Either install from the repo, or use the custom conda channel - it also has a figure_report package

# Standard (user) install

If cython and pysam are installed into the environment, the program
will attempt to compile the extension modules. If they are not installed,
the precompiled extension modules will be used. In this case, the required
pysam version will be automatically installed with the package.

Try to install tested cython and pysam version.
If you use untested versions, a warning will be displayed,
but compilation will still be attempted. The currently recommended cython and pysam
versions are in requirements.txt and requirements.yml

# Develop install

General
- Install cython and pysam before attempting install of mqc

If you get mysterious errors, try removing
mqc/src/mqc/pileup/bsseq_pileup_read.c
and
mqc/src/mqc/pileup/bsseq_pileup_read.*.so

to test the code, you need to install the develop dependencies (see setup.py)


# test the installation

```
pytest -v -n 12 -m 'not acceptance_test' ~/projects/mqc/tests
# expect 5 failures in tests/test_mbias.py::TestBinomPvalueBasedCuttingSiteDetermination::test_finds_gap_repair_nucleotides_in_flen_smaller_than_read_length, and 37 warnings
pytest -v -n 12 -m 'acceptance_test' ~/projects/mqc/tests
# expect no errors and 12 warnings
```
