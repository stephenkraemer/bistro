"""An extensible and flexible WGBS data parser built upon pysam and htslib"""

"""
fix develop install
-------------------

- if cython and/or pysam are missing, and the C files are not there (i.e. we are in the distribution): abort and say: for develop install. pysam and cython must be installed before calling setup.py
- i would not make cython and pysam build requirements. rather, by default, the 'recommended pysam version' should be used on the already cythonized ext. modules
- if cython and pysam are already present, we compile against them, because the user obviously wants to use them
- if only cython is present, I would just use the precompiled files, because the user does not enforce a pysam version - so can just take ours, and then we can take the C extension modules from us
- common case: download repo, develop. Dev. wants to user 'working' pysam and cython versions.
  - requirements.txt is associated with 'deployment config'. This is what we need here: record the config of the deployment done with this tagged version. An extern dev checking out the repo can then use the deployment config we used when testing or deploying a commit

ToDo
####

- implement these cases
	case1: pysam only installed, not cython
	if pysam version == requirements.txt version:
	    log message
            from precompiled_ext_modules
	else:
	    raise: cython needed to comply with pysam version

	case2: pysam and cython installed:    
	   re-cythonize
	   log message

	case3: no pysam version given
	   use precompiled extension modules
	   log message

- implement these cases
   case 1: pysam installed, not cython
	if pysam version == requirements.txt version:
	    log message
            from precompiled_ext_modules
	else:
	    raise
   case 2: pysam not installed
	   use precompiled extension modules
	   log message
   case 3: pysam and cython
       # need to distinguish this case from pysam presence only to allow for recythonizing code during development when cython version changes, but not pysam version
       force recythonize

- function use precompiled extension modules
  - raise if extension modules are not there.  write warning message when installing from repo: pysam and cython need to be installed before. To get versions of last release, use requirements.txt or use new versions.
- recythonize function:
  - raise if no cython available
- check logic
- commit
- test local
- test on cluster
- push to github
- push to phabricator

Done
####
- transfer cython and pysam versions to requirements.txt
- re-establish reading pysam version from requirements.txt
"""

import sys
import os.path as op

import os
from pathlib import Path

from setuptools import setup, find_packages
from setuptools.extension import Extension
from setuptools.command.build_ext import build_ext
from distutils import sysconfig
from functools import partial
from pkg_resources import get_distribution, DistributionNotFound

TESTED_PYSAM_VERSIONS=['0.11.2.2']
TESTED_CYTHON_VERSIONS=['0.26']

def no_cythonize(extensions, **_ignore):
    for extension in extensions:
        sources = []
        for sfile in extension.sources:
            path, ext = op.splitext(sfile)
            if ext in ('.pyx', '.py'):
                if extension.language == 'c++':
                    ext = '.cpp'
                else:
                    ext = '.c'
                sfile = path + ext
            sources.append(sfile)
        extension.sources[:] = sources
    return extensions


try:
    # get_distribution raises DistributionNotFound if package not installed
    if (get_distribution('cython').version not in TESTED_CYTHON_VERSIONS
        or get_distribution('pysam').version not in TESTED_PYSAM_VERSIONS):
        print('\033[1;101mYou are using untested cython and pysam version,'
              ' trying to cythonize without guarantee \033[0m')

    from Cython.Build import cythonize

    prepare_extensions = partial(cythonize, force=True)
    required_pysam = []
    print("Cython and pysam are available, cythonizing where necessary.")

except (ModuleNotFoundError, DistributionNotFound):
    print("Cython or pysam not available,"
          " using precompiled C extension modules.")
    prepare_extensions = no_cythonize
    requirements_path = Path(__file__).parent / 'requirements.txt'
    with requirements_path.open() as fin:
        required_pysam = [s.strip()
                          for s in fin.readlines()
                          if s.startswith('pysam')]
    # Test whether we are in a sdist obtained from pypi, where we have c files
    # or whether this is code retrieved from github, where we don't have c files
    # and can't install without cython
    # print(os.listdir("mqc/pileup"))
    # print(os.getcwd())
    # print(__file__)
    if not op.exists("src/mqc/pileup/bsseq_pileup_read.c"):
        raise ImportError("Can't install because cython is not available"
                          "and neither are C extension modules."
                          "Are you trying to install from the git repo without"
                          "cython? Please install through pypi if you don't"
                          "want to install cython")


python_lib = sysconfig.get_python_lib()
extensions = [
    Extension(
        "mqc.pileup.bsseq_pileup_read",
        ["src/mqc/pileup/bsseq_pileup_read.pyx"],
        include_dirs = [sysconfig.get_python_inc(prefix=sys.prefix),
                        op.join(python_lib, 'pysam'),
                        op.join(python_lib, 'pysam/include/htslib')],
    ),
]


setup(name='mqc',
      version='0.2',
      description='WGBS parser with highly configurable QC functions',
      long_description=__doc__,
      url='http://github.com/sjkv/mqc',
      author='Stephen Kraemer',
      author_email='stephenkraemer@gmail.com',
      license='MIT',

      packages=find_packages(where='src'),
      package_dir={'': 'src'},

      entry_points = {
          'console_scripts': ['mqc=mqc.cli:mqc'],
      },

      install_requires=[
          'pandas',
          'numpy',
          'seaborn',
          'pytoml==0.1.11',
          'click',
          'joblib',
          'python-magic>=0.4.13',
          'ipywidgets',
          'matplotlib',
          'plotnine',
          'rpy2',
      ] + required_pysam,

      extras_require={
          'dev': [
              'pytest',
              'pytest-mock',
              # 'cython (>=0.25)',
          ]
      },

      cmdclass = {'build_ext': build_ext},
      ext_modules=prepare_extensions(extensions),

      package_data={'mqc.resources': ['*']},
      )
