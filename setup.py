"""Install mqc"""

import sys
import os.path as op

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
      zip_safe=False,

      entry_points = {
          'console_scripts': ['mqc=mqc.cli:mqc'],
      },

      install_requires=[
          # 'rpy2',  only for mbias plots if plotting with R ggplot2
          'altair',
          'click',
          'dataclasses',
          # 'dpcontracts',  (github)
          'feather-format',
          'figure_report',
          'joblib',
          'matplotlib',
          'more_itertools',
          'numpy',
          'pandas',
          'plotnine',
          'python-magic>=0.4.13',
          'scikit-learn',
          'scipy',
          'seaborn',
          'statsmodels',
          'toml',
          'toolz',
      ] + required_pysam,

      extras_require={
          'dev': [
              'pytest',
              'pytest-mock',
              'pytest-xdist',
              # 'cython (>=0.25)',
              'mypy==0.610',
              'numpy-stubs',
          ],
      },

      dependency_links=[
          'git+https://github.com/numpy/numpy-stubs.git#egg=numpy-stubs-0.01',
      ],

      cmdclass = {'build_ext': build_ext},
      ext_modules=prepare_extensions(extensions),

      package_data={'mqc.resources': ['*'],
                    'mqc': ['py.typed']},
      )
