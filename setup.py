"""An extensible and flexible WGBS data parser built upon pysam and htslib"""

import sys
import os.path as op

import os
from pathlib import Path

from setuptools import setup, find_packages
from setuptools.extension import Extension
from setuptools.command.build_ext import build_ext
from distutils import sysconfig


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
    from Cython.Build import cythonize
    print("Cython is available, cythonizing where necessary.")
    prepare_extensions = cythonize
    # will proceed to compile pyx to C
    required_pysam_str = 'pysam'
except ModuleNotFoundError:
    print("Cython is not available, using precompiled C extension modules.")
    prepare_extensions = no_cythonize
    requirements_path = Path(__file__).parent / 'requirements.txt'
    with requirements_path.open() as fin:
        required_pysam_str = [s.strip()
                              for s in fin.readlines()
                              if s.startswith('pysam')][0]
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
          required_pysam_str,
          'numpy',
          'seaborn',
          'pytoml==0.1.11',
          'click',
          'joblib',
          'python-magic (>=0.4.13)',
          'pysam',
          'ipywidgets',
          'matplotlib',
      ],

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
