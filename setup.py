from os.path import join
from setuptools import setup, find_packages
from setuptools.extension import Extension
from Cython.Build import cythonize
from distutils import sysconfig
import sys

python_lib = sysconfig.get_python_lib()
extensions = [
    Extension(
        "mqc.pileup.bsseq_pileup_read",
        ["mqc/pileup/bsseq_pileup_read.pyx"],
        include_dirs = [sysconfig.get_python_inc(prefix=sys.prefix),
                        join(python_lib, 'pysam'),
                        join(python_lib, 'pysam/include/htslib')],
        # libraries=['fftw3', 'fftw3f', 'fftw3l', 'fftw3_threads', 'fftw3f_threads', 'fftw3l_threads'],
        # library_dirs=['/some/path/to/include/'], # not needed for fftw unless it is installed in an unusual place
    ),
]

setup(name='mqc',
      version='0.2',
      description='WGBS parser with highly configurable QC functions',
      url='http://github.com/sjkv/mqc',
      author='Stephen Kraemer',
      author_email='stephenkraemer@gmail.com',
      license='MIT',
      packages=find_packages(),

      entry_points = {
          'console_scripts': ['mqc=mqc.cli:mqc'],
      },
      install_requires=[
          'cython (>=0.25)',
          'pandas',
          'numpy',
          'seaborn',
          'pytoml',
          'click',
          'joblib',
          'python-magic (>=0.4.13)',
          'pytest',
          'pytest-mock',
          'pysam',
          'ipywidgets',
          'matplotlib',
          'python-magic'
      ],
      ext_modules=cythonize(extensions)
      )



# setup(
    # packages = find_packages(),
    # ext_modules = cythonize(extensions)
# )
