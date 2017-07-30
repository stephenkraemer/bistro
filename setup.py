from setuptools import setup, find_packages
from setuptools.extension import Extension
from Cython.Build import cythonize
from distutils import sysconfig


extensions = [
    Extension(
        "mqc.pileup.bsseq_pileup_read",
        ["mqc/pileup/bsseq_pileup_read.pyx"],
        include_dirs = ['/home/kraemers/programs/python_virtualenvs/python_3.6.0/include/python3.6m',
                        '/home/kraemers/programs/python_virtualenvs/python_3.6.0/lib/python3.6/site-packages/pysam/',
                        '/home/kraemers/programs/python_virtualenvs/python_3.6.0/lib/python3.6/site-packages/pysam/include/htslib'],
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
      packages=['mqc'],

      entry_points = {
          'console_scripts': ['mqc=mqc.cli:mqc'],
      },
      install_requires=[
          'pandas',
          'numpy',
          'seaborn',
          'pytoml',
          'pytest',
          'click',
          'joblib',
          'cython (>=0.25)',
          'python-magic (>=0.4.13)',
          'pytest-mock',
      ],
      ext_modules=cythonize(extensions)
      )



# setup(
    # packages = find_packages(),
    # ext_modules = cythonize(extensions)
# )
