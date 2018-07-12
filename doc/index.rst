mqc
###

Contents:

.. toctree::
   :maxdepth: 2

   intro
   callme
   developers
   design_decisions
   plateau_detection
   api


Installation
============

First, create or activate a virtual environment (python>=3.6), where cython and pysam are
installed. Conda or virtualenv is both ok, with conda you would create an
environment as follows: ::

    conda create -n mqc python=3.6
    source activate mqc
    conda install cython pysam

Then retrieve the code and install with develop install mode: ::

    cd $where_you_want
    git clone https://eilslabs-phabricator.dkfz.de/diffusion/177/mqc.git
    cd mqc
    git checkout develop
    pip install -e .

To run the tests, you may want to separate unit tests and acceptance tests
(which could take a couple of minutes): ::

    echo -e 'Running unit tests\n*****************************************'
    pytest -m "not acceptance_test" mqc/tests/

    echo -e 'Running acceptance tests\n*****************************************'
    pytest -m "acceptance_test" mqc/tests/


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
