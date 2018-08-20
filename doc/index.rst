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

Install release
---------------

    custom_conda_channel=/icgc/dkfzlsdf/analysis/B080/kraemers/projects/bistro/channel
    conda create -n bistro_release python=3.6
    source activate bistro_release
    conda install -c file:///$custom_conda_channel mqc=0.2

Develop install
---------------

    conda create -n bistro_dev python=3.6
    source activate bistro_dev
    conda install pysam cython

    cd $where_you_want

    git clone $figure_report_repo_TODO
    cd figure_report
    pip install -e .

    git clone https://eilslabs-phabricator.dkfz.de/diffusion/177/mqc.git
    cd mqc
    git checkout develop
    pip install -e .[dev]

Run tests
---------
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
