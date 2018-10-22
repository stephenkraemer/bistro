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

    *With pip*
    A develop install should pull all dependencies as specified in the setup.py. Currently,
    figure_report must be preinstalled by hand, as well as cython and pysam.

    As I have not yet tested the package across a range of pysam and cython versions, the safest thing
    is to use the versions I used for development.

    pip install cython=0.26
    pip install pysam=0.11.2.2

    git clone $figure_report_repo  # still need to publish this - for now ask me for it
    cd figure_report
    pip install -e .

    git clone https://eilslabs-phabricator.dkfz.de/diffusion/177/mqc.git
    cd mqc
    git checkout develop
    pip install -e .[dev]


    *With conda*
    Basically, you'll want to create an environment with conda and install all dependencies
    with conda. Then install the figure_report and mqc packages as develop install using pip -e.
    As I have not yet tested the package across a range of pysam and cython versions, the safest thing
    is to use the versions I used for development.

    conda create -n bistro_dev python=3.6
    source activate bistro_dev
    conda install pysam=0.11.2.2 cython=0.26
    # install other dependencies through conda
    # until this is possible by utilizing the meta.yaml / recipe
    # (https://github.com/conda/conda/issues/6788), the fastest way
    # is to copy the meta.yaml requirements in an environment file

    cd $where_you_want

    git clone $figure_report_repo  # still need to publish this - for now ask me for it
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
