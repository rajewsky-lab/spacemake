Installation
============

Step 1: create conda environment
--------------------------------

The most straightforward way to install spacemake, is to first create conda environment with the above packages.

To do this, we highly recommend using `mamba <https://github.com/mamba-org/mamba>`_, a much faster conda package manager than conda itself.

After mamba is installed, download the `environment.yaml <https://raw.githubusercontent.com/rajewsky-lab/spacemake/master/environment.yaml>`_. This environment.yaml, contains all the dependencies required by spacemake.

Once downloaded, to install all spacemake dependencies type::

    mamba env create -f environment.yaml

This will create a conda environment called ``spacemake``.

Too activate the newly created environment type::

   conda activate spacemake

Step 2: download Dropseq-tools
------------------------------

To work with spacemake, currently it is needed to download `Dropseq-tools from here <https://github.com/broadinstitute/Drop-seq>`_.
This packages is a collection of processing tools originally written for `Drop-seq <https://www.cell.com/cell/fulltext/S0092-8674(15)00549-8>`_. Spacemake uses several functions from this package during pre-processing and processing and without it it is impossible to run spacemake.

Simply download one of the releases (we recommend using `2.5.1 <https://github.com/broadinstitute/Drop-seq/releases/download/v2.5.1/Drop-seq_tools-2.5.1.zip>`_) and place it somewhere in your filesystem.


Step 3: install spacemake
-------------------------

**After creating the conda environment and downloading Dropseq-tools** (as described above) spacemake can be installed via ``pip``::

   pip install spacemake

This will install spacemake, you should be good to go :)

.. warning::
    Make sure to first create the conda environment as described above.

    Although it is also possible to install the required packages independently, and then
    to install spacemake, this option has not been tested, and one can quickly run into
    dependency issues and errors.

To install the developmental version of spacemake (``dev`` branch from github) type the following command::

   pip install git+https://github.com/rajewsky-lab/spacemake.git@dev

