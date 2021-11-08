Installation
============

.. _installation:

Create a conda environment
--------------------------

The most straightforward way of installing spacemake, is first creating a conda environment with the above packages.

To do this, we highly recommend using `mamba <https://github.com/mamba-org/mamba>`_, a much faster conda package manager than conda itself.

After mamba is installed, download the `environment.yaml <https://raw.githubusercontent.com/rajewsky-lab/spacemake/dev/environment.yaml>`_. This environment.yaml, contains all the dependencies required by spacemake.

Once downloaded, to install all spacemake dependencies type::

    mamba env create -f environment.yaml

This will create a conda environment called ``spacemake``.

Too activate the newly created environment type::

   conda activate spacemake

Installing spacemake
--------------------

After creating the environment (as described above) spacemake can be installed via ``pip``::

   pip install spacemake

This will install spacemake, you should be good to go :)

To install the developmental version of spacemake (``dev`` branch from github) type the following command::

   pip install git+https://github.com/rajewsky-lab/spacemake.git@dev

