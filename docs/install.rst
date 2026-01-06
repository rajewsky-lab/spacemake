Installation
============


Using apptainer (recommended)
=============================
The easiest way to run spacemake is to use `apptainer <https://apptainer.org/>`_ (formerly known as Singularity). apptainer allows to run containerized applications without requiring root access.
This is especially useful when working on HPC clusters where users do not have root access.

If you have not installed apptainer yet, follow the installation instructions on the `apptainer website <https://apptainer.org/docs/admin/latest/installation.html>`_.
For a non-provileged installation, you can simply:

.. code-block:: console

    curl -s https://raw.githubusercontent.com/apptainer/apptainer/main/tools/install-unprivileged.sh | \
      bash -s - install-dir

where `install-dir` is the directory you want to install it into. Once installed, you can run it with `install-dir/bin/apptainer`. For the rest of this guide, we will assume that `apptainer` is in your PATH.

To run spacemake using apptainer, first download the latest spacemake apptainer image from the Rajewsky lab: `smmain.sif <https://bimsbstatic.mdc-berlin.de/rajewsky/containers/spacemake/smmain.sif>`_:

.. code-block:: console

    wget https://bimsbstatic.mdc-berlin.de/rajewsky/containers/spacemake/smmain.sif

You can then run spacemake using apptainer like this:

.. code-block:: console

    apptainer run -B /tmp --writable-tmpfs smmain.sif --help


For your convenience, you can create an alias for spacemake like this (we recommend to also add this line to your ``.bashrc`` or ``.zshrc`` file to make it more permanent):

.. code-block:: console

    alias spacemake="apptainer run -B /tmp --writable-tmpfs /path/to/smmain.sif"

This will allow you to run spacemake commands directly using simply ``spacemake``.

.. note:: 
   
    The apptainer version of spacemake comes pre-packaged with `Dropseq-tools`, so you do not need to download it separately.
    When starting a new project, simply omit the ``--dropseq-tools`` argument, and spacemake will automatically use the Dropseq-tools version included in the container.


.. note::

    Apptainer may not by default expose the file systems you want to work with inside the container.
    You can use the ``-B`` option to bind mount directories from your host system into the container.
    In the above example, we bind mount ``/tmp`` as an example (and to speed things up). You can bind mount any other directory you want 
    to work with by adding additional ``-B /path/on/host`` options to the apptainer command.

Installing via conda and pip
============================


Step 1: create conda environment
--------------------------------

The most straightforward way to install spacemake, is to first create conda environment with the above packages.
We highly recommend using `mamba <https://github.com/mamba-org/mamba>`_, a much faster conda package manager than conda itself.
After mamba is installed, download the `environment.yaml <https://raw.githubusercontent.com/rajewsky-lab/spacemake/master/environment.yaml>`_.
This file contains all dependencies required by spacemake.

Once downloaded, to install all spacemake dependencies type::

    mamba env create -f environment.yaml

This will create a conda environment called ``spacemake``. To activate the newly created environment type::

   conda activate spacemake

Step 2: download Dropseq-tools
------------------------------

Spacemake currently requires to download `Dropseq-tools <https://github.com/broadinstitute/Drop-seq>`_.
This package is a collection of processing tools originally written for `Drop-seq <https://www.cell.com/cell/fulltext/S0092-8674(15)00549-8>`_. Spacemake uses several functions from this package during pre-processing and processing, and without it it is impossible to run spacemake.

Simply download one of the releases (we recommend using `2.5.1 <https://github.com/broadinstitute/Drop-seq/releases/download/v2.5.1/Drop-seq_tools-2.5.1.zip>`_) and place it somewhere in your filesystem.


Step 3: install spacemake
-------------------------

**After creating the conda environment and downloading Dropseq-tools** (as described above),
spacemake can be installed via ``pip``::

   pip install spacemake

This will install spacemake, you should be good to go :)

.. warning::
    Make sure to first create the conda environment as described above.

    Although it is also possible to install the required packages independently, and then
    to install spacemake, this option has not been tested, and one can quickly run into
    dependency issues and errors.

To make sure spacemake has been properly installed, run::

   spacemake --version

This should output the latest spacemake version available in ``pip``.