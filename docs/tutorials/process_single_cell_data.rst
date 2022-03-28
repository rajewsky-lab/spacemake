Processing a custom single-cell sample
======================================

In this tutorial we will process a custom single cell sample. 

As an example we will be using 1 million reads from `this Visium dataset <https://www.10xgenomics.com/resources/datasets/mouse-brain-section-coronal-1-standard-1-0-0>`_.

Step 1: inistall and initialize spacemake
-----------------------------------------

To install spacemake follow the :ref:`installation guide here <installation>`.

To initialize spacemake follow the :ref:`initialization guide here <initialization>`.

Step 2: add a new species
-------------------------

.. note::

    If you initialized spacemake with the `--download-species` flag, you can
    omit this step, as spacemake will automatically download and configure
    mm10 mouse genome.fa and annotation.gtf files for you.

The sample we are working with here is a mouse brain sample, so we have to add a new species:

.. code-bloc:: console

    
