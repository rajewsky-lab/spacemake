.. include:: ../links.rst

Processing a custom single-cell sample
======================================

In this tutorial we will process a custom single cell sample. 

As an example we will be using 1 million reads from `this Visium dataset <https://www.10xgenomics.com/resources/datasets/mouse-brain-section-coronal-1-standard-1-0-0>`_.

.. note::
    
    Firstly, the example data used here is a 10X `Visium`_ dataset, hence it is spatial.
    However, for the sake of this tutorial, we will be treating it as a single-cell sample.
    
    Secondly, for many methods (such as `Visium`_, `10X Chromium`_ `Slide-seq`_ or `Seq-scope`_)
    spacemake provides pre-defined variables. If you are using
    one of these methods follow our :ref:`Quick start guide <quick start guide>` instead.

Step 1: install and initialize spacemake
-----------------------------------------

To install spacemake follow the :ref:`installation guide here <installation>`.

To initialize spacemake follow the :ref:`initialization guide here <initialization>`.

Step 2: download test data
--------------------------

For the sake of this tutorial we will work with a test dataset: 1 million Read1 and 1 million Read2 reads from a `Visium`_ adult mouse brain.

To download the test data:

.. code-block::

    wget -nv http://bimsbstatic.mdc-berlin.de/rajewsky/spacemake-test-data/visium/test_fastq/visium_public_lane_joined_1m_R1.fastq.gz
    wget -nv http://bimsbstatic.mdc-berlin.de/rajewsky/spacemake-test-data/visium/test_fastq/visium_public_lane_joined_1m_R2.fastq.gz

.. note:: 

    If there is already data available, to be processed and analyzed, this step can be omitted.

Step 3: add a new species
-------------------------

.. note::

    If you initialized spacemake with the ``--download-species`` flag, you can
    omit this step, as spacemake will automatically download and configure
    mm10 mouse genome.fa and annotation.gtf files for you.

The sample we are working with here is a mouse brain sample, so we have to add a new species:

.. code-block:: console

   spacemake config add_species --name mouse \
   --annotation /path/to/mouse/annotation.gtf \
   --genome /path/to/mouse/genome.fa


Step 4: add a new barcode\_flavor
---------------------------------

The ``barcode_flavor`` will decide which nucletodies of Read1/Read2 extract the UMIs and cell-barcodes from.

In this perticular test sample, the first 16 nucleotides of Read1 are the cell-barcode, and the following 12 nucleotides are the UMIs.

Consequently, we create a new ``barcode_flavor`` like this:

.. code-block:: console

    spacemake config add_barcode_flavor --name test_barcode_flavor \
    --cell_barcode r1[0:16] \
    --umi r1[16:28]

.. note:: 

    There are several ``barcode_flavors`` provided by spacemake out of the box,
    such as ``visium`` for 10X `Visium`_ or ``sc_10x_v2`` for `10X Chromium`_ v2 
    kits. The ``default`` flavor is identical to a `Drop-seq`_ library, with 12
    nucleotide cell-barcode and 8 nucleotide UMI. 

    :ref:`More info about provided flavors here <provided barcode\\_flavors>`.

    If you want to use one of these, there is no need to add your own flavor.

Step 5: add a new run\_mode
---------------------------

A ``run_mode`` in spacemake defines how a sample should processed downstream. 
In this tutorial, we will trim the PolyA stretches from the 3' end of Read2,
count both exonic and intronic reads, expect 5000 cells, and analyze the data,
turn off multi-mapper counting (so only unique reads are counted),
using 50, 100 and 300 UMI cutoffs. To set these parameters, we define a 
``test_run_mode`` like this:

.. code-block:: console

    spacemake config add_run_mode --name test_run_mode \
    --polyA_adapter_trimming True \
    --count_mm_reads False \
    --n_beads 5000 \
    --count_intronic_reads True \
    --umi_cutoff 50 100 300

.. note:: 

    As with ``barcode_flavors``, spacemake provides several ``run_modes`` out
    of the box. For more info :ref:`check out a more detailed guide here <configure run\\_modes>`.

Step 6: add the sample
----------------------

After configuring all the steps above, we are ready to add our (test) sample:

.. code-block:: console

    spacemake projects add_sample --project_id test_project \
    --sample_id test_sample \
    --R1 visium_public_lane_joined_1m_R1.fastq.gz \
    --R2 visium_public_lane_joined_1m_R1.fastq.gz \
    --species mouse \
    --barcode_flavor test_barcode_flavor \
    --run_mode test_run_mode

.. note::

    If there is already data available, here the Read1 and Read2 ``.fastq.gz`` files should be added,
    instead of the test files.

Step 7: runn spacemake
----------------------

Now we can process our samples with spacemake. Since we added only one sample, only one sample will be processed
and analyzed. To start spacemake, simply write:

.. code-block:: console
    
    spacemake run --cores 16

.. note::
    
    The number of cores used should be suited for the machine on which spacemake is ran.
    When processing more than one samle, we recommend using spacemake with at least 8 cores.
    In order to achieve maximum parallelism.

Step 8: results 
---------------

The results of the analysis for this sample will be under ``projects/test_project/processed_data/test_sample/illumina/complete_data/``

Under this directory, there are several files and directories which are important:

* ``final.polyA_adapter_trimmed.bam``: final, mapped, tagged ``.bam`` file. ``CB`` tag contains the cell barcode, and the ``MI`` contains the UMI-s. 

* ``qc_sheet_test_sample_no_spatial_data.html``: the QC-sheet for this sample, as a self-contained ``.html`` file.

* ``dge/``: a directory containing the Digital Expression Matrices (DGEs)

    * ``dge.all.polyA_adapter_trimmed.5000_beads.txt.gz``: a compressed, text based DGE

    * ``dge.all.polyA_adapter_trimmed.5000_beads.h5ad``: the same DGE but stored in ``.h5ad`` format (`used by the anndata python package <https://github.com/theislab/anndata/issues/180>`_). This matrix is stored as a Compressed Sparse Column matrix (using `scipy.sparse.csc_matrix <https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.csc_matrix.html>`_).

    * ``dge.all.polyA_adapter_trimmed.5000_beads.summary.txt``: the summary of the DGE, one line per cell.

    * ``dge.all.polyA_adapter_trimmed.5000_beads.obs.csv``: the observation table of the matrix. Similar to the previous file, more detailed.

* ``automated_analysis/test_run_mode/umi_cutoff_50/``: In this directory the results of the automated analysis can be found. As it can be seen under the ``automated_analysis`` directory there are two further levels, one for ``run_mode`` and one for ``umi_cutoff``. This is because one sample can have several ``run_modes`` and in the same way one ``run_mode`` can have several UMI cutoffs.

    * ``results.h5ad``: the result of the automated analysis, stored in an anndata object. Same as the DGE before, but containing processed data.

    * ``test_sample_no_spatial_data_illumina_automated_report.html``: automated analysis self-contained ``.html`` report.

.. note::

    If the ``test_project`` had more samples, than those would be automatically placed under ``projects/test_project``. Similarily, under one spacemake
    directory there can be several projects in parallel, and each will have their own directory structure under the ``projects/`` folder.

