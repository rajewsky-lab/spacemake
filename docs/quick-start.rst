Quick start guide
=================

The examples here are minimal code pieces how to start spacemake. They assume that you have initialised spacemake
with ``spacemake init`` as described :ref:`here <installation>`.

As spacemake comes with no ``default`` value for ``species``, before you can start anything, you need to add one
using ``spacemake config add_species`` as described :ref:`here <configure-species>`.

Visium quick start
------------------

To add a visium project, simply type::

   spacemake projects add_sample \
      --project_id <project_id> \
      --sample_id <sample_id> \
      --R1 <path_to_R1.fastq.gz> \
      --R2 <path_to_R2.fastq.gz> \
      --species <species> \
      --puck visium \
      --run_mode visium \
      --barcode_flavor visium

Above we add a new visium project with ``puck, run_mode, barcode_flavor`` all set to ``visium``.

This is possible as spacemake comes with pre-defined variables, all suited for visium. The visium ``run_mode`` will process the 
sample in the same way as `spaceranger <https://support.10xgenomics.com/spatial-gene-expression/software/pipelines/latest/what-is-space-ranger>`_ would: intronic reads will not be counted, multi-mappers (where the multi-mapping read maps only to one CDS or UTR region) will be counted,
3' polyA stretches will not be trimmed from Read2.

To see the values of these predefined variables checkout the :ref:`configuration:Configuration`.

**To add several visium projects at once, follow** :ref:`the tutorial here <add-several-projects>`

Slide-seq quick start
---------------------

To add a Slide-seq project, simply type::

   spacemake projects add_sample \
      --project_id <project_id> \
      --sample_id <sample_id> \
      --R1 <path_to_R1.fastq.gz> \
      --R2 <path_to_R2.fastq.gz> \
      --species <species> \
      --puck slideseq \
      --run_mode slideseq \
      --barcode_flavor slideseq_14bc \
      --puck_barcode_file <path_to_puck_barcode_file>

Above we add a new Slide-seq project with the ``puck, run_mode`` will be set to ``slideseq``
which are pre-defined settings for Slide-seq samples.

.. note::
   For spatial samples other than visium - such as Slide-seq - we need to provide a
   ``puck_barcode_file`` (since each puck has different barcodes, unlike for visium samples).
   This file should be a comma or tab separated, containing column names as first row. Acceptable column names are:

   - ``cell_bc``, ``barcodes``  or ``barcode`` for cell-barcode
   - ``xcoord`` or ``x_pos`` for x-positions
   - ``ycoord`` or ``y_pos`` for y-positions

In this example ``barcode_flavor`` will be set to ``slideseq_14bc``,
a pre-defined ``barcode_flavor`` in spacemake, where the ``cell_barcode`` comes from the first 14nt of Read1, and the ``UMI`` comes from nt 13-22 (remaining 9 nt). 
The other pre-defined ``barcode_flavor`` for Slide-seq is ``slideseq_15bc``: here ``cell_barcode`` again comes from the first 14nt of Read1, but the ``UMI`` comes from nt 14-22 (remaining 8) of Read1.

To see the values of these predefined variables checkout the :ref:`configuration:Configuration`.

**To add several slideseq projects at once, follow** :ref:`the tutorial here <add-several-projects>`

Drop-seq quick start
-------------------

spacemake was written as a spatial-transcriptomics pipeline, however it will also work for
single-cell experiments, where there is no spatial information available. 

To add a drop-seq project, simply type::

   spacemake projects add_sample \
      --project_id <project_id> \
      --sample_id <sample_id> \
      --R1 <path_to_R1.fastq.gz> \
      --R2 <path_to_R2.fastq.gz> \
      --species <species> \
      --run_mode scRNAseq

As seen above, we define fewer variables as before: only ``species`` and ``run_mode`` are needed.

Behind the scenes, the ``barcode_flavor`` will be set to ``default``, which by default contains 
the Drop-seq barcoding strategy (``cell_barcode`` is 1-12nt of Read1, ``UMI`` is 13-20 nt of Read1).

For 10X samples either the ``sc_10x_v2`` (10X Chromium Single Cell 3' V2)
or ``visium`` (10X Chromium Single Cell 3' V3, same as visium) would had to be used as
``barcode_flavor``. Both are pre-defined in spacemake.

By setting ``run_mode`` to ``scRNAseq`` we used the pre-defined ``run_mode`` settings tailored for single-cell experiments: expected number of beads will be 10k, introns will be counted, UMI cutoff will be at 500, multi-mappers will not be counted and polyA and adapter sequences will be trimmed from Read2. 

Of course, running single-cell samples with other ``run_mode`` settings is also possible. 

To see the values of these predefined variables checkout the :ref:`configuration:Configuration`.

**To add several single-cell projects at once, follow** :ref:`the tutorial here <add-several-projects>`

Running spacemake
-----------------

After a sample is added, like in one of the examples above, spacemake can be run with::

   spacemake run --cores <n_cores> --keep-going

The ``--keep-going`` flag is optional, however it will ensure that spacemake runs all
the jobs it can, even if one job fails (this logic is directly taken from snakemake).

Listing projects
----------------

To list projects which were added so far, and their information, type::

   spacemake projects list

