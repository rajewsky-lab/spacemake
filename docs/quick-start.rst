Quick-start
===========

The examples here are minimal code pieces how to start spacemake. They assume that you have initialised spacemake
with ``spacemake init`` as described :ref:`here <_installation>`.

As spacemake comes with no ``default`` value for ``species``, before you can start anything, you need to add one
using ``spacemake config add_species`` as described :ref:`here <_configure-species>`.

Visium quick-start
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

As you can see, above we add a new visium project with ``puck, run_mode, barcode_flavor`` all set to ``visium``.

This is possible as spacemake comes with pre-defined variables, all suited for visium. The visium ``run_mode`` will process the 
sample in the same way as `spaceranger <https://support.10xgenomics.com/spatial-gene-expression/software/pipelines/latest/what-is-space-ranger>`_ would: intronic reads will not be counted, multi-mappers (where the multi-mapping read maps only to one CDS or UTR region) will be counted,
3' polyA stretches will not be trimmed from Read2.

**To add several visium projects at once, follow :ref:`the tutorial here <_add-several-projects>`**
