Managing projects and samples
=============================

In spacemake each sample, and it's settings, are stored in the ``project_df.csv`` under the root
directory of the spacemake project.

Each sample will have exactly one row in this ``project_df.csv`` file. In the back-end, spacemake uses a ``pandas.DataFrame`` to load, and save this ``.csv`` file on disk. This data-frame
will be indexed by key ``(project_id, sample_id)``

The spacemake class responsible for this back-end logic is the :ref:`ProjectDF<ProjectDF>` class.

Add a single sample
-------------------

In spacemake each sample can have the folloing variables:

``project_id``
   ``project_id`` of a sample

``sample_id``
   ``sample_id`` of a sample

``R1``
   .fastq.gz file path to Read1 reads

``R2``
   .fastq.gz file path to R2 reads

``longreads`` (optional)
   fastq(.gz)|fq(.gz)|bam file path to pacbio long reads for library debugging

``longread-signature`` (optional)
   identify the expected longread signature (see longread.yaml)

``barcode_flavor`` (optional)
   ``barcode_flavor`` of the sample. If not provided, ``default`` will be used (Drop-seq).

``species``
   ``species`` of the sample

``puck`` (optional)
   name of the ``puck`` for this sample. if puck contains a ``barcodes`` variable, with a path
   to a coordinate file, those coordinates will be used when processing this sample.
   If not provided, a ``default`` puck will be used with ``width_um=3000``,
   ``spot_diameter_um=10``.

``puck_id`` (optional)
   ``puck_id`` of a sample

``puck_barcode_file`` (optional)
    the path to the file contining (x,y) positions of the barcodes. If the ``puck`` for this
    sample has a ``barcodes`` variable, it will be ignored, and ``puck_barcode_file`` will
    be used.

``investigator`` (optional)
   name investigator(s) responsible for this sample

``experiment`` (optional)
   description of the experiment

``sequencing_date``
   sequencing date of the sample

``run_mode`` (optional)
   A list of ``run_mode`` names for this sample. The sample will be processed as defined in 
   the ``run_mode``-s provided. If not provided, the ``default`` ``run_mode`` will be used.


To add a single sample, we can use the following command::

   spacemake projects add_sample \
      --project_id PROJECT_ID \                 # required
      --sample_id SAMPLE_ID \                   # required
      --R1 R1 \                                 # required, if no longreads
      --R2 R2 \                                 # required, if no longreads
      --longreads LONGREADS \                   # required, if no R1 & R2
      --longread-signature LONGREAD_SIGNATURE \ # optional
      --barcode_flavor BARCODE_FLAVOR \         # optional
      --species SPECIES \                       # required
      --puck PUCK \                             # optional
      --puck_id PUCK_ID \                       # optional
      --puck_barcode_file PUCK_BARCODE_FILE \   # optional
      --investigator INVESTIGATOR \             # optional
      --experiment EXPERIMENT \                 # optional
      --sequencing_date SEQUENCING_DATE \       # optional
      --run_mode RUN_MODE [RUN_MODE ...] \      # optional


.. warning::

   A sample is spatial only if: either a ``puck_barcode_file`` is provided, or the sample's
   ``puck`` has a ``barcodes`` variable pointing to a barcode position file.
   If this is not the case, spacemake won't be able to find the spatial barcodes for
   this sample, and the sampe will be processed as a single-cell sample.

   In case both the ``puck_barcode_file`` is provided and the sample's ``puck`` has the
   ``barcodes`` variable set, ``puck_barcode_file`` will be used for the spatial coordinates.

Add several samples at once
---------------------------

.. _add-several-samples:

It is possible to add several samples in just one command. First, the sample variables have
to be defined in a ``samples.yaml`` file, then we can run the following command::

   spacemake projects add_samples_from_yaml --samples_yaml samples.yaml

The ``samples.yaml`` should have the following structure:

.. code-block:: yaml

   additional_projects:
      - project_id: visium
        sample_id: visium_1
        R1: <path_to_visium_1_R1.fastq.gz>
        R2: <path_to_visium_1_R2.fastq.gz>
        species: mouse
        puck: visium
        barcode_flavor: visium
        run_mode: [visium]
      - project_id: visium
        sample_id: visium_2
        R1: <path_to_visium_2_R1.fastq.gz>
        R2: <path_to_visium_2_R2.fastq.gz>
        species: human
        puck: visium
        barcode_flavor: visium
        run_mode: [visium]
      - project_id: slideseq
        sample_id: slideseq_1
        R1: <path_to_slideseq_1_R1.fastq.gz>
        R2: <path_to_slideseq_1_R2.fastq.gz>
        species: mouse
        puck: slideseq
        barcode_flavor: slideseq_14bc
        run_mode: [default, slideseq]
        puck_barcode_file: <path_to_slideseq_puck_barcode_file>

Under ``additional_projects`` we define a list where each element will be a key:value pair, to be inserted in the ``project_df.csv``

.. note::
   When using the above command, if a sample is already present in the ``project_df.csv`` rather than added again, it will be updated with the new values, if they changed.
   
   If someone runs ``spacemake projects add_samples_from_yaml --samples yaml samples.yaml`` and
   then modifies something in the ``samples.yaml``, and runs the command again, the ``project_df.csv``
   will contain the updated version of the settings.

Add samples from illumina sample-sheet
--------------------------------------

Coming soon...
