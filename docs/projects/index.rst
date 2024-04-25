Manage projects and samples
===========================

In spacemake, each sample and its settings are stored in the ``project_df.csv`` under the root
directory of the spacemake project.

Each sample has exactly one row in this ``project_df.csv`` file. In the back-end, spacemake uses 
a ``pandas.DataFrame`` to load and save this ``.csv`` file on disk. This data frame
is indexed by key ``(project_id, sample_id)``

The spacemake class responsible for this back-end logic is the :ref:`ProjectDF<ProjectDF>` class.

Add a single sample
-------------------

Sample parameters
^^^^^^^^^^^^^^^^^

In spacemake each sample can have the folloing variables:

``project_id``
   ``project_id`` of a sample

``sample_id``
   ``sample_id`` of a sample

``R1``
   ``.fastq.gz`` file path(s) to Read1 read file(s). Can be either a single file, or a space separated list of consecutive files. If a list provided, the files will be merged together and the merged ``R1.fastq.gz`` will be processed downstream.

``R2``
    same as before, but for Read2 read file(s).    

``longreads`` (optional)
   fastq(.gz)|fq(.gz)|bam file path to pacbio long reads for library debugging

``longread-signature`` (optional)
   identify the expected longread signature (see longread.yaml)

``dge`` (optional)
    Since the ``0.1`` version of spacemake, it is possible to only provide the count matrix as input data for spacemake.
    Note: a raw count matrix is expected, if a non count matrix is provided, spacemake will raise an error. 

``barcode_flavor`` (optional)
   ``barcode_flavor`` of the sample. If not provided, ``default`` will be used (Drop-seq).

``species``
   ``species`` of the sample

``map_strategy``
    As of version ``0.7`` you can provide a mapping strategy as a string which gets converted into a 
    series of map rules. These rules translate into BAM names and their dependencies. 

    *Example:*

        `bowtie2:rRNA->STAR:genome:final`

    The input to the first mappings are always the pre-processed
    but unmapped reads.
    map-rules are composed of two or three paramters, separated by ':'.
    The first two parameters for the mapping are <mapper>:<reference>. The target BAM will have 
    the name <reference>.<mapper>.bam. Supported mappers are currently `bowtie2` and `STAR`
    The reference names must be defined for the species you assign to the sample.

    Optionally, a triplet can be used <mapper>:<reference>:<symlink> where the presence of <symlink> 
    indicates that the BAM file should additionally be made accessible under the name <symlink>.bam, a useful 
    shorthand or common hook expected by other stages of SPACEMAKE. A special symlink name is "final"
    which is required by downstream stages of SPACEMAKE. If no "final" is specified, the last map-rule
    automatically is selected and symlinked as "final".
    Note, due to the special importance of "final", it may get modified to contain other flags of the run-mode
    that are presently essential for SPACEMAKE to have in the file name ("final.bam" may in fact be 
    "final.polyA_adapter_trimmed.bam")

    The example above is going to create
        
        (1) rRNA.bowtie2.bam
        
            using bowtie2 and the index associated with the "rRNA" reference


        (2) genome.STAR.bam 
        
            using STAR on the *unmapped* reads from BAM (1)


        (3) final.bam
        
            a symlink pointing to the actual BAM created in (2).

    Note that one BAM must be designated 'final.bam', or the last BAM file created will be selected as final.
    (used as input to downstream processing rules for DGE creation, etc.)

    .. note:: 
      
        Parallel mappings can be implemented by using commata:

        `bowtie2:rRNA:rRNA,STAR:genome:final`

        This rule differs from the first example because it will align the unmapped reads from the uBAM
        in parallel to the rRNA reference and to the genome. In this way the same reads can match to both
        indices.

    .. note:: 
      
      Gene tagging will be applied automatically if annotation data were provided for the associated 
      reference index (by using 'spacemake config add_species --annotation=... ')

    .. warning::
      
       Please be sure to escape your `map_strategy` with double-quotes to prevent bash from
       interpreting the '>' as a redirect instruction.



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

``sequencing_date`` (optional)
   sequencing date of the sample

``run_mode`` (optional)
   A list of ``run_mode`` names for this sample. The sample will be processed as defined in 
   the ``run_mode``-s provided. If not provided, the ``default`` ``run_mode`` will be used.


To add a single sample, we can use the following command::

   spacemake projects add_sample \
      --project_id PROJECT_ID \                 # required
      --sample_id SAMPLE_ID \                   # required
      --R1 R1 [R1 R1 ...] \                     # required, if no longreads
      --R2 R2 [R2 R2 ...] \                     # required, if no longreads
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

Add a Visium/Seq-scope/Slide-seq sample
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Currently spacemake works out of the box with three spatial methods: `Visium <https://www.10xgenomics.com/products/spatial-gene-expression>`_, `Seq-scope <https://www.sciencedirect.com/science/article/abs/pii/S0092867421006279>`_ and `Slide-seq <https://pubmed.ncbi.nlm.nih.gov/33288904/>`_.

* To add a Visium sample, follow the :ref:`quick start guide here <step 1: add a visium sample>`.
* To add a Seq-scope sample, follow the :ref:`quick start guide here <step 1: add a seq-scope sample>`.
* To add a Slide-seq sample, follow the :ref:`quick start guide here <step 1: add a slide-seq sample>`.

Add a custom spatial sample
^^^^^^^^^^^^^^^^^^^^^^^^^^^

In order to process a custom spatial sample with spacemake follow the step by step guide below.

Step 1: specifying a puck
"""""""""""""""""""""""""

Each spatial sample will need a so-called puck to be configured first. By 'puck' we mean the physical properties of the underlying methods.
Visium for instance works with 6.5mm by 6.5mm sized capture areas, where each spot has 55 microns diameter. To configure a custom puck :ref:`follow the guide here <configure pucks>`.

.. warning::

    If a puck is not specified, spacemake will still run but will use the ``default`` puck as specified :ref:`here <provided pucks>`.

Step 2: formatting a custom puck_barcode_file
"""""""""""""""""""""""""""""""""""""""""""""

For all spatial samples we need to provide a ``puck_barcode_file``. This file needs to be a comma or tab separated, and it needs to have the following three (named) columns:

   - ``cell_bc``, ``barcodes``  or ``barcode`` for cell-barcode
   - ``xcoord`` or ``x_pos`` for x-positions
   - ``ycoord`` or ``y_pos`` for y-positions

Step 3: configure run\_mode(s), barcode\_flavor and species
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Before a custom sample is added the run\_mode(s), barcode\_flavor and species should be configured. The guides on how to do this can be found :ref:`here for run-modes <configure run\\_modes>`, :ref:`here for <configure barcode\\_flavors>` and :ref:`here for species <configure species>`.

The configured run\_mode(s) will specify how a sample is processed downstream, and the barcode\_flavor will specify the barcoding strategy used (ie how many nucleotides are used for UMI, which nucleotides are used for the spot barcodes).

.. warning::

    If no run\_mode(s) are provided spacemake will use the ``default`` run\_mode as specified :ref:`here <provided run\\_mode(s)>`.

    Similarily if there is no barcode\_flavor specified spacemake will use the ``default`` barcode\_flavor as specified :ref:`here <provided barcode\\_flavors>`.

Step 4: add your sample
"""""""""""""""""""""""

Once everything is configured you can add your custom spatial sample with the following command::

    spacemake projects add_sample \
        # your sample's project_id \
        --project_id PROJECT_ID \
        # your sample's sample_id \
        --sample_id SAMPLE_ID \
        # one or more R1.fastq.gz files
        --R1 R1 [R1 R1 ...] \
        # one or more R2.fastq.gz files
        --R2 R2 [R2 R2 ...] \
        # name of the barcode\_flavor, configured in Step 3 \
        --barcode_flavor BARCODE_FLAVOR \
        # name of the species, configured in Step 3 \
        --species SPECIES \
        # name of the puck, configured in Step 1 \
        --puck PUCK \
        # path to your custom barcode file, configured in Step 2 \
        --puck_barcode_file PUCK_BARCODE_FILE \
        # name of the run\_mode(s), configured in Step 3 \
        --run_mode RUN_MODE [RUN_MODE ...]

Add a single-cell sample
^^^^^^^^^^^^^^^^^^^^^^^^

To add a single-cell sample follow the :ref:`quick start guide here <step 1: add a single-cell rna-seq sample>`.

Add a pre-processed count-matrix
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Coming soon!

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
   When using the above command, if a sample is already present in the ``project_df.csv`` rather than adding it again, spacemake will update it.
   
   If someone runs ``spacemake projects add_samples_from_yaml --samples yaml samples.yaml`` and
   then modifies something in the ``samples.yaml``, and runs the command again, the ``project_df.csv``
   will contain the updated version of the settings.

Add samples from an Illumina sample-sheet
--------------------------------------

You can add samples directly from an Illumina sample-sheet, assuming the sample-sheet is configured appropriately and a basecalls folder is available. Spacemake will then automatically process the sample-sheet, create the appropriate directories, and begin demultiplexing the data. Once the data is demultiplexed, spacemake will continue with the processing as described above.

To use this functionality, type::

   spacemake projects add_sample_sheet \
       --sample_sheet <path_to_sample_sheet> \
       --basecalls_dir <path_to_basecalls_folder>

The sample-sheet columns have to obey certain conventions for spacemake to parse it properly:

* ``Sample_ID`` contains the ``sample_id`` in the project.
* ``Sample_Project`` contains the ``project_id`` in the project.
* ``Description`` must end with ``_species``, where species is the one configured for the samples in the project, e.g. ``HEK293_wt_human``.

Spacemake will also parse the fields ``Investigator``, ``Date``, and ``Experiment`` from the sample-sheet and add them to the project metadata.

Listing projects
----------------

To list projects, which are already configured and added, simply type::
    
    spacemake projects list

It will show the main variables for each project in the ``project_df.csv``. 

To view extra variables which are not shown, use the ``--variables`` option 
to specify which extra variables to show.

Merging samples
----------------

Spacemake can merge samples that have been resequenced to increase the number of quantified molecules in the data. To merge samples, first configure, add, and process the individual samples as they are. Make sure that the samples belong in the same project, e.g. have the same ``project_id``. Then merge them by typing::

   spacemake projects merge_samples \
       --merge_project_id <project_id> \
       --merged_sample_id <sample_merged> \
       --sample_id_list <sample_a> <sample_b>

The above command will merge the two samples by creating a new sample with the same variables. Spacemake performs the merging at the level of the ``bam`` files, thus properly processing the merged sample by collapsing PCR duplicates. Processing will automatically run until the creation of the ``qc_sheets`` and the automated analyses.
