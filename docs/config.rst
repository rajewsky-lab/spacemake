Configuration
=============

Once installed and initialized, spacemake needs to be configured.

.. include:: shared/shared_sample_variables.rst

Configure species
-----------------

.. _configure-species:

To add species, the following command can be used::

   spacemake config add_species \
       --name NAME \         # name of the species to be added
       --reference REF \     # name of the reference sequence 
                             # ('genome', 'rRNA', 'spike_in', ...)
                             # if omitted defaults to 'genome'
       --sequence SEQUENCE \ # path to the reference sequence file 
                             # (.fa) to be added
       --genome SEQUENCE \   # DEPRECATED! Please use --sequence instead.
       --annotation ANNOTATION \
                             # path to the annotation (.gtf) file for the species
                             # to be added

The ``spacemake config update-species`` takes the same arguments as above, while ``spacemake config delete-species`` takes only ``--name``.

As of version ``0.7`` you can add multiple reference sequences per species. For that, 
simply execute ``add-species`` multiple times, varying ``--reference ...`` but keeping ``--name`` constant.


To list the currently available ``species``, type::
   
   spacemake config list-species

Configure adapter-flavors and pre-processing
--------------------------------------------

.. _configure-adapter-flavor:

Spacemake allows to pre-process raw reads based on adapter-flavors. An adapter-flavor describes how adapters and
polyA stretches should be trimmed from the cDNA read (usually read2). The complete set of operations that can be performed are:
- trim polyA stretches
- trim adapters
- clip low-quality bases
- clip fixed number of bases from either end of read2

Access to these operations is provided through the ``adapter-flavors`` section of ``config.yaml`` only. Here is an example of an adapter-flavor:

.. code-block:: yaml

   adapter_flavors:
      example:
         - nextseq_quality:
               cutoff: 25
         - polyA:
         - adapter:
               name: SMART
               seq: AAGCAGTGGTATCAACGCAGAGTGAATGGG
               where: left
               min_overlap: 10
               max_errors: 0.1


Below follows a list of each operation and the supported parameters and default values. 

quality
^^^^^^^

Trim low-quality bases from 3' end and/or 5' end of read. Functionality is provided by `cutadapt <https://cutadapt.readthedocs.io/en/stable/guide.html#quality-trimming>`_.
Two parameters are supported: ``left`` and ``right``, which define the quality threshold below which bases will be trimmed from the 5' and 3' end of read2, respectively. 
Default is ``left: 0`` and ``right: 25``.

nextseq_quality
^^^^^^^^^^^^^^^

Trim low-quality bases from 3' end of read. Functionality is provided by `cutadapt <https://cutadapt.readthedocs.io/en/stable/guide.html#quality-trimming>`_.
The sole parameter is ``cutoff``, which defines the quality threshold below which bases will be trimmed. Analogous to quality with ``right=cutoff``, except that terminal `G` nucleotides
are always treated as below cutoff quality. Default is ``cutoff: 25``.

.. note::

    Before version ``0.9.1`` there was no quality trimming of bases at all, which led to issues on some runs. Between versions ``0.9.1`` and ``0.9.5``, the default was set 
    to ``nextseq_quality`` with ``cutoff: 32``, which is a common default for quality trimming, but relatively strict. In version ``0.9.5`` the default was changed to ``cutoff: 25``, 
    which is in our experience a good compromise, because low quality bases may still be soft-clipped in the mapping stage. However, if you experience a drop in UMI counts between pre ``0.9.1`` and 
    current versions, you can try lowering the quality cutoff further (or even set it to 0) and rerun your samples, to restore pre ``0.9.1`` behavior.

clip
^^^^

Clip bases from either end of read2. Two parameters are supported: ``left`` and ``right``, which define how many bases should be clipped from the 5' and 3' end of read2, 
respectively. Default is ``left: 0`` and ``right: 0``.

polyA
^^^^^

Trim polyA stretches from 3' end of read. Functionality is provided by `cutadapt <https://cutadapt.readthedocs.io/en/stable/guide.html#quality-trimming>`_.
The only supported parameter is ``revcomp``, which if set to ``True`` will trim polyT stretches instead of polyA. Default is ``revcomp: False``.

adapter
^^^^^^^

Trim adapters from either end of read. Functionality is provided by `cutadapt <https://cutadapt.readthedocs.io/en/stable/guide.html#search-parameters>`_.
Paraneters are:

  - ``name``: name of the adapter. Only for logging purposes.
  - ``seq``: sequence of the adapter to be trimmed.
  - ``min_overlap``: minimum overlap between read and adapter for a successful trimming. Default is ``3``.
  - ``max_errors``: maximum error rate allowed for a successful trimming. Default is ``0.1``.
  - ``where``: where to search for the adapter. Possible values are ``'left'``, and ``'right'``. Default is ``'right'`` (3 prime end of cDNA).

.. note::
   Internally, spacemake uses the cutadapt python module to perform all trimming operations. If ``where == 'left'`` we use ``cutadapt.adapters.NonInternalFrontAdapter``,
   for ``where == 'right'`` we use ``cutadapt.adapters.BackAdapter``.

   For more information about the parameters and their meaning, please refer to the `cutadapt source code <https://github.com/marcelm/cutadapt/blob/main/src/cutadapt/adapters.py>`_.

Each adapter-flavor in the ``config.yaml`` is a list of operations to be performed in the given order. If needed, you can chain multiple operations of the same type (for 
example to remove multiple adapters).


CRAM/BAM tags with pre-processing info
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Note that `spacemake` keeps a record of pre-processing steps
for each read in the CRAM/BAM file tags, so it is always possible to track which operations were performed:
- A3: comma-separated list of adapters detected and trimmed from the 3'end. May also contain "polyA" if polyA trimming was performed and/or "Q" if quality trimming was performed.
- A5: comma-separated list of adapters detected and trimmed from the 5'end.
- T3: comma-separated list of number of bases trimmed from the 3'end (synced with A3).
- T5: comma-separated list of number of bases trimmed from the 5'end (synced with A5).

Here is an example:

.. code-block:: text

   read_name 163 chr1 1000 60 50M = 1050 100 ACGT... * NM:i:0 A3:Z:Q,polyA A5:Z:SMART T3:Z:5,10 T5:22

The tags indicate that from the 3' end of the read, first 5 bases were trimmed due to quality (Q), then 10 bases of a polyA stretch.
From the 5' end, 22 bases of SMART adapter were trimmed.

Configure barcode-flavors
--------------------------

.. _configure-barcode-flavor:

This sample-variable describes how the cell-barcode and the UMI should be extracted from Read1 and Read2.
The ``default`` value for barcode\_flavor will be dropseq: ``cell = r1[0:12]`` (cell-barcode comes from first 12nt of Read1) and
``UMI = r1[12:20]`` (UMI comes from the 13-20 nt of Read1). 

**If a sample has no barcode\_flavor provided, the default barcode\_flavor will be used**

Barcode correction
^^^^^^^^^^^^^^^^^^

As of version ``0.9.3``, spacemake performs spatial barcode correction with edit distance 1, which boosts counts by ~5-15% for many samples.
For performance reasons, this employs some heuristics:
- all ``N`` bases are replaced with ``A``, in the reference (flowcell) catalog, as well as in the samples.
- a capture-area catalog of reference barcodes is built for each samples, based on exact match counts alone.
- exact matches to the capture-area catalog are searched first and preferred. Unmatched barcodes go on to a second stage of potential error correction.
- spacemake looks all edit distance 1 variants of an unmatched sample barcode in the capture-area catalog in a defined order. 
  The first match is reported and no further matches are considered. 
  The order is as follows: (1) substitutions, (2) insertions, (3) deletions. This means that if a barcode has no exact matches, but **multiple** edit 1 matches, the 
  correction will be deterministic, but is not guaranteed to be correct. In practice, however, the fraction of barcodes with multiple edit 1 matches is extremely low 
  and dwarfed by other sources of experimental and technical noise.

.. note::
   Barcode correction requires to configure ``--puck-barcode-files`` for your sample. Otherwise it will not be treated as a spatial sample and no capture-area catalog 
   can be built.
   
Provided barcode-flavors
^^^^^^^^^^^^^^^^^^^^^^^^^

.. note::

   Future versions of spacemake will merge barcode-flavors into adapter-flavors (which arguably become pre-processing flavors at that point)
   by defining ``barcode`` as a pre-processing step with ``cell`` and ``UMI`` as parameters.
   In the current implementation, barcode-flavors are kept separate for backwards compatibility. The new implementation will give additional
   flexibity, for example to remove additional adapters/primers, or clip the read further, after barcode extraction. Currently, if ``barcode``
   is not in the list of pre-processing steps, it is taken to be implied as the last step and its parameters are loaded from the barcode-flavor.

Spacemake provides the following barcode-flavors out of the box:

.. code-block:: yaml

    default:
        cell: "r1[0:12]"
        UMI: "r1[12:20]"
    openst:
        cell: "r1[2:27]"
        UMI: "r2[0:9]"
    sc_10x_v2:
        cell: "r1[0:16]"
        UMI: "r1[16:26]"
    seq_scope:
        UMI: "r2[0:9]"
        cell: "r1[0:20]"
    slide_seq_14bc:
        cell: "r1[0:14]"
        UMI: "r1[14:23]"
    slide_seq_15bc:
        cell: "r1[0:14]"
        UMI: "r1[15:23]"
    visium:
        cell: "r1[0:16]"
        UMI: "r1[16:28]"

To list the currently available ``barcode-flavor``-s, type::
   
   spacemake config list_barcode-flavors

.. warning::

    The command line interface for adding, updating, and deleting barcode-flavors will be deprecated in future versions of spacemake.
    Please consider editing the ``config.yaml`` file directly to manage barcode-flavors.

Add a new barcode\_flavor
^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block::

   spacemake config add_barcode-flavor \
      --name NAME \
         # name of the barcode flavor

      --umi UMI \
         # structure of UMI, using python's list syntax.
         # Example: to set UMI to 13-20 NT of Read1, use --umi r1[12:20].
         # It is also possible to use the first 8nt of Read2 as UMI: --umi r2[0:8].

      --cell-barcode CELL-BARCODE
         # structure of CELL BARCODE, using python's list syntax.
         # Example: to set the cell-barcode to 1-12 nt of Read1, use --cell-barcode r1[0:12].
         # It is also possible to reverse the CELL BARCODE, for instance with r1[0:12][::-1]. 


Update/delete a barcode-flavor
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The ``spacemake config update-barcode-flavor`` takes the same arguments as above, while ``spacemake config delete-barcode-flavor`` takes only ``--name``.

Configure run-modes
--------------------

.. _configure-run-mode:

Specifying a "run mode" is an essential flexibity that spacemake offers.
Through setting a ``run-mode``, a sample can be processed and analysed downstream in various fashions.

Each ``run-mode`` can have the following variables:

``n_beads``
   number of cell-barcode expected

``umi_cutoff``
   a list of integers. downstream the analysis will be run using these UMI cutoffs,
   that is cell-barcodes with less UMIs will be discarded

``clean_dge``
   whether to clean cell-barcodes from overhang primers, before creating the DGE.

``detect_tissue`` (spatial only)
   if ``True``, apart from UMI cutoff spacemake will try to detect the tissue *in-silico*.

``polyA_adapter_trimming``
   if ``True`` 3' polyA stretches and apaters will be trimmed from Read2.

``count_intronic_reads``
   if ``True`` intronic reads will be counted when creating the DGE.

``count_mm_reads``
   if ``True`` multi-mappers will be counted. Only those multi-mapping reads will be
   counted this way, which map to exactly one CDS or UTR segment of a gene.

``mesh_data`` (spatial only)
   if ``True`` a mesh will be created when running this ``run-mode``.

``mesh_type`` (spatial only)
   spacemake currently offers two types of meshes: (1) ``circle``, where circles with a given
   ``mesh_spot_diameter_um`` will be placed in a hexagonal grid, ``mesh_spot_distance_um``
   distance apart; (2) a hexagonal grid, where equal hexagons with ``mesh_spot_diameter_um``
   sides will be placed in a full mesh grid, such that the whole area is covered.

``mesh_spot_diameter_um`` (spatial only)
   the diameter of the mesh spatial-unit, in microns.

``mesh_spot_distance_um`` (spatial only, only for circle mesh)
   distance between the meshed circles, in microns.

``spatial_barcode_min_matches`` (spatial only)
   ratio spatial barcode matches, expressed as 0-1 interval, used as a minimum threshold to
   filter out pucks from DGE creation and subsequent steps of the pipeline. If set to 0, 
   no pucks are excluded.

``parent_run-mode``
   Each ``run-mode`` can have a parent, to which it will fall back.
   If a one of the ``run-mode`` variables is missing, the variable of the parent will be used.
   If parent is not provided, the ``default`` ``run-mode`` will be the parent. 

Provided run-modes
^^^^^^^^^^^^^^^^^^^^^

.. code-block:: yaml

  default:
    clean_dge: false
    count_intronic_reads: true
    count_mm_reads: false
    detect_tissue: false
    mesh_data: false
    mesh_spot_diameter_um: 55
    mesh_spot_distance_um: 100
    mesh_type: circle
    n_beads: 100000
    polyA_adapter_trimming: true
    spatial_barcode_min_matches: 0
    umi_cutoff:
    - 100
    - 300
    - 500
  openst:
    clean_dge: false
    count_intronic_reads: true
    count_mm_reads: true
    detect_tissue: false
    mesh_data: true
    mesh_spot_diameter_um: 7
    mesh_spot_distance_um: 7
    mesh_type: hexagon
    n_beads: 100000
    polyA_adapter_trimming: true
    spatial_barcode_min_matches: 0.1
    umi_cutoff:
    - 100
    - 250
    - 500
  scRNA_seq:
    count_intronic_reads: true
    count_mm_reads: false
    detect_tissue: false
    n_beads: 10000
    umi_cutoff:
    - 500
  seq_scope:
    clean_dge: false
    count_intronic_reads: false
    count_mm_reads: false
    detect_tissue: false
    mesh_data: true
    mesh_spot_diameter_um: 10
    mesh_spot_distance_um: 15
    mesh_type: hexagon
    n_beads: 1000
    umi_cutoff:
    - 100
    - 300
  slide_seq:
    clean_dge: false
    detect_tissue: false
    n_beads: 100000
    umi_cutoff:
    - 50
  visium:
    clean_dge: false
    count_intronic_reads: false
    count_mm_reads: true
    detect_tissue: true
    n_beads: 10000
    umi_cutoff:
    - 1000

.. note::
   If a sample has no ``run-mode`` provided, the ``default`` will be used

.. note:: 
   If a ``run-mode`` variable is not provided, the variable of the default ``run-mode`` will be used

To list the currently available ``run-mode``-s, type::
   
   spacemake config list_run-modes

.. warning::

    The command line interface for adding, updating, and deleting run\_modes will be deprecated in future versions of spacemake.
    Please consider editing the ``config.yaml`` file directly to manage run-modes.

Add a new run\_mode
^^^^^^^^^^^^^^^^^^^

See the :ref:`variable descriptions <configure-run-mode>` above.

.. code-block::

   spacemake config add_run-mode \
      --name NAME \ 
      --parent_run_mode PARENT_RUN_MODE \
      --umi_cutoff UMI_CUTOFF [UMI_CUTOFF ...] \
      --n_beads N_BEADS \
      --clean_dge {True,true,False,false} \
      --detect_tissue {True,true,False,false} \
      --polyA_adapter_trimming {True,true,False,false} \
      --count_intronic_reads {True,true,False,false} \
      --count_mm_reads {True,true,False,false} \
      --mesh_data {True,true,False,false} \
      --mesh_type {circle,hexagon} \
      --mesh_spot_diameter_um MESH_SPOT_DIAMETER_UM \
      --mesh_spot_distance_um MESH_SPOT_DISTANCE_UM

Update/delete a run-mode
^^^^^^^^^^^^^^^^^^^^^^^^^

The ``spacemake config update-run-mode`` takes the same arguments as above, while ``spacemake config delete-run-mode`` takes only ``--name``.


Configure pucks
---------------

.. _configure-puck:

Each spatial sample is associated with a ``puck``. The ``puck`` variable defines the 
dimensionality of the underlying spatial structure, which spacemake uses
during the automated analysis and plotting, as well as the binning (meshing) of
the data when selected in the ``run-mode``.

Each puck has the following variables:

- ``width_um``: the width of the puck, in microns
- ``spot_diameter_um``: the diameter of bead on this puck, in microns.
- ``barcodes`` (optional): the path to the barcode file, containing the cell\_barcode
  and (x,y) position for each. This is handy when several pucks have the same barcodes,
  such as for 10x Visium.
- ``coordinate_system`` (optional): the path to the coordinate system file, containing puck
  IDs and the (x,y,z) position for each, in global coordinates. This coordinate system is analogous
  to the global coordinate system for image stitching. When specified, this 'stitching' is
  automatically performed on ``puck``-s with spatial information.


Provided pucks
^^^^^^^^^^^^^^

.. code-block:: yaml

  default:
    coordinate_system: ''
    spot_diameter_um: 10
    width_um: 3000
  openst:
    coordinate_system: puck_data/openst_coordinate_system.csv
    spot_diameter_um: 0.6
    width_um: 1200
  seq_scope:
    spot_diameter_um: 1
    width_um: 1000
  slide_seq:
    spot_diameter_um: 10
    width_um: 3000
  visium:
    barcodes: puck_data/visium_barcode_positions.csv
    spot_diameter_um: 55
    width_um: 6500

The ``visium`` puck comes with a ``barcodes`` variable, which points to
``puck_data/visium_barcode_positions.csv``. Similarly, the ``openst`` puck comes with
a ``coordinate_system`` variable, pointing to ``puck_data/openst_coordinate_system.csv``.

Upon initiation, these files will automatically placed there by spacemake

To list the currently available ``puck``-s, type::
   
   spacemake config list_pucks


.. warning::

    The command line interface for adding, updating, and deleting pucks will be deprecated in future versions of spacemake.
    Please consider editing the ``config.yaml`` file directly to manage pucks.


Add a new puck
^^^^^^^^^^^^^^

.. code-block::

   spacemake config add_puck \
      --name NAME \        # name of the puck
      --width_um WIDTH_UM \
      --spot_diameter_um SPOT_DIAMETER_UM \
      --barcodes BARCODES \ # path to the barcode file, optional 
      --coordinate_system COORDINATE_SYSTEM # path to the coordinate system file, optional 


Custom snakemake rules
----------------------

As of version ``0.7`` it is now possible to add custom snakemake rules to your spacemake workflow.
Simply add the following line to the ``config.yaml`` in your spacemake root folder:

.. code-block:: yaml

   custom_rules: /path/to/my_own_custom_snakefile.smk

Within your custom code, you can import spacemake modules and have access to internal variables.
If you need to make spacemake aware of new top-level targets that have to be made, 
you can register a callback

.. code-block:: python

   register_module_output_hook(get_my_custom_targets, "my_own_custom_snakefile.smk")

The function `get_my_custom_targets()` will be called once all other, internal spacemake code has been executed
and is expected to return a list of files that will be appended to the ``input:`` dependencies of the top-level 
rule. Providing rules to make these files is up to your custom rules.

The second parameter is more for logging purposes and allows to track which module or part of the 
code injected which dependencies. By default, it is good practive to use the filename.

