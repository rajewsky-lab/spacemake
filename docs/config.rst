Configuration
=============

Once installed, spacemake configured before running.

Initialise
----------

To initialise a spacemake project, ``cd`` into the directory in which you want to start spacemake, this will be your ``project_root``.

Then simply type::
   
   spacemake init \
      --dropseq_tools <path_to_dropseq_tools_dir> \
      --picard_tools <path_to_picardtools_picard.jar>

Optionally, you can also provide the ``--download_species`` flag, which will download Gencode genomes and
annotations for ``mouse`` and ``human``, and place them under ``project\_root/species\_data/<species>``,
where <species> is either mouse or human.

Sample-variables
----------------

In spacemake, each sample has to have the following sample-variables:

species
   a collection of genome, annotation and rRNA\_genome. There is no default species, and each sample can have exactly one species.

barcode\_flavor
   the variable which specifies the structure of Read1 and Read2, namely how the cell\_barcode and UMI should be extracted. If no value provided for a sample, the default will be used.

run\_mode
   each sample can have several ``run_mode``-s, all of which are user definable. If no ``run_mode``-s are specified, a sample will be processed using ``default`` ``run_mode`` settings.

puck (for spatial samples only)
   if a sample is spatial, it has to have a puck variable. If no puck is specified, a default puck will be used.  


To add, update, delete or list a sample-variable, you can use the following commands::

   spacemake config add_<sample-variable>
   spacemake config update_<sample-variable>
   spacemake config delete_<sample-variable>
   spacemake config list_<sample-variable>

where ``<sample-variable>`` can be ``species, barcode_flavor, run_mode or puck``

Configure species
-----------------

.. _configure-species:

To add species, the following command can be used::

   spacemake config add_species \
       --name NAME \         # name of the species to be added
       --genome GENOME \     # path to the genome (.fa) file for the species to
                             # be added
       --annotation ANNOTATION \
                             # path to the annotation (.gtf) file for the species
                             # to be added
       --rRNA_genome RRNA_GENOME
                             # (optional) path to the ribosomal-RNA genome (.fa)
                             # file for the species to be added

The ``spacemake config update_species`` takes the same arguments as above, while ``spacemake config delete_species`` takes only ``--name``.

To list the currently available ``species``, type::
   
   spacemake config list_species

Configure barcode\_flavors
--------------------------

.. _configure-barcode_flavor:

This sample-variable describes how the cell-barcode and the UMI should be extracted from Read1 and Read2.
The ``default`` value for barcode\_flavor will be dropseq: ``cell_barcode = r1[0:12]`` (cell-barcode comes from first 12nt of Read1) and
``UMI = r1[12:20]`` (UMI comes from the 13-20 nt of Read1). 

**If a sample has no barcode\_flavor provided, the default run\_mode will be used**

Provided barcode\_flavors
^^^^^^^^^^^^^^^^^^^^^^^^^

Spacemake provides the following barcode\_flavors out of the box:

.. code-block:: yaml

    default:
        cell: "r1[0:12]"
        UMI: "r1[12:20]"
    slideseq_14bc:
        cell: "r1[0:14]"
        UMI: "r1[14:23]"
    slideseq_15bc:
        cell: "r1[0:14]"
        UMI: "r1[15:23]"
    visium:
        cell: "r1[0:16]"
        UMI: "r1[16:28]"
    sc_10x_v2:
        cell: "r1[0:16]"
        UMI: "r1[16:26]"

To list the currently available ``barcode_flavor``-s, type::
   
   spacemake config list_barcode_flavors

Add a new barcode\_flavor
^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block::

   spacemake config add_barcode_flavor \
      --name NAME \
         # name of the barcode flavor

      --umi UMI \
         # structure of UMI, using python's list syntax.
         # Example: to set UMI to 13-20 NT of Read1, use --umi r1[12:20].
         # It is also possible to use the first 8nt of Read2 as UMI: --umi r2[0:8].

      --cell_barcode CELL_BARCODE
         # structure of CELL BARCODE, using python's list syntax.
         # Example: to set the cell_barcode to 1-12 nt of Read1, use --cell_barcode r1[0:12].
         # It is also possible to reverse the CELL BARCODE, for instance with r1[0:12][::-1]. 


Update/delete a barcode\_flavor
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The ``spacemake config update_barcode_flavor`` takes the same arguments as above, while ``spacemake config delete_barcode_flavor`` takes only ``--name``.

Configure run\_modes
--------------------

.. _configure-run_mode:

Specifying a "run mode" is an essential flexibity that spacemake offers.
Through setting a ``run_mode``, a sample can be processed and analysed downstream in various fashions.

Each ``run_mode`` can have the following variables:

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
   if ``True`` a mesh will be created when running this ``run_mode``.

``mesh_type`` (spatial only)
   spacemake currently offers two types of meshes: (1) ``circle``, where circles with a given
   ``mesh_spot_diameter_um`` will be placed in a hexagonal grid, ``mesh_spot_distance_um``
   distance apart; (2) a hexagonal grid, where equal hexagons with ``mesh_spot_diameter_um``
   sides will be placed in a full mesh grid, such that the whole area is covered.

``mesh_spot_diameter_um`` (spatial only)
   the diameter of the mesh spatial-unit, in microns.

``mesh_spot_distance_um`` (spatial only, only for circle mesh)
   distance between the meshed circles, in microns.

``parent_run_mode``
   Each ``run_mode`` can have a parent, to which it will fall back.
   If a one of the ``run_mode`` variables is missing, the variable of the parent will be used.
   If parent is not provided, the ``default`` ``run_mode`` will be the parent. 

Provided run\_mode(s)
^^^^^^^^^^^^^^^^^^^^^

.. code-block:: yaml

    default:
        n_beads: 100000
        umi_cutoff: [100, 300, 500]
        clean_dge: False
        detect_tissue: False
        polyA_adapter_trimming: True
        count_intronic_reads: True
        count_mm_reads: False
        mesh_data: False
        mesh_type: 'circle'
        mesh_spot_diameter_um: 55
        mesh_spot_distance_um: 100
    visium:
        n_beads: 5000
        umi_cutoff: [1000]
        clean_dge: False
        detect_tissue: True
        polyA_adapter_trimming: False
        count_intronic_reads: False
        count_mm_reads: True
    slideseq:
        n_beads: 100000
        umi_cutoff: [50]
        clean_dge: False
        detect_tissue: False
    scRNAseq:
        n_beads: 10000
        umi_cutoff: [500]
        detect_tissue: False
        polyA_adapter_trimming: True
        count_intronic_reads: True
        count_mm_reads: False

.. note::
   If a sample has no ``run_mode`` provided, the ``default`` will be used

.. note:: 
   If a ``run_mode`` variable is not provided, the variable of the default ``run_mode`` will be used

To list the currently available ``run_mode``-s, type::
   
   spacemake config list_run_modes

Add a new run\_mode
^^^^^^^^^^^^^^^^^^^

.. code-block::

   spacemake config add_run_mode \
      --name NAME \        # name of the run_mode to be added
      --parent_run_mode PARENT_RUN_MODE \
                           # Name of the parent run_mode.
                           # All run_modes will fall back to 'default'
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

Update/delete a run\_mode
^^^^^^^^^^^^^^^^^^^^^^^^^

The ``spacemake config update_run_mode`` takes the same arguments as above, while ``spacemake config delete_run_mode`` takes only ``--name``.


Configure pucks
---------------

.. _configure-puck:

Each spatial sample, needs to have a ``puck``. The ``puck`` sample-variable will define the 
dimensionality of the underlying spatial structure, which then spacemake will use
during the autmated analysis and plotting. 

Each puck has the following variables:

- ``width_um``: the width of the puck, in microns
- ``spot_diameter_um``: the diameter of bead on this puck, in microns.
- ``barcodes`` (optional): the path to the barcode file, containing the cell\_barcode
  and (x,y) position for each. This is handy, when several pucks have the same barcodes,
  such as for 10x visium.


Provided pucks
^^^^^^^^^^^^^^

.. code-block:: yaml

    default:
        width_um: 3000
        spot_diameter_um: 10
    visium:
        barcodes: 'puck_data/visium_barcode_positions.csv'
        width_um: 6500
        spot_diameter_um: 55
    seqscope:
        width_um: 1000
        spot_diameter_um: 1

as you can see, the ``visium`` puck comes with a ``barcodes`` variable, which points to
``puck_data/visium_barcode_positions.csv``. Upon initiation, this file will automatically placed 
there by spacemake

To list the currently available ``puck``-s, type::
   
   spacemake config list_pucks


Add a new puck
^^^^^^^^^^^^^^

.. code-block::

   spacemake config add_puck \
      --name NAME \        # name of the puck
      --width_um WIDTH_UM \
      --spot_diameter_um SPOT_DIAMETER_UM \
      --barcodes BARCODES # path to the barcode file, optional 

