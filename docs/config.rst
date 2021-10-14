Configuration
=============

Once installed, spacemake configured before running.

Initialise
----------

To initialise a spacemake project, ``cd`` into the directory, in which you want to start spacemake, this will be your project\_root.

Then simply type::
   
   spacemake init --dropseq_tools <path_to_dropseq_tools_dir> --picard_tools <path_to_picardtools_picard.jar>

Optionally, you can also provide the ``--download_species`` flag, which will download Gencode genomes and annotations for mouse and human, and place them under project\_root/species\_data/<species>,
where <species> is either mouse or human.

Important variables to set
--------------------------

In spacemake, each has to have the following variables:

species
   a collection of genome, annotation and rRNA\_genome. There is no default species, and each sample can have exactly one species.

barcode\_flavor
   the variable which specifies the structure of Read1 and Read2, namely how the cell\_barcode and UMI should be extracted. If no value provided for a sample, the default will be used.

run\_mode
   each sample can have several run\_modes, all of which are user definable. If no run\_mode is specified, a sample will be processed using the default settings.

puck (for spatial samples only)
   if a sample is spatial, it has to have a puck variable. If no puck is specified, a default puck will be used.  


Each of these variables can be added, updated and deleted with::

   spacemake config add_<variable>
   spacemake config update_<variable>
   spacemaek config delete_<variable>

where ``<variable`` can be ``species, barcode_flavor, run_mode or puck``

Configure species
-----------------

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

Configure barcode\_flavor
-------------------------

The ``default`` value for barcode\_flavor will be dropseq: ``cell_barcode = r1[0:12]`` (cell-barcode comes from first 12nt of Read1) and ``UMI = r1[12:20]`` (UMI comes from the 13-20 nt of Read1).

To add a new barcode\_flavor::

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


The ``spacemake config update_barcode_flavor_species`` takes the same arguments as above, while ``spacemake config delete_barcode_flavor`` takes only ``--name``.

Provided barcode\_flavors
^^^^^^^^^^^^^^^^^^^^^^^^^

Spacemake provides the following barcode flavors out of the box::

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
