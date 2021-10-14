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
   the variable which specifies the structure of R1 and R2, namely how the cell\_barcode and UMI should be extracted. If no value provided for a sample, the default will be used.

run\_mode
   each sample can have several run\_modes, all of which are user definable. If no run\_mode is specified, a sample will be processed using the default settings.

puck (for spatial samples only)
   if a sample is spatial, it has to have a puck variable. If no puck is specified, a default puck will be used.  


Each of these variables can be added, updated and deleted with::

   spacemake config add_<variable>
   spacemake config update_<variable>
   spacemaek config delete_<variable>

where ``<variable`` can be ``species, barcode_flavor, run_mode or puck``

Adding species
--------------

To add species, the following command can be used::

   spacemake config add_species \
   --name NAME           # name of the species to be added \
   --genome GENOME       # path to the genome (.fa) file for the species to be added \
   --annotation ANNOTATION \
                         # path to the annotation (.gtf) file for the species to be added \
   --rRNA_genome RRNA_GENOME \
                         # (optional) path to the ribosomal-RNA genome (.fa) file for the species to be added
