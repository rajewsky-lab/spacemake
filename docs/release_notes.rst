Release notes
=============

Version 0.9
-----------
0.9.0
~~~~~
* Uses ``cram`` instead of ``bam`` to significantly reduce disk space.
* The ``migrate`` tool enable updating of v0.8x processed data into 0.9.
* Fixed bugs.

Version 0.8
-----------
0.8.1
~~~~~
* Fixed all automated analyses errors.
* Fixed a bug that occurred with toy data during the ``filter_mm_reads`` step.

0.8.0
~~~~~
* Faster cmdline interface
* improvement in speed for ``fastq_to_BAM``
* Added ``mrfifo`` implementation for the preprocessing steps
* Lower memory footprint for BamTagHistogram and for ``n_intersect_sequences``
* Added several unit tests with pytest
* Added kebab-case functionality
* Fixed bugs

Version 0.7
-----------
0.7.9
~~~~~
* Fixed issues with STAR shared memory indexes

0.7.7
~~~~~
This release provides better RAM usage and faster processing through the following features:
* Lower RAM consumption during the fastq to BAM step.
* Lower RAM consumption and increased ~5x faster processing during the ``h5ad`` creation step.
* STAR now loads the genome index only once in memory and shares it across workers.

0.7.6
~~~~~
* Fixed a critical bug when processing merging samples
* Updated documentation to describe the Open-ST workflow

0.7.5
~~~~~
* Fixed several bugs.
* Added support for open-ST ``run_mode``, and ``barcode_flavor``.

0.7.2
~~~~~
* Updated the ``environment.yaml`` to pin the pandas version.
* Fixed bugs and errors.

0.7.0
~~~~~
* Added support for multiple references per species (genome, rRNA, mirgenedb etc).
* The multiple species references are controlled through the ``map_strategy`` variable (e.g. rRNA->genome).
* Implemented an optimized algorithm for binning (meshing) coordinate systems.
* Added support for easy addition of custom rules into the workflow.
* Added unit tests.
* Updated the ``environment.yaml`` to pin the OpenJDK version.
* Fixed several bugs and errors.

Version 0.5
-----------

0.5.5
~~~~~
* Fixed issue with ``--version`` introduced in v0.5.4
* Fixed command line error where updating a sample would remove ``puck_ids``

0.5.4
~~~~~
* Each sample can have multiple barcode files
* Added a novosparc tutorial
* Fixed bug on read-only ``project_df.csv``` during spacemake run
* Fixed bug in adding new puck command

0.5.3
~~~~~
* Further bug fixes
* Completed docs with longread tutorial

0.5.2
~~~~~
* Fixed errors during initialization and species download
* Updated docs

0.5.1
~~~~~
* Added H&E integration
* Updated documentation
* Bug fixes

Version 0.4
-----------

0.4.3
~~~~~
* Bug fixes

0.4.2
~~~~~
* Added novosparc command line tools
* Bug fixes

0.4.1
~~~~~
Initial release of spacemake containing

* base spacemake for processing sequencing based spatial transcriptomics data
* default values for methods such as: Visium, Slide-seq, Seq-scope and general scRNA-seq
* modules: technical replicate merging, downsampling, long-read analysis
