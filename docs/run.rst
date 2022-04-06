.. _Running spacemake general:

Running spacemake
=================

Main modules
------------

After spacemake in configured with the ``spacemake config`` command, and projects/samples
are added with the ``spacemake projects`` command, spacemake can be run with the 
``spacemake run`` command. It takes the following parameters::

    spacemake run \ 
        --cores CORES \     # number of cores to be used in total
        --dryrun, -n  \     # invokes a dry snakemake run, printing only commands
        --rerun-incomplete, --ri \
                            # forces snakemake to rerun incompletely generated files
        --keep-going  \     # if a job fails, keep executing independent jobs.
                            # we recommend to always set this when running spacemake
                            # overnight
        --printshellcmds, -p \
                            # print shell commands for each rule, if exist
        --touch, -t   \     # rather than running the rules, just touch each file
        --with_fastqc, -wfqc
                            # Run also fastqc as part of the spacemake run

Downsampling
------------

To run a downsampling (or saturation) analysis, one can use the following command::

    spacemake run downsample \
        --project_id_list [PROJECT_ID_LIST ...] \
        --sample_id_list [SAMPLE_ID_LIST ...]

In the ``project_id_list`` and ``sample_id_list`` arguments one can specify which a 
list of ``project_id``-s and ``sample_id``-s respectively, for which the downsampling
should be run. It is possible to set only one, or both of these arguments. If both are
set the downsampling will be run on samples for which the ``project_id`` and the ``sample_id`` are in both lists (intersection).

.. note::

    In addition to the list arguments specified above, the downsample command also
    takes the same arguments as the simple ``spacemake run`` command.
