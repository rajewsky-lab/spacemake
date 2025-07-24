.. _Migrate general:

Migrate
=======

Starting with v0.9, spacemake adopts the ``cram`` format to reduce disk requirements.

.. note::

    In addition to the ``.fa`` genome sequence, the ``cram`` format requires an index file, 
    ``.fai``. If the file does not exist, it will be created during the ``spacemake run``.
    In that case, ensure you have write permissions for the directory containing the sequence file.

To smoothly migrate existing project structures created with earlier versions (v0.8x)
and make them compatible with v0.9+, spacemake offers the ``migrate`` tool.

If you attempt to run a newer spacemake version inside a project directory created
with an older one, you'll get the following warning::

    ###############################################################
    #                                                             #
    #      WARNING: incompatible version detected                 #
    #                                                             #
    #   You are using a newer version of spacemake which is       #
    #   incompatible with your current project setup.             #
    #                                                             #
    #     Please migrate your projects using:                     #
    #      $ spacemake migrate                                    #
    #                                                             #
    #   For more info, see:                                       #
    #   https://spacemake.readthedocs.io/en/latest/migrate.html   #
    #                                                             #
    ###############################################################

To reshape the existing processed data and make it compatible with the latest
spacemake version, use ``spacemake migrate``::

    spacemake migrate \
        --project-id                # the project-id of the sample to migrate
        --sample-id                 # the sample-id of the sample to migrate
        --all                       # if provided, spacemake will migrate all samples in the
                                    # current folder
        --cores CORES               # number of cores to be used in total
        --dryrun, -n                # invokes a dry snakemake run, printing only commands
        --rerun-incomplete, --ri    # forces snakemake to rerun incompletely generated files
        --keep-going                # if a job fails, keep executing independent jobs.
                                    # we recommend to always set this when running spacemake
                                    # overnight
        --printshellcmds, -p        # print shell commands for each rule, if exist
        --touch, -t                 # rather than running the rules, just touch each file

.. note::

    The ``migrate`` command requires the user to provide either a pair of
    (``project-id``, ``sample-id``), or the ``--all`` flag. If the ``--all`` flag is
    provided, spacemake will migrate all samples within the current directory. 

The migration will convert all ``.bam`` files into ``.cram``, remove the ``.bam`` files,
and update the ``config.yaml`` to register spacemake's current version. The list of
.bam files removed during migration will be stored inside a file called ``bam_files_removed.txt``,
together with the total disk space saved.

.. warning::
    
    Migration will permanently remove all ``.bam`` files from the project directory.
    Make sure to adjust your workflows for working with ``.cram`` files and/or keep copies
    of ``.bam`` if it is necessary.