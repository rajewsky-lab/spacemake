After a sample is added spacemake can be run with::

   spacemake run --cores <n_cores> --keep-going

The ``--keep-going`` flag is optional, however it will ensure that spacemake runs all
the jobs it can, even if one job fails (this logic is directly taken from snakemake).

For a complete explanation on the `spacemake run` command :ref:`check out the documentation here <Running spacemake general>`.
