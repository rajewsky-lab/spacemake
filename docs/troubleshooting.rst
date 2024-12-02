Troubleshooting
===============
Below is a list of known issues you may encounter when running spacemake.
As spacemake depends on several external libraries and tools, it is not
always possible to resolve some of these issues.

Bumped on another problem which is not documented here? Feel free to `open
an issue on Github. <https://github.com/rajewsky-lab/spacemake/issues>`_


GLIBCXX_xxx not found
^^^^^^^^^^^^^^^^^^^^^
In certain environments you might run into the following error:

.. code-block:: console

    ImportError: /lib/x86_64-linux-gnu/libstdc++.so.6: version `GLIBCXX_3.4.29' not found (required by /.../envs/spacemake/lib/python3.10/site-packages/matplotlib/_path.cpython-310-x86_64-linux-gnu.so)

Certain dependencies (e.g. scipy) might affect this. To overcome it, try the following:

.. code-block:: console

   export LD_LIBRARY_PATH=/conda_installation_folder/envs/spacemake/lib:$LD_LIBRARY_PATH

For more details on this and further troubleshooting visit `stackoverflow <https://stackoverflow.com/questions/72540359/glibcxx-3-4-30-not-found-for-librosa-in-conda-virtual-environment-after-tryin>`_.


Issues with memory usage
^^^^^^^^^^^^^^^^^^^^^^^^
Certain steps of the spacemake workflow might currently result in excessive memory
usage. These occur in samples where a large number of barcodes exist in the data,
such as for ``stero-seq`` or ``open-ST``. The excessive memory usage is due to using
the ``Drop-seq`` tools, where a specific memory size is allocated for ``java``.

We are working on removing the depndency to the ``Drop-seq`` tools altogether,
which will also speed up several steps. If you run into memory errors, however, you
can solve them by modifying the ``main.smk`` file inside your spacemake installation,
which should be somewhere in

.. code-block:: console

    /path_to_conda/envs/spacemake/lib/python3.10/site-packages/spacemake/snakemake

inside your conda installation folder. Simply modify the following lines

.. code-block:: console

    {dropseq_tools}/BamTagHistogram -m 32g
    {dropseq_tools}/DigitalExpression -m 16g

by increasing the value of ``-m`` accordingly.


Issues with STAR
^^^^^^^^^^^^^^^^
To reduce memory usage when running several samples at the same time, 
spacemake uses STAR's shared memory capability. This currently has
the following limitations:

1. It is not possible for one user to run two distinct spacemake instances with the same genome index. Multiple spacemake instances (each processing several samples) can run at the same time if different species indexes are used.
2. Similarly, it is not possible that two users can run spacemake with the same genome index loaded at the same time.

In addition to the above, you might run into STAR-related errors if the spacemake
instance was killed before finishing. This occurs when the genome index is still loaded
into memory, and STAR will either throw an error and exit, or just stall. In that case,
try to run:

.. code-block:: console

    STAR --genomeLoad Remove --genomeDir <genome_dir>

In case the shared memory cannot be released (Linux), try the following:

.. code-block:: console

    ipcs -m | grep `whoami` | awk '{ print $2 }' | xargs -n1 ipcrm -m
