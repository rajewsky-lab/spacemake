Troubleshooting
===============
Below is a list of known issues you may encounter when running spacemake.
As spacemake depends on several external libraries and tools, it is not
always possible to resolve some of these issues.

Bumped on another problem which is not documented here? Feel free to `open
an issue on Github. <https://github.com/rajewsky-lab/spacemake/issues>`_

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


GLIBCXX_xxx not found
^^^^^^^^^^^^^^^^^^
In certain environments you might run into the following error:

.. code-block:: console

    ImportError: /lib/x86_64-linux-gnu/libstdc++.so.6: version `GLIBCXX_3.4.29' not found (required by /.../envs/spacemake/lib/python3.10/site-packages/matplotlib/_path.cpython-310-x86_64-linux-gnu.so)

Certain dependencies (e.g. scipy) might affect this. To overcome it, try the following:

.. code-block:: console

   export LD_LIBRARY_PATH=/conda_installation_folder/envs/spacemake/lib:$LD_LIBRARY_PATH

For more details on this and further troubleshooting visit `stackoverflow <https://stackoverflow.com/questions/72540359/glibcxx-3-4-30-not-found-for-librosa-in-conda-virtual-environment-after-tryin>`_.