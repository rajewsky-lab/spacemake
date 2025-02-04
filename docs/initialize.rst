Initialization
==============

Initializing using required arguments
-------------------------------------

.. include:: shared/spacemake_init.rst

Optional arguments
------------------

The `spacemake init` command takes the following optional arguments:

``root-dir``
    The ``root-dir`` for the spacemake instance. Defaults to ``.``, the directory in which `spacemake init` is ran.

``temp-dir``
    Path to the temporary directory, defaults to ``/tmp``.

``download-species``
    If set, spacemake will download the genome (.fa) and annotation (.gtf) files for mouse and 
    human from gencode, as specified `here <https://github.com/rajewsky-lab/spacemake/blob/master/spacemake/data/config/species_data_url.yaml>`_.

Hence, the complete `spacemake init` command looks like this::
    
    spacemake init \
      --root-dir ROOT-DIR \             # optional
      --temp-dir TEMP-DIR \             # optional
      --download-species \              # optional
      --dropseq-tools DROPSEQ-TOOLS     # required
