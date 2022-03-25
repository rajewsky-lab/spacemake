Initialization
==============

Initializing using required arguments
-------------------------------------

.. include:: shared/spacemake_init.rst

Optional arguments
------------------

The `spacemake init` command takes the following optional arguments:

``root_dir``
    The ``root_dir`` for the spacemake instance. Defaults to ``.``, the directory in which `spacemake init` is ran.

``temp_dir``
    Path to the temporary directory, defaults to ``/tmp``.

``download_species``
    If set, spacemake will download the genome (.fa) and annotation (.gtf) files for mouse and human (from gencode, as specified `here <https://github.com/rajewsky-lab/spacemake/blob/master/spacemake/data/config/species_data_url.yaml>`_.

Hence, the complete `spacemake init` command looks like this::
    
    spacemake init \
      --root_dir ROOT_DIR \             # optional
      --temp_dir TEMP_DIR \             # optional
      --download_species \              # optional
      --dropseq_tools DROPSEQ_TOOLS     # required
