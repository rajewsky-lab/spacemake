import pytest
import sys
import os
from spacemake.cmdline import *

from fixtures import (
    initialized_root,
    with_species,
    with_tile_test_data,
    configured_root,
    sm,
    spacemake_dir,
)


def test_from_yaml(with_species):
    import spacemake.config as smc

    os.chdir(with_species.as_posix())
    config = smc.ConfigFile.from_yaml()

    print(config.get_variable("default_mapper_settings", name="STAR"))
