# __version__ = 1.0
# import matplotlib._path
# from . import preprocess as pp
# from . import spatial as sp

# from .smk import Spacemake
import importlib

__all__ = ['pl']

def __getattr__(name):
    if name == "pl":
        mod = importlib.import_module(".plotting", __name__)
        globals()[name] = mod
        return mod
    raise AttributeError(f"module {__name__} has no attribute {name}")