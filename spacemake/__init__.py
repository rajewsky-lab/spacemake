# __version__ = 1.0
# import matplotlib._path
# from . import preprocess as pp
# from . import spatial as sp

# from .smk import Spacemake
import importlib

# implement lazy loading/importing to prevent slow commandline interface
# the 'mod' returned by import_module() is only really loading the module
# upon first access to its content. It is also cached through the normal
# mechanisms of module imports.

__all__ = ['pl', 'pp', 'sp']

_lazy_modules = {
  'pl' : '.plotting',
  'pp' : '.preprocess',
  'sp' : '.spatial',
}

def __getattr__(name):
    if name in _lazy_modules:
        mod = importlib.import_module(_lazy_modules[name], __name__)
        globals()[name] = mod
        return mod
    else:
        raise AttributeError(f"module {__name__} has no attribute {name}")
