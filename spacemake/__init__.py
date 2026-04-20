import importlib
from .contrib import __version__, __author__, __email__, __license__

# implement lazy loading/importing to prevent slow commandline interface.
# we only really load a sub-module upon first access to it.
# >>> import spacemake # <- fast
# >>> import spacemake.pl # <- slow (matplotlib pulls in a lot of weight)
# >>> import spacemake.pl as pl2 # <- fast (cached, returns the same module as above import)

__all__ = ["pl", "pp", "sp", "__version__", "__author__", "__email__", "__license__"]

_lazy_modules = {
    "pl": ".plotting",
    "pp": ".preprocess",
    "sp": ".spatial",
}


def __getattr__(name):
    g = globals()
    if name in g:
        # we already loaded this before. Don't even go through importlib at all
        return g[name]

    elif name in _lazy_modules:

        # first access. Need to load
        mod = importlib.import_module(_lazy_modules[name], __name__)
        globals()[name] = mod
        return mod

    else:
        # re-create normal error behavior
        raise AttributeError(f"module {__name__} has no attribute {name}")
