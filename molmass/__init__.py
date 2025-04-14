# molmass/__init__.py

from . import elements
from .elements import *
from .molmass import *
from .molmass import __all__, __doc__, __version__

__all__ = __all__ + elements.__all__
del elements


def _set_module() -> None:
    """Set __module__ attribute for all public objects."""
    globs = globals()
    module = globs['__name__']
    for item in __all__:
        obj = globs[item]
        if hasattr(obj, '__module__'):
            obj.__module__ = module


_set_module()
