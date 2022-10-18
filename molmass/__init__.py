# molmass/__init__.py

from .molmass import __doc__, __all__, __version__
from .molmass import *
from .elements import *

from . import elements

__all__ = __all__ + elements.__all__
del elements
