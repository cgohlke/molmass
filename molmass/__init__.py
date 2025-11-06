# molmass/__init__.py

from . import elements
from .elements import *
from .molmass import *
from .molmass import __all__, __doc__, __version__

__all__ = __all__ + elements.__all__
del elements

# constants are repeated for documentation

__version__ = __version__
"""Molmass version string."""

ELEMENTS = ELEMENTS
"""Collection of chemical elements with lookup by number, symbol, and name.

An ordered collection of Element instances with lookup by
atomic number (int), chemical symbol (str), or element name (str).
"""

ELEMENTARY_CHARGE = ELEMENTARY_CHARGE
"""Elementary charge in coulombs."""

ELECTRON = ELECTRON
"""Electron particle."""

PROTON = PROTON
"""Proton particle."""

NEUTRON = NEUTRON
"""Neutron particle."""

POSITRON = POSITRON
"""Positron or antielectron particle."""
