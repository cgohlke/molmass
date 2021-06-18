# elements.py

# Copyright (c) 2005-2021, Christoph Gohlke
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice,
#    this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its
#    contributors may be used to endorse or promote products derived from
#    this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.

"""Properties of the chemical elements.

Each chemical element is represented as an object instance. Physicochemical
and descriptive properties of the elements are stored as instance attributes.

:Author: `Christoph Gohlke <https://www.lfd.uci.edu/~gohlke/>`_

:License: BSD 3-Clause

:Version: 2021.6.18

Requirements
------------
* `CPython >= 3.7 <https://www.python.org>`_

Revisions
---------
2021.6.18
    Remove support for Python 3.6 (NEP 29).
    Add Particle types (#5).
2020.1.1
    Update atomic weights and isotopic compositions from NIST.
    Move element descriptions into separate module.
    Remove support for Python 2.7 and 3.5.
    Update copyright.
2018.8.15
    Move modules into molmass package.
2018.5.25
    Style and docstring fixes.
2016.2.25
    Fixed some ionization energies.

References
----------
1. https://www.nist.gov/pml/
   atomic-weights-and-isotopic-compositions-relative-atomic-masses
2. https://en.wikipedia.org/wiki/{element.name}

Examples
--------
>>> from molmass import ELEMENTS
>>> ele = ELEMENTS['C']
>>> ele.number
6
>>> ele.symbol
'C'
>>> ele.name
'Carbon'
>>> ele.description[:21]
'Carbon is a member of'
>>> ele.eleconfig
'[He] 2s2 2p2'
>>> ele.eleconfig_dict
{(1, 's'): 2, (2, 's'): 2, (2, 'p'): 2}
>>> str(ELEMENTS[6])
'Carbon'
>>> len(ELEMENTS)
109
>>> sum(ele.mass for ele in ELEMENTS)
14693.181589001004
>>> for ele in ELEMENTS:
...     ele.validate()
...     ele = eval(repr(ele))

"""

__version__ = '2021.6.18'

__all__ = (
    'ELEMENTS',
    'Element',
    'Isotope',
    'Particle',
    'ELEMENTARY_CHARGE',
    'ELECTRON',
    'PROTON',
    'NEUTRON',
    'POSITRON',
    'sqlite_script',
)


class lazyattr:
    """Attribute whose value is computed on first access.

    Lazyattrs are not thread-safe.

    """

    # TODO: replace with functools.cached_property? requires Python >= 3.8
    __slots__ = ('func', '__dict__')

    def __init__(self, func):
        """Initialize instance from decorated function."""
        self.func = func
        self.__doc__ = func.__doc__
        self.__module__ = func.__module__
        self.__name__ = func.__name__
        self.__qualname__ = func.__qualname__
        # self.lock = threading.RLock()

    def __get__(self, instance, owner):
        # with self.lock:
        if instance is None:
            return self
        try:
            value = self.func(instance)
        except AttributeError as exc:
            raise RuntimeError(exc)
        if value is NotImplemented:
            return getattr(super(owner, instance), self.func.__name__)
        setattr(instance, self.func.__name__, value)
        return value


class Particle:
    """Particle, e.g. electron, proton, or neutron.

    Attributes
    ----------
    name : str
        Name in English.
    mass : float
        Relative mass. Ratio of the average mass of atoms of the particle
        to 1/12 of the mass of an atom of 12C.
    charge : float
        Electric charge in coulomb.

    """

    __slots__ = ('name', 'mass', 'charge')

    def __init__(self, name, mass, charge):
        self.name = name
        self.mass = float(mass)
        self.charge = float(charge)


class Element:
    """Chemical element.

    Attributes
    ----------
    number : int
        Atomic number.
    symbol : str of length 1 or 2
        Chemical symbol.
    name : str
        Name in English.
    group : int
        Group in periodic table.
    period : int
        Period in periodic table.
    block : int
        Block in periodic table.
    series : int
        Index to chemical series.
    protons : int
        Number of protons.
    neutrons : int
        Number of neutrons in the most abundant naturally occurring stable.
        isotope
    nominalmass : int
        Mass number of the most abundant naturally occurring stable isotope.
    electrons : int
        Number of electrons.
    mass : float
        Relative atomic mass. Ratio of the average mass of atoms
        of the element to 1/12 of the mass of an atom of 12C.
    exactmass : float
        Relative atomic mass calculated from the isotopic composition.
    eleneg : float
        Electronegativity (Pauling scale).
    covrad : float
        Covalent radius in Angstrom.
    atmrad : float
        Atomic radius in Angstrom.
    vdwrad : float
        Van der Waals radius in Angstrom.
    tboil : float
        Boiling temperature in K.
    tmelt : float
        Melting temperature in K.
    density : float
        Density at 295K in g/cm3 respectively g/L.
    oxistates : str
        Oxidation states.
    eleaffin : float
        Electron affinity in eV.
    eleconfig : str
        Ground state electron configuration.
    eleconfig_dict : dict
        Ground state electron configuration (shell, subshell): electrons.
    eleshells : int
        Number of electrons per shell.
    ionenergy : tuple
        Ionization energies in eV
    isotopes : dict
        Isotopic composition.
        keys: isotope mass number
        values: Isotope(relative atomic mass, abundance)

    """

    def __init__(self, number, symbol, name, **kwargs):
        self.number = number
        self.symbol = symbol
        self.name = name
        self.electrons = number
        self.protons = number
        self.group = 0
        self.period = 0
        self.block = ''
        self.series = 0
        self.mass = 0.0
        self.eleneg = 0.0
        self.eleaffin = 0.0
        self.covrad = 0.0
        self.atmrad = 0.0
        self.vdwrad = 0.0
        self.tboil = 0.0
        self.tmelt = 0.0
        self.density = 0.0
        self.eleconfig = ''
        self.oxistates = ''
        self.ionenergy = ()
        self.isotopes = {}
        self.__dict__.update(kwargs)

    def __str__(self):
        return self.name

    def __repr__(self):
        ionenergy = []
        for i, j in enumerate(self.ionenergy):
            if i and (i % 5 == 0):
                ionenergy.append(f'\n        {j}')
            else:
                ionenergy.append(f'{j}')
        ionenergy = ', '.join(ionenergy)
        if len(self.ionenergy) > 5:
            ionenergy = f'(\n        {ionenergy},\n    )'
        elif len(self.ionenergy) == 1:
            ionenergy = f'({ionenergy},)'
        else:
            ionenergy = f'({ionenergy})'

        isotopes = []
        for massnum in sorted(self.isotopes):
            iso = self.isotopes[massnum]
            isotopes.append(
                '{0}: Isotope({1}, {2}, {0})'.format(
                    massnum, iso.mass, iso.abundance
                )
            )
        isotopes = ',\n        '.join(isotopes)
        if len(self.isotopes) > 1:
            isotopes = f'{{\n        {isotopes},\n    }},'
        else:
            isotopes = f'{{{isotopes}}},'

        return ',\n    '.join(
            (
                f"Element(\n    {self.number}, '{self.symbol}', '{self.name}'",
                f"group={self.group}, period={self.period},"
                f" block='{self.block}', series={self.series}",
                f"mass={self.mass}, eleneg={self.eleneg},"
                f" eleaffin={self.eleaffin}",
                f"covrad={self.covrad}, atmrad={self.atmrad},"
                f" vdwrad={self.vdwrad}",
                f"tboil={self.tboil}, tmelt={self.tmelt}, density={self.density}",
                f"eleconfig='{self.eleconfig}'",
                f"oxistates='{self.oxistates}'",
                f"ionenergy={ionenergy}",
                f"isotopes={isotopes}\n)",
            )
        )

    @lazyattr
    def nominalmass(self):
        """Return mass number of most abundant natural stable isotope."""
        nominalmass = 0
        maxabundance = 0
        for massnum, iso in self.isotopes.items():
            if iso.abundance > maxabundance:
                maxabundance = iso.abundance
                nominalmass = massnum
        return nominalmass

    @lazyattr
    def neutrons(self):
        """Return number neutrons in most abundant natural stable isotope."""
        return self.nominalmass - self.protons

    @lazyattr
    def exactmass(self):
        """Return relative atomic mass calculated from isotopic composition."""
        return sum(iso.mass * iso.abundance for iso in self.isotopes.values())

    @lazyattr
    def eleconfig_dict(self):
        """Return electron configuration as dict."""
        adict = {}
        if self.eleconfig.startswith('['):
            base = self.eleconfig.split(' ', 1)[0][1:-1]
            adict.update(ELEMENTS[base].eleconfig_dict)
        for e in self.eleconfig.split()[bool(adict) :]:
            adict[(int(e[0]), e[1])] = int(e[2:]) if len(e) > 2 else 1
        return adict

    @lazyattr
    def eleshells(self):
        """Return number of electrons in shell as tuple."""
        eleshells = [0, 0, 0, 0, 0, 0, 0]
        for key, val in self.eleconfig_dict.items():
            eleshells[key[0] - 1] += val
        return tuple(ele for ele in eleshells if ele)

    @lazyattr
    def description(self):
        """Return text description of element."""
        try:
            from .elements_descriptions import elements_descriptions
        except ImportError:
            try:
                from elements_descriptions import elements_descriptions
            except ImportError:
                return ''

        return elements_descriptions(ELEMENTS, self.symbol)

    def validate(self):
        """Check consistency of data. Raise Error on failure."""
        assert self.period in PERIODS
        assert self.group in GROUPS
        assert self.block in BLOCKS
        assert self.series in SERIES

        if self.number != self.protons:
            raise ValueError(
                f'{self.symbol} - atomic number must equal proton number'
            )
        if self.protons != sum(self.eleshells):
            raise ValueError(
                f'{self.symbol} - number of protons must equal electrons'
            )
        if len(self.ionenergy) > 1:
            ionev_ = self.ionenergy[0]
            for ionev in self.ionenergy[1:]:
                if ionev <= ionev_:
                    raise ValueError(
                        f'{self.symbol} - ionenergy not increasing'
                    )
                ionev_ = ionev

        mass = 0.0
        frac = 0.0
        for iso in self.isotopes.values():
            mass += iso.abundance * iso.mass
            frac += iso.abundance
        if abs(mass - self.mass) > 0.03:
            raise ValueError(
                f'{self.symbol} - average of isotope masses '
                f'({mass:.4f}) != mass ({self.mass:.4f})'
            )
        if abs(frac - 1.0) > 1e-9:
            raise ValueError(
                f'{self.symbol} - sum of isotope abundances != 1.0'
            )


class Isotope:
    """Isotope massnumber, relative atomic mass, and abundance."""

    __slots__ = ('massnumber', 'mass', 'abundance')

    def __init__(self, mass=0.0, abundance=1.0, massnumber=0):
        self.mass = mass
        self.abundance = abundance
        self.massnumber = massnumber

    def __str__(self):
        return '{}, {:.4f}, {:.6f}%'.format(
            self.massnumber, self.mass, self.abundance * 100
        )

    def __repr__(self):
        return 'Isotope({}, {}, {})'.format(
            repr(self.mass), repr(self.abundance), repr(self.massnumber)
        )


class Elements:
    """Ordered dict of Elements with lookup by number, symbol, and name."""

    def __init__(self, *elements):
        self._list = []
        self._dict = {}
        for element in elements:
            if element.number > len(self._list) + 1:
                raise ValueError('Elements must be added in order')
            if element.number <= len(self._list):
                self._list[element.number - 1] = element
            else:
                self._list.append(element)
            self._dict[element.number] = element
            self._dict[element.symbol] = element
            self._dict[element.name] = element

    def __str__(self):
        return '[{}]'.format(', '.join(ele.symbol for ele in self._list))

    def __repr__(self):
        elements = ',\n    '.join(
            '\n    '.join(line for line in repr(element).splitlines())
            for element in self._list
        )
        elements = f'Elements(\n    {elements},\n)'
        return elements

    def __contains__(self, item):
        return item in self._dict

    def __iter__(self):
        return iter(self._list)

    def __len__(self):
        return len(self._list)

    def __getitem__(self, key):
        try:
            return self._dict[key]
        except KeyError:
            try:
                start, stop, step = key.indices(len(self._list))
                return self._list[slice(start - 1, stop - 1, step)]
            except Exception:
                raise KeyError


# fmt: off
ELEMENTS = Elements(
    Element(
        1, 'H', 'Hydrogen',
        group=1, period=1, block='s', series=1,
        mass=1.007941, eleneg=2.2, eleaffin=0.75420375,
        covrad=0.32, atmrad=0.79, vdwrad=1.2,
        tboil=20.28, tmelt=13.81, density=0.084,
        eleconfig='1s',
        oxistates='1*, -1',
        ionenergy=(13.5984,),
        isotopes={
            1: Isotope(1.00782503223, 0.999885, 1),
            2: Isotope(2.01410177812, 0.000115, 2),
        },
    ),
    Element(
        2, 'He', 'Helium',
        group=18, period=1, block='s', series=2,
        mass=4.002602, eleneg=0.0, eleaffin=0.0,
        covrad=0.93, atmrad=0.49, vdwrad=1.4,
        tboil=4.216, tmelt=0.95, density=0.1785,
        eleconfig='1s2',
        oxistates='*',
        ionenergy=(24.5874, 54.416),
        isotopes={
            3: Isotope(3.0160293201, 1.34e-06, 3),
            4: Isotope(4.00260325413, 0.99999866, 4),
        },
    ),
    Element(
        3, 'Li', 'Lithium',
        group=1, period=2, block='s', series=3,
        mass=6.94, eleneg=0.98, eleaffin=0.618049,
        covrad=1.23, atmrad=2.05, vdwrad=1.82,
        tboil=1615.0, tmelt=453.7, density=0.53,
        eleconfig='[He] 2s',
        oxistates='1*',
        ionenergy=(5.3917, 75.638, 122.451),
        isotopes={
            6: Isotope(6.0151228874, 0.0759, 6),
            7: Isotope(7.0160034366, 0.9241, 7),
        },
    ),
    Element(
        4, 'Be', 'Beryllium',
        group=2, period=2, block='s', series=4,
        mass=9.0121831, eleneg=1.57, eleaffin=0.0,
        covrad=0.9, atmrad=1.4, vdwrad=0.0,
        tboil=3243.0, tmelt=1560.0, density=1.85,
        eleconfig='[He] 2s2',
        oxistates='2*',
        ionenergy=(9.3227, 18.211, 153.893, 217.713),
        isotopes={9: Isotope(9.012183065, 1.0, 9)},
    ),
    Element(
        5, 'B', 'Boron',
        group=13, period=2, block='p', series=5,
        mass=10.811, eleneg=2.04, eleaffin=0.279723,
        covrad=0.82, atmrad=1.17, vdwrad=0.0,
        tboil=4275.0, tmelt=2365.0, density=2.46,
        eleconfig='[He] 2s2 2p',
        oxistates='3*',
        ionenergy=(8.298, 25.154, 37.93, 59.368, 340.217),
        isotopes={
            10: Isotope(10.01293695, 0.199, 10),
            11: Isotope(11.00930536, 0.801, 11),
        },
    ),
    Element(
        6, 'C', 'Carbon',
        group=14, period=2, block='p', series=1,
        mass=12.01074, eleneg=2.55, eleaffin=1.262118,
        covrad=0.77, atmrad=0.91, vdwrad=1.7,
        tboil=5100.0, tmelt=3825.0, density=3.51,
        eleconfig='[He] 2s2 2p2',
        oxistates='4*, 2, -4*',
        ionenergy=(
            11.2603, 24.383, 47.877, 64.492, 392.077,
            489.981,
        ),
        isotopes={
            12: Isotope(12.0, 0.9893, 12),
            13: Isotope(13.00335483507, 0.0107, 13),
        },
    ),
    Element(
        7, 'N', 'Nitrogen',
        group=15, period=2, block='p', series=1,
        mass=14.006703, eleneg=3.04, eleaffin=-0.07,
        covrad=0.75, atmrad=0.75, vdwrad=1.55,
        tboil=77.344, tmelt=63.15, density=1.17,
        eleconfig='[He] 2s2 2p3',
        oxistates='5, 4, 3, 2, -3*',
        ionenergy=(
            14.5341, 39.601, 47.488, 77.472, 97.888,
            522.057, 667.029,
        ),
        isotopes={
            14: Isotope(14.00307400443, 0.99636, 14),
            15: Isotope(15.00010889888, 0.00364, 15),
        },
    ),
    Element(
        8, 'O', 'Oxygen',
        group=16, period=2, block='p', series=1,
        mass=15.999405, eleneg=3.44, eleaffin=1.461112,
        covrad=0.73, atmrad=0.65, vdwrad=1.52,
        tboil=90.188, tmelt=54.8, density=1.33,
        eleconfig='[He] 2s2 2p4',
        oxistates='-2*, -1',
        ionenergy=(
            13.6181, 35.116, 54.934, 77.412, 113.896,
            138.116, 739.315, 871.387,
        ),
        isotopes={
            16: Isotope(15.99491461957, 0.99757, 16),
            17: Isotope(16.9991317565, 0.00038, 17),
            18: Isotope(17.99915961286, 0.00205, 18),
        },
    ),
    Element(
        9, 'F', 'Fluorine',
        group=17, period=2, block='p', series=6,
        mass=18.998403163, eleneg=3.98, eleaffin=3.4011887,
        covrad=0.72, atmrad=0.57, vdwrad=1.47,
        tboil=85.0, tmelt=53.55, density=1.58,
        eleconfig='[He] 2s2 2p5',
        oxistates='-1*',
        ionenergy=(
            17.4228, 34.97, 62.707, 87.138, 114.24,
            157.161, 185.182, 953.886, 1103.089,
        ),
        isotopes={19: Isotope(18.99840316273, 1.0, 19)},
    ),
    Element(
        10, 'Ne', 'Neon',
        group=18, period=2, block='p', series=2,
        mass=20.1797, eleneg=0.0, eleaffin=0.0,
        covrad=0.71, atmrad=0.51, vdwrad=1.54,
        tboil=27.1, tmelt=24.55, density=0.8999,
        eleconfig='[He] 2s2 2p6',
        oxistates='*',
        ionenergy=(
            21.5645, 40.962, 63.45, 97.11, 126.21,
            157.93, 207.27, 239.09, 1195.797, 1362.164,
        ),
        isotopes={
            20: Isotope(19.9924401762, 0.9048, 20),
            21: Isotope(20.993846685, 0.0027, 21),
            22: Isotope(21.991385114, 0.0925, 22),
        },
    ),
    Element(
        11, 'Na', 'Sodium',
        group=1, period=3, block='s', series=3,
        mass=22.98976928, eleneg=0.93, eleaffin=0.547926,
        covrad=1.54, atmrad=2.23, vdwrad=2.27,
        tboil=1156.0, tmelt=371.0, density=0.97,
        eleconfig='[Ne] 3s',
        oxistates='1*',
        ionenergy=(
            5.1391, 47.286, 71.64, 98.91, 138.39,
            172.15, 208.47, 264.18, 299.87, 1465.091,
            1648.659,
        ),
        isotopes={23: Isotope(22.989769282, 1.0, 23)},
    ),
    Element(
        12, 'Mg', 'Magnesium',
        group=2, period=3, block='s', series=4,
        mass=24.3051, eleneg=1.31, eleaffin=0.0,
        covrad=1.36, atmrad=1.72, vdwrad=1.73,
        tboil=1380.0, tmelt=922.0, density=1.74,
        eleconfig='[Ne] 3s2',
        oxistates='2*',
        ionenergy=(
            7.6462, 15.035, 80.143, 109.24, 141.26,
            186.5, 224.94, 265.9, 327.95, 367.53,
            1761.802, 1962.613,
        ),
        isotopes={
            24: Isotope(23.985041697, 0.7899, 24),
            25: Isotope(24.985836976, 0.1, 25),
            26: Isotope(25.982592968, 0.1101, 26),
        },
    ),
    Element(
        13, 'Al', 'Aluminium',
        group=13, period=3, block='p', series=7,
        mass=26.9815385, eleneg=1.61, eleaffin=0.43283,
        covrad=1.18, atmrad=1.82, vdwrad=0.0,
        tboil=2740.0, tmelt=933.5, density=2.7,
        eleconfig='[Ne] 3s2 3p',
        oxistates='3*',
        ionenergy=(
            5.9858, 18.828, 28.447, 119.99, 153.71,
            190.47, 241.43, 284.59, 330.21, 398.57,
            442.07, 2085.983, 2304.08,
        ),
        isotopes={27: Isotope(26.98153853, 1.0, 27)},
    ),
    Element(
        14, 'Si', 'Silicon',
        group=14, period=3, block='p', series=5,
        mass=28.0855, eleneg=1.9, eleaffin=1.389521,
        covrad=1.11, atmrad=1.46, vdwrad=2.1,
        tboil=2630.0, tmelt=1683.0, density=2.33,
        eleconfig='[Ne] 3s2 3p2',
        oxistates='4*, -4',
        ionenergy=(
            8.1517, 16.345, 33.492, 45.141, 166.77,
            205.05, 246.52, 303.17, 351.1, 401.43,
            476.06, 523.5, 2437.676, 2673.108,
        ),
        isotopes={
            28: Isotope(27.97692653465, 0.92223, 28),
            29: Isotope(28.9764946649, 0.04685, 29),
            30: Isotope(29.973770136, 0.03092, 30),
        },
    ),
    Element(
        15, 'P', 'Phosphorus',
        group=15, period=3, block='p', series=1,
        mass=30.973761998, eleneg=2.19, eleaffin=0.7465,
        covrad=1.06, atmrad=1.23, vdwrad=1.8,
        tboil=553.0, tmelt=317.3, density=1.82,
        eleconfig='[Ne] 3s2 3p3',
        oxistates='5*, 3, -3',
        ionenergy=(
            10.4867, 19.725, 30.18, 51.37, 65.023,
            220.43, 263.22, 309.41, 371.73, 424.5,
            479.57, 560.41, 611.85, 2816.943, 3069.762,
        ),
        isotopes={31: Isotope(30.97376199842, 1.0, 31)},
    ),
    Element(
        16, 'S', 'Sulfur',
        group=16, period=3, block='p', series=1,
        mass=32.0648, eleneg=2.58, eleaffin=2.0771029,
        covrad=1.02, atmrad=1.09, vdwrad=1.8,
        tboil=717.82, tmelt=392.2, density=2.06,
        eleconfig='[Ne] 3s2 3p4',
        oxistates='6*, 4, 2, -2',
        ionenergy=(
            10.36, 23.33, 34.83, 47.3, 72.68,
            88.049, 280.93, 328.23, 379.1, 447.09,
            504.78, 564.65, 651.63, 707.14, 3223.836,
            3494.099,
        ),
        isotopes={
            32: Isotope(31.9720711744, 0.9499, 32),
            33: Isotope(32.9714589098, 0.0075, 33),
            34: Isotope(33.967867004, 0.0425, 34),
            36: Isotope(35.96708071, 0.0001, 36),
        },
    ),
    Element(
        17, 'Cl', 'Chlorine',
        group=17, period=3, block='p', series=6,
        mass=35.4529, eleneg=3.16, eleaffin=3.612724,
        covrad=0.99, atmrad=0.97, vdwrad=1.75,
        tboil=239.18, tmelt=172.17, density=2.95,
        eleconfig='[Ne] 3s2 3p5',
        oxistates='7, 5, 3, 1, -1*',
        ionenergy=(
            12.9676, 23.81, 39.61, 53.46, 67.8,
            98.03, 114.193, 348.28, 400.05, 455.62,
            529.97, 591.97, 656.69, 749.75, 809.39,
            3658.425, 3946.193,
        ),
        isotopes={
            35: Isotope(34.968852682, 0.7576, 35),
            37: Isotope(36.965902602, 0.2424, 37),
        },
    ),
    Element(
        18, 'Ar', 'Argon',
        group=18, period=3, block='p', series=2,
        mass=39.948, eleneg=0.0, eleaffin=0.0,
        covrad=0.98, atmrad=0.88, vdwrad=1.88,
        tboil=87.45, tmelt=83.95, density=1.66,
        eleconfig='[Ne] 3s2 3p6',
        oxistates='*',
        ionenergy=(
            15.7596, 27.629, 40.74, 59.81, 75.02,
            91.007, 124.319, 143.456, 422.44, 478.68,
            538.95, 618.24, 686.09, 755.73, 854.75,
            918.0, 4120.778, 4426.114,
        ),
        isotopes={
            36: Isotope(35.967545105, 0.003336, 36),
            38: Isotope(37.96273211, 0.000629, 38),
            40: Isotope(39.9623831237, 0.996035, 40),
        },
    ),
    Element(
        19, 'K', 'Potassium',
        group=1, period=4, block='s', series=3,
        mass=39.0983, eleneg=0.82, eleaffin=0.501459,
        covrad=2.03, atmrad=2.77, vdwrad=2.75,
        tboil=1033.0, tmelt=336.8, density=0.86,
        eleconfig='[Ar] 4s',
        oxistates='1*',
        ionenergy=(
            4.3407, 31.625, 45.72, 60.91, 82.66,
            100.0, 117.56, 154.86, 175.814, 503.44,
            564.13, 629.09, 714.02, 787.13, 861.77,
            968.0, 1034.0, 4610.955, 4933.931,
        ),
        isotopes={
            39: Isotope(38.9637064864, 0.932581, 39),
            40: Isotope(39.963998166, 0.000117, 40),
            41: Isotope(40.9618252579, 0.067302, 41),
        },
    ),
    Element(
        20, 'Ca', 'Calcium',
        group=2, period=4, block='s', series=4,
        mass=40.078, eleneg=1.0, eleaffin=0.02455,
        covrad=1.74, atmrad=2.23, vdwrad=0.0,
        tboil=1757.0, tmelt=1112.0, density=1.54,
        eleconfig='[Ar] 4s2',
        oxistates='2*',
        ionenergy=(
            6.1132, 11.71, 50.908, 67.1, 84.41,
            108.78, 127.7, 147.24, 188.54, 211.27,
            591.25, 656.39, 726.03, 816.61, 895.12,
            974.0, 1087.0, 1157.0, 5129.045, 5469.738,
        ),
        isotopes={
            40: Isotope(39.962590863, 0.96941, 40),
            42: Isotope(41.95861783, 0.00647, 42),
            43: Isotope(42.95876644, 0.00135, 43),
            44: Isotope(43.95548156, 0.02086, 44),
            46: Isotope(45.953689, 4e-05, 46),
            48: Isotope(47.95252276, 0.00187, 48),
        },
    ),
    Element(
        21, 'Sc', 'Scandium',
        group=3, period=4, block='d', series=8,
        mass=44.955908, eleneg=1.36, eleaffin=0.188,
        covrad=1.44, atmrad=2.09, vdwrad=0.0,
        tboil=3109.0, tmelt=1814.0, density=2.99,
        eleconfig='[Ar] 3d 4s2',
        oxistates='3*',
        ionenergy=(
            6.5615, 12.8, 24.76, 73.47, 91.66,
            110.1, 138.0, 158.7, 180.02, 225.32,
            249.8, 685.89, 755.47, 829.79, 926.0,
        ),
        isotopes={45: Isotope(44.95590828, 1.0, 45)},
    ),
    Element(
        22, 'Ti', 'Titanium',
        group=4, period=4, block='d', series=8,
        mass=47.867, eleneg=1.54, eleaffin=0.084,
        covrad=1.32, atmrad=2.0, vdwrad=0.0,
        tboil=3560.0, tmelt=1935.0, density=4.51,
        eleconfig='[Ar] 3d2 4s2',
        oxistates='4*, 3',
        ionenergy=(
            6.8281, 13.58, 27.491, 43.266, 99.22,
            119.36, 140.8, 168.5, 193.5, 215.91,
            265.23, 291.497, 787.33, 861.33,
        ),
        isotopes={
            46: Isotope(45.95262772, 0.0825, 46),
            47: Isotope(46.95175879, 0.0744, 47),
            48: Isotope(47.94794198, 0.7372, 48),
            49: Isotope(48.94786568, 0.0541, 49),
            50: Isotope(49.94478689, 0.0518, 50),
        },
    ),
    Element(
        23, 'V', 'Vanadium',
        group=5, period=4, block='d', series=8,
        mass=50.9415, eleneg=1.63, eleaffin=0.525,
        covrad=1.22, atmrad=1.92, vdwrad=0.0,
        tboil=3650.0, tmelt=2163.0, density=6.09,
        eleconfig='[Ar] 3d3 4s2',
        oxistates='5*, 4, 3, 2, 0',
        ionenergy=(
            6.7462, 14.65, 29.31, 46.707, 65.23,
            128.12, 150.17, 173.7, 205.8, 230.5,
            255.04, 308.25, 336.267, 895.58, 974.02,
        ),
        isotopes={
            50: Isotope(49.94715601, 0.0025, 50),
            51: Isotope(50.94395704, 0.9975, 51),
        },
    ),
    Element(
        24, 'Cr', 'Chromium',
        group=6, period=4, block='d', series=8,
        mass=51.9961, eleneg=1.66, eleaffin=0.67584,
        covrad=1.18, atmrad=1.85, vdwrad=0.0,
        tboil=2945.0, tmelt=2130.0, density=7.14,
        eleconfig='[Ar] 3d5 4s',
        oxistates='6, 3*, 2, 0',
        ionenergy=(
            6.7665, 16.5, 30.96, 49.1, 69.3,
            90.56, 161.1, 184.7, 209.3, 244.4,
            270.8, 298.0, 355.0, 384.3, 1010.64,
        ),
        isotopes={
            50: Isotope(49.94604183, 0.04345, 50),
            52: Isotope(51.94050623, 0.83789, 52),
            53: Isotope(52.94064815, 0.09501, 53),
            54: Isotope(53.93887916, 0.02365, 54),
        },
    ),
    Element(
        25, 'Mn', 'Manganese',
        group=7, period=4, block='d', series=8,
        mass=54.938044, eleneg=1.55, eleaffin=0.0,
        covrad=1.17, atmrad=1.79, vdwrad=0.0,
        tboil=2235.0, tmelt=1518.0, density=7.44,
        eleconfig='[Ar] 3d5 4s2',
        oxistates='7, 6, 4, 3, 2*, 0, -1',
        ionenergy=(
            7.434, 15.64, 33.667, 51.2, 72.4,
            95.0, 119.27, 196.46, 221.8, 248.3,
            286.0, 314.4, 343.6, 404.0, 435.3,
            1136.2,
        ),
        isotopes={55: Isotope(54.93804391, 1.0, 55)},
    ),
    Element(
        26, 'Fe', 'Iron',
        group=8, period=4, block='d', series=8,
        mass=55.845, eleneg=1.83, eleaffin=0.151,
        covrad=1.17, atmrad=1.72, vdwrad=0.0,
        tboil=3023.0, tmelt=1808.0, density=7.874,
        eleconfig='[Ar] 3d6 4s2',
        oxistates='6, 3*, 2, 0, -2',
        ionenergy=(
            7.9024, 16.18, 30.651, 54.8, 75.0,
            99.0, 125.0, 151.06, 235.04, 262.1,
            290.4, 330.8, 361.0, 392.2, 457.0,
            485.5, 1266.1,
        ),
        isotopes={
            54: Isotope(53.93960899, 0.05845, 54),
            56: Isotope(55.93493633, 0.91754, 56),
            57: Isotope(56.93539284, 0.02119, 57),
            58: Isotope(57.93327443, 0.00282, 58),
        },
    ),
    Element(
        27, 'Co', 'Cobalt',
        group=9, period=4, block='d', series=8,
        mass=58.933194, eleneg=1.88, eleaffin=0.6633,
        covrad=1.16, atmrad=1.67, vdwrad=0.0,
        tboil=3143.0, tmelt=1768.0, density=8.89,
        eleconfig='[Ar] 3d7 4s2',
        oxistates='3, 2*, 0, -1',
        ionenergy=(
            7.881, 17.06, 33.5, 51.3, 79.5,
            102.0, 129.0, 157.0, 186.13, 276.0,
            305.0, 336.0, 376.0, 411.0, 444.0,
            512.0, 546.8, 1403.0,
        ),
        isotopes={59: Isotope(58.93319429, 1.0, 59)},
    ),
    Element(
        28, 'Ni', 'Nickel',
        group=10, period=4, block='d', series=8,
        mass=58.6934, eleneg=1.91, eleaffin=1.15716,
        covrad=1.15, atmrad=1.62, vdwrad=1.63,
        tboil=3005.0, tmelt=1726.0, density=8.91,
        eleconfig='[Ar] 3d8 4s2',
        oxistates='3, 2*, 0',
        ionenergy=(
            7.6398, 18.168, 35.17, 54.9, 75.5,
            108.0, 133.0, 162.0, 193.0, 224.5,
            321.2, 352.0, 384.0, 430.0, 464.0,
            499.0, 571.0, 607.2, 1547.0,
        ),
        isotopes={
            58: Isotope(57.93534241, 0.68077, 58),
            60: Isotope(59.93078588, 0.26223, 60),
            61: Isotope(60.93105557, 0.011399, 61),
            62: Isotope(61.92834537, 0.036346, 62),
            64: Isotope(63.92796682, 0.009255, 64),
        },
    ),
    Element(
        29, 'Cu', 'Copper',
        group=11, period=4, block='d', series=8,
        mass=63.546, eleneg=1.9, eleaffin=1.23578,
        covrad=1.17, atmrad=1.57, vdwrad=1.4,
        tboil=2840.0, tmelt=1356.6, density=8.92,
        eleconfig='[Ar] 3d10 4s',
        oxistates='2*, 1',
        ionenergy=(
            7.7264, 20.292, 26.83, 55.2, 79.9,
            103.0, 139.0, 166.0, 199.0, 232.0,
            266.0, 368.8, 401.0, 435.0, 484.0,
            520.0, 557.0, 633.0, 671.0, 1698.0,
        ),
        isotopes={
            63: Isotope(62.92959772, 0.6915, 63),
            65: Isotope(64.9277897, 0.3085, 65),
        },
    ),
    Element(
        30, 'Zn', 'Zinc',
        group=12, period=4, block='d', series=8,
        mass=65.38, eleneg=1.65, eleaffin=0.0,
        covrad=1.25, atmrad=1.53, vdwrad=1.39,
        tboil=1180.0, tmelt=692.73, density=7.14,
        eleconfig='[Ar] 3d10 4s2',
        oxistates='2*',
        ionenergy=(
            9.3942, 17.964, 39.722, 59.4, 82.6,
            108.0, 134.0, 174.0, 203.0, 238.0,
            274.0, 310.8, 419.7, 454.0, 490.0,
            542.0, 579.0, 619.0, 698.8, 738.0,
            1856.0,
        ),
        isotopes={
            64: Isotope(63.92914201, 0.4917, 64),
            66: Isotope(65.92603381, 0.2773, 66),
            67: Isotope(66.92712775, 0.0404, 67),
            68: Isotope(67.92484455, 0.1845, 68),
            70: Isotope(69.9253192, 0.0061, 70),
        },
    ),
    Element(
        31, 'Ga', 'Gallium',
        group=13, period=4, block='p', series=7,
        mass=69.723, eleneg=1.81, eleaffin=0.41,
        covrad=1.26, atmrad=1.81, vdwrad=1.87,
        tboil=2478.0, tmelt=302.92, density=5.91,
        eleconfig='[Ar] 3d10 4s2 4p',
        oxistates='3*',
        ionenergy=(5.9993, 20.51, 30.71, 64.0),
        isotopes={
            69: Isotope(68.9255735, 0.60108, 69),
            71: Isotope(70.92470258, 0.39892, 71),
        },
    ),
    Element(
        32, 'Ge', 'Germanium',
        group=14, period=4, block='p', series=5,
        mass=72.63, eleneg=2.01, eleaffin=1.232712,
        covrad=1.22, atmrad=1.52, vdwrad=0.0,
        tboil=3107.0, tmelt=1211.5, density=5.32,
        eleconfig='[Ar] 3d10 4s2 4p2',
        oxistates='4*',
        ionenergy=(7.8994, 15.934, 34.22, 45.71, 93.5),
        isotopes={
            70: Isotope(69.92424875, 0.2057, 70),
            72: Isotope(71.922075826, 0.2745, 72),
            73: Isotope(72.923458956, 0.0775, 73),
            74: Isotope(73.921177761, 0.365, 74),
            76: Isotope(75.921402726, 0.0773, 76),
        },
    ),
    Element(
        33, 'As', 'Arsenic',
        group=15, period=4, block='p', series=5,
        mass=74.921595, eleneg=2.18, eleaffin=0.814,
        covrad=1.2, atmrad=1.33, vdwrad=1.85,
        tboil=876.0, tmelt=1090.0, density=5.72,
        eleconfig='[Ar] 3d10 4s2 4p3',
        oxistates='5, 3*, -3',
        ionenergy=(
            9.7886, 18.633, 28.351, 50.13, 62.63,
            127.6,
        ),
        isotopes={75: Isotope(74.92159457, 1.0, 75)},
    ),
    Element(
        34, 'Se', 'Selenium',
        group=16, period=4, block='p', series=1,
        mass=78.971, eleneg=2.55, eleaffin=2.02067,
        covrad=1.16, atmrad=1.22, vdwrad=1.9,
        tboil=958.0, tmelt=494.0, density=4.82,
        eleconfig='[Ar] 3d10 4s2 4p4',
        oxistates='6, 4*, -2',
        ionenergy=(
            9.7524, 21.9, 30.82, 42.944, 68.3,
            81.7, 155.4,
        ),
        isotopes={
            74: Isotope(73.922475934, 0.0089, 74),
            76: Isotope(75.919213704, 0.0937, 76),
            77: Isotope(76.919914154, 0.0763, 77),
            78: Isotope(77.91730928, 0.2377, 78),
            80: Isotope(79.9165218, 0.4961, 80),
            82: Isotope(81.9166995, 0.0873, 82),
        },
    ),
    Element(
        35, 'Br', 'Bromine',
        group=17, period=4, block='p', series=6,
        mass=79.9035, eleneg=2.96, eleaffin=3.363588,
        covrad=1.14, atmrad=1.12, vdwrad=1.85,
        tboil=331.85, tmelt=265.95, density=3.14,
        eleconfig='[Ar] 3d10 4s2 4p5',
        oxistates='7, 5, 3, 1, -1*',
        ionenergy=(
            11.8138, 21.8, 36.0, 47.3, 59.7,
            88.6, 103.0, 192.8,
        ),
        isotopes={
            79: Isotope(78.9183376, 0.5069, 79),
            81: Isotope(80.9162897, 0.4931, 81),
        },
    ),
    Element(
        36, 'Kr', 'Krypton',
        group=18, period=4, block='p', series=2,
        mass=83.798, eleneg=0.0, eleaffin=0.0,
        covrad=1.12, atmrad=1.03, vdwrad=2.02,
        tboil=120.85, tmelt=116.0, density=4.48,
        eleconfig='[Ar] 3d10 4s2 4p6',
        oxistates='2*',
        ionenergy=(
            13.9996, 24.359, 36.95, 52.5, 64.7,
            78.5, 110.0, 126.0, 230.39,
        ),
        isotopes={
            78: Isotope(77.92036494, 0.00355, 78),
            80: Isotope(79.91637808, 0.02286, 80),
            82: Isotope(81.91348273, 0.11593, 82),
            83: Isotope(82.91412716, 0.115, 83),
            84: Isotope(83.9114977282, 0.56987, 84),
            86: Isotope(85.9106106269, 0.17279, 86),
        },
    ),
    Element(
        37, 'Rb', 'Rubidium',
        group=1, period=5, block='s', series=3,
        mass=85.4678, eleneg=0.82, eleaffin=0.485916,
        covrad=2.16, atmrad=2.98, vdwrad=0.0,
        tboil=961.0, tmelt=312.63, density=1.53,
        eleconfig='[Kr] 5s',
        oxistates='1*',
        ionenergy=(
            4.1771, 27.28, 40.0, 52.6, 71.0,
            84.4, 99.2, 136.0, 150.0, 277.1,
        ),
        isotopes={
            85: Isotope(84.9117897379, 0.7217, 85),
            87: Isotope(86.909180531, 0.2783, 87),
        },
    ),
    Element(
        38, 'Sr', 'Strontium',
        group=2, period=5, block='s', series=4,
        mass=87.62, eleneg=0.95, eleaffin=0.05206,
        covrad=1.91, atmrad=2.45, vdwrad=0.0,
        tboil=1655.0, tmelt=1042.0, density=2.63,
        eleconfig='[Kr] 5s2',
        oxistates='2*',
        ionenergy=(
            5.6949, 11.03, 43.6, 57.0, 71.6,
            90.8, 106.0, 122.3, 162.0, 177.0,
            324.1,
        ),
        isotopes={
            84: Isotope(83.9134191, 0.0056, 84),
            86: Isotope(85.9092606, 0.0986, 86),
            87: Isotope(86.9088775, 0.07, 87),
            88: Isotope(87.9056125, 0.8258, 88),
        },
    ),
    Element(
        39, 'Y', 'Yttrium',
        group=3, period=5, block='d', series=8,
        mass=88.90584, eleneg=1.22, eleaffin=0.307,
        covrad=1.62, atmrad=2.27, vdwrad=0.0,
        tboil=3611.0, tmelt=1795.0, density=4.47,
        eleconfig='[Kr] 4d 5s2',
        oxistates='3*',
        ionenergy=(
            6.2173, 12.24, 20.52, 61.8, 77.0,
            93.0, 116.0, 129.0, 146.52, 191.0,
            206.0, 374.0,
        ),
        isotopes={89: Isotope(88.9058403, 1.0, 89)},
    ),
    Element(
        40, 'Zr', 'Zirconium',
        group=4, period=5, block='d', series=8,
        mass=91.224, eleneg=1.33, eleaffin=0.426,
        covrad=1.45, atmrad=2.16, vdwrad=0.0,
        tboil=4682.0, tmelt=2128.0, density=6.51,
        eleconfig='[Kr] 4d2 5s2',
        oxistates='4*',
        ionenergy=(6.6339, 13.13, 22.99, 34.34, 81.5),
        isotopes={
            90: Isotope(89.9046977, 0.5145, 90),
            91: Isotope(90.9056396, 0.1122, 91),
            92: Isotope(91.9050347, 0.1715, 92),
            94: Isotope(93.9063108, 0.1738, 94),
            96: Isotope(95.9082714, 0.028, 96),
        },
    ),
    Element(
        41, 'Nb', 'Niobium',
        group=5, period=5, block='d', series=8,
        mass=92.90637, eleneg=1.6, eleaffin=0.893,
        covrad=1.34, atmrad=2.08, vdwrad=0.0,
        tboil=5015.0, tmelt=2742.0, density=8.58,
        eleconfig='[Kr] 4d4 5s',
        oxistates='5*, 3',
        ionenergy=(
            6.7589, 14.32, 25.04, 38.3, 50.55,
            102.6, 125.0,
        ),
        isotopes={93: Isotope(92.906373, 1.0, 93)},
    ),
    Element(
        42, 'Mo', 'Molybdenum',
        group=6, period=5, block='d', series=8,
        mass=95.95, eleneg=2.16, eleaffin=0.7472,
        covrad=1.3, atmrad=2.01, vdwrad=0.0,
        tboil=4912.0, tmelt=2896.0, density=10.28,
        eleconfig='[Kr] 4d5 5s',
        oxistates='6*, 5, 4, 3, 2, 0',
        ionenergy=(
            7.0924, 16.15, 27.16, 46.4, 61.2,
            68.0, 126.8, 153.0,
        ),
        isotopes={
            92: Isotope(91.90680796, 0.1453, 92),
            94: Isotope(93.9050849, 0.0915, 94),
            95: Isotope(94.90583877, 0.1584, 95),
            96: Isotope(95.90467612, 0.1667, 96),
            97: Isotope(96.90601812, 0.096, 97),
            98: Isotope(97.90540482, 0.2439, 98),
            100: Isotope(99.9074718, 0.0982, 100),
        },
    ),
    Element(
        43, 'Tc', 'Technetium',
        group=7, period=5, block='d', series=8,
        mass=97.9072, eleneg=1.9, eleaffin=0.55,
        covrad=1.27, atmrad=1.95, vdwrad=0.0,
        tboil=4538.0, tmelt=2477.0, density=11.49,
        eleconfig='[Kr] 4d5 5s2',
        oxistates='7*',
        ionenergy=(7.28, 15.26, 29.54),
        isotopes={98: Isotope(97.9072124, 1.0, 98)},
    ),
    Element(
        44, 'Ru', 'Ruthenium',
        group=8, period=5, block='d', series=8,
        mass=101.07, eleneg=2.2, eleaffin=1.04638,
        covrad=1.25, atmrad=1.89, vdwrad=0.0,
        tboil=4425.0, tmelt=2610.0, density=12.45,
        eleconfig='[Kr] 4d7 5s',
        oxistates='8, 6, 4*, 3*, 2, 0, -2',
        ionenergy=(7.3605, 16.76, 28.47),
        isotopes={
            96: Isotope(95.90759025, 0.0554, 96),
            98: Isotope(97.9052868, 0.0187, 98),
            99: Isotope(98.9059341, 0.1276, 99),
            100: Isotope(99.9042143, 0.126, 100),
            101: Isotope(100.9055769, 0.1706, 101),
            102: Isotope(101.9043441, 0.3155, 102),
            104: Isotope(103.9054275, 0.1862, 104),
        },
    ),
    Element(
        45, 'Rh', 'Rhodium',
        group=9, period=5, block='d', series=8,
        mass=102.9055, eleneg=2.28, eleaffin=1.14289,
        covrad=1.25, atmrad=1.83, vdwrad=0.0,
        tboil=3970.0, tmelt=2236.0, density=12.41,
        eleconfig='[Kr] 4d8 5s',
        oxistates='5, 4, 3*, 1*, 2, 0',
        ionenergy=(7.4589, 18.08, 31.06),
        isotopes={103: Isotope(102.905498, 1.0, 103)},
    ),
    Element(
        46, 'Pd', 'Palladium',
        group=10, period=5, block='d', series=8,
        mass=106.42, eleneg=2.2, eleaffin=0.56214,
        covrad=1.28, atmrad=1.79, vdwrad=1.63,
        tboil=3240.0, tmelt=1825.0, density=12.02,
        eleconfig='[Kr] 4d10',
        oxistates='4, 2*, 0',
        ionenergy=(8.3369, 19.43, 32.93),
        isotopes={
            102: Isotope(101.9056022, 0.0102, 102),
            104: Isotope(103.9040305, 0.1114, 104),
            105: Isotope(104.9050796, 0.2233, 105),
            106: Isotope(105.9034804, 0.2733, 106),
            108: Isotope(107.9038916, 0.2646, 108),
            110: Isotope(109.9051722, 0.1172, 110),
        },
    ),
    Element(
        47, 'Ag', 'Silver',
        group=11, period=5, block='d', series=8,
        mass=107.8682, eleneg=1.93, eleaffin=1.30447,
        covrad=1.34, atmrad=1.75, vdwrad=1.72,
        tboil=2436.0, tmelt=1235.1, density=10.49,
        eleconfig='[Kr] 4d10 5s',
        oxistates='2, 1*',
        ionenergy=(7.5762, 21.49, 34.83),
        isotopes={
            107: Isotope(106.9050916, 0.51839, 107),
            109: Isotope(108.9047553, 0.48161, 109),
        },
    ),
    Element(
        48, 'Cd', 'Cadmium',
        group=12, period=5, block='d', series=8,
        mass=112.414, eleneg=1.69, eleaffin=0.0,
        covrad=1.48, atmrad=1.71, vdwrad=1.58,
        tboil=1040.0, tmelt=594.26, density=8.64,
        eleconfig='[Kr] 4d10 5s2',
        oxistates='2*',
        ionenergy=(8.9938, 16.908, 37.48),
        isotopes={
            106: Isotope(105.9064599, 0.0125, 106),
            108: Isotope(107.9041834, 0.0089, 108),
            110: Isotope(109.90300661, 0.1249, 110),
            111: Isotope(110.90418287, 0.128, 111),
            112: Isotope(111.90276287, 0.2413, 112),
            113: Isotope(112.90440813, 0.1222, 113),
            114: Isotope(113.90336509, 0.2873, 114),
            116: Isotope(115.90476315, 0.0749, 116),
        },
    ),
    Element(
        49, 'In', 'Indium',
        group=13, period=5, block='p', series=7,
        mass=114.818, eleneg=1.78, eleaffin=0.404,
        covrad=1.44, atmrad=2.0, vdwrad=1.93,
        tboil=2350.0, tmelt=429.78, density=7.31,
        eleconfig='[Kr] 4d10 5s2 5p',
        oxistates='3*',
        ionenergy=(5.7864, 18.869, 28.03, 55.45),
        isotopes={
            113: Isotope(112.90406184, 0.0429, 113),
            115: Isotope(114.903878776, 0.9571, 115),
        },
    ),
    Element(
        50, 'Sn', 'Tin',
        group=14, period=5, block='p', series=7,
        mass=118.71, eleneg=1.96, eleaffin=1.112066,
        covrad=1.41, atmrad=1.72, vdwrad=2.17,
        tboil=2876.0, tmelt=505.12, density=7.29,
        eleconfig='[Kr] 4d10 5s2 5p2',
        oxistates='4*, 2*',
        ionenergy=(7.3439, 14.632, 30.502, 40.734, 72.28),
        isotopes={
            112: Isotope(111.90482387, 0.0097, 112),
            114: Isotope(113.9027827, 0.0066, 114),
            115: Isotope(114.903344699, 0.0034, 115),
            116: Isotope(115.9017428, 0.1454, 116),
            117: Isotope(116.90295398, 0.0768, 117),
            118: Isotope(117.90160657, 0.2422, 118),
            119: Isotope(118.90331117, 0.0859, 119),
            120: Isotope(119.90220163, 0.3258, 120),
            122: Isotope(121.9034438, 0.0463, 122),
            124: Isotope(123.9052766, 0.0579, 124),
        },
    ),
    Element(
        51, 'Sb', 'Antimony',
        group=15, period=5, block='p', series=5,
        mass=121.76, eleneg=2.05, eleaffin=1.047401,
        covrad=1.4, atmrad=1.53, vdwrad=0.0,
        tboil=1860.0, tmelt=903.91, density=6.69,
        eleconfig='[Kr] 4d10 5s2 5p3',
        oxistates='5, 3*, -3',
        ionenergy=(
            8.6084, 16.53, 25.3, 44.2, 56.0,
            108.0,
        ),
        isotopes={
            121: Isotope(120.903812, 0.5721, 121),
            123: Isotope(122.9042132, 0.4279, 123),
        },
    ),
    Element(
        52, 'Te', 'Tellurium',
        group=16, period=5, block='p', series=5,
        mass=127.6, eleneg=2.1, eleaffin=1.970875,
        covrad=1.36, atmrad=1.42, vdwrad=2.06,
        tboil=1261.0, tmelt=722.72, density=6.25,
        eleconfig='[Kr] 4d10 5s2 5p4',
        oxistates='6, 4*, -2',
        ionenergy=(
            9.0096, 18.6, 27.96, 37.41, 58.75,
            70.7, 137.0,
        ),
        isotopes={
            120: Isotope(119.9040593, 0.0009, 120),
            122: Isotope(121.9030435, 0.0255, 122),
            123: Isotope(122.9042698, 0.0089, 123),
            124: Isotope(123.9028171, 0.0474, 124),
            125: Isotope(124.9044299, 0.0707, 125),
            126: Isotope(125.9033109, 0.1884, 126),
            128: Isotope(127.90446128, 0.3174, 128),
            130: Isotope(129.906222748, 0.3408, 130),
        },
    ),
    Element(
        53, 'I', 'Iodine',
        group=17, period=5, block='p', series=6,
        mass=126.90447, eleneg=2.66, eleaffin=3.059038,
        covrad=1.33, atmrad=1.32, vdwrad=1.98,
        tboil=457.5, tmelt=386.7, density=4.94,
        eleconfig='[Kr] 4d10 5s2 5p5',
        oxistates='7, 5, 1, -1*',
        ionenergy=(10.4513, 19.131, 33.0),
        isotopes={127: Isotope(126.9044719, 1.0, 127)},
    ),
    Element(
        54, 'Xe', 'Xenon',
        group=18, period=5, block='p', series=2,
        mass=131.293, eleneg=0.0, eleaffin=0.0,
        covrad=1.31, atmrad=1.24, vdwrad=2.16,
        tboil=165.1, tmelt=161.39, density=4.49,
        eleconfig='[Kr] 4d10 5s2 5p6',
        oxistates='2, 4, 6',
        ionenergy=(12.1298, 21.21, 32.1),
        isotopes={
            124: Isotope(123.905892, 0.000952, 124),
            126: Isotope(125.9042983, 0.00089, 126),
            128: Isotope(127.903531, 0.019102, 128),
            129: Isotope(128.9047808611, 0.264006, 129),
            130: Isotope(129.903509349, 0.04071, 130),
            131: Isotope(130.90508406, 0.212324, 131),
            132: Isotope(131.9041550856, 0.269086, 132),
            134: Isotope(133.90539466, 0.104357, 134),
            136: Isotope(135.907214484, 0.088573, 136),
        },
    ),
    Element(
        55, 'Cs', 'Caesium',
        group=1, period=6, block='s', series=3,
        mass=132.90545196, eleneg=0.79, eleaffin=0.471626,
        covrad=2.35, atmrad=3.34, vdwrad=0.0,
        tboil=944.0, tmelt=301.54, density=1.9,
        eleconfig='[Xe] 6s',
        oxistates='1*',
        ionenergy=(3.8939, 25.1),
        isotopes={133: Isotope(132.905451961, 1.0, 133)},
    ),
    Element(
        56, 'Ba', 'Barium',
        group=2, period=6, block='s', series=4,
        mass=137.327, eleneg=0.89, eleaffin=0.14462,
        covrad=1.98, atmrad=2.78, vdwrad=0.0,
        tboil=2078.0, tmelt=1002.0, density=3.65,
        eleconfig='[Xe] 6s2',
        oxistates='2*',
        ionenergy=(5.2117, 100.004),
        isotopes={
            130: Isotope(129.9063207, 0.00106, 130),
            132: Isotope(131.9050611, 0.00101, 132),
            134: Isotope(133.90450818, 0.02417, 134),
            135: Isotope(134.90568838, 0.06592, 135),
            136: Isotope(135.90457573, 0.07854, 136),
            137: Isotope(136.90582714, 0.11232, 137),
            138: Isotope(137.905247, 0.71698, 138),
        },
    ),
    Element(
        57, 'La', 'Lanthanum',
        group=3, period=6, block='f', series=9,
        mass=138.90547, eleneg=1.1, eleaffin=0.47,
        covrad=1.69, atmrad=2.74, vdwrad=0.0,
        tboil=3737.0, tmelt=1191.0, density=6.16,
        eleconfig='[Xe] 5d 6s2',
        oxistates='3*',
        ionenergy=(5.5769, 11.06, 19.175),
        isotopes={
            138: Isotope(137.9071149, 0.0008881, 138),
            139: Isotope(138.9063563, 0.9991119, 139),
        },
    ),
    Element(
        58, 'Ce', 'Cerium',
        group=3, period=6, block='f', series=9,
        mass=140.116, eleneg=1.12, eleaffin=0.5,
        covrad=1.65, atmrad=2.7, vdwrad=0.0,
        tboil=3715.0, tmelt=1071.0, density=6.77,
        eleconfig='[Xe] 4f 5d 6s2',
        oxistates='4, 3*',
        ionenergy=(5.5387, 10.85, 20.2, 36.72),
        isotopes={
            136: Isotope(135.90712921, 0.00185, 136),
            138: Isotope(137.905991, 0.00251, 138),
            140: Isotope(139.9054431, 0.8845, 140),
            142: Isotope(141.9092504, 0.11114, 142),
        },
    ),
    Element(
        59, 'Pr', 'Praseodymium',
        group=3, period=6, block='f', series=9,
        mass=140.90766, eleneg=1.13, eleaffin=0.5,
        covrad=1.65, atmrad=2.67, vdwrad=0.0,
        tboil=3785.0, tmelt=1204.0, density=6.48,
        eleconfig='[Xe] 4f3 6s2',
        oxistates='4, 3*',
        ionenergy=(5.473, 10.55, 21.62, 38.95, 57.45),
        isotopes={141: Isotope(140.9076576, 1.0, 141)},
    ),
    Element(
        60, 'Nd', 'Neodymium',
        group=3, period=6, block='f', series=9,
        mass=144.242, eleneg=1.14, eleaffin=0.5,
        covrad=1.64, atmrad=2.64, vdwrad=0.0,
        tboil=3347.0, tmelt=1294.0, density=7.0,
        eleconfig='[Xe] 4f4 6s2',
        oxistates='3*',
        ionenergy=(5.525, 10.72),
        isotopes={
            142: Isotope(141.907729, 0.27152, 142),
            143: Isotope(142.90982, 0.12174, 143),
            144: Isotope(143.910093, 0.23798, 144),
            145: Isotope(144.9125793, 0.08293, 145),
            146: Isotope(145.9131226, 0.17189, 146),
            148: Isotope(147.9168993, 0.05756, 148),
            150: Isotope(149.9209022, 0.05638, 150),
        },
    ),
    Element(
        61, 'Pm', 'Promethium',
        group=3, period=6, block='f', series=9,
        mass=144.9128, eleneg=1.13, eleaffin=0.5,
        covrad=1.63, atmrad=2.62, vdwrad=0.0,
        tboil=3273.0, tmelt=1315.0, density=7.22,
        eleconfig='[Xe] 4f5 6s2',
        oxistates='3*',
        ionenergy=(5.582, 10.9),
        isotopes={145: Isotope(144.9127559, 1.0, 145)},
    ),
    Element(
        62, 'Sm', 'Samarium',
        group=3, period=6, block='f', series=9,
        mass=150.36, eleneg=1.17, eleaffin=0.5,
        covrad=1.62, atmrad=2.59, vdwrad=0.0,
        tboil=2067.0, tmelt=1347.0, density=7.54,
        eleconfig='[Xe] 4f6 6s2',
        oxistates='3*, 2',
        ionenergy=(5.6437, 11.07),
        isotopes={
            144: Isotope(143.9120065, 0.0307, 144),
            147: Isotope(146.9149044, 0.1499, 147),
            148: Isotope(147.9148292, 0.1124, 148),
            149: Isotope(148.9171921, 0.1382, 149),
            150: Isotope(149.9172829, 0.0738, 150),
            152: Isotope(151.9197397, 0.2675, 152),
            154: Isotope(153.9222169, 0.2275, 154),
        },
    ),
    Element(
        63, 'Eu', 'Europium',
        group=3, period=6, block='f', series=9,
        mass=151.964, eleneg=1.2, eleaffin=0.5,
        covrad=1.85, atmrad=2.56, vdwrad=0.0,
        tboil=1800.0, tmelt=1095.0, density=5.25,
        eleconfig='[Xe] 4f7 6s2',
        oxistates='3*, 2',
        ionenergy=(5.6704, 11.25),
        isotopes={
            151: Isotope(150.9198578, 0.4781, 151),
            153: Isotope(152.921238, 0.5219, 153),
        },
    ),
    Element(
        64, 'Gd', 'Gadolinium',
        group=3, period=6, block='f', series=9,
        mass=157.25, eleneg=1.2, eleaffin=0.5,
        covrad=1.61, atmrad=2.54, vdwrad=0.0,
        tboil=3545.0, tmelt=1585.0, density=7.89,
        eleconfig='[Xe] 4f7 5d 6s2',
        oxistates='3*',
        ionenergy=(6.1498, 12.1),
        isotopes={
            152: Isotope(151.9197995, 0.002, 152),
            154: Isotope(153.9208741, 0.0218, 154),
            155: Isotope(154.9226305, 0.148, 155),
            156: Isotope(155.9221312, 0.2047, 156),
            157: Isotope(156.9239686, 0.1565, 157),
            158: Isotope(157.9241123, 0.2484, 158),
            160: Isotope(159.9270624, 0.2186, 160),
        },
    ),
    Element(
        65, 'Tb', 'Terbium',
        group=3, period=6, block='f', series=9,
        mass=158.92535, eleneg=1.2, eleaffin=0.5,
        covrad=1.59, atmrad=2.51, vdwrad=0.0,
        tboil=3500.0, tmelt=1629.0, density=8.25,
        eleconfig='[Xe] 4f9 6s2',
        oxistates='4, 3*',
        ionenergy=(5.8638, 11.52),
        isotopes={159: Isotope(158.9253547, 1.0, 159)},
    ),
    Element(
        66, 'Dy', 'Dysprosium',
        group=3, period=6, block='f', series=9,
        mass=162.5, eleneg=1.22, eleaffin=0.5,
        covrad=1.59, atmrad=2.49, vdwrad=0.0,
        tboil=2840.0, tmelt=1685.0, density=8.56,
        eleconfig='[Xe] 4f10 6s2',
        oxistates='3*',
        ionenergy=(5.9389, 11.67),
        isotopes={
            156: Isotope(155.9242847, 0.00056, 156),
            158: Isotope(157.9244159, 0.00095, 158),
            160: Isotope(159.9252046, 0.02329, 160),
            161: Isotope(160.9269405, 0.18889, 161),
            162: Isotope(161.9268056, 0.25475, 162),
            163: Isotope(162.9287383, 0.24896, 163),
            164: Isotope(163.9291819, 0.2826, 164),
        },
    ),
    Element(
        67, 'Ho', 'Holmium',
        group=3, period=6, block='f', series=9,
        mass=164.93033, eleneg=1.23, eleaffin=0.5,
        covrad=1.58, atmrad=2.47, vdwrad=0.0,
        tboil=2968.0, tmelt=1747.0, density=8.78,
        eleconfig='[Xe] 4f11 6s2',
        oxistates='3*',
        ionenergy=(6.0215, 11.8),
        isotopes={165: Isotope(164.9303288, 1.0, 165)},
    ),
    Element(
        68, 'Er', 'Erbium',
        group=3, period=6, block='f', series=9,
        mass=167.259, eleneg=1.24, eleaffin=0.5,
        covrad=1.57, atmrad=2.45, vdwrad=0.0,
        tboil=3140.0, tmelt=1802.0, density=9.05,
        eleconfig='[Xe] 4f12 6s2',
        oxistates='3*',
        ionenergy=(6.1077, 11.93),
        isotopes={
            162: Isotope(161.9287884, 0.00139, 162),
            164: Isotope(163.9292088, 0.01601, 164),
            166: Isotope(165.9302995, 0.33503, 166),
            167: Isotope(166.9320546, 0.22869, 167),
            168: Isotope(167.9323767, 0.26978, 168),
            170: Isotope(169.9354702, 0.1491, 170),
        },
    ),
    Element(
        69, 'Tm', 'Thulium',
        group=3, period=6, block='f', series=9,
        mass=168.93422, eleneg=1.25, eleaffin=0.5,
        covrad=1.56, atmrad=2.42, vdwrad=0.0,
        tboil=2223.0, tmelt=1818.0, density=9.32,
        eleconfig='[Xe] 4f13 6s2',
        oxistates='3*, 2',
        ionenergy=(6.1843, 12.05, 23.71),
        isotopes={169: Isotope(168.9342179, 1.0, 169)},
    ),
    Element(
        70, 'Yb', 'Ytterbium',
        group=3, period=6, block='f', series=9,
        mass=173.054, eleneg=1.1, eleaffin=0.5,
        covrad=1.74, atmrad=2.4, vdwrad=0.0,
        tboil=1469.0, tmelt=1092.0, density=9.32,
        eleconfig='[Xe] 4f14 6s2',
        oxistates='3*, 2',
        ionenergy=(6.2542, 12.17, 25.2),
        isotopes={
            168: Isotope(167.9338896, 0.00123, 168),
            170: Isotope(169.9347664, 0.02982, 170),
            171: Isotope(170.9363302, 0.1409, 171),
            172: Isotope(171.9363859, 0.2168, 172),
            173: Isotope(172.9382151, 0.16103, 173),
            174: Isotope(173.9388664, 0.32026, 174),
            176: Isotope(175.9425764, 0.12996, 176),
        },
    ),
    Element(
        71, 'Lu', 'Lutetium',
        group=3, period=6, block='d', series=9,
        mass=174.9668, eleneg=1.27, eleaffin=0.5,
        covrad=1.56, atmrad=2.25, vdwrad=0.0,
        tboil=3668.0, tmelt=1936.0, density=9.84,
        eleconfig='[Xe] 4f14 5d 6s2',
        oxistates='3*',
        ionenergy=(5.4259, 13.9),
        isotopes={
            175: Isotope(174.9407752, 0.97401, 175),
            176: Isotope(175.9426897, 0.02599, 176),
        },
    ),
    Element(
        72, 'Hf', 'Hafnium',
        group=4, period=6, block='d', series=8,
        mass=178.49, eleneg=1.3, eleaffin=0.0,
        covrad=1.44, atmrad=2.16, vdwrad=0.0,
        tboil=4875.0, tmelt=2504.0, density=13.31,
        eleconfig='[Xe] 4f14 5d2 6s2',
        oxistates='4*',
        ionenergy=(6.8251, 14.9, 23.3, 33.3),
        isotopes={
            174: Isotope(173.9400461, 0.0016, 174),
            176: Isotope(175.9414076, 0.0526, 176),
            177: Isotope(176.9432277, 0.186, 177),
            178: Isotope(177.9437058, 0.2728, 178),
            179: Isotope(178.9458232, 0.1362, 179),
            180: Isotope(179.946557, 0.3508, 180),
        },
    ),
    Element(
        73, 'Ta', 'Tantalum',
        group=5, period=6, block='d', series=8,
        mass=180.94788, eleneg=1.5, eleaffin=0.322,
        covrad=1.34, atmrad=2.09, vdwrad=0.0,
        tboil=5730.0, tmelt=3293.0, density=16.68,
        eleconfig='[Xe] 4f14 5d3 6s2',
        oxistates='5*',
        ionenergy=(7.5496,),
        isotopes={
            180: Isotope(179.9474648, 0.0001201, 180),
            181: Isotope(180.9479958, 0.9998799, 181),
        },
    ),
    Element(
        74, 'W', 'Tungsten',
        group=6, period=6, block='d', series=8,
        mass=183.84, eleneg=2.36, eleaffin=0.815,
        covrad=1.3, atmrad=2.02, vdwrad=0.0,
        tboil=5825.0, tmelt=3695.0, density=19.26,
        eleconfig='[Xe] 4f14 5d4 6s2',
        oxistates='6*, 5, 4, 3, 2, 0',
        ionenergy=(7.864,),
        isotopes={
            180: Isotope(179.9467108, 0.0012, 180),
            182: Isotope(181.94820394, 0.265, 182),
            183: Isotope(182.95022275, 0.1431, 183),
            184: Isotope(183.95093092, 0.3064, 184),
            186: Isotope(185.9543628, 0.2843, 186),
        },
    ),
    Element(
        75, 'Re', 'Rhenium',
        group=7, period=6, block='d', series=8,
        mass=186.207, eleneg=1.9, eleaffin=0.15,
        covrad=1.28, atmrad=1.97, vdwrad=0.0,
        tboil=5870.0, tmelt=3455.0, density=21.03,
        eleconfig='[Xe] 4f14 5d5 6s2',
        oxistates='7, 6, 4, 2, -1',
        ionenergy=(7.8335,),
        isotopes={
            185: Isotope(184.9529545, 0.374, 185),
            187: Isotope(186.9557501, 0.626, 187),
        },
    ),
    Element(
        76, 'Os', 'Osmium',
        group=8, period=6, block='d', series=8,
        mass=190.23, eleneg=2.2, eleaffin=1.0778,
        covrad=1.26, atmrad=1.92, vdwrad=0.0,
        tboil=5300.0, tmelt=3300.0, density=22.61,
        eleconfig='[Xe] 4f14 5d6 6s2',
        oxistates='8, 6, 4*, 3, 2, 0, -2',
        ionenergy=(8.4382,),
        isotopes={
            184: Isotope(183.9524885, 0.0002, 184),
            186: Isotope(185.953835, 0.0159, 186),
            187: Isotope(186.9557474, 0.0196, 187),
            188: Isotope(187.9558352, 0.1324, 188),
            189: Isotope(188.9581442, 0.1615, 189),
            190: Isotope(189.9584437, 0.2626, 190),
            192: Isotope(191.961477, 0.4078, 192),
        },
    ),
    Element(
        77, 'Ir', 'Iridium',
        group=9, period=6, block='d', series=8,
        mass=192.217, eleneg=2.2, eleaffin=1.56436,
        covrad=1.27, atmrad=1.87, vdwrad=0.0,
        tboil=4700.0, tmelt=2720.0, density=22.65,
        eleconfig='[Xe] 4f14 5d7 6s2',
        oxistates='6, 4*, 3, 2, 1*, 0, -1',
        ionenergy=(8.967,),
        isotopes={
            191: Isotope(190.9605893, 0.373, 191),
            193: Isotope(192.9629216, 0.627, 193),
        },
    ),
    Element(
        78, 'Pt', 'Platinum',
        group=10, period=6, block='d', series=8,
        mass=195.084, eleneg=2.28, eleaffin=2.1251,
        covrad=1.3, atmrad=1.83, vdwrad=1.75,
        tboil=4100.0, tmelt=2042.1, density=21.45,
        eleconfig='[Xe] 4f14 5d9 6s',
        oxistates='4*, 2*, 0',
        ionenergy=(8.9588, 18.563),
        isotopes={
            190: Isotope(189.9599297, 0.00012, 190),
            192: Isotope(191.9610387, 0.00782, 192),
            194: Isotope(193.9626809, 0.3286, 194),
            195: Isotope(194.9647917, 0.3378, 195),
            196: Isotope(195.96495209, 0.2521, 196),
            198: Isotope(197.9678949, 0.07356, 198),
        },
    ),
    Element(
        79, 'Au', 'Gold',
        group=11, period=6, block='d', series=8,
        mass=196.966569, eleneg=2.54, eleaffin=2.30861,
        covrad=1.34, atmrad=1.79, vdwrad=1.66,
        tboil=3130.0, tmelt=1337.58, density=19.32,
        eleconfig='[Xe] 4f14 5d10 6s',
        oxistates='3*, 1',
        ionenergy=(9.2255, 20.5),
        isotopes={197: Isotope(196.96656879, 1.0, 197)},
    ),
    Element(
        80, 'Hg', 'Mercury',
        group=12, period=6, block='d', series=8,
        mass=200.592, eleneg=2.0, eleaffin=0.0,
        covrad=1.49, atmrad=1.76, vdwrad=0.0,
        tboil=629.88, tmelt=234.31, density=13.55,
        eleconfig='[Xe] 4f14 5d10 6s2',
        oxistates='2*, 1',
        ionenergy=(10.4375, 18.756, 34.2),
        isotopes={
            196: Isotope(195.9658326, 0.0015, 196),
            198: Isotope(197.9667686, 0.0997, 198),
            199: Isotope(198.96828064, 0.1687, 199),
            200: Isotope(199.96832659, 0.231, 200),
            201: Isotope(200.97030284, 0.1318, 201),
            202: Isotope(201.9706434, 0.2986, 202),
            204: Isotope(203.97349398, 0.0687, 204),
        },
    ),
    Element(
        81, 'Tl', 'Thallium',
        group=13, period=6, block='p', series=7,
        mass=204.3834, eleneg=2.04, eleaffin=0.377,
        covrad=1.48, atmrad=2.08, vdwrad=1.96,
        tboil=1746.0, tmelt=577.0, density=11.85,
        eleconfig='[Xe] 4f14 5d10 6s2 6p',
        oxistates='3, 1*',
        ionenergy=(6.1082, 20.428, 29.83),
        isotopes={
            203: Isotope(202.9723446, 0.2952, 203),
            205: Isotope(204.9744278, 0.7048, 205),
        },
    ),
    Element(
        82, 'Pb', 'Lead',
        group=14, period=6, block='p', series=7,
        mass=207.2, eleneg=2.33, eleaffin=0.364,
        covrad=1.47, atmrad=1.81, vdwrad=2.02,
        tboil=2023.0, tmelt=600.65, density=11.34,
        eleconfig='[Xe] 4f14 5d10 6s2 6p2',
        oxistates='4, 2*',
        ionenergy=(7.4167, 15.032, 31.937, 42.32, 68.8),
        isotopes={
            204: Isotope(203.973044, 0.014, 204),
            206: Isotope(205.9744657, 0.241, 206),
            207: Isotope(206.9758973, 0.221, 207),
            208: Isotope(207.9766525, 0.524, 208),
        },
    ),
    Element(
        83, 'Bi', 'Bismuth',
        group=15, period=6, block='p', series=7,
        mass=208.9804, eleneg=2.02, eleaffin=0.942363,
        covrad=1.46, atmrad=1.63, vdwrad=0.0,
        tboil=1837.0, tmelt=544.59, density=9.8,
        eleconfig='[Xe] 4f14 5d10 6s2 6p3',
        oxistates='5, 3*',
        ionenergy=(
            7.2855, 16.69, 25.56, 45.3, 56.0,
            88.3,
        ),
        isotopes={209: Isotope(208.9803991, 1.0, 209)},
    ),
    Element(
        84, 'Po', 'Polonium',
        group=16, period=6, block='p', series=5,
        mass=208.9824, eleneg=2.0, eleaffin=1.9,
        covrad=1.46, atmrad=1.53, vdwrad=0.0,
        tboil=0.0, tmelt=527.0, density=9.2,
        eleconfig='[Xe] 4f14 5d10 6s2 6p4',
        oxistates='6, 4*, 2',
        ionenergy=(8.414,),
        isotopes={209: Isotope(208.9824308, 1.0, 209)},
    ),
    Element(
        85, 'At', 'Astatine',
        group=17, period=6, block='p', series=6,
        mass=209.9871, eleneg=2.2, eleaffin=2.8,
        covrad=1.45, atmrad=1.43, vdwrad=0.0,
        tboil=610.0, tmelt=575.0, density=0.0,
        eleconfig='[Xe] 4f14 5d10 6s2 6p5',
        oxistates='7, 5, 3, 1, -1*',
        ionenergy=(),
        isotopes={210: Isotope(209.9871479, 1.0, 210)},
    ),
    Element(
        86, 'Rn', 'Radon',
        group=18, period=6, block='p', series=2,
        mass=222.0176, eleneg=0.0, eleaffin=0.0,
        covrad=0.0, atmrad=1.34, vdwrad=0.0,
        tboil=211.4, tmelt=202.0, density=9.23,
        eleconfig='[Xe] 4f14 5d10 6s2 6p6',
        oxistates='2*',
        ionenergy=(10.7485,),
        isotopes={222: Isotope(222.0175782, 1.0, 222)},
    ),
    Element(
        87, 'Fr', 'Francium',
        group=1, period=7, block='s', series=3,
        mass=223.0197, eleneg=0.7, eleaffin=0.0,
        covrad=0.0, atmrad=0.0, vdwrad=0.0,
        tboil=950.0, tmelt=300.0, density=0.0,
        eleconfig='[Rn] 7s',
        oxistates='1*',
        ionenergy=(4.0727,),
        isotopes={223: Isotope(223.019736, 1.0, 223)},
    ),
    Element(
        88, 'Ra', 'Radium',
        group=2, period=7, block='s', series=4,
        mass=226.0254, eleneg=0.9, eleaffin=0.0,
        covrad=0.0, atmrad=0.0, vdwrad=0.0,
        tboil=1413.0, tmelt=973.0, density=5.5,
        eleconfig='[Rn] 7s2',
        oxistates='2*',
        ionenergy=(5.2784, 10.147),
        isotopes={226: Isotope(226.0254103, 1.0, 226)},
    ),
    Element(
        89, 'Ac', 'Actinium',
        group=3, period=7, block='f', series=10,
        mass=227.0278, eleneg=1.1, eleaffin=0.0,
        covrad=0.0, atmrad=0.0, vdwrad=0.0,
        tboil=3470.0, tmelt=1324.0, density=10.07,
        eleconfig='[Rn] 6d 7s2',
        oxistates='3*',
        ionenergy=(5.17, 12.1),
        isotopes={227: Isotope(227.0277523, 1.0, 227)},
    ),
    Element(
        90, 'Th', 'Thorium',
        group=3, period=7, block='f', series=10,
        mass=232.0377, eleneg=1.3, eleaffin=0.0,
        covrad=1.65, atmrad=0.0, vdwrad=0.0,
        tboil=5060.0, tmelt=2028.0, density=11.72,
        eleconfig='[Rn] 6d2 7s2',
        oxistates='4*',
        ionenergy=(6.3067, 11.5, 20.0, 28.8),
        isotopes={232: Isotope(232.0380558, 1.0, 232)},
    ),
    Element(
        91, 'Pa', 'Protactinium',
        group=3, period=7, block='f', series=10,
        mass=231.03588, eleneg=1.5, eleaffin=0.0,
        covrad=0.0, atmrad=0.0, vdwrad=0.0,
        tboil=4300.0, tmelt=1845.0, density=15.37,
        eleconfig='[Rn] 5f2 6d 7s2',
        oxistates='5*, 4',
        ionenergy=(5.89,),
        isotopes={231: Isotope(231.0358842, 1.0, 231)},
    ),
    Element(
        92, 'U', 'Uranium',
        group=3, period=7, block='f', series=10,
        mass=238.02891, eleneg=1.38, eleaffin=0.0,
        covrad=1.42, atmrad=0.0, vdwrad=1.86,
        tboil=4407.0, tmelt=1408.0, density=18.97,
        eleconfig='[Rn] 5f3 6d 7s2',
        oxistates='6*, 5, 4, 3',
        ionenergy=(6.1941,),
        isotopes={
            234: Isotope(234.0409523, 5.4e-05, 234),
            235: Isotope(235.0439301, 0.007204, 235),
            238: Isotope(238.0507884, 0.992742, 238),
        },
    ),
    Element(
        93, 'Np', 'Neptunium',
        group=3, period=7, block='f', series=10,
        mass=237.0482, eleneg=1.36, eleaffin=0.0,
        covrad=0.0, atmrad=0.0, vdwrad=0.0,
        tboil=4175.0, tmelt=912.0, density=20.48,
        eleconfig='[Rn] 5f4 6d 7s2',
        oxistates='6, 5*, 4, 3',
        ionenergy=(6.2657,),
        isotopes={237: Isotope(237.0481736, 1.0, 237)},
    ),
    Element(
        94, 'Pu', 'Plutonium',
        group=3, period=7, block='f', series=10,
        mass=244.0642, eleneg=1.28, eleaffin=0.0,
        covrad=0.0, atmrad=0.0, vdwrad=0.0,
        tboil=3505.0, tmelt=913.0, density=19.74,
        eleconfig='[Rn] 5f6 7s2',
        oxistates='6, 5, 4*, 3',
        ionenergy=(6.026,),
        isotopes={244: Isotope(244.0642053, 1.0, 244)},
    ),
    Element(
        95, 'Am', 'Americium',
        group=3, period=7, block='f', series=10,
        mass=243.0614, eleneg=1.3, eleaffin=0.0,
        covrad=0.0, atmrad=0.0, vdwrad=0.0,
        tboil=2880.0, tmelt=1449.0, density=13.67,
        eleconfig='[Rn] 5f7 7s2',
        oxistates='6, 5, 4, 3*',
        ionenergy=(5.9738,),
        isotopes={243: Isotope(243.0613813, 1.0, 243)},
    ),
    Element(
        96, 'Cm', 'Curium',
        group=3, period=7, block='f', series=10,
        mass=247.0704, eleneg=1.3, eleaffin=0.0,
        covrad=0.0, atmrad=0.0, vdwrad=0.0,
        tboil=0.0, tmelt=1620.0, density=13.51,
        eleconfig='[Rn] 5f7 6d 7s2',
        oxistates='4, 3*',
        ionenergy=(5.9914,),
        isotopes={247: Isotope(247.0703541, 1.0, 247)},
    ),
    Element(
        97, 'Bk', 'Berkelium',
        group=3, period=7, block='f', series=10,
        mass=247.0703, eleneg=1.3, eleaffin=0.0,
        covrad=0.0, atmrad=0.0, vdwrad=0.0,
        tboil=0.0, tmelt=1258.0, density=13.25,
        eleconfig='[Rn] 5f9 7s2',
        oxistates='4, 3*',
        ionenergy=(6.1979,),
        isotopes={247: Isotope(247.0703073, 1.0, 247)},
    ),
    Element(
        98, 'Cf', 'Californium',
        group=3, period=7, block='f', series=10,
        mass=251.0796, eleneg=1.3, eleaffin=0.0,
        covrad=0.0, atmrad=0.0, vdwrad=0.0,
        tboil=0.0, tmelt=1172.0, density=15.1,
        eleconfig='[Rn] 5f10 7s2',
        oxistates='4, 3*',
        ionenergy=(6.2817,),
        isotopes={251: Isotope(251.0795886, 1.0, 251)},
    ),
    Element(
        99, 'Es', 'Einsteinium',
        group=3, period=7, block='f', series=10,
        mass=252.083, eleneg=1.3, eleaffin=0.0,
        covrad=0.0, atmrad=0.0, vdwrad=0.0,
        tboil=0.0, tmelt=1130.0, density=0.0,
        eleconfig='[Rn] 5f11 7s2',
        oxistates='3*',
        ionenergy=(6.42,),
        isotopes={252: Isotope(252.08298, 1.0, 252)},
    ),
    Element(
        100, 'Fm', 'Fermium',
        group=3, period=7, block='f', series=10,
        mass=257.0951, eleneg=1.3, eleaffin=0.0,
        covrad=0.0, atmrad=0.0, vdwrad=0.0,
        tboil=0.0, tmelt=1800.0, density=0.0,
        eleconfig='[Rn] 5f12 7s2',
        oxistates='3*',
        ionenergy=(6.5,),
        isotopes={257: Isotope(257.0951061, 1.0, 257)},
    ),
    Element(
        101, 'Md', 'Mendelevium',
        group=3, period=7, block='f', series=10,
        mass=258.0984, eleneg=1.3, eleaffin=0.0,
        covrad=0.0, atmrad=0.0, vdwrad=0.0,
        tboil=0.0, tmelt=1100.0, density=0.0,
        eleconfig='[Rn] 5f13 7s2',
        oxistates='3*',
        ionenergy=(6.58,),
        isotopes={258: Isotope(258.0984315, 1.0, 258)},
    ),
    Element(
        102, 'No', 'Nobelium',
        group=3, period=7, block='f', series=10,
        mass=259.101, eleneg=1.3, eleaffin=0.0,
        covrad=0.0, atmrad=0.0, vdwrad=0.0,
        tboil=0.0, tmelt=1100.0, density=0.0,
        eleconfig='[Rn] 5f14 7s2',
        oxistates='3, 2*',
        ionenergy=(6.65,),
        isotopes={259: Isotope(259.10103, 1.0, 259)},
    ),
    Element(
        103, 'Lr', 'Lawrencium',
        group=3, period=7, block='d', series=10,
        mass=262.1096, eleneg=1.3, eleaffin=0.0,
        covrad=0.0, atmrad=0.0, vdwrad=0.0,
        tboil=0.0, tmelt=1900.0, density=0.0,
        eleconfig='[Rn] 5f14 6d 7s2',
        oxistates='3*',
        ionenergy=(4.9,),
        isotopes={262: Isotope(262.10961, 1.0, 262)},
    ),
    Element(
        104, 'Rf', 'Rutherfordium',
        group=4, period=7, block='d', series=8,
        mass=267.1218, eleneg=0.0, eleaffin=0.0,
        covrad=0.0, atmrad=0.0, vdwrad=0.0,
        tboil=0.0, tmelt=0.0, density=0.0,
        eleconfig='[Rn] 5f14 6d2 7s2',
        oxistates='*',
        ionenergy=(6.0,),
        isotopes={267: Isotope(267.12179, 1.0, 267)},
    ),
    Element(
        105, 'Db', 'Dubnium',
        group=5, period=7, block='d', series=8,
        mass=268.1257, eleneg=0.0, eleaffin=0.0,
        covrad=0.0, atmrad=0.0, vdwrad=0.0,
        tboil=0.0, tmelt=0.0, density=0.0,
        eleconfig='[Rn] 5f14 6d3 7s2',
        oxistates='*',
        ionenergy=(),
        isotopes={268: Isotope(268.12567, 1.0, 268)},
    ),
    Element(
        106, 'Sg', 'Seaborgium',
        group=6, period=7, block='d', series=8,
        mass=271.1339, eleneg=0.0, eleaffin=0.0,
        covrad=0.0, atmrad=0.0, vdwrad=0.0,
        tboil=0.0, tmelt=0.0, density=0.0,
        eleconfig='[Rn] 5f14 6d4 7s2',
        oxistates='*',
        ionenergy=(),
        isotopes={271: Isotope(271.13393, 1.0, 271)},
    ),
    Element(
        107, 'Bh', 'Bohrium',
        group=7, period=7, block='d', series=8,
        mass=272.1383, eleneg=0.0, eleaffin=0.0,
        covrad=0.0, atmrad=0.0, vdwrad=0.0,
        tboil=0.0, tmelt=0.0, density=0.0,
        eleconfig='[Rn] 5f14 6d5 7s2',
        oxistates='*',
        ionenergy=(),
        isotopes={272: Isotope(272.13826, 1.0, 272)},
    ),
    Element(
        108, 'Hs', 'Hassium',
        group=8, period=7, block='d', series=8,
        mass=270.1343, eleneg=0.0, eleaffin=0.0,
        covrad=0.0, atmrad=0.0, vdwrad=0.0,
        tboil=0.0, tmelt=0.0, density=0.0,
        eleconfig='[Rn] 5f14 6d6 7s2',
        oxistates='*',
        ionenergy=(),
        isotopes={270: Isotope(270.13429, 1.0, 270)},
    ),
    Element(
        109, 'Mt', 'Meitnerium',
        group=9, period=7, block='d', series=8,
        mass=276.1516, eleneg=0.0, eleaffin=0.0,
        covrad=0.0, atmrad=0.0, vdwrad=0.0,
        tboil=0.0, tmelt=0.0, density=0.0,
        eleconfig='[Rn] 5f14 6d7 7s2',
        oxistates='*',
        ionenergy=(),
        isotopes={276: Isotope(276.15159, 1.0, 276)},
    ),
)
# fmt: on

ELEMENTARY_CHARGE = 1.602176634e-19

ELECTRON = Particle('Electron', 5.48579909065e-4, -ELEMENTARY_CHARGE)
PROTON = Particle('Proton', 1.007276466621, ELEMENTARY_CHARGE)
NEUTRON = Particle('Neutron', 1.00866491595, 0.0)
POSITRON = Particle('Positron', 5.48579909065e-4, ELEMENTARY_CHARGE)

PERIODS = {1: 'K', 2: 'L', 3: 'M', 4: 'N', 5: 'O', 6: 'P', 7: 'Q'}

BLOCKS = {'s': '', 'g': '', 'f': '', 'd': '', 'p': ''}

GROUPS = {
    1: ('IA', 'Alkali metals'),
    2: ('IIA', 'Alkaline earths'),
    3: ('IIIB', ''),
    4: ('IVB', ''),
    5: ('VB', ''),
    6: ('VIB', ''),
    7: ('VIIB', ''),
    8: ('VIIIB', ''),
    9: ('VIIIB', ''),
    10: ('VIIIB', ''),
    11: ('IB', 'Coinage metals'),
    12: ('IIB', ''),
    13: ('IIIA', 'Boron group'),
    14: ('IVA', 'Carbon group'),
    15: ('VA', 'Pnictogens'),
    16: ('VIA', 'Chalcogens'),
    17: ('VIIA', 'Halogens'),
    18: ('VIIIA', 'Noble gases'),
}

SERIES = {
    1: 'Nonmetals',
    2: 'Noble gases',
    3: 'Alkali metals',
    4: 'Alkaline earth metals',
    5: 'Metalloids',
    6: 'Halogens',
    7: 'Poor metals',
    8: 'Transition metals',
    9: 'Lanthanides',
    10: 'Actinides',
}


def sqlite_script():
    """Return SQL script to create sqlite database of elements.

    Examples
    --------
    >>> import sqlite3
    >>> con = sqlite3.connect(':memory:')
    >>> cur = con.executescript(sqlite_script())
    >>> con.commit()
    >>> for r in cur.execute("SELECT name FROM element WHERE number=6"):
    ...     str(r[0])
    'Carbon'
    >>> con.close()

    """
    sql = [
        """
        CREATE TABLE "period" (
            "number" TINYINT NOT NULL PRIMARY KEY,
            "label" CHAR NOT NULL UNIQUE,
            "description" VARCHAR(64)
        );
        CREATE TABLE "group" (
            "number" TINYINT NOT NULL PRIMARY KEY,
            "label" VARCHAR(8) NOT NULL,
            "description" VARCHAR(64)
        );
        CREATE TABLE "block" (
            "label" CHAR NOT NULL PRIMARY KEY,
            "description" VARCHAR(64)
        );
        CREATE TABLE "series" (
            "id" TINYINT NOT NULL PRIMARY KEY,
            "label"  VARCHAR(32) NOT NULL,
            "description" VARCHAR(256)
        );
        CREATE TABLE "element" (
            "number" TINYINT NOT NULL PRIMARY KEY,
            "symbol" VARCHAR(2) UNIQUE NOT NULL,
            "name" VARCHAR(16) UNIQUE NOT NULL,
            "period" TINYINT NOT NULL,
            --FOREIGN KEY("period") REFERENCES "period"(number),
            "group" TINYINT NOT NULL,
            --FOREIGN KEY("group") REFERENCES "group"(number),
            "block" CHAR NOT NULL,
            --FOREIGN KEY("block") REFERENCES "block"(label),
            "series" TINYINT NOT NULL,
            --FOREIGN KEY("series") REFERENCES "series"(id),
            "mass" REAL NOT NULL,
            "eleneg" REAL,
            "covrad" REAL,
            "atmrad" REAL,
            "vdwrad" REAL,
            "tboil" REAL,
            "tmelt" REAL,
            "density" REAL,
            "eleaffin" REAL,
            "eleconfig" VARCHAR(32),
            "oxistates" VARCHAR(32),
            "description" VARCHAR(2048)
        );
        CREATE TABLE "isotope" (
            "element" TINYINT NOT NULL,
            --FOREIGN KEY ("element") REFERENCES "element"("number"),
            "massnum" TINYINT NOT NULL,
            "mass" REAL NOT NULL,
            "abundance" REAL NOT NULL,
            PRIMARY KEY ("element", "massnum")
        );
        CREATE TABLE "eleconfig" (
            "element" TINYINT NOT NULL,
            --FOREIGN KEY ("element") REFERENCES "element"("number"),
            "shell" TINYINT NOT NULL,
            --FOREIGN KEY ("shell") REFERENCES "period"("number"),
            "subshell" CHAR NOT NULL,
            --FOREIGN KEY ("subshell") REFERENCES "block"("label"),
            "count" TINYINT,
            PRIMARY KEY ("element", "shell", "subshell")
        );
        CREATE TABLE "ionenergy" (
            "element" TINYINT NOT NULL,
            --FOREIGN KEY ("element") REFERENCES "element"("number"),
            "number" TINYINT NOT NULL,
            "energy" REAL NOT NULL,
            PRIMARY KEY ("element", "number")
        );
    """
    ]

    for key, label in PERIODS.items():
        sql.append(
            f"""INSERT INTO "period" VALUES ({key}, '{label}', NULL);"""
        )

    for key, (label, descr) in GROUPS.items():
        sql.append(
            f"""INSERT INTO "group" VALUES ({key}, '{label}', '{descr}');"""
        )

    for data in BLOCKS.items():
        sql.append(
            f"""INSERT INTO "block" VALUES ('{data[0]}', '{data[1]}');"""
        )

    for series in sorted(SERIES):
        sql.append(
            f"""INSERT INTO "series" VALUES (
            {series}, '{SERIES[series]}', ''\n);"""
        )

    for ele in ELEMENTS:
        descr = word_wrap(
            ele.description.replace("'", "\'\'").replace("\"", "\"\""),
            linelen=74,
            indent=0,
            joinstr='\n ',
        )
        sql.append(
            f"""INSERT INTO "element" VALUES (
            {ele.number}, '{ele.symbol}', '{ele.name}', {ele.period},
            {ele.group}, '{ele.block}', {ele.series}, {ele.mass:.10f},
            {ele.eleneg:.4f}, {ele.covrad:.4f}, {ele.atmrad:.4f},
            {ele.vdwrad:.4f}, {ele.tboil:4f}, {ele.tmelt:.4f},
            {ele.density:.4f}, {ele.eleaffin:.8f}, '{ele.eleconfig}',
            '{ele.oxistates}', '{descr}'\n);"""
        )

    for ele in ELEMENTS:
        for iso in ele.isotopes.values():
            sql.append(
                f"""INSERT INTO "isotope" VALUES (
                {ele.number}, {iso.massnumber},
                {iso.mass:.10f}, {iso.abundance:.8f}\n);"""
            )

    for ele in ELEMENTS:
        for (shell, subshell), count in ele.eleconfig_dict.items():
            sql.append(
                f"""INSERT INTO "eleconfig" VALUES (
                {ele.number}, {shell}, '{subshell}', {count}\n);"""
            )

    for ele in ELEMENTS:
        for i, ionenergy in enumerate(ele.ionenergy):
            sql.append(
                f"""INSERT INTO "ionenergy" VALUES (
                {ele.number}, {i + 1}, {ionenergy:.4f}\n);"""
            )

    sql = '\n'.join(sql)
    sql = sql.replace('                ', '            ')
    sql = sql.replace('        ', '')
    return sql


def word_wrap(text, linelen=79, indent=0, joinstr='\n'):
    """Return string, word wrapped at linelen."""
    if len(text) < linelen:
        return text
    result = []
    line = []
    llen = -indent
    for word in text.split():
        llen += len(word) + 1
        if llen < linelen:
            line.append(word)
        else:
            result.append(' '.join(line))
            line = [word]
            llen = len(word)
    if line:
        result.append(' '.join(line))
    return joinstr.join(result)


if __name__ == '__main__':
    import doctest

    print(f'ELEMENTS = {repr(ELEMENTS)}')
    # print(sqlite_script())
    doctest.testmod(verbose=False)
