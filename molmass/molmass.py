# molmass.py

# Copyright (c) 1990-2022, Christoph Gohlke
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

"""Molecular mass calculations.

Molmass is a Python library, console script, and web application to calculate
the molecular mass (average, nominal, and isotopic pure), the elemental
composition, and the mass distribution spectrum of a molecule given by its
chemical formula, relative element weights, or sequence.

Calculations are based on the isotopic composition of the elements. Mass
deficiency due to chemical bonding is not taken into account.

The library includes a database of physicochemical and descriptive properties
of the chemical elements.

:Author: `Christoph Gohlke <https://www.cgohlke.com>`_
:License: BSD 3-Clause
:Version: 2022.12.8
:DOI: 10.5281/zenodo.7135495

Quickstart
----------

Install the molmass package and all dependencies from the
Python Package Index::

    python -m pip install -U molmass[all]

Print the console script usage::

    python -m molmass --help

Run the web application::

    python -m molmass --web

The molmass library is documented via docstrings.

See `Examples`_ for using the programming interface.

Source code and support are available on
`GitHub <https://github.com/cgohlke/molmass>`_.

Requirements
------------

This release has been tested with the following requirements and dependencies
(other versions may work):

- `CPython 3.8.10, 3.9.13, 3.10.9, 3.11.1 <https://www.python.org>`_
- `Flask 2.2.2 <https://pypi.org/project/Flask/>`_ (optional)
- `Pandas 1.5.2 <https://pypi.org/project/pandas/>`_ (optional)
- `wxPython 4.2.0 <https://pypi.org/project/wxPython/>`_ (optional)

Revisions
---------

2022.12.8

- Fix split_charge formula with trailing ]] (#11).

2022.10.18

- Several breaking changes.
- Add experimental support for ion charges (#5).
- Change Element, Isotope, and Particle to dataclass (breaking).
- Change types of Spectrum and Composition (breaking).
- Add functions to export Spectrum and Composition as Pandas DataFrames.
- Replace lazyattr with functools.cached_property.
- Rename molmass_web to web (breaking).
- Change output of web application (breaking).
- Run web application using Flask if installed.
- Add options to specify URL of web application and not opening web browser.
- Convert to Google style docstrings.
- Add type hints.
- Remove support for Python 3.7.

2021.6.18

- Add Particle types to elements (#5).
- Fix molmass_web failure on WSL2 (#9).
- Fix elements_gui layout issue.
- Remove support for Python 3.6.

2020.6.10

- Fix elements_gui symbol size on WSL2.
- Support wxPython 4.1.

2020.1.1

- Update elements atomic weights and isotopic compositions from NIST.
- Move element descriptions into separate module.
- Remove support for Python 2.7 and 3.5.

2018.8.15

- Move modules into molmass package.

2018.5.29

- Add option to start web interface from console.
- Separate styles from content and use CSS flex layout in molmass_web.

2018.5.25

- Style and docstring fixes.
- Make from_fractions output deterministic.
- Accept Flask request.args in molmass_web.
- Style and template changes in molmass_web.

2016.2.25

- Fix some elements ionization energies.

2005.x.x

- Initial release.

Examples
--------

Calculate the molecular mass, elemental composition, and mass distribution of
a molecule from its chemical formula:

>>> from molmass import Formula
>>> f = Formula('C8H10N4O2')  # Caffeine
>>> f
Formula('C8H10N4O2')
>>> f.formula  # hill notation
'C8H10N4O2'
>>> f.empirical
'C4H5N2O'
>>> f.mass  # average mass
194.1909...
>>> f.nominal_mass  # == f.isotope.massnumber
194
>>> f.monoisotopic_mass  # == f.isotope.mass
194.0803...
>>> f.atoms
24
>>> f.charge
0
>>> f.composition().dataframe()
         Count  Relative mass  Fraction
Element...
C            8      96.085920  0.494801
H           10      10.079410  0.051905
N            4      56.026812  0.288514
O            2      31.998810  0.164780
>>> f.spectrum(min_intensity=0.01).dataframe()
             Relative mass  Fraction  Intensity %         m/z
Mass number...
194             194.080376  0.898828   100.000000  194.080376
195             195.082873  0.092625    10.305100  195.082873
196             196.084968  0.008022     0.892492  196.084968
197             197.087214  0.000500     0.055681  197.087214

Access physicochemical and descriptive properties of the chemical elements:

>>> from molmass import ELEMENTS, Element
>>> e = ELEMENTS['C']
>>> e
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
)
>>> e.number
6
>>> e.symbol
'C'
>>> e.name
'Carbon'
>>> e.description
'Carbon is a member of group 14 of the periodic table...'
>>> e.eleconfig
'[He] 2s2 2p2'
>>> e.eleconfig_dict
{(1, 's'): 2, (2, 's'): 2, (2, 'p'): 2}
>>> str(ELEMENTS[6])
'Carbon'
>>> len(ELEMENTS)
109
>>> sum(e.mass for e in ELEMENTS)
14693.181589001004
>>> for e in ELEMENTS:
...     e.validate()
...     e = eval(repr(e))

"""

from __future__ import annotations

__version__ = '2022.12.8'

__all__ = [
    'Composition',
    'CompositionItem',
    'Formula',
    'FormulaError',
    'Spectrum',
    'SpectrumEntry',
    'analyze',
    'format_charge',
    'from_elements',
    'from_fractions',
    'from_oligo',
    'from_peptide',
    'from_sequence',
    'from_string',
    'hill_sorted',
    'join_charge',
    'main',
    'split_charge',
    'test',
    'AMINOACIDS',
    'DEOXYNUCLEOTIDES',
    'GROUPS',
    'NUCLEOTIDES',
    'PREPROCESSORS',
]

import sys
import re
import math
import copy
from dataclasses import dataclass
from functools import reduce, cached_property

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from typing import Callable, Iterable, Iterator, Sequence
    import pandas

try:
    from .elements import ELEMENTS, ELECTRON, Isotope
except ImportError:
    from elements import ELEMENTS, ELECTRON, Isotope  # type: ignore


def analyze(
    formula: str,
    /,
    *,
    maxatoms: int = 512,
    min_intensity: float = 1e-4,
    debug: bool = False,
) -> str:
    """Return analysis of chemical formula as string.

    Parameters:
        formula: Chemical formula.
        maxatoms: Number of atoms below which to calculate spectrum.
        min_intensity: Minimum intensity to include in spectrum.

    Examples:
        >>> print(analyze('C8H10N4O2', min_intensity=0.01))
        Formula: C8H10N4O2
        Empirical formula: C4H5N2O
        <BLANKLINE>
        Nominal mass: 194
        Average mass: 194.19095
        Monoisotopic mass: 194.08038 (89.883%)
        Most abundant mass: 194.08038 (89.883%)
        Mean of distribution: 194.18604
        <BLANKLINE>
        Number of atoms: 24
        <BLANKLINE>
        Elemental Composition
        <BLANKLINE>
        Element  Count  Relative mass  Fraction %
        C            8       96.08592     49.4801
        H           10       10.07941      5.1905
        N            4       56.02681     28.8514
        O            2       31.99881     16.4780
        <BLANKLINE>
        Mass Distribution
        <BLANKLINE>
        A    Relative mass  Fraction %  Intensity %
        194      194.08038   89.882781   100.000000
        195      195.08287    9.262511    10.305100
        196      196.08497    0.802196     0.892492
        197      197.08721    0.050048     0.055681

    """
    result: list[str] = []
    try:
        f = Formula(formula)
        c = f.composition()
        if f.atoms < maxatoms:
            s = f.spectrum(min_intensity=min_intensity)
        else:
            s = None

        if len(str(f)) <= 50:
            result.append(f'Formula: {f}')
        if formula != f.formula:
            result.append(f'Hill notation: {f.formula}')
        if f.formula != f.empirical:
            result.append(f'Empirical formula: {f.empirical}')

        prec = precision_digits(f.mass, 9)
        result.append(f'\nNominal mass: {f.nominal_mass}')
        result.append(f'Average mass: {f.mass:.{prec}f}')
        result.append(
            'Monoisotopic mass: '
            f'{f.isotope.mass:.{prec}f} ({f.isotope.abundance * 100:.3f}%)'
        )
        if s is not None:
            result.append(
                'Most abundant mass: '
                f'{s.peak.mass:.{prec}f} ({s.peak.fraction * 100:.3f}%)'
            )
            result.append(f'Mean of distribution: {s.mean:.{prec}f}')

        result.append(f'\nNumber of atoms: {f.atoms}')

        if len(c) > 1:
            result.extend(('\nElemental Composition\n', str(c)))

        if s is not None:
            if len(s) > 1:
                result.append(f'\nMass Distribution\n\n{s}')

    except Exception as exc:
        if debug:
            raise
        result.append(f'Error: {exc}')

    return '\n'.join(result)


class Formula:
    """Chemical formula.

    Calculate various properties from formula string, such as hill notation,
    empirical formula, mass, elemental composition, and mass distribution.

    Parameters:
        formula:
            Chemical formula. May contain only symbols of chemical elements,
            groups, isotopes, parentheses, numbers, and trailing charge.
            The formula is not validated until it is parsed on demand
            during instance attribute access or method calls.
        groups:
            Mapping of chemical group name to formula in Hill notation.
            The default is :py:attr:`GROUPS`.

    Examples:
        Elements and counts:

        >>> Formula('H2O')
        Formula('H2O')

        Isotopes:

        >>> Formula('D2O')
        Formula('[2H]2O')

        >>> Formula('[30Si]3O2')
        Formula('[30Si]3O2')

        Ion charges:

        >>> Formula('[AsO4]3-')
        Formula('[AsO4]3-')

        Abbreviations of chemical groups:

        >>> Formula('EtOH')
        Formula('(C2H5)OH')

        Simple arithmetic:

        >>> Formula('(COOH)2')
        Formula('(COOH)2')

        >>> Formula('CuSO4.5H2O')
        Formula('CuSO4(H2O)5')

        Relative element weights:

        >>> Formula('O: 0.26, 30Si: 0.74')
        Formula('O2[30Si]3')

        Nucleotide sequences:

        >>> Formula('CGCGAATTCGCG')
        Formula('((C10H12N5O5P)2(C9H12N3O6P)4(C10H12N5O6P)4(C10H13N2O7P)2H2O)')

        >>> Formula('dsrna(CCUU)')
        Formula('((C10H12N5O6P)2(C9H12N3O7P)2...(C9H11N2O8P)2(H2O)2)')

        Peptide sequences:

        >>> Formula('MDRGEQGLLK')
        Formula('((C4H5NO3)(C5H7NO3)(C2H3NO)2...(C6H12N4O)H2O)')

        >>> Formula('peptide(CPK)')
        Formula('((C3H5NOS)(C6H12N2O)(C5H7NO)H2O)')

    """

    _charge: int
    _formula: str
    _formula_nocharge: str

    def __init__(
        self,
        formula: str = '',
        groups: dict[str, str] | None = None,
        /,
    ) -> None:
        self._formula = from_string(formula, groups)
        self._formula_nocharge, self._charge = split_charge(self._formula)

    @cached_property
    def _elements(self) -> dict[str, dict[int, int]]:
        """Number of atoms and isotopes by element.

        A ``dict[symbol: dict[massnumber: count]]``, where `massnumber` is
        either an isotope mass number or zero for a natural distribution
        of isotopes.

        Raises :py:class:`FormulaError` if formula is invalid.

        >>> Formula('H')._elements
        {'H': {0: 1}}
        >>> Formula('[2H]2O')._elements
        {'O': {0: 1}, 'H': {2: 2}}

        :meta public:

        """
        formula = self._formula_nocharge
        if not formula:
            raise FormulaError('empty formula', formula, 0)

        validchars = set('([{<123456789ABCDEFGHIKLMNOPRSTUVWXYZ')

        if formula[0] not in validchars:
            raise FormulaError(
                f'unexpected character {formula[0]!r}', formula, 0
            )

        validchars |= set(']})>0abcdefghiklmnoprstuy')

        elements: dict[str, dict[int, int]] = {}
        ele = ''  # parsed element
        num = 0  # number
        level = 0  # parenthesis level
        counts = [1]  # parenthesis level multiplication
        i = len(formula)
        while i:
            i -= 1
            char = formula[i]
            if char not in validchars:
                raise FormulaError(
                    f'unexpected character {char!r}', formula, i
                )
            if char in '([{<':
                level -= 1
                if level < 0 or num != 0:
                    raise FormulaError(
                        "missing closing parenthesis ')]}>'", formula, i
                    )
            elif char in ')]}>':
                if num == 0:
                    num = 1
                level += 1
                if level > len(counts) - 1:
                    counts.append(0)
                counts[level] = num * counts[level - 1]
                num = 0
            elif char.isdigit():
                j = i
                while i and formula[i - 1].isdigit():
                    i -= 1
                num = int(formula[i : j + 1])
                if num == 0:
                    raise FormulaError('count is zero', formula, i)
            elif char.islower():
                if not formula[i - 1].isupper():
                    raise FormulaError(
                        f'unexpected character {char!r}', formula, i
                    )
                ele = char
            elif char.isupper():
                ele = char + ele
                if num == 0:
                    num = 1
                if ele not in ELEMENTS:
                    raise FormulaError(f'unknown symbol {ele!r}', formula, i)
                iso_str = ''
                j = i
                while i and formula[i - 1].isdigit():
                    i -= 1
                    iso_str = formula[i] + iso_str
                if iso_str and i and formula[i - 1] not in '([{<':
                    i = j
                    iso_str = ''
                if iso_str:
                    iso = int(iso_str)
                    if iso not in ELEMENTS[ele].isotopes:
                        raise FormulaError(
                            f"unknown isotope '{iso}{ele}'", formula, i
                        )
                else:
                    iso = 0
                number = num * counts[level]
                if ele in elements:
                    item = elements[ele]
                    if iso in item:
                        item[iso] += number
                    else:
                        item[iso] = number
                else:
                    elements[ele] = {iso: number}
                ele = ''
                num = 0

        if num != 0:
            raise FormulaError('number preceding formula', formula, 0)

        if level != 0:
            raise FormulaError(
                "missing opening parenthesis '([{<'", formula, 0
            )

        if not elements:
            raise FormulaError('invalid formula', formula, 0)

        return elements

    @cached_property
    def formula(self) -> str:
        """Formula string in Hill notation.

        >>> Formula('BrC2H5').formula
        'C2H5Br'
        >>> Formula('[(CH3)3Si2]2NNa').formula
        'C6H18NNaSi4'

        """
        return from_elements(self._elements, charge=self._charge)

    @cached_property
    def empirical(self) -> str:
        """Empirical formula in Hill notation.

        The empirical formula has the simplest whole number ratio of atoms
        of each element present in the formula.

        >>> Formula('H2O').empirical
        'H2O'
        >>> Formula('C6H12O6').empirical
        'CH2O'
        >>> Formula('[SO4]2_4-').empirical
        '[O4S]2-'

        """
        return from_elements(
            self._elements, charge=self._charge, divisor=self.gcd
        )

    @cached_property
    def atoms(self) -> int:
        """Number of atoms.

        >>> Formula('CH3COOH').atoms
        8

        """
        return sum(sum(i.values()) for i in self._elements.values())

    @property
    def charge(self) -> int:
        """Charge number in units of elementary charge.

        >>> Formula('SO4_2-').charge
        -2

        """
        return self._charge

    @cached_property
    def gcd(self) -> int:
        """Greatest common divisor of element counts.

        >>> Formula('H2').gcd
        2
        >>> Formula('H2O').gcd
        1
        >>> Formula('C6H12O6').gcd
        6

        """
        values: list[int] = []
        for j in self._elements.values():
            values.extend(j.values())
        if abs(self._charge) > 0:
            values.append(abs(self._charge))
        return gcd(values)

    @cached_property
    def mass(self) -> float:
        """Average relative molecular mass.

        The sum of the relative atomic masses of all atoms and charges in
        the formula.
        Equals the molar mass in g/mol, i.e., the mass of one mole of
        substance.

        >>> Formula('H').mass
        1.007941
        >>> Formula('H+').mass
        1.007392...
        >>> Formula('SO4_2-').mass
        96.06351...
        >>> Formula('12C').mass
        12.0
        >>> Formula('C8H10N4O2').mass
        194.1909...
        >>> Formula('C48H32AgCuO12P2Ru4').mass
        1438.404...

        """
        result = 0.0
        for symbol in self._elements:
            ele = ELEMENTS[symbol]
            for massnumber, count in self._elements[symbol].items():
                if massnumber:
                    result += ele.isotopes[massnumber].mass * count
                else:
                    result += ele.mass * count
        return result - ELECTRON.mass * self._charge

    @property
    def monoisotopic_mass(self) -> float:
        """Mass of isotope composed of most abundant elemental isotopes.

        >>> Formula('C8H10N4O2').monoisotopic_mass
        194.08037...

        """
        return self.isotope.mass

    @property
    def nominal_mass(self) -> int:
        """Monoisotopic mass number.

        The number of protons and neutrons in the isotope composed of the most
        abundant elemental isotopes.

        >>> Formula('C8H10N4O2').nominal_mass
        194

        """
        return self.isotope.massnumber

    @property
    def mz(self) -> float:
        """Mass-to-charge ratio.

        >>> Formula('H').mz
        1.007941
        >>> Formula('H+').mz
        1.007392...
        >>> Formula('SO4_2-').mz
        48.03175...

        """
        if self._charge == 0:
            return self.mass
        return self.mass / abs(self._charge)

    @cached_property
    def isotope(self) -> Isotope:
        """Isotope composed of most abundant elemental isotopes.

        >>> Formula('C').isotope.mass
        12.0
        >>> Formula('13C').isotope.massnumber
        13
        >>> Formula('C48H32AgCuO12P2Ru4').isotope
        Isotope(mass=1439.588..., abundance=0.00205..., massnumber=1440...)

        """
        result = Isotope(-ELECTRON.mass * self._charge, 1.0, 0, self._charge)
        for symbol in self._elements:
            ele = ELEMENTS[symbol]
            for massnumber, count in self._elements[symbol].items():
                if massnumber != 0:
                    isotope = ele.isotopes[massnumber]
                else:
                    isotope = ele.isotopes[ele.nominalmass]
                result.mass += isotope.mass * count
                result.massnumber += isotope.massnumber * count
                result.abundance *= isotope.abundance**count
        return result

    def composition(self, isotopic: bool = True) -> Composition:
        """Return elemental composition.

        Parameters:
            isotopic:
                List isotopes separately as opposed to part of an element.

        Examples:
            >>> print(Formula('[12C]C').composition())
            Element  Count  Relative mass  Fraction %
            C            1      12.010740     50.0224
            12C          1      12.000000     49.9776

            >>> print(Formula('[12C]C').composition(False))
            Element  Count  Relative mass  Fraction %
            C            2      24.010740    100.0000

        """
        elements = self._elements
        result = []
        if isotopic:
            for symbol in hill_sorted(elements):
                ele = ELEMENTS[symbol]
                iso = elements[symbol]
                for massnumber in sorted(iso):
                    count = iso[massnumber]
                    if massnumber:
                        mass = ele.isotopes[massnumber].mass * count
                        symbol = f'{massnumber}{symbol}'
                    else:
                        mass = ele.mass * count
                    result.append((symbol, count, mass, mass / self.mass))
        else:
            for symbol in hill_sorted(elements):
                ele = ELEMENTS[symbol]
                mass = 0.0
                counter = 0
                for massnumber, count in elements[symbol].items():
                    counter += count
                    if massnumber:
                        mass += ele.isotopes[massnumber].mass * count
                    else:
                        mass += ele.mass * count
                result.append((symbol, counter, mass, mass / self.mass))
        if self._charge != 0:
            mass = -ELECTRON.mass * self._charge
            result.append(('e-', -self._charge, mass, mass / self.mass))
        return Composition(result)

    def spectrum(
        self,
        *,
        min_fraction: float = 1e-16,
        min_intensity: float | None = None,
    ) -> Spectrum:
        """Return low resolution mass spectrum.

        Calculated by combining the mass numbers of the elemental isotopes.

        Parameters:
            min_fraction: Minimum of fraction to return.
            min_intensity: Minimum intensity to return.

        Returns:
            Mapping of massnumber to mass, fraction, and intensity.

        Examples:
            >>> print(Formula('D').spectrum())
            A  Relative mass  Fraction %  Intensity %
            2      2.0141018  100.000000   100.000000

            >>> print(Formula('H').spectrum())
            A  Relative mass  Fraction %  Intensity %
            1      1.0078250   99.988500   100.000000
            2      2.0141018    0.011500     0.011501

            >>> print(Formula('H+').spectrum())
            A  Relative mass  Fraction %  Intensity %
            1      1.0072765   99.988500   100.000000
            2      2.0135532    0.011500     0.011501

            >>> print(Formula('D2').spectrum())
            A  Relative mass  Fraction %  Intensity %
            4      4.0282036  100.000000   100.000000

            >>> print(Formula('DH').spectrum())
            A  Relative mass  Fraction %  Intensity %
            3      3.0219268   99.988500   100.000000
            4      4.0282036    0.011500     0.011501

            >>> print(Formula('H2').spectrum())
            A  Relative mass  Fraction %  Intensity %
            2      2.0156501   99.977001   100.000000
            3      3.0219268    0.022997     0.023003
            4      4.0282036    0.000001     0.000001

            >>> print(Formula('DHO').spectrum())
            A   Relative mass  Fraction %  Intensity %
            19      19.016841   99.745528   100.000000
            20      20.021536    0.049468     0.049594
            21      21.021087    0.204981     0.205504
            22      22.027363    0.000024     0.000024

            >>> print(Formula('SO4_2-').spectrum(min_intensity=0.1))
            A   Relative mass  Fraction %  Intensity %  m/z
            96      95.952827   94.070057   100.000000  47.976413
            97      96.952996    0.886071     0.941927  48.476498
            98      97.949936    4.983307     5.297443  48.974968

        """
        spectrum: dict[int, list[float]] = {0: [0.0, 1.0, 0.0]}
        elements = self._elements

        for symbol in elements:
            ele = ELEMENTS[symbol]
            for massnumber, count in elements[symbol].items():
                if massnumber:
                    # specific isotope
                    iso = ele.isotopes[massnumber]
                    for key in reversed(sorted(spectrum)):
                        t = spectrum[key]
                        del spectrum[key]
                        if t[1] < min_fraction:
                            continue
                        f = t[1]
                        m = t[0] + iso.mass * count
                        k = key + iso.massnumber * count
                        if k in spectrum:
                            s = spectrum[k]
                            s[0] += (s[1] * s[0] + f * m) / (s[1] + f)
                            s[1] += f
                        else:
                            spectrum[k] = [m, f]
                else:
                    # mixture of isotopes
                    isotopes = ele.isotopes.values()
                    for _ in range(count):
                        for key in reversed(sorted(spectrum)):
                            t = spectrum[key]
                            del spectrum[key]
                            if t[1] < min_fraction:
                                continue
                            for iso in isotopes:
                                f = t[1] * iso.abundance
                                m = t[0] + iso.mass
                                k = key + iso.massnumber
                                if k in spectrum:
                                    s = spectrum[k]
                                    s[0] = (s[1] * s[0] + f * m) / (s[1] + f)
                                    s[1] += f
                                else:
                                    spectrum[k] = [m, f]

        # filter low intensities
        if min_intensity is not None:
            norm = 100 / max(v[1] for v in spectrum.values())
            for massnumber, value in spectrum.copy().items():
                if value[1] * norm < min_intensity:
                    del spectrum[massnumber]

        if self._charge != 0:
            # correct for electron mass
            mc = ELECTRON.mass * self._charge
            for value in spectrum.values():
                value[0] -= mc

        return Spectrum(spectrum, self._charge)

    def __mul__(self, number: int, /) -> Formula:
        """Return formula repeated number times.

        >>> Formula('HO-') * 2
        Formula('[(HO)2]2-')

        """
        if not isinstance(number, int) or number < 1:
            raise TypeError('can only multipy with positive number')
        return Formula(
            join_charge(
                f'({self._formula_nocharge}){number}', self._charge * number
            )
        )

    def __rmul__(self, number: int, /) -> Formula:
        """Return formula repeated number times.

        >>> 2 * Formula('H2O')
        Formula('(H2O)2')

        """
        return self.__mul__(number)

    def __add__(self, other: Formula, /) -> Formula:
        """Add this and other formula.

        >>> Formula('H2O') + Formula('HO-')
        Formula('[(H2O)(HO)]-')

        """
        if not isinstance(other, Formula):
            raise TypeError('can only add Formula instance')
        return Formula(
            join_charge(
                f'({self._formula_nocharge})({other._formula_nocharge})',
                self._charge + other.charge,
            )
        )

    def __sub__(self, other: Formula, /) -> Formula:
        """Subtract elements of other formula.

        >>> Formula('H2O') - Formula('O')
        Formula('H2')

        """
        if not isinstance(other, Formula):
            raise TypeError('can only subtract Formula instance')
        _elements = copy.deepcopy(self._elements)
        for symbol, isotopes in other._elements.items():
            if symbol not in _elements:
                raise ValueError(f'element {symbol} not in {self}')
            element = _elements[symbol]
            for massnumber, count in isotopes.items():
                if massnumber not in element:
                    raise ValueError(
                        f'element {massnumber}{symbol} not in {self}'
                    )
                element[massnumber] -= count
                if element[massnumber] < 0:
                    raise ValueError(
                        f'negative number of element {massnumber}{symbol}'
                    )
                if element[massnumber] == 0:
                    del element[massnumber]
                if not element:
                    del _elements[symbol]
        return Formula(
            from_elements(_elements, charge=self._charge - other.charge)
        )

    def __str__(self) -> str:
        return self._formula

    def __repr__(self) -> str:
        return f'Formula({self._formula!r})'


@dataclass
class CompositionItem:
    """Item of elemental composition."""

    symbol: str
    """Chemical symbol of element or electron."""

    count: int
    """Number of elements or electrons in composition."""

    mass: float
    """Relative mass."""

    fraction: float
    """Mass fraction."""

    def astuple(self) -> tuple[str, int, float, float]:
        """Return (symbol, count, mass, fraction)."""
        return self.symbol, self.count, self.mass, self.fraction


class Composition:
    """Elemental composition.

    The basic interface is an ordered ``dict[symbol, CompositionItem]]``.

    Parameters:
        items:
            Symbol, count, mass, and fraction for all elements in formula.
            Symbols must not be repeated.

    Examples:
        >>> c = Composition((('2H', 2, 4.028, 0.201), ('O', 1, 15.999, 0.799)))
        >>> c
        <Composition([('2H', 2, 4.028, 0.201), ...])>
        >>> print(c)
        Element  Count  Relative mass  Fraction %
        2H           2       4.028000     20.1000
        O            1      15.999000     79.9000
        >>> c['O']
        CompositionItem(symbol='O', count=1, mass=15.999, fraction=0.799)

    """

    _items: dict[str, CompositionItem]

    def __init__(self, items: Sequence[tuple[str, int, float, float]]) -> None:
        self._items = {
            e[0]: CompositionItem(e[0], e[1], e[2], e[3]) for e in items
        }

    @cached_property
    def total(self) -> CompositionItem:
        """Sums of counts (excluding electrons), masses, and fractions.

        >>> Formula('2HO').composition().total
        CompositionItem(symbol='Total', count=2, mass=18.0135..., fraction=1.0)

        """
        result = CompositionItem('Total', 0, 0.0, 0.0)
        for item in self._items.values():
            if item.symbol != 'e-':
                result.count += item.count
            result.mass += item.mass
            result.fraction += item.fraction
        return result

    def dataframe(self) -> pandas.DataFrame:
        """Return composition as pandas DataFrame.

        >>> Formula('2HO').composition().dataframe()
                 Count  Relative mass  Fraction
        Element...
        2H           1       2.014102  0.111811
        O            1      15.999405  0.888189

        """
        from pandas import DataFrame

        data = {
            'Element': [i.symbol for i in self._items.values()],
            'Count': [i.count for i in self._items.values()],
            'Relative mass': [i.mass for i in self._items.values()],
            'Fraction': [i.fraction for i in self._items.values()],
        }
        return DataFrame(data).set_index('Element')

    def asdict(self) -> dict[str, tuple[int, float, float]]:
        """Return dict[symbol, tuple[symbol, count, mass, fraction]].

        >>> Formula('2HO').composition().asdict()
        {'2H': (1, 2.0141..., 0.1118...), 'O': (1, 15.9994..., 0.8881...)}

        """
        return {key: value.astuple()[1:] for key, value in self._items.items()}

    def astuple(self) -> tuple[tuple[str, int, float, float], ...]:
        """Return tuple[tuple[symbol, count, mass, fraction], ...].

        >>> Formula('2HO').composition().astuple()
        (('2H', 1, 2.0141..., 0.1118...), ('O', 1, 15.9994..., 0.8881...))

        """
        return tuple(value.astuple() for value in self._items.values())

    def keys(self) -> Iterator[str]:
        return iter(self._items.keys())

    def values(self) -> Iterator[CompositionItem]:
        return iter(self._items.values())

    def items(self) -> Iterator[tuple[str, CompositionItem]]:
        return iter(self._items.items())

    def __getitem__(self, key: str) -> CompositionItem:
        return self._items[key]

    def __len__(self) -> int:
        return len(self._items)

    def __iter__(self) -> Iterator[str]:
        return iter(self._items)

    def __repr__(self) -> str:
        if len(self._items) > 1:
            return f'<Composition([{self.astuple()[0]}, ...])>'
        return f'Composition({list(self.astuple())})'

    def __str__(self) -> str:
        if len(self) == 0:
            return ''
        precision = precision_digits(self.total.mass, 9)
        result = ['Element  Count  Relative mass  Fraction %']
        for item in self._items.values():
            result.append(
                f'{item.symbol:<7s}'
                f'  {item.count:5}'
                f'  {item.mass:13.{precision}f}'
                f'  {item.fraction * 100:10.4f}'
            )
        return '\n'.join(result)


@dataclass
class SpectrumEntry:
    """Entry of mass distribution spectrum."""

    massnumber: int
    """Mass number, bin of spectrum."""

    mass: float
    """Relative mass."""

    fraction: float
    """Mass fraction or abundance."""

    intensity: float
    """Fraction normalized by maximum fraction in %."""

    mz: float
    """Mass-to-charge ratio."""

    def astuple(self) -> tuple[int, float, float, float, float]:
        """Return attributes (massnumber, mass, fraction, intensity, mz)."""
        return (
            self.massnumber,
            self.mass,
            self.fraction,
            self.intensity,
            self.mz,
        )


class Spectrum:
    """Mass distribution.

    The basic interface is a sorted ``dict[massnumber, SpectrumEntry]``.

    Parameters:
        spectrum: Mapping of mass number to relative mass and fraction.
        charge: Charge number.

    Examples:
        >>> s = Spectrum({1: [1.078, 0.9999], 2: [2.014, 0.0001]})
        >>> s
        <Spectrum({1: [1, 1.078, 0.9999, 100.0, 1.078], ...})>
        >>> print(s)
        A  Relative mass  Fraction %  Intensity %
        1      1.0780000   99.990000   100.000000
        2      2.0140000    0.010000     0.010001
        >>> s[2]
        SpectrumEntry(massnumber=2, mass=2.014, fraction=0.0001, ..., mz=2.014)

    """

    _spectrum: dict[int, SpectrumEntry]
    _charge: int

    def __init__(
        self, spectrum: dict[int, list[float]], /, charge: int = 0
    ) -> None:
        mz = 1.0 / (1 if charge == 0 else abs(charge))
        intensity = 100 / max(v[1] for v in spectrum.values())
        self._spectrum = {
            massnumber: SpectrumEntry(
                massnumber=massnumber,
                mass=entry[0],
                fraction=entry[1],
                intensity=entry[1] * intensity,
                mz=entry[0] * mz,
            )
            for massnumber, entry in sorted(spectrum.items())
        }
        self._charge = charge

    @cached_property
    def peak(self) -> SpectrumEntry:
        """Most abundant entry.

        >>> Formula('C8H10N4O2').spectrum().peak.mass
        194.0803...

        """
        return max(self._spectrum.values(), key=lambda x: x.fraction)

    @cached_property
    def mean(self) -> float:
        """Mean of all masses.

        >>> Formula('C8H10N4O2').spectrum().mean
        194.1909...

        """
        return sum((v.mass * v.fraction) for v in self._spectrum.values())

    @cached_property
    def range(self) -> tuple[int, int]:
        """Smallest and largest massnumbers.

        >>> Formula('C8H10N4O2').spectrum().range
        (194, 205)

        """
        return min(self._spectrum.keys()), max(self._spectrum.keys())

    def dataframe(self) -> pandas.DataFrame:
        """Return composition as pandas DataFrame.

        >>> Formula('SO4_2-').spectrum(min_intensity=0.1).dataframe()
                     Relative mass  Fraction  Intensity %        m/z
        Mass number...
        96               95.952827  0.940701   100.000000  47.976413
        97               96.952996  0.008861     0.941927  48.476498
        98               97.949936  0.049833     5.297443  48.974968

        """
        from pandas import DataFrame

        data = {
            'Mass number': [i.massnumber for i in self._spectrum.values()],
            'Relative mass': [i.mass for i in self._spectrum.values()],
            'Fraction': [i.fraction for i in self._spectrum.values()],
            'Intensity %': [i.intensity for i in self._spectrum.values()],
            'm/z': [i.mz for i in self._spectrum.values()],
        }
        return DataFrame(data).set_index('Mass number')

    def asdict(self) -> dict[int, tuple[float, float, float, float]]:
        """Return dict[massnumber, tuple[mass, fraction, intensity, mz]].

        >>> Formula('DH').spectrum().asdict()
        {3: (3.0..., 0.9..., 100.0, 3...), 4: (4.0..., 0.0..., 0.0..., 4...)}

        """
        return {
            key: value.astuple()[1:] for key, value in self._spectrum.items()
        }

    def astuple(self) -> tuple[tuple[int, float, float, float, float], ...]:
        """Return tuple[tuple[massnumber, mass, fraction, intensity, mz], ...].

        >>> Formula('DH').spectrum().astuple()
        ((3, 3.0..., 0.9..., 100.0, 3...), (4, 4.0..., 0.0..., 0.0..., 4...))

        """
        return tuple(value.astuple() for value in self._spectrum.values())

    def keys(self) -> Iterator[int]:
        return iter(self._spectrum.keys())

    def values(self) -> Iterator[SpectrumEntry]:
        return iter(self._spectrum.values())

    def items(self) -> Iterator[tuple[int, SpectrumEntry]]:
        return iter(self._spectrum.items())

    def __getitem__(self, key: int) -> SpectrumEntry:
        return self._spectrum[key]

    def __len__(self) -> int:
        return len(self._spectrum)

    def __iter__(self) -> Iterator[int]:
        return iter(self._spectrum)

    def __repr__(self) -> str:
        item = list(list(self._spectrum.values())[0].astuple())
        return '<Spectrum({' + f'{item[0]!r}: {item!r}, ...' + '})>'

    def __str__(self) -> str:
        if not self._spectrum:
            return ''
        a = len(str(self.range[-1]))
        mz = '  m/z' if abs(self._charge) > 1 else ''
        result = [f'A{" " * a} Relative mass  Fraction %  Intensity %{mz}']
        precision = precision_digits(self.peak.mass, 9)
        for entry in self.values():
            mz = f'  {entry.mz:<13.{precision}f}' if mz else ''
            result.append(
                f'{entry.massnumber:<{a}d}'
                f'  {entry.mass:>13.{precision}f}'
                f' {entry.fraction*100:11.6f}'
                f' {entry.intensity:12.6f}'
                f'{mz}'.strip()
            )
        return '\n'.join(result)


class FormulaError(Exception):
    """Error in chemical formula.

    Parameters:
        message: Error message.
        formula: Chemical formula in which error occurred.
        position: Position in formula where error occurred.

    Examples:
        >>> Formula('abc').formula
        Traceback (most recent call last):
         ...
        molmass.molmass.FormulaError: unexpected character 'a'
        abc
        ^

        >>> Formula('(H2O)2-H2O').formula
        Traceback (most recent call last):
         ...
        molmass.molmass.FormulaError: unexpected character '-'
        (H2O)2-H2O
        ......^

        >>> Formula('[11C]').formula
        Traceback (most recent call last):
         ...
        molmass.molmass.FormulaError: unknown isotope '11C'
        [11C]
        .^

    """

    message: str
    """Error message."""

    formula: str
    """Formula in which error occurred."""

    position: int
    """Position in formula where error occurred."""

    def __init__(
        self, message: str, formula: str = '', position: int = -1
    ) -> None:
        self.position = position
        self.message = message
        self.formula = formula
        Exception.__init__(self, message, formula, position)

    def __str__(self) -> str:
        if self.position < 0:
            return str(self.message)
        return f"{self.message}\n{self.formula}\n{'.' * self.position}^"


def from_string(formula: str, groups: dict[str, str] | None = None, /) -> str:
    r"""Return formula string from user input string.

    Parameters:
        formula:
            Chemical formula.
            Supports simple, non-nested arithmetic (\+ and \*), abbreviations
            of common chemical groups, peptides, oligos, and mass fractions.
        groups:
            Mapping of chemical group name to formula in Hill notation.
            The default is :py:attr:`GROUPS`.

    Returns:
        Chemical formula composed only of element and isotope symbols,
        parentheses, numbers, and trailing ion charges.

    Raises:
        FormulaError: Invalid formula.

    Examples:
        >>> from_string('Valohp')
        '(C5H8NO2)'
        >>> from_string('HLeu2OH')
        'H(C6H11NO)2OH'
        >>> from_string('D2O')
        '[2H]2O'
        >>> from_string('O: 0.26, 30Si: 0.74')
        'O2[30Si]3'
        >>> from_string('PhNH2.HCl')
        '(C6H5)NH2HCl'
        >>> from_string('CuSO4.5H2O')
        'CuSO4(H2O)5'
        >>> from_string('CuSO4+5*H2O')
        'CuSO4(H2O)5'
        >>> from_string('ssdna(AC)')
        '((C10H12N5O5P)(C9H12N3O6P)H2O)'
        >>> from_string('peptide(GG)')
        '((C2H3NO)2H2O)'

    """
    try:
        formula = formula.strip().replace(' ', '')
    except AttributeError as exc:
        raise ValueError('formula must be a string') from exc

    # abbreviations of common chemical groups
    if groups is None:
        groups = GROUPS
    if groups:
        for grp in reversed(sorted(groups)):
            formula = formula.replace(grp, f'({groups[grp]})')

    # list of mass fractions
    if ':' in formula and ',' in formula:
        fractions: dict[str, float] = {}
        try:
            for item in formula.split(','):
                items = item.split(':')
                fractions[items[0].strip()] = float(items[1].strip())
        except Exception as exc:
            raise FormulaError(
                'invalid list of mass fractions', formula
            ) from exc
        return from_fractions(fractions)

    # oligos and peptides
    if len(formula) > 1:
        fset = set(formula)
        if fset <= set('ATCG') and fset & set('ATG'):
            return from_oligo(formula, 'ssdna')
        if fset <= set('AUCG') and fset & set('AG'):
            return from_oligo(formula, 'ssrna')
        if fset <= set(AMINOACIDS.keys()) and fset & set('AEGMLQRT'):
            return from_peptide(formula)
        for dtype, func in PREPROCESSORS.items():
            for match in re.findall(dtype + r'\((.*?)\)', formula):
                formula = formula.replace(f'{dtype}({match})', func(match))

    # Deuterium
    formula = re.sub('(D)(?![a-z])', '[2H]', formula)

    # charge
    formula, charge = split_charge(formula)

    # arithmetic
    formula = formula.replace('.', '+')
    if '+' in formula:
        for match in re.findall(
            r'(?:\+|^)((\d+)\*?(.*?))(?:(?=\+)|$)', formula
        ):
            formula = formula.replace(match[0], f'({match[2]}){match[1]}')
        formula = formula.replace('+', '')
    if '-' in formula:
        FormulaError('subtraction not allowed', formula, formula.index('-'))

    if charge != 0:
        formula = f'[{formula}]{format_charge(charge)}'
    return formula


def from_elements(
    elements: dict[str, dict[int, int]],
    /,
    *,
    divisor: int = 1,
    charge: int = 0,
    fmt: tuple[str, str, str, str] | None = None,
) -> str:
    """Return formula string in Hill notation from elements dict.

    Parameters:
        elements:
            Number of atoms and isotopes by element.
            See :py:meth:`Formula._elements`.
        divisor:
            Number by which to divide element counts and charge.
        charge:
            Charge number.
        fmt:
            Format strings.

    Examples:
        >>> from_elements({'C': {0: 4, 12: 2}}, divisor=2)
        'C2[12C]'

        >>> from_elements(
        ...     {'C': {0: 4, 12: 2}}, divisor=2, fmt=('{}', '{}<sub>{}</sub>',
        ...     '<sup>{}</sup>{}', '<sup>{}</sup>{}<sub>{}</sub>'))
        'C<sub>2</sub><sup>12</sup>C'

    """
    if fmt is None:
        fmt = ('{}', '{}{}', '[{}{}]', '[{}{}]{}')
    formula = []
    for symbol in hill_sorted(elements):
        isotopes = elements[symbol]
        for massnumber in sorted(isotopes):
            count = isotopes[massnumber] // divisor
            if massnumber:
                if count == 1:
                    formula.append(fmt[2].format(massnumber, symbol))
                else:
                    formula.append(fmt[3].format(massnumber, symbol, count))
            else:
                if count == 1:
                    formula.append(fmt[0].format(symbol))
                else:
                    formula.append(fmt[1].format(symbol, count))
    return join_charge(''.join(formula), charge // divisor)


def from_fractions(
    fractions: dict[str, float],
    /,
    *,
    maxcount: int = 10,
    precision: float = 1e-4,
) -> str:
    """Return formula string from elemental mass fractions.

    Parameters:
        fractions: Mapping of element symbols to abundances.
        maxcount: Maximum count of single element in formula.
        precision: Threshold for finding smallest acceptable element count.

    Examples:
        >>> from_fractions({'H': 0.112, 'O': 0.888})
        'H2O'
        >>> from_fractions({'D': 0.2, 'O': 0.8})
        'O[2H]2'
        >>> from_fractions({'H': 8.97, 'C': 59.39, 'O': 31.64})
        'C5H9O2'
        >>> from_fractions({'O': 0.26, '30Si': 0.74})
        'O2[30Si]3'

    """
    if not fractions:
        return ''
    # divide normalized fractions by element/isotope mass
    numbers = {}
    sumfractions = sum(fractions.values())
    for symbol, fraction in fractions.items():
        if symbol == 'D':  # Deuterium
            symbol = '2H'
        if symbol[0].isupper():
            try:
                mass = ELEMENTS[symbol].mass
            except KeyError as exc:
                raise FormulaError(f'unknown element {symbol!r}') from exc
        else:
            if symbol.startswith('['):
                symbol = symbol[1:-1]
            i = 0
            while symbol[i].isdigit():
                i += 1
            massnum = int(symbol[:i])
            symbol = symbol[i:]
            try:
                mass = ELEMENTS[symbol].isotopes[massnum].mass
            except KeyError as exc:
                raise FormulaError(
                    f"unknown isotope '[{massnum}{symbol}]'"
                ) from exc
            symbol = f'[{massnum}{symbol}]'
        numbers[symbol] = fraction / (sumfractions * mass)

    # divide numbers by smallest number
    smallest = min(numbers.values())
    for symbol in numbers:
        numbers[symbol] /= smallest

    # find smallest factor that turns all numbers into integers
    precision *= len(numbers)
    best = 1e6
    factor = 1
    for i in range(1, maxcount):
        x = sum(abs((i * n) - round(i * n)) for n in numbers.values())
        if x < best:
            best = x
            factor = i
            if best < i * precision:
                break

    formula = []
    for symbol, number in sorted(numbers.items()):
        count = int(round(factor * number))
        if count > 1:
            formula.append(f'{symbol}{count}')
        else:
            formula.append(symbol)
    return ''.join(formula)


def from_sequence(sequence: str, groups: dict[str, str], /) -> str:
    """Return sequence as chemical formula.

    Parameters:
        sequence:
            DNA, RNA, or peptide sequence.
        groups:
            Mapping of sequence item to chemical formula in Hill notation.
            One of :py:attr:`DEOXYNUCLEOTIDES`, :py:attr:`NUCLEOTIDES`, or
            :py:attr:`AMINOACIDS`.

    Examples:
        >>> from_sequence('A', {'A': 'B'})
        '(B)'
        >>> from_sequence('AA', {'A': 'B'})
        '(B)2'

    """
    sequence, charge = split_charge(sequence)
    counts = {key: 0 for key in groups}
    for item in sequence:
        counts[item] += 1
    formula = []
    for key in sorted(groups):
        num = counts[key]
        if num == 1:
            formula.append(f'({groups[key]})')
        elif num:
            formula.append(f'({groups[key]}){num}')
    return join_charge(''.join(formula), charge)


def from_peptide(sequence: str, /) -> str:
    """Return chemical formula for polymer of unmodified amino acids.

    Parameters:
        sequence: Perptide sequence.

    Examples:
        >>> from_peptide('GG')
        '((C2H3NO)2H2O)'

        >>> f = Formula(from_peptide('GPAVL IMCFY WHKRQ NEDST_2+'))
        >>> print(f.formula, f.atoms, f.mass)
        [C107H159N29O30S2]2+ 327 2395.7168...

    """
    sequence = sequence.replace(' ', '')
    sequence, charge = split_charge(sequence)
    formula = f'({from_sequence(sequence, AMINOACIDS)}H2O)'
    return join_charge(formula, charge)


def from_oligo(sequence: str, dtype: str = 'ssdna', /) -> str:
    """Return chemical formula for polymer of unmodified (deoxy)nucleotides.

    Each strand includes a 5' monophosphate.

    Parameters:
        sequence:
            DNA or RNA sequence.
        dtype:
            Nucleic acid sequence type.
            One of 'ssdna', 'dsdna', 'ssrna', or 'dsrna'.

    Examples:
        >>> from_oligo('AC', 'ssdna')
        '((C10H12N5O5P)(C9H12N3O6P)H2O)'

        >>> from_oligo('AU', 'dsrna')
        '((C10H12N5O6P)2(C9H11N2O8P)2(H2O)2)'

        >>> f = Formula(from_oligo('ATC G', 'dsdna'))
        >>> print(f.formula, f.atoms, f.mass)
        C78H102N30O50P8 268 2507.60913...

        >>> f = Formula(from_oligo('AUC G_2+', 'ssrna'))
        >>> print(f.formula, f.atoms, f.mass)
        [C38H49N15O29P4]2+ 135 1303.7744...

    """
    sequence, charge = split_charge(sequence)
    sequence = sequence.replace(' ', '')
    dtype_str = dtype.lower()
    if 'rna' in dtype:
        items = NUCLEOTIDES
        complements = NUCLEOTIDE_COMPLEMENTS
    else:
        items = DEOXYNUCLEOTIDES
        complements = DEOXYNUCLEOTIDE_COMPLEMENTS

    if dtype_str.startswith('ds'):
        t = ''.join(complements[i] for i in sequence)
        formula = from_sequence(sequence + t, items)
        formula = f'({formula}(H2O)2)'
    else:
        formula = from_sequence(sequence, items)
        formula = f'({formula}H2O)'
    return join_charge(formula, charge)


def hill_sorted(symbols: Iterable[str], /) -> Iterator[str]:
    """Return iterator over element symbols in order of Hill notation.

    Parameters:
        symbols: Element symbols.

    Yields:
        Element symbols in order of Hill notation.

    Examples:
        >>> tuple(hill_sorted('HCO'))
        ('C', 'H', 'O')

    """
    symbols_set = set(symbols)
    if 'C' in symbols:
        symbols_set.remove('C')
        yield 'C'
        if 'H' in symbols_set:
            symbols_set.remove('H')
            yield 'H'
    yield from sorted(symbols_set)


def gcd(numbers: Iterable[int], /) -> int:
    """Return greatest common divisor of integer numbers.

    Using Euclid's algorithm.

    Parameters:
        numbers: Integer numbers.

    Examples:
        >>> gcd([4])
        4
        >>> gcd([3, 6])
        3
        >>> gcd([6, 7])
        1

    """

    def _gcd(a: int, b: int, /) -> int:
        # return greatest common divisor of two integer numbers
        while b:
            a, b = b, a % b
        return a

    return reduce(_gcd, set(numbers))


def precision_digits(f: float, width: int, /) -> int:
    """Return number of digits after decimal point to print f in width chars.

    Parameters:
        f: Floating point number to print in `width` characters.
        width: Maximum length of printed number.

    Examples:
        >>> precision_digits(-0.12345678, 5)
        2
        >>> precision_digits(1.23456789, 5)
        3
        >>> precision_digits(12.3456789, 5)
        2
        >>> precision_digits(12345.6789, 5)
        1

    """
    precision = math.log(abs(f), 10)
    if precision < 0:
        precision = 0
    precision = width - int(math.floor(precision))
    precision -= 3 if f < 0 else 2  # sign and decimal point
    if precision < 1:
        precision = 1
    return precision


def split_charge(formula: str, /) -> tuple[str, int]:
    """Return formula stripped from charge, and charge.

    Parameters:
        formula: Chemical formula with appended charge.

    Returns:
        Chemical formula without charge, and charge number.

    Examples:
        >>> split_charge('Formula')
        ('Formula', 0)
        >>> split_charge('Formula+')
        ('Formula', 1)
        >>> split_charge('Formula++')
        ('Formula', 2)
        >>> split_charge('[Formula]2+')
        ('Formula', 2)
        >>> split_charge('[[Formula]]2-')
        ('[Formula]', -2)
        >>> split_charge('[Formula]_2-')
        ('Formula', -2)
        >>> split_charge('Formula_2-')
        ('Formula', -2)
        >>> split_charge('Formula_+')
        ('Formula', 1)
        >>> split_charge('Formula+-')
        ('Formula', 0)
        >>> split_charge('Formul+a+')
        ('Formul+a', 1)

    """
    charge = 0
    m = re.search(r'([\]_])([0-9]{1,})([+-]{1,})$', formula)
    if m:
        m_delim, m_count, m_sign = m.groups()
        if m_count == '':
            charge = int(f'{m_sign}1')
        else:
            charge = int(f'{m_sign}{m_count}')
    else:
        m = re.search(r'([\]_]?)([+-]{1,})$', formula)
        charge = 0
        if m:
            m_delim, m_sign = m.groups()
            for char in m_sign:
                if char == '+':
                    charge += 1
                elif char == '-':
                    charge -= 1
    if m:
        if m_delim == '_':
            formula = formula.rsplit('_', 1)[0]
        elif m_delim == '':
            formula = formula.strip(m_sign)
        elif m_delim == ']':
            formula = formula.rsplit(']', 1)[0] + ']'
        if formula.startswith('[') and formula.endswith(']'):
            formula = formula[1:-1]
    return formula, charge


def join_charge(formula: str, charge: int, separator: str = '', /) -> str:
    """Return formula with charge appended.

    Parameters:
        formula: Chemical formula without charge.
        charge: Charge number.
        separator: Character separating formula from charge. One of '' or '_'.

    Returns:
        Formula with charge != 0 appended.
        If `separator` is empty, the formula is wrapped in square brackets.

    Examples:
        >>> join_charge('Formula', 0)
        'Formula'
        >>> join_charge('Formula', 1)
        '[Formula]+'
        >>> join_charge('Formula', 2)
        '[Formula]2+'
        >>> join_charge('Formula', -2, '_')
        'Formula_2-'

    """
    if charge != 0:
        if separator:
            formula = f'{formula}{separator}{format_charge(charge, "")}'
        else:
            formula = f'[{formula}]{format_charge(charge, "")}'
    return formula


def format_charge(charge: int, prefix='', /) -> str:
    """Return string representation of charge.

    Parameters:
        charge: Charge number.
        prefix: Character placed in front of charge in case ``abs(charge)>1``.

    Examples:
        >>> format_charge(-2, '_')
        '_2-'
        >>> format_charge(1, '_')
        '+'

    """
    if charge == 0:
        return '0'
    if prefix and abs(charge) <= 1:
        prefix = ''
    sign = '+' if charge > 0 else '-'
    count = f'{abs(charge)}' if abs(charge) > 1 else ''
    return f'{prefix}{count}{sign}'


# Common chemical groups
GROUPS: dict[str, str] = {
    'Abu': 'C4H7NO',
    'Acet': 'C2H3O',
    'Acm': 'C3H6NO',
    'Adao': 'C10H15O',
    'Ade': 'C5H5N5',
    'Aib': 'C4H7NO',
    'Ala': 'C3H5NO',
    'Arg': 'C6H12N4O',
    'Argp': 'C6H11N4O',
    'Asn': 'C4H6N2O2',
    'Asnp': 'C4H5N2O2',
    'Asp': 'C4H5NO3',
    'Aspp': 'C4H4NO3',
    'Asu': 'C8H13NO3',
    'Asup': 'C8H12NO3',
    'Boc': 'C5H9O2',
    'Bom': 'C8H9O',
    'Bpy': 'C10H8N2',  # Bipyridine
    'Brz': 'C8H6BrO2',
    'Bu': 'C4H9',
    'Bum': 'C5H11O',
    'Bz': 'C7H5O',
    'Bzl': 'C7H7',
    'Bzlo': 'C7H7O',
    'Cha': 'C9H15NO',
    'Chxo': 'C6H11O',
    'Cit': 'C6H11N3O2',
    'Citp': 'C6H10N3O2',
    'Clz': 'C8H6ClO2',
    'Cp': 'C5H5',
    'Cy': 'C6H11',
    'Cys': 'C3H5NOS',
    'Cysp': 'C3H4NOS',
    'Cyt': 'C4H5N3O',
    'Dde': 'C10H13O2',
    'Dnp': 'C6H3N2O4',
    'Et': 'C2H5',
    'Fmoc': 'C15H11O2',
    'For': 'CHO',
    'Gln': 'C5H8N2O2',
    'Glnp': 'C5H7N2O2',
    'Glp': 'C5H5NO2',
    'Glu': 'C5H7NO3',
    'Glup': 'C5H6NO3',
    'Gly': 'C2H3NO',
    'Gua': 'C5H5N5O',
    'Hci': 'C7H13N3O2',
    'Hcip': 'C7H12N3O2',
    'His': 'C6H7N3O',
    'Hisp': 'C6H6N3O',
    'Hser': 'C4H7NO2',
    'Hserp': 'C4H6NO2',
    'Hx': 'C6H11',
    'Hyp': 'C5H7NO2',
    'Hypp': 'C5H6NO2',
    'Ile': 'C6H11NO',
    'Ivdde': 'C14H21O2',
    'Leu': 'C6H11NO',
    'Lys': 'C6H12N2O',
    'Lysp': 'C6H11N2O',
    'Mbh': 'C15H15O2',
    'Me': 'CH3',
    'Mebzl': 'C8H9',
    'Meobzl': 'C8H9O',
    'Met': 'C5H9NOS',
    'Mmt': 'C20H17O',
    'Mtc': 'C14H19O3S',
    'Mtr': 'C10H13O3S',
    'Mts': 'C9H11O2S',
    'Mtt': 'C20H17',
    'Nle': 'C6H11NO',
    'Npys': 'C5H3N2O2S',
    'Nva': 'C5H9NO',
    'Odmab': 'C20H26NO3',
    'Orn': 'C5H10N2O',
    'Ornp': 'C5H9N2O',
    'Pbf': 'C13H17O3S',
    'Pen': 'C5H9NOS',
    'Penp': 'C5H8NOS',
    'Ph': 'C6H5',
    'Phe': 'C9H9NO',
    'Phepcl': 'C9H8ClNO',
    'Phg': 'C8H7NO',
    'Pmc': 'C14H19O3S',
    'Ppa': 'C8H7O2',
    'Pro': 'C5H7NO',
    'Prop': 'C3H7',
    'Py': 'C5H5N',
    'Pyr': 'C5H5NO2',
    'Sar': 'C3H5NO',
    'Ser': 'C3H5NO2',
    'Serp': 'C3H4NO2',
    'Sta': 'C8H15NO2',
    'Stap': 'C8H14NO2',
    'Tacm': 'C6H12NO',
    'Tbdms': 'C6H15Si',
    'Tbu': 'C4H9',
    'Tbuo': 'C4H9O',
    'Tbuthio': 'C4H9S',
    'Tfa': 'C2F3O',
    'Thi': 'C7H7NOS',
    'Thr': 'C4H7NO2',
    'Thrp': 'C4H6NO2',
    'Thy': 'C5H6N2O2',
    'Tips': 'C9H21Si',
    'Tms': 'C3H9Si',
    'Tos': 'C7H7O2S',
    'Trp': 'C11H10N2O',
    'Trpp': 'C11H9N2O',
    'Trt': 'C19H15',
    'Tyr': 'C9H9NO2',
    'Tyrp': 'C9H8NO2',
    'Ura': 'C4H4N2O2',
    'Val': 'C5H9NO',
    'Valoh': 'C5H9NO2',
    'Valohp': 'C5H8NO2',
    'Xan': 'C13H9O',
}

# Amino acids - H2O
AMINOACIDS: dict[str, str] = {
    'G': 'C2H3NO',  # Glycine, Gly
    'P': 'C5H7NO',  # Proline, Pro
    'A': 'C3H5NO',  # Alanine, Ala
    'V': 'C5H9NO',  # Valine, Val
    'L': 'C6H11NO',  # Leucine, Leu
    'I': 'C6H11NO',  # Isoleucine, Ile
    'M': 'C5H9NOS',  # Methionine, Met
    'C': 'C3H5NOS',  # Cysteine, Cys
    'F': 'C9H9NO',  # Phenylalanine, Phe
    'Y': 'C9H9NO2',  # Tyrosine, Tyr
    'W': 'C11H10N2O',  # Tryptophan, Trp
    'H': 'C6H7N3O',  # Histidine, His
    'K': 'C6H12N2O',  # Lysine, Lys
    'R': 'C6H12N4O',  # Arginine, Arg
    'Q': 'C5H8N2O2',  # Glutamine, Gln
    'N': 'C4H6N2O2',  # Asparagine, Asn
    'E': 'C5H7NO3',  # Glutamic Acid, Glu
    'D': 'C4H5NO3',  # Aspartic Acid, Asp
    'S': 'C3H5NO2',  # Serine, Ser
    'T': 'C4H7NO2',  # Threonine, Thr
}

# Deoxynucleoside monophosphates - H2O
DEOXYNUCLEOTIDES: dict[str, str] = {
    'A': 'C10H12N5O5P',
    'T': 'C10H13N2O7P',
    'C': 'C9H12N3O6P',
    'G': 'C10H12N5O6P',
}

# Nucleoside monophosphates - H2O
NUCLEOTIDES: dict[str, str] = {
    'A': 'C10H12N5O6P',
    'U': 'C9H11N2O8P',
    'C': 'C9H12N3O7P',
    'G': 'C10H12N5O7P',
}

DEOXYNUCLEOTIDE_COMPLEMENTS: dict[str, str] = {
    'A': 'T',
    'T': 'A',
    'C': 'G',
    'G': 'C',
}

NUCLEOTIDE_COMPLEMENTS: dict[str, str] = {
    'A': 'U',
    'U': 'A',
    'C': 'G',
    'G': 'C',
}

# Formula preprocessors
PREPROCESSORS: dict[str, Callable[[str], str]] = {
    'peptide': from_peptide,
    'ssdna': lambda x: from_oligo(x, 'ssdna'),
    'dsdna': lambda x: from_oligo(x, 'dsdna'),
    'ssrna': lambda x: from_oligo(x, 'ssrna'),
    'dsrna': lambda x: from_oligo(x, 'dsrna'),
}


def test(verbose: bool = False) -> None:
    """Test module.

    Parameters:
        verbose: Print status of testing.

    """
    formula: str
    empirical: str
    mass: float

    f = Formula('EtOH')
    assert f.formula == 'C2H6O'
    assert f.empirical == 'C2H6O'
    assert f.mass == 46.068531
    assert f.atoms == 9
    assert f.gcd == 1
    assert f.isotope.mass == 46.041864812949996
    assert f.isotope.massnumber == 46
    assert f.isotope.abundance == 0.9756627354527866
    assert f.composition().astuple() == (
        ('C', 2, 24.02148, 0.5214292593788155),
        ('H', 6, 6.047646, 0.13127499116479316),
        ('O', 1, 15.999405, 0.34729574945639136),
    )
    spectrum = f.spectrum(min_fraction=1e-9, min_intensity=1e-9)
    assert spectrum.asdict() == {
        46: (46.04186481295, 0.9756627354527866, 100.0, 46.04186481295),
        47: (
            47.04532293299781,
            0.022149945780234485,
            2.2702461593918635,
            47.04532293299781,
        ),
        48: (
            48.04629172951169,
            0.002142167346900965,
            0.21956023009393977,
            48.04629172951169,
        ),
        49: (
            49.04956894435918,
            4.4886295288353736e-05,
            0.004600595437061851,
            49.04956894435918,
        ),
        50: (
            50.05315154462678,
            2.6452603656753454e-07,
            2.711244643824313e-05,
            50.05315154462678,
        ),
        51: (
            51.05909655531064,
            1.6186163397759485e-10,
            1.6589916586542372e-08,
            51.05909655531064,
        ),
    }, spectrum.asdict()

    assert spectrum.mean == 46.06852122027406, spectrum.mean
    assert spectrum.peak.astuple() == (
        46,
        46.04186481295,
        0.9756627354527866,
        100.0,
        46.04186481295,
    ), spectrum.peak.astuple()

    # these formulas should pass
    for formula, empirical, mass in [
        (''.join(e.symbol for e in ELEMENTS), '', 14693.181589000998),
        ('12C', '[12C]', 12.0),
        ('12CC', 'C[12C]', 24.0107),
        ('SO4_2-', '[O4S]2-', 96.06351715981813),
        ('[CHNOP[13C]]2-', '[C[13C]HNOP]2-', 87.003003),
        ('[CHNOP[13C]]_2-', '[C[13C]HNOP]2-', 87.003003),
        ('Co(Bpy)(CO)4', '', 327.16),
        ('CH3CH2Cl', 'C2H5Cl', 64.5147),
        ('C1000H1000', 'CH', 13018.68),
        ('Ru2(CO)8', 'C4O4Ru', 426.2232),
        ('RuClH(CO)(PPh3)3', 'C55H46ClOP3Ru', 952.41392),
        ('PhSiMe3', 'C9H14Si', 150.29566),
        ('Ph(CO)C(CH3)3', 'C11H14O', 162.23156),
        ('HGlyGluTyrOH', 'C16H21N3O7', 367.35864),
        ('HCysTyrIleGlnAsnCysProLeuNH2', 'C41H65N11O11S2', 952.1519),
        ('CGCGAATTCGCG', 'C116H148N46O73P12', 3726.4),
        ('MDRGEQGLLK', 'C47H83N15O16S', 1146.3),
        ('CDCl3', 'C[2H]Cl3', 120.384),
        ('[13C]Cl4', '[13C]Cl4', 154.8153),
        ('C5(PhBu(EtCHBr)2)3', 'C53H78Br6', 1194.626),
        ('AgCuRu4(H)2[CO]12{PPh3}2', 'C48H32AgCuO12P2Ru4', 1438.4022),
        ('PhNH2.HCl', 'C6H8ClN', 129.5892),
        ('NH3.BF3', 'BF3H3N', 84.8357),
        ('CuSO4.5H2O', 'CuH10O9S', 249.68),
        (
            'HCysp(Trt)Tyrp(Tbu)IleGlnp(Trt)Asnp(Trt)ProLeuGlyNH2',
            'C101H113N11O11S',
            1689.13532,
        ),
    ]:
        if verbose:
            print(f'Trying Formula({formula!r}) ...', end='')
        try:
            f = Formula(formula)
            f.empirical
            f.mass
            f.spectrum
        except FormulaError as exc:
            print('Error:', exc)
            continue
        if empirical and f.empirical != empirical:
            print(
                f'Failure for {formula!r}:\n    Expected {empirical!r}, '
                f'got {f.empirical!r}:'
            )
            continue
        if mass and abs(f.mass - mass) > 0.1:
            print(
                f'Failure for {formula!r}:\n    Expected {mass}, got {f.mass}'
            )
            continue
        if verbose:
            print('ok')

    # these formulas are expected to fail
    for formula in [
        '',
        '()',
        '2',
        'a',
        '(a)',
        'C:H',
        'H:',
        'C[H',
        'H)2',
        'A',
        'Aa',
        '2lC',
        '1C',
        '[11C]',
        'H0',
        '()0',
        '(H)0C',
        'Ox: 0.26, 30Si: 0.74',
        'H^++',
        '[CHNOP[13C]]__2-',
    ]:
        if verbose:
            print(f'Trying Formula({formula!r}) ...', end='')
        try:
            empirical = Formula(formula).empirical
        except FormulaError as exc:
            if verbose:
                print('ok\nExpected error:', exc)
        else:
            print(
                f'Failure expected for {formula!r}, got '
                f'{Formula(formula).formula!r}'
            )

    from molmass import PROTON

    assert abs(Formula('[1H]+').mass - PROTON.mass) < 1e-7


def main(argv: list[str] | None = None, /) -> int:
    """Command line usage main function.

    Parameters:
        argv: Command line arguments.

    """
    if argv is None:
        argv = sys.argv

    import optparse

    def search_doc(r, d):
        return re.search(r, __doc__).group(1) if __doc__ else d

    parser = optparse.OptionParser(
        usage='usage: %prog [options] formula',
        description=search_doc('\n\n([^|]*?)\n\n', ''),
        version='%prog {}'.format(search_doc(':Version: (.*)', 'Unknown')),
        prog='molmass',
    )
    opt = parser.add_option
    opt(
        '--web',
        dest='web',
        action='store_true',
        default=False,
        help='start web application and open it in a web browser',
    )
    opt(
        '--url',
        dest='url',
        default=None,
        help='URL to run web application',
    )
    opt(
        '--nobrowser',
        dest='nobrowser',
        action='store_true',
        default=False,
        help='do not open web browser',
    )
    opt(
        '--test',
        dest='test',
        action='store_true',
        default=False,
        help='test the module',
    )
    opt(
        '--doctest',
        dest='doctest',
        action='store_true',
        default=False,
        help='run the internal tests',
    )
    opt('-v', '--verbose', dest='verbose', action='store_true', default=False)

    settings, formula_list = parser.parse_args()

    if settings.web:
        try:
            from .web import main as web_main
        except ImportError:
            from web import main as web_main  # type: ignore

        if formula_list:
            form = {'q': ''.join(formula_list)}
        else:
            form = None

        return web_main(
            url=settings.url, open_browser=not settings.nobrowser, form=form
        )
    if settings.doctest:
        import doctest

        try:
            import molmass.molmass as m
        except ImportError:
            m = None  # type: ignore
        doctest.testmod(m, optionflags=doctest.ELLIPSIS)
        return 0
    if settings.test:
        test(settings.verbose)
        return 0
    if formula_list:
        formula = ''.join(formula_list)
    else:
        parser.error('no formula specified')

    try:
        results = analyze(formula)
    except Exception as exc:
        print('\nError: \n  ', exc, sep='')
        raise exc
    else:
        print('\n', results, sep='')

    return 0


if '--doctest' in sys.argv:
    # enable doctest for @cached_property
    __test__ = {  # noqa
        'Formula._elements': Formula._elements,
        'Formula.formula': Formula.formula,
        'Formula.empirical': Formula.empirical,
        'Formula.atoms': Formula.atoms,
        'Formula.gcd': Formula.gcd,
        'Formula.mass': Formula.mass,
        'Formula.isotope': Formula.isotope,
        'Spectrum.range': Spectrum.range,
        'Spectrum.peak': Spectrum.peak,
        'Spectrum.mean': Spectrum.mean,
        'Composition.total': Composition.total,
    }

if __name__ == '__main__':
    sys.exit(main())
