# molmass.py

# Copyright (c) 1990-2026, Christoph Gohlke
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
:License: BSD-3-Clause
:Version: 2026.11.1
:DOI: `10.5281/zenodo.7135495 <https://doi.org/10.5281/zenodo.7135495>`_

Quickstart
----------

Install the molmass package and all dependencies from the
`Python Package Index <https://pypi.org/project/molmass/>`_::

    python -m pip install -U "molmass[all]"

See `Examples`_ and `Documentation <https://www.cgohlke.com/docs/molmass/>`_
for using the programming interface.

Print the console script usage::

    python -m molmass --help

Run the web application::

    python -m molmass --web

Try the `Molecular Mass Calculator <https://www.cgohlke.com/molmass/>`_ online.

Source code and support are available on
`GitHub <https://github.com/cgohlke/molmass>`_.

Requirements
------------

This revision was tested with the following requirements and dependencies
(other versions may work):

- `CPython <https://www.python.org>`_ 3.12.10, 3.13.16, 3.14.10, 3.15.2
- `Pandas <https://pypi.org/project/pandas/>`_ 3.1.0 (optional)
- `wxPython <https://pypi.org/project/wxPython/>`_ 4.3.2 (optional)

Revisions
---------

2026.11.1

- Reject empty brackets, for example C() or C1.5 (#24, breaking).
- Fix elements_gui rendering issues on Windows and macOS.

2026.10.31

- Update atomic weights and isotopic compositions (CIAAW 2024, AME2020).
- Update Dalton and particle masses to CODATA 2022.
- Reject charges inside brackets, for example [Na+] (#20, breaking).
- Reject negative and non-finite mass fractions (#22, breaking).
- Reject more invalid formulas and arguments (breaking).
- Make Element hashable and compare by identity (breaking).
- Compare Formula by Hill notation and support sum (breaking).
- Use abundance 1 for specified isotopes in Formula.isotope (breaking).
- Fix and improve from_fractions (breaking).
- Fix several formula parsing bugs.
- Fix Formula.isotope, arithmetic, spectrum, and concentration validation.
- Fix isotope, eleconfig, and ionenergy tables in elements.sqlite_script.
- Add density_unit argument to ppm conversion methods.
- Add Element.mass_uncertainty property.
- Fix and improve web and elements_gui applications.
- Improve docstrings and examples.
- Replace optparse with argparse.
- Run web application in http.server instead of Flask or CGI.

2026.8.15

- Add HTML application based on Pyodide.

2026.6.9

- Add convenience methods for concentration calculations.
- Use HTML Living Standard in web application.
- Drop support for Python 3.11, support Python 3.15.

2026.1.8

- …

Refer to the CHANGES file for older revisions.

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
    mass=12.01074, mass_uncertainty=0.001,
    eleneg=2.55, eleaffin=1.262118,
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
        13: Isotope(13.00335483534, 0.0107, 13),
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

"""

from __future__ import annotations

__version__ = '2026.11.1'

__all__ = [
    'AMINOACIDS',
    'DEOXYNUCLEOTIDES',
    'GROUPS',
    'NUCLEOTIDES',
    'PREPROCESSORS',
    'Composition',
    'CompositionItem',
    'DensityUnit',
    'Formula',
    'FormulaError',
    'MassUnit',
    'Spectrum',
    'SpectrumEntry',
    'VolumeUnit',
    '__version__',
    'analyze',
    'dilution',
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
    'mass_charge_ratio',
    'split_charge',
]

import copy
import math
import re
import sys
from dataclasses import dataclass
from functools import cached_property, reduce
from typing import TYPE_CHECKING, Literal

if TYPE_CHECKING:
    from collections.abc import Callable, Iterable, Iterator, Sequence

    import pandas

try:
    from .elements import ELECTRON, ELEMENTS, Isotope
except ImportError:
    from elements import ELECTRON, ELEMENTS, Isotope  # type: ignore[no-redef]


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
        formula:
            Chemical formula.
        maxatoms:
            Number of atoms below which to calculate spectrum.
        min_intensity:
            Minimum intensity percentage to include in spectrum.
        debug:
            Raise exceptions instead of printing error message.

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
            # avoid expensive spectrum calculation for many atoms
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
        if s is not None and len(s) > 0:
            result.append(
                'Most abundant mass: '
                f'{s.peak.mass:.{prec}f} ({s.peak.fraction * 100:.3f}%)'
            )
            result.append(f'Mean of distribution: {s.mean:.{prec}f}')

        result.append(f'\nNumber of atoms: {f.atoms}')

        if len(c) > 1:
            result.extend(('\nElemental Composition\n', str(c)))

        if s is not None and len(s) > 1:
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
            Some errors are raised during initialization, others only
            when the formula is parsed on demand during instance attribute
            access or method calls.
        groups:
            Mapping of chemical group name to formula in Hill notation.
            The default is :py:attr:`GROUPS`.
        parse_groups:
            Parse chemical group names. Enabled by default.
        parse_oligos:
            Parse pure peptide and nucleotide sequences. Enabled by default.
        parse_fractions:
            Parse list of mass fractions. Enabled by default.
        parse_arithmetic:
            Parse simple arithmetic operators. Enabled by default.
        allow_empty:
            Allow empty formulas.
            If False, raise :py:class:`FormulaError` for empty formulas.

    Examples:
        Elements and counts:

        >>> Formula('')
        Formula('')

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

        >>> Formula('O2-2')
        Formula('[O2]2-')

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
        /,
        groups: dict[str, str] | None = None,
        *,
        parse_groups: bool = True,
        parse_oligos: bool = True,
        parse_fractions: bool = True,
        parse_arithmetic: bool = True,
        allow_empty: bool = True,
    ) -> None:
        self._formula = from_string(
            formula,
            groups,
            parse_groups=parse_groups,
            parse_oligos=parse_oligos,
            parse_fractions=parse_fractions,
            parse_arithmetic=parse_arithmetic,
        )
        self._formula_nocharge, self._charge = split_charge(self._formula)
        if not self._formula_nocharge and self._charge != 0:
            msg = 'charge without formula'
            raise FormulaError(msg, formula, 0)
        if not allow_empty and not self._formula_nocharge:
            msg = 'empty formula'
            raise FormulaError(msg, formula, 0)

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
            return {}

        validchars = set('([{<123456789ABCDEFGHIKLMNOPRSTUVWXYZ')

        if formula[0] not in validchars:
            msg = f'unexpected character {formula[0]!r}'
            raise FormulaError(msg, formula, 0)

        validchars |= set(']})>0abcdefghiklmnoprstuy')

        elements: dict[str, dict[int, int]] = {}
        ele = ''  # parsed element
        num = 0  # number
        level = 0  # parenthesis level
        counts = [1]  # parenthesis level multiplication
        closing: list[str] = []  # stack of closing parentheses
        i = len(formula)
        while i:
            i -= 1
            char = formula[i]
            if char not in validchars:
                msg = f'unexpected character {char!r}'
                raise FormulaError(msg, formula, i)
            if char in '([{<':
                level -= 1
                if level < 0 or num != 0:
                    msg = "missing closing parenthesis ')]}>'"
                    raise FormulaError(msg, formula, i)
                if closing.pop() != ')]}>'['([{<'.index(char)]:
                    msg = 'mismatched parentheses'
                    raise FormulaError(msg, formula, i)
                if formula[i + 1] in ')]}>':
                    msg = 'empty brackets'
                    raise FormulaError(msg, formula, i)
            elif char in ')]}>':
                if num == 0:
                    num = 1
                level += 1
                if level > len(counts) - 1:
                    counts.append(0)
                counts[level] = num * counts[level - 1]
                num = 0
                closing.append(char)
            elif char.isdigit():
                j = i
                while i and formula[i - 1].isdigit():
                    i -= 1
                num = int(formula[i : j + 1])
                if num == 0:
                    msg = 'count is zero'
                    raise FormulaError(msg, formula, i)
            elif char.islower():
                if not formula[i - 1].isupper():
                    msg = f'unexpected character {char!r}'
                    raise FormulaError(msg, formula, i)
                ele = char
            elif char.isupper():
                ele = char + ele
                if num == 0:
                    num = 1
                if ele not in ELEMENTS:
                    msg = f'unknown symbol {ele!r}'
                    raise FormulaError(msg, formula, i)
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
                        msg = f"unknown isotope '{iso}{ele}'"
                        raise FormulaError(msg, formula, i)
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
            msg = 'number preceding formula'
            raise FormulaError(msg, formula, 0)

        if level != 0:
            msg = "missing opening parenthesis '([{<'"
            raise FormulaError(msg, formula, 0)

        if not elements:
            msg = 'invalid formula'
            raise FormulaError(msg, formula, 0)

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

    @property
    def expanded(self) -> str:
        """Formula string with expanded groups, charges, and sequences.

        The input formula parsed by :py:func:`from_string`.

        >>> Formula('WQ').expanded  # peptide sequence
        '((C5H8N2O2)(C11H10N2O)H2O)'

        """
        return self._formula

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
        Equals the molar mass in g/mol, that is, the mass of one mole of
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

        """
        result = 0.0
        for symbol, massnumber_counts in self._elements.items():
            ele = ELEMENTS[symbol]
            for massnumber, count in massnumber_counts.items():
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

        Using the average relative molecular mass.

        >>> Formula('H').mz
        1.007941
        >>> Formula('SO4_2-').mz
        48.03175...

        """
        return mass_charge_ratio(self.mass, self._charge)

    @property
    def isotope(self) -> Isotope:
        """Isotope composed of most abundant elemental isotopes.

        Specified isotopes are considered pure, with an abundance of 1.

        >>> Formula('C').isotope.mass
        12.0
        >>> Formula('13C').isotope
        Isotope(mass=13.0033..., abundance=1.0, massnumber=13...)
        >>> Formula('C48H32AgCuO12P2Ru4').isotope
        Isotope(mass=1439.588..., abundance=0.00205..., massnumber=1440...)

        """
        result = Isotope(-ELECTRON.mass * self._charge, 1.0, 0, self._charge)
        for symbol, massnumber_counts in self._elements.items():
            ele = ELEMENTS[symbol]
            for massnumber, count in massnumber_counts.items():
                if massnumber != 0:
                    isotope = ele.isotopes[massnumber]
                else:
                    isotope = ele.isotopes[ele.nominalmass]
                    result.abundance *= isotope.abundance**count
                result.mass += isotope.mass * count
                result.massnumber += isotope.massnumber * count
        return result

    def composition(self, *, isotopic: bool = True) -> Composition:
        """Return elemental composition.

        Parameters:
            isotopic:
                List isotopes separately as opposed to part of an element.

        Returns:
            Elemental composition.

        Examples:
            >>> print(Formula('[12C][13C]C').composition())
            Element  Count  Relative mass  Fraction %
            C            1      12.010740     32.4491
            12C          1      12.000000     32.4201
            13C          1      13.003355     35.1308

            >>> print(Formula('[12C][13C]C').composition(isotopic=False))
            Element  Count  Relative mass  Fraction %
            C            3      37.014095    100.0000

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
                        symbol_iso = f'{massnumber}{symbol}'
                    else:
                        mass = ele.mass * count
                        symbol_iso = symbol
                    result.append((symbol_iso, count, mass, mass / self.mass))
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
            min_fraction:
                Minimum fraction for pruning intermediate paths.
                Prune low-probability paths during combinatorial isotope
                expansion; does not filter final spectrum entries.
                Must be in range [0, 1].
            min_intensity:
                Minimum intensity percentage for final entries.
                Remove entries below this threshold, relative to the most
                abundant peak, from the returned spectrum.
                When set, ``Spectrum.mean`` will underestimate the true
                distribution mean.

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

            >>> print(Formula('SO4_2-').spectrum(min_intensity=0.1))
            A   Relative mass  Fraction %  Intensity %  m/z
            96      95.952827   94.070057   100.000000  47.976413
            97      96.952996    0.886071     0.941927  48.476498
            98      97.949936    4.983307     5.297443  48.974968

        """
        if not 0.0 <= min_fraction <= 1.0:
            msg = f'min_fraction must be in range [0, 1], not {min_fraction}'
            raise ValueError(msg)
        spectrum: dict[int, list[float]] = {0: [0.0, 1.0]}
        elements: dict[str, dict[int, int]] = self._elements

        if not elements:
            return Spectrum({}, self._charge)

        for symbol, massnumber_counts in elements.items():
            ele = ELEMENTS[symbol]
            for massnumber, count in massnumber_counts.items():
                if massnumber:
                    # specific isotope
                    iso = ele.isotopes[massnumber]
                    for key in sorted(spectrum, reverse=True):
                        t = spectrum[key]
                        del spectrum[key]
                        if t[1] < min_fraction:
                            continue
                        f = t[1]
                        m = t[0] + iso.mass * count
                        k = key + iso.massnumber * count
                        if k in spectrum:
                            s = spectrum[k]
                            s[0] = (s[1] * s[0] + f * m) / (s[1] + f)
                            s[1] += f
                        else:
                            spectrum[k] = [m, f]
                else:
                    # mixture of isotopes
                    isotopes = ele.isotopes.values()
                    for _ in range(count):
                        for key in sorted(spectrum, reverse=True):
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
        if min_intensity is not None and spectrum:
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
        if isinstance(number, bool) or not isinstance(number, int):
            return NotImplemented
        if number < 1:
            msg = 'can only multiply with positive integer'
            raise ValueError(msg)
        if not self._formula_nocharge:
            return Formula()
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
            return NotImplemented
        formula = ''.join(
            f'({f})'
            for f in (self._formula_nocharge, other._formula_nocharge)
            if f
        )
        return Formula(join_charge(formula, self._charge + other.charge))

    def __radd__(self, other: int, /) -> Formula:
        """Return formula for ``0 + formula``, enabling use with ``sum()``.

        >>> sum([Formula('H2O'), Formula('HO-')])
        Formula('[(H2O)(HO)]-')

        """
        if isinstance(other, bool) or other != 0:
            return NotImplemented
        return self

    def __sub__(self, other: Formula, /) -> Formula:
        """Subtract elements of other formula.

        >>> Formula('H2O') - Formula('O')
        Formula('H2')

        """
        if not isinstance(other, Formula):
            return NotImplemented
        _elements = copy.deepcopy(self._elements)
        for symbol, isotopes in other._elements.items():
            if symbol not in _elements:
                msg = f'element {symbol} not in {self.formula}'
                raise FormulaError(msg)
            element = _elements[symbol]
            for massnumber, count in isotopes.items():
                if massnumber not in element:
                    msg = f'element {massnumber}{symbol} not in {self.formula}'
                    raise FormulaError(msg)
                element[massnumber] -= count
                if element[massnumber] < 0:
                    msg = f'negative number of element {massnumber}{symbol}'
                    raise FormulaError(msg)
                if element[massnumber] == 0:
                    del element[massnumber]
                if not element:
                    del _elements[symbol]
        return Formula(
            from_elements(_elements, charge=self._charge - other.charge)
        )

    def moles(self, mass: float, /, unit: MassUnit = 'g') -> float:
        """Return number of moles from mass.

        Use to determine amount of substance in moles given mass.
        Calculate number of moles as mass / molar mass.

        Parameters:
            mass:
                Mass value.
            unit:
                Unit of mass
                {'g', 'mg', 'μg', 'ug', 'ng', 'kg', 'Da', 'kDa', 'MDa'}.

        Returns:
            Number of moles.

        Examples:
            >>> f = Formula('NaCl')
            >>> f.moles(5.844)
            0.09999...
            >>> f.moles(500, 'mg')
            0.008555...
            >>> f.moles(1.0, 'Da')
            2.841...e-26

        """
        if self.mass == 0.0:
            msg = 'cannot calculate moles for formula with zero mass'
            raise ValueError(msg)
        return float(mass) * unit_factor(MASS_UNIT, unit) / self.mass

    def mass_from_moles(self, moles: float, /, unit: MassUnit = 'g') -> float:
        """Return mass from number of moles.

        Use to determine mass of a substance given amount in moles.
        Calculate mass as number of moles * molar mass.

        Parameters:
            moles:
                Number of moles.
            unit:
                Desired output unit
                {'g', 'mg', 'μg', 'ug', 'ng', 'kg', 'Da', 'kDa', 'MDa'}.

        Returns:
            Mass in specified unit.

        Examples:
            >>> f = Formula('NaCl')
            >>> f.mass_from_moles(0.1)
            5.8442...
            >>> f.mass_from_moles(0.1, 'mg')
            5844.26...
            >>> f.mass_from_moles(1.0, 'Da')
            3.519...e+25

        """
        if self.mass == 0.0:
            msg = 'cannot calculate mass for formula with zero mass'
            raise ValueError(msg)
        return (float(moles) * self.mass) / unit_factor(MASS_UNIT, unit)

    def molarity(
        self,
        mass: float,
        volume: float,
        *,
        mass_unit: MassUnit = 'g',
        volume_unit: VolumeUnit = 'L',
    ) -> float:
        """Return molarity (mol/L) from mass and volume.

        Use to determine concentration of a solution.
        Calculate molarity as mass / molar mass / volume.

        Parameters:
            mass:
                Mass of solute.
            volume:
                Volume of solution.
            mass_unit:
                Unit of mass
                {'g', 'mg', 'μg', 'ug', 'ng', 'kg', 'Da', 'kDa', 'MDa'}.
            volume_unit:
                Unit of volume
                {'L', 'mL', 'μL', 'uL', 'nL', 'dL', 'm3', 'cm3', 'dm3'}.

        Returns:
            Molarity in mol/L.

        Examples:
            >>> f = Formula('NaCl')
            >>> f.molarity(5.844, 0.1)
            0.9999...
            >>> f.molarity(292.2, 50, mass_unit='mg', volume_unit='mL')
            0.09999...

        """
        volume = float(volume) * unit_factor(VOLUME_UNIT, volume_unit)
        if not volume > 0.0:
            msg = 'volume must be positive'
            raise ValueError(msg)
        if self.mass == 0.0:
            msg = 'cannot calculate molarity for formula with zero mass'
            raise ValueError(msg)
        mass = float(mass) * unit_factor(MASS_UNIT, mass_unit)
        return mass / self.mass / volume

    def mass_for_molarity(
        self,
        molarity: float,
        volume: float,
        *,
        mass_unit: MassUnit = 'g',
        volume_unit: VolumeUnit = 'L',
    ) -> float:
        """Return mass needed for desired molarity.

        Use to prepare solutions of desired concentration.
        Calculate mass as molarity * volume * molar mass.

        Parameters:
            molarity:
                Target concentration (mol/L).
            volume:
                Volume of solution.
            mass_unit:
                Desired output unit
                {'g', 'mg', 'μg', 'ug', 'ng', 'kg', 'Da', 'kDa', 'MDa'}.
            volume_unit:
                Unit of volume
                {'L', 'mL', 'μL', 'uL', 'nL', 'dL', 'm3', 'cm3', 'dm3'}.

        Returns:
            Mass in specified unit.

        Examples:
            >>> f = Formula('NaCl')
            >>> f.mass_for_molarity(0.15, 1.0)
            8.7664...
            >>> f.mass_for_molarity(0.15, 250, volume_unit='mL')
            2.1916...
            >>> f.mass_for_molarity(0.15, 0.25, mass_unit='mg')
            2191.6...

        """
        if self.mass == 0.0:
            msg = 'cannot calculate mass for formula with zero mass'
            raise ValueError(msg)
        volume = float(volume) * unit_factor(VOLUME_UNIT, volume_unit)
        if not volume >= 0.0:
            msg = 'volume must not be negative'
            raise ValueError(msg)
        return (float(molarity) * volume * self.mass) / unit_factor(
            MASS_UNIT, mass_unit
        )

    def ppm_to_molarity(
        self,
        ppm: float,
        /,
        density: float = 1.0,
        *,
        density_unit: DensityUnit = 'g/mL',
    ) -> float:
        """Return molarity from parts per million.

        Use to convert concentration in ppm (mg/kg) to molarity (mol/L).

        Parameters:
            ppm:
                Concentration in ppm (mg/kg).
            density:
                Solution density (default 1.0 g/mL for water).
            density_unit:
                Unit of density
                {'g/mL', 'g/cm3', 'kg/L', 'kg/dm3', 'g/L', 'kg/m3', 'mg/mL'}.

        Returns:
            Molarity in mol/L.

        Examples:
            >>> f = Formula('NaCl')
            >>> f.ppm_to_molarity(1000)
            0.01711...
            >>> f.ppm_to_molarity(1000, 1.2)
            0.02053...
            >>> f.ppm_to_molarity(1000, 1200, density_unit='kg/m3')
            0.02053...

        """
        density = float(density) * unit_factor(DENSITY_UNIT, density_unit)
        # ppm = mg solute / kg solution
        # convert to g/L: ppm (mg/kg) * density (kg/L) / 1000 (mg/g) = g/L
        # then divide by molar mass to get mol/L
        if self.mass == 0.0:
            msg = 'cannot calculate molarity for formula with zero mass'
            raise ValueError(msg)
        if not density > 0.0:
            msg = 'density must be positive'
            raise ValueError(msg)
        return ppm * density / self.mass / 1000.0

    def molarity_to_ppm(
        self,
        molarity: float,
        /,
        density: float = 1.0,
        *,
        density_unit: DensityUnit = 'g/mL',
    ) -> float:
        """Return parts per million from molarity.

        Use to convert molarity (mol/L) to concentration in ppm (mg/kg).

        Parameters:
            molarity:
                Concentration in mol/L.
            density:
                Solution density (default 1.0 g/mL for water).
            density_unit:
                Unit of density
                {'g/mL', 'g/cm3', 'kg/L', 'kg/dm3', 'g/L', 'kg/m3', 'mg/mL'}.

        Returns:
            Concentration in ppm (mg/kg).

        Examples:
            >>> f = Formula('NaCl')
            >>> f.molarity_to_ppm(0.1)
            5844.26...
            >>> f.molarity_to_ppm(0.1, 1.2)
            4870.22...
            >>> f.molarity_to_ppm(0.1, 1200, density_unit='g/L')
            4870.22...

        """
        density = float(density) * unit_factor(DENSITY_UNIT, density_unit)
        # molarity (mol/L) * molar mass (g/mol) = g/L
        # convert to ppm: g/L * 1000 (mg/g) / density (kg/L) = mg/kg
        if self.mass == 0.0:
            msg = 'cannot calculate ppm for formula with zero mass'
            raise ValueError(msg)
        if not density > 0.0:
            msg = 'density must be positive'
            raise ValueError(msg)
        return self.mass * molarity * 1000.0 / density

    def __eq__(self, other: object, /) -> bool:
        """Return whether formulas have same Hill notation.

        >>> Formula('HOH') == Formula('H2O')
        True

        """
        if not isinstance(other, Formula):
            return NotImplemented
        return self.formula == other.formula

    def __hash__(self) -> int:
        return hash(self.formula)

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
        >>> c = Formula('D2O').composition()
        >>> print(c)
        Element  Count  Relative mass  Fraction %
        2H           2       4.028204     20.1133
        O            1      15.999405     79.8867
        >>> c['O']
        CompositionItem(symbol='O', count=1, mass=15.999405, fraction=0.79...)

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
        CompositionItem(symbol='Total', count=2, mass=18.0135..., fraction=...)

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
        """Return iterator of symbols."""
        return iter(self._items.keys())

    def values(self) -> Iterator[CompositionItem]:
        """Return iterator of CompositionItems."""
        return iter(self._items.values())

    def items(self) -> Iterator[tuple[str, CompositionItem]]:
        """Return iterator of (symbol, CompositionItem) pairs."""
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
        result = [
            'Element  Count  Relative mass  Fraction %',
            *(
                f'{item.symbol:<7s}'
                f'  {item.count:5}'
                f'  {item.mass:13.{precision}f}'
                f'  {item.fraction * 100:10.4f}'
                for item in self._items.values()
            ),
        ]
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
        spectrum:
            Mapping of mass number to relative mass and fraction.
        charge:
            Charge number.

    Examples:
        >>> s = Formula('DH').spectrum()
        >>> print(s)
        A  Relative mass  Fraction %  Intensity %
        3      3.0219268   99.988500   100.000000
        4      4.0282036    0.011500     0.011501
        >>> s[4]
        SpectrumEntry(massnumber=4, mass=4.0282..., fraction=0.000115, ...)

    """

    _spectrum: dict[int, SpectrumEntry]
    _charge: int

    def __init__(
        self, spectrum: dict[int, list[float]], /, charge: int = 0
    ) -> None:
        self._charge = charge
        if not spectrum:
            self._spectrum = {}
            return
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

    @cached_property
    def peak(self) -> SpectrumEntry:
        """Most abundant entry.

        Raise ValueError if spectrum is empty.

        >>> Formula('C8H10N4O2').spectrum().peak.mass
        194.0803...

        """
        if not self._spectrum:
            msg = 'spectrum is empty'
            raise ValueError(msg)
        return max(self._spectrum.values(), key=lambda x: x.fraction)

    @cached_property
    def mean(self) -> float:
        """Weighted mean mass of spectrum entries.

        Calculated as the sum of each entry's mass weighted by its fraction.
        Accurate only for unfiltered spectra (where fractions sum to 1).
        When the spectrum was created with ``min_intensity`` filtering,
        low-abundance entries are excluded and the result is an underestimate
        of the true distribution mean.

        >>> Formula('C8H10N4O2').spectrum().mean
        194.1909...

        """
        return sum((v.mass * v.fraction) for v in self._spectrum.values())

    @cached_property
    def range(self) -> tuple[int, int]:
        """Smallest and largest massnumbers.

        Raise ValueError if spectrum is empty.

        >>> Formula('C8H10N4O2').spectrum().range
        (194, 205)

        """
        if not self._spectrum:
            msg = 'spectrum is empty'
            raise ValueError(msg)
        return min(self._spectrum.keys()), max(self._spectrum.keys())

    def dataframe(self) -> pandas.DataFrame:
        """Return spectrum as pandas DataFrame.

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

        >>> Formula('DH').spectrum().asdict()[4]
        (4.0282..., 0.000115, 0.0115..., 4.0282...)

        """
        return {
            key: value.astuple()[1:] for key, value in self._spectrum.items()
        }

    def astuple(self) -> tuple[tuple[int, float, float, float, float], ...]:
        """Return tuple[tuple[massnumber, mass, fraction, intensity, mz], ...].

        >>> Formula('DH').spectrum().astuple()[0]
        (3, 3.0219..., 0.999885, 100.0, 3.0219...)

        """
        return tuple(value.astuple() for value in self._spectrum.values())

    def keys(self) -> Iterator[int]:
        """Return iterator of massnumbers."""
        return iter(self._spectrum.keys())

    def values(self) -> Iterator[SpectrumEntry]:
        """Return iterator of SpectrumEntry objects."""
        return iter(self._spectrum.values())

    def items(self) -> Iterator[tuple[int, SpectrumEntry]]:
        """Return iterator of (massnumber, SpectrumEntry) pairs."""
        return iter(self._spectrum.items())

    def __getitem__(self, key: int) -> SpectrumEntry:
        return self._spectrum[key]

    def __len__(self) -> int:
        return len(self._spectrum)

    def __iter__(self) -> Iterator[int]:
        return iter(self._spectrum)

    def __repr__(self) -> str:
        charge = f', {self._charge}' if self._charge != 0 else ''
        if len(self._spectrum) == 0:
            return f'Spectrum({{}}{charge})'
        item = list(next(iter(self._spectrum.values())).astuple())
        return f'<Spectrum({{{item[0]!r}: {item!r}, ...}}{charge})>'

    def __str__(self) -> str:
        if not self._spectrum:
            return ''
        a = len(str(self.range[-1]))
        show_mz = abs(self._charge) > 1
        mz_header = '  m/z' if show_mz else ''
        result = [
            f'A{" " * a} Relative mass  Fraction %  Intensity %{mz_header}'
        ]
        precision = precision_digits(self.peak.mass, 9)
        for entry in self.values():
            mz = f'  {entry.mz:<13.{precision}f}' if show_mz else ''
            result.append(
                f'{entry.massnumber:<{a}d}'
                f'  {entry.mass:>13.{precision}f}'
                f' {entry.fraction * 100:11.6f}'
                f' {entry.intensity:12.6f}'
                f'{mz}'.strip()
            )
        return '\n'.join(result)


class FormulaError(ValueError):
    """Error in chemical formula.

    Parameters:
        message:
            Error message.
        formula:
            Chemical formula in which error occurred.
        position:
            Position in formula where error occurred.

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
        molmass.molmass.FormulaError: subtraction not allowed
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
        ValueError.__init__(self, message, formula, position)

    def __str__(self) -> str:
        if self.position < 0:
            return str(self.message)
        return f"{self.message}\n{self.formula}\n{'.' * self.position}^"


def from_string(
    formula: str,
    /,
    groups: dict[str, str] | None = None,
    *,
    parse_groups: bool = True,
    parse_oligos: bool = True,
    parse_fractions: bool = True,
    parse_arithmetic: bool = True,
) -> str:
    r"""Return formula string from user input string.

    Parameters:
        formula:
            Chemical formula.
            Supports simple, non-nested arithmetic (\+ and \*), abbreviations
            of common chemical groups, peptides, oligos, and mass fractions.
        groups:
            Mapping of chemical group name to formula in Hill notation.
            The default is :py:attr:`GROUPS`.
        parse_groups:
            Parse chemical group names.
        parse_oligos:
            Parse pure peptide and nucleotide sequences.
        parse_fractions:
            Parse list of mass fractions.
        parse_arithmetic:
            Parse simple arithmetic operators.

    Returns:
        Chemical formula composed only of element and isotope symbols,
        parentheses, numbers, and trailing ion charges.

    Raises:
        FormulaError: Invalid formula.
        TypeError: Non-string formula.

    Examples:
        >>> from_string('Valohp')
        '(C5H8NO2)'
        >>> from_string('Valohp', parse_groups=False)
        'Valohp'
        >>> from_string('HLeu2OH')
        'H(C6H11NO)2OH'
        >>> from_string('O: 0.26, 30Si: 0.74')
        'O2[30Si]3'
        >>> from_string('PhNH2.HCl')
        '(C6H5)NH2HCl'
        >>> from_string('PhNH2.HCl', parse_arithmetic=False)
        '(C6H5)NH2.HCl'
        >>> from_string('CuSO4.5H2O')
        'CuSO4(H2O)5'
        >>> from_string('CuSO4+5*H2O')
        'CuSO4(H2O)5'
        >>> from_string('5*H2O')
        '(H2O)5'
        >>> from_string('O2+12C')  # not [12C]O2
        'O2(C)12'
        >>> from_string('ssdna(AC)')
        '((C10H12N5O5P)(C9H12N3O6P)H2O)'
        >>> from_string('peptide(GG)')
        '((C2H3NO)2H2O)'
        >>> from_string('WQ')
        '((C5H8N2O2)(C11H10N2O)H2O)'
        >>> from_string('WQ', parse_oligos=False)
        'WQ'

    """
    try:
        formula = formula.strip().replace(' ', '')
    except AttributeError:
        msg = 'formula must be a string'
        raise TypeError(msg) from None

    if parse_groups:
        if groups is None:
            groups = GROUPS
        if groups:
            for grp in sorted(groups, reverse=True):
                formula = formula.replace(grp, f'({groups[grp]})')

    if parse_fractions and ':' in formula and ',' in formula:
        fractions: dict[str, float] = {}
        try:
            for item in formula.split(','):
                items = item.split(':')
                key = items[0].strip()
                value = float(items[1].strip())
                fractions[key] = fractions.get(key, 0.0) + value
        except Exception as exc:
            msg = 'invalid list of mass fractions'
            raise FormulaError(msg, formula) from exc
        return from_fractions(fractions)

    sequence, _ = split_charge(formula)
    if parse_oligos and len(sequence) > 1:
        fset = set(sequence)
        if fset <= set('ATCG') and fset & set('ATG'):
            return from_oligo(formula, 'ssdna')
        if fset <= set('AUCG') and fset & set('AG'):
            return from_oligo(formula, 'ssrna')
        if fset <= set(AMINOACIDS.keys()) and fset & set('AEGMLQRT'):
            return from_peptide(formula)
        for dtype, func in PREPROCESSORS.items():
            for match in re.findall(dtype + r'\((.*?)\)', formula):
                formula = formula.replace(f'{dtype}({match})', func(match))

    # replace Deuterium
    formula = _DEUTERIUM_PATTERN.sub('[2H]', formula)

    level = 0
    for i, char in enumerate(formula):
        if char in '([{<':
            level += 1
        elif char in ')]}>':
            level -= 1
        elif level > 0 and char in '+-':
            msg = (
                f'{char!r} inside brackets, '
                'charge must follow complete formula, for example [Na]+'
            )
            raise FormulaError(msg, formula, i)

    # charge
    formula, charge = split_charge(formula)

    if parse_arithmetic:
        formula = formula.replace('.', '+')
        if '+' in formula or '*' in formula:
            # TODO: document O2+12C is C12O2, not [12C]O2
            formula = _ARITHMETIC_PATTERN.sub(r'+(\3)\2', formula)
            # invalid lower case after operator
            m = _ARITHMETIC_INVALID_PATTERN.search(formula)
            if m:
                msg = 'invalid symbol'
                raise FormulaError(msg, formula, formula.index(m[0]) + 1)
            formula = formula.replace('+', '')
        if '-' in formula:
            msg = 'subtraction not allowed'
            raise FormulaError(msg, formula, formula.index('-'))

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
        ...     {'C': {0: 4, 12: 2}},
        ...     divisor=2,
        ...     fmt=(
        ...         '{}',
        ...         '{}<sub>{}</sub>',
        ...         '<sup>{}</sup>{}',
        ...         '<sup>{}</sup>{}<sub>{}</sub>',
        ...     ),
        ... )
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
            elif count == 1:
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
    tolerance: float | None = None,
) -> str:
    """Return empirical formula from elemental mass fractions.

    The result is the simplest whole number ratio of atoms in Hill order.

    Parameters:
        fractions:
            Mapping of element or isotope symbols to relative mass fractions.
            Fractions need not sum to 1. Fractions of symbols denoting the
            same element or isotope, for example 'D' and '[2H]', are added.
        maxcount:
            Largest multiplier tried to turn atom ratios into integers.
            Atom counts of individual elements may exceed `maxcount`.
        precision:
            Accept the smallest multiplier for which the summed deviation of
            atom counts from integers is below
            ``precision * number of symbols * multiplier``.
        tolerance:
            Maximum deviation of any atom count from an integer.
            By default, return the best fit regardless of deviation.

    Raises:
        FormulaError:
            Invalid symbol or fraction, or deviation exceeds `tolerance`.

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

    # number of atoms relative to each other
    numbers: dict[str, float] = {}
    for symbol, value in fractions.items():
        try:
            fraction = float(value)
        except (TypeError, ValueError):
            msg = f'invalid fraction {value!r}'
            raise FormulaError(msg) from None
        if not (math.isfinite(fraction) and fraction > 0.0):
            msg = 'fractions must be positive and finite'
            raise FormulaError(msg)

        bracketed = symbol.startswith('[')
        if bracketed != symbol.endswith(']'):
            msg = f'invalid symbol {symbol!r}'
            raise FormulaError(msg)
        name = symbol[1:-1] if bracketed else symbol
        if name == 'D':  # Deuterium
            name = '2H'
        i = len(name) - len(name.lstrip('0123456789'))
        element = name[i:]
        # ELEMENTS also maps element names, which are not valid symbols
        if element not in ELEMENTS or ELEMENTS[element].symbol != element:
            if element[:1].isupper():
                msg = f'unknown element {element!r}'
            else:
                msg = f'invalid symbol {symbol!r}'
            raise FormulaError(msg)
        ele = ELEMENTS[element]
        if i:
            massnum = int(name[:i])
            if massnum not in ele.isotopes:
                msg = f"unknown isotope '[{massnum}{element}]'"
                raise FormulaError(msg)
            key = f'[{massnum}{element}]'
            mass = ele.isotopes[massnum].mass
        else:
            key = element
            mass = ele.mass
        numbers[key] = numbers.get(key, 0.0) + fraction / mass

    smallest = min(numbers.values())
    numbers = {key: value / smallest for key, value in numbers.items()}

    # find smallest multiplier that turns all numbers into integers
    precision *= len(numbers)
    best = math.inf
    factor = 1
    for i in range(1, maxcount + 1):
        x = sum(abs((i * n) - round(i * n)) for n in numbers.values())
        if x < best:
            best = x
            factor = i
            if best < i * precision:
                break

    counts = {key: math.floor(factor * n + 0.5) for key, n in numbers.items()}
    if tolerance is not None:
        deviation = max(
            abs(factor * n - counts[key]) for key, n in numbers.items()
        )
        if deviation > tolerance:
            msg = f'atom counts deviate from integers by {deviation:.3g}'
            raise FormulaError(msg)

    return ''.join(
        key if counts[key] == 1 else f'{key}{counts[key]}'
        for key in hill_sorted(counts)
    )


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
        >>> from_sequence('GAG', AMINOACIDS)
        '(C3H5NO)(C2H3NO)2'

    """
    sequence, charge = split_charge(sequence)
    counts = dict.fromkeys(groups, 0)
    for i, item in enumerate(sequence):
        if item not in counts:
            msg = f'unknown sequence item {item!r}'
            raise FormulaError(msg, sequence, i)
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
        sequence: Peptide sequence.

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
    dtype = dtype.lower()
    if dtype not in {'ssdna', 'dsdna', 'ssrna', 'dsrna'}:
        msg = f'unknown nucleotide sequence type {dtype!r}'
        raise ValueError(msg)

    if 'rna' in dtype:
        items = NUCLEOTIDES
        complements = NUCLEOTIDE_COMPLEMENTS
    else:
        items = DEOXYNUCLEOTIDES
        complements = DEOXYNUCLEOTIDE_COMPLEMENTS

    if dtype.startswith('ds'):
        # unknown items are passed through and rejected by from_sequence
        t = ''.join(complements.get(i, i) for i in sequence)
        formula = from_sequence(sequence + t, items)
        formula = f'({formula}(H2O)2)'
    else:
        formula = from_sequence(sequence, items)
        formula = f'({formula}H2O)'
    return join_charge(formula, charge)


def mass_charge_ratio(mass: float, charge: int, /) -> float:
    """Return mass-to-charge ratio.

    Parameters:
        mass:
            Mass to divide by charge.
        charge:
            Charge number in units of elementary charge.

    Examples:
        >>> mass_charge_ratio(1.0, -2)
        0.5
        >>> mass_charge_ratio(1.0, 0)
        1.0
        >>> f = Formula('SO4_2-')
        >>> mass_charge_ratio(f.monoisotopic_mass, f.charge)
        47.976...

    """
    return mass if charge == 0 else mass / abs(charge)


def dilution(
    molarity: float,
    volume: float,
    *,
    final_molarity: float | None = None,
    final_volume: float | None = None,
) -> tuple[float, float]:
    """Return dilution parameters using M1V1 = M2V2.

    Use to calculate final concentration or volume after dilution.
    Calculate using the formula M1V1 = M2V2.

    Volumes must be in the same unit; the unit does not affect the
    calculation.

    Must provide either final_molarity or final_volume.

    Parameters:
        molarity:
            Initial concentration (M).
        volume:
            Initial volume (any unit).
        final_molarity:
            Final concentration.
        final_volume:
            Final volume (same unit).

    Returns:
        Tuple of (final_molarity, final_volume).

    Examples:
        >>> dilution(1.0, 10, final_volume=100)
        (0.1, 100)
        >>> dilution(12.0, 10, final_molarity=1.0)
        (1.0, 120.0)

    """
    if final_volume is None and final_molarity is None:
        msg = 'must provide either final_volume or final_molarity'
        raise ValueError(msg)
    if final_volume is not None and final_molarity is not None:
        msg = 'provide only one: final_volume or final_molarity'
        raise ValueError(msg)

    if final_volume is not None:
        if final_volume == 0.0:
            msg = 'final_volume must be non-zero'
            raise ValueError(msg)
        return (molarity * volume) / final_volume, final_volume

    assert final_molarity is not None
    if final_molarity == 0.0:
        msg = 'final_molarity must be non-zero'
        raise ValueError(msg)
    return final_molarity, (molarity * volume) / final_molarity


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
    if 'C' in symbols_set:
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
        >>> gcd([4, 6])
        2
        >>> gcd([])
        1

    """

    def _gcd(a: int, b: int, /) -> int:
        # return greatest common divisor of two integer numbers
        while b:
            a, b = b, a % b
        return a

    try:
        return reduce(_gcd, set(numbers))
    except TypeError:
        return 1


def precision_digits(f: float, width: int, /) -> int:
    """Return number of digits after decimal point to print f in width chars.

    Parameters:
        f:
            Floating point number to print in `width` characters.
        width:
            Maximum length of printed number.

    Examples:
        >>> precision_digits(1.23456789, 5)
        3
        >>> precision_digits(12345.6789, 5)
        1

    """
    try:
        precision = math.log10(abs(f))
    except ValueError:
        # math domain error
        precision = 0.0
    precision = max(precision, 0)
    precision = width - math.floor(precision)
    precision -= 3 if f < 0 else 2  # sign and decimal point
    return max(precision, 1)


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
        >>> split_charge('Formula_2-')
        ('Formula', -2)
        >>> split_charge('Formula-2')
        ('Formula', -2)

    """
    charge = 0
    # matches ]2+, _2-
    m = _CHARGE_PATTERN.search(formula)
    if m:
        m_delim, m_count, m_sign = m.groups()
        if m_count == '':
            charge = int(f'{m_sign}1')  # unreachable?
        else:
            charge = int(f'{m_sign}{m_count}')
    else:
        # matches +2, -2
        m = _CHARGE_NO_DELIM_PATTERN.search(formula)
        if m:
            m_sign, m_count = m.groups()
            m_delim = f'{m_sign}{m_count}'
            charge = int(m_delim)
        else:
            # matches _+, _--, +, --
            m = _CHARGE_NO_COUNT_PATTERN.search(formula)
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
            formula = formula[: -len(m_sign)]
        elif m_delim == ']':
            formula = formula.rsplit(']', 1)[0] + ']'
        else:
            formula = formula.rsplit(m_delim, 1)[0]
        if formula.startswith('[') and formula.endswith(']'):
            # only strip if first bracket closes at the end
            level = 0
            for i, char in enumerate(formula):
                if char in '([{<':
                    level += 1
                elif char in ')]}>':
                    level -= 1
                    if level == 0 and i < len(formula) - 1:
                        break
            else:
                formula = formula[1:-1]
    return formula, charge


def join_charge(formula: str, charge: int, separator: str = '', /) -> str:
    """Return formula with charge appended.

    Parameters:
        formula:
            Chemical formula without charge.
        charge:
            Charge number.
        separator:
            Character separating formula from charge. One of '' or '_'.

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


def format_charge(charge: int, prefix: str = '', /) -> str:
    """Return string representation of charge.

    Parameters:
        charge:
            Charge number.
        prefix:
            Character placed in front of charge in case ``abs(charge)>1``.

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


def unit_factor[U: str](units: dict[U, float], unit: U, /) -> float:
    """Return conversion factor of unit or raise ValueError."""
    try:
        return units[unit]
    except KeyError:
        msg = f'unknown unit {unit!r}, must be one of {tuple(units)}'
        raise ValueError(msg) from None


MassUnit = Literal[
    'g',
    'mg',
    'μg',
    'ug',
    'ng',
    'kg',
    'Da',
    'kDa',
    'MDa',
]
"""Mass unit symbols."""

VolumeUnit = Literal[
    'L',
    'mL',
    'μL',
    'uL',
    'nL',
    'dL',
    'm3',
    'cm3',
    'dm3',
]
"""Volume unit symbols."""

DensityUnit = Literal[
    'g/mL',
    'g/cm3',
    'kg/L',
    'kg/dm3',
    'g/L',
    'kg/m3',
    'mg/mL',
]
"""Density unit symbols."""

MASS_UNIT: dict[MassUnit, float] = {
    'g': 1.0,
    'mg': 1e-3,
    'μg': 1e-6,
    'ug': 1e-6,  # support 'ug' as alternative to 'μg'
    'ng': 1e-9,
    'kg': 1e3,
    'Da': 1.66053906892e-24,  # atomic mass unit to grams, CODATA 2022
    'kDa': 1.66053906892e-21,
    'MDa': 1.66053906892e-18,
}
"""Mapping of mass unit symbols to conversion factors to grams."""

VOLUME_UNIT: dict[VolumeUnit, float] = {
    'L': 1.0,
    'mL': 1e-3,
    'μL': 1e-6,
    'uL': 1e-6,  # support 'uL' as alternative to 'μL'
    'nL': 1e-9,
    'dL': 1e-1,
    'm3': 1e3,
    'cm3': 1e-3,
    'dm3': 1.0,
}
"""Mapping of volume unit symbols to conversion factors to liters."""

DENSITY_UNIT: dict[DensityUnit, float] = {
    'g/mL': 1.0,
    'g/cm3': 1.0,
    'kg/L': 1.0,
    'kg/dm3': 1.0,
    'g/L': 1e-3,
    'kg/m3': 1e-3,
    'mg/mL': 1e-3,
}
"""Mapping of density unit symbols to conversion factors to g/mL."""

# precompile regex patterns
_DEUTERIUM_PATTERN = re.compile(r'(D)(?![a-z])')
_ARITHMETIC_PATTERN = re.compile(r'(?:\+|^)((\d+)\*?(.*?))(?:(?=\+)|$)')
_ARITHMETIC_INVALID_PATTERN = re.compile(r'(\+[a-z])')
_CHARGE_PATTERN = re.compile(r'([\]_])([0-9]{1,})([+-]{1,})$')
_CHARGE_NO_DELIM_PATTERN = re.compile(r'([+-]{1})([0-9]{1,})$')
_CHARGE_NO_COUNT_PATTERN = re.compile(r'([\]_]?)([+-]{1,})$')

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
"""Mapping of common chemical group names to Hill notation."""

ADDUCTS: dict[str, str] = {
    'ACN': 'C2H3N',  # Acetonitrile
    'DMSO': 'C2H6OS',  # Dimethyl sulfoxide
    'FA': 'CH2O2',  # Formic Acid
    'HAc': 'CH4CO2',  # Acetate
    'IsoProp': 'C3H8O',  # Isopropyl alcohol
    'TFA': 'C2HF3O2',  # Trifluoroacetic acid
}
"""Mapping of common adduct names to Hill notation."""

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
"""Mapping of amino acid one-letter codes to Hill notation.

H2O is removed from formulas to represent polymerization.
"""

DEOXYNUCLEOTIDES: dict[str, str] = {
    'A': 'C10H12N5O5P',
    'T': 'C10H13N2O7P',
    'C': 'C9H12N3O6P',
    'G': 'C10H12N5O6P',
}
"""Mapping of Deoxynucleoside monophosphates one-letter codes to Hill notation.

H2O is removed from formulas to represent polymerization.
"""

NUCLEOTIDES: dict[str, str] = {
    'A': 'C10H12N5O6P',
    'U': 'C9H11N2O8P',
    'C': 'C9H12N3O7P',
    'G': 'C10H12N5O7P',
}
"""Mapping of Nucleoside monophosphate one-letter codes to Hill notation.

H2O is removed from formulas to represent polymerization.
"""

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

PREPROCESSORS: dict[str, Callable[[str], str]] = {
    'peptide': from_peptide,
    'ssdna': lambda x: from_oligo(x, 'ssdna'),
    'dsdna': lambda x: from_oligo(x, 'dsdna'),
    'ssrna': lambda x: from_oligo(x, 'ssrna'),
    'dsrna': lambda x: from_oligo(x, 'dsrna'),
}
"""Mapping of formula preprocessors."""


def main(argv: list[str] | None = None, /) -> int:
    """Command line usage main function.

    Parameters:
        argv: Command line arguments.

    """
    if argv is None:
        argv = sys.argv[1:]

    import argparse

    def search_doc(r: str, d: str) -> str:
        if not __doc__:
            return d
        match = re.search(r, __doc__)
        if match is None:
            return d
        return match.group(1)

    parser = argparse.ArgumentParser(
        usage='%(prog)s [options] formula',
        description=search_doc('\n\n([^|]*?)\n\n', ''),
        prog='molmass',
    )
    opt = parser.add_argument
    opt(
        '--version',
        action='version',
        version='%(prog)s ' + search_doc(':Version: (.*)', 'Unknown'),
    )
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
    opt('-v', '--verbose', dest='verbose', action='store_true', default=False)
    opt('formula', nargs='*', help=argparse.SUPPRESS)

    settings = parser.parse_intermixed_args(argv)
    formula_list = settings.formula

    if settings.web:
        try:
            from .web import main as web_main
        except ImportError:
            from web import main as web_main  # type: ignore[no-redef]
        form = {'q': ''.join(formula_list)} if formula_list else None

        return web_main(
            url=settings.url, open_browser=not settings.nobrowser, form=form
        )
    if formula_list:
        formula = ''.join(formula_list)
    else:
        parser.error('no formula specified')

    results = analyze(formula)
    print('\n', results, sep='')  # noqa: T201

    return 0


if __name__ == '__main__':
    sys.exit(main())
