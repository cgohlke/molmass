# molmass.py

# Copyright (c) 1990-2020, Christoph Gohlke
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

"""Molecular Mass Calculations.

Molmass is a Python library and console script to calculate the molecular mass
(average, nominal, and isotopic pure), the elemental composition, and the
mass distribution spectrum of a molecule given by its chemical formula,
relative element weights, or sequence.

Calculations are based on the isotopic composition of the elements. Mass
deficiency due to chemical bonding is not taken into account.

Examples of valid formulas are ``H2O``, ``[2H]2O``, ``CH3COOH``, ``EtOH``,
``CuSO4.5H2O``, ``(COOH)2``, ``AgCuRu4(H)2[CO]12{PPh3}2``, ``CGCGAATTCGCG``,
and ``MDRGEQGLLK``.

Formulas are case sensitive and ``+`` denotes the arithmetic operator,
not an ion charge.

For command line usage run ``python -m molmass --help``

:Author: `Christoph Gohlke <https://www.lfd.uci.edu/~gohlke/>`_

:License: BSD 3-Clause

:Version: 2020.1.1

Requirements
------------
* `CPython >= 3.6 <https://www.python.org>`_

Revisions
---------
2020.1.1
    Update elements.py to version 2020.1.1.
    Remove support for Python 2.7 and 3.5.
    Update copyright.
2018.8.15
    Move modules into molmass package.
2018.5.29
    Add option to start web interface from console.
2018.5.25
    Style and docstring fixes.
    Make 'from_fractions' output deterministic.
2005.x.x
    Initial release.

Examples
--------
>>> from molmass import Formula
>>> f = Formula('D2O')  # heavy water
>>> f.formula  # hill notation
'[2H]2O'
>>> f.empirical
'[2H]2O'
>>> f.mass  # average mass
20.02760855624
>>> f.isotope.massnumber  # nominal mass
20
>>> f.isotope.mass  # monoisotopic mass
20.02311817581
>>> f.atoms
3
>>> print(f.composition())
Element  Number  Relative mass  Fraction %
2H            2       4.028204     20.1133
O             1      15.999405     79.8867
Total:        3      20.027609    100.0000
>>> print(f.spectrum())
Relative mass    Fraction %      Intensity
20.023118         99.757000     100.000000
21.027335          0.038000       0.038093
22.027363          0.205000       0.205499

"""

__version__ = '2020.1.1'

__all__ = (
    'analyze', 'Formula', 'FormulaError', 'Spectrum', 'Composition',
    'test', 'main', 'from_fractions', 'from_elements', 'from_oligo',
    'from_peptide', 'from_sequence', 'from_string', 'hill_sorted',
    'GROUPS', 'AMINOACIDS', 'DEOXYNUCLEOTIDES', 'PREPROCESSORS'
)

import sys
import re
import math
import copy
from functools import reduce

try:
    from .elements import ELEMENTS, Isotope
except ImportError:
    from elements import ELEMENTS, Isotope


def analyze(formula, maxatoms=250):
    """Return analysis of chemical formula as string.

    Calculate mass spectrum if number of atoms is smaller than maxatoms.

    """
    result = []
    try:
        f = Formula(formula)

        if len(str(f)) <= 50:
            result.append(f'Formula: {f}')
        if formula != f.formula:
            result.append(f'Hill notation: {f.formula}')
        if f.formula != f.empirical:
            result.append(f'Empirical formula: {f.empirical}')

        prec = precision_digits(f.mass, 9)
        if f.mass != f.isotope.mass:
            result.append(f'\nAverage mass: {f.mass:.{prec}f}')
        result.extend((
            f'Monoisotopic mass: {f.isotope.mass:.{prec}f}',
            f'Nominal mass: {f.isotope.massnumber}'))

        c = f.composition()
        if len(c) > 1:
            result.extend(('\nElemental Composition\n', str(c)))

        if f.atoms < maxatoms:
            s = f.spectrum()
            if len(s) > 1:
                result.extend((
                    '\nMass Distribution',
                    f'\nMost abundant mass: {s.peak[0]:.{prec}f} '
                    f'({s.peak[1] * 100:.3f}%)',
                    f'Mean mass: {s.mean:.{prec}f}\n',
                    str(s)
                ))

    except Exception as exc:
        result.append(f'Error: {exc}')

    return '\n'.join(result)


class lazyattr:
    """Lazy object attribute whose value is computed on first access."""

    def __init__(self, func):
        self.func = func
        # crude hack to keep docstrings and allow doctests
        if func.__doc__:
            self.__doc__ = func.__doc__
            lazyattr.docstrings.__doc__ += '\n\n' + func.__doc__

    def __get__(self, instance, owner):
        if instance is None:
            return self
        result = self.func(instance)
        if result is NotImplemented:
            return getattr(super(owner, instance), self.func.__name__)
        setattr(instance, self.func.__name__, result)
        return result

    @staticmethod
    def docstrings():
        """Docstrings of lazy attributes get appended to this string."""


class Formula:
    """Chemical formula.

    Input string may contain only symbols of chemical elements and isotopes,
    parentheses, and numbers.
    Calculate various properties from formula string, such as hill notation,
    empirical formula, mass, elemental composition, and mass distribution.
    Raise FormulaError on errors.

    Examples
    --------
    >>> Formula('H2O')
    Formula('H2O')

    """

    def __init__(self, formula='', groups=None):
        self._formula = from_string(formula, groups)

    def __str__(self):
        return self._formula

    def __repr__(self):
        return f"Formula('{self._formula}')"

    def __mul__(self, number):
        """Return this formula repeated number times as new Formula.

        Examples
        --------
        >>> Formula('H2O') * 2
        Formula('(H2O)2')

        """
        if not isinstance(number, int) or number < 1:
            raise TypeError('can only multipy with positive number')
        return Formula(f'({self._formula}){number}')

    def __rmul__(self, number):
        """Return this formula repeated number times as new Formula.

        Examples
        --------
        >>> 2 * Formula('H2O')
        Formula('(H2O)2')

        """
        return self.__mul__(number)

    def __add__(self, other):
        """Add this and other formula and return as new Formula.

        Examples
        --------
        >>> Formula("H2O") + Formula("H2O")
        Formula('(H2O)(H2O)')

        """
        if not isinstance(other, Formula):
            raise TypeError('can only add Formula instance')
        return Formula(f'({self})({other})')

    def __sub__(self, other):
        """Subtract elements of other formula and return as new Formula.

        Examples
        --------
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
        return Formula(from_elements(_elements))

    @lazyattr
    def _elements(self):
        """Return number of atoms and isotopes by element as dict.

        Return type is dict{element symbol: dict{mass number: count}}

        Examples
        --------
        >>> Formula('H')._elements
        {'H': {0: 1}}
        >>> pprint(Formula('[2H]2O')._elements)  # pprint sorts dict
        {'H': {2: 2}, 'O': {0: 1}}

        """
        formula = self._formula
        if not formula:
            raise FormulaError('empty formula', formula, 0)

        validchars = set('([{<123456789ABCDEFGHIKLMNOPRSTUVWXYZ')

        if not formula[0] in validchars:
            raise FormulaError(
                f"unexpected character '{formula[0]}'", formula, 0
            )

        validchars |= set(']})>0abcdefghiklmnoprstuy')

        elements = {}
        ele = ''      # parsed element
        num = 0       # number
        level = 0     # parenthesis level
        counts = [1]  # parenthesis level multiplication
        i = len(formula)
        while i:
            i -= 1
            char = formula[i]
            if char not in validchars:
                raise FormulaError(f"unexpected character {char}'", formula, i)
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
                num = int(formula[i:j + 1])
                if num == 0:
                    raise FormulaError('count is zero', formula, i)
            elif char.islower():
                if not formula[i - 1].isupper():
                    raise FormulaError(
                        f"unexpected character '{char}'", formula, i)
                ele = char
            elif char.isupper():
                ele = char + ele
                if num == 0:
                    num = 1
                if ele not in ELEMENTS:
                    raise FormulaError(f"unknown symbol '{ele}'", formula, i)
                iso = ''
                j = i
                while i and formula[i - 1].isdigit():
                    i -= 1
                    iso = formula[i] + iso
                if iso and i and not formula[i - 1] in '([{<':
                    i = j
                    iso = ''
                if iso:
                    iso = int(iso)
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
                "missing opening parenthesis '([{<'", formula, 0)

        if not elements:
            raise FormulaError('invalid formula', formula, 0)

        return elements

    @lazyattr
    def formula(self):
        """Return formula string in Hill notation.

        Examples
        --------
        >>> Formula('BrC2H5').formula
        'C2H5Br'
        >>> Formula('HBr').formula
        'BrH'
        >>> Formula('[(CH3)3Si2]2NNa').formula
        'C6H18NNaSi4'

        """
        return from_elements(self._elements)

    @lazyattr
    def empirical(self):
        """Return empirical formula in Hill notation.

        The empirical formula has the simplest whole number ratio of atoms
        of each element present in formula.

        Examples
        --------
        >>> Formula('H2O').empirical
        'H2O'
        >>> Formula('S4').empirical
        'S'
        >>> Formula('C6H12O6').empirical
        'CH2O'

        """
        return from_elements(self._elements, self.gcd)

    @lazyattr
    def atoms(self):
        """Return number of atoms.

        Examples
        --------
        >>> Formula('CH3COOH').atoms
        8

        """
        return sum(sum(i.values()) for i in self._elements.values())

    @lazyattr
    def gcd(self):
        """Return greatest common divisor of element counts.

        Examples
        --------
        >>> Formula('H2').gcd
        2
        >>> Formula('H2O').gcd
        1
        >>> Formula('C6H12O6').gcd
        6

        """
        return gcd(
            {list(i)[0] for i in (j.values() for j in self._elements.values())}
        )

    @lazyattr
    def mass(self):
        """Return average relative molecular mass.

        Sums the relative atomic masses of all atoms in molecule.
        Equals the molar mass in g/mol, i.e. the mass of one mole of substance.

        Examples
        --------
        >>> Formula('C').mass
        12.01074
        >>> Formula('12C').mass
        12.0
        >>> print('{:.2f}'.format(Formula('C48H32AgCuO12P2Ru4').mass))
        1438.40

        """
        result = 0.0
        for symbol in self._elements:
            ele = ELEMENTS[symbol]
            for massnumber, count in self._elements[symbol].items():
                if massnumber:
                    result += ele.isotopes[massnumber].mass * count
                else:
                    result += ele.mass * count
        return result

    @lazyattr
    def isotope(self):
        """Return isotope composed of most abundant elemental isotopes.

        Examples
        --------
        >>> print(Formula('C').isotope.mass)
        12.0
        >>> Formula('13C').isotope.massnumber
        13
        >>> print(Formula('C48H32AgCuO12P2Ru4').isotope)
        1440, 1439.5890, 0.205075%

        """
        result = Isotope()
        for symbol in self._elements:
            ele = ELEMENTS[symbol]
            for massnumber, count in self._elements[symbol].items():
                if massnumber:
                    isotope = ele.isotopes[massnumber]
                else:
                    isotope = ele.isotopes[ele.nominalmass]
                result.mass += isotope.mass * count
                result.massnumber += isotope.massnumber * count
                result.abundance *= isotope.abundance ** count
        return result

    def composition(self, isotopic=True):
        """Return elemental composition as Composition instance.

        Return type is tuple(tuple(symbol, count, mass, fraction), ).

        If isotopic is True, isotopes specified in the formula are listed
        separately, otherwise they are listed as part of an element.

        Examples
        --------
        >>> Formula('[12C]C').composition(False)
        (('C', 2, 24.01074, 1.0),)
        >>> for i in Formula('[12C]C').composition(True): print(i)
        ('C', 1, 12.01074, 0.5002236499166623)
        ('12C', 1, 12.0, 0.49977635008333776)

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
        return Composition(result)

    def spectrum(self, minfract=1e-9):
        """Return low resolution mass spectrum as Spectrum instance.

        Return type is dict{massnumber: list[mass, fraction]}.

        Calculated by combining the mass numbers of the elemental isotopes.

        Examples
        --------
        >>> def _(spectrum):
        ...     for key, val in spectrum.items():
        ...         print(f'{key}, {val[0]:.4f}, {val[1]*100:.6f}%')
        >>> _(Formula('D').spectrum())
        2, 2.0141, 100.000000%
        >>> _(Formula('H').spectrum())
        1, 1.0078, 99.988500%
        2, 2.0141, 0.011500%
        >>> _(Formula('D2').spectrum())
        4, 4.0282, 100.000000%
        >>> _(Formula('DH').spectrum())
        3, 3.0219, 99.988500%
        4, 4.0282, 0.011500%
        >>> _(Formula('H2').spectrum())
        2, 2.0157, 99.977001%
        3, 3.0219, 0.022997%
        4, 4.0282, 0.000001%
        >>> _(Formula('DHO').spectrum())
        19, 19.0168, 99.745528%
        20, 20.0215, 0.049468%
        21, 21.0211, 0.204981%
        22, 22.0274, 0.000024%

        """
        spectrum = {0: [0.0, 1.0]}
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
                        if t[1] < minfract:
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
                            if t[1] < minfract:
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
        return Spectrum(spectrum)


class Spectrum(dict):
    """Mass distribution.

    Basic type is dict{bin: list[mass, fraction]}.
    Iterators over the dict are sorted by bin/mass.

    Examples
    --------
    >>> print(Spectrum({1: [1.078, 0.9999], 2: [2.014, 0.0001]}))
    Relative mass    Fraction %      Intensity
    1.0780000         99.990000     100.000000
    2.0140000          0.010000       0.010001

    """

    def __init__(self, *args, **kwds):
        dict.__init__(self, *args, **kwds)
        self._sorted_keys = sorted(dict.keys(self))

    def __iter__(self):
        return iter(self._sorted_keys)

    def keys(self):
        return iter(self._sorted_keys)

    def values(self):
        return (self[key] for key in self._sorted_keys)

    def items(self):
        return ((key, self[key]) for key in self._sorted_keys)

    @lazyattr
    def range(self):
        """Return smallest and largest bin number."""
        return min(dict.keys(self)), max(dict.keys(self))

    @lazyattr
    def peak(self):
        """Return most abundant mass and fraction."""
        mass = 0.0
        fraction = 0.0
        for m, f in dict.values(self):
            if f > fraction:
                fraction = f
                mass = m
        return mass, fraction

    @lazyattr
    def mean(self):
        """Return mean of all masses in spectrum."""
        return sum((mass * fraction) for mass, fraction in dict.values(self))

    def __str__(self):
        if len(self) == 0:
            return ''
        result = ['Relative mass    Fraction %      Intensity']
        prec = precision_digits(self.peak[0], 9)
        norm = 100.0 / self.peak[1]
        for mass, fraction in self.values():
            result.append(
                '{:<13.{}f}   {:11.6f}   {:12.6f}'.format(
                    mass, prec, fraction * 100.0, fraction * norm
                )
            )
        return '\n'.join(result)


class Composition(tuple):
    """Elemental composition.

    Basic type is tuple(tuple(symbol, count, mass, fraction), ).

    Examples
    --------
    >>> print(Composition((('2H', 2, 4.028, 0.201), ('O', 1, 15.999, 0.799))))
    Element  Number  Relative mass  Fraction %
    2H            2       4.028000     20.1000
    O             1      15.999000     79.9000
    Total:        3      20.027000    100.0000

    """

    @lazyattr
    def total(self):
        """Return sums of counts, masses, and fractions."""
        result = [0, 0.0, 0.0]
        for item in self:
            result[0] += item[1]
            result[1] += item[2]
            result[2] += item[3]
        return tuple(result)

    def __str__(self):
        if len(self) == 0:
            return ''
        prec = precision_digits(self.total[1], 9)
        result = ['Element  Number  Relative mass  Fraction %']
        for symbol, count, mass, fraction in self:
            result.append(
                '{:<6s} {:8}  {:13.{}f} {:11.4f}'.format(
                    symbol, count, mass, prec, fraction * 100
                )
            )
        if len(self) > 1:
            count, mass, fraction = self.total
            result.append(
                '{:<6s} {:8}  {:13.{}f} {:11.4f}'.format(
                    'Total:', count, mass, prec, fraction * 100
                )
            )
        return '\n'.join(result)


class FormulaError(Exception):
    """Custom exception to report errors in the Formula object."""

    def __init__(self, msg, formula='', pos=-1):
        self.position = pos
        self.message = msg
        self.formula = formula
        Exception.__init__(self, msg, formula, pos)

    def __str__(self):
        if self.position < 0:
            return str(self.message)
        return f"{self.message}\n{self.formula}\n{'.' * self.position}^"


def from_string(formula, groups=None):
    """Return formula string from user input string.

    Return string should be composed of chemical elements/isotopes,
    parentheses, and numbers only. Raise FormulaError on errors.

    Supports simple, non-nested, arithmetic (+ and *), abbreviations of
    common chemical groups, peptides, oligos, and mass fractions.

    Examples
    --------
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
        raise FormulaError('not a string') from exc

    # abbreviations of common chemical groups
    if groups is None:
        groups = GROUPS
    if groups:
        for grp in reversed(sorted(groups)):
            formula = formula.replace(grp, f'({groups[grp]})')

    # list of mass fractions
    if ':' in formula and ',' in formula:
        fractions = {}
        try:
            for item in formula.split(','):
                item = item.split(':')
                fractions[item[0].strip()] = float(item[1].strip())
        except Exception as exc:
            raise FormulaError(
                'invalid list of mass fractions', formula) from exc
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

    # arithmetic
    formula = formula.replace('.', '+')
    if '+' in formula:
        for match in re.findall(
            r'(?:\+|^)((\d+)\*?(.*?))(?:(?=\+)|$)', formula
        ):
            formula = formula.replace(match[0], f'({match[2]}){match[1]}')
        formula = formula.replace('+', '')
    if '-' in formula:
        FormulaError('subtraction not supported yet')

    return formula


def from_elements(elements, divisor=1, *fmt):
    """Return formula string in Hill notation from elements dict.

    Examples
    --------
    >>> from_elements({'C': {0: 4, 12: 2}}, 2)
    'C2[12C]'
    >>> from_elements({'C': {0: 4, 12: 2}}, 2, '{}', '{}<sub>{}</sub>',
    ...     '<sup>{}</sup>{}', '<sup>{}</sup>{}<sub>{}</sub>')
    'C<sub>2</sub><sup>12</sup>C'

    """
    if not fmt:
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
    return ''.join(formula)


def from_fractions(fractions, maxcount=10, precision=1e-4):
    """Return formula string from elemental mass fractions.

    Examples
    --------
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
                raise FormulaError(f"unknown element '{symbol}'") from exc
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
                    f"unknown isotope '[{massnum}{symbol}]'") from exc
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


def from_sequence(sequence, items):
    """Translate sequence using items dict and return histogram of items in
    translated sequence as formula string.

    Raise KeyError if a sequence item is not in items.

    Examples
    --------
    >>> from_sequence('A', {'A': 'B'})
    '(B)'
    >>> from_sequence('AA', {'A': 'B'})
    '(B)2'

    """
    counts = {key: 0 for key in items}
    for item in sequence:
        counts[item] += 1
    formula = []
    for key in sorted(items):
        num = counts[key]
        if num == 1:
            formula.append(f'({items[key]})')
        elif num:
            formula.append(f'({items[key]}){num}')
    return ''.join(formula)


def from_peptide(sequence):
    """Return chemical formula for polymer of unmodified amino acids.

    >>> from_peptide('GG')
    '((C2H3NO)2H2O)'
    >>> f = Formula(from_peptide('GPAVL IMCFY WHKRQ NEDST'))
    >>> print(f.formula, f.atoms, round(f.mass, 6))
    C107H159N29O30S2 327 2395.717936

    """
    sequence = sequence.replace(' ', '')
    return f'({from_sequence(sequence, AMINOACIDS)}H2O)'


def from_oligo(sequence, dtype='ssdna'):
    """Return chemical formula for polymer of unmodified (deoxy)nucleotides.

    Dtype is 'ssdna', 'dsdna', 'ssrna', or 'dsrna'.
    Each strand includes a 5' monophosphate.

    Examples
    --------
    >>> from_oligo('AC', 'ssdna')
    '((C10H12N5O5P)(C9H12N3O6P)H2O)'
    >>> from_oligo('AU', 'dsrna')
    '((C10H12N5O6P)2(C9H11N2O8P)2(H2O)2)'

    >>> f = Formula(from_oligo('ATC G', 'dsdna'))
    >>> print(f.formula, f.atoms, round(f.mass, 6))
    C78H102N30O50P8 268 2507.609138
    >>> f = Formula(from_oligo('AUC G', 'ssrna'))
    >>> print(f.formula, f.atoms, round(f.mass, 6))
    C38H49N15O29P4 135 1303.775567

    """
    dtype = dtype.lower()
    sequence = sequence.replace(' ', '')
    items = NUCLEOTIDES if 'rna' in dtype else DEOXYNUCLEOTIDES
    if dtype.startswith('ds'):
        complements = items['complements']
        t = ''.join(complements[i] for i in sequence)
        formula = from_sequence(sequence + t, items)
        formula = f'({formula}(H2O)2)'
    else:
        formula = from_sequence(sequence, items)
        formula = f'({formula}H2O)'
    return formula


def hill_sorted(symbols):
    """Return iterator over element symbols in order of Hill notation.

    Examples
    --------
    >>> for i in hill_sorted('HCO'): print(i, end='')
    CHO

    """
    symbols = set(symbols)
    if 'C' in symbols:
        symbols.remove('C')
        yield 'C'
        if 'H' in symbols:
            symbols.remove('H')
            yield 'H'
    yield from sorted(symbols)


def gcd(numbers):
    """Return greatest common divisor of integer numbers.

    Using Euclid's algorithm.

    Examples
    --------
    >>> gcd([4])
    4
    >>> gcd([3, 6])
    3
    >>> gcd([6, 7])
    1

    """

    def _gcd(a, b):
        """Return greatest common divisor of two integer numbers."""
        while b:
            a, b = b, a % b
        return a

    return reduce(_gcd, numbers)


def precision_digits(f, width):
    """Return number of digits after decimal point to print f in width chars.

    Examples
    --------
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


# Common chemical groups
GROUPS = {
    'Abu': 'C4H7NO',
    'Acet': 'C2H3O',
    'Acm': 'C3H6NO',
    'Adao': 'C10H15O',
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
    'Tips': 'C9H21Si',
    'Tms': 'C3H9Si',
    'Tos': 'C7H7O2S',
    'Trp': 'C11H10N2O',
    'Trpp': 'C11H9N2O',
    'Trt': 'C19H15',
    'Tyr': 'C9H9NO2',
    'Tyrp': 'C9H8NO2',
    'Val': 'C5H9NO',
    'Valoh': 'C5H9NO2',
    'Valohp': 'C5H8NO2',
    'Xan': 'C13H9O',
}

# Amino acids - H2O
AMINOACIDS = {
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

# Deoxynucleotide monophosphates - H2O
DEOXYNUCLEOTIDES = {
    'A': 'C10H12N5O5P',
    'T': 'C10H13N2O7P',
    'C': 'C9H12N3O6P',
    'G': 'C10H12N5O6P',
    'complements': {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'},
}

# Nucleotide monophosphates - H2O
NUCLEOTIDES = {
    'A': 'C10H12N5O6P',
    'U': 'C9H11N2O8P',
    'C': 'C9H12N3O7P',
    'G': 'C10H12N5O7P',
    'complements': {'A': 'U', 'U': 'A', 'C': 'G', 'G': 'C'},
}

# Formula preprocessors
PREPROCESSORS = {
    'peptide': from_peptide,
    'ssdna': lambda x: from_oligo(x, 'ssdna'),
    'dsdna': lambda x: from_oligo(x, 'dsdna'),
    'ssrna': lambda x: from_oligo(x, 'ssrna'),
    'dsrna': lambda x: from_oligo(x, 'dsrna'),
}


def test(verbose=False):
    """Test the module and the examples in docstrings."""
    import doctest

    doctest.testmod(verbose=verbose)

    # these formulas should pass
    for data in [
        (''.join(e.symbol for e in ELEMENTS), '', 14693.181589000998),
        ('12C', '[12C]', 12.0),
        ('12CC', 'C[12C]', 24.0107),
        ('Co(Bpy)(CO)4', '', 327.16),
        ('CH3CH2Cl', 'C2H5Cl', 64.5147),
        ('C1000H1000', 'CH', 13018.68),
        ('Ru2(CO)8', 'C4O4Ru', 426.2232),
        ('RuClH(CO)(PPh3)3', 'C55H46ClOP3Ru', 952.41392),
        ('PhSiMe3', 'C9H14Si', 150.29566),
        ('Ph(CO)C(CH3)3', 'C11H14O', 162.23156),
        ('HGlyGluTyrOH', 'C16H21N3O7', 367.35864),
        ('HCysTyrIleGlnAsnCysProLeuNH2', 'C41H65N11O11S2', 952.1519),
        ('HCysp(Trt)Tyrp(Tbu)IleGlnp(Trt)Asnp(Trt)ProLeuGlyNH2',
         'C101H113N11O11S', 1689.13532),
        ('CGCGAATTCGCG', 'C116H148N46O73P12', 3726.4),
        ('MDRGEQGLLK', 'C47H83N15O16S', 1146.3),
        ('CDCl3', 'C[2H]Cl3', 120.384),
        ('[13C]Cl4', '[13C]Cl4', 154.8153),
        ('C5(PhBu(EtCHBr)2)3', 'C53H78Br6', 1194.626),
        ('AgCuRu4(H)2[CO]12{PPh3}2', 'C48H32AgCuO12P2Ru4', 1438.4022),
        ('PhNH2.HCl', 'C6H8ClN', 129.5892),
        ('NH3.BF3', 'BF3H3N', 84.8357),
        ('CuSO4.5H2O', 'CuH10O9S', 249.68),
    ]:
        if verbose:
            print(f"Trying Formula('{data[0]}') ...", end='')
        try:
            f = Formula(data[0])
            f.empirical
            f.mass
            f.spectrum
        except FormulaError as exc:
            print('Error:', exc)
            continue
        if data[1] and f.empirical != data[1]:
            print(
                "Failure for {}:\n    Expected '{}', got '{}':".format(
                    data[0], data[1], f.empirical
                )
            )
            continue
        if data[2] and abs(f.mass - data[2]) > 0.1:
            print(
                'Failure for {}:\n    Expected {}, got {}'.format(
                    data[0], data[2], f.mass
                )
            )
            continue
        if verbose:
            print('ok')

    # these formulas are expected to fail
    for data in [
        '', '()', '2', 'a', '(a)', 'C:H', 'H:', 'C[H', 'H)2',
        'A', 'Aa', '2lC', '1C', '[11C]', 'H0', '()0', '(H)0C',
        'Ox: 0.26, 30Si: 0.74',
    ]:
        if verbose:
            print(f"Trying Formula('{data}') ...", end='')
        try:
            f = Formula(data).empirical
        except FormulaError as exc:
            if verbose:
                print('ok\nExpected error:', exc)
        else:
            print(
                f"Failure expected for '{data}', got '{Formula(data).formula}'"
            )


def main(argv=None):
    """Command line usage main function."""
    if argv is None:
        argv = sys.argv

    import optparse

    def search_doc(r, d):
        return re.search(r, __doc__).group(1) if __doc__ else d

    parser = optparse.OptionParser(
        usage='usage: %prog [options] formula',
        description=search_doc('\n\n([^|]*?)\n\n', ''),
        version='%prog {}'.format(search_doc(':Version: (.*)', 'Unknown')),
        prog='molmass'
    )
    opt = parser.add_option
    opt('--web', dest='web', action='store_true', default=False,
        help='start web application and open it in a web browser')
    opt('--test', dest='test', action='store_true', default=False,
        help='test the module')
    opt('-v', '--verbose', dest='verbose', action='store_true', default=False)

    settings, formula = parser.parse_args()

    if settings.web:
        try:
            from . import molmass_web  # noqa: delayed import
        except ImportError:
            import molmass_web  # noqa: delayed import

        return molmass_web.main()
    if settings.test:
        test(settings.verbose)
        return 0
    if formula:
        formula = ''.join(formula)
    else:
        parser.error('no formula specified')

    try:
        results = analyze(formula)
    except Exception as exc:
        print('\nError: \n  ', exc, sep='')
        raise exc
    else:
        print('\n', results, sep='')


if __name__ == '__main__':
    from pprint import pprint  # noqa: used by doctest

    sys.exit(main())
