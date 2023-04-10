Molecular mass calculations
===========================

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
:Version: 2023.4.10
:DOI: `10.5281/zenodo.7135495 <https://doi.org/10.5281/zenodo.7135495>`_

Quickstart
----------

Install the molmass package and all dependencies from the
`Python Package Index <https://pypi.org/project/molmass/>`_::

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

This revision was tested with the following requirements and dependencies
(other versions may work):

- `CPython <https://www.python.org>`_ 3.8.10, 3.9.13, 3.10.11, 3.11.3
- `Flask <https://pypi.org/project/Flask/>`_ 2.2.3 (optional)
- `Pandas <https://pypi.org/project/pandas/>`_ 1.5.3 (optional)
- `wxPython <https://pypi.org/project/wxPython/>`_ 4.2.0 (optional)

Revisions
---------

2023.4.10

- Support rdkit-style ionic charges (#11, #12).
- Enable multiplication without addition in from_string.

2022.12.9

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
