..
  This file is updated by setup.py. Do not edit.

Revisions
=========

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

- Improve code quality.

2025.12.12

- Make boolean arguments keyword-only (breaking).

2025.11.11

- Allow empty formulas (breaking).
- Derive FormulaError from ValueError.
- Move tests to separate test module.

2025.9.4

- Precompile regex patterns.
- Remove doctest command line option.
- Drop support for Python 3.10, support Python 3.14.

2025.4.14

- Add mass_charge_ratio helper function (#17).
- Drop support for Python 3.9.

2024.10.25

- Fix composition of formula with multiple isotopes of same element (#16).

2024.5.24

- Fix docstring examples not correctly rendered on GitHub.

2024.5.10

- Add options to disable parsing groups, oligos, fractions, arithmetic (#14).
- Add Formula.expanded property.

2023.8.30

- Fix linting issues.
- Add py.typed marker.
- Drop support for Python 3.8.

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
- Drop support for Python 3.7.

2021.6.18

- Add Particle types to elements (#5).
- Fix molmass_web failure on WSL2 (#9).
- Fix elements_gui layout issue.
- Drop support for Python 3.6.

2020.6.10

- Fix elements_gui symbol size on WSL2.
- Support wxPython 4.1.

2020.1.1

- Update elements atomic weights and isotopic compositions from NIST.
- Move element descriptions into separate module.
- Drop support for Python 2.7 and 3.5.

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
