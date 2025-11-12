# test_molmass.py

# Copyright (c) 1990-2025, Christoph Gohlke
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# * Redistributions of source code must retain the above copyright
#   notice, this list of conditions and the following disclaimer.
# * Redistributions in binary form must reproduce the above copyright
#   notice, this list of conditions and the following disclaimer in the
#   documentation and/or other materials provided with the distribution.
# * Neither the name of the copyright holders nor the names of any
#   contributors may be used to endorse or promote products derived
#   from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.

"""Unittests for the molmass package.

:Version: 2025.11.11

"""

import pytest

import molmass
from molmass import (
    AMINOACIDS,
    DEOXYNUCLEOTIDES,
    ELEMENTS,
    GROUPS,
    NUCLEOTIDES,
    PREPROCESSORS,
    PROTON,
    Composition,
    CompositionItem,
    Element,
    Formula,
    FormulaError,
    Isotope,
    Spectrum,
    SpectrumEntry,
    __version__,
    analyze,
    format_charge,
    from_elements,
    from_fractions,
    from_oligo,
    from_peptide,
    from_sequence,
    from_string,
    hill_sorted,
    join_charge,
    main,
    mass_charge_ratio,
    split_charge,
)
from molmass.molmass import gcd, precision_digits


@pytest.mark.skipif(__doc__ is None, reason='__doc__ is None')
def test_version():
    """Assert molmass versions match docstrings."""
    ver = ':Version: ' + __version__
    assert ver in __doc__
    assert ver in molmass.__doc__


def test_empty():
    """Test Formula('')."""
    with pytest.raises(FormulaError):
        Formula('', allow_empty=False)

    with pytest.raises(FormulaError):
        Formula(' ', allow_empty=False)

    f = Formula('', allow_empty=True)
    assert f.formula == ''
    assert f.expanded == ''
    assert f.empirical == ''
    assert f.atoms == 0
    assert f.gcd == 1
    assert f.charge == 0
    assert f.mass == 0.0
    assert f.monoisotopic_mass == 0.0
    assert f.nominal_mass == 0
    assert f.mz == 0
    assert f.isotope.mass == 0.0
    assert f.isotope.massnumber == 0
    assert f.isotope.abundance == 1.0

    composition = f.composition()
    assert composition.astuple() == ()
    assert composition.asdict() == {}
    assert len(composition) == 0
    assert composition.dataframe().empty
    assert str(composition) == ''
    assert repr(composition) == 'Composition([])'

    spectrum = f.spectrum(min_fraction=1e-9, min_intensity=1e-9)
    assert len(spectrum) == 0
    assert spectrum.asdict() == {}
    assert spectrum.mean == 0.0
    assert spectrum.dataframe().empty
    assert str(spectrum) == ''
    assert repr(spectrum) == 'Spectrum({})'
    with pytest.raises(ValueError):
        spectrum.range
    with pytest.raises(ValueError):
        spectrum.peak


def test_etoh():
    """Test Formula('EtOH')."""
    f = Formula('EtOH')
    assert f.formula == 'C2H6O'
    assert f.expanded == '(C2H5)OH'
    assert f.empirical == 'C2H6O'
    assert f.atoms == 9
    assert f.gcd == 1
    assert f.charge == 0
    assert f.mass == 46.068531
    assert f.monoisotopic_mass == 46.041864812949996
    assert f.nominal_mass == 46
    assert f.mz == 46.068531
    assert f.isotope.mass == 46.041864812949996
    assert f.isotope.massnumber == 46
    assert f.isotope.abundance == 0.9756627354527866

    composition = f.composition()
    assert str(composition).startswith(
        'Element  Count  Relative mass  Fraction %'
    )
    assert repr(composition).startswith('<Composition(')
    assert len(composition) == 3
    assert composition.astuple() == (
        ('C', 2, 24.02148, 0.5214292593788155),
        ('H', 6, 6.047646, 0.13127499116479316),
        ('O', 1, 15.999405, 0.34729574945639136),
    )
    assert not composition.dataframe().empty

    spectrum = f.spectrum(min_fraction=1e-9, min_intensity=1e-9)
    assert isinstance(spectrum[46], SpectrumEntry)
    assert str(spectrum).startswith(
        'A   Relative mass  Fraction %  Intensity %'
    )
    assert repr(spectrum).startswith('<Spectrum({46: [46, 46.041')
    assert len(spectrum) == 6
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
    assert not spectrum.dataframe().empty

    assert abs(spectrum.mean - 46.06852122027406) < 1e-9, spectrum.mean
    assert spectrum.peak.astuple() == (
        46,
        46.04186481295,
        0.9756627354527866,
        100.0,
        46.04186481295,
    ), spectrum.peak.astuple()


@pytest.mark.parametrize(
    'formula, empirical, mass',
    [
        (
            ''.join(e.symbol for e in ELEMENTS),
            'CHAcAgAlAmArAsAtAuBBaBeBhBiBkBrCaCdCeCfClCmCoCrCsCuDbDyErEsEuFFe'
            'FmFrGaGdGeHeHfHgHoHsIInIrKKrLaLiLrLuMdMgMnMoMtNNaNbNdNeNiNoNpOOs'
            'PPaPbPdPmPoPrPtPuRaRbReRfRhRnRuSSbScSeSgSiSmSnSrTaTbTcTeThTiTlTm'
            'UVWXeYYbZnZr',
            14693.181589000998,
        ),
        ('', '', 0.0),
        ('1H+', '[[1H]]+', PROTON.mass),
        ('12C', '[12C]', 12.0),
        ('12CC', 'C[12C]', 24.01074),
        ('SO4_2-', '[O4S]2-', 96.06351715981813),
        ('[SO4]2_4-', '[O4S]2-', 192.12703431963627),
        ('[CHNOP[13C]]2-', '[C[13C]HNOP]2-', 87.003),
        ('[CHNOP[13C]]_2-', '[C[13C]HNOP]2-', 87.003),
        ('CHNOP[13C]-2', '[C[13C]HNOP]2-', 87.003),
        ('Co(Bpy)(CO)4', 'C14H8CoN2O4', 327.158108),
        ('CH3CH2Cl', 'C2H5Cl', 64.514085),
        ('C1000H1000', 'CH', 13018.68),
        ('Ru2(CO)8', 'C4O4Ru', 426.22116),
        ('RuClH(CO)(PPh3)3', 'C55H46ClOP3Ru', 952.399576994),
        ('PhSiMe3', 'C9H14Si', 150.293334),
        ('Ph(CO)C(CH3)3', 'C11H14O', 162.228719),
        ('HGlyGluTyrOH', 'C16H21N3O7', 367.354545),
        ('HCysTyrIleGlnAsnCysProLeuNH2', 'C41H65N11O11S2', 952.153293),
        ('CGCGAATTCGCG', 'C116H148N46O73P12', 3726.371154976),
        ('MDRGEQGLLK', 'C47H83N15O16S', 1146.319708),
        ('CDCl3', 'C[2H]Cl3', 120.38354177811999),
        ('[13C]Cl4', '[13C]Cl4', 154.81495483507),
        ('C5(PhBu(EtCHBr)2)3', 'C53H78Br6', 1194.609618),
        ('AgCuRu4(H)2[CO]12{PPh3}2', 'C48H32AgCuO12P2Ru4', 1438.404215996),
        ('PhNH2.HCl', 'C6H8ClN', 129.587571),
        ('NH3.BF3', 'BF3H3N', 84.8367355),
        ('CuSO4.5H2O', 'CuH10O9S', 249.68485),
        ('5*H2O+CuSO4', 'CuH10O9S', 249.68485),
        ('5*H2O', 'H2O', 90.076435),
        (
            'HCysp(Trt)Tyrp(Tbu)IleGlnp(Trt)Asnp(Trt)ProLeuGlyNH2',
            'C101H113N11O11S',
            1689.114061,
        ),
    ],
)
def test_formulas_mass(formula, empirical, mass):
    """Test valid formulas."""
    f = Formula(formula)
    assert f.empirical == empirical
    assert f.mass == pytest.approx(mass)
    f.spectrum()
    str(f)


@pytest.mark.parametrize(
    'formula',
    [
        # '',
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
        'O2_-2',
        'C+a',
    ],
)
def test_formulas_invalid(formula):
    """Test invalid formulas."""
    with pytest.raises(FormulaError):
        Formula(formula).empirical


@pytest.mark.parametrize(
    'formula, message, detail',
    [
        ('abc', 'unexpected character', 'abc\n^'),
        ('(H2O)2-H2O', 'subtraction not allowed', '(H2O)2-H2O\n......^'),
        ('[11C]', 'unknown isotope', '[11C]\n.^'),
        ('O: 0;26, 30Si: 0.74', 'invalid list of mass fractions', ''),
    ],
)
def test_formula_error(formula, message, detail):
    """Test Formula errors."""
    with pytest.raises(FormulaError) as excinfo:
        Formula(formula).formula
    assert message in str(excinfo.value)
    assert detail in str(excinfo.value)


@pytest.mark.parametrize(
    'formula, expected',
    [
        # elements and counts
        ('H2O', 'H2O'),
        ('D2O', '[2H]2O'),
        # isotopes
        ('[30Si]3O2', '[30Si]3O2'),
        # ion charges
        ('[AsO4]3-', '[AsO4]3-'),
        ('O2-2', '[O2]2-'),
        # abbreviations of chemical groups
        ('EtOH', '(C2H5)OH'),
        # simple arithmetic
        ('(COOH)2', '(COOH)2'),
        ('CuSO4.5H2O', 'CuSO4(H2O)5'),
        # relative element weights
        ('O: 0.26, 30Si: 0.74', 'O2[30Si]3'),
        # nucleotide sequences
        (
            'CGCGAATTCGCG',
            '((C10H12N5O5P)2(C9H12N3O6P)4(C10H12N5O6P)4(C10H13N2O7P)2H2O)',
        ),
        (
            'dsrna(CCUU)',
            '((C10H12N5O6P)2(C9H12N3O7P)2(C10H12N5O7P)2(C9H11N2O8P)2(H2O)2)',
        ),
        # peptide sequences
        (
            'MDRGEQGLLK',
            '((C4H5NO3)(C5H7NO3)(C2H3NO)2(C6H12N2O)'
            '(C6H11NO)2(C5H9NOS)(C5H8N2O2)(C6H12N4O)H2O)',
        ),
        ('peptide(CPK)', '((C3H5NOS)(C6H12N2O)(C5H7NO)H2O)'),
    ],
)
def test_formula_repr(formula, expected):
    """Test repr(Formula)."""
    assert repr(Formula(formula)) == f'Formula({expected!r})'


@pytest.mark.parametrize(
    'formula, expected',
    [
        ('H', {'H': {0: 1}}),
        ('[2H]2O', {'O': {0: 1}, 'H': {2: 2}}),
    ],
)
def test_formula_elements(formula, expected):
    """Test Formula._elements property."""
    assert Formula(formula)._elements == expected


@pytest.mark.parametrize(
    'formula, expected',
    [
        ('BrC2H5', 'C2H5Br'),
        ('[(CH3)3Si2]2NNa', 'C6H18NNaSi4'),
    ],
)
def test_formula_formula(formula, expected):
    """Test Formula.formula property."""
    assert Formula(formula).formula == expected


@pytest.mark.parametrize(
    'formula, expected',
    [
        ('H2O', 'H2O'),
        ('C6H12O6', 'CH2O'),
        ('[SO4]2_4-', '[O4S]2-'),
    ],
)
def test_formula_empirical(formula, expected):
    """Test Formula.empirical property."""
    assert Formula(formula).empirical == expected


@pytest.mark.parametrize(
    'formula, expected',
    [
        ('EtOH', '(C2H5)OH'),
        ('CuSO4.5H2O', 'CuSO4(H2O)5'),
        ('WQ', '((C5H8N2O2)(C11H10N2O)H2O)'),
    ],
)
def test_formula_expanded(formula, expected):
    """Test Formula.expanded property."""
    assert Formula(formula).expanded == expected


@pytest.mark.parametrize(
    'formula, expected',
    [
        ('CH3COOH', 8),
        ('WQ', 44),
    ],
)
def test_formula_atoms(formula, expected):
    """Test Formula.atoms property."""
    assert Formula(formula).atoms == expected


@pytest.mark.parametrize(
    'formula, expected',
    [
        ('H2O', 0),
        ('SO4_2-', -2),
    ],
)
def test_formula_charge(formula, expected):
    """Test Formula.charge property."""
    assert Formula(formula).charge == expected


@pytest.mark.parametrize(
    'formula, expected',
    [
        ('H2O', 1),
        ('H2', 2),
        ('C6H12O6', 6),
    ],
)
def test_formula_gcd(formula, expected):
    """Test Formula.gcd property."""
    assert Formula(formula).gcd == expected


@pytest.mark.parametrize(
    'formula, expected',
    [
        ('H', 1.007941),
        ('H+', 1.007392),
        ('SO4_2-', 96.06351),
        ('12C', 12.0),
        ('C8H10N4O2', 194.1909),
        ('C48H32AgCuO12P2Ru4', 1438.404216),
    ],
)
def test_formula_mass(formula, expected):
    """Test Formula.mass property."""
    assert Formula(formula).mass == pytest.approx(expected)


@pytest.mark.parametrize(
    'formula, expected',
    [
        ('C8H10N4O2', 194.08037),
    ],
)
def test_formula_monoisotopic_mass(formula, expected):
    """Test Formula.monoisotopic_mass property."""
    assert Formula(formula).monoisotopic_mass == pytest.approx(expected)


@pytest.mark.parametrize(
    'formula, expected',
    [
        ('C8H10N4O2', 194),
    ],
)
def test_formula_nominal_mass(formula, expected):
    """Test Formula.nominal_mass property."""
    assert Formula(formula).nominal_mass == expected


@pytest.mark.parametrize(
    'formula, expected',
    [
        ('H', 1.007941),
        ('H+', 1.007392),
        ('SO4_2-', 48.03175),
    ],
)
def test_formula_mz(formula, expected):
    """Test Formula.mz property."""
    assert Formula(formula).mz == pytest.approx(expected)


@pytest.mark.parametrize(
    'formula, expected',
    [
        ('C', Isotope(12.0, 0.9893, 12)),
        ('13C', Isotope(13.003355, 0.0107, 13)),
        ('SO4_2-', Isotope(95.952826812, 0.9407, 96, -2)),
        ('C48H32AgCuO12P2Ru4', Isotope(1439.588966, 0.0020507511, 1440)),
    ],
)
def test_formula_isotope(formula, expected):
    """Test Formula.isotope property."""
    isotope = Formula(formula).isotope
    assert isinstance(isotope, Isotope)
    str(isotope)
    repr(isotope)
    assert isotope.mass == pytest.approx(expected.mass)
    assert isotope.abundance == pytest.approx(expected.abundance)
    assert isotope.massnumber == expected.massnumber
    assert isotope.charge == expected.charge


def test_formula_composition():
    """Test Formula.composition method."""
    composition = Formula('[12C][13C]C').composition()

    composition.astuple()
    composition.asdict()
    composition.dataframe()

    assert 'Element  Count  Relative mass  Fraction %' in str(composition)
    assert '12C' in str(composition)
    assert '13C' in str(composition)
    assert isinstance(composition, Composition)
    assert repr(composition).startswith("<Composition([('C', 1, 12.0")
    assert len(composition) == 3
    assert len(list(composition.items())) == 3
    assert len(list(composition.keys())) == 3
    assert len(list(composition.values())) == 3

    item = composition['C']
    assert isinstance(item, CompositionItem)
    repr(item)
    str(item)
    assert item.symbol == 'C'
    assert item.count == 1
    assert item.mass == pytest.approx(12.010740)
    assert item.fraction == pytest.approx(0.324490982)
    assert item.astuple() == (
        item.symbol,
        item.count,
        item.mass,
        item.fraction,
    )
    assert composition.astuple()[0] == item.astuple()
    assert composition.asdict()['C'] == item.astuple()[1:]

    item = composition['12C']
    assert item.symbol == '12C'
    assert item.count == 1
    assert item.mass == 12.0
    assert item.fraction == pytest.approx(0.324201)
    assert composition.astuple()[1] == item.astuple()
    assert composition.asdict()['12C'] == item.astuple()[1:]

    composition = Formula('[12C][13C]C+').composition(False)
    assert '12C' not in str(composition)
    item = composition['e-']
    assert item.symbol == 'e-'
    assert item.count == -1
    item = composition['C']
    assert item.symbol == 'C'
    assert item.count == 3
    assert item.mass == pytest.approx(37.014095)
    assert item.fraction == pytest.approx(1.000014821057817)
    assert composition.astuple()[0] == item.astuple()
    assert composition.asdict()['C'] == item.astuple()[1:]


@pytest.mark.parametrize(
    'formula, mean, expected',
    [
        ('D', 2.0141018, ((2, 2.0141018, 1.0, 100.0, 2.0141018),)),
        ('D2', 4.0282036, ((4, 4.0282036, 1.0, 100.0, 4.0282036),)),
        (
            'H',
            1.0079408,
            (
                (1, 1.0078250, 0.999885, 100.0, 1.0078250),
                (2, 2.0141017, 0.000115, 0.0115013, 2.0141018),
            ),
        ),
        (
            'H+',
            1.0073922,
            (
                (1, 1.0072765, 0.999885, 100.0, 1.0072764),
                (2, 2.0135532, 0.000115, 0.0115013, 2.0135532),
            ),
        ),
        (
            'DH',
            3.0220425,
            (
                (3, 3.0219268, 0.999885, 100.0, 3.0219268),
                (4, 4.0282036, 0.000115, 0.0115013, 4.0282036),
            ),
        ),
        (
            'DHO',
            19.0214475,
            (
                (19, 19.0168414, 0.9974553, 100.0, 19.0168414),
                (20, 20.0215362, 0.00049468, 0.04959389, 20.0215362),
                (21, 21.0210866, 0.00204981, 0.20550374, 21.0210866),
                (22, 22.0273632, 2.35750000e-07, 2.36351448e-05, 22.0273632),
            ),
        ),
        (
            'SO4_2-',
            96.0030982,
            (
                (96, 95.9528268, 0.9407006, 100.0, 47.9764134),
                (97, 96.9529958, 0.0088607, 0.9419271, 48.4764979),
                (98, 97.9499357, 0.0498331, 5.2974427, 48.9749678),
            ),
        ),
    ],
)
def test_formula_spectrum(formula, mean, expected):
    """Test Formula.spectrum method."""
    min_intensity = 0.1 if len(formula) > 3 else 1e-16

    spectrum = Formula(formula).spectrum(min_intensity=min_intensity)
    assert isinstance(spectrum, Spectrum)
    repr(spectrum)
    str(spectrum)

    assert spectrum.mean == pytest.approx(mean, abs=1e-6)

    for (k, v), e in zip(spectrum.items(), expected):
        assert isinstance(v, SpectrumEntry)
        repr(v)
        str(v)
        assert k == e[0]
        assert v.massnumber == e[0]
        assert v.mass == pytest.approx(e[1], abs=1e-6)
        assert v.fraction == pytest.approx(e[2], abs=1e-6)
        assert v.intensity == pytest.approx(e[3], abs=1e-6)
        assert v.mz == pytest.approx(e[4], abs=1e-6)

    for (k, v), e in zip(spectrum.asdict().items(), expected):
        assert k == e[0]
        assert v == pytest.approx(e[1:], abs=1e-6)

    for a, b in zip(spectrum.astuple(), expected):
        assert a == pytest.approx(b, abs=1e-6)

    expected_peak = SpectrumEntry(*expected[0]).astuple()
    assert spectrum.peak.astuple() == pytest.approx(expected_peak, abs=1e-6)

    assert spectrum.range == (expected[0][0], expected[-1][0])

    spectrum.dataframe()


def test_formula_mul():
    """Test Formula multiplication."""
    formula = Formula('HO-') * 2
    assert formula.expanded == '[(HO)2]2-'


def test_formula_rmul():
    """Test Formula right multiplication."""
    formula = 2 * Formula('HO-')
    assert formula.expanded == '[(HO)2]2-'


def test_formula_add():
    """Test Formula addition."""
    formula = Formula('H2O') + Formula('HO-')
    assert formula.expanded == '[(H2O)(HO)]-'

    formula = Formula('H2O') + Formula('')
    assert formula.expanded == '(H2O)()'

    with pytest.raises(TypeError):
        Formula('H2O') + 'O'


def test_formula_sub():
    """Test Formula subtraction."""
    formula = Formula('H2O') - Formula('O')
    assert formula.expanded == 'H2'

    formula = Formula('H2O') - Formula('H2O')
    assert formula.expanded == ''

    formula = Formula('H2O') - Formula('')
    assert formula.expanded == 'H2O'

    with pytest.raises(ValueError):
        Formula('H2O') - Formula('O2')

    with pytest.raises(TypeError):
        Formula('H2O') - 'O'


@pytest.mark.parametrize(
    'numbers, expected',
    [
        ([], 1),
        ([4], 4),
        ([3, 6], 3),
        ([6, 7], 1),
    ],
)
def test_gcd(numbers, expected):
    """Test gcd function."""
    assert gcd(numbers) == expected


@pytest.mark.parametrize(
    'value, digits, expected',
    [
        (0.0, 5, 3),
        (-0.12345678, 5, 2),
        (1.23456789, 5, 3),
        (12.3456789, 5, 2),
        (12345.6789, 5, 1),
    ],
)
def test_precision_digits(value, digits, expected):
    """Test precision_digits function."""
    assert precision_digits(value, digits) == expected


@pytest.mark.parametrize(
    'value, charge, expected',
    [
        (1.0, -2, 0.5),
        (1.0, 1, 1.0),
        (1.0, 0, 1.0),
    ],
)
def test_mass_charge_ratio(value, charge, expected):
    """Test mass_charge_ratio function."""
    assert mass_charge_ratio(value, charge) == expected


@pytest.mark.parametrize(
    'formula, expected',
    [
        ('Formula', ('Formula', 0)),
        ('Formula+', ('Formula', 1)),
        ('Formula+1', ('Formula', 1)),
        ('Formula++', ('Formula', 2)),
        ('[Formula]2+', ('Formula', 2)),
        ('Formula+2', ('Formula', 2)),
        ('[[Formula]]2-', ('[Formula]', -2)),
        ('[Formula]-2', ('Formula', -2)),
        ('[Formula]_2-', ('Formula', -2)),
        ('Formula_2-', ('Formula', -2)),
        ('Formula-2', ('Formula', -2)),
        ('Formula_+', ('Formula', 1)),
        ('Formula+-', ('Formula', 0)),
        ('Formul+a+', ('Formul+a', 1)),
    ],
)
def test_split_charge(formula, expected):
    """Test split_charge function."""
    assert split_charge(formula) == expected


@pytest.mark.parametrize(
    'formula, kwargs, expected',
    [
        ('Valohp', {}, '(C5H8NO2)'),
        ('Valohp', {'parse_groups': False}, 'Valohp'),
        ('HLeu2OH', {}, 'H(C6H11NO)2OH'),
        ('D2O', {}, '[2H]2O'),
        ('PhNH2.HCl', {}, '(C6H5)NH2HCl'),
        ('PhNH2.HCl', {'parse_arithmetic': False}, '(C6H5)NH2.HCl'),
        ('CuSO4.5H2O', {}, 'CuSO4(H2O)5'),
        ('CuSO4+5*H2O', {}, 'CuSO4(H2O)5'),
        ('5*H2O', {}, '(H2O)5'),
        ('O2+12C', {}, 'O2(C)12'),
        ('ATC', {}, '((C10H12N5O5P)(C9H12N3O6P)(C10H13N2O7P)H2O)'),
        ('AUC', {}, '((C10H12N5O6P)(C9H12N3O7P)(C9H11N2O8P)H2O)'),
        ('ssdna(AC)', {}, '((C10H12N5O5P)(C9H12N3O6P)H2O)'),
        ('peptide(GG)', {}, '((C2H3NO)2H2O)'),
        ('WQ', {}, '((C5H8N2O2)(C11H10N2O)H2O)'),
        ('WQ', {'parse_oligos': False}, 'WQ'),
        ('O: 0.26, 30Si: 0.74', {}, 'O2[30Si]3'),
        (
            'O: 0.26, 30Si: 0.74',
            {'parse_fractions': False},
            'O:0(,30Si:0)26()74',
        ),
    ],
)
def test_from_string(formula, kwargs, expected):
    """Test from_string function."""
    assert from_string(formula, **kwargs) == expected


def test_from_string_error():
    """Test from_string function error."""
    with pytest.raises(TypeError):
        from_string(1)


@pytest.mark.parametrize(
    'formula, kwargs, expected',
    [
        ({}, {}, ''),
        (
            {'C': {0: 4, 12: 2}},
            {
                'divisor': 2,
                'charge': 2,
                'fmt': (
                    '{}',
                    '{}<sub>{}</sub>',
                    '<sup>{}</sup>{}',
                    '<sup>{}</sup>{}<sub>{}</sub>',
                ),
            },
            '[C<sub>2</sub><sup>12</sup>C]+',
        ),
        (
            {'C': {0: 4, 12: 2}},
            {
                'charge': 2,
                'fmt': (
                    '{}',
                    '{}<sub>{}</sub>',
                    '<sup>{}</sup>{}',
                    '<sup>{}</sup>{}<sub>{}</sub>',
                ),
            },
            '[C<sub>4</sub><sup>12</sup>C<sub>2</sub>]2+',
        ),
    ],
)
def test_from_elements(formula, kwargs, expected):
    """Test from_elements function."""
    assert from_elements(formula, **kwargs) == expected


@pytest.mark.parametrize(
    'formula, kwargs, expected',
    [
        ({}, {}, ''),
        ({'H': 0.5}, {'maxcount': 2, 'precision': 1.0}, 'H'),
        ({'H': 0.112, 'O': 0.888}, {}, 'H2O'),
        ({'D': 0.2, 'O': 0.8}, {}, 'O[2H]2'),
        ({'[2H]': 0.2, 'O': 0.8}, {}, 'O[2H]2'),
        ({'H': 8.97, 'C': 59.39, 'O': 31.64}, {}, 'C5H9O2'),
        ({'O': 0.26, '30Si': 0.74}, {}, 'O2[30Si]3'),
    ],
)
def test_from_fractions(formula, kwargs, expected):
    """Test from_fractions function."""
    assert from_fractions(formula, **kwargs) == expected


def test_from_fractions_error():
    """Test from_fractions function."""
    with pytest.raises(FormulaError):
        from_fractions({'[12O]': 1.0})


@pytest.mark.parametrize(
    'sequence, groups, expected',
    [
        ('A', {'A': 'B'}, '(B)'),
        ('AA', {'A': 'B'}, '(B)2'),
    ],
)
def test_from_sequence(sequence, groups, expected):
    """Test from_sequence function."""
    assert from_sequence(sequence, groups) == expected


@pytest.mark.parametrize(
    'sequence, expected',
    [
        ('GG', '((C2H3NO)2H2O)'),
        (
            'GPAVL IMCFY WHKRQ NEDST_2+',
            '[((C3H5NO)(C3H5NOS)(C4H5NO3)(C5H7NO3)(C9H9NO)(C2H3NO)(C6H7N3O)'
            '(C6H11NO)(C6H12N2O)(C6H11NO)(C5H9NOS)(C4H6N2O2)(C5H7NO)(C5H8N2O2)'
            '(C6H12N4O)(C3H5NO2)(C4H7NO2)(C5H9NO)(C11H10N2O)(C9H9NO2)H2O)]2+',
        ),
    ],
)
def test_from_peptide(sequence, expected):
    """Test from_peptide function."""
    assert from_peptide(sequence) == expected


@pytest.mark.parametrize(
    'sequence, dtype, expected',
    [
        ('AC', 'ssdna', '((C10H12N5O5P)(C9H12N3O6P)H2O)'),
        ('AU', 'dsrna', '((C10H12N5O6P)2(C9H11N2O8P)2(H2O)2)'),
        (
            'ATC G',
            'dsdna',
            '((C10H12N5O5P)2(C9H12N3O6P)2(C10H12N5O6P)2(C10H13N2O7P)2(H2O)2)',
        ),
        (
            'AUC G_2+',
            'ssrna',
            '[((C10H12N5O6P)(C9H12N3O7P)(C10H12N5O7P)(C9H11N2O8P)H2O)]2+',
        ),
    ],
)
def test_from_oligo(sequence, dtype, expected):
    """Test from_oligo function."""
    assert from_oligo(sequence, dtype) == expected


@pytest.mark.parametrize(
    'formula, expected',
    [
        (('H', 'C', 'O'), ('C', 'H', 'O')),
        (('O', 'H'), ('H', 'O')),
        (('Na', 'Cl'), ('Cl', 'Na')),
    ],
)
def test_hill_sorted(formula, expected):
    """Test hill_sorted function."""
    assert tuple(hill_sorted(formula)) == expected


@pytest.mark.parametrize(
    'args, expected',
    [
        ((0,), 'Formula'),
        ((1,), '[Formula]+'),
        ((2,), '[Formula]2+'),
        ((-2, '_'), 'Formula_2-'),
    ],
)
def test_join_charge(args, expected):
    """Test join_charge function."""
    assert join_charge('Formula', *args) == expected


@pytest.mark.parametrize(
    'args, expected',
    [
        ((0,), '0'),
        ((1,), '+'),
        ((-1,), '-'),
        ((2,), '2+'),
        ((-2,), '2-'),
        ((0, '_'), '0'),
        ((1, '_'), '+'),
        ((-1, '_'), '-'),
        ((2, '_'), '_2+'),
        ((-2, '_'), '_2-'),
    ],
)
def test_format_charge(args, expected):
    """Test format_charge function."""
    assert format_charge(*args) == expected


def test_elements():
    """Test element properties."""
    for e in ELEMENTS:
        assert isinstance(e, Element)
        e.validate()
        eval(repr(e))


def test_analyze():
    """Test analyze function."""
    result = analyze('C8H10N4O2', min_intensity=0.01)
    assert 'Formula: C8H10N4O2' in result
    assert 'Hill notation' not in result
    assert 'Empirical formula: C4H5N2O' in result
    assert 'O            2       31.99881     16.4780' in result
    assert 'Mass Distribution' in result
    assert '197      197.08721    0.050048     0.055681' in result

    result = analyze('H10N4O2C8', min_intensity=0.01, maxatoms=10)
    assert 'Formula: H10N4O2C8' in result
    assert 'Hill notation: C8H10N4O2' in result
    assert 'Empirical formula: C4H5N2O' in result
    assert 'Mass Distribution' not in result

    result = analyze('', debug=True)
    assert 'Nominal mass: 0' in result
    assert 'Average mass: 0.000000' in result

    result = analyze('Xy')
    assert 'Error: ' in result


def test_main(capsys):
    """Test main function."""
    assert main(['C8H10N4O2', '-v']) == 0
    result = capsys.readouterr().out
    assert 'Formula: C8H10N4O2' in result
    assert 'Hill notation' not in result
    assert 'Empirical formula: C4H5N2O' in result
    assert 'O            2       31.99881     16.4780' in result
    assert 'Mass Distribution' in result
    assert '197      197.08721    0.050048     0.055681' in result

    assert main(['Xy']) == 0
    result = capsys.readouterr().out
    assert 'Error: unknown symbol' in result

    try:
        assert main(['--help']) == 0
    except SystemExit:
        pass
    result = capsys.readouterr().out
    assert 'Usage: molmass [options] formula' in result


if __name__ == '__main__':
    import sys
    import warnings

    # warnings.simplefilter('always')
    warnings.filterwarnings('ignore', category=ImportWarning)
    argv = sys.argv
    argv.append('--cov-report=html')
    argv.append('--cov=molmass')
    argv.append('--verbose')
    sys.exit(pytest.main(argv))

# mypy: allow-untyped-defs
# mypy: check-untyped-defs=False
