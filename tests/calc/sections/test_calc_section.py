"""Tests for module calc_section of viperleed.calc.sections."""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-11-24'
__license__ = 'GPLv3+'

import pytest
from pytest_cases import parametrize

from viperleed.calc.sections.calc_section import CalcSection as Section

_ALIASES = {
    Section.INITIALIZATION: ('0', 'ini', 'init', 'initialize',
                             'initialise', 'initialization',
                             'initialisation', 'INITIALIZATION'),
    Section.REFCALC: ('1', 'ref', 'fd', 'refcalc', 'ref_calc',
                      'reference', 'reference_calculation', 'ref-calc',
                      'fd_calc', 'full_dynamic', 'fulldynamic_calculation'),
    Section.RFACTOR_REFCALC: ('11', 'rfac', 'r_fac', 'rfactor', 'r_factor',
                              'RFACTOR_REFCALC'),
    Section.RFACTOR_SUPERPOS: ('12', 'rfacsuper', 'rfac_super',
                               'r_factor_superpos', 'RFACTOR_SUPERPOS'),
    Section.DELTAS: ('2', 'del', 'delta', 'deltas', 'delta_calculation',
                     'deltascalc', 'deltaamplitudes', 'delta_amplitudes',
                     'delta-amplitudes', 'DELTA-AMPLITUDES'),
    Section.SEARCH: ('3', 'search', 'structure_search', 'struct_search',
                     'structsearch'),
    Section.SUPERPOS: ('31', 'sup', 'super', 'superpos'),
    Section.DOMAINS: ('4', 'dom', 'domain', 'domains',
                      'domain_calculation', 'domainscalc', 'domains-calc'),
    Section.ERRORCALC: ('5', 'err', 'error', 'errors', 'errorcalc',
                        'error_calculation', 'errors_calc',
                        'error-calculation'),
    Section.FD_OPTIMIZATION: ('6', 'opt', 'optimize', 'optimisation',
                              'fdopt', 'fd_opt', 'fd_optimization',
                              'full_dynamic_optimization',
                              'fulldynamic_optimise', 'fd-opt',
                              'fd-optimization'),
    }


@parametrize('section,alias', ((section, alias)
                               for section, aliases in _ALIASES.items()
                               for alias in aliases))
def test_alias(section, alias):
    """Check that `alias` is interpreted as `section`."""
    assert Section(alias) is section


@parametrize(alias=('Ref', 'RefCalc', 'SuperPos', 'RFAC_Super', 'FD_Opt'))
def test_alias_case_insensitive(alias):
    """Check that aliases are not case sensitive."""
    assert Section(alias) is Section(alias.lower())


_STARTS_WITH_ALIAS = {
    'superposition': Section.SUPERPOS,
    'super_pos': Section.SUPERPOS,
    'rfactors': Section.RFACTOR_REFCALC,
    'rfac_refcalc': Section.RFACTOR_REFCALC,
    'rfactor_superposition': Section.RFACTOR_SUPERPOS,
    'refcalcx': Section.REFCALC,
    'fd_': Section.REFCALC,
    'fd_optimizer': Section.FD_OPTIMIZATION,
    'optim': Section.FD_OPTIMIZATION,
    'initt': Section.INITIALIZATION,
    'deltaz': Section.DELTAS,
    'Superposition': Section.SUPERPOS,
    }


@parametrize('name,section', _STARTS_WITH_ALIAS.items())
def test_starts_with_alias(name, section):
    """Check that names starting with an alias give the longest one."""
    assert Section(name) is section


_INVALID = ('unknown', 'xref', 'ref calc', 'ref-calcx', 'delta-', '',
            'ref:search', 'search:superpos', 'super:pos')


@parametrize(alias=_INVALID)
def test_invalid_alias(alias):
    """Check that an unknown section name is not accepted."""
    with pytest.raises(ValueError, match='is not a valid CalcSection'):
        Section(alias)


@parametrize('string_,expect', (
    ('ref-search', (Section.REFCALC, Section.DELTAS, Section.SEARCH)),
    ('init-ref', (Section.INITIALIZATION, Section.REFCALC)),
    ('delta-search', (Section.DELTAS, Section.SEARCH)),
    ('ref:search', (Section.REFCALC, Section.DELTAS, Section.SEARCH)),
    ))
def test_range_of_aliases(string_, expect):
    """Check that a hyphen between aliases gives a range of sections."""
    assert Section.sequence_from_string(string_) == expect


@parametrize(string_=('search:superpos', 'ref-calc:search', 'fd-opt-x'))
def test_invalid_range(string_):
    """Check that a string that is no section or range is refused."""
    with pytest.raises(ValueError, match='CalcSection'):
        Section.sequence_from_string(string_)
