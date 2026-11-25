"""Tests of the stored eeasisss outputs and of how they are checked.

The run_phaseshift fixture replaces eeasisss by stored raw outputs for
the slabs of MOCKED_CASES (see stored_eeasisss.py). The release check
writes the stored data of each of these slabs anew with the real
eeasisss, as by the regeneration command, and compares it with the
stored one. This takes long (Al2O3/NiAl(111) about 13 min on one
core), so it only runs with --release. The other tests, in the
default tier, check the tolerance of the comparison, the tree in
which the release check runs eeasisss, and that the tier tables list
the mocked slabs.
"""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-11-25'
__license__ = 'GPLv3+'

from types import SimpleNamespace

import pytest
from pytest_cases import parametrize

from ...harness import tiers
from . import stored_eeasisss

_FIXTURE = 'run_phaseshift'  # The fixture that uses the stored data
_WINDOWS_EXECUTABLE = 'eeasisss.exe'  # What psgen runs on Windows

# Excerpts of the formats of PS.r.* files and of PHASESHIFTS
_STORED = """\
 17  -9.70   -2.99 PS.r.14.100   stored-input 260101-000000
    20.0 -0.1947  0.2808  1.2E-05  0.0000
-0.6222-0.1719-0.1963 0.0150
"""
_WITHIN_TOLERANCE = {
    'identical': _STORED,
    'other whitespace': _STORED.replace('    20.0', ' 20.0'),
    'one unit up': _STORED.replace('0.2808', '0.2809'),
    'one unit down': _STORED.replace('0.2808', '0.2807'),
    'one unit in packed fields': _STORED.replace('-0.1719', '-0.1720'),
    'one unit with exponent': _STORED.replace('1.2E-05', '1.3E-05'),
    'Fortran exponent': _STORED.replace('1.2E-05', '1.2D-05'),
    'sign of zero': _STORED.replace(' 0.0000', '-0.0000'),
    }
_BEYOND_TOLERANCE = {
    'two units': _STORED.replace('0.2808', '0.2810'),
    'few units in packed fields': _STORED.replace('-0.1719', '-0.1722'),
    'few units with exponent': _STORED.replace('1.2E-05', '1.5E-05'),
    'other name': _STORED.replace('PS.r.14.100', 'PS.r.14.101'),
    'other integer': _STORED.replace(' 17 ', ' 18 '),
    'number for name': _STORED.replace('stored-input', '0.1234'),
    'token missing': _STORED.replace(' 0.0150', ''),
    'token added': _STORED.replace('0.0150', '0.0150 0.0000'),
    'packed field missing': _STORED.replace('-0.1719-0.1963', '-0.1719'),
    'fewer decimals': _STORED.replace('0.2808', '0.28'),
    'more decimals': _STORED.replace('0.0150', '0.01509'),
    }


@pytest.mark.release
@parametrize('case_id', stored_eeasisss.MOCKED_CASES)
def test_stored_data_up_to_date(case_id, tensorleed_path, tmp_path):
    """Check that the real eeasisss gives the stored data of a slab."""
    stored_eeasisss.check_stored_data(case_id, tensorleed_path, tmp_path)


def test_real_tree_on_windows(tmp_path, mocker):
    """Check that the copied tree has the executable psgen runs."""
    mocker.patch.object(stored_eeasisss, 'os',
                        SimpleNamespace(name=stored_eeasisss.WINDOWS))
    real_parent = tmp_path / 'real'
    real_parent.mkdir()
    real_tree = stored_eeasisss.make_stub_tree(real_parent)
    tree = stored_eeasisss.copy_real_tree(tmp_path, real_tree)
    assert (tree / _WINDOWS_EXECUTABLE).is_file()


def test_tier_tables_list_mocked_cases():
    """Check that tiers.py has the mocked and real ids of the slabs."""
    case_ids = set(stored_eeasisss.MOCKED_CASES)
    mocked_ids = {case_id + stored_eeasisss.MOCKED_SUFFIX
                  for case_id in case_ids}
    assert set(tiers.MOCKED_CASES[_FIXTURE]) == mocked_ids
    assert set(tiers.RELEASE_CASES[_FIXTURE]) == case_ids


@parametrize('written',
             _WITHIN_TOLERANCE.values(),
             ids=tuple(_WITHIN_TOLERANCE))
def test_mismatch_within_tolerance(written):
    """Check that numbers off by one unit in the last digit match."""
    assert stored_eeasisss.first_mismatch(_STORED, written) is None


@parametrize('written',
             _BEYOND_TOLERANCE.values(),
             ids=tuple(_BEYOND_TOLERANCE))
def test_mismatch_beyond_tolerance(written):
    """Check that larger or non-numeric differences are reported."""
    assert stored_eeasisss.first_mismatch(_STORED, written)
