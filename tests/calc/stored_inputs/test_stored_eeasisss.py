"""Tests of the stored eeasisss outputs and of how they are checked.

The run_phaseshift fixture replaces eeasisss by stored raw outputs for
the slabs of MOCKED_CASES (see stored_eeasisss.py). The release check
writes the stored data of each of these slabs anew with the real
eeasisss, as by the regeneration command, and compares it with the
stored one. This takes long (Al2O3/NiAl(111) about 13 min on one
core), so it only runs with --release. The other tests, in the
default tier, check the tree in which the release check runs
eeasisss, and that the tier tables list the mocked slabs.
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
