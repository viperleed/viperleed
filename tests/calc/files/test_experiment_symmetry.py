"""Tests for module viperleed.calc.files.experiment_symmetry."""

__authors__ = (
    'Alexander M. Imre (@amimre)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2024-07-01'
__license__ = 'GPLv3+'

import pytest
from pytest_cases import fixture
from pytest_cases import parametrize_with_cases

from viperleed.calc.classes.rparams.special.energy_range import TheoEnergies
from viperleed.calc.classes.slab import MissingBulkSlabError
from viperleed.calc.files import experiment_symmetry
from viperleed.calc.lib.context import execute_in_dir

from pathlib import Path

from .. import poscar_slabs
from ..tags import CaseTag as Tag

CasePOSCARSlabs = poscar_slabs.CasePOSCARSlabs

with_bulk_info = parametrize_with_cases('args', cases=CasePOSCARSlabs,
                                        has_tag=Tag.BULK_PROPERTIES)


@fixture(name='mock_energies')
def fixture_mock_energies():
    """Give an example THEO_ENERGIES to rpars."""
    def _mock(rpars):
        rpars.THEO_ENERGIES = TheoEnergies(10, 100, 2)
    return _mock


@with_bulk_info
def test_success(args, mock_energies, tmp_path):
    """Check the expected result of writing experiment_symmetry.ini."""
    slab, rpars, *_ = args
    slab.make_bulk_slab(rpars)
    mock_energies(rpars)
    with execute_in_dir(tmp_path):
        experiment_symmetry.write(slab, rpars)
    written = tmp_path / 'experiment_symmetry.ini'
    assert written.exists()
    content = written.read_text()
    assert 'superlattice' in content
    assert 'eMax = 100.00' in content


@parametrize_with_cases('args', cases=CasePOSCARSlabs,
                        has_tag=Tag.NO_INFO)
def test_raises_without_bulk_info(args, mock_energies):
    """Check complaints if no bulk-symmetry information is present."""
    slab, rpars, *_ = args
    mock_energies(rpars)
    with pytest.raises(MissingBulkSlabError):
        experiment_symmetry.write(slab, rpars)


@with_bulk_info
def test_write_failing(args, mock_energies, mocker):
    """Check complaints when writing to file fails."""
    slab, rpars, *_ = args
    slab.make_bulk_slab(rpars)
    mock_energies(rpars)
    mocker.patch('builtins.open', side_effect=OSError)
    with pytest.raises(OSError):
        experiment_symmetry.write(slab, rpars)
