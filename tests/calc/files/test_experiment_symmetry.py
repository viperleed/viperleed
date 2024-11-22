"""Tests for module viperleed.calc.files.experiment_symmetry."""

__authors__ = (
    'Alexander M. Imre (@amimre)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2024-07-01'
__license__ = 'GPLv3+'

import os

import pytest
from pytest_cases import parametrize_with_cases
from viperleed.calc.classes.slab import MissingBulkSlabError
from viperleed.calc.files import experiment_symmetry

from pathlib import Path

from .. import poscar_slabs
from ..tags import CaseTag as Tag

CasePOSCARSlabs = poscar_slabs.CasePOSCARSlabs

@parametrize_with_cases('args', cases=CasePOSCARSlabs,
                        has_tag=Tag.BULK_PROPERTIES)
def test_write_experiment_symmetry(args, tmp_path):
    slab, rpars, info = args

    os.chdir(tmp_path)
    path = Path.cwd()

    # create bulkslab
    slab.make_bulk_slab(rpars)

    # put some fake THEO_ENERGIES
    rpars.THEO_ENERGIES = rpars.THEO_ENERGIES.from_value([10, 100, 2])

    experiment_symmetry.write(slab, rpars)
    assert (path / 'experiment_symmetry.ini').exists()

    content = (path / 'experiment_symmetry.ini').read_text()
    assert 'superlattice' in content
    assert 'eMax = 100.00' in content


@parametrize_with_cases('args', cases=CasePOSCARSlabs,
                        has_tag=Tag.NO_INFO)
def test_write_fails_without_bulk_info(args, tmp_path):
    slab, rpars, info = args

    os.chdir(tmp_path)

    # put some fake THEO_ENERGIES
    rpars.THEO_ENERGIES = rpars.THEO_ENERGIES.from_value([10, 100, 2])

    with pytest.raises(MissingBulkSlabError):
        experiment_symmetry.write(slab, rpars)

@parametrize_with_cases('args', cases=CasePOSCARSlabs,
                        has_tag=Tag.BULK_PROPERTIES)
def test_write_failing(args, monkeypatch):
    slab, rpars, info = args

    # create bulkslab
    slab.make_bulk_slab(rpars)

    # Monkey patch the `open` function to raise OSError
    def mock_open(*args, **kwargs):
        raise OSError("Mocked OSError")

    # put some fake THEO_ENERGIES
    rpars.THEO_ENERGIES = rpars.THEO_ENERGIES.from_value([10, 100, 2])

    monkeypatch.setattr('builtins.open', mock_open)
    with pytest.raises(OSError):
        experiment_symmetry.write(slab, rpars)
