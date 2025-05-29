"""Tests for module viperleed.calc.files.experiment_symmetry."""

__authors__ = (
    'Alexander M. Imre (@amimre)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2024-07-01'
__license__ = 'GPLv3+'

import re

import pytest
from pytest_cases import fixture
from pytest_cases import parametrize_with_cases

from viperleed.calc import symmetry
from viperleed.calc.classes.rparams.special.energy_range import TheoEnergies
from viperleed.calc.classes.slab import MissingBulkSlabError
from viperleed.calc.files import experiment_symmetry
from viperleed.calc.lib.context import execute_in_dir

from ...helpers import exclude_tags
from ..poscar_slabs import CasePOSCARSlabs
from ..symmetry.test_symmetry import TestPlaneGroupFinding as _TestPG
from ..symmetry.test_symmetry import _reconstruct_case_id
from ..tags import CaseTag as Tag

_MODULE = 'viperleed.calc.files.experiment_symmetry'
SYMM_FAILS = _TestPG.known_incorrect_groups
NO_GROUP = 'unknown'
_INT_RE = r'-?\d+'
_FLOAT_RE = r'-?\d+(.\d*)?'
_MATRIX_RE = r'\[\[{0},\s*{0}\],\s*\[{0},\s*{0}\]\]'


def mock_energies(rpars):
    """Give an example THEO_ENERGIES to rpars."""
    rpars.THEO_ENERGIES = TheoEnergies(10, 100, 2)


def prepare_slab_and_rpars(slab, rpars, make_bulk=True):
    """Prepare slab the same way as initialization would."""
    mock_energies(rpars)
    if slab.foundplanegroup == NO_GROUP:
        symmetry.findSymmetry(slab, rpars, forceFindOri=True)
    if make_bulk:
        bulk = slab.make_bulk_slab(rpars)
        symmetry.findBulkSymmetry(bulk, rpars)
        symmetry.findSymmetry(bulk, rpars, bulk=True)


@fixture(name='with_info')
@parametrize_with_cases('args',
                        cases=CasePOSCARSlabs,
                        filter=exclude_tags(Tag.NO_INFO))
def fixture_with_info(args):
    """Return a prepared slab and rpars."""
    slab, rpars, *_ = args
    prepare_slab_and_rpars(slab, rpars)
    return args


def test_success(with_info, tmp_path, subtests, first_case):
    """Check the expected result of writing experiment_symmetry.ini."""
    slab, rpars, info = with_info
    with execute_in_dir(tmp_path):
        experiment_symmetry.write(slab, rpars)
    written = tmp_path / 'experiment_symmetry.ini'
    assert written.exists()
    contents = [line.strip() for line in written.read_text().splitlines()]

    this_case = _reconstruct_case_id(first_case)
    expected = (
        r'\[\w+\]',  # Header
        rf'SUPERLATTICE = {_MATRIX_RE.format(_INT_RE)}',
        rf'surfBasis = {_MATRIX_RE.format(_FLOAT_RE)}',
        'eMax = 100(.00)?',  # GUI doesn't care if it's float
        rf'screenAperture = {_FLOAT_RE}',
        rf'beamIncidence = \({_FLOAT_RE}, {_FLOAT_RE}\)',
        )
    if info.symmetry.hermann and this_case not in SYMM_FAILS:
        expected += (f'surfGroup = {info.symmetry.hermann}.*',)
    for pattern in expected:
        with subtests.test(pattern):
            assert any(re.match(pattern, line) for line in contents)


def test_raises_without_bulk_info(mocker):
    """Check complaints if no bulk-symmetry information is present."""
    slab, rpars = (mocker.MagicMock() for _ in range(2))
    slab.bulkslab = None
    with pytest.raises(MissingBulkSlabError):
        experiment_symmetry.write(slab, rpars)


def test_write_failing(mocker):
    """Check complaints when writing to file fails."""
    mock_dict = {'SUPERLATTICE': mocker.MagicMock(),
                 'surfBasis': mocker.MagicMock()}
    mocker.patch(f'{_MODULE}.getLEEDdict', return_value=mock_dict)
    mocker.patch('builtins.open', side_effect=OSError)
    args = (mocker.MagicMock() for _ in range(2))
    with pytest.raises(OSError):
        experiment_symmetry.write(*args)


def test_get_dict_fails(mocker):
    """Check complaints when fetching a dict of LEED parameters fails."""
    mocker.patch(f'{_MODULE}.getLEEDdict', return_value=None)
    args = (mocker.MagicMock() for _ in range(2))
    with pytest.raises(ValueError):
        experiment_symmetry.write(*args)
