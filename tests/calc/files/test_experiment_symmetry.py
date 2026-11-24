"""Tests for module viperleed.calc.files.experiment_symmetry."""

__authors__ = (
    'Alexandra Mia Imre (@alexmiame)',
    'Michele Riva (@michele-riva)',
)
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2024-07-01'
__license__ = 'GPLv3+'

import re

import pytest
from pytest_cases import fixture
from pytest_cases import parametrize
from pytest_cases import parametrize_with_cases

from viperleed.calc import symmetry
from viperleed.calc.classes.rparams.special.energy_range import TheoEnergies
from viperleed.calc.classes.slab import MissingBulkSlabError
from viperleed.calc.files import experiment_symmetry
from viperleed.calc.files.parameters import ParameterInterpreter
from viperleed.calc.files.parameters.utils import Assignment
from viperleed.calc.lib.context import execute_in_dir
from viperleed.calc.lib.version import Version
from viperleed.calc.sections.initialization import _inherit_from_main

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


@fixture(name='with_info')
@parametrize_with_cases('args',
                        cases=CasePOSCARSlabs,
                        filter=exclude_tags(Tag.NO_INFO))
def fixture_with_info(args):
    """Return a prepared slab and rpars."""
    slab, rpars, *_ = args
    # Add some energies to check that eMax is set correctly
    rpars.THEO_ENERGIES = TheoEnergies(10, 100, 2)

    # Find plane groups of slab and bulk to ensure symmetry
    # info is written correctly.
    if slab.foundplanegroup == NO_GROUP:
        symmetry.findSymmetry(slab, rpars, forceFindOri=True)
    bulk = slab.make_bulk_slab(rpars)
    symmetry.findBulkSymmetry(bulk, rpars)
    symmetry.findSymmetry(bulk, rpars, bulk=True)
    return args


def _prepare_for_writing(slab, rpars):
    """Prepare `slab` and `rpars` for writing experiment_symmetry."""
    rpars.THEO_ENERGIES = TheoEnergies(10, 100, 2)
    symmetry.findSymmetry(slab, rpars, forceFindOri=True)
    bulk = slab.make_bulk_slab(rpars)
    symmetry.findBulkSymmetry(bulk, rpars)
    symmetry.findSymmetry(bulk, rpars, bulk=True)


def _written_beam_incidence(slab, rpars, tmp_path):
    """Return the beamIncidence line written for `slab` and `rpars`."""
    with execute_in_dir(tmp_path):
        experiment_symmetry.write(slab, rpars)
    written = (tmp_path / experiment_symmetry.FILENAME).read_text()
    return next(line for line in written.splitlines()
                if line.startswith('beamIncidence'))


_BEAM_INCIDENCE_LINES = {
    'default': (None, 'beamIncidence = (0, 0)'),
    'normal': ('THETA 0, PHI 0', 'beamIncidence = (0.0, 0.0)'),
    'off-normal': ('THETA 20.5, PHI 30', 'beamIncidence = (20.5, 30.0)'),
    'negative theta': ('THETA -20, PHI 30',
                       'beamIncidence = (20.0, 210.0)'),
    }


@fixture(name='ag100_main')
@parametrize('incidence,expect',
             _BEAM_INCIDENCE_LINES.values(),
             ids=_BEAM_INCIDENCE_LINES)
def fixture_ag100_main(incidence, expect):
    """Return a Ag(100) slab and an Rparams with BEAM_INCIDENCE read."""
    slab, rpars, *_ = CasePOSCARSlabs().case_poscar_ag100()
    if incidence is not None:
        interpreter = ParameterInterpreter(rpars)
        interpreter.interpret_beam_incidence(
            Assignment(incidence, 'BEAM_INCIDENCE')
            )
    return slab, rpars, expect


def test_beam_incidence_line(ag100_main, tmp_path):
    """Check the exact beamIncidence line written."""
    slab, rpars, expect = ag100_main
    _prepare_for_writing(slab, rpars)
    assert _written_beam_incidence(slab, rpars, tmp_path) == expect


def test_beam_incidence_line_domain(ag100_main, tmp_path, mocker):
    """Check the exact beamIncidence line written for a domain."""
    _, main_rpars, expect = ag100_main
    main_rpars.TL_VERSION = Version('2.0.0')
    slab, rpars, *_ = CasePOSCARSlabs().case_poscar_ag100()
    domain = mocker.MagicMock(rpars=rpars, workdir=tmp_path)
    mocker.patch('shutil.copy2')  # IVBEAMS
    mocker.patch('viperleed.calc.files.parameters.modify')
    with execute_in_dir(tmp_path):
        _inherit_from_main((domain,), main_rpars, {domain: ()})
    _prepare_for_writing(slab, rpars)
    assert _written_beam_incidence(slab, rpars, tmp_path) == expect


def test_get_dict_fails(mocker):
    """Check complaints when fetching a dict of LEED parameters fails."""
    mocker.patch(f'{_MODULE}.getLEEDdict', return_value=None)
    args = (mocker.MagicMock() for _ in range(2))
    with pytest.raises(ValueError):
        experiment_symmetry.write(*args)


def test_raises_without_bulk_info(mocker):
    """Check complaints if no bulk-symmetry information is present."""
    slab, rpars = (mocker.MagicMock() for _ in range(2))
    slab.bulkslab = None
    with pytest.raises(MissingBulkSlabError):
        experiment_symmetry.write(slab, rpars)


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


def test_write_fails(mocker):
    """Check complaints when writing to file fails."""
    mock_dict = {'SUPERLATTICE': mocker.MagicMock(),
                 'surfBasis': mocker.MagicMock()}
    mocker.patch(f'{_MODULE}.getLEEDdict', return_value=mock_dict)
    mocker.patch('builtins.open', side_effect=OSError)
    args = (mocker.MagicMock() for _ in range(2))
    with pytest.raises(OSError):
        experiment_symmetry.write(*args)
