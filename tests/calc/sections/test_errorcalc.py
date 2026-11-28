"""Tests for module errorcalc of viperleed.calc.sections."""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-10-08'
__license__ = 'GPLv3+'

import logging

import numpy as np
import pytest
from pytest_cases import parametrize

from viperleed.calc.classes.rparams.special.r_factor_type import RFactorType
from viperleed.calc.lib.version import Version
from viperleed.calc.sections import errorcalc

_MODULE = 'viperleed.calc.sections.errorcalc'
_V0I = 4.5
_ENERGY_RANGE = 500.0
_MODE = 'geo'


# R_FACTOR_TYPE, R_FACTOR_LEGACY, TL_VERSION, and the R factor computed
_COMPUTED_R_FACTOR = (
    pytest.param('pendry', True, '2.1.0', 'pendry', id='pendry'),
    pytest.param('smooth', True, '2.1.0', 'smooth', id='smooth'),
    pytest.param('smooth', False, '2.1.0', 'smooth', id='smooth-python'),
    pytest.param('smooth', True, '2.0.0', 'pendry', id='smooth-tl2.0'),
    pytest.param('smooth', False, '2.0.0', 'smooth',
                 id='smooth-tl2.0-python'),
    pytest.param('r2', True, '2.0.0', 'r2', id='r2-tl2.0'),
    )


@parametrize('r_factor_type,legacy,tl_version,expected', _COMPUTED_R_FACTOR)
# pylint: disable-next=too-many-arguments,too-many-positional-arguments
def test_figures_get_r_factor_type(r_factor_type, legacy, tl_version,
                                   expected, mocker, tmp_path, monkeypatch):
    """Check that error plots are labelled with the R computed."""
    monkeypatch.chdir(tmp_path)
    # Fake Rparams with one atom varied once per mode
    atom = mocker.MagicMock()
    search_par = mocker.MagicMock(atom=atom, steps=2, el='vac')
    rpars = mocker.MagicMock(
        domainParams=[],
        best_v0r=0.,
        search_index=0,
        disp_blocks=[None],
        STOP=False,
        halt=0,
        HALTING=2,
        search_atlist=[atom],
        searchpars=[search_par],
        R_FACTOR_TYPE=RFactorType(r_factor_type),
        R_FACTOR_LEGACY=legacy,
        TL_VERSION=Version(tl_version),
        )
    for name in ('R_Error',
                 'section_deltas',
                 'section_superpos'):
        mocker.patch(f'{_MODULE}.{name}')
    mocker.patch(f'{_MODULE}.readDISPLACEMENTS_block', return_value=True)
    mocker.patch(f'{_MODULE}.section_rfactor', return_value=[0.1])
    io_mock = mocker.patch(f'{_MODULE}.ioerrorcalc')
    io_mock.generate_errors_csv.return_value = ('', {})

    errorcalc.errorcalc(mocker.MagicMock(), rpars)

    make_figs = io_mock.make_errors_figs
    make_figs.assert_called_once()
    _, kwargs = make_figs.call_args  # .kwargs needs Python 3.8
    assert kwargs['r_factor_type'] is RFactorType(expected)


@parametrize(r_factor_type=('pendry', 'smooth'))
def test_var_r_reaches_output(r_factor_type, mocker, tmp_path,
                              monkeypatch, caplog):
    """Check that var(R) from a real R_Error reaches log and outputs."""
    monkeypatch.chdir(tmp_path)
    caplog.set_level(logging.INFO)
    # Fake Rparams with one atom varied in five geo steps
    displacements = [-0.2, -0.1, 0.0, 0.1, 0.2]
    rfactors = [0.5, 0.3, 0.2, 0.25, 0.45]
    atom = mocker.MagicMock(num=1, el='O', disp_geo={'all': displacements})
    atom.displist[0].disp_labels = {_MODE: 'z'}
    atom.displist[0].disp_lin_steps = {_MODE: displacements}
    search_par = mocker.MagicMock(atom=atom, steps=5, el='O', mode=_MODE)
    rpars = mocker.MagicMock(
        domainParams=[],
        best_v0r=0.,
        search_index=0,
        disp_blocks=[None],
        STOP=False,
        halt=0,
        HALTING=2,
        search_atlist=[atom],
        searchpars=[search_par],
        R_FACTOR_TYPE=RFactorType(r_factor_type),
        R_FACTOR_LEGACY=False,
        V0_IMAG=_V0I,
        )
    rpars.total_energy_range.return_value = _ENERGY_RANGE
    for name in ('section_deltas',
                 'section_superpos'):
        mocker.patch(f'{_MODULE}.{name}')
    mocker.patch(f'{_MODULE}.readDISPLACEMENTS_block',
                 side_effect=lambda *_, only_mode: only_mode == _MODE)
    mocker.patch(f'{_MODULE}.section_rfactor', return_value=rfactors)
    generate_csv = mocker.patch(f'{_MODULE}.ioerrorcalc.generate_errors_csv',
                                return_value=('', {}))
    for name in ('write_errors_archive',
                 'write_errors_summary_csv',
                 'make_errors_figs',
                 'write_errors_pdf'):
        mocker.patch(f'{_MODULE}.ioerrorcalc.{name}')

    errorcalc.errorcalc(mocker.MagicMock(), rpars)

    var_r = np.sqrt(8 * _V0I / _ENERGY_RANGE) * min(rfactors)
    (errors,), _ = generate_csv.call_args  # .args needs Python 3.8
    assert errors[0].r_type is RFactorType(r_factor_type)
    assert errors[0].var_r == pytest.approx(var_r)
    assert f'Found values for var(R): {_MODE}: {var_r:.4f}' in caplog.text
