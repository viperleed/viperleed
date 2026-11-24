"""Tests for module fd_optimization of viperleed.calc.sections."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2025-03-25'
__license__ = 'GPLv3+'

import logging
from types import SimpleNamespace

import pytest
from pytest_cases import fixture
from pytest_cases import parametrize

from viperleed.calc.classes.rparams.rparams import Rparams
from viperleed.calc.classes.rparams.special.beam_incidence import (
    BeamIncidence,
    )
from viperleed.calc.files.parameters import ParameterInterpreter
from viperleed.calc.files.parameters.utils import Assignment
from viperleed.calc.files.parameters.write import ModifiedParameterValue
from viperleed.calc.lib.context import execute_in_dir
from viperleed.calc.sections.fd_optimization import _replace_incidence_angle
from viperleed.calc.sections.fd_optimization import fd_optimization
from viperleed.calc.sections.fd_optimization import get_fd_r

_MODULE = 'viperleed.calc.sections.fd_optimization'


class TestReplaceIncidenceAngle:
    """Tests for the _replace_incidence_angle helper."""

    _valid = {
        'theta': ('theta', 20, (20.0, 50.0)),
        'theta negative': ('theta', -30, (30.0, 230.0)),
        'theta beyond grazing': ('theta', 95, (90.0, 50.0)),
        'theta beyond negative grazing': ('theta', -95, (90.0, 230.0)),
        'phi': ('phi', 70, (10.0, 70.0)),
        'phi above 360': ('phi', 370, (10.0, 10.0)),
        'phi 360': ('phi', 360, (10.0, 0.0)),
        'phi negative': ('phi', -20, (10.0, 340.0)),
        }

    @parametrize('which,value,expect', _valid.values(), ids=_valid)
    def test_replace(self, which, value, expect):
        """Check the incidence obtained by replacing one angle."""
        rpars = Rparams()
        rpars.BEAM_INCIDENCE = BeamIncidence(10, 50)
        incidence = _replace_incidence_angle(rpars, which, value)
        assert (incidence.theta, incidence.phi) == expect
        assert rpars.BEAM_INCIDENCE == BeamIncidence(10, 50)  # Unchanged

    def test_replace_theta_wraps_phi(self):
        """Check that a phi of 360 degrees becomes zero."""
        rpars = Rparams()
        rpars.BEAM_INCIDENCE = BeamIncidence(10, 360)
        incidence = _replace_incidence_angle(rpars, 'theta', 20)
        assert (incidence.theta, incidence.phi) == (20, 0)


class TestGetFdR:
    """Tests for the get_fd_r function."""

    @fixture(name='mock_args')
    def fixture_mock_args(self, mocker, tmp_path):
        """Return arguments for get_fd_r."""
        return (mocker.MagicMock(name='slab'),
                Rparams(),
                tmp_path/'fd_op_partial')

    @fixture(name='mocks')
    def fixture_mock_implementation(self, mocker):
        """Replace implementation details with mocks."""
        return {
            'refcalc': mocker.patch(f'{_MODULE}.section_refcalc'),
            'r_factor': mocker.patch(f'{_MODULE}.section_rfactor',
                                     return_value=mocker.MagicMock()),
            'copy': mocker.patch('shutil.copy2'),
            }

    def test_success(self, mock_args, mocks, tmp_path, caplog):
        """Check results of a successful call."""
        slab, rpars, workdir = mock_args
        with execute_in_dir(tmp_path), caplog.at_level(logging.WARNING):
            rfactor, rfactor_per_beam = get_fd_r(*mock_args)
        assert rfactor is rpars.last_R
        assert rfactor_per_beam is mocks['r_factor'].return_value
        assert rpars.TENSOR_OUTPUT == [0]  # All off
        assert workdir.is_dir()
        assert (workdir/rpars.paths.compile_logs.name).is_dir()

        mocks['refcalc'].assert_called_once_with(slab,
                                                 rpars,
                                                 parent_dir=tmp_path)
        mocks['r_factor'].assert_called_once_with(slab, rpars, 11)
        mocks['copy'].assert_any_call(tmp_path/'BEAMLIST', workdir)
        mocks['copy'].assert_any_call(tmp_path/'PHASESHIFTS', workdir)
        assert not caplog.text

    _phi = {
        'phi 360': ((10, 360), (10, 0)),
        'phi below 360': ((10, 50), (10, 50)),
        'normal incidence': ((0, 0), (0, 0)),
        }

    @parametrize('incidence,expect', _phi.values(), ids=_phi)
    # pylint: disable-next=too-many-arguments,too-many-positional-arguments
    def test_phi_wrapped(self, incidence, expect,
                         mock_args, mocks, tmp_path, caplog):
        """Check that the calculation runs with phi in [0, 360)."""
        _, rpars, _ = mock_args
        rpars.BEAM_INCIDENCE = BeamIncidence(*incidence)
        used_in_refcalc = []

        def _store_incidence(_, rpars, **__):
            incidence = rpars.BEAM_INCIDENCE
            used_in_refcalc.append((incidence.theta, incidence.phi))

        mocks['refcalc'].side_effect = _store_incidence
        self.test_success(mock_args, mocks, tmp_path, caplog)
        assert used_in_refcalc == [expect]

    _raises = {
        'copy': OSError,
        'refcalc': Exception,
        'r_factor': Exception,
        }

    @parametrize('callee,exc', _raises.items(), ids=_raises)
    # pylint: disable-next=too-many-arguments  # 4/7 are fixtures
    def test_raises(self, callee, exc, mock_args, mocks, tmp_path, caplog):
        """Check complaints when failing to execute a `callee`."""
        mocks[callee].side_effect = exc
        with pytest.raises(exc):
            with execute_in_dir(tmp_path):
                get_fd_r(*mock_args)
        assert caplog.text


class TestOptimizeIncidence:
    """Tests for fd_optimization of the angles of BEAM_INCIDENCE."""

    @fixture(name='optimize')
    def fixture_optimize(self, mocker, tmp_path):
        """Return a function that optimizes an incidence angle."""
        def _fake_r(_, rpars, workdir):
            angle = float(workdir.name.split('_')[-1])
            return (angle - rpars.fake_minimum)**2, []

        def _run(which, incidence, minimum):
            rpars = Rparams()
            ParameterInterpreter(rpars).interpret_beam_incidence(
                Assignment(incidence, 'BEAM_INCIDENCE')
                )
            rpars.OPTIMIZE['which'] = which
            rpars.FORTRAN_COMP = ['gfortran', '']
            rpars.N_CORES = 1
            rpars.fake_minimum = minimum
            mocker.patch(f'{_MODULE}.get_fd_r', side_effect=_fake_r)
            mocker.patch(f'{_MODULE}.iofdopt')
            written = []

            def _store_written(rpars, param, **_):
                modified = ModifiedParameterValue(param, rpars, original=None)
                written.append(modified.fmt_value)

            mocker.patch(f'{_MODULE}.parameters.modify',
                         side_effect=_store_written)
            slab = SimpleNamespace(check_atom_collisions=lambda *_: None)
            with execute_in_dir(tmp_path):
                fd_optimization(slab, rpars)
            return written

        return _run

    _written = {
        'theta, phi 360': (('theta', 'THETA 10, PHI 360', 11.0),
                           'THETA 11.0000, PHI 0.0000'),
        'theta beyond grazing': (('theta', 'THETA 89, PHI 30', 90.3),
                                 'THETA 90.0000, PHI 30.0000'),
        }

    @parametrize('args,expect', _written.values(), ids=_written)
    def test_written(self, args, expect, optimize):
        """Check the BEAM_INCIDENCE written to PARAMETERS."""
        assert optimize(*args) == [expect]
