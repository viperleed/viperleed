"""Tests for module checker of viperleed.calc.files.parameters."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2023-10-25'
__license__ = 'GPLv3+'

import pytest
from pytest_cases import fixture, parametrize

from viperleed.calc.classes.rparams import Rparams
from viperleed.calc.files.parameters.checker import ParametersChecker
from viperleed.calc.files.parameters.errors import ParameterConflictError

from ....helpers import not_raises

_MODULE = 'viperleed.calc.files.parameters.checker'


@fixture(name='rpars_with_attrs', scope='session')
def factory_rpars_with_attrs():
    """Return an Rparams with specific attributes set."""
    def _make(**attribute_defs):
        rpars = Rparams()
        for name, value in attribute_defs.items():
            setattr(rpars, name, value)
        return rpars
    return _make


@fixture(name='checker')
def fixture_checker():
    """Return a fresh ParametersChecker object."""
    return ParametersChecker()


class TestElementConflicts:
    """Tests for element-name conflict checks."""

    def test_no_collision(self, checker, rpars_with_attrs):
        """Check that no exceptions are raised if there's no collision."""
        rpars = rpars_with_attrs(ELEMENT_RENAME={'A': 'B'},
                                 ELEMENT_MIX={'C': 'D'})
        with not_raises(Exception):
            checker.check_parameter_conflicts(rpars)

    def test_collision(self, checker, rpars_with_attrs):
        """Ensure exceptions are raised if element name collisions exist."""
        rpars = rpars_with_attrs(ELEMENT_RENAME={'A': 'B', 'C': 'D'},
                                 ELEMENT_MIX={'C': 'E'})
        with pytest.raises(ParameterConflictError):
            checker.check_parameter_conflicts(rpars)

    missing = {
        'no mix': {'ELEMENT_RENAME': {'A': 'B'}},
        'no rename': {'ELEMENT_MIX': {'C': 'D'}},
        'neither': {},
        }

    @parametrize(rpars_kwargs=missing.values(), ids=missing)
    def test_missing_parameters(self, rpars_kwargs, checker, rpars_with_attrs):
        """Check no complaints when some values are at their default."""
        rpars = rpars_with_attrs(**rpars_kwargs)
        with not_raises(Exception):
            checker.check_parameter_conflicts(rpars)


class TestFortranCompUpdated:
    """Check that FORTRAN_COMP is updated by the checker."""

    _unchanged = {
        'ifort post': {'FORTRAN_COMP': ('ifort', 'cli flags')},
        'gfortran post': {'FORTRAN_COMP': ('gfortran', 'cli flags')},
        'mpiifort post': {'FORTRAN_COMP_MPI': ('mpiifort', 'cli flags')},
        'mpifort post': {'FORTRAN_COMP_MPI': ('mpifort', 'cli flags')},
        'user': {'FORTRAN_COMP': ('user compiler', '')},
        'user mpi': {'FORTRAN_COMP_MPI': ('user compiler', '')},
        }

    @parametrize(fortran_comp=_unchanged.values())
    def test_not_called(self, fortran_comp, checker, rpars_with_attrs):
        """Check that fortran_comp is not messed with."""
        rpars = rpars_with_attrs(**fortran_comp)
        before = tuple(getattr(rpars, attr) for attr in fortran_comp)
        checker.check_parameter_conflicts(rpars)
        after = tuple(getattr(rpars, attr) for attr in fortran_comp)
        assert after == before

    _verify = {
        'intel': {'FORTRAN_COMP': ('ifort', ''),
                  'FORTRAN_COMP_MPI': ('mpiifort', '')},
        'gnu': {'FORTRAN_COMP': ('gfortran', ''),
                'FORTRAN_COMP_MPI': ('mpifort', '')},
        }

    @parametrize(comp=_verify.values(), ids=_verify)
    def test_verified(self, comp, checker, rpars_with_attrs, mocker):
        """Check that an automatic comp is verified."""
        rpars = rpars_with_attrs(**comp)
        mock_mpi = mocker.patch.object(rpars, 'getFortranMpiComp')
        mock_non_mpi = mocker.patch.object(rpars, 'getFortranComp')
        checker.check_parameter_conflicts(rpars)
        mock_mpi.assert_called_once()
        mock_non_mpi.assert_called_once()


class TestNCoresChecked:
    """Ensure that N_CORES is checked against the number of CPUs."""

    def test_not_defined(self, checker, rpars_with_attrs, caplog):
        """Test that no check is performed when N_CORES is undefined."""
        rpars = rpars_with_attrs()
        checker.check_parameter_conflicts(rpars)
        assert not caplog.text

    @parametrize(cpus=(0, -1))
    # pylint: disable-next=too-many-arguments  # 4/6 are fixtures
    def test_fail_to_detect(self, cpus, checker, rpars_with_attrs,
                            caplog, mocker):
        """Test that no check is performed when failing to detect CPUs."""
        rpars = rpars_with_attrs(N_CORES=5)
        mocker.patch(f'{_MODULE}.available_cpu_count', return_value=cpus)
        checker.check_parameter_conflicts(rpars)
        assert not caplog.text

    def test_enough_cpus(self, checker, rpars_with_attrs, caplog, mocker):
        """Check no warnings with N_CORES < nr. CPUs."""
        rpars = rpars_with_attrs(N_CORES=5)
        mocker.patch(f'{_MODULE}.available_cpu_count', return_value=10)
        checker.check_parameter_conflicts(rpars)
        assert not caplog.text

    def test_too_few_cpus(self, checker, rpars_with_attrs, caplog, mocker):
        """Check warnings are emitted when N_CORES > nr. CPUs."""
        rpars = rpars_with_attrs(N_CORES=5)
        mocker.patch(f'{_MODULE}.available_cpu_count', return_value=1)
        checker.check_parameter_conflicts(rpars)
        msg = 'Consider reducing N_CORES'
        assert msg in caplog.text
