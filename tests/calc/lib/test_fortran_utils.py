"""Tests for module fortran_utils of viperleed.calc.lib."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2024-08-26'
__license__ = 'GPLv3+'
import subprocess

import pytest
from pytest_cases import parametrize

from viperleed.calc.lib.fortran_utils import wrap_fortran_line
from viperleed.calc.lib.fortran_utils import get_mpifort_version
from viperleed.calc.lib.fortran_utils import MpifortNotFoundError
from viperleed.calc.lib.fortran_utils import CouldNotDeterminMpifortVersionError


class TestWrapFortranLine:  # pylint: disable=too-few-public-methods
    """Tests for the wrap_fortran_line function."""

    _valid = {
        'at limit': ('x' * 72, 'x' * 72),
        'empty': ('', ''),
        'non ASCII': (
            'é' * 156,
            'é' * 72 + '&\n     &' + 'é' * 66 + '&\n     &' + 'é' * 18,
            ),
        'short': ('12345', '12345'),
        'spaces': ('123456' + ' ' * 76,
                   '123456' + ' ' * 66 + '&\n     &' + ' ' * 10),
        'wrap once': ('x' * 78, 'x' * 72 + '&\n     &' + 'x' * 6),
        'wrap twice exactly': (
            'x' * 72 + 'y' * 66,
            'x' * 72 + '&\n     &' + 'y' * 66,
            ),
        'wrap multiple': (
            'x' * 72 + 'y' * 66 + 'z' * 12,
            'x' * 72 + '&\n     &' + 'y' * 66 + '&\n     &' + 'z' * 12,
            ),
        'code line': (
            '      DATA NLMS(1), NLMS(2), NLMS(3), NLMS(4), NLMS(5), NLMS(6), '
            'NLMS(7), NLMS(8), NLMS(9), NLMS(10), NLMS(11), NLMS(12), '
            'NLMS(13), NLMS(14), NLMS(15), NLMS(16), NLMS(17), NLMS(18)/ 0, '
            '76, 284, 809, 1925, 4032, 7680, 13593, 22693, 36124, 55276, '
            '81809, 117677, 165152, 226848, 305745, 405213, 529036/',
            '''\
      DATA NLMS(1), NLMS(2), NLMS(3), NLMS(4), NLMS(5), NLMS(6), NLMS(7)&
     &, NLMS(8), NLMS(9), NLMS(10), NLMS(11), NLMS(12), NLMS(13), NLMS(1&
     &4), NLMS(15), NLMS(16), NLMS(17), NLMS(18)/ 0, 76, 284, 809, 1925,&
     & 4032, 7680, 13593, 22693, 36124, 55276, 81809, 117677, 165152, 22&
     &6848, 305745, 405213, 529036/'''
            )
        }

    @parametrize('string,expect', _valid.values(), ids=_valid)
    def test_valid(self, string, expect):
        """Check expected outcome with valid arguments."""
        assert wrap_fortran_line(string) == expect

class TestGetMpifortVersion:
    def test_get_mpifort_version_success(self, monkeypatch):
        """Test successful retrieval of mpifort version."""

        # Mock function to simulate subprocess.run behavior
        def mock_run_success(cmd, shell, check, stdout=None, stderr=None):
            if "mpifort --version | sed" in " ".join(cmd):
                # Simulate the expected output for version extraction
                return subprocess.CompletedProcess(
                    args=cmd, returncode=0, stdout=b"13.2.0\n"
                )
            elif "mpifort --version" in " ".join(cmd):
                # Simulate the initial check call
                return subprocess.CompletedProcess(args=cmd, returncode=0)
            raise ValueError("Unexpected command")  # Catch unexpected calls

        monkeypatch.setattr(subprocess, "run", mock_run_success)

        # Run the function and check the version output
        version = get_mpifort_version()
        assert str(version) == "13.2.0"


    def test_get_mpifort_version_not_installed(self, monkeypatch):
        """Test behavior when mpifort is not installed."""

        def mock_run_not_installed(cmd, shell, check, stdout=None, stderr=None):
            raise subprocess.CalledProcessError(1, cmd)

        monkeypatch.setattr(subprocess, "run", mock_run_not_installed)
        with pytest.raises(MpifortNotFoundError):
            get_mpifort_version()


    def test_get_mpifort_version_could_not_determine(self, monkeypatch):
        """Test behavior when mpifort version cannot be determined."""

        def mock_run_could_not_determine(
            cmd, shell, check, stdout=None, stderr=None
        ):
            if "mpifort --version" in cmd:
                return subprocess.CompletedProcess(
                    args=cmd, returncode=0
                )  # Simulate successful initial check
            raise subprocess.CalledProcessError(1, cmd)

        monkeypatch.setattr(subprocess, "run", mock_run_could_not_determine)
        with pytest.raises(CouldNotDeterminMpifortVersionError):
            get_mpifort_version()
