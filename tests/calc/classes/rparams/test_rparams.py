"""Tests for module rparams of viperleed.calc.classes.rparams."""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
)
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-11-23'
__license__ = 'GPLv3+'

import logging

from pytest_cases import parametrize

from viperleed.calc.classes.rparams.rparams import Rparams

_PRE = 'gfortran -O2'
_GENERIC_POST = '-llapack -lpthread -lblas -fbacktrace'
_OPENBLAS_POST = '-lopenblas -fbacktrace'
_PROBE = 'viperleed.calc.lib.fortran_utils.gfortran_finds_openblas'


class TestGetFortranCompGfortran:  # pylint: disable=too-few-public-methods
    """Tests for the gfortran flags of Rparams.getFortranComp."""

    _platforms = {  # id: (os_name, probe_result, expect_post, probed)
        'nt, OpenBLAS found': ('nt', True, _OPENBLAS_POST, True),
        'nt, OpenBLAS not found': ('nt', False, _GENERIC_POST, True),
        'posix': ('posix', True, _GENERIC_POST, False),
        }

    @parametrize(platform=_platforms.values(), ids=_platforms)
    def test_post(self, platform, mocker, caplog):
        """Check and log the link flags chosen for gfortran."""
        os_name, found, expect, probed = platform
        rpars = Rparams()  # Before patching os.name: it makes Paths
        mocker.patch('shutil.which', return_value='gfortran')
        mocker.patch('os.name', os_name)
        probe = mocker.patch(_PROBE, return_value=found)
        caplog.set_level(logging.DEBUG)
        rpars.getFortranComp(comp='gfortran')
        pre, post = rpars.FORTRAN_COMP
        assert pre == _PRE
        assert post == expect
        assert probe.called == probed
        blas_flags = expect.replace(' -fbacktrace', '')
        assert f'BLAS flags {blas_flags}' in caplog.text
