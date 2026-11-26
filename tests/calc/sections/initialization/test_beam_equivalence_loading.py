"""Tests for loading BEAM_EQUIVALENCE during initialization."""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-11-25'
__license__ = 'GPLv3+'

import pytest
from pytest_cases import parametrize

from viperleed.calc.sections.initialization import init_domains
from viperleed.calc.sections.initialization import initialization


class _StopError(Exception):
    """Raised to interrupt initialization after the file loading."""


class TestLoadingInInitialization:
    """Check which Rparams objects load BEAM_EQUIVALENCE."""

    @parametrize('subdomain', [True, False])
    def test_loaded_only_for_main(self, subdomain, mocker):
        """Check that per-domain Rparams do not load the file."""
        rpars = mocker.MagicMock()
        rpars.initTheoEnergies.side_effect = _StopError
        with pytest.raises(_StopError):
            initialization(mocker.MagicMock(), rpars, subdomain=subdomain)
        loader = rpars.try_loading_beam_equivalence_file
        assert loader.called is (not subdomain)

    def test_loaded_in_init_domains(self, mocker):
        """Check that the domain initialization loads the file."""
        rpars = mocker.MagicMock(DOMAINS=[])
        rpars.initTheoEnergies.side_effect = _StopError
        with pytest.raises(_StopError):
            init_domains(rpars)
        rpars.try_loading_beam_equivalence_file.assert_called_once()
