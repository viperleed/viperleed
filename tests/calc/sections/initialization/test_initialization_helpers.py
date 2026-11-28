"""Tests for helper functions of the initialization section."""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-11-24'
__license__ = 'GPLv3+'

import logging

import pytest
from pytest_cases import parametrize

from viperleed.calc.classes.rparams.special.layer_cuts import LayerCuts
from viperleed.calc.lib.context import execute_in_dir
from viperleed.calc.sections.initialization import _log_auto_layer_cuts
from viperleed.calc.sections.initialization import initialization


_MODULE = 'viperleed.calc.sections.initialization'


class TestLogAutoLayerCuts:
    """Tests for the _log_auto_layer_cuts helper."""

    _cuts = {
        'dz': 'dz(1.2)',
        'dc': 'dc(1.2)',
        'mixed': '0.35 < dz(1.2)',
        }

    @parametrize(layer_cuts=_cuts.values(), ids=_cuts)
    def test_auto_cuts_logged(self, layer_cuts, ag100, caplog):
        """Check that automatically detected cuts are logged."""
        slab, rpars, *_ = ag100
        rpars.LAYER_CUTS = LayerCuts.from_string(layer_cuts)
        slab.create_layers(rpars)
        caplog.set_level(logging.INFO)
        _log_auto_layer_cuts(rpars, slab, user_auto_cuts=True)
        expect = '0.3500 0.4500 0.5500 0.6500 0.7500'
        assert expect in caplog.text

    def test_explicit_cuts_logged_as_given(self, ag100, caplog):
        """Check that explicit cuts are logged as given."""
        slab, rpars, *_ = ag100
        # 0.32 is not a midpoint between Ag(100) layers
        rpars.LAYER_CUTS = LayerCuts.from_string('0.32 < dz(1.2)')
        slab.create_layers(rpars)
        caplog.set_level(logging.INFO)
        _log_auto_layer_cuts(rpars, slab, user_auto_cuts=True)
        expect = '0.3200 0.4500 0.5500 0.6500 0.7500'
        assert expect in caplog.text

    def test_explicit_cuts_not_logged(self, ag100, caplog):
        """Check that nothing is logged for purely explicit cuts."""
        slab, rpars, *_ = ag100
        rpars.LAYER_CUTS = LayerCuts.from_string('0.35 0.45 0.55 0.65 0.75')
        slab.create_layers(rpars)
        caplog.set_level(logging.INFO)
        _log_auto_layer_cuts(rpars, slab, user_auto_cuts=False)
        assert not caplog.text


class _StopInitializationError(Exception):
    """Raised to stop the initialization after logging layer cuts."""


def test_auto_cuts_decided_before_bulk_like_below(ag100, tmp_path, mocker):
    """Check that initialization decides from the user's LAYER_CUTS."""
    slab, rpars, *_ = ag100
    rpars.LAYER_CUTS = LayerCuts.from_string('0.4 < dz(1.3)')
    rpars.BULK_LIKE_BELOW = 0.65
    rpars.BULK_REPEAT = None
    mocker.patch(f'{_MODULE}.parameters')  # Do not edit PARAMETERS
    log_cuts = mocker.patch(f'{_MODULE}._log_auto_layer_cuts',
                            side_effect=_StopInitializationError)
    with execute_in_dir(tmp_path), pytest.raises(_StopInitializationError):
        initialization(slab, rpars)
    assert not any(t.is_auto_cut for t in rpars.LAYER_CUTS)
    log_cuts.assert_called_once_with(rpars, slab, user_auto_cuts=True)
