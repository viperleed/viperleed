"""Tests for module _write of viperleed.tleedmlib.files.parameters.

Created on 2023-10-18

@author: Michele Riva (@michele-riva)
"""

import numpy as np
import pytest
from pytest_cases import parametrize

from viperleed.tleedmlib.base import strip_comments
from viperleed.tleedmlib.classes.rparams import Rparams
from viperleed.tleedmlib.files.parameters._write import ModifiedParameterValue


class TestModifiedParameterValue:
    """collection of tests for the (internal) ModifiedParameterValue class."""

    def test_comment_out_only(self):
        """Check attributes when parameter is to be commented out."""
        mod_param = ModifiedParameterValue('param1', 'value1',
                                           only_comment_out=True)
        assert mod_param.only_comment_out
        assert mod_param.param == 'param1'
        assert not strip_comments(mod_param.line)

    def test_param_comment(self):
        """Check attributes when a comment is requested for a parameter."""
        mod_param = ModifiedParameterValue('N_BULK_LAYERS', 'value1',
                                           comment='Test Comment')
        assert not mod_param.only_comment_out
        assert mod_param.param == 'N_BULK_LAYERS'
        assert mod_param.comment == 'Test Comment'
        assign, _comment = (s.strip() for s in mod_param.line.split('!'))
        assert _comment == 'Test Comment'
        assert assign == 'N_BULK_LAYERS = value1'

    _fmt_values = {
        'string': ('N_BULK_LAYERS', 3, '3'),
        'simple': ('V0_IMAG', 2.345, '2.3450'),
        'vector': ('THEO_ENERGIES', (1.2, 3.5, 0.28), '1.2 3.5 0.28'),
        'woods': ('SUPERLATTICE', np.array(((1, 0), (0, 2))),
                  'SUPERLATTICE = (1x2)'),
        'matrix': ('SYMMETRY_CELL_TRANSFORM', np.array(((1, -4), (8, 2))),
                   'SYMMETRY_CELL_TRANSFORM M = 1 -4, 8 2'),
        'BULK_REPEAT float': ('BULK_REPEAT', 9.234, '9.23400'),
        'BULK_REPEAT vector': ('BULK_REPEAT', (1.3, -2.38, 4.997),
                               '[1.30000 -2.38000 4.99700]'),
        'LMAX single': ('LMAX', (10, 10), '10'),
        'LMAX range': ('LMAX', (6, 15), '6-15'),
        'LAYER_CUTS list': ('LAYER_CUTS', (0.1, 0.3, 0.87),
                            '0.1000 0.3000 0.8700'),
        'LAYER_CUTS dz': ('LAYER_CUTS', '0.1 < dz(1.3)', '0.1 < dz(1.3)'),
        }

    @parametrize('attr,value,expected', _fmt_values.values(), ids=_fmt_values)
    def test_format_value(self, attr, value, expected):
        """Check the expected formatting of example values."""
        if attr == 'LAYER_CUTS' and 'dz' in value:
            pytest.xfail(reason='Known bug for formatting LAYER_CUTS with dz')
        rpars = Rparams()
        setattr(rpars, attr, value)
        mod_param = ModifiedParameterValue(attr, rpars)
        assert mod_param.fmt_value == expected

    _fmt_special_values = {
        'BEAM_INCIDENCE': ({'PHI': 0, 'THETA': 12},
                           'THETA 12.0000, PHI 0.0000')
        }

    @parametrize(args=_fmt_special_values.items(), ids=_fmt_special_values)
    def test_format_special_getter(self, args):
        """Check formatting of values that require special getters."""
        param, (rpars_settings, expected) = args
        rpars = Rparams()
        for attr, value in rpars_settings.items():
            setattr(rpars, attr, value)
        mod_param = ModifiedParameterValue(param, rpars)
        assert mod_param.fmt_value == expected

    _invalid = {
        'not a param': 'WHO_KNOWS',
        'special_characters': 's@mething_#wrONg?',
        'no getter': 'V0_IMAG',  # valid, but no special getter method
        }

    @parametrize(attr=_invalid.values(), ids=_invalid)
    def test_raises(self, attr):
        """Check that unknown inputs raise ValueError."""
        rpars = Rparams()
        try:
            delattr(rpars, attr)
        except AttributeError:
            pass
        with pytest.raises(ValueError):
            ModifiedParameterValue(attr, rpars)

    def test_raises_no_formatter(self):
        """Check complaints when asking for a parameter without a formatter."""
        with pytest.raises(ValueError):
            ModifiedParameterValue('PARAM_WITHOUT_FORMATTER', 1)
