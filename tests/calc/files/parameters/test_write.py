"""Tests for module write of viperleed.calc.files.parameters."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2023-10-18'
__license__ = 'GPLv3+'

import numpy as np
import pytest
from pytest_cases import parametrize

from viperleed.calc.classes.rparams import LMax
from viperleed.calc.classes.rparams import LayerCuts
from viperleed.calc.classes.rparams import Rparams
from viperleed.calc.files.parameters.write import ModifiedParameterValue
from viperleed.calc.files.parameters.write import ParametersFileEditor
from viperleed.calc.files.parameters.write import comment_out, modify
from viperleed.calc.lib.base import strip_comments

from ....helpers import execute_in_dir


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
        'LMAX single': ('LMAX', LMax(10, 10), '10'),
        'LMAX range': ('LMAX', LMax(6, 15), '6-15'),
        }

    @parametrize('attr,value,expected', _fmt_values.values(), ids=_fmt_values)
    def test_format_value(self, attr, value, expected):
        """Check the expected formatting of example values."""
        rpars = Rparams()
        setattr(rpars, attr, value)
        mod_param = ModifiedParameterValue(attr, rpars)
        assert mod_param.fmt_value == expected

    _fmt_layer_cuts = {
        'LAYER_CUTS list': ((0.1, 0.3, 0.87), '0.1000 0.3000 0.8700'),
        'LAYER_CUTS dz': ('0.1 < dz(1.3)', '0.1000 < dz(1.3)'),
        }

    @parametrize('val,expected', _fmt_layer_cuts.values(), ids=_fmt_layer_cuts)
    def test_format_layer_cuts_value(self, val, expected):
        """Check the expected formatting of example values."""
        attr = 'LAYER_CUTS'
        rpars = Rparams()
        setattr(rpars, attr, LayerCuts.as_layer_cuts(val))
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


class TestParametersEditor:
    """Collection of tests for the ParametersFileEditor class."""

    def test_comment_out_parameter(self):
        """Test successful execution of comment_out_parameter method."""
        rpars = Rparams()
        editor = ParametersFileEditor(rpars)
        modpar, comment = 'PARAM1', 'Test Comment'
        modified = editor.comment_out_parameter(modpar, comment=comment)
        assert modified.only_comment_out
        assert modified.param == modpar
        assert modified.comment == comment

    def test_modify_parameter_explicit_value(self):
        """Test successful execution of modify_param method."""
        rpars = Rparams()
        editor = ParametersFileEditor(rpars)
        modpar, new_value, comment = 'BULK_LIKE_BELOW', 0.23, 'Test Comment'
        modified = editor.modify_param(modpar, new_value, comment=comment)
        assert not modified.only_comment_out
        assert modified.param == modpar
        assert modified.comment == comment
        assert modified.fmt_value == '0.2300'

    def test_dont_save_nonexisting_file(self):
        """Check no duplicate file is created if one did not exist."""
        rpars = Rparams()
        editor = ParametersFileEditor(rpars, path='_a_non_exis_ting_pa_th__',
                                      save_existing_parameters_file=True)
        editor.save_existing_parameters_file()
        # pylint: disable-next=protected-access
        assert not any(editor._path.glob('*'))

    def test_save_existing_parameters_file(self, read_one_param_file):
        """Check that original file is duplicated correctly."""
        fpath, rpars = read_one_param_file
        rpars.timestamp = timestamp = 'new_timestamp_123456'
        editor = ParametersFileEditor(rpars, path=fpath.parent)
        editor.save_existing_parameters_file()
        new_param = next(fpath.parent.glob(f'PARAMETERS*{timestamp}'), None)
        assert new_param is not None
        with fpath.open('r', encoding='utf-8') as old_file:
            with new_param.open('r', encoding='utf-8') as new_file:
                assert new_file.read() == old_file.read()

    def test_write_modified_nothing(self, read_one_param_file):
        """Check that file is unchanged with no edits."""
        fpath, rpars = read_one_param_file
        editor = ParametersFileEditor(rpars, path=fpath.parent)
        editor.write_modified_parameters()
        ori_path = next(fpath.parent.glob(f'PARAMETERS*{rpars.timestamp}'))
        with fpath.open('r', encoding='utf-8') as mod_file:
            with ori_path.open('r', encoding='utf-8') as ori_file:
                assert mod_file.read() == ori_file.read()

    @staticmethod
    def _modify_one_param(editor, rpars):
        """Edit one parameter, return a comment."""
        rpars.BULK_LIKE_BELOW = 0.9876
        comment = 'Decreasing digits'
        editor.modify_param('BULK_LIKE_BELOW', comment=comment)
        return comment, 'BULK_LIKE_BELOW = 0.9876'

    def test_write_modified_called(self, read_one_param_file):
        """Check editing for explicit call to write_modified_parameters."""
        fpath, rpars = read_one_param_file
        with execute_in_dir(fpath.parent):
            fpath.unlink()
            editor = ParametersFileEditor(rpars)
            comment, value = self._modify_one_param(editor, rpars)
            editor.write_modified_parameters()
            check_file_modified(fpath.name, value, comment)

    def test_write_modified_context(self, read_one_param_file):
        """Check editing behaviour when used as a context manager."""
        fpath, rpars = read_one_param_file
        with ParametersFileEditor(rpars, path=fpath.parent) as editor:
            comment, value = self._modify_one_param(editor, rpars)
            try:
                check_file_modified(fpath, value, comment)
            except AssertionError:
                pass
            else:
                pytest.fail(reason='File was modified too early')
        check_file_modified(fpath, value, comment)

    def test_write_commented_param(self, read_one_param_file):
        """check correct commenting-out of one parameter."""
        fpath, rpars = read_one_param_file
        with ParametersFileEditor(rpars, path=fpath.parent) as editor:
            editor.comment_out_parameter('SITE_DEF')  # Two of them!
        assert all_commented_out(fpath, 'SITE_DEF')


def all_commented_out(fpath, param):
    """Return whether all fpath lines with param are commented out."""
    with fpath.open('r', encoding='utf-8') as _file:
        return not any(strip_comments(line) for line in _file if param in line)


def check_file_modified(fpath, assign_str, comment=''):
    """Check that the expected modifications were carried out."""
    with open(fpath, 'r', encoding='utf-8') as mod_file:
        lines = (line for line in mod_file if assign_str in line)
        line = next(lines, None)
        assert line is not None
        if comment:
            assert line.strip().endswith(comment)


class TestCommentOutAndModifyFunctions:
    """Collection of tests for the public functions of .write."""

    def test_comment_out_param(self, read_one_param_file):
        """Check effective commenting-out of one parameter."""
        fpath, rpars = read_one_param_file
        comment_out(rpars, 'BULK_LIKE_BELOW', path=fpath.parent)
        assert all_commented_out(fpath, 'BULK_LIKE_BELOW')

    def test_modify_param(self, read_one_param_file):
        """Check effective modification of one parameter."""
        fpath, rpars = read_one_param_file
        rpars.LMAX = LMax(3, 16)
        with execute_in_dir(fpath.parent):
            modify(rpars, 'LMAX')
        check_file_modified(fpath, 'LMAX = 3-16')
