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
from viperleed.calc.files import parameters
from viperleed.calc.files.parameters.utils import Assignment
from viperleed.calc.files.parameters.write import ModifiedParameterValue
from viperleed.calc.files.parameters.write import ParametersFileEditor
from viperleed.calc.files.parameters.write import comment_out
from viperleed.calc.files.parameters.write import modify
from viperleed.calc.lib.context import execute_in_dir
from viperleed.calc.lib.string_utils import strip_comments


class TestModifiedParameterValue:
    """collection of tests for the (internal) ModifiedParameterValue class."""

    def test_comment_out_only(self):
        """Check attributes when parameter is to be commented out."""
        mod_param = ModifiedParameterValue('param1', 'value1',
                                           original=None,
                                           only_comment_out=True)
        assert mod_param.only_comment_out
        assert mod_param.param == 'param1'
        assert not strip_comments(mod_param.line)
        assert not mod_param.original
        assert mod_param.to_assignment() is None

    def test_param_comment(self):
        """Check attributes when a comment is requested for a parameter."""
        mod_param = ModifiedParameterValue('N_BULK_LAYERS', 'value1',
                                           original=None,
                                           comment='Test Comment')
        assert not mod_param.only_comment_out
        assert mod_param.param == 'N_BULK_LAYERS'
        assert mod_param.comment == 'Test Comment'
        assign, _comment = (s.strip() for s in mod_param.line.split('!'))
        assert _comment == 'Test Comment'
        assert assign == 'N_BULK_LAYERS = value1'
        new_assignment = mod_param.to_assignment()
        assert new_assignment == Assignment('value1', 'N_BULK_LAYERS', assign)

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
        mod_param = ModifiedParameterValue(attr, rpars, original=None)
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
        mod_param = ModifiedParameterValue(attr, rpars, original=None)
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
        mod_param = ModifiedParameterValue(param, rpars, original=None)
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
            ModifiedParameterValue(attr, rpars, original=None)

    def test_raises_no_formatter(self):
        """Check complaints when asking for a parameter without a formatter."""
        with pytest.raises(ValueError):
            ModifiedParameterValue('PARAM_WITHOUT_FORMATTER', 1, original=None)


class TestParametersEditor:
    """Collection of tests for the ParametersFileEditor class."""

    def test_comment_out_parameter(self):
        """Test successful execution of comment_out_parameter method."""
        rpars = Rparams()
        editor = ParametersFileEditor(rpars)
        modpar, comment = 'PARAM1', 'Test Comment'
        modified = editor.comment_out_parameter(modpar, comment=comment)
        assert len(modified) == 1
        assert modified[0].only_comment_out
        assert modified[0].param == modpar
        assert modified[0].comment == comment

    def test_modify_parameter_explicit_value(self):
        """Test successful execution of modify_param method."""
        rpars = Rparams()
        editor = ParametersFileEditor(rpars)
        modpar, new_value, comment = 'BULK_LIKE_BELOW', 0.23, 'Test Comment'
        modified = editor.modify_param(modpar, new_value, comment=comment)
        assert len(modified) == 1
        assert not modified[0].only_comment_out
        assert modified[0].param == modpar
        assert modified[0].comment == comment
        assert modified[0].fmt_value == '0.2300'

    def test_modify_twice(self):
        """Test correct subsequent modification of the same parameter."""
        rpars = Rparams()
        editor = ParametersFileEditor(rpars)
        modpar = 'BULK_LIKE_BELOW'
        once = {'new_value': 0.23,
                'comment': 'First edit'}
        twice = {'new_value': 0.99,
                 'comment': 'Second edit'}
        editor.modify_param(modpar, **once)
        modified = editor.modify_param(modpar, **twice)
        assert len(modified) == 1
        assert not modified[0].only_comment_out
        assert modified[0].param == modpar
        assert modified[0].comment == twice['comment']
        assert modified[0].fmt_value == f'{twice["new_value"]:.4f}'

    def test_write_comment_out_then_edit(self, read_one_param_file):
        """Check successive commenting and editing of a parameter."""
        fpath, rpars = read_one_param_file
        with ParametersFileEditor(rpars, path=fpath.parent) as editor:
            editor.comment_out_parameter('LMAX')
            rpars.LMAX = LMax(3, 16)
            editor.modify_param('LMAX')
        old_value = '! LMAX = 8-14'
        comment = '! line automatically changed to:'
        check_file_modified(fpath, old_value, comment)

    def test_write_edit_then_comment_out(self, read_one_param_file):
        """Check successive editing and commenting of a parameter."""
        fpath, rpars = read_one_param_file
        with ParametersFileEditor(rpars, path=fpath.parent) as editor:
            rpars.LMAX = LMax(3, 16)
            editor.modify_param('LMAX')
            editor.comment_out_parameter('LMAX')
        old_value = '! LMAX = 8-14'
        comment = '! line commented out automatically'
        check_file_modified(fpath, old_value, comment)
        for not_present in ('LMAX = 3-16', '! LMAX = 3-16'):
            try:
                check_file_modified(fpath, not_present)
            except AssertionError:
                pass
            else:
                pytest.fail(reason=f'Unexpectedly wrote {not_present!r}')

    def test_write_modified_nothing(self, read_one_param_file):
        """Check that file is unchanged with no edits."""
        fpath, rpars = read_one_param_file
        rpars.files_to_out.clear()  # In case STOP was commented out
        editor = ParametersFileEditor(rpars, path=fpath.parent)
        editor.write_modified_parameters()
        ori_path = next(fpath.parent.glob(f'PARAMETERS*{rpars.timestamp}'))
        with fpath.open('r', encoding='utf-8') as mod_file:
            with ori_path.open('r', encoding='utf-8') as ori_file:
                assert mod_file.read() == ori_file.read()
        assert not rpars.files_to_out

    @staticmethod
    def _modify_one_param(editor, rpars, new_value=0.9876):
        """Edit one parameter, return a comment."""
        rpars.BULK_LIKE_BELOW = new_value
        comment = f'Value changed to {new_value}'
        editor.modify_param('BULK_LIKE_BELOW', comment=comment)
        return comment, f'BULK_LIKE_BELOW = {new_value:.4f}'

    def test_write_modified_called(self, read_one_param_file):
        """Check editing for explicit call to write_modified_parameters."""
        fpath, rpars = read_one_param_file
        with execute_in_dir(fpath.parent):
            fpath.unlink()
            editor = ParametersFileEditor(rpars)
            comment, value = self._modify_one_param(editor, rpars)
            editor.write_modified_parameters()
            check_file_modified(fpath.name, value, comment)
        check_marked_as_edited(rpars)

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
        check_marked_as_edited(rpars)

    def test_write_modified_twice(self, read_one_param_file):
        """Check editing behaviour when used as a context manager."""
        fpath, rpars = read_one_param_file
        with ParametersFileEditor(rpars, path=fpath.parent) as editor:
            comment_1, value_1 = self._modify_one_param(editor, rpars, 0.1234)
            comment_2, value_2 = self._modify_one_param(editor, rpars, 0.7865)
        try:
            check_file_modified(fpath, value_1, comment_1)
        except AssertionError:
            pass
        else:
            pytest.fail(reason='Only the last modification should be written')
        check_file_modified(fpath, value_2, comment_2)
        check_marked_as_edited(rpars)

    def test_write_commented_param(self, read_one_param_file):
        """check correct commenting-out of one parameter."""
        fpath, rpars = read_one_param_file
        with ParametersFileEditor(rpars, path=fpath.parent) as editor:
            editor.comment_out_parameter('BULK_LIKE_BELOW')  # 4 times!
        assert all_commented_out(fpath, 'BULK_LIKE_BELOW')
        check_marked_as_edited(rpars)


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


def check_marked_as_edited(rpars):
    """Check that rpars has PARAMETERS among the edited files."""
    # pylint: disable-next=magic-value-comparison
    assert 'PARAMETERS' in rpars.files_to_out


class TestCommentOutAndModifyFunctions:
    """Collection of tests for the public functions of .write."""

    def test_comment_out_param(self, read_one_param_file):
        """Check effective commenting-out of one parameter."""
        fpath, rpars = read_one_param_file
        comment_out(rpars, 'BULK_LIKE_BELOW', path=fpath.parent)
        assert all_commented_out(fpath, 'BULK_LIKE_BELOW')
        check_marked_as_edited(rpars)

    def test_modify_param(self, read_one_param_file):
        """Check effective modification of one parameter."""
        fpath, rpars = read_one_param_file
        rpars.LMAX = LMax(3, 16)
        with execute_in_dir(fpath.parent):
            modify(rpars, 'LMAX')
        check_file_modified(fpath, 'LMAX = 3-16')
        check_marked_as_edited(rpars)

    def test_add_new_param(self, read_one_param_file):
        """Check correct addition of a new PARAMETER."""
        fpath, rpars = read_one_param_file
        with execute_in_dir(fpath.parent):
            modify(rpars, 'N_BULK_LAYERS')
            rpars_again = parameters.read()
        check_file_modified(fpath, 'N_BULK_LAYERS = 1')
        check_marked_as_edited(rpars)
        assert rpars.readParams == rpars_again.readParams

    def test_modify_param_twice(self, read_one_param_file):
        """Check effective modification of one parameter."""
        fpath, rpars = read_one_param_file
        once, twice = LMax(3, 16), LMax(9, 11)
        with execute_in_dir(fpath.parent):
            for new_value in (once, twice):
                rpars.LMAX = new_value
                modify(rpars, 'LMAX')
        check_file_modified(fpath, 'LMAX = 9-11')
        check_file_modified(fpath,
                            '! LMAX = 3-16',
                            comment='! line automatically changed to:')
        check_marked_as_edited(rpars)

    def test_read_params_updated_when_modified(self, read_one_param_file):
        """Check that rpars.readParams is updated upon modification."""
        fpath, rpars = read_one_param_file
        rpars.LMAX = LMax(3, 16)
        with execute_in_dir(fpath.parent):
            modify(rpars, 'LMAX')
            rpars_again = parameters.read()
        assert rpars.readParams == rpars_again.readParams

    def test_read_params_updated_when_commented(self, read_one_param_file):
        """Check that rpars.readParams is updated upon modification."""
        fpath, rpars = read_one_param_file
        with execute_in_dir(fpath.parent):
            comment_out(rpars, 'BULK_LIKE_BELOW')  # 4 of them!
            rpars_again = parameters.read()
        assert rpars.readParams == rpars_again.readParams
