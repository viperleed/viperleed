"""Tests for module utils of viperleed.calc.bookkeeper."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2024-10-14'
__license__ = 'GPLv3+'

from pathlib import Path
from unittest.mock import MagicMock
from unittest.mock import patch

import pytest
from pytest_cases import fixture

from viperleed.calc.bookkeeper.utils import discard_files
from viperleed.calc.bookkeeper.utils import make_property
from viperleed.calc.bookkeeper.utils import needs_update_for_attr
from viperleed.calc.bookkeeper.utils import _get_attr_or_dict_item


_UPDATED_VALUE = 'updated_value'
_OTHER_UPDATED_VALUE = 'other_update_value'
_DICT_KEY = 'key'
_TEST = 'test_value'


class MockClass:
    """A class for testing utils functions."""

    def __init__(self, attr_value=None, dict_value=None):
        """Initialize instance with an attribute and a dictionary value."""
        self.attr_value = attr_value
        self.dict_value = (
            {_DICT_KEY: dict_value} if dict_value is not None
            else {}
            )

    def update_from_cwd(self):
        """Set attr_value to _UPDATED_VALUE."""
        self.attr_value = _UPDATED_VALUE

    def other_update(self):
        """Set attr_value to _OTHER_UPDATED_VALUE."""
        self.attr_value = _OTHER_UPDATED_VALUE


class TestDiscardFiles:
    """Tests for the discard_files function."""

    @fixture(name='make_fake_path')
    def factory_make_fake_path(self):
        """Return a fake, existing path."""
        def _make():
            mock_path = MagicMock(spec=Path)
            mock_path.exists.return_value = True
            return mock_path
        return _make

    @fixture(name='mock_dir')
    def fixture_mock_dir(self, make_fake_path):
        """Return a fake path to a directory."""
        mock_dir = make_fake_path()
        mock_dir.is_file.return_value = False
        mock_dir.is_dir.return_value = True
        return mock_dir

    @fixture(name='mock_file')
    def fixture_mock_file(self, make_fake_path):
        """Return a fake path to a file."""
        mock_file = make_fake_path()
        mock_file.is_file.return_value = True
        return mock_file

    @fixture(name='mock_log_error')
    def fixture_mock_log_error(self):
        """Return a fake LOGGER.error."""
        patch_ = patch('viperleed.calc.bookkeeper.utils.LOGGER.error')
        with patch_ as mock_log_error:
            yield mock_log_error

    @fixture(name='mock_rmtree')
    def fixture_mock_rmtree(self):
        """Return a fake shutil.rmtree."""
        patch_ = patch('viperleed.calc.bookkeeper.utils.shutil.rmtree')
        with patch_ as mock_rmtree:
            yield mock_rmtree

    def test_deletion_not_exists(self, mock_file, mock_dir, mock_rmtree):
        """Test deletion of non-existing files/directories."""
        for mock_ in (mock_file, mock_dir):
            mock_.exists.return_value = False
        discard_files(mock_file, mock_dir)
        mock_file.unlink.assert_not_called()
        mock_rmtree.assert_not_called()

    def test_directory_deletion(self, mock_dir, mock_rmtree):
        """Check successful deletion of a directory."""
        discard_files(mock_dir)
        mock_rmtree.assert_called_once_with(mock_dir)

    def test_file_deletion(self, mock_file):
        """Check successful deletion of a single file."""
        discard_files(mock_file)
        mock_file.unlink.assert_called_once()

    def test_file_deletion_fails(self, mock_file, mock_log_error):
        """Test a failing file deletion."""
        # Test failed file deletion with logging
        mock_file.unlink.side_effect = OSError
        discard_files(mock_file)
        mock_log_error.assert_called_once()

    def test_directory_deletion_fails(self, mock_dir,
                                      mock_log_error,
                                      mock_rmtree):
        """Test failure to discard a directory."""
        mock_rmtree.side_effect = OSError
        discard_files(mock_dir)
        mock_log_error.assert_called_once()

    def test_multiple_deletion(self, mock_file, mock_dir, mock_rmtree):
        """Test deletion of multiple files/directories."""
        discard_files(mock_file, mock_dir)
        mock_file.unlink.assert_called_once()
        mock_rmtree.assert_called_once_with(mock_dir)


class TestMakeProperty:
    """Tests for the make_property descriptor creator."""

    def test_without_update(self):
        """Test make_property when needs_update is False."""
        mock_obj = MockClass(attr_value=_TEST)
        prop = make_property('attr_value')
        assert prop.fget(mock_obj) is _TEST

    def test_with_update(self):
        """Test make_property when needs_update is True."""
        mock_obj = MockClass()
        prop = make_property('attr_value',
                             needs_update=True,
                             updater='update_from_cwd')
        with pytest.raises(AttributeError, match=r'update_from_cwd\(\)'):
            prop.fget(mock_obj)
        mock_obj.update_from_cwd()
        assert prop.fget(mock_obj) is _UPDATED_VALUE

    def test_dict_access(self):
        """Test make_property for dictionary attribute access."""
        mock_obj = MockClass(dict_value=_TEST)
        prop = make_property('dict_value[key]')
        assert prop.fget(mock_obj) is _TEST


class TestNeedsUpdateForAttr:
    """Tests for the needs_update_for_attr decorator."""

    def test_decorator(self):
        """Test needs_update_for_attr decorator."""
        mock_obj = MockClass()

        @needs_update_for_attr('attr_value', updater='other_update')
        def some_method(obj):
            return obj.attr_value

        # If attr_value is None, the decorator raises AttributeError
        with pytest.raises(AttributeError, match=r'other_update\(\)'):
            some_method(mock_obj)

        # After update, the method should return the correct value
        mock_obj.other_update()
        assert some_method(mock_obj) is _OTHER_UPDATED_VALUE

    def test_dict_access(self):
        """Test needs_update_for_attr decorator with dictionary access."""
        mock_obj = MockClass()

        @needs_update_for_attr(f'dict_value[{_DICT_KEY}]',
                               updater='update_from_cwd')
        def some_method(obj):
            return obj.dict_value[_DICT_KEY]

        # If dict_value['key'] is None, the decorator raises
        with pytest.raises(AttributeError, match=r'update_from_cwd\(\)'):
            some_method(mock_obj)

        # With updated dictionary, the method returns the correct value
        mock_obj.dict_value[_DICT_KEY] = _TEST
        assert some_method(mock_obj) is _TEST


class TestGetAttrOrDictItem:
    """Tests for the _get_attr_or_dict_item function."""

    def test_attrgetter(self):
        """Test _get_attr_or_dict_item for normal attribute access."""
        _attr_name = 'attr_value'
        kwargs = {_attr_name: _TEST}
        mock_obj = MockClass(**kwargs)
        getter, attr_name = _get_attr_or_dict_item(_attr_name)

        assert getter(mock_obj) is _TEST
        assert attr_name == _attr_name

    def test_dict_access(self):
        """Test _get_attr_or_dict_item for dictionary item access."""
        mock_obj = MockClass(dict_value=_TEST)
        getter, attr_name = _get_attr_or_dict_item(f'dict_value[{_DICT_KEY}]')

        assert getter(mock_obj) is _TEST
        assert attr_name == _DICT_KEY
