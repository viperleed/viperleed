"""Tests for module utils of viperleed.calc.bookkeeper."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2024-10-14'
__license__ = 'GPLv3+'

from pathlib import Path

import pytest
from pytest_cases import fixture
from pytest_cases import parametrize

from viperleed.calc.bookkeeper.errors import NotAnInteractiveShellError
from viperleed.calc.bookkeeper.utils import ask_user_confirmation
from viperleed.calc.bookkeeper.utils import discard_files
from viperleed.calc.bookkeeper.utils import file_contents_identical
from viperleed.calc.bookkeeper.utils import make_property
from viperleed.calc.bookkeeper.utils import needs_update_for_attr
from viperleed.calc.bookkeeper.utils import _get_attr_or_dict_item


_MODULE = 'viperleed.calc.bookkeeper.utils'
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


class MockInput:  # pylint: disable=too-few-public-methods
    """Fake replacement for the input built-in function."""

    def __init__(self, *responses):
        """Initialize with some expected user responses."""
        self._responses = iter(responses)

    def __call__(self, *_):
        """Return a user response."""
        return next(self._responses, 'yes')


class TestAskUserConfirmation:
    """Tests for the ask_user_confirmation function."""

    _user_replies = {
        'no reply': ('', False),  # No by default
        'invalid reply, then no': (
            'please do not',
            'NoPe',  # This is the one that is used
            False,
            ),
        'confirmed': ('YES please', True),
        'invalid reply, then yes': (
            'maybe',
            'y',     # This is the one that is used
            True,
            ),
        }

    @parametrize(replies_and_expect=_user_replies.values(), ids=_user_replies)
    def test_user_confirmed(self, replies_and_expect, mocker):
        """Check the result of asking user confirmation to proceed."""
        *replies, expect = replies_and_expect
        mocker.patch('builtins.input', new=MockInput(*replies))
        mock_stdin = mocker.patch('sys.stdin')
        mock_stdin.isatty.return_value = True
        assert ask_user_confirmation(mode=None) == expect

    def test_not_interactive(self, mocker, caplog):
        """Check complaints when sys.stdin is not an interactive shell."""
        mock_stdin = mocker.patch('sys.stdin')
        mock_stdin.isatty.return_value = False
        mock_mode = mocker.MagicMock(long_flag='--mock-mode')
        with pytest.raises(NotAnInteractiveShellError):
            ask_user_confirmation(mock_mode)
        # pylint: disable-next=magic-value-comparison
        assert 'interactive' in caplog.text


class TestDiscardFiles:
    """Tests for the discard_files function."""

    @fixture(name='make_fake_path')
    def factory_make_fake_path(self, mocker):
        """Return a fake, existing path."""
        def _make():
            mock_path = mocker.MagicMock(spec=Path)
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
    def fixture_mock_log_error(self, mocker):
        """Return a fake LOGGER.error."""
        yield mocker.patch(f'{_MODULE}.LOGGER.error')

    @fixture(name='mock_rmtree')
    def fixture_mock_rmtree(self, mocker):
        """Return a fake shutil.rmtree."""
        yield mocker.patch('shutil.rmtree')

    def test_deletion_not_exists(self, mock_file, mock_dir, mock_rmtree):
        """Test deletion of non-existing files/directories."""
        for mock_ in (mock_file, mock_dir):
            mock_.exists.return_value = False
        n_discarded = discard_files(mock_file, mock_dir)
        mock_file.unlink.assert_not_called()
        mock_rmtree.assert_not_called()
        assert not n_discarded

    def test_directory_deletion(self, mock_dir, mock_rmtree):
        """Check successful deletion of a directory."""
        n_discarded = discard_files(mock_dir)
        mock_rmtree.assert_called_once_with(mock_dir)
        assert n_discarded == 1

    def test_file_deletion(self, mock_file):
        """Check successful deletion of a single file."""
        n_discarded = discard_files(mock_file)
        mock_file.unlink.assert_called_once()
        assert n_discarded == 1

    def test_file_deletion_fails(self, mock_file, mock_log_error):
        """Test a failing file deletion."""
        # Test failed file deletion with logging
        mock_file.unlink.side_effect = OSError
        n_discarded = discard_files(mock_file)
        mock_log_error.assert_called_once()
        assert not n_discarded

    def test_directory_deletion_fails(self, mock_dir,
                                      mock_log_error,
                                      mock_rmtree):
        """Test failure to discard a directory."""
        mock_rmtree.side_effect = OSError
        n_discarded = discard_files(mock_dir)
        mock_log_error.assert_called_once()
        assert not n_discarded

    def test_multiple_deletion(self, mock_file, mock_dir, mock_rmtree):
        """Test deletion of multiple files/directories."""
        n_discarded = discard_files(mock_file, mock_dir)
        mock_file.unlink.assert_called_once()
        mock_rmtree.assert_called_once_with(mock_dir)
        # pylint: disable-next=magic-value-comparison
        assert n_discarded == 2


class TestFileContentsIdentical:
    """Tests for the file_contents_identical function."""

    @fixture(name='mock_files')
    def factory_mock_files(self, tmp_path):
        """Create files with contents in a temporary directory."""
        def _make(*contents, n_files=2):
            n_files = max(n_files, len(contents))
            files = [tmp_path / f'file_{i}' for i in range(n_files)]
            for file, file_contents in zip(files, contents):
                file.write_text(file_contents)
            return files
        return _make

    _valid = {
        'same': (('same', 'same'), True),
        'different': (('contents A', 'contents B'), False),
        'one missing': (('first file',), False),
        'both missing': ((), False),
        }

    @parametrize('contents,expect', _valid.values(), ids=_valid)
    def test_file_compare(self, contents, expect, mock_files):
        """Check the expected return of file_contents_identical."""
        files = mock_files(*contents, n_files=2)
        result = file_contents_identical(*files)
        assert result == expect


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
