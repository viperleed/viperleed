"""Tests for module fs_util of viperleed.calc.lib."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2025-03-11'
__license__ = 'GPLv3+'

from pathlib import Path
import shutil
import sys

import pytest

from viperleed.calc.lib.fs_util import copytree_exists_ok, move

from ...helpers import filesystem_from_dict
from ...helpers import filesystem_to_dict


class TestCopytreeExistsOk:
    """Tests for the copytree_exists_ok function."""

    def test_after_py38(self, mocker):
        """Check that for python>=3.8 the shutil function is used."""
        mocker.patch('sys.version_info', (3, 8))
        sh_copy = mocker.patch('shutil.copytree')
        args = mocker.MagicMock(), mocker.MagicMock()
        copytree_exists_ok(*args)
        sh_copy.assert_called_once_with(*args, dirs_exist_ok=True)

    def test_before_py38(self, mocker, tmp_path):
        """Check usage of the backward-compatible version for python<=3.8."""
        mocker.patch('sys.version_info', (3, 7))
        self.test_existing_directory(tmp_path)

    def test_existing_directory(self, tmp_path):
        """Test copying a directory when destination already exists."""
        before = {'source': {'file.txt': 'content',
                             'subfolder': {'subfile': ''}},
                  'destination': {}}
        expect = before.copy()
        expect['destination'] = before['source']
        filesystem_from_dict(before, tmp_path)
        src, dst = (tmp_path/f for f in expect)
        copytree_exists_ok(src, dst)
        after = filesystem_to_dict(tmp_path)
        assert after == expect

    def test_source_does_not_exist(self, tmp_path):
        """Test complaints if the source is not a directory."""
        src = tmp_path / 'source'
        dst = tmp_path / 'destination'

        with pytest.raises(FileNotFoundError):
            copytree_exists_ok(src, dst)

    def test_new_directory(self, tmp_path):
        """Test copying a directory when destination does not exist."""
        before = {'source': {'file.txt': 'content',
                             'subfolder': {'subfile': ''}},
                 }
        expect = before.copy()
        expect['destination'] = before['source']
        filesystem_from_dict(before, tmp_path)
        src, dst = (tmp_path/f for f in expect)
        copytree_exists_ok(src, dst)
        after = filesystem_to_dict(tmp_path)
        assert after == expect


class TestMove:
    """Tests for the move function."""

    def test_after_py39(self, mocker):
        """Test that on Python >=3.9 Paths are passed along unaltered."""
        mocker.patch('sys.version_info', (3, 9))
        mocker.patch('pathlib.Path.exists', return_value=False)
        sh_move = mocker.patch('shutil.move')
        args = mocker.MagicMock(spec=Path), mocker.MagicMock(spec=Path)
        move(*args)
        sh_move.assert_called_once_with(*args, copy_function=shutil.copy2)

    def test_before_py39(self, mocker, tmp_path):
        """Test move function for Python <3.9."""
        mocker.patch('sys.version_info', (3, 7))
        self.test_directory(tmp_path)

    def test_directory(self, tmp_path):
        """Test successfully moving a directory."""
        before = {'src': {'file': 'contents', 'dir': {'file': 'contents'}}}
        expect = {'dst': before['src']}
        filesystem_from_dict(before, tmp_path)
        args = tmp_path/'src', tmp_path/'dst'
        move(*args)
        after = filesystem_to_dict(tmp_path)
        assert after == expect

    def test_file(self, mocker, tmp_path):
        """Test successfully moving a file."""
        before = {'file.txt': 'content'}
        expect = {'moved.txt': before['file.txt']}
        filesystem_from_dict(before, tmp_path)
        args = tmp_path/'file.txt', tmp_path/'moved.txt'
        move(*args)
        after = filesystem_to_dict(tmp_path)
        assert after == expect

    def test_dst_exists(self, tmp_path):
        """Test moving a file when the destination already exists."""
        tree = {'src': {'file': 'contents', 'dir': {'file': 'contents'}},
                'dst': {'src': {'already_there.txt': ''}}}
        filesystem_from_dict(tree, tmp_path)
        args = (tmp_path/f for f in tree)
        with pytest.raises(FileExistsError):
            move(*args)
