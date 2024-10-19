"""Tests for module meta of viperleed.calc.bookkeeper.history."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2024-10-13'
__license__ = 'GPLv3+'

from configparser import ConfigParser
import hashlib
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest
from pytest_cases import fixture
from pytest_cases import parametrize
from pytest_cases import parametrize_with_cases

from viperleed.calc.bookkeeper.history.meta import BookkeeperMetaFile
from viperleed.calc.bookkeeper.history.meta import _read_chunked
from viperleed.calc.bookkeeper.history.meta import _HEADER
from viperleed.calc.bookkeeper.history.meta import _METADATA_NAME
from viperleed.calc.bookkeeper.history.meta import _SECTIONS


_MOCK_HASH = 'this is a fake hash'
_MOCK_OTHER_HASH = 'this is another fake hash'


@fixture(name='metadata_file')
def fixture_metadata_file(tmp_path):
    """Return a metadata file at tmp_path."""
    parser, file = ConfigParser(), tmp_path/_METADATA_NAME
    parser.read_dict(_SECTIONS)
    parser['archived']['hash'] = _MOCK_HASH
    parser['archived']['with'] = _MOCK_OTHER_HASH
    with file.open('w', encoding='utf-8') as metadata:
        metadata.write(_HEADER)
        parser.write(metadata)
    return file


class CasesHistoryTree:
    """Collection of test cases that mimic a file-system tree."""

    def _tree_from_dict(self, tree_as_dict, root):
        """Create files and directories from a dictionary version of a tree."""
        for entry, contents in tree_as_dict.items():
            path = root / entry
            if isinstance(contents, dict):
                path.mkdir()
                self._tree_from_dict(contents, path)
                continue
            path.write_text(contents)

    def case_empty_tree(self, tmp_path):
        """Return a tree without any contents."""
        return tmp_path, None

    def case_only_metadata_file(self, tmp_path):
        """Return a tree containing only the metadata file."""
        self._tree_from_dict({_METADATA_NAME: ''}, tmp_path)
        return tmp_path, None

    def case_sample_tree(self, tmp_path):
        """Return the root of a fake history folder and its hash."""
        # Create a bunch of files and folders
        sample_file_contents = 'These are the test contents of a file'
        tree = {
            _METADATA_NAME : sample_file_contents,  # Not used for hash
            'file_1': sample_file_contents,
            'file_2': sample_file_contents,
            'folder_1': {
                'subfile1': sample_file_contents,
                'empty subfolder1': {},
                'subfolder2': {'subfile1.1': sample_file_contents},
                },
            }
        self._tree_from_dict(tree, tmp_path)
        return tmp_path, '9d46c0aae42777a813780a6ba49befec'


@fixture(name='meta')
def fixture_bookkeeper_meta(mock_path):
    """Return a BookkeeperMetaFile instance at `mock_path`."""
    return BookkeeperMetaFile(mock_path)


@fixture(name='meta_read')
def fixture_bookkeeper_meta_read(metadata_file):
    """Return a BookkeeperMetaFile instance read from file."""
    meta = BookkeeperMetaFile(metadata_file.parent)
    meta.read()
    return meta


class TestBookkeeperMetaFile:
    """Tests for the BookkeeperMetaFile class."""

    def test_init(self, meta, mock_path):
        """Test initialization of BookkeeperMetaFile."""
        assert meta.path == Path(mock_path)
        # pylint: disable-next=protected-access           # OK in tests
        parser = meta._parser
        assert parser is not None
        assert isinstance(parser, ConfigParser)
        assert parser.has_section('archived')
        assert meta.parent is None

    def test_file_property(self, meta, mock_path):
        """Test file property."""
        expected_file = Path(mock_path) / _METADATA_NAME
        assert meta.file == expected_file

    def test_hash_property(self, meta):
        """Test accessing hash_ after compute_hash is called."""
        hash_ = hashlib.md5(b'test')
        # pylint: disable-next=protected-access           # OK in tests
        meta._hash = hash_
        assert meta.hash_ == hash_.hexdigest()

    def test_collect_from(self, meta):
        """Test collect_from method."""
        other_meta = MagicMock()
        other_meta.hash_ = _MOCK_HASH
        meta.collect_from(other_meta)
        assert meta.parent == _MOCK_HASH

    module = 'viperleed.calc.bookkeeper.history.meta'

    @patch(f'{module}.hashlib.md5')
    @patch(f'{module}.BookkeeperMetaFile._update_hash_from_folder')
    def test_compute_hash_fake(self, mock_update_folder, mock_md5, meta):
        """Test compute_hash method by mocking methods."""
        mock_md5.return_value = _MOCK_OTHER_HASH
        meta.compute_hash()
        mock_update_folder.assert_called_once_with(meta.path)
        assert meta.hash_ == _MOCK_OTHER_HASH
        # pylint: disable-next=protected-access           # OK in tests
        assert meta._parser['archived']['hash'] == _MOCK_OTHER_HASH

    @parametrize_with_cases('tree', cases=CasesHistoryTree)
    def test_compute_hash(self, tree):
        """Test computation of hash of a real file-system tree."""
        root, expect_hash = tree
        meta = BookkeeperMetaFile(root)
        meta.compute_hash()
        try:
            assert meta.hash_ == expect_hash
        except AttributeError:
            # pylint: disable-next=protected-access       # OK in tests
            hash_ = meta._hash
            assert hash_ is None
            assert hash_ is expect_hash

    def test_read(self, meta_read):
        """Test read method when the file exists."""
        assert meta_read.hash_ == _MOCK_HASH
        assert meta_read.parent == _MOCK_OTHER_HASH

    def test_write(self, meta_read):
        """Check correct writing to file."""
        contents_before = meta_read.file.read_text(encoding='utf-8')
        meta_read.file.unlink()
        meta_read.write()
        contents_after = meta_read.file.read_text(encoding='utf-8')
        assert contents_before == contents_after


class TestBookkeeperMetaFileRaises:
    """Tests for complaints raised by BookkeeperMetaFile."""

    def test_hash_property_before_computation(self, meta):
        """Test accessing hash_ before compute_hash is called."""
        with pytest.raises(AttributeError):
            _ = meta.hash_

    def test_read_file_not_found(self, meta):
        """Test read method when the metadata file is missing."""
        # NB: Path instances have immutable .is_file. Setting it
        # onto the class rather than the instance circumvents this.
        mock_is_file = patch.object(type(meta.path),
                                    'is_file',
                                    return_value=False)
        with mock_is_file, pytest.raises(FileNotFoundError):
            meta.read()

    def test_write_without_hash(self, meta):
        """Test complaints when writing without a hash."""
        with pytest.raises(TypeError):
            meta.write()


_read_chunked_bytes = (
    (b'data1', b'data2'),
    (b'some longer data1', b'some longer data2', b'shorter'),
    )

@parametrize(chunks=_read_chunked_bytes)
def test_read_chunked(chunks):
    """Test _read_chunked function."""
    mock_file = MagicMock()
    mock_file.read.side_effect = (*chunks, b'')
    result = tuple(_read_chunked(mock_file, len(chunks[0])))
    assert result == chunks
    mock_file.read.assert_called()
