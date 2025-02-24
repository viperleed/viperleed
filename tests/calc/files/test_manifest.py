"""Tests for module manifest of viperleed.calc.files."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2025-02-24'
__license__ = 'GPLv3+'

from pathlib import Path

import pytest
from pytest_cases import fixture
from pytest_cases import parametrize

from viperleed.calc.files.manifest import ManifestFile
from viperleed.calc.files.manifest import ManifestFileError
from viperleed.calc.files.manifest import InconsistenPathsError
from viperleed.calc.lib.context import execute_in_dir


_SAMPLE_CONTENTS = 'file1.txt', 'file2.log', 'some_root_folder'


@fixture(name='manifest')
def fixture_manifest(tmp_path):
    """Return a ManifestFile at tmp_path."""
    with execute_in_dir(tmp_path):
        return ManifestFile(*_SAMPLE_CONTENTS)


@fixture(name='mock_open')
def fixure_mock_open(mocker):
    """Fake the results of opening a Path."""
    def _mock(*lines, is_file=True):
        if is_file:
            mocker.patch('pathlib.Path.is_file', return_value=True)
        mock_open = mocker.MagicMock()
        mock_open.__enter__.return_value = lines
        mocker.patch('pathlib.Path.open', return_value=mock_open)
        return mock_open
    return _mock


class TestManifestFile:
    """Tests for the ManifestFile class."""

    @staticmethod
    def check_contents(manifest, path, *items):
        """Ensure manifest has items at path."""
        # pylint: disable-next=protected-access           # OK in tests
        section = manifest._sections[path]
        assert all(item in section for item in items)

    @staticmethod
    def _translate_abs_paths(root, headers_and_lines):
        """Edit absolute paths in headers and contents depending on OS."""
        for header in headers_and_lines.copy():
            if not header:
                continue
            path_str = header.split('at ')[1].replace(']', '')
            path = Path(path_str)
            try:
                path.relative_to(root)
            except ValueError:
                pass  # Handled below
            else:     # Nothing to do, it's relative
                continue
            abs_path_str = str(path.resolve().as_posix())
            lines = headers_and_lines.pop(header)
            header = header.replace(path_str, abs_path_str)
            headers_and_lines[header] = {
                line.replace(path_str, abs_path_str, 1)
                for line in lines
                }

    def check_written_contents(self, manifest, expect):
        """Check that the manifest file at path has expected contents."""
        # Do OS-dependent path changes:
        self._translate_abs_paths(manifest.path, expect)
        written = (manifest.path/'manifest').read_text().splitlines()
        try:
            indices = [written.index(h) for h in expect if h]
        except ValueError as exc:
            raise AssertionError(f'Header not found in {written}') from exc
        ranges = zip((None, *indices), (*indices, None))
        for (start, stop), (header, contents) in zip(ranges, expect.items()):
            written_section = written[start:stop]
            if header:
                written_contents = written_section[1:]
                assert written_section[0] == header
            else:
                written_contents = written_section
            assert len(set(written_contents)) == len(written_contents)
            assert set(written_contents) == contents
            # Either nothing was written, or there
            # should be an empty line at the end
            assert not written_section or not written_section[-1]

    def test_init_no_contents(self):
        """Test ManifestFile initialization."""
        manifest = ManifestFile()
        assert manifest.path == Path.cwd()
        assert len(manifest.paths) == 1
        assert manifest.name == 'manifest'

    def test_init_with_contents(self, manifest):
        """Check correct initialization of a manifest with contents."""
        self.check_contents(manifest, manifest.path, *_SAMPLE_CONTENTS)

    def test_init_with_path(self):
        """Test initialization with a specific path."""
        path_name = '/some/test/path'
        manifest = ManifestFile(path=path_name)
        assert manifest.path == Path(path_name).resolve()

    def test_init_with_other_manifest(self, manifest):
        """Check correct initialization of a manifest with contents."""
        other_contents = 'a.txt', 'other_folder'
        other_manifest = ManifestFile(*other_contents)
        manifest = ManifestFile(other_manifest)
        self.check_contents(manifest, other_manifest.path, *other_contents)
        # The two manifests have the same path
        self.check_contents(manifest, manifest.path, *other_contents)

    def test_add_root_item(self, manifest):
        """Test adding an item to the root of a manifest."""
        manifest.add('newfile.cfg')
        self.check_contents(manifest, manifest.path, 'newfile.cfg')

    def test_add_manifest_other_path(self, manifest):
        """Test adding another ManifestFile."""
        label = 'extra'
        other_manifest = ManifestFile('b.log', path='other')
        manifest.add_manifest(other_manifest, label=label)
        self.check_contents(manifest, other_manifest.path, 'b.log')
        # pylint: disable-next=protected-access           # OK in tests
        assert manifest._labels[other_manifest.path] == label

    def test_add_manifest_same_path(self, manifest, tmp_path):
        """Test adding another ManifestFile."""
        label = 'extra'
        other_manifest = ManifestFile('b.log', path=tmp_path)
        manifest.add_manifest(other_manifest, label=label)
        self.check_contents(manifest, other_manifest.path, 'b.log')
        # pylint: disable-next=protected-access           # OK in tests
        assert manifest._labels[other_manifest.path] != label

    def test_contains(self, manifest):
        """Check the correct outcome of the "in" operation."""
        assert _SAMPLE_CONTENTS[0] in manifest
        assert 'not_in_manifest' not in manifest

    def test_inplace_add_same_path(self, manifest):
        """Test in-place addition."""
        other = ManifestFile('b.log', path=manifest.path)
        manifest += other
        self.check_contents(manifest, manifest.path, *_SAMPLE_CONTENTS)
        self.check_contents(manifest, manifest.path, 'b.log')

    def test_inplace_other_path(self, manifest):
        """Test in-place addition of a sub-manifest."""
        other = ManifestFile('b.log', path=manifest.path/'subpath')
        manifest += other
        self.check_contents(manifest, manifest.path, *_SAMPLE_CONTENTS)
        self.check_contents(manifest, other.path, 'b.log')
        # Ensure the contents of the other have been added identically
        other.add('some_new_file')
        self.check_contents(manifest, other.path, 'some_new_file')

    def test_iter_sections(self, manifest):
        """Test iter_sections() method."""
        sections = dict(manifest.iter_sections())
        assert manifest.path in sections
        for item in _SAMPLE_CONTENTS:
            assert item in sections[manifest.path]

    def test_iter_sections_relative(self, manifest):
        """Test iter_sections() method."""
        sections = dict(manifest.iter_sections(relative=True))
        assert Path() in sections
        for item in _SAMPLE_CONTENTS:
            assert item in sections[Path()]

    def test_iter_sections_relative_raises(self, manifest):
        """Test iter_sections() method."""
        manifest.add_manifest(ManifestFile(path='/some/other/path'))
        with pytest.raises(ManifestFileError):
            dict(manifest.iter_sections(relative=True))

    def test_read_contents_line_raises(self, manifest):
        """Check complaints for an invalid content path."""
        wrong = 'wrong_path/contents'
        with pytest.raises(InconsistenPathsError):
            manifest._read_contents_line(wrong, None, 'right_path')

    _read_line = {  # line, relative_path, expect
        'root content': ('content', '', 'content'),
        'relative content': ('path/to/content', 'path', 'to/content'),
        'abs content': ('/abs/content', '/abs', 'content'),
        }

    @parametrize('line,rel_path,expect', _read_line.values(), ids=_read_line)
    def test_read_contents_line_success(self, line, rel_path, expect):
        """Check complaints for an invalid content path."""
        manifest = ManifestFile()
        # pylint: disable=protected-access                # OK in tests
        manifest._sections['SECTION'] = set()
        manifest._read_contents_line(line, 'SECTION', rel_path)
        assert manifest._sections['SECTION'] == {expect}

    def test_read_header_line_trims_label(self):
        """Ensure _read_header_line() trims spaces in labels."""
        manifest = ManifestFile()
        label, path = 'My Label', '/test/path'
        line = f'[  {label}  at  {path}  ]'
        # pylint: disable-next=protected-access           # OK in tests
        read_path, _ = manifest._read_header_line(line)
        # pylint: disable-next=protected-access           # OK in tests
        assert manifest._labels[read_path] == label
        assert read_path == Path(path).resolve()

    def test_read_missing_file(self, tmp_path):
        """Test read() when the manifest file does not exist."""
        with execute_in_dir(tmp_path):
            manifest = ManifestFile()
            manifest.read()
        assert not any(c for _, c in manifest.iter_sections())

    def test_read_valid_file(self, mock_open):
        """Test reading from a valid manifest file."""
        mock_data = (
            'cwd_file.txt\n',
            '[some root folder at /test]\n',
            '/test/file1.txt\n',
            '/test/file2.log\n',
            '\n',
            '[domain ABC at ./folder]\n',
            'folder/domain_file\n',
            '[at other]\n',
            )
        mock_open(*mock_data)
        manifest = ManifestFile()
        manifest.read()
        expect_contents = {
            Path.cwd(): {'cwd_file.txt'},
            Path('/test').resolve(): {'file1.txt', 'file2.log'},
            Path('./folder').resolve(): {'domain_file'},
            Path('other').resolve(): set(),
            }
        expect_labels = {
            Path.cwd(): 'root',
            Path('/test').resolve(): 'some root folder',
            Path('./folder').resolve(): 'domain ABC',
            Path('other').resolve(): '',
            }
        assert manifest.paths == tuple(expect_contents)
        for path, contents in expect_contents.items():
            self.check_contents(manifest, path, *contents)
        # pylint: disable-next=protected-access           # OK in tests
        assert manifest._labels == expect_labels

    def test_write_children(self, manifest):
        """Check expected results of writing a manifest without children."""
        with execute_in_dir(manifest.path):
            Path('path_to_empty').mkdir()
            Path('path_to_other').mkdir()
            empty_manifest = ManifestFile(path='path_to_empty')
            other_manifest = ManifestFile('content', path='path_to_other')
            manifest.add_manifest(empty_manifest, label='empty')
            manifest.add_manifest(other_manifest, label='other')
            manifest.write()
        assert not manifest.has_absolute_paths
        expect = {
            None: {*_SAMPLE_CONTENTS, ''},
            '[empty at path_to_empty]': {''},
            '[other at path_to_other]': {'path_to_other/content', ''}
            }
        self.check_written_contents(manifest, expect)

    def test_write_children_root_empty(self, tmp_path):
        """Check expected results of writing a manifest without children."""
        with execute_in_dir(tmp_path):
            manifest = ManifestFile()
            empty_manifest = ManifestFile(path='path_to_empty')
            other_manifest = ManifestFile('content', path='path_to_other')
            manifest.add_manifest(empty_manifest, label='empty')
            manifest.add_manifest(other_manifest, label='other')
            manifest.write()
        expect = {
            None: set(),  # No lines for empty root manifest
            '[empty at path_to_empty]': {''},
            '[other at path_to_other]': {'path_to_other/content', ''}
            }
        self.check_written_contents(manifest, expect)

    def test_write_no_children(self, manifest):
        """Check expected results of writing a manifest without children."""
        with execute_in_dir(manifest.path):
            manifest.write()
        assert not manifest.has_absolute_paths
        written = (manifest.path/'manifest').read_text().splitlines()
        assert len(set(written)) == len(written)  # No lines missed
        expect = {*_SAMPLE_CONTENTS, ''}
        assert set(written) == expect

    def test_write_absolute_paths(self, manifest, tmp_path):
        """Test writing when subfolder paths are absolute."""
        with execute_in_dir(tmp_path):
            manifest += ManifestFile('external.txt', path='/external')
            manifest.write()
        assert manifest.has_absolute_paths
        expect = {
            None: {*_SAMPLE_CONTENTS, ''},
            '[at /external]': {'/external/external.txt', ''},
            }
        self.check_written_contents(manifest, expect)

    def test_write_then_read(self, manifest):
        """Check that reading back a written manifest changes nothing."""
        with execute_in_dir(manifest.path):
            Path('path_to_empty').mkdir()
            Path('path_to_other').mkdir()
            empty_manifest = ManifestFile(path='path_to_empty')
            other_manifest = ManifestFile('content', path='path_to_other')
            manifest.add_manifest(empty_manifest, label='empty')
            manifest.add_manifest(other_manifest, label='other')
            manifest.write()
            read_back = ManifestFile()
            read_back.read()
        written = dict(manifest.iter_sections())
        read = dict(read_back.iter_sections())
        assert written == read


class TestManifestFileRaises:
    """Tests for situations that cause exceptions in a ManifestFile."""

    _cant_add = 'abc', None, tuple()

    @parametrize(other=_cant_add)
    def test_add_unsupported(self, manifest, other):
        """Check complaints when trying to add unsupported contents."""
        with pytest.raises(TypeError):
            manifest += other

    def test_read_invalid_header(self, manifest):
        """Test reading an invalid header line."""
        manifest = ManifestFile()
        with pytest.raises(ValueError):
            # pylint: disable-next=protected-access       # OK in tests
            manifest._read_header_line('invalid header\n')
