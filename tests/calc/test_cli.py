"""Tests for the viperleed.calc command-line interface."""

__authors__ = (
    'Alexander M. Imre (@amimre)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2023-08-02'
__license__ = 'GPLv3+'

import os

import pytest
from pytest_cases import fixture
from pytest_cases import parametrize

from viperleed.calc.bookkeeper.mode import BookkeeperMode
from viperleed.calc.cli import LOG_VERBOSE
from viperleed.calc.cli import LOG_VERY_VERBOSE
from viperleed.calc.cli import ViPErLEEDCalcCLI
from viperleed.calc.cli import _copy_files_from_manifest
from viperleed.calc.cli import _copy_input_files_to_work
from viperleed.calc.cli import _copy_tensors_and_deltas_to_work
from viperleed.calc.cli import _verbosity_to_log_level
from viperleed.calc.constants import DEFAULT_WORK
from viperleed.calc.files.manifest import ManifestFileError
from viperleed.calc.lib.context import execute_in_dir
from viperleed.calc.sections.calc_section import ALL_INPUT_FILES

from ..helpers import filesystem_from_dict
from ..helpers import filesystem_to_dict

_MODULE = 'viperleed.calc.cli'


class TestCalcCliCall:
    """Tests for the __call__ method of ViPErLEEDCalcCLI."""

    @fixture(name='mock_implementation')
    def fixture_mock_implementation(self, mocker):
        """Replace implementation details with mocks."""
        def _mock(exit_code):
            return {
                'bookie': mocker.patch(f'{_MODULE}.Bookkeeper'),
                'inputs': mocker.patch(f'{_MODULE}._copy_input_files_to_work'),
                'log_level': mocker.patch(
                    f'{_MODULE}._verbosity_to_log_level'
                    ),
                'manifest': mocker.patch(
                    f'{_MODULE}._copy_files_from_manifest',
                    ),
                'multiprocess': mocker.patch('multiprocessing.freeze_support'),
                'rmtree': mocker.patch('shutil.rmtree'),
                'run': mocker.patch(f'{_MODULE}.run_calc',
                                    return_value=(exit_code, None)),
                'tensors_deltas': mocker.patch(
                    f'{_MODULE}._copy_tensors_and_deltas_to_work',
                    ),
                }
        return _mock

    def test_delete_workdir(self, tmp_path, mock_implementation):
        """Check the successful removal of the work directory."""
        cli = ViPErLEEDCalcCLI()
        mocks = mock_implementation(exit_code=0)
        with execute_in_dir(tmp_path):
            error = cli(['--delete-workdir'])
        assert not error
        mocks['rmtree'].assert_called_once_with(tmp_path/DEFAULT_WORK)

    def test_delete_workdir_fails(self,
                                  tmp_path,
                                  mock_implementation,
                                  capsys):
        """Check complaints are printed when removing work fails."""
        cli = ViPErLEEDCalcCLI()
        mocks = mock_implementation(exit_code=0)
        mocks['rmtree'].side_effect = OSError
        with execute_in_dir(tmp_path):
            error = cli(['--delete-workdir'])
        assert not error
        captured = capsys.readouterr().out
        # pylint: disable-next=magic-value-comparison
        assert 'Error deleting' in captured

    def test_exit_code(self, tmp_path, mock_implementation, mocker):
        """Check the result of a failed execution."""
        cli = ViPErLEEDCalcCLI()
        exit_code = mocker.MagicMock()
        mocks = mock_implementation(exit_code)
        with execute_in_dir(tmp_path):
            result = cli([])
        assert result is exit_code
        assert bool(result)    # Should fake an error condition
        mocks['rmtree'].assert_not_called()  # work should stay

    def test_success(self, tmp_path, mock_implementation, mocker):
        """Check the successful call to ViPErLEEDCalcCLI."""
        cli = ViPErLEEDCalcCLI()
        mocks = mock_implementation(exit_code=0)
        with execute_in_dir(tmp_path):
            error = cli([])
        assert not error

        # Check bookkeeper calls
        mocks['bookie'].assert_called_once()  # Made instance
        bookkeeper = mocks['bookie'].return_value
        assert bookkeeper.run.mock_calls == [
            mocker.call(mode=BookkeeperMode.CLEAR),
            mocker.call(mode=BookkeeperMode.ARCHIVE),
            ]
        # As well as calls to other functions. Don't bother specifying
        # the arguments for those functions that require the parsed
        # command-line arguments.
        expect_calls = {
            'inputs': mocker.call(tmp_path/DEFAULT_WORK),
            'log_level': None,
            'manifest': mocker.call(tmp_path),
            'multiprocess': mocker.call(),
            'run': None,
            'tensors_deltas': None,
            }
        for mock_name, call in expect_calls.items():
            mock = mocks[mock_name]
            if call is None:
                mock.assert_called_once()
            else:
                assert mock.mock_calls == [call]
        mocks['rmtree'].assert_not_called()


class TestCalcParser:
    """Tests for parsing CLI arguments of viperleed.calc."""

    @fixture(name='calc_parser')
    def fixture_calc_parser(self):
        """Return a CLI argument parser for viperleed.calc."""
        return ViPErLEEDCalcCLI().parser

    def test_parse_version(self, calc_parser):
        """Check that requesting the version exits afterwards."""
        with pytest.raises(SystemExit):
            calc_parser.parse_args(['--version'])

    def test_parse_work(self, calc_parser, tmp_path):
        """Check interpretation of -w flag."""
        parsed = calc_parser.parse_args(['-w', str(tmp_path)])
        assert parsed.work == str(tmp_path)

    @parametrize(v_flag=('-v', '--verbose'))
    def test_parse_verbose(self, calc_parser, v_flag):
        """Check interpretation of -v flag."""
        assert calc_parser.parse_args([v_flag,]).verbose

    @parametrize(v_flag=('-vv', '--very-verbose'))
    def test_parse_very_verbose(self, calc_parser, v_flag):
        """Check interpretation of -vv flag."""
        assert calc_parser.parse_args([v_flag,]).very_verbose


class TestCopyFilesFromManifest:
    """Tests for the _copy_files_from_manifest helper function."""

    @fixture(name='manifest')
    def fixture_manifest(self, tmp_path):
        """Create a manifest file and its contents at `tmp_path`."""
        manifest = 'file1.txt \nfile2  \n  \n\n  folder\n'
        manifest += '''
[folder at one]
one/subfile.txt

[Domain 3 at two]
two/domain_file
'''
        copied = {
            'file1.txt': 'Test file 1',
            'file2': 'Test file 2',
            'folder': {},
            'one': {'subfile.txt': 'Subfile contents'},
            'two': {'domain_file': 'domain file contents'},
            }
        stay = {
            'manifest': manifest,
            'file_not_in_manifest': None,
            'folder_not_in_manifest': {}
            }
        tree = {**copied, **stay}
        filesystem_from_dict(tree, tmp_path)
        return copied, stay

    def test_no_manifest_file(self, tmp_path):
        """Check that no resources are copied if manifest does not exist."""
        dest = tmp_path/'dest'
        dest.mkdir()
        with execute_in_dir(tmp_path):
            _copy_files_from_manifest(dest)
        copied = filesystem_to_dict(dest)
        assert not copied

    def test_copy_successful(self, manifest, tmp_path):
        """Check the successful copy of files/folders."""
        dest = tmp_path/'dest'
        dest.mkdir()
        with execute_in_dir(tmp_path):
            _copy_files_from_manifest(dest)
        expect_copy, expect_stay = manifest
        copied = filesystem_to_dict(dest)
        assert copied == expect_copy
        assert not any(s in copied for s in expect_stay)

    def test_copy_fails(self, manifest, tmp_path, capsys):
        """Check complaints are printed if copying resources fails."""
        manifest_file = tmp_path/'manifest'
        with manifest_file.open('a', encoding='utf-8') as file:
            file.write('two/this_does_not_exist\n')
        self.test_copy_successful(manifest, tmp_path)
        expect_print = f'Error copying two{os.sep}this_does_not_exist'
        stdout = capsys.readouterr().out
        assert expect_print in stdout

    def test_manifest_absolute_paths_raises(self, mocker):
        """Check complaints are printed if copying resources fails."""
        manifest = mocker.patch(f'{_MODULE}.ManifestFile')
        manifest.has_absolute_paths = True
        with pytest.raises(ManifestFileError):
            _copy_files_from_manifest(manifest)


class TestCopyInputFiles:
    """Tests for the _copy_input_files_to_work helper."""

    def test_success(self, mocker):
        """Check the successful copy of input files."""
        copy = mocker.patch('shutil.copy2')
        _copy_input_files_to_work(mocker.MagicMock())
        assert copy.call_count == len(ALL_INPUT_FILES)

    def test_copy_fails(self, mocker, caplog):
        """Check no complaints when files are missing."""
        caplog.set_level(0)  # All messages
        copy = mocker.patch('shutil.copy2', side_effect=FileNotFoundError)
        _copy_input_files_to_work(mocker.MagicMock())
        assert copy.call_count == len(ALL_INPUT_FILES)
        assert not caplog.text


class TestCopyTensorsDeltas:
    """Tests for the _copy_tensors_and_deltas_to_work helper."""

    def test_copy_all(self, mocker):
        """Check correct copying of all Tensors/Deltas."""
        copy = mocker.patch(f'{_MODULE}.copytree_exists_ok')
        _copy_tensors_and_deltas_to_work(mocker.MagicMock(), all_tensors=True)
        n_folders = 2  # Tensors and Deltas folders
        assert copy.call_count == n_folders

    def test_copy_all_not_found(self, mocker, caplog):
        """Check correct copying of all Tensors/Deltas."""
        caplog.set_level(0)  # All messages
        copy = mocker.patch(f'{_MODULE}.copytree_exists_ok',
                            side_effect=FileNotFoundError)
        _copy_tensors_and_deltas_to_work(mocker.MagicMock(), all_tensors=True)
        n_folders_tried = 2  # Tried Tensors and Deltas folders
        assert copy.call_count == n_folders_tried
        assert not caplog.text

    def test_copy_most_recent(self, mocker):
        """Check no complaints when no Tensors exist."""
        mocker.patch(f'{_MODULE}.getMaxTensorIndex', return_value=123)
        mocker.patch('pathlib.Path.is_file', return_value=True)
        copy = mocker.patch('shutil.copy2')
        _copy_tensors_and_deltas_to_work(mocker.MagicMock(), all_tensors=False)
        n_files = 2  # Tensors.zip and Deltas.zip
        assert copy.call_count == n_files

    def test_copy_most_recent_missing(self, mocker):
        """Check no complaints when no Tensors exist."""
        def _is_file(path):
            # pylint: disable-next=magic-value-comparison
            return 'Tensors' in path.name
        mocker.patch(f'{_MODULE}.getMaxTensorIndex', return_value=123)
        mocker.patch('pathlib.Path.is_file', _is_file)
        copy = mocker.patch('shutil.copy2')
        _copy_tensors_and_deltas_to_work(mocker.MagicMock(), all_tensors=False)
        n_files = 1  # Only Tensors.zip, no Deltas.zip
        assert copy.call_count == n_files

    def test_copy_most_recent_none(self, mocker, caplog):
        """Check no complaints when no Tensors exist."""
        caplog.set_level(0)  # All messages
        mocker.patch(f'{_MODULE}.getMaxTensorIndex', return_value=None)
        copy = mocker.patch('shutil.copy2')
        _copy_tensors_and_deltas_to_work(mocker.MagicMock(), all_tensors=False)
        copy.assert_not_called()
        assert not caplog.text


_presets = {  # cli_args: expected_presets
    (): {},
    ('verbose',): {'LOG_LEVEL': LOG_VERBOSE},
    ('very_verbose',): {'LOG_LEVEL': LOG_VERY_VERBOSE},
    ('verbose', 'very_verbose'): {'LOG_LEVEL': LOG_VERY_VERBOSE},
    }

@parametrize('args,expect', _presets.items())
def test_verbosity_to_log(args, expect, mocker):
    """Tests for the _verbosity_to_log_level helper."""
    cli_args = mocker.MagicMock(verbose=False, very_verbose=False)
    for attr_name in args:
        setattr(cli_args, attr_name, True)
    presets = {}
    _verbosity_to_log_level(cli_args, presets)
    assert presets == expect
