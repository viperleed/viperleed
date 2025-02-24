"""Tests for cleanup.move_oldruns of viperleed.calc.sections."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2025-02-03'
__license__ = 'GPLv3+'

from copy import deepcopy
from pathlib import Path
import re

import pytest
from pytest_cases import fixture
from pytest_cases import parametrize

from viperleed.calc.constants import DEFAULT_DELTAS
from viperleed.calc.constants import DEFAULT_OUT
from viperleed.calc.constants import DEFAULT_SUPP
from viperleed.calc.constants import DEFAULT_TENSORS
from viperleed.calc.constants import DEFAULT_WORK_HISTORY
from viperleed.calc.constants import LOG_PREFIX
from viperleed.calc.lib.context import execute_in_dir
from viperleed.calc.sections.cleanup import PREVIOUS_LABEL
from viperleed.calc.sections.cleanup import move_oldruns
from viperleed.calc.sections.cleanup import _collect_worhistory_contents
from viperleed.calc.sections.cleanup import _find_next_workistory_contents
from viperleed.calc.sections.cleanup import _find_next_workistory_dir_name
from viperleed.calc.sections.cleanup import _find_next_workistory_run_number
from viperleed.calc.sections.cleanup import _make_new_workhistory_subfolder

from ....helpers import raises_test_exception
from .conftest import _MODULE


all_prerun = parametrize(prerun=(True, False))


class TestCollectWorhistoryContents:
    """Tests for the _collect_worhistory_contents helper."""

    @fixture(name='mock_implementation')
    def fixture_mock_implementation(self, mocker):
        """Replace implementation details with mocks."""
        def _mock(files_and_dirs=None):
            return {
                'contents': mocker.patch(
                    f'{_MODULE}._find_next_workistory_contents',
                    return_value=files_and_dirs or ((), ()),
                    ),
                }
        return _mock

    @fixture(name='run')
    def fixture_run(self, rpars, mock_implementation, tmp_path):
        """Call _collect_worhistory_contents in tmp_path."""
        subfolder = tmp_path/'subfolder'
        def _run(prerun, **kwargs):
            subfolder.mkdir()
            with execute_in_dir(tmp_path):
                mocks = mock_implementation(**kwargs)
                _collect_worhistory_contents(rpars, prerun, subfolder)
                return mocks
        return _run

    _directories = {  # (dirname, prerun) : called
        ('some_dir', True): 'shutil.move',
        ('some_dir', False): 'shutil.copytree',
        (DEFAULT_SUPP, True): 'shutil.move',
        (DEFAULT_SUPP, False): 'shutil.copytree',
        }
    _files = {  # (file, prerun): called
        ('some_file', True): 'shutil.move',
        ('some_file', False): 'shutil.copy2',
        ('control.chem', True): 'shutil.copy2',
        ('control.chem', False): 'shutil.copy2',
        }

    @parametrize('args,raises', _directories.items())
    # pylint: disable-next=too-many-arguments  # 3/6 fixtures
    def test_copy_directory_raises(self, args, raises, run, caplog, mocker):
        """Check that failure to move/copy a directory emits warnings."""
        directory, prerun = args
        mocker.patch(raises, side_effect=OSError)
        run(prerun, files_and_dirs=((), (directory,)))
        expect_log = re.compile(rf'.*Error copying {directory} to .*'
                                r'Files in directory may .*overwritten.\n')
        assert expect_log.fullmatch(caplog.text)

    @parametrize('args,raises', _files.items())
    # pylint: disable-next=too-many-arguments  # 3/6 fixtures
    def test_copy_file_raises(self, args, raises, run, caplog, mocker):
        """Check that failure to move/copy a directory emits warnings."""
        file, prerun = args
        mocker.patch(raises, side_effect=OSError)
        run(prerun, files_and_dirs=((file,), ()))
        expect_log = re.compile(rf'.*Error copying {file} to .*'
                                r'File may .*overwritten.\n')
        assert expect_log.fullmatch(caplog.text)

    @parametrize('args,expect', _directories.items())
    def test_copy_or_move_directory(self, args, expect, run, mocker):
        """Check that a file is moved or copied depending on prerun."""
        directory, prerun = args
        mocks = {'shutil.move': mocker.patch('shutil.move'),
                 'shutil.copytree': mocker.patch('shutil.copytree')}
        called = mocks.pop(expect)
        run(prerun, files_and_dirs=((), (directory,)))
        called.assert_called_once()
        for not_called in mocks.values():
            not_called.assert_not_called()

    @parametrize('args,expect', _files.items())
    def test_copy_or_move_file(self, args, expect, run, mocker):
        """Check that a file is moved or copied depending on prerun."""
        file, prerun = args
        mocks = {'shutil.move': mocker.patch('shutil.move'),
                 'shutil.copy2': mocker.patch('shutil.copy2')}
        called = mocks.pop(expect)
        run(prerun, files_and_dirs=((file,), ()))
        called.assert_called_once()
        for not_called in mocks.values():
            not_called.assert_not_called()


class TestFindNextWorkistoryContents:
    """Tests for the _find_next_workistory_contents helper."""

    @all_prerun
    def test_no_files(self, prerun, rpars, mocker):
        """Check outcome when there are no files/directories."""
        mocker.patch('pathlib.Path.iterdir', return_value=())
        files, directories = _find_next_workistory_contents(rpars, prerun)
        assert not files
        assert not directories

    def test_no_prerun(self, rpars, mocker):
        """Check expected results for a non-prerun call."""
        calc_log = f'{LOG_PREFIX}-timestamp.log'
        collected_files = {'file1.log',
                           'output1.dat',
                           'supp1.txt'}
        collected_dirs = {'some_directory'}
        skipped_files = {  # Not in manifest
            'another_file.txt',
            'io_file.ext',
            calc_log,
            }
        skipped_dirs = {   # Even if in manifest
            DEFAULT_DELTAS,
            DEFAULT_TENSORS,
            DEFAULT_WORK_HISTORY,
            }
        all_files = skipped_files | collected_files
        all_dirs = skipped_dirs | collected_dirs
        for item in (*collected_files, *all_dirs, calc_log):
            rpars.manifest.add(item)
        contents = tuple(Path(f) for f in (*all_files, *all_dirs))
        mocker.patch('pathlib.Path.iterdir', return_value=contents)
        mocker.patch('pathlib.Path.is_file', lambda f: f.name in all_files)
        mocker.patch('pathlib.Path.is_dir', lambda d: d.name in all_dirs)
        files, directories = _find_next_workistory_contents(rpars,
                                                            prerun=False)
        assert set(files) == collected_files
        assert set(directories) == collected_dirs

    def test_prerun(self, rpars, mocker):
        """Check expected results for a prerun call."""
        mocker.patch(f'{_MODULE}._OUT_FILES', ('output1.dat',))
        mocker.patch(f'{_MODULE}._SUPP_FILES', ('supp1.txt',))
        mocker.patch(f'{_MODULE}._IOFILES', ('io_file.ext',))
        rpars.manifest.add('ignore_me.dat')
        rpars.files_to_out.add('this_file_was_generated')
        collected_files = {
            'file1.log',
            'output1.dat',  # OUT
            'supp1.txt',    # SUPP
            }
        skipped_files = {
            'another_file.txt',         # Not a known file
            'io_file.ext'               # IO
            'ignore_me.dat',            # in manifest
            'this_file_was_generated',  # in files_to_out
            }
        collected_dirs = {
            DEFAULT_OUT,
            DEFAULT_SUPP,
            }
        skipped_dirs = {
            'some_directory',
            DEFAULT_TENSORS,
            DEFAULT_DELTAS,
            DEFAULT_WORK_HISTORY,
            }
        all_files = collected_files | skipped_files
        all_dirs = collected_dirs | skipped_dirs
        contents = tuple(Path(f) for f in (*all_files, *all_dirs))
        mocker.patch('pathlib.Path.iterdir', return_value=contents)
        mocker.patch('pathlib.Path.is_file', lambda f: f.name in all_files)
        mocker.patch('pathlib.Path.is_dir', lambda d: d.name in all_dirs)
        files, directories = _find_next_workistory_contents(rpars, prerun=True)
        assert set(files) == collected_files
        assert set(directories) == collected_dirs


class TestFindNextWorkistoryDirName:
    """Tests for the _find_next_workistory_dir_name helper."""

    @fixture(name='mock_rpars')
    def fixture_mock_rpars(self, rpars):
        """Fill an Rparams with some default values."""
        rpars.TENSOR_INDEX = 1
        rpars.timestamp = '20240217'
        return rpars

    @fixture(name='run_number')
    def fixture_run_number(self, mocker):
        """Return a fake run index."""
        mocker.patch(f'{_MODULE}._find_next_workistory_run_number',
                     return_value=1234)

    fake_run = pytest.mark.usefixtures('run_number')

    @fake_run
    def test_logs_no_prerun(self, mock_rpars, mocker):
        """Check the expected result when there are log files."""
        mock_glob = mocker.patch('pathlib.Path.glob', return_value=())
        dirname = _find_next_workistory_dir_name(mock_rpars, prerun=False)
        expected = f't001.r1234_{mock_rpars.timestamp}'
        assert dirname == expected
        mock_glob.assert_not_called()

    @fake_run
    def test_logs_prerun(self, mock_rpars, mocker):
        """Check the expected result when there are log files."""
        logs = (
            Path(f'{LOG_PREFIX}_240101-120000.log'),
            Path(f'{LOG_PREFIX}_240216-235959.log'),  # Most recent
            Path(f'{LOG_PREFIX}_250216-000000.log'),  # In manifest
            )
        mock_rpars.manifest.add(logs[-1].name)
        mocker.patch('pathlib.Path.glob', return_value=logs)
        mocker.patch('pathlib.Path.is_file', return_value=True)
        old_runs_before = mock_rpars.lastOldruns
        old_runs_before_contents = deepcopy(old_runs_before)
        dirname = _find_next_workistory_dir_name(mock_rpars, prerun=True)
        expected = f't001.r1234_{PREVIOUS_LABEL}_240216-235959'
        assert dirname == expected
        assert mock_rpars.lastOldruns is old_runs_before
        assert mock_rpars.lastOldruns == old_runs_before_contents

    _old_runs = {
        'same history': ([1, 2, 3], [1, 2, 3], ''),
        'new runs': ([1, 2, 3, 1, 2], [1, 2, 3], 'RD'),
        'more old runs': ([1, ], [1, 2, 3], ''),
        }

    @fake_run
    @parametrize('history,old,expect', _old_runs.values(), ids=_old_runs)
    def test_old_runs(self, mock_rpars, history, old, expect):
        """Check expected results with both runHistory and lastOldruns."""
        mock_rpars.runHistory = history
        mock_rpars.lastOldruns = old
        dirname = _find_next_workistory_dir_name(mock_rpars, prerun=False)
        if expect:
            expect = '_' + expect
        expected_dirname = f't001.r1234{expect}_{mock_rpars.timestamp}'
        assert dirname == expected_dirname
        assert mock_rpars.lastOldruns is not old
        assert mock_rpars.lastOldruns is not history
        assert mock_rpars.lastOldruns == history

    @fake_run
    def test_no_logs_no_prerun(self, mock_rpars, mocker):
        """Check the expected result when there are no log files."""
        mocker.patch('pathlib.Path.glob', return_value=())
        dirname = _find_next_workistory_dir_name(mock_rpars, prerun=False)
        expected = f't001.r1234_{mock_rpars.timestamp}'
        assert dirname == expected

    @fake_run
    def test_no_logs_prerun(self, mock_rpars, mocker):
        """Check the expected result when there are no log files."""
        mocker.patch('pathlib.Path.glob', return_value=())
        old_runs_before = mock_rpars.lastOldruns
        old_runs_before_contents = deepcopy(old_runs_before)
        dirname = _find_next_workistory_dir_name(mock_rpars, prerun=True)
        expected = f't001.r1234_{PREVIOUS_LABEL}_moved-{mock_rpars.timestamp}'
        assert dirname == expected
        assert mock_rpars.lastOldruns is old_runs_before
        assert mock_rpars.lastOldruns == old_runs_before_contents

    _run_history = {
        'refcalc, delta, search': ([1, 2, 3], 'RDS'),
        'refcalc only': ([1,], 'R'),
        'no abbreviation': ([1, 11, 3, 31, 4, 95], 'RS'),
        }

    @fake_run
    @parametrize('history,expect', _run_history.values(), ids=_run_history)
    def test_run_history_no_oldruns(self, mock_rpars, history, expect):
        """Check the expected result when there is a runHistory."""
        mock_rpars.runHistory = history
        dirname = _find_next_workistory_dir_name(mock_rpars, prerun=False)
        expected_dirname = f't001.r1234_{expect}_{mock_rpars.timestamp}'
        assert dirname == expected_dirname
        assert mock_rpars.lastOldruns is not history
        assert mock_rpars.lastOldruns == history


class TestFindNextWorkistoryRunNumber:
    """Tests for the _find_next_workistory_run_number helper."""

    @all_prerun
    def test_existing_directories(self, rpars, prerun, mocker):
        """Check expected result when worhistory subfolders are present."""
        rpars.TENSOR_INDEX = 234
        existing = (
            Path('t234.r001_'),
            Path('t234.r002_'),
            Path('t234.r099_more_text_with.r578'),
            Path('t999.r999_'),
            )
        mocker.patch('pathlib.Path.iterdir', return_value=existing)
        mocker.patch('pathlib.Path.is_dir', return_value=True)
        number =_find_next_workistory_run_number(rpars, prerun)
        expect = 100  # There's one folder for a different TENSOR_INDEX
        assert number == expect

    @all_prerun
    def test_more_than_1000(self, rpars, prerun, mocker):
        """Check expected result if there are more than 1000 runs."""
        rpars.TENSOR_INDEX = 2
        existing = 't002.r1001_', 't099.r213_'
        existing = tuple(Path(p) for p in existing)
        mocker.patch('pathlib.Path.iterdir', return_value=existing)
        mocker.patch('pathlib.Path.is_dir', return_value=True)
        number =_find_next_workistory_run_number(rpars, prerun)
        expect = 1002
        assert number == expect

    @all_prerun
    def test_not_a_folder(self, rpars, prerun, mocker):
        """Check that only folders are considered."""
        rpars.TENSOR_INDEX = 2
        files = 't002.r012_', 't099.r000_'
        directories = 't002.r001_', 't099.r213_'
        def _mock_isdir(path):
            return path.name in directories
        existing = *files, *directories
        existing = tuple(Path(p) for p in existing)
        mocker.patch('pathlib.Path.iterdir', return_value=existing)
        mocker.patch('pathlib.Path.is_dir', _mock_isdir)
        number =_find_next_workistory_run_number(rpars, prerun)
        expect = 2
        assert number == expect

    @all_prerun
    def test_no_directories(self, rpars, prerun, mocker):
        """Check expected outcome when no workhistory subfolders exist."""
        rpars.TENSOR_INDEX = 123
        mocker.patch('pathlib.Path.iterdir', return_value=())
        number = _find_next_workistory_run_number(rpars, prerun)
        assert number == (0 if prerun else 1)

    @all_prerun
    def test_no_workhistory(self, rpars, prerun, mocker):
        """Check expected outcome when workhistory does not exist."""
        rpars.TENSOR_INDEX = -55
        mocker.patch('pathlib.Path.iterdir', side_effect=FileNotFoundError)
        number = _find_next_workistory_run_number(rpars, prerun)
        assert number == (0 if prerun else 1)

    @all_prerun
    def test_other_subfolders(self, rpars, prerun, mocker):
        """Check expected outcome when stray subfolders are present."""
        rpars.TENSOR_INDEX = 99
        existing = 'some_folder', 't099.rXYZ_', 't099.r008_', 't999.r999_'
        existing = tuple(Path(p) for p in existing)
        mocker.patch('pathlib.Path.iterdir', return_value=existing)
        mocker.patch('pathlib.Path.is_dir', return_value=True)
        number =_find_next_workistory_run_number(rpars, prerun)
        expect = 9
        assert number == expect

    @all_prerun
    def test_unsorted_dirs(self, rpars, prerun, mocker):
        """Check expected outcome when subfolders are not sorted."""
        rpars.TENSOR_INDEX = 2
        existing = (
            Path('test'),
            Path('t002.r005_'),
            Path('t000.r008_'),
            Path('t002.r002_'),
            Path('t002.r004_'),
            )
        mocker.patch('pathlib.Path.iterdir', return_value=existing)
        mocker.patch('pathlib.Path.is_dir', return_value=True)
        number =_find_next_workistory_run_number(rpars, prerun)
        expect = 6
        assert number == expect


class TestMakeNewWorkhistorySubfolder:
    """Tests for the _make_new_workhistory_subfolder helper."""

    @fixture(name='mock_implementation')
    def fixture_mock_implementation(self, mocker):
        """Replace implementation details with mocks."""
        def _mock(dirname='tst_workhistory_dir'):
            return {
                'dirname': mocker.patch(
                    f'{_MODULE}._find_next_workistory_dir_name',
                    return_value=dirname,
                    ),
                }
        return _mock

    @fixture(name='run')
    def fixture_run(self, rpars, mock_implementation, tmp_path):
        """Call _make_new_workhistory_subfolder in tmp_path."""
        def _run(prerun, **kwargs):
            with execute_in_dir(tmp_path):
                mocks = mock_implementation(**kwargs)
                _make_new_workhistory_subfolder(rpars, prerun)
                return mocks
        return _run

    @all_prerun
    def test_create_workhistory_oserror(self, prerun, run, caplog, mocker):
        """Check complaints when creating workhistory fails."""
        _mkdir = Path.mkdir
        def _mkdir_fails(path, **kwargs):
            if path.name == DEFAULT_WORK_HISTORY:
                raise OSError
            _mkdir(path, **kwargs)
        mocker.patch('pathlib.Path.mkdir', _mkdir_fails)
        with pytest.raises(OSError):
            run(prerun)
        expect_log = f'Error creating {DEFAULT_WORK_HISTORY} folder'
        assert expect_log in caplog.text

    @all_prerun
    def test_fails_to_create_subfolder(self, prerun, run, caplog, mocker):
        """Check complaints when creating workhistory fails."""
        _mkdir = Path.mkdir
        dirname = 'test_subfolder'
        def _mkdir_fails(path, **kwargs):
            if path.name == dirname:
                raise OSError
            _mkdir(path, **kwargs)
        mocker.patch('pathlib.Path.mkdir', _mkdir_fails)
        with pytest.raises(OSError):
            run(prerun, dirname=dirname)
        expect_log = 'Error creating'
        assert expect_log in caplog.text

    @all_prerun
    def test_fails_to_create_workhistory(self, prerun, run):
        """Check complaints when creating workhistory fails."""
        with raises_test_exception('pathlib.Path.mkdir'):
            run(prerun)

    @all_prerun
    def test_implementation_called(self, prerun, run, rpars, mocker):
        """Check that helper methods are called as expected."""
        mocks = run(prerun)
        expect_calls = (mocker.call(rpars, prerun),)
        for mock in mocks.values():
            assert mock.call_count == len(expect_calls)
            mock.assert_has_calls(expect_calls)

    @all_prerun
    def test_workhistory_created(self, prerun, run, tmp_path):
        """Check that workhistory exists after execution."""
        workhistory = tmp_path/DEFAULT_WORK_HISTORY
        workhistory_subfolder = workhistory/'test'
        run(prerun, dirname=workhistory_subfolder.name)
        assert workhistory.is_dir()
        assert workhistory_subfolder.is_dir()

    @all_prerun
    def test_workhistory_exists(self, prerun, run, tmp_path):
        """Check that workhistory exists after execution."""
        workhistory = tmp_path/DEFAULT_WORK_HISTORY
        workhistory.mkdir()
        workhistory_subfolder = workhistory/'test'
        run(prerun, dirname=workhistory_subfolder.name)  # No complaint
        assert workhistory_subfolder.is_dir()


class TestMoveOldruns:
    """Tests for the move_oldruns function."""

    @fixture(name='mock_implementation')
    def fixture_mock_implementation(self, tmp_path, mocker):
        """Replace implementation details of move_oldruns with mocks."""
        subfolder = tmp_path/'tst_workhistory_dir'
        def _mock():
            return {
                'make_subfolder': mocker.patch(
                    f'{_MODULE}._make_new_workhistory_subfolder',
                    return_value=subfolder,
                    ),
                'contents': mocker.patch(
                    f'{_MODULE}._collect_worhistory_contents',
                    ),
                'organize_workdir': mocker.patch(
                    f'{_MODULE}.organize_workdir'
                    ),
                }
        return _mock

    @fixture(name='run')
    def fixture_run(self, rpars, mock_implementation, tmp_path):
        """Call move_oldruns at tmp_path, optionally with mocks."""
        def _run(prerun):
            with execute_in_dir(tmp_path):
                mocked = mock_implementation()
                move_oldruns(rpars, prerun)
                return mocked
        return _run

    @all_prerun
    def test_implementation_called(self, prerun, run, rpars, mocker):
        """Check that helper methods are called as expected."""
        mocks = run(prerun)
        expect_calls = {
            'make_subfolder': mocker.call(rpars, prerun),
            'contents': mocker.call(
                rpars,
                prerun,
                mocks['make_subfolder'].return_value,
                ),
            }
        for called_func, call in expect_calls.items():
            mock = mocks[called_func]
            mock.assert_called_once()
            mock.assert_has_calls((call,))

    @all_prerun
    def test_manifest_updated(self, rpars, prerun, run):
        """Check that rpars.manifest is updated as expected."""
        run(prerun)
        should_contain_workhistory = not prerun
        contains_workhistory = DEFAULT_WORK_HISTORY in rpars.manifest
        assert contains_workhistory == should_contain_workhistory

    def test_organize_workdir_calls(self, rpars, run, mocker):
        """Check expected calls to organize_workdir."""
        common_args = {'delete_unzipped': False,
                       'tensors': False,
                       'deltas': False}
        domain_rp = mocker.MagicMock()
        domain_wrk = mocker.MagicMock()
        rpars.domainParams = [
            mocker.MagicMock(rp=domain_rp, workdir=domain_wrk),
            ]
        expect_calls = [
            # One for the main
            mocker.call(rpars, path='', **common_args),
            # and one for each domain
            mocker.call(domain_rp, path=domain_wrk, **common_args),
            ]
        mocks = run(prerun=False)
        organize_workdir = mocks['organize_workdir']
        assert organize_workdir.call_count == len(expect_calls)
        organize_workdir.has_calls(expect_calls)

    def test_prerun_does_not_touch_workdir(self, run):
        """Check that organize_workdir is called only when running."""
        mocks = run(prerun=True)
        organize_workdir = mocks['organize_workdir']
        organize_workdir.assert_not_called()
