"""Tests for module root_explorer of viperleed.calc.bookkeeper."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2024-10-16'
__license__ = 'GPLv3+'

from collections import namedtuple
from itertools import combinations
from itertools import product
from operator import attrgetter
from pathlib import Path
import time
from unittest.mock import patch

import pytest
from pytest_cases import fixture
from pytest_cases import parametrize

from viperleed.calc.bookkeeper.constants import STATE_FILES
from viperleed.calc.bookkeeper.constants import EDITED_SUFFIX
from viperleed.calc.bookkeeper.constants import ORI_SUFFIX
from viperleed.calc.bookkeeper.errors import FileOperationFailedError
from viperleed.calc.bookkeeper.root_explorer import DomainRootExplorer
from viperleed.calc.bookkeeper.root_explorer import RootExplorer
from viperleed.calc.bookkeeper.root_explorer import TensorAndDeltaInfo
from viperleed.calc.constants import DEFAULT_DELTAS
from viperleed.calc.constants import DEFAULT_OUT
from viperleed.calc.constants import DEFAULT_SUPP
from viperleed.calc.constants import DEFAULT_TENSORS
from viperleed.calc.constants import LOG_PREFIX
from viperleed.calc.lib.time_utils import DateTimeFormat

from ...helpers import filesystem_from_dict
from ...helpers import filesystem_to_dict
from ...helpers import make_obj_raise
from ...helpers import not_raises
from ...helpers import raises_exception
from .conftest import MOCK_STATE_FILES


_MOCK_ROOT = '/mock/root'
_MODULE = 'viperleed.calc.bookkeeper.root_explorer'
_TEST_STRING = 'this is a test value'
patch_discard = patch(f'{_MODULE}.discard_files')


@fixture(name='explorer')
def fixture_explorer(mocker):
    """Return a RootExplorer with a fake bookkeeper."""
    mock_bookkeeper = mocker.MagicMock()
    root_path = Path(_MOCK_ROOT)
    mocker.patch('pathlib.Path.iterdir', return_value=())
    return RootExplorer(root_path, mock_bookkeeper)


@fixture(name='tensors')
def fixture_tensor_and_delta_info():
    """Return a TensorAndDeltaInfo with a mock root path."""
    return TensorAndDeltaInfo(Path('/fake/root'))


@fixture(name='mock_history')
def fixture_mock_history(mocker):
    """Return a fake HistoryExplorer."""
    def _mock(max_run_per_tensor):
        return mocker.MagicMock(
            max_run_per_tensor=max_run_per_tensor,
            last_folder_and_siblings=[mocker.MagicMock(tensor_num=t)
                                      for t in max_run_per_tensor],
            )
    return _mock


@fixture(name='patched_path')
def fixture_patched_path(explorer, mocker):
    """Return a version of pathlib.Path with fake methods."""
    MockedPath = namedtuple('MockedPath', ('read', 'write', 'glob'))
    notes_file = explorer.path / 'notes.txt'
    mock_glob = mocker.patch('pathlib.Path.glob',
                             return_value=iter([notes_file]))
    mock_read = mocker.patch('pathlib.Path.read_text',
                             return_value=_TEST_STRING)
    mock_write = mocker.patch('pathlib.Path.write_text')
    yield MockedPath(mock_read, mock_write, mock_glob)


def called_or_not(condition):
    """Return the name of an assertion method depending on `condition`."""
    return 'assert_called' if condition else 'assert_not_called'


def mock_exists(collection):
    """Return a callable that can replace pathlib.Path.exists."""
    def _exists(path):
        return path in collection
    return _exists


class TestRootExplorer:
    """Tests for the RootExplorer class."""

    _recent = {
        'no log': None,
        'log file': _TEST_STRING,
        }

    @parametrize(expect=_recent.values(), ids=_recent)
    def test_calc_timestamp(self, explorer, expect, mocker):
        """Test calc_timestamp property."""
        most_recent = (expect if expect is None
                       else mocker.MagicMock(timestamp=expect))
        # pylint: disable-next=protected-access           # OK in tests
        explorer._logs = mocker.MagicMock(most_recent=most_recent)
        assert explorer.calc_timestamp == expect

    _called = {
        # outer_call: {mocked_attr: inner_call or None}
        # if None, take mocked_attr as the one that is called
        'clear_for_next_calc_run' : {
            '_logs': 'logs.discard',
            '_remove_out_and_supp': None,
            '_remove_ori_files': None,
            'collect_info': None,
            },
        'revert_to_previous_calc_run' : {
            '_logs': 'logs.discard',
            '_remove_out_and_supp': None,
            '_replace_state_files_from_ori': None,
            'collect_info': None,
            },
        }

    @parametrize('method_name,called', _called.items(), ids=_called)
    @patch_discard
    def test_calls_methods(self, _, explorer, method_name, called,
                           check_methods_called):
        """Check calling of expected methods."""
        check_methods_called(explorer, method_name, **called)

    _clean_up_root_internals = tuple(product((True, False), repeat=3))

    @parametrize(inner_returns=_clean_up_root_internals)
    def test_clean_up_root_return_value(self, inner_returns, explorer, mocker):
        """Check the expected return value of _clean_up_root."""
        logs_return, supp_out_return, ori_return = inner_returns
        explorer.collect_info()  # So .logs is available
        mocker.patch.object(explorer.logs, 'discard', return_value=logs_return)
        mocker.patch.object(explorer, '_remove_out_and_supp',
                            return_value=supp_out_return)
        # pylint: disable-next=protected-access           # OK in tests
        result = explorer._clean_up_root(lambda: ori_return)
        assert result == (logs_return or supp_out_return or ori_return)

    @parametrize(method=('revert_to_previous_calc_run',
                         'clear_for_next_calc_run'))
    def test_clean_up_root_callers_return(self, method, explorer, mocker):
        """Check the expected return value of callers of _clean_up_root."""
        expect = mocker.MagicMock()
        mocker.patch.object(explorer, '_clean_up_root', return_value=expect)
        result = getattr(explorer, method)()
        assert result is expect

    def test_collect_info(self, explorer, mock_attributes, mocker):
        """Check calling of expected method when collecting info."""
        mock_logs = mocker.patch(f'{_MODULE}.LogFiles')
        mock_attributes(mock_logs, 'collect')
        mock_attributes(explorer, '_collect_files_to_archive')
        mock_attributes(explorer, '_find_potential_domain_subfolders')

        explorer.collect_info()

        explorer.logs.collect.assert_called_once()
        called = {
            '_collect_files_to_archive',
            'collect_info',
            '_find_potential_domain_subfolders',
            }
        for method in dir(explorer):
            try:
                called_ok = getattr(explorer, called_or_not(method in called))
            except AttributeError:
                pass
            else:
                called_ok()

    def test_collect_files_to_archive(self, explorer, mocker):
        """Test the _collect_files_to_archive method."""
        # pylint: disable-next=protected-access           # OK in tests
        explorer._logs = mocker.MagicMock(calc=(explorer.path / 'calc.log',))
        mocker.patch.object(explorer, 'workhistory')
        explorer.workhistory.find_current_directories.return_value = (
            explorer.path / 'workhist_folder',
            )
        expected_files = (
            DEFAULT_OUT,
            DEFAULT_SUPP,
            'calc.log',
            'workhist_folder',
            )
        mocker.patch('pathlib.Path.is_file', return_value=True)
        # pylint: disable-next=protected-access           # OK in tests
        explorer._collect_files_to_archive()
        # pylint: disable-next=protected-access           # OK in tests
        to_archive = explorer._files_to_archive
        assert to_archive == tuple(explorer.path / f for f in expected_files)

    @parametrize(has_domains=(True, False))
    def test_find_domains(self, has_domains, explorer, mocker):
        """Test the _find_potential_domain_subfolders method."""
        mock_finder = mocker.MagicMock()
        mock_domains = (mocker.MagicMock(),) if has_domains else ()
        mock_finder.find_potential_domains.return_value = mock_domains
        mocker.patch(f'{_MODULE}.DomainFinder', return_value=mock_finder)
        # pylint: disable-next=protected-access           # OK in tests
        explorer._find_potential_domain_subfolders()
        mock_finder.find_potential_domains.assert_called_once()
        assert explorer.has_domains == has_domains

    _remove_files = {
        '_remove_ori_files': [f'{file}{ORI_SUFFIX}' for file in STATE_FILES],
        '_remove_out_and_supp': (DEFAULT_OUT, DEFAULT_SUPP),
        }

    @parametrize('method_name,files', _remove_files.items(), ids=_remove_files)
    @patch_discard
    # pylint: disable-next=too-many-arguments  # 1/6 patch, 2/6 fixture
    def test_files_discarded(self, discard_files,
                             method_name, files,
                             explorer, mocker):
        """Check that files and folders are removed."""
        method = getattr(explorer, method_name)
        discard_files.return_value = expect = mocker.MagicMock()
        returned = method()
        removed = (explorer.path / f for f in files)
        discard_files.assert_called_once_with(*removed)
        assert returned is expect

    def test_infer_run_info(self, explorer, mocker):
        """Check correct result of inferring info from log files."""
        mocker.patch.object(explorer, '_logs')
        mocker.patch.object(explorer.logs,
                            'infer_run_info',
                            return_value=_TEST_STRING)
        assert explorer.infer_run_info() is _TEST_STRING

    _archive_files = {
        False: (),
        True: ('file1', 'file2'),
        }

    @parametrize('expect,files', _archive_files.items(), ids=_archive_files)
    def test_needs_archiving(self, explorer, files, expect):
        """Check the needs_archiving property."""
        # pylint: disable-next=protected-access           # OK in tests
        explorer._files_to_archive = files
        assert explorer.needs_archiving == expect

    def test_read_and_clear_notes_file(self, explorer, patched_path):
        """Test successful reading and clearing of a notes file."""
        notes = explorer.read_and_clear_notes_file()
        assert notes == _TEST_STRING
        patched_path.read.assert_called_once_with(encoding='utf-8')
        patched_path.write.assert_called_once_with('', encoding='utf-8')

    def test_read_and_clear_notes_file_no_notes(self, explorer, patched_path):
        """Test result of reading notes whn no file exists."""
        patched_path.glob.return_value = iter(())
        notes = explorer.read_and_clear_notes_file()
        assert not notes
        patched_path.read.assert_not_called()
        patched_path.write.assert_not_called()

    def test_remove_tensors_and_deltas(self, explorer, mocker):
        """Check delegation of calls in remove_tensors_and_deltas."""
        mock_discard = mocker.patch.object(explorer.tensors, 'discard')
        mock_folders = mocker.MagicMock()
        explorer.remove_tensors_and_deltas(mock_folders)
        mock_discard.assert_called_once_with(explorer.history, mock_folders)

    def test_replace_state_files_from_ori_no_ori(self, explorer):
        """Check there are no complaints if an _ori file is not there."""
        with make_obj_raise('pathlib.Path.replace', FileNotFoundError):
            with not_raises(FileNotFoundError):
                # pylint: disable-next=protected-access   # OK in tests
                explorer._replace_state_files_from_ori()

    _fake_ori_files = {
        'no files': (),
        'one file': ((Path('dst'), Path('src')),)
        }

    @parametrize(files=_fake_ori_files.values(), ids=_fake_ori_files)
    def test_replace_state_files_from_ori_return(self, files,
                                                 explorer, mocker):
        """Check the expected return value of _replace_state_files_from_ori."""
        mocker.patch.object(explorer, 'list_files_to_replace',
                            return_value=files)
        mocker.patch('pathlib.Path.replace')
        # pylint: disable-next=protected-access           # OK in tests
        returned = explorer._replace_state_files_from_ori()
        assert returned == bool(files)


class TestRootExplorerCopyStateFilesFrom:
    """Tests for the _copy_state_files_from method."""

    copy_info = namedtuple('copy_info', ('missing', 'fail'))
    copy_info_domain = namedtuple('copy_info',
                                  ('missing', 'fail', 'no_complain'))

    @fixture(name='call_copy')
    def fixture_call_copy(self, src, explorer, mock):
        """Call _copy_state_files_from, return expected failures."""
        def _call(missing, fail, only_files):
            failures = mock(missing, fail, only_files)
            context = pytest.raises if failures else not_raises
            with context(FileOperationFailedError) as exc_info:
                # pylint: disable-next=protected-access
                explorer._copy_state_files_from(src, '{}', '{}_fmt',
                                                only_files=only_files)
            return failures, exc_info
        return _call

    @fixture(name='make_expected_failures')
    def fixture_make_expected_failures(self, src, oserror):
        """Return a dict of failures for copying state files."""
        def _make(missing, fail, only_files):
            failures = {}
            for file in (*missing, *fail):
                if only_files and file not in only_files:
                    continue
                if not only_files and file not in MOCK_STATE_FILES:
                    continue
                failures[file] = (oserror if file in fail
                                  else [src/file, src/f'{file}_fmt'])
            return failures
        return _make

    @fixture(name='mock')
    def fixture_mock_implementation(self,
                                    oserror,
                                    make_expected_failures,
                                    mocker):
        """Replace implementation details with mocks."""
        def _mock(missing, fail, only_files):
            def is_file(path):
                return not path.name.startswith(missing)
            def copy2(_, dst):
                if dst.name in fail:
                    raise oserror
            def _glob(path, pattern):
                globbed = next(f for f in only_files or MOCK_STATE_FILES
                               if pattern.startswith(f))
                return (path/globbed,)
            mocker.patch('pathlib.Path.is_file', is_file)
            mocker.patch('pathlib.Path.glob', _glob)
            mocker.patch('shutil.copy2', copy2)
            return make_expected_failures(missing, fail, only_files)
        return _mock

    @fixture(name='oserror')
    def fixture_oserror(self):
        """Return an OSError instance."""
        return OSError()

    @fixture(name='src')
    def fixture_src(self):
        """Return a non-existing path."""
        return Path.cwd()/'does_not_exist'

    _copy_from = {  # missing, copy failures
        'successful': copy_info((), ()),
        'not found': copy_info(('fake', 'PARAMETERS'), ()),
        'fail to copy': copy_info((), ('fake',)),
        }

    @parametrize(info=_copy_from.values(), ids=_copy_from)
    @parametrize(only_files=(None, ('fake',)))
    def test_copy_state_files_from(self, info, only_files, call_copy):
        """Test the _copy_state_files_from method."""
        failures, exc_info = call_copy(*info, only_files)
        if failures:
            assert exc_info.value.failures == failures

    _copy_from_domain_main = {  # missing, copy failures, no complaints
        'successful': copy_info_domain((), (), ()),
        'not found': copy_info_domain(('PARAMETERS', 'POSCAR'),
                                      (),
                                      ('POSCAR',)),
        # VIBROCC should not be there, but if it is and
        # we fail to copy it, we should still complain
        'fail to copy': copy_info_domain((), ('VIBROCC',), ()),
        }

    @parametrize(info=_copy_from_domain_main.values(),
                 ids=_copy_from_domain_main)
    def test_copy_state_files_from_domains(self,
                                           explorer,
                                           info,
                                           call_copy):
        """Test the _copy_state_files_from method for the main domain root."""
        explorer.has_domains = True
        failures, exc_info = call_copy(info.missing, info.fail, None)
        for file in info.no_complain:
            try:
                del failures[file]
            except KeyError:
                pass
        if failures:
            assert exc_info.value.failures == failures

    def test_log_failures_exception(self, explorer, caplog):
        """Check logging emission when file copy raises exceptions."""
        exc_txt = 'test exception'
        failures = {'file': Exception(exc_txt)}
        # pylint: disable-next=protected-access           # OK in tests
        explorer._log_failures_when_copying_input_files(failures)
        assert exc_txt in caplog.text

    def test_log_failures_not_found(self, explorer, caplog):
        """Check logging emission when file copy fails because not found."""
        rel_paths = ['try_one', 'try_folder/try_file']
        failures = {'file': [explorer.path/f for f in rel_paths]}
        # pylint: disable-next=protected-access           # OK in tests
        explorer._log_failures_when_copying_input_files(failures)
        log = caplog.text
        assert all(p in log for p in rel_paths)


class TestRootExplorerNextCalc:
    """Tests for prepare_for_next_calc_run and related methods."""

    # The methods that are called in prepare_for_next_calc_run
    called = (
        '_mark_state_files_as_ori',
        '_copy_state_files_from_out_or_original_inputs',
        'complain_about_edited_files',
        )

    def test_methods_called(self, explorer, mocker):
        """Check that prepare_for_next_calc_run calls the expected methods."""
        mock_called = [mocker.patch.object(explorer, m) for m in self.called]
        explorer.prepare_for_next_calc_run()
        for method in mock_called:
            method.assert_called_once()

    # Methods whose exceptions are delayed to
    # the end of prepare_for_next_calc_run
    errors_delayed = (
        '_mark_state_files_as_ori',
        '_copy_state_files_from_out_or_original_inputs',
        )

    @parametrize(method=errors_delayed)
    def test_errors_delayed(self, method, explorer, mocker):
        """Check failure when one method fails."""
        mock_called = [mocker.patch.object(explorer, m) for m in self.called]
        exc = FileOperationFailedError({'fail_file': 'fail_reason'})
        getattr(explorer, method).side_effect = exc
        with pytest.raises(OSError) as exc_info:
            explorer.prepare_for_next_calc_run()
        assert exc_info.match('fail_file: fail_reason')
        # All methods called even if there was an exception
        for called_method in mock_called:
            called_method.assert_called()

    @parametrize(missing=(True,False))
    def test_mark_state_files_as_ori(self, missing, explorer, mocker):
        """Check correct suffixing of root inputs as _ori."""
        mock_replace = mocker.patch(
            'pathlib.Path.replace',
            side_effect=FileNotFoundError if missing else None
            )
        # pylint: disable-next=protected-access           # OK in tests
        explorer._mark_state_files_as_ori()
        assert mock_replace.call_count == len(MOCK_STATE_FILES)

    def test_mark_state_files_as_ori_fails(self, explorer, mocker, caplog):
        """Check no complaints if root files to rename are missing."""
        exc = OSError('test')
        mock_replace = mocker.patch('pathlib.Path.replace', side_effect=exc)
        with pytest.raises(FileOperationFailedError) as exc_info:
            # pylint: disable-next=protected-access       # OK in tests
            explorer._mark_state_files_as_ori()
        assert mock_replace.call_count == len(MOCK_STATE_FILES)
        assert exc_info.value.failures == dict.fromkeys(MOCK_STATE_FILES, exc)
        assert all(f in caplog.text for f in MOCK_STATE_FILES)

    class EqualExc(Exception):
        """An exception with equality."""

        def __eq__(self, other):
            return type(self) is type(other) and self.args == other.args

    _copy_excs = {
        # Keys are: Exceptions on attempts from OUT/original_inputs
        # (or None) and the "combined" failures from both
        'OUT exists': (None, None, None),
        'original_inputs exists': (
            FileOperationFailedError({'fake_file': ['not found in OUT']}),
            None,
            None,
            ),
        'file missing': (
            FileOperationFailedError({'fake_file': ['not found in OUT']}),
            FileOperationFailedError({'fake_file': ['not found in orig']}),
            {'fake_file': ['not found in OUT', 'not found in orig']},
            ),
        'exc in OUT': (
            FileOperationFailedError({'fake_file': EqualExc('failed in OUT')}),
            FileOperationFailedError({'fake_file': ['not found in orig']}),
            {'fake_file': EqualExc('failed in OUT')},
            ),
        'exc in original_inputs': (
            FileOperationFailedError({'fake_file': ['not found in OUT']}),
            FileOperationFailedError({'fake_file': EqualExc('failed in ori')}),
            {'fake_file': EqualExc('failed in ori')},
            ),
        }

    @parametrize('out_exc,ori_exc,fails', _copy_excs.values(), ids=_copy_excs)
    # pylint: disable-next=too-many-arguments  # 2 fixtures
    def test_copy_state_files(self, out_exc, ori_exc, fails,
                              explorer, mocker):
        """Test failure to copy new inputs from OUT/original_inputs."""
        def file_missing(src, *_, **__):
            """Fake a failed copy from src."""
            if out_exc and src.name == DEFAULT_OUT:
                raise out_exc
            if ori_exc:
                raise ori_exc
        mock_copy_from = mocker.patch.object(explorer,
                                             '_copy_state_files_from',
                                             side_effect=file_missing)
        # Patch out logging, as it requires paths in the failures
        mock_log = mocker.patch.object(
            explorer,
            '_log_failures_when_copying_input_files'
            )
        context = pytest.raises if fails else not_raises
        with context(FileOperationFailedError) as exc_info:
            # pylint: disable-next=protected-access       # OK in tests
            explorer._copy_state_files_from_out_or_original_inputs()
        expect_calls = [
            mocker.call(explorer.path / DEFAULT_OUT, '{}', '{}_OUT_[0-9]*'),
            mocker.call(explorer.orig_inputs_dir,
                        only_files=out_exc.failures if out_exc else {}),
            ]
        assert mock_copy_from.mock_calls == expect_calls
        if fails:
            assert exc_info.value.failures == fails
            mock_log.assert_called_once_with(fails)


class TestRootExplorerListFiles:
    """Tests for listing files/directories to be removed."""

    def test_to_discard_no_files(self, explorer, mocker):
        """Test list_paths_to_discard when no discardable files exist."""
        explorer.collect_info()
        mocker.patch.object(Path, 'exists', return_value=False)
        paths_to_discard = explorer.list_paths_to_discard()
        assert not any(paths_to_discard)

    def test_to_discard_files_exist(self, explorer, mocker, mock_history):
        """Test list_paths_to_discard when some discardable files exist."""
        mock_log_path = explorer.path / 'log_file.log'
        mock_files = (
            mock_log_path,
            explorer.path / DEFAULT_TENSORS / f'{DEFAULT_TENSORS}_001.zip',
            )
        mock_paths = (
            explorer.path / DEFAULT_OUT,  # No SUPP
            *mock_files
            )
        explorer.collect_info()
        mocker.patch.object(Path, 'exists', new=mock_exists(mock_paths))
        mocker.patch.object(Path, 'is_file', new=mock_exists(mock_files))
        mocker.patch.object(type(explorer.logs), 'files', [mock_log_path])
        mocker.patch.object(explorer, 'history', mock_history({1: 1}))
        paths_to_discard = explorer.list_paths_to_discard()
        assert paths_to_discard == mock_paths

    def test_to_replace_no_ori(self, explorer, mocker):
        """Test list_files_to_replace when no '_ori' files are present."""
        mocker.patch.object(Path, 'exists', return_value=False)
        files_to_replace = explorer.list_files_to_replace()
        assert not any(files_to_replace)

    _ori_files = {}
    for repeats in range(len(STATE_FILES)):
        _ori_files.update({'+'.join(files): files
                           for files in combinations(STATE_FILES, repeats+1)})

    @parametrize(state_files=_ori_files.values(), ids=_ori_files)
    def test_to_replace_ori_files_exist(self, state_files, explorer, mocker):
        """Test list_files_to_replace when '_ori' files are present."""
        mock_ori = [explorer.path / f'{f}{ORI_SUFFIX}' for f in state_files]
        mocker.patch.object(Path, 'exists', new=mock_exists(mock_ori))
        files_to_replace = explorer.list_files_to_replace()
        expected_files_to_replace = {
            explorer.path / f: ori_f
            for f, ori_f in zip(state_files, mock_ori)
            }
        assert dict(files_to_replace) == expected_files_to_replace


class TestRootExplorerMarkEditedFiles:
    """Tests for the mark_edited_files method of RootExplorer."""

    @fixture(name='make_calc_inputs')
    def factory_calc_inputs(self):
        """Fill explorer.path with calc input files."""
        def _make(explorer):
            filesystem_from_dict(dict.fromkeys(MOCK_STATE_FILES), explorer.path)
        return _make

    @fixture(name='make_calc_log')
    def factory_calc_log(self):
        """Add a calc log file to explorer.path."""
        def _make(explorer):
            timestamp = DateTimeFormat.FILE_SUFFIX.now()
            (explorer.path / f'{LOG_PREFIX}-{timestamp}.log').touch()
        return _make

    @fixture(name='mock_edited_after_calc')
    def factory_edited_after_calc(self,
                                  make_calc_log,
                                  make_calc_inputs,
                                  explorer):
        """Prepare a root tree with input files edited after calc started."""
        make_calc_log(explorer)
        time.sleep(1)  # NB: need >= 1 s, as timestamp has no micros
        make_calc_inputs(explorer)

    @fixture(name='explorer')  # NB: replaces the module fixture!
    def fixture_explorer(self, tmp_path, mocker):
        """Return a RootExplorer instance for a real temporary path."""
        return RootExplorer(tmp_path, bookkeeper=mocker.MagicMock())

    @staticmethod
    def check_has_edited_files(explorer, expected):
        """Check that the expected _edited files exist in explorer.path."""
        edited = set(f.name.replace(EDITED_SUFFIX, '')
                     for f in explorer.path.glob(f'*{EDITED_SUFFIX}*'))
        assert edited == set(expected)

    @staticmethod
    def check_warns_edited(log_text, expected):
        """Check that all expected edited files are named in log_text."""
        assert all(f+EDITED_SUFFIX in log_text for f in expected)

    def test_complain_edited(self, explorer, mocker, caplog):
        """Check warnings when edited files are present."""
        mocker.patch('pathlib.Path.glob',
                     return_value=(explorer.path/f'file{EDITED_SUFFIX}',))
        explorer.complain_about_edited_files()
        assert caplog.text

    def test_no_log_file(self, explorer, make_calc_inputs):
        """Check that no files are _edited if there's no calc log file."""
        make_calc_inputs(explorer)
        explorer.collect_info()
        explorer.mark_edited_files()
        self.check_has_edited_files(explorer, ())

    def test_earlier(self, explorer, make_calc_inputs, make_calc_log):
        """Check that no files are _edited when they predate calc."""
        # Prepare contents of the root directory
        make_calc_inputs(explorer)
        time.sleep(1)  # NB: need >= 1 s, as timestamp has no micros
        make_calc_log(explorer)

        # Collect info, then ensure no files are _edited due to sleep
        explorer.collect_info()
        explorer.mark_edited_files()
        self.check_has_edited_files(explorer, ())

    @pytest.mark.usefixtures('mock_edited_after_calc')
    def test_later_no_ori(self, explorer, caplog, mocker):
        """Check _edited files when they follow calc and there's no ori."""
        # Patch pathlib.Path.replace to test a replacement failure
        failed, *expect_edited = MOCK_STATE_FILES
        old_replace = Path.replace
        def _mock_replace(path_, other_):
            """A callable to replace pathlib.Path.replace."""
            if path_.name != failed:
                old_replace(path_, other_)
            raise OSError('failed to replace')
        mocker.patch('pathlib.Path.replace', _mock_replace)

        # Collect info, then ensure expect_edited files are _edited
        explorer.collect_info()
        explorer.mark_edited_files()
        self.check_has_edited_files(explorer, expect_edited)
        self.check_warns_edited(caplog.text, MOCK_STATE_FILES)
        assert f'Failed to rename {failed}' in caplog.text

    @pytest.mark.usefixtures('mock_edited_after_calc')
    def test_later_with_ori(self, explorer, mocker, caplog):
        """Check _edited files when they follow calc."""
        # Add original_inputs files
        same_contents, *expect_edited = MOCK_STATE_FILES

        def _mock_same_contents(_, ori_file):
            """A replacement for file_contents_identical."""
            return ori_file.name == same_contents

        mocker.patch(f'{_MODULE}.file_contents_identical', _mock_same_contents)

        # Ensure only files with different contents are _edited
        explorer.collect_info()
        explorer.mark_edited_files()
        self.check_has_edited_files(explorer, expect_edited)
        self.check_warns_edited(caplog.text, expect_edited)
        assert same_contents not in caplog.text  # No warning if same

    @pytest.mark.usefixtures('mock_edited_after_calc')
    def test_root_file_missing(self, explorer, caplog):
        """Check no _edited version is created for a missing input file."""
        # Delete one file
        missing, *expect_edited = MOCK_STATE_FILES
        (explorer.path / missing).unlink()
        explorer.collect_info()
        explorer.mark_edited_files()
        self.check_has_edited_files(explorer, expect_edited)
        self.check_warns_edited(caplog.text, expect_edited)
        assert missing not in caplog.text  # No warning if missing


class TestRootExplorerRaises:
    """Tests for conditions that causes exceptions or complaints."""

    def test_cannot_replace_state_files_from_ori(self, explorer, mocker):
        """Check complaints when a root file cannot be replaced with _ori."""
        mock_error = mocker.patch(f'{_MODULE}.LOGGER.error')
        with raises_exception('pathlib.Path.replace', OSError):
            mocker.patch('pathlib.Path.exists', return_value=True)
            # pylint: disable-next=protected-access       # OK in tests
            explorer._replace_state_files_from_ori()
            mock_error.assert_called_once()

    _read_notes_raises = {
        'read': '',
        'write': _TEST_STRING,
        }

    @parametrize('method,expect',
                 _read_notes_raises.items(),
                 ids=_read_notes_raises)
    # pylint: disable-next=too-many-arguments  # 3/6 are fixtures
    def test_read_and_clear_notes_file(self, method, expect, explorer,
                                       patched_path, mocker):
        """Test that complaints by read/write_text emit errors."""
        complaining = getattr(patched_path, method)
        complaining.side_effect = OSError
        mock_error = mocker.patch(f'{_MODULE}.LOGGER.error')
        notes = explorer.read_and_clear_notes_file()
        assert notes == expect
        mock_error.assert_called_once()

    _attr_needs_update = (
        'calc_timestamp',
        'needs_archiving',
        )
    _method_needs_update = (  # Only those without args
        'clear_for_next_calc_run',
        'infer_run_info',
        'list_paths_to_discard',
        'mark_edited_files',
        'revert_to_previous_calc_run',
        '_collect_files_to_archive',
        )

    @parametrize(attr=_attr_needs_update)
    def test_too_early_attribute_access(self, explorer, attr):
        """Check that accessing attributes before collection fails."""
        with pytest.raises(AttributeError,
                           match=r'.*collect_(info|subfolders).*'):
            attrgetter(attr)(explorer)

    @parametrize(method_name=_method_needs_update)
    def test_too_early_method_call(self, explorer, method_name):
        """Check that calling methods before collection fails."""
        method = attrgetter(method_name)(explorer)
        with pytest.raises(AttributeError,
                           match=r'.*collect_(info|subfolders).*'):
            method()


class TestRootExplorerUnlabledFiles:
    """Tests for the ensure_has_unlabled_inputs method."""

    def test_no_ori_files(self, explorer, mocker):
        """Check no file operations happen if no _ori files exist."""
        mock_copy = mocker.patch('shutil.copy2')
        missing = explorer.ensure_has_unlabled_inputs()
        assert not missing
        mock_copy.assert_not_called()

    def test_oserror(self, explorer, mocker, caplog):
        """Check that failure to copy files emits log messages."""
        mocker.patch('pathlib.Path.glob',
                     return_value=(explorer.path/'input_file_ori',))
        mocker.patch('pathlib.Path.is_file', return_value=False)
        mock_copy = mocker.patch('shutil.copy2', side_effect=OSError)
        missing = explorer.ensure_has_unlabled_inputs()
        assert missing
        mock_copy.assert_called_once()  # One "fake" _ori file
        assert caplog.text

    def test_expected_files(self, tmp_path, mocker):
        """Ensure the expected unlabeled files are created."""
        tree_before = {
            'file_1': 'contents of file_1',
            'file_1_ori': 'file_1 exists already ',
            'file_2_ori': 'file_2 instead is missing',
            # NB: given the current implementation of bookkeeper there
            # should never be both _ori and _edited files. However,
            # better be sure that we handle correctly this case too.
            'file_3_edited': 'contents of the edited version of file_3',
            'file_3_ori': 'file_3 has been edited instead',
            }
        expect_tree = {
            **tree_before,
            'file_2': tree_before['file_2_ori'],  # Only file_2 created
            }
        filesystem_from_dict(tree_before, tmp_path)
        explorer = RootExplorer(tmp_path, mocker.MagicMock())
        missing = explorer.ensure_has_unlabled_inputs()
        assert missing
        tree_after = filesystem_to_dict(tmp_path)
        assert tree_after == expect_tree


class TestDomainRootExplorer:
    """Tests for the DomainRootExplorer subclass of RootExplorer."""

    @fixture(name='domain')
    def fixture_domain(self, explorer):
        """Return an initialized DomainRootExplorer."""
        return DomainRootExplorer(explorer.path/'domain',
                                  explorer.workhistory.bookkeeper,
                                  explorer)

    def test_calc_timestamp_from_domain(self, domain, explorer, mocker):
        """Check pulling of the calc_timestamp property from root logs."""
        most_recent_log_in_domain = mocker.MagicMock()
        # pylint: disable-next=protected-access           # OK in tests
        explorer._logs = mocker.MagicMock(most_recent=None)
        # pylint: disable-next=protected-access           # OK in tests
        domain._logs = mocker.MagicMock(most_recent=most_recent_log_in_domain)
        assert domain.calc_timestamp is most_recent_log_in_domain.timestamp

    def test_calc_timestamp_from_main(self, domain, explorer, mocker):
        """Check pulling of the calc_timestamp property from main."""
        most_recent_log_in_main = mocker.MagicMock()
        # pylint: disable-next=protected-access           # OK in tests
        explorer._logs = mocker.MagicMock(most_recent=most_recent_log_in_main)
        # pylint: disable-next=protected-access           # OK in tests
        domain._logs = mocker.MagicMock(most_recent=None)
        assert domain.calc_timestamp is most_recent_log_in_main.timestamp

    def test_infer_run_info_from_domain(self, domain, explorer, mocker):
        """Check inferring run info from the log file in the domain folder."""
        for root in (domain, explorer):
            mocker.patch.object(root, '_logs')
        domain_info = mocker.MagicMock()
        explorer.logs.infer_run_info.return_value = mocker.MagicMock()
        domain.logs.infer_run_info.return_value = domain_info
        assert domain.infer_run_info() is domain_info

    def test_infer_run_info_from_main(self, domain, explorer, mocker):
        """Check inferring run info from the main log file."""
        for root in (domain, explorer):
            mocker.patch.object(root, '_logs')
        main_info = mocker.MagicMock()
        explorer.logs.infer_run_info.return_value = main_info
        domain.logs.infer_run_info.return_value = {}
        assert domain.infer_run_info() is main_info

    def test_relative_path_from_domain(self, domain):
        """Check the expected return of _relative_path."""
        # pylint: disable-next=protected-access           # OK in tests
        domain._path = Path('/some/other/path')
        # pylint: disable-next=protected-access           # OK in tests
        result = domain._relative_path(domain.path/'abc')
        assert result == Path('abc').as_posix()

    def test_relative_path_from_main(self, domain):
        """Check the expected return of _relative_path."""
        # pylint: disable-next=protected-access           # OK in tests
        result = domain._relative_path(domain.path/'abc')
        assert result == Path('domain/abc').as_posix()


class TestTensorsAndDeltasInfo:
    """Tests for the TensorsAndDeltasInfo class."""

    tensor_index = 5  # Example index to use for tensor/delta files

    @fixture(name='mock_get_tensors', autouse=True)
    def patch_get_tensors(self, mocker):
        """Replace getMaxTensorIndex with a fake one."""
        return mocker.patch(f'{_MODULE}.getMaxTensorIndex',
                            return_value=self.tensor_index)

    @fixture(name='simple_history')
    def fixture_simple_history(self, mock_history):
        """Return a fake HistoryExplorer with one tensor and one run."""
        return mock_history({self.tensor_index: 1})

    def test_collect(self, tensors, mock_get_tensors):
        """Test that collect sets the correct tensor index."""
        tensors.collect()
        # pylint: disable-next=protected-access       # OK in tests
        mock_get_tensors.assert_called_once_with(home=tensors._root,
                                                 zip_only=True)
        assert tensors.most_recent == self.tensor_index

    def test_to_discard(self, tensors, simple_history, mocker):
        """Test list_paths_to_discard when tensor and delta files exist."""
        history = simple_history
        history.last_folder_and_siblings.extend(
            # Simulate some siblings too
            mocker.MagicMock(tensor_num=self.tensor_index)
            for _ in range(4)
            )
        expected = (
            f'{DEFAULT_DELTAS}/{DEFAULT_DELTAS}_{self.tensor_index:03d}.zip',
            f'{DEFAULT_TENSORS}/{DEFAULT_TENSORS}_{self.tensor_index:03d}.zip',
            )
        # pylint: disable-next=protected-access           # OK in tests
        expected = tuple(tensors._root / f for f in expected)
        tensors.collect()
        mocker.patch.object(Path, 'is_file', new=mock_exists(expected))
        discard_paths = tensors.list_paths_to_discard(history)
        assert discard_paths == expected

    def test_to_discard_none(self, tensors, simple_history, mocker):
        """Test list_paths_to_discard when no tensor/delta files exist."""
        tensors.collect()
        mocker.patch.object(Path, 'is_file', return_value=False)
        assert not any(tensors.list_paths_to_discard(simple_history))

    def test_to_discard_history_empty(self, tensors, mock_history, mocker):
        """Check that no Tensors/Deltas are discarded with empty history."""
        tensors.collect()
        mocker.patch.object(Path, 'is_file', return_value=True)
        history = mock_history({})
        assert not any(tensors.list_paths_to_discard(history))

    def test_to_discard_older_runs(self, tensors, mock_history, mocker):
        """Check that no Tensors/Deltas are discarded when more runs exist."""
        tensors.collect()
        mocker.patch.object(Path, 'is_file', return_value=True)
        history = mock_history({self.tensor_index: 3})
        assert not any(tensors.list_paths_to_discard(history))

    @patch_discard
    def test_discard(self, mock_discard_files,
                     tensors, simple_history, mocker):
        """Test removal of tensors and deltas."""
        removed = (
            f'{DEFAULT_DELTAS}/{DEFAULT_DELTAS}_{self.tensor_index:03d}.zip',
            f'{DEFAULT_TENSORS}/{DEFAULT_TENSORS}_{self.tensor_index:03d}.zip',
            )
        # pylint: disable-next=protected-access           # OK in tests
        removed = tuple(tensors._root/f for f in removed)
        tensors.collect()
        mocker.patch.object(Path, 'is_file', return_value=True)
        tensors.discard(simple_history,
                        simple_history.last_folder_and_siblings)
        mock_discard_files.assert_called_once_with(*removed)

    _attr_needs_update = (
        'most_recent',
        )
    _method_needs_update = (
        # None for now
        )

    @parametrize(attr=_attr_needs_update)
    def test_too_early_attribute_access(self, tensors, attr):
        """Check that accessing attributes before collection fails."""
        with pytest.raises(AttributeError, match=r'.*collect.*'):
            attrgetter(attr)(tensors)

    @parametrize(method_name=_method_needs_update)
    def test_too_early_method_call(self, tensors, method_name):
        """Check that calling methods before collection fails."""
        method = attrgetter(method_name)(tensors)
        with pytest.raises(AttributeError, match=r'.*collect.*'):
            method()
