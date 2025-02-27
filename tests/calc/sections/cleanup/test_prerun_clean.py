"""Tests for cleanup.prerun_clean (and helpers) of viperleed.calc.section."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2025-02-03'
__license__ = 'GPLv3+'

from pathlib import Path
import shutil

import pytest
from pytest_cases import fixture
from pytest_cases import parametrize

from viperleed.calc.constants import DEFAULT_OUT
from viperleed.calc.constants import DEFAULT_SUPP
from viperleed.calc.constants import DEFAULT_WORK_HISTORY
from viperleed.calc.classes.rparams.domain_params import DomainParameters
from viperleed.calc.classes.rparams.rparams import Rparams
from viperleed.calc.lib.context import execute_in_dir
from viperleed.calc.sections.calc_section import ALL_INPUT_FILES
from viperleed.calc.sections.cleanup import move_oldruns
from viperleed.calc.sections.cleanup import prerun_clean

from ....helpers import filesystem_from_dict
from ....helpers import filesystem_to_dict
from .conftest import _MODULE


class TestPrerunClean:
    """Tests for the prerun_clean function."""

    executables = {  # Some compiled executables
        'refcalc-010203-040506': 'refcalc exe contents',
        'rfactor-999999-888888': 'rfactor exe contents',
        'search-000001-000002': 'search exe contents',
        'superpos-123456-123456': 'superpos exe contents',
        # And also their windows versions
        'refcalc-010203-040506.exe': 'refcalc exe windows contents',
        'rfactor-999999-888888.exe': 'rfactor exe windows contents',
        'search-000001-000002.exe': 'search exe windows contents',
        'superpos-123456-123456.exe': 'superpos windows exe contents',
        }
    old_out_files = {
        # Some old-style _OUT files
        'POSCAR_OUT': 'POSCAR_OUT contents',
        'PARAMETERS_OUT': 'PARAMETERS_OUT contents',
        'VIBROCC_OUT': 'VIBROCC_OUT contents',
        'R_OUT_refcalc_R=0.123': 'R_OUT contents',
        }
    silent_fail = {  # Those for which we don't complain
        'fortran-compile.log': '',
        }
    calc_tree = {
        DEFAULT_WORK_HISTORY: {'subfolder':{}},
        DEFAULT_OUT: {'out_file': 'out file contents'},
        DEFAULT_SUPP: {'supp_file': 'supp file contents'},
        **old_out_files,
        **executables,
        **silent_fail,
        }

    def check_fails_to_remove_file(self, file, run, mocker):
        """Check failure to remove file."""
        _path_unlink = Path.unlink
        def _fails_to_remove(path):
            if path.name == file:
                raise OSError
            _path_unlink(path)
        mocker.patch('pathlib.Path.unlink', new=_fails_to_remove)
        expect, clean = run()
        # The failing file should still be there
        expect[file] = self.calc_tree[file]
        assert clean == expect

    @fixture(name='run')
    def fixture_run_prerun_clean(self, tmp_path, rpars, work_tree):
        """Execute prerun_clean at tmp_path."""
        def _run(nothing_to_move=False):
            expect, log_name, mock_move_oldruns = work_tree
            prerun_clean(rpars, logname=log_name)
            clean = filesystem_to_dict(tmp_path)
            if nothing_to_move:
                mock_move_oldruns.assert_not_called()
            else:
                kwargs = {'prerun': True, }#'recurse_domains': False}
                mock_move_oldruns.assert_called_once_with(rpars, **kwargs)
            return expect, clean
        return _run

    @fixture(name='work_tree')
    def fixture_work_tree(self, tmp_path, mocker):
        """Prepare a directory with leftovers from a calc run."""
        # Mock move_oldruns to avoid having
        # to deal with its implementation
        mock_move = mocker.patch(f'{_MODULE}.move_oldruns')
        # Input files should be untouched, as is the current log
        surviving = {f: f'{f} contents' for f in ALL_INPUT_FILES}
        log_name = 'current_log.log'
        surviving[log_name] = 'current log contents'
        # Since we mock out move_oldruns, the old log
        # is not moved to a new workhistory subfolder
        surviving['old_log.log'] = 'old log contents'
        tree = {**self.calc_tree, **surviving}
        filesystem_from_dict(tree, tmp_path)
        with execute_in_dir(tmp_path):
            yield surviving, log_name, mock_move

    def test_dont_catch_move_oldruns_bugs(self, run, work_tree):
        """Check that no generic exception is masked."""
        *_, mock_move = work_tree
        mock_move.side_effect = Exception
        with pytest.raises(Exception):
            run()

    def test_fails_to_move_oldruns(self, run, work_tree, caplog):
        """Check warnings when failing to sort previous executions."""
        *_, mock_move = work_tree
        mock_move.side_effect = OSError
        expect, clean = run()
        expect_log = 'old files may be lost'
        assert clean == expect
        assert expect_log in caplog.text

    @parametrize(file=executables)
    def test_fails_to_remove_executable(self, file, run, caplog, mocker):
        """Check warnings when failing to remove a compiled executable."""
        caplog.set_level(0)  # We emit debug
        self.check_fails_to_remove_file(file, run, mocker)
        expect_log = f'Failed to delete file {file}'
        assert expect_log in caplog.text

    @parametrize(file=silent_fail)
    def test_fails_to_remove_log(self, file, run, caplog, mocker):
        """Check there no complaints when failing to remove some log files."""
        caplog.set_level(0)  # Collect all messages
        self.check_fails_to_remove_file(file, run, mocker)
        assert not caplog.text

    @parametrize(file=old_out_files)
    def test_fails_to_remove_old_out_file(self, file, run, caplog, mocker):
        """Check warnings when failing to remove an old-style _OUT file."""
        self.check_fails_to_remove_file(file, run, mocker)
        expect_log = f'Failed to delete previous {file}'
        assert expect_log in caplog.text

    @parametrize(folder=(DEFAULT_WORK_HISTORY, DEFAULT_OUT, DEFAULT_SUPP))
    def test_fails_to_rmdir(self, folder, run, caplog, mocker):
        """Check complaints when deleting a folder fails."""
        _rmtree = shutil.rmtree
        def _rmtree_raises(path):
            if Path(path).name == folder:
                raise OSError
            _rmtree(path)
        mocker.patch('shutil.rmtree', new=_rmtree_raises)
        expect, clean = run()
        # The stuff that failed to be removed should stay
        expect[folder] = self.calc_tree[folder]
        assert clean == expect
        expect_log = ('Failed to clear', folder)
        assert all(m in caplog.text for m in expect_log)

    def test_no_need_to_move_oldruns(self, run, work_tree, tmp_path, caplog):
        """Check no complaints if there's no old stuff to archive."""
        expect, *_ = work_tree
        old_logs = ('old_log.log', *self.silent_fail)
        for log in old_logs:
            (tmp_path/log).unlink()
            try:
                del expect[log]
            except KeyError:
                pass
        _, clean = run(nothing_to_move=True)
        assert clean == expect
        assert not caplog.text

    def test_not_an_executable(self, run, tmp_path, caplog):
        """Check that non-executable files are retained."""
        not_an_exe = tmp_path/'refcalc-not_a_time'
        not_an_exe.touch()
        expect, clean = run()
        expect[not_an_exe.name] = ''
        assert clean == expect
        assert not caplog.text

    @parametrize(folder=(DEFAULT_WORK_HISTORY, DEFAULT_OUT, DEFAULT_SUPP))
    # pylint: disable-next=too-many-arguments  # 4/6 are fixtures
    def test_skips_missing_folder(self, folder, work_tree,
                                  tmp_path, run, caplog):
        """Check no complaints when folder is not found."""
        expect, *_ = work_tree
        shutil.rmtree(tmp_path/folder)
        _, clean = run()
        assert expect == clean
        assert not caplog.text

    @parametrize(file=silent_fail)
    # pylint: disable-next=too-many-arguments  # 4/6 are fixtures
    def test_skips_missing_log(self, file, work_tree, tmp_path, run, caplog):
        """Check no complaints when folder is not found."""
        expect, *_ = work_tree
        (tmp_path/file).unlink()
        _, clean = run()
        assert expect == clean
        assert not caplog.text

    def test_success(self, run, caplog):
        """Check the successful cleanup of a work directory."""
        expect, clean = run()
        assert clean == expect
        assert not caplog.text


class TestPropagateToDomains:
    """Tests for prerun_clean propagation to nested domains."""

    def make_domains(self, rpars, top_level_domains, root_path, mocker):
        """Create nested domains at root_path, return calls."""
        filesystem_from_dict(top_level_domains, root_path)
        calls = {}
        for domain, subdomains in top_level_domains.items():
            if isinstance(subdomains, str):  # A file
                continue
            workdir = root_path/domain
            domain = DomainParameters(workdir, domain)
            domain.rpars = Rparams()
            domain.rpars.TENSOR_INDEX = 1
            rpars.domainParams.append(domain)
            has_subdomains = any(isinstance(s, dict)
                                 for s in subdomains.values())
            call = mocker.call(domain.rpars,
                               logname='',
                               )#recurse_domains=has_subdomains)
            calls[str(workdir)] = call
            if has_subdomains:
                subcalls = self.make_domains(domain.rpars,
                                             subdomains,
                                             workdir,
                                             mocker)
                calls.update(subcalls)
        return calls

    @fixture(name='domains')
    def fixture_domains(self, rpars, tmp_path, mocker):
        """Create sample nested domains at tmp_path."""
        domain_contents = {'log.log': ''}  # To check workhistory
        nested_domains = {
            'simple_domain': domain_contents,
            'nested_domain': {'subdomain_1': domain_contents,
                              'subdomain_2': domain_contents,
                              **domain_contents},
            }
        nested_calls = self.make_domains(rpars,
                                         nested_domains,
                                         tmp_path,
                                         mocker)
        return nested_calls, nested_domains, tmp_path

    @fixture(name='register_calls')
    def fixture_register_calls(self, mocker):
        """Return all calls to the patched function."""
        def _register_calls(func, calls):
            def _wrap(*args, **kwargs):
                calls.append(mocker.call(*args, **kwargs))
                func(*args, **kwargs)
            return _wrap

        def _patch(func):
            calls = []
            mocker.patch(f'{_MODULE}.{func.__name__}',
                         _register_calls(func, calls))
            return calls
        return _patch

    def _check_workhistory_created(self, root_tree, domains):
        """Check that all domain subfolders contain a non-empty workhistory."""
        for domain, subdomains in domains.items():
            if isinstance(subdomains, str):     # A file
                assert domain not in root_tree  # Should be moved
                continue
            assert domain in root_tree
            assert DEFAULT_WORK_HISTORY in root_tree[domain]
            assert root_tree[domain][DEFAULT_WORK_HISTORY]
            has_subdomains = any(isinstance(s, dict)
                                 for s in subdomains.values())
            if has_subdomains:
                self._check_workhistory_created(root_tree[domain], subdomains)

    def test_propagate_to_domains(self, rpars, register_calls, domains):
        """Check that prerun_clean is correctly propagated to all domains."""
        move_calls = register_calls(move_oldruns)
        _, nested_domains, root = domains
        prerun_clean(rpars)

        # Check that move_oldruns was called correctly
        n_domains = (1     # simple_domain
                     + 1   # nested_domain
                     + 2)  # subdomains of nested_domain
        assert len(move_calls) == n_domains  # Only once per domain
        domain_tree = filesystem_to_dict(root)
        self._check_workhistory_created(domain_tree, nested_domains)
