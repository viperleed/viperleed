"""Module run_bookkeeper_base of tests/calc/bokkeeper/bookkeeper.

Defines the _TestBookkeeperRunBase class that collects basic
functionality for testing execution of the bookkeeper in
various modes.
"""

__authors__ = (
    'Alexander M. Imre (@amimre)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2024-11-15'
__license__ = 'GPLv3+'

import logging
import re

from viperleed.calc.bookkeeper.constants import EDITED_SUFFIX
from viperleed.calc.bookkeeper.constants import ORI_SUFFIX
from viperleed.calc.bookkeeper.domain_finder import DomainFinder
from viperleed.calc.bookkeeper.exit_code import BookkeeperExitCode
from viperleed.calc.bookkeeper.history.constants import HISTORY_INFO_NAME
from viperleed.calc.bookkeeper.history.info import HistoryInfoFile
from viperleed.calc.bookkeeper.history.meta import BookkeeperMetaFile
from viperleed.calc.bookkeeper.history.meta import _METADATA_NAME
from viperleed.calc.bookkeeper.log import BOOKIE_LOGFILE
from viperleed.calc.bookkeeper.mode import BookkeeperMode as Mode
from viperleed.calc.constants import DEFAULT_HISTORY
from viperleed.calc.constants import DEFAULT_OUT
from viperleed.calc.constants import DEFAULT_SUPP

from ....helpers import filesystem_to_dict
from ..conftest import MOCK_INPUT_CONTENT
from ..conftest import MOCK_ORIG_CONTENT
from ..conftest import MOCK_OUT_CONTENT
from ..conftest import MOCK_OUT_SUFFIXED_CONTENT
from ..conftest import MOCK_STATE_FILES
from ..conftest import MOCK_WORKHISTORY


class _TestCheckers:
    """Collection of useful methods for checking bookkeeper behavior."""

    @staticmethod
    def _check_domains_history_info(main_run, domains, skip=()):
        """Ensure consistent information has been saved to history.info."""
        main_bookie, *_ = main_run
        entries = {}
        for root_path in (main_bookie.cwd, *domains):
            info = HistoryInfoFile(root_path)
            info.read()
            entries[root_path] = info.last_entry
        assert all(entries.values())  # All should have a valid entry

        identical_fields = ('run_info', 'timestamp', 'is_discarded')
        for attr in identical_fields:
            attr_values = {getattr(e, attr)
                           for p, e in entries.items()
                           if p not in skip}
            assert len(attr_values) == 1, f'Found inconsistent {attr}'

    @staticmethod
    def _check_domains_registered(main_run, domains_info, skip):
        """Ensure domain information has been registered in metadata."""
        main_bookie, main_history_path, *_ = main_run
        main_meta = BookkeeperMetaFile(main_history_path)
        main_meta.read()

        # Collect hashes of the history folders of the domains. While
        # doing so, check that the main has been correctly registered.
        expect_domain_register = (str(main_bookie.cwd), main_meta.hash_)
        expect_main_register = []
        for domain_path, domain_info in domains_info.items():
            if domain_path in skip:
                continue
            domain_history_path = domain_info['history_run']
            domain_meta = BookkeeperMetaFile(domain_history_path)
            domain_meta.read()
            registered_in_domain = domain_meta.domains['main']
            assert registered_in_domain == expect_domain_register
            expect_main_register.append((domain_path.name, domain_meta.hash_))

        # Finally, check the registered info in the main
        registered_in_main = main_meta.domains['domains']
        assert registered_in_main == tuple(expect_main_register)

    @staticmethod
    def _check_file_contents(path, *expected_contents):
        """Check that the file at `path` has any of the `expected_contents`."""
        assert path.is_file()
        contents = path.read_text()
        assert any(e in contents for e in expected_contents)

    def _check_history_exists(self, bookkeeper, history_run_path, *_):
        """Test that `history_run_path` and history.info exist."""
        assert history_run_path.is_dir()
        assert (bookkeeper.cwd / HISTORY_INFO_NAME).exists()
        self.check_metadata_exists(history_run_path)

    def _check_workhistory_archived(self, bookkeeper, *_):
        """Ensure the workhistory folders have been stored correctly."""
        history = bookkeeper.history.path
        # pylint: disable-next=protected-access           # OK in tests
        workhistory = bookkeeper._workhistory.path
        assert not workhistory.is_dir()
        for ori_name, hist_name in MOCK_WORKHISTORY.items():
            if hist_name is None:  # File should be deleted
                continue
            moved_dir = history/hist_name
            moved_file = moved_dir/'file'
            assert moved_dir.is_dir()
            self._check_file_contents(moved_file, ori_name)

    def _check_input_files_in_history(self, *run):
        """Make sure that input files were stored in history."""
        *_, history_run_path, _ = run
        for file in MOCK_STATE_FILES:
            archived_input = history_run_path / file
            assert archived_input.is_file()
            self._check_file_contents(archived_input, MOCK_ORIG_CONTENT)

    def _check_out_files_in_history(self, *run, **kwargs):
        """Check that the expected state files are stored in 'OUT'."""
        *_, history_run_path, _ = run
        self.check_out_files_at_path(history_run_path, **kwargs)

    @staticmethod
    def check_complained_about_edited(caplog):
        """Make sure there was a log message about finding _edited files."""
        # pylint: disable=magic-value-comparison
        assert 'user-edited' in caplog.text

    @staticmethod
    def check_exit_code_ok(exit_code):
        """Make sure bookkeeper exited with a non-failure condition."""
        assert exit_code is not BookkeeperExitCode.FAIL

    def check_has_archived(self, run, **kwargs):
        """Check that the root directory has been archived to history."""
        self._check_history_exists(*run)
        kwargs['out_suffixed'] = kwargs.pop('has_out_suffixed', False)
        kwargs['missing'] = kwargs.get('missing_out_files', {})
        self._check_out_files_in_history(*run, **kwargs)
        self._check_input_files_in_history(*run)
        if kwargs['mode'] is Mode.ARCHIVE:
            self.check_root_after_archive(*run, **kwargs)
        else:
            self.check_root_is_clean(*run)
        # Check that the workhistory directories are
        # also where they should be
        bookkeeper, *_ = run
        self._check_workhistory_archived(bookkeeper)

    @staticmethod
    def check_history_folder_empty(bookkeeper, *_):
        """Test that an empty (except bookkeeper log) history folder exists."""
        cwd = bookkeeper.cwd
        assert (cwd / DEFAULT_HISTORY).exists()
        history_contents = (f for f in (cwd / DEFAULT_HISTORY).iterdir()
                            if f.name != BOOKIE_LOGFILE)
        assert not any(history_contents)

    @staticmethod
    def check_last_log_message(caplog, mode, exit_code):
        """Check that bookkeeper emits one meaningful last log message."""
        records = ((r, r.getMessage()) for r in reversed(caplog.records))
        skip_msg_contents = (EDITED_SUFFIX, 'domain folders')
        last_record, msg = next(
            (r, m)
            for r, m in records
            if m
            and not any(c in m for c in skip_msg_contents)
            )
        expect_level = {
            BookkeeperExitCode.FAIL: logging.ERROR,
            BookkeeperExitCode.NOTHING_TO_DO: logging.INFO,
            BookkeeperExitCode.SUCCESS: logging.INFO,
            }
        expect_msg = {
            BookkeeperExitCode.FAIL: 'Fail.*',
            BookkeeperExitCode.NOTHING_TO_DO: (
                '.*Exiting without doing anything.' if mode is Mode.ARCHIVE
                else 'Found nothing to do. Exiting...'
                ),
            BookkeeperExitCode.SUCCESS: '(Successfully|Done).*',
            }
        assert last_record.levelno == expect_level[exit_code]
        assert re.fullmatch(expect_msg[exit_code], msg)

    @staticmethod
    def check_metadata_exists(history_folder):
        """Test that the metadata file is present in `history_folder`."""
        assert (history_folder / _METADATA_NAME).is_file()

    @staticmethod
    def check_no_duplicate_logs(caplog):
        """Ensure all log messages appear only once."""
        all_messages = (r.getMessage() for r in caplog.records)
        messages = [m for m in all_messages if m.strip()]
        assert len(set(messages)) == len(messages)

    @staticmethod
    def check_no_edited_files(bookkeeper, *_):
        """Ensure no _edited files exist in bookkeeper's cwd."""
        assert not any((bookkeeper.cwd / f'{f}{EDITED_SUFFIX}').exists()
                       for f in MOCK_STATE_FILES)

    @staticmethod
    def check_no_warnings(caplog, exclude_msgs=()):
        """Check that there are no warnings or errors."""
        def _record_faulty(record):
            """Return whether record is not an excluded one."""
            if record.levelno < logging.WARNING:
                return False
            msg = str(record)
            return not any(exclude in msg for exclude in exclude_msgs)
        faulty = [rec for rec in caplog.records if _record_faulty(rec)]
        assert not any(faulty), f'Found: {faulty[0].getMessage()!r}'

    def check_out_files_at_path(self,
                                path,
                                out_suffixed=False,
                                missing=(),
                                **_):
        """Check that all output files are found at `path`/DEFAULT_OUT."""
        expected_contents = (MOCK_OUT_SUFFIXED_CONTENT if out_suffixed
                             else MOCK_OUT_CONTENT)
        out_path = path / DEFAULT_OUT
        for file in MOCK_STATE_FILES:
            if file in missing:
                assert not (out_path/file).exists()
                continue
            try:
                out_suffixed_file = next(out_path.glob(f'{file}_OUT_*'))
            except StopIteration:
                # pylint: disable-next=magic-value-comparison
                if out_suffixed and file == 'PARAMETERS':
                    # There was no PARAMETERS in OUT for < v0.13.0
                    continue
                if out_suffixed:
                    raise
                out_suffixed_file = out_path/f'{file}_OUT'  # A dummy
            expect_file = out_suffixed_file if out_suffixed else out_path/file
            other_file = out_path/file if out_suffixed else out_suffixed_file
            self._check_file_contents(expect_file, expected_contents)
            # Only one between _OUT-suffixed and non-suffixed
            assert not other_file.is_file()

    def check_out_files_untouched(self, bookkeeper, *_, **kwargs):
        """Ensure all expected files are found in OUT."""
        self.check_out_files_at_path(bookkeeper.cwd, **kwargs)

    def check_root_after_archive(self,
                                 *after_archive,
                                 out_suffixed=False,
                                 missing_out_files=(),
                                 **_):
        """Make sure the root is structured as expected after archiving."""
        self.check_root_inputs_renamed_to_ori(*after_archive)
        self.check_out_files_untouched(*after_archive,
                                       out_suffixed=out_suffixed,
                                       missing=missing_out_files)
        self.check_root_inputs_replaced_by_out_or_ori(
            *after_archive,
            out_suffixed=out_suffixed,
            )

    def check_root_inputs_renamed_to_ori(self, bookkeeper, *_):
        """Check that the input files have now a _ori suffix."""
        cwd = bookkeeper.cwd
        for file in MOCK_STATE_FILES:
            self._check_file_contents(cwd/f'{file}{ORI_SUFFIX}',
                                      MOCK_INPUT_CONTENT)

    def check_root_inputs_replaced_by_out_or_ori(self, bookkeeper, *_,
                                                 out_suffixed=False):
        """Check that input files in root come from OUT or original_inputs."""
        cwd = bookkeeper.cwd
        out_contents = (MOCK_OUT_SUFFIXED_CONTENT if out_suffixed
                        else MOCK_OUT_CONTENT)
        for file in MOCK_STATE_FILES:
            self._check_file_contents(cwd/file,
                                      out_contents,
                                      MOCK_ORIG_CONTENT)
            if any((cwd/DEFAULT_OUT).glob(f'{file}*')):
                self._check_file_contents(cwd/file, out_contents)

    def check_root_inputs_untouched(self, bookkeeper, *_):
        """Check the the original state files have not been moved."""
        cwd = bookkeeper.cwd
        for file in MOCK_STATE_FILES:
            self._check_file_contents(cwd/file, MOCK_INPUT_CONTENT)

    @staticmethod
    def check_root_is_clean(bookkeeper, *_):
        """Check that no calc output is present in the main directory."""
        cwd = bookkeeper.cwd
        for file in MOCK_STATE_FILES:
            assert not (cwd / f'{file}{ORI_SUFFIX}').is_file()
        assert not (cwd / DEFAULT_SUPP).exists()
        assert not (cwd / DEFAULT_OUT).exists()
        assert not any(cwd.glob('*.log'))

    def check_root_reverted_to_previous_calc_run(self, bookkeeper, *_):
        """Ensure root contains the same as after the previous calc run."""
        self.check_root_is_clean(bookkeeper)
        self.check_root_inputs_untouched(bookkeeper)

    def has_out_suffixed(self, bookkeeper, *_):
        """Return whether there are any old-style _OUT files in OUT.

        It is critical to execute this function BEFORE bokkeeper runs.

        Parameters
        ----------
        bookkeeper : Bookkeeper
            A bookkeeper instance. Its .cwd/OUT folder is checked
            for _OUT files.
        *_ : object
            Other unused positional arguments.

        Returns
        -------
        bool
        """
        return any((bookkeeper.cwd/DEFAULT_OUT).glob('*_OUT*'))


class _TestBookkeeperRunBase(_TestCheckers):
    """Base class for checking correct execution of bookkeeper."""

    mode = None

    def check_domains_archived(self, main_run, domains, skip=()):
        """Check that all subdomains have been archived correctly."""
        *_, mocker = main_run
        # We call the same methods as run_archive_after_calc_and_check
        # does, but on all subdomain root paths.
        for domain_root, domain_info in domains.items():
            mock_run = (mocker.MagicMock(cwd=domain_root),
                        domain_info['history_run'],
                        mocker)
            self.check_has_archived(mock_run, mode=self.mode)

        # Make sure the domains have been correctly registered,
        # and that the correct information has been written to
        # history.info files.
        self._check_domains_registered(main_run, domains, skip)
        self._check_domains_history_info(main_run, domains, skip)

    def collect_root_contents(self, bookkeeper):
        """Return a dictionary of the current contents of bookkeeper's CWD."""
        skip = {  # Files whose contents are complex to keep track of
            BOOKIE_LOGFILE,
            HISTORY_INFO_NAME,
            _METADATA_NAME,
            }
        return filesystem_to_dict(bookkeeper.cwd, skip=skip)

    def patch_for_domains(self, mocker):
        """Exclude checks that are knowingly going to fail in DOMAINS."""
        # Patch away checkers that we know will fail:
        # _check_workhistory_archived: main has no workhistory folder,
        #   while, for the subdomains, the target names of folders in
        #   history differ from those in MOCK_WORKHISTORY.
        mocker.patch.object(self, '_check_workhistory_archived')
        # check_no_duplicate_logs: messages repeated for each domain
        mocker.patch.object(self, 'check_no_duplicate_logs')

    def run_after_archive_and_check(self, after_archive, caplog, **kwargs):
        """Check that running bookkeeper after ARCHIVE does basic stuff."""
        bookkeeper, *_ = after_archive
        # See if we have legacy _OUT-suffixed files:
        has_out_suffixed = self.has_out_suffixed(bookkeeper)
        # bookkeeper should not think that it needs archiving
        assert not bookkeeper.archiving_required
        self._run_bookkeeper(bookkeeper, kwargs, caplog)
        bookkeeper.update_from_cwd(silent=True)
        kwargs['has_out_suffixed'] = has_out_suffixed
        self.check_has_archived(after_archive, **kwargs)

    def run_after_calc_exec_and_check(self,
                                      after_calc_execution,
                                      caplog,
                                      check_archiving_required=True,
                                      missing_out_files=(),
                                      **kwargs):
        """Check that running bookkeeper after calc does some basic stuff."""
        bookkeeper, *_, mocker = after_calc_execution
        # See if we have legacy _OUT-suffixed files:
        has_out_suffixed = self.has_out_suffixed(bookkeeper)
        if check_archiving_required:
            # bookkeeper should think that it needs archiving
            assert bookkeeper.archiving_required
        # Wrap the check for edited files to ensure it is executed
        check_edited_files = mocker.patch.object(
            # pylint: disable-next=protected-access       # OK in tests
            bookkeeper._root,
            'complain_about_edited_files',
            # pylint: disable-next=protected-access       # OK in tests
            wraps=bookkeeper._root.complain_about_edited_files
            )
        self._run_bookkeeper(bookkeeper, kwargs, caplog)
        bookkeeper.update_from_cwd(silent=True)
        kwargs['has_out_suffixed'] = has_out_suffixed
        kwargs['missing_out_files'] = missing_out_files
        self.check_has_archived(after_calc_execution, **kwargs)
        check_edited_files.assert_called()

    @staticmethod
    def _edit_input_file_contents(at_path):
        """Modify the contents of inputs `at_path` with some sentinel text."""
        sentinel_text = ('This is some other text to ensure file contents '
                         'are not modified when running twice')
        for file in MOCK_STATE_FILES:
            (at_path / file).write_text(sentinel_text)

    def run_again_and_check_nothing_changed(self, run, caplog,
                                            acceptable_warnings=()):
        """Ensure running bookkeeper again with the same mode does nothing."""
        bookkeeper, *_ = run
        # We will collect the current state of the root folder (except
        # for the bookkeeper.log that changes) and ensure that running
        # doesn't change it.
        kwargs, skip = {}, {BOOKIE_LOGFILE}

        # To ensure the input files are not fiddled with, replace
        # their contents with some other text.
        cwd = bookkeeper.cwd
        self._edit_input_file_contents(cwd)

        # Do the same for any domain
        finder = DomainFinder(bookkeeper)
        # About the disable: the public find_potential_domains method
        # raises NotImplementedError since the fix for #344. However,
        # here it is still OK to use the version that is problematic
        # with respect to the work folder, as we do not have a work
        # folder.
        # pylint: disable-next=protected-access
        domains = finder._find_potential_domains_buggy()
        for domain in domains:
            self._edit_input_file_contents(domain)
        before_run = filesystem_to_dict(cwd, skip=skip)
        exit_code = self._run_bookkeeper(bookkeeper, kwargs, caplog)
        after_run = filesystem_to_dict(cwd, skip=skip)

        # Exit code should be SUCCESS or NOTHING_TO_DO
        self.check_exit_code_ok(exit_code)
        assert after_run == before_run
        self.check_no_warnings(caplog, exclude_msgs=acceptable_warnings)

    def run_and_check_prerun_archiving(self, after_calc_execution, caplog,
                                       exclude_warnings=()):
        """Execute bookkeeper, and verify it has archived a calc run."""
        self.run_after_calc_exec_and_check(after_calc_execution, caplog)
        self.check_no_warnings(caplog,
                               exclude_msgs=('metadata',)+exclude_warnings)
        self.check_root_is_clean(*after_calc_execution)

        # Original SHOULD NOT be replaced by output:
        # ARCHIVE does not run only if the run crashed,
        # in which case we don't want to overwrite
        self.check_root_inputs_untouched(*after_calc_execution)

        # There should be no file marked as _edited
        self.check_no_edited_files(*after_calc_execution)

    def run_and_check_prerun_archiving_domains(self,
                                               domains_after_calc_execution,
                                               caplog,
                                               exclude_warnings=(),
                                               already_processed=()):
        """Execute bookkeeper, and verify it has archived a domain run."""
        *main_run, domains = domains_after_calc_execution
        *_, mocker = main_run
        self.patch_for_domains(mocker)
        self.run_and_check_prerun_archiving(main_run,
                                            caplog,
                                            exclude_warnings=exclude_warnings)

        # Repeat the same checks as in run_and_check_prerun_archiving
        # for all subdomains.
        self.check_domains_archived(main_run, domains, skip=already_processed)
        for domain_path in domains:
            if domain_path in already_processed:
                continue
            mock_bookie = mocker.MagicMock(cwd=domain_path)
            self.check_root_is_clean(mock_bookie)
            self.check_root_inputs_untouched(mock_bookie)
            self.check_no_edited_files(mock_bookie)

    def run_archive_after_calc_and_check(self, after_calc_execution,
                                         caplog,
                                         check_archiving_required=True,
                                         **kwargs):
        """Check correct storage of history files in ARCHIVE mode."""
        # See if we have legacy _OUT-suffixed files:
        has_out_suffixed = self.has_out_suffixed(*after_calc_execution)
        kwargs.update({
            'check_archiving_required': check_archiving_required,
            'mode': 'archive',
            })
        self.run_after_calc_exec_and_check(after_calc_execution,
                                           caplog,
                                           **kwargs)
        kwargs['out_suffixed'] = has_out_suffixed
        self.check_root_after_archive(*after_calc_execution, **kwargs)
        self.check_no_warnings(caplog, exclude_msgs=('metadata',))

    def run_before_calc_exec_and_check(self,
                                       before_calc_execution,
                                       caplog,
                                       check_archiving_required=True,
                                       **kwargs):
        """Check that running bookkeeper before calc does almost nothing."""
        bookkeeper = before_calc_execution
        if check_archiving_required:
            # bookkeeper should not think that it needs archiving
            assert not bookkeeper.archiving_required
        # NOTICE: we do not update_from_cwd on purpose. This way,
        # the bookkeeper is left in the exact same state it has
        # after an execution as part of the viperleed.calc CLI.
        exit_code = self._run_bookkeeper(bookkeeper, kwargs, caplog)
        # Bookkeeper should not do anything (except for logging)
        self.check_history_folder_empty(before_calc_execution)
        self.check_root_inputs_untouched(before_calc_execution)
        self.check_exit_code_ok(exit_code)

    def _run_bookkeeper(self, bookkeeper, kwargs, caplog):
        """Execute a bookkeeper run and return its exit code."""
        caplog.set_level(5)  # < DEBUG
        kwargs.setdefault('mode', self.mode)
        exit_code = bookkeeper.run(**kwargs)
        kwargs['mode'] = Mode(kwargs['mode'])
        # NOTICE: we do not update_from_cwd on purpose. This way,
        # the bookkeeper is left in the exact same state it has
        # after an execution as part of the viperleed.calc CLI.
        self.check_no_duplicate_logs(caplog)
        self.check_last_log_message(caplog, kwargs['mode'], exit_code)
        return exit_code
