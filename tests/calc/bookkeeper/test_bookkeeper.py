"""Tests for module viperleed.calc.bookkeeper."""

__authors__ = (
    'Alexander M. Imre (@amimre)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2023-08-02'
__license__ = 'GPLv3+'

import logging

from pytest_cases import fixture

from viperleed.calc import DEFAULT_HISTORY
from viperleed.calc import ORIGINAL_INPUTS_DIR_NAME
from viperleed.calc.bookkeeper.bookkeeper import Bookkeeper
from viperleed.calc.bookkeeper.constants import HISTORY_INFO_NAME
from viperleed.calc.bookkeeper.history import _DISCARDED
from viperleed.calc.bookkeeper.mode import BookkeeperMode
from viperleed.calc.sections.cleanup import DEFAULT_OUT
from viperleed.calc.sections.cleanup import DEFAULT_SUPP

from ...helpers import execute_in_dir
from .conftest import MOCK_INPUT_CONTENT
from .conftest import MOCK_ORIG_CONTENT
from .conftest import MOCK_OUT_CONTENT
from .conftest import MOCK_TIMESTAMP
from .conftest import MOCK_STATE_FILES
from .conftest import NOTES_TEST_CONTENT


@fixture(name='after_archive')
def fixture_after_archive(after_calc_run):
    """Yield a temporary directory for testing the bookkeeper."""
    bookkeeper, *_ = after_calc_run
    bookkeeper.run(mode=BookkeeperMode.ARCHIVE)
    bookkeeper.update_from_cwd(silent=True)
    return after_calc_run


@fixture(name='before_calc_run')
def fixture_before_calc_run(tmp_path):
    """Yield a temporary directory for testing the bookkeeper.

    This represents a new calculation, i.e., before any viperleed calc or
    bookkeeper run.

    Parameters
    ----------
    tmp_path : fixture
        Path to a temporary directory.

    Yields
    ------
    bookkeeper : Bookkeeper
        A bookkeeper instance ready for running in tmp_path.
    """
    # create mock input files
    for file in MOCK_STATE_FILES:
        (tmp_path / file).write_text(MOCK_INPUT_CONTENT)

    bookkeeper = Bookkeeper(cwd=tmp_path)
    bookkeeper.update_from_cwd(silent=True)
    with execute_in_dir(tmp_path):
        yield bookkeeper
    # It would be nice to clean up, but the following line causes
    # a PermissionError. Likely because of logging keeping a hold
    # of the bookkeeper.log file
    # shutil.rmtree(tmp_path)


@fixture(name='history_path')
def fixture_history_path(mock_tree_after_calc_run):
    """Return the path to a history subfolder of `mock_tree_after_calc_run`."""
    return mock_tree_after_calc_run / DEFAULT_HISTORY


@fixture(name='history_run_path')
def fixture_history_run_path(history_path):
    """Return the path to a history run subfolder of `history_path`."""
    return history_path / f't000.r001_{MOCK_TIMESTAMP}'


class _TestBookkeeperRunBase:
    """Base class for checking correct execution of bookkeeper."""

    mode = None

    def check_history_exists(self, bookkeeper, history_run_path):
        """Test that history_path and directory/history.info exist."""
        assert history_run_path.is_dir()
        assert (bookkeeper.cwd / HISTORY_INFO_NAME).exists()

    def check_history_folder_empty(self, bookkeeper, *_):
        """Test that an empty history folder exists, except bookkeeper.log."""
        cwd = bookkeeper.cwd
        assert (cwd / DEFAULT_HISTORY).exists()
        history_contents = (f for f in (cwd / DEFAULT_HISTORY).iterdir()
                            if f.name != 'bookkeeper.log')
        assert not any(history_contents)

    def check_no_warnings(self, caplog, contain=None, exclude_msgs=()):
        """Check that there are no warnings or errors."""
        if contain is None:
            def _record_faulty(record):
                return record.levelno >= logging.WARNING
        else:
            def _record_faulty(record):
                return (record.levelno >= logging.WARNING
                        and contain in str(record))
        assert not any(
            _record_faulty(rec)
            for rec in caplog.records
            if not any(exclude in str(rec) for exclude in exclude_msgs)
            )

    def check_out_files_in_history(self, *run):
        """Check that the expected state files are stored in 'OUT'."""
        *_, history_run_path = run
        for file in MOCK_STATE_FILES:
            hist_file = history_run_path / DEFAULT_OUT / f'{file}_OUT'
            assert hist_file.is_file()
            hist_content = hist_file.read_text()
            assert MOCK_OUT_CONTENT in hist_content

    def check_root_inputs_replaced_by_out(self, bookkeeper, *_):
        """Check that the input files in root come from OUT."""
        cwd = bookkeeper.cwd
        for file in MOCK_STATE_FILES:
            out_content = (cwd / file).read_text()
            assert MOCK_OUT_CONTENT in out_content

    def check_root_inputs_untouched(self, bookkeeper, *_):
        """Check the the original state files have not been moved."""
        cwd = bookkeeper.cwd
        for file in MOCK_STATE_FILES:
            assert (cwd / file).is_file()
            input_content = (cwd / file).read_text()
            assert MOCK_INPUT_CONTENT in input_content

    def check_root_is_clean(self, bookkeeper, *_):
        """Check that no calc output is present in the main directory."""
        cwd = bookkeeper.cwd
        for file in MOCK_STATE_FILES:
            assert not (cwd / f'{file}_ori').is_file()
        assert not (cwd / DEFAULT_SUPP).exists()
        assert not (cwd / DEFAULT_OUT).exists()
        assert not any(cwd.glob('*.log'))

    def _test_after_archive_base(self, after_archive):
        """Check that running bookkeeper after ARCHIVE does basic stuff."""
        bookkeeper, *_ = after_archive
        # bookkeeper should not think that it needs archiving
        assert not bookkeeper.archiving_required
        bookkeeper.run(mode=self.mode)
        bookkeeper.update_from_cwd(silent=True)
        self.check_history_exists(*after_archive)
        self.check_out_files_in_history(*after_archive)
        self.check_root_is_clean(*after_archive)

    def _test_after_calc_run_base(self, after_calc_run):
        """Check that running bookkeeper after calc does some basic stuff."""
        bookkeeper, *_ = after_calc_run
        # bookkeeper should think that it needs archiving
        assert bookkeeper.archiving_required
        bookkeeper.run(mode=self.mode)
        bookkeeper.update_from_cwd(silent=True)
        self.check_history_exists(*after_calc_run)
        self.check_out_files_in_history(*after_calc_run)

    def test_before_calc_run(self, before_calc_run):
        """Check that running bookkeeper before calc does almost nothing."""
        bookkeeper = before_calc_run
        # bookkeeper should not think that it needs archiving
        assert not bookkeeper.archiving_required
        bookkeeper.run(mode=self.mode)
        # Bookkeeper should not do anything (except for logging)
        self.check_history_folder_empty(before_calc_run)
        self.check_root_inputs_untouched(before_calc_run)


class TestBookkeeperArchive(_TestBookkeeperRunBase):
    """Tests for correct behavior of ARCHIVE bookkeeper runs."""

    mode = BookkeeperMode.ARCHIVE

    def _check_root_inputs_renamed_to_ori(self, bookkeeper, *_):
        """Check that the input files have now a _ori suffix."""
        cwd = bookkeeper.cwd
        for file in MOCK_STATE_FILES:
            ori_file = cwd / f'{file}_ori'
            assert ori_file.is_file()
            input_content = (cwd / file).read_text()
            assert MOCK_OUT_CONTENT in input_content

    def test_archive_after_calc_run(self, after_calc_run, caplog):
        """Check correct storage of history files in ARCHIVE mode."""
        self._test_after_calc_run_base(after_calc_run)
        self._check_root_inputs_renamed_to_ori(*after_calc_run)
        self.check_root_inputs_replaced_by_out(*after_calc_run)
        self.check_no_warnings(caplog)

    def test_archive_again(self, after_archive, caplog):
        """Bookkeeper ARCHIVE after ARCHIVE should not do anything."""
        bookkeeper, *_ = after_archive
        assert not bookkeeper.archiving_required
        # write stuff to files to check they are not overwritten
        cwd = bookkeeper.cwd
        for file in MOCK_STATE_FILES:
            (cwd / file).write_text('something else')
        bookkeeper.run(mode=self.mode)
        for file in MOCK_STATE_FILES:
            assert (cwd / file).read_text() == 'something else'
        self.check_no_warnings(caplog)

    # pylint: disable-next=arguments-differ
    def test_before_calc_run(self, before_calc_run, caplog):
        """Check no archiving happens before calc runs."""
        super().test_before_calc_run(before_calc_run)
        self.check_no_warnings(caplog)


class TestBookkeeperClear(_TestBookkeeperRunBase):
    """Tests for correct behavior of CLEAR bookkeeper runs."""

    mode = BookkeeperMode.CLEAR

    # pylint: disable-next=arguments-differ
    def test_before_calc_run(self, before_calc_run, caplog):
        """Check correct overwriting of input files in CLEAR mode."""
        super().test_before_calc_run(before_calc_run)
        self.check_no_warnings(caplog)

    def test_clear_after_archive(self, after_archive, caplog):
        """Check behavior of CLEAR after ARCHIVE (e.g., manual call)."""
        self._test_after_archive_base(after_archive)
        self.check_root_inputs_replaced_by_out(*after_archive)
        self.check_no_warnings(caplog)

    def test_clear_after_calc_run(self, after_calc_run, caplog):
        """Check behavior of CLEAR after a non-ARCHIVEd calc run.

        This may happen, for example, if the previous (calc or
        bookkeeper) execution crashed.

        Parameters
        ----------
        after_calc_run: fixture
        caplog: fixture

        Returns
        -------
        None.
        """
        self._test_after_calc_run_base(after_calc_run)
        self.check_no_warnings(caplog)
        self.check_root_is_clean(*after_calc_run)

        # Original SHOULD NOT be replaced by output:
        # ARCHIVE does not run only if the run crashed,
        # in which case we don't want to overwrite
        self.check_root_inputs_untouched(*after_calc_run)


class TestBookkeeperDiscard(_TestBookkeeperRunBase):
    """Tests for correct behavior of DISCARD bookkeeper runs."""

    mode = BookkeeperMode.DISCARD

    # pylint: disable-next=arguments-differ
    def test_before_calc_run(self, before_calc_run, caplog):
        """Check correct overwriting of input files in CLEAR mode."""
        super().test_before_calc_run(before_calc_run)
        self.check_no_warnings(
            caplog,
            exclude_msgs=('Failed to mark last entry as discarded',),
            )

    def test_discard_after_archive(self, after_archive, caplog):
        """Check reverting of state when DISCARDing an ARCHIVEd calc run."""
        self._test_after_archive_base(after_archive)

        # Original be replaced by output                                        # TODO: this does something else!
        bookkeeper, *_ = after_archive
        cwd = bookkeeper.cwd
        for file in MOCK_STATE_FILES:
            out_content = (cwd / file).read_text()
            assert MOCK_INPUT_CONTENT in out_content
        # A 'DISCARDED' note should be in history.info
        assert _DISCARDED in bookkeeper.history_info.path.read_text()
        # Some fields are knowingly faulty,
        # but we can still DISCARD them.
        faulty_entry_logs = (
            'Found entry with',
            'Could not understand',
            'Faulty entry is',
            )
        self.check_no_warnings(caplog, exclude_msgs=faulty_entry_logs)

    def test_discard_after_calc_run(self, after_calc_run, caplog):
        """Check behavior of DISCARD after a non-ARCHIVEd calc run.

        This may happen, for example, if the previous (calc or
        bookkeeper) execution crashed.

        Parameters
        ----------
        after_calc_run: fixture
        caplog: fixture

        Returns
        -------
        None.
        """
        self._test_after_calc_run_base(after_calc_run)
        self.check_no_warnings(caplog, exclude_msgs=('discarded',))
        self.check_root_is_clean(*after_calc_run)

        # Original SHOULD NOT be replaced by output:
        # ARCHIVE does not run only if the run crashed,
        # in which case we don't want to overwrite
        self.check_root_inputs_untouched(*after_calc_run)

        # A 'DISCARDED' note should be in history.info
        bookkeeper, *_ = after_calc_run
        assert bookkeeper.history_info.last_entry_was_discarded
        assert _DISCARDED in bookkeeper.history_info.path.read_text()


class TestBookkeeperDiscardFull(_TestBookkeeperRunBase):
    """Tests for correct behavior of DISCARD_FULL bookkeeper runs."""

    mode = BookkeeperMode.DISCARD_FULL

    # pylint: disable-next=arguments-differ
    def test_before_calc_run(self, before_calc_run, caplog):
        """Check correct overwriting of input files in DISCARD_FULL mode.

        This should do the same as a normal DISCARD.

        Parameters
        ----------
        before_calc_run: fixture
        caplog: fixture

        Returns
        -------
        None.
        """
        super().test_before_calc_run(before_calc_run)
        self.check_no_warnings(
            caplog,
            exclude_msgs=('remove',)   # Catches two expected warnings
            )

    def test_discard_full_after_archive(self, after_archive, caplog):
        """Check correct behavior of DISCARD_FULL on an ARCHIVEd calc run.

        Should revert state to before the run and purge the history.

        Parameters
        ----------
        after_archive: fixture
        caplog: fixture

        Returns
        -------
        None.
        """
        bookkeeper, *_, history_run_path = after_archive
        last_entry = bookkeeper.history_info.last_entry
        bookkeeper.run(mode=self.mode)
        if not last_entry.can_be_removed:
            # Should prevent the removal of the history
            assert history_run_path.is_dir()
            return
        assert not history_run_path.is_dir()
        self.check_root_is_clean(*after_archive)
        self.check_root_inputs_untouched(*after_archive)
        self.check_no_warnings(caplog)

    def test_discard_full_after_calc_run(self, after_calc_run, caplog):
        """Check behavior of DISCARD_FULL on a non-ARCHIVEd calc run."""
        bookkeeper, *_, history_run_path = after_calc_run
        cwd = bookkeeper.cwd
        notes_in_history_info = (
            NOTES_TEST_CONTENT in (cwd / HISTORY_INFO_NAME).read_text()
            )
        last_entry = bookkeeper.history_info.last_entry
        bookkeeper.run(mode=self.mode)
        # Since the run was not archived, the history should be empty
        assert not history_run_path.is_dir()
        if notes_in_history_info:
            assert 'last entry in history.info has user notes' in caplog.text
        elif last_entry and not last_entry.can_be_removed:
            expected_logs = (
                'contains invalid fields that could not be interpreted',
                'contains fields with non-standard format',
                'some expected fields were deleted',
                )
            assert any(msg in caplog.text for msg in expected_logs)
        else:
            expected_logs = (
                'could not identify directory to remove',
                'No entries to remove',
                )
            assert any(msg in caplog.text for msg in expected_logs)
