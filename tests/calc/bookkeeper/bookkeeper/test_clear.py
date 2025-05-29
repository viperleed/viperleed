"""Tests for module bookkeeper of viperleed.calc.bookkeeper.

Collects tests for running bookkeeper in CLEAR mode.
"""

__authors__ = (
    'Alexander M. Imre (@amimre)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2023-08-02'
__license__ = 'GPLv3+'

from copy import deepcopy

from pytest_cases import fixture
from pytest_cases import parametrize

from viperleed.calc.bookkeeper.constants import ORI_SUFFIX
from viperleed.calc.bookkeeper.constants import STATE_FILES
from viperleed.calc.bookkeeper.mode import BookkeeperMode as Mode
from viperleed.calc.constants import DEFAULT_OUT
from viperleed.calc.constants import DEFAULT_SUPP

from .run_bookkeeper_base import _TestBookkeeperRunBase


@fixture(name='after_clear')
def fixture_after_clear(after_archive, after_bookkeper_run):
    """Prepare a directory like the one after CLEAR was executed."""
    return after_bookkeper_run(after_archive, Mode.CLEAR)


class TestBookkeeperClear(_TestBookkeeperRunBase):
    """Tests for correct behavior of CLEAR bookkeeper runs."""

    mode = Mode.CLEAR

    def test_run_before_calc_exec(self, before_calc_execution, caplog):
        """Check correct overwriting of input files in CLEAR mode."""
        self.run_before_calc_exec_and_check(before_calc_execution, caplog)
        self.check_no_warnings(caplog)

    def test_clear_after_archive(self, after_archive, caplog):
        """Check behavior of CLEAR after ARCHIVE (e.g., manual call)."""
        has_out_suffixed = self.has_out_suffixed(*after_archive)
        self.run_after_archive_and_check(after_archive, caplog)
        self.check_root_inputs_replaced_by_out_or_ori(
            *after_archive,
            out_suffixed=has_out_suffixed,
            )
        self.check_no_warnings(caplog, exclude_msgs=('metadata',))

    def test_clear_after_archive_domains(self, archived_domains, caplog):
        """Check behavior of CLEAR after ARCHIVE (e.g., manual call)."""
        *main_run, domains = archived_domains
        *_, mocker = main_run
        self.patch_for_domains(mocker)
        self.run_after_archive_and_check(main_run, caplog)
        self.check_root_inputs_replaced_by_out_or_ori(*main_run)
        # Repeat the check for all subdomains
        for domain_path in domains:
            mock_bookie = mocker.MagicMock(cwd=domain_path)
            self.check_root_inputs_replaced_by_out_or_ori(mock_bookie)
        self.check_no_warnings(caplog, exclude_msgs=('metadata',))

    def test_clear_after_calc_exec(self, after_calc_execution, caplog):
        """Check behavior of CLEAR after a non-ARCHIVEd calc run.

        This may happen, for example, if the previous (calc or
        bookkeeper) execution crashed.

        Parameters
        ----------
        after_calc_execution : fixture
            A bookkeeper and information on a root directory right
            after viperleed.calc has run, and before any bookkeeper
            execution (even the default --archive).
        caplog : fixture
            The pytest.caplog fixture.

        Returns
        -------
        None.
        """
        self.run_and_check_prerun_archiving(after_calc_execution, caplog)

    def get_tree_from_archived(self, archived):
        """Return a dictionary of the root after CLEAR from the ARCHIVE one."""
        cleared = deepcopy(archived)
        _deleted = (
            DEFAULT_OUT,
            DEFAULT_SUPP,
            *(f'{file}{ORI_SUFFIX}' for file in STATE_FILES),
            *(file for file in cleared if file.endswith('.log')),
            )
        for file in _deleted:
            try:
                del cleared[file]
            except KeyError:
                pass
        return cleared

    def test_clear_after_calc_with_edited(self,
                                          after_calc_with_edited_file,
                                          caplog):
        """Check root after CLEARing a non-ARCHIVEd run with user edits."""
        bookkeeper, before, archived = after_calc_with_edited_file
        self._run_bookkeeper(bookkeeper, {}, caplog)

        expect = self.get_tree_from_archived(archived)
        # While the above is identical to the ARCHIVE-then-CLEAR
        # case, we DO NOT pull an _edited file from original_inputs,
        # as we want to avoid running with the "wrong" file.
        # See viperleed/pull/198#issuecomment-2508005204.
        del expect['POSCAR']

        # Also, differently from the ARCHIVE-then-CLEAR case, we DO NOT
        # pull OUT files as new inputs. We expect users to have done so
        # when they realized that ARCHIVE failed.
        expect['PARAMETERS'] = before['PARAMETERS']
        assert self.collect_root_contents(bookkeeper) == expect

    def test_clear_archived_with_edited(self,
                                        after_calc_with_edited_file,
                                        after_bookkeper_run,
                                        caplog):
        """Check root after CLEARing an ARCHIVEd run with user edits."""
        # Run in archive mode first to mimic the root structure
        bookkeeper, _, archived = after_bookkeper_run(
            after_calc_with_edited_file,
            'archive',
            )
        caplog.clear()  # The ARCHIVE logs
        self._run_bookkeeper(bookkeeper, {}, caplog)

        # See viperleed/pull/198#issuecomment-2506549827
        expect = self.get_tree_from_archived(archived)
        assert self.collect_root_contents(bookkeeper) == expect

    def test_clear_domains(self, domains_after_calc_execution, caplog):
        """Check behavior of CLEAR after a non-ARCHIVEd domain run."""
        self.run_and_check_prerun_archiving_domains(
            domains_after_calc_execution,
            caplog,
            )

    @parametrize(mode=(Mode.ARCHIVE, Mode.CLEAR, Mode.DISCARD))
    def test_one_domain_already_processed(self,
                                          mode,
                                          manual_run_one_domain,
                                          caplog):
        """Check CLEAR when one domain subfolder was manually processed."""
        domain_run, already_processed = manual_run_one_domain(mode)
        pre_archived = mode is Mode.ARCHIVE
        *_, mocker, _ = domain_run
        caplog.clear()
        self.run_and_check_prerun_archiving_domains(
            domain_run,
            caplog,
            already_processed=(already_processed,),
            )
        # Check also the one that was manually processed.
        mock_bookie = mocker.MagicMock(cwd=already_processed)
        # The root folder should be free of calc results, irrespective
        # of whether the domain was processed or not (since we run in
        # clear mode in the root).
        self.check_root_is_clean(mock_bookie)
        check_inputs = (
            self.check_root_inputs_replaced_by_out_or_ori if pre_archived
            else self.check_root_inputs_untouched
            )
        check_inputs(mock_bookie)
        self.check_no_edited_files(mock_bookie)

    def test_twice(self, after_clear, caplog):
        """Ensure running CLEAR twice does nothing."""
        warnings = ('metadata',)
        self.run_again_and_check_nothing_changed(after_clear, caplog,
                                                 acceptable_warnings=warnings)

    def test_twice_domains(self,
                           archived_domains,
                           after_bookkeper_run,
                           caplog,
                           mocker):
        """Ensure running CLEAR twice does nothing."""
        after_clear = after_bookkeper_run(archived_domains, self.mode)
        caplog.clear()
        self.patch_for_domains(mocker)
        self.test_twice(after_clear, caplog)
