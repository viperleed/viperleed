"""Tests for module deltas of viperleed.calc.sections.

This module defines tests for the deltas_domains function, the main
entry point of a delta-amplitudes calculation for a multi-domain
system.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2025-07-05'
__license__ = 'GPLv3+'

import logging
from pathlib import Path

import pytest
from pytest_cases import fixture

from viperleed.calc.sections.deltas import DeltaCompileTask
from viperleed.calc.sections.deltas import DeltaRunTask
from viperleed.calc.sections.deltas import deltas_domains
from viperleed.calc.sections.deltas import run_delta

from .conftest import _MODULE

use = pytest.mark.usefixtures


@fixture(name='make_domain')
def factory_make_domain(mocker, tmp_path):
    """Return a fake domain."""
    def _make(name):
        work = tmp_path/name
        mock_attrs = {
            'slab': mocker.MagicMock(name='slab'),
            'rpars': mocker.MagicMock(name='rpars'),
            'workdir': work,
            }
        work.mkdir()
        return mocker.MagicMock(**mock_attrs)
    return _make


@fixture(name='make_rpars')
def factory_rpars(make_domain, rpars):
    """Return an Rparams with domains."""
    def _make(n_domains):
        rpars.domainParams = [make_domain(name=f'domain {i}')
                              for i in range(n_domains)]
        return rpars
    return _make


@fixture(name='mocks')
def fixture_mock_implementation(mocker):
    """Replace implementation details with mocks."""
    return {
        'one_domain': mocker.patch(f'{_MODULE}.deltas', return_value=None),
        'pool': mocker.patch(f'{_MODULE}.parallelization.monitoredPool'),
        'copy_log': mocker.patch(f'{_MODULE}.leedbase.copy_compile_log'),
        'rmtree': mocker.patch('shutil.rmtree'),
        'checksums': mocker.patch(f'{_MODULE}.validate_multiple_files'),
        }


@fixture(name='mock_tasks')
def fixture_mock_domain_tasks(mocks, rpars):
    """Make the single-domain implementation return given tasks."""
    def _make(n_comptasks, n_runtasks):
        if len(n_comptasks) != len(n_runtasks):
            raise ValueError('Inconsistent number of domains from '
                             'comp_tasks and run_tasks.')
        src_dir = rpars.get_tenserleed_directory.return_value.path
        comp_tasks = (
            [DeltaCompileTask(f'param for domain {d}', 'hash', src_dir, i)
             for i in range(n_comptask_domain)]
            for d, n_comptask_domain in enumerate(n_comptasks)
            )
        run_tasks = (
            [DeltaRunTask(f'comptask {i} of domain {d}')
             for i in range(n_runtask_domain)]
            for d, n_runtask_domain in enumerate(n_runtasks)
            )
        all_tasks = tuple(zip(comp_tasks, run_tasks))
        mocks['one_domain'].side_effect = all_tasks
        return [t
                for domain_comp_tasks, _ in all_tasks
                for t in domain_comp_tasks]
    return _make


class TestDeltasDomains:
    """Tests for the deltas_domains function."""

    def test_compile_and_run(self, make_rpars, mocks, mock_tasks, caplog):
        """Check the expected calls for a successful full execution."""
        caplog.set_level(logging.INFO)
        n_domains = 2
        rpars = make_rpars(n_domains)
        n_comptasks = (2, 3)
        n_runtasks = (3, 3)
        comp_tasks = mock_tasks(n_comptasks, n_runtasks)
        deltas_domains(rpars)
        n_parallel = 2  # compile, run
        n_comptasks_total = sum(n_comptasks)

        assert mocks['one_domain'].call_count == n_domains
        assert mocks['pool'].call_count == n_parallel
        assert mocks['copy_log'].call_count == n_comptasks_total
        assert mocks['rmtree'].call_count == n_comptasks_total
        mocks['checksums'].assert_called_once()
        rpars.updateCores.assert_called_once()
        assert all(c.fortran_comp is rpars.FORTRAN_COMP for c in comp_tasks)

        expect_log = (
            'Getting input for delta calculations:',
            'Getting input for delta calculations:',
            'Compiling fortran files...',
            'Running delta calculations...',
            'Delta calculations finished.',
            )
        log_info = [r.getMessage()
                    for r in caplog.records
                    if r.levelno == logging.INFO]
        assert len(log_info) == len(expect_log)
        assert all(msg.startswith(expect)
                   for msg, expect in zip(log_info, expect_log))

    @use('mocks')
    def test_fortran_comp_known(self, make_rpars, mock_tasks):
        """Check that no compiler is searched if already available."""
        n_domains = 5
        rpars = make_rpars(n_domains)
        rpars.FORTRAN_COMP[0] = 'some compiler'
        mock_tasks(n_comptasks=(1,)*n_domains, n_runtasks=(1,)*n_domains)
        deltas_domains(rpars)
        rpars.getFortranComp.assert_not_called()

    @use('mocks')
    def test_fortran_comp_unknown(self, make_rpars, mock_tasks):
        """Check that no compiler is searched if already available."""
        n_domains = 5
        rpars = make_rpars(n_domains=5)
        rpars.FORTRAN_COMP[0] = ''
        mock_tasks(n_comptasks=(1,)*n_domains, n_runtasks=(1,)*n_domains)
        deltas_domains(rpars)
        rpars.getFortranComp.assert_called_once()

    def test_execute_in_domain_work(self, make_rpars, mocks):
        """Check that single-domain inputs are created in the correct paths."""
        rpars = make_rpars(4)
        expect_directories = [d.workdir for d in rpars.domainParams]
        directories = []
        def _register_directory(*_, **__):
            directories.append(Path.cwd())
        mocks['one_domain'].side_effect = _register_directory
        deltas_domains(rpars)
        assert directories == expect_directories

    def test_execution_suppressed(self, make_rpars, mocks, mock_tasks, caplog):
        """Check expected behavior when only inputs are collected."""
        caplog.set_level(logging.INFO)
        rpars = make_rpars(1)
        rpars.SUPPRESS_EXECUTION = True
        mock_tasks(n_comptasks=(10,), n_runtasks=(35,))
        deltas_domains(rpars)
        rpars.updateCores.assert_not_called()
        rpars.setHaltingLevel.assert_called_once_with(3)
        mocks['pool'].assert_not_called()
        mocks['copy_log'].assert_not_called()
        mocks['rmtree'].assert_not_called()
        mocks['checksums'].assert_not_called()
        expect_log = 'Getting input for delta calculations:'
        assert expect_log in caplog.text
        not_in_log = (
            'Compiling fortran files...',
            'Running delta calculations...',
            'Delta calculations finished.',
            )
        assert not any(msg in caplog.text for msg in not_in_log)

    def test_nothing_to_compile(self, make_rpars, mocks, mock_tasks, caplog):
        """Check that no compilation occurs if no domain requires it."""
        caplog.set_level(logging.INFO)
        n_domains = 3
        rpars = make_rpars(n_domains)
        mock_tasks(n_comptasks=(0,)*n_domains, n_runtasks=(1,)*n_domains)
        deltas_domains(rpars)
        rpars.updateCores.assert_called_once()
        assert mocks['one_domain'].call_count == n_domains
        n_parallel = 1  # only run
        assert mocks['pool'].call_count == n_parallel
        mocks['copy_log'].assert_not_called()
        mocks['rmtree'].assert_not_called()
        mocks['checksums'].assert_not_called()
        not_in_log = (
            'Compiling fortran files...',
            )
        assert not any(msg in caplog.text for msg in not_in_log)

    def test_nothing_to_execute(self, make_rpars, mocks, mock_tasks, caplog):
        """Check that no calculation occurs if no domain requires it."""
        caplog.set_level(logging.INFO)
        n_domains = 95
        rpars = make_rpars(n_domains)
        mock_tasks(n_comptasks=(1,)*n_domains, n_runtasks=(0,)*n_domains)
        deltas_domains(rpars)
        rpars.updateCores.assert_called_once()
        assert mocks['one_domain'].call_count == n_domains
        n_parallel = 1  # only compilation
        assert mocks['pool'].call_count == n_parallel
        assert mocks['copy_log'].call_count == n_domains  # 1 task each
        assert mocks['rmtree'].call_count == n_domains    # 1 task each
        mocks['checksums'].assert_called_once()
        not_in_log = (
            'Running delta calculations...',
            'Delta calculations finished.',
            )
        assert not any(msg in caplog.text for msg in not_in_log)

    def test_single_domain_returns_early(self, make_rpars, mocks):
        """Check no execution occurs if single-domain runs return early."""
        rpars = make_rpars(3)
        deltas_domains(rpars)
        mocks['pool'].assert_not_called()
        mocks['copy_log'].assert_not_called()
        mocks['rmtree'].assert_not_called()
        mocks['checksums'].assert_not_called()

    def test_stop_after_compile(self, make_rpars, mocks, mock_tasks, caplog):
        """Check effect of stopping right after compilation."""
        caplog.set_level(logging.INFO)
        n_domains = 12
        rpars = make_rpars(n_domains)
        mock_tasks(n_comptasks=(1,)*n_domains, n_runtasks=(3,)*n_domains)
        def _set_stop(*_, **__):
            rpars.STOP = True
        mocks['pool'].side_effect = _set_stop
        deltas_domains(rpars)
        assert mocks['pool'].call_count == 1  # Only compilation
        mocks['checksums'].assert_called_once()
        mocks['copy_log'].assert_not_called()
        mocks['rmtree'].assert_not_called()
        expect_log = 'Compiling fortran files...'
        not_in_log = (
            'Running delta calculations...',
            'Delta calculations finished.',
            )
        assert expect_log in caplog.text
        assert not any(msg in caplog.text for msg in not_in_log)

    def test_stop_after_run(self, make_rpars, mocks, mock_tasks, caplog):
        """Check effect of stopping right after execution."""
        caplog.set_level(logging.INFO)
        n_domains = 2
        rpars = make_rpars(n_domains)
        mock_tasks(n_comptasks=(1,)*n_domains, n_runtasks=(3,)*n_domains)
        def _set_stop(*args, **__):
            if args[2] is run_delta:
                rpars.STOP = True
        mocks['pool'].side_effect = _set_stop
        deltas_domains(rpars)
        n_parallel = 2
        assert mocks['pool'].call_count == n_parallel
        mocks['checksums'].assert_called_once()
        mocks['copy_log'].assert_not_called()
        mocks['rmtree'].assert_not_called()
        expect_log = (
            'Compiling fortran files...',
            'Running delta calculations...',
            )
        not_in_log = (
            'Delta calculations finished.',
            )
        assert all(msg in caplog.text for msg in expect_log)
        assert not any(msg in caplog.text for msg in not_in_log)


class TestDeltasDomainsRaises:
    """Tests for conditions that cause the deltas_domains function to fail."""

    @use('mocks')
    def test_find_compiler_fails(self, make_rpars, mock_tasks, caplog):
        """Check complaints if no compiler is available."""
        n_domains = 5
        rpars = make_rpars(n_domains)
        rpars.getFortranComp.side_effect = Exception('find compiler failed')
        mock_tasks(n_comptasks=(1,)*n_domains, n_runtasks=(2,)*n_domains)
        exc_msg = 'No Fortran compiler'
        with pytest.raises(RuntimeError, match=exc_msg):
            deltas_domains(rpars)
        expect_log = 'No fortran compiler found, cancelling...'
        assert expect_log in caplog.text

    def test_cleanup_fails(self, make_rpars, mocks, mock_tasks, caplog):
        """Check logging if cleaning up compile folders fails."""
        n_domains = 4
        rpars = make_rpars(n_domains)
        n_comptasks = (3,)*n_domains
        mock_tasks(n_comptasks, n_runtasks=(1,)*n_domains)
        mocks['rmtree'].side_effect = Exception
        deltas_domains(rpars)
        n_tasks_total = sum(n_comptasks)
        assert mocks['copy_log'].call_count == n_tasks_total
        assert mocks['rmtree'].call_count == n_tasks_total
        expect_log = 'Error deleting delta compile folder'
        assert caplog.text.count(expect_log) == n_tasks_total

    def test_inner_call_raises(self, make_rpars, mocks, caplog):
        """Check complaints if single-domain runs fail unexpectedly."""
        rpars = make_rpars(5)
        exc_msg = 'First domain fails'
        mocks['one_domain'].side_effect = Exception(exc_msg)
        with pytest.raises(Exception, match=exc_msg):
            deltas_domains(rpars)
        expect_log = 'Error while creating delta input for'
        assert expect_log in caplog.text

    def test_inner_call_unexpected_return(self, make_rpars, mocks):
        """Check complaints if single-domain runs fail unexpectedly."""
        rpars = make_rpars(4)
        mocks['one_domain'].return_value = 'bad value'
        exc_msg = 'Unknown error while creating delta input'
        with pytest.raises(RuntimeError, match=exc_msg):
            deltas_domains(rpars)
