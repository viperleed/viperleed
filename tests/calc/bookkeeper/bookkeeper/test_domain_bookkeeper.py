"""Tests for DomainBookkeeper of viperleed.calc.bookkeeper.bookkeeper."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2025-04-11'
__license__ = 'GPLv3+'

from operator import attrgetter
from pathlib import Path

import pytest
from pytest_cases import fixture
from pytest_cases import parametrize

from viperleed.calc.bookkeeper.bookkeeper import DomainBookkeeper

from .conftest import _MODULE


@fixture(name='make_domain_bookie')
def fixture_domain_bookie(mocker):
    """Return a DomainBookkeeper with a fake path and main root."""
    def _make():
        mock_cwd = mocker.MagicMock(spec=Path)
        mock_domain_root_instance = mocker.MagicMock(path=mock_cwd)
        mock_domain_root_cls = mocker.patch(
            f'{_MODULE}.DomainRootExplorer',
            return_value=mock_domain_root_instance,
            )
        mock_main_root = mocker.MagicMock()
        bookie = DomainBookkeeper(mock_main_root, cwd=mock_cwd)
        return bookie, mock_main_root, mock_domain_root_cls
    return _make


class TestDomainBookkeeper:
    """Tests for the DomainBookkeeper class."""

    def test_check_may_discard_full(self, make_domain_bookie, mocker):
        """Check private calls in check_may_discard_full."""
        domain_bookie, *_ = make_domain_bookie()
        private = mocker.patch.object(domain_bookie, '_check_may_discard_full')
        domain_bookie.check_may_discard_full()
        private.assert_called_once_with()

    def test_clean_state(self, make_domain_bookie, mocker):
        """Test the _clean_state private method."""
        domain_bookie, main_root, mock_domain_root_cls = make_domain_bookie()
        # pylint: disable=protected-access                # OK in tests
        domain_bookie._state_info['logger_prepared'] = log = mocker.MagicMock()
        domain_bookie._clean_state()

        # Check unchanged attributes
        assert domain_bookie._main_root is main_root
        assert domain_bookie._state_info['logger_prepared'] is log
        # DomainRootExplorer should have been called twice, once
        # at creation of the domain_bookie instance above, once
        # as part of the _clean_state call.
        n_domain_root_calls = 2
        assert mock_domain_root_cls.call_count == n_domain_root_calls

    def test_init(self, make_domain_bookie, mocker):
        """Test initialization."""
        mock_super_root_instance = mocker.MagicMock()
        mock_super_root_cls = mocker.patch(
            f'{_MODULE}.RootExplorer',
            return_value=mock_super_root_instance,
            )
        domain_bookie, main_root, mock_domain_root_cls = make_domain_bookie()

        # Check delegation to super().__init__
        mock_super_root_cls.assert_called_once_with(
            path=domain_bookie.cwd,
            bookkeeper=domain_bookie,
            )
        # As well as replacement of the _root attribute with
        # a DomainRootExplorer instead of a RootExplorer.
        mock_domain_root_cls.assert_called_once_with(
            mock_super_root_instance.path,
            domain_bookie,
            main_root,
            )
        mock_domain_root_instance = mock_domain_root_cls()
        # pylint: disable=protected-access                # OK in tests
        assert domain_bookie._root is mock_domain_root_instance
        assert domain_bookie._main_root is main_root

    @parametrize(archives=(True, False))
    # Indeed, there are quite a few local variables in the next test,
    # but splitting it up into separate bits does not make much sense,
    # as we would not really reuse them.
    # pylint: disable-next=too-many-locals
    def test_run_in_subdomain(self, archives, make_domain_bookie, mocker):
        """Check delegation to private methods of run_in_subdomain."""
        domain_bookie, mock_main_root, *_ = make_domain_bookie()
        mock_return = (
            mocker.MagicMock(),
            mocker.MagicMock() if archives else None,
            )
        mock_args = [mocker.MagicMock() for _ in range(10)]
        mock_kwargs = {f'kwarg_{k}': mocker.MagicMock() for k in range(10)}

        mocks = {
            'domain_folder': mock_return[1],
            'main_folder': mocker.MagicMock(),
            'run_one': mocker.patch.object(domain_bookie,
                                           '_run_one_domain',
                                           return_value=mock_return),
            'log': mocker.patch(f'{_MODULE}.log.remove_bookkeeper_logfile'),
            }
        calls = {
            'run_one': mocker.call(*mock_args, **mock_kwargs),
            'log': mocker.call(domain_bookie.history.path),
            }
        if archives:
            calls.update({
            'domain_folder.mark_as_domain': mocker.call(
                mock_main_root.path,
                mocks['main_folder'],
                ),
            'domain_folder.metadata.write': mocker.call(),
            })
        result = domain_bookie.run_in_subdomain(mocks['main_folder'],
                                                *mock_args,
                                                **mock_kwargs)
        assert result == mock_return
        for callee, call in calls.items():
            mock_name, *method = callee.split('.', maxsplit=1)
            mock = mocks[mock_name]
            if method:
                mock = attrgetter(method[0])(mock)
            assert mock.mock_calls == [call]

    def test_run_in_subdomain_removes_log(self, make_domain_bookie, mocker):
        """Check removal of log-file handler in case of failure."""
        domain_bookie, *_ = make_domain_bookie()
        mocker.patch.object(domain_bookie,
                            '_run_one_domain',
                            side_effect=BaseException)
        mock_log = mocker.patch(f'{_MODULE}.log.remove_bookkeeper_logfile')

        with pytest.raises(BaseException):
            domain_bookie.run_in_subdomain(None)
        mock_log.assert_called_once()
