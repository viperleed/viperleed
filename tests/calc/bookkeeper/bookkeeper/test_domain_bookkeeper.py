"""Tests for DomainBookkeeper of viperleed.calc.bookkeeper.bookkeeper."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2025-04-11'
__license__ = 'GPLv3+'

from pathlib import Path

from viperleed.calc.bookkeeper.bookkeeper import DomainBookkeeper


_MODULE = 'viperleed.calc.bookkeeper.bookkeeper'


class TestDomainBookkeeper:
    """Tests for the DomainBookkeeper class."""

    def test_clean_state(self, mocker):
        """Test the _clean_state private method."""
        mock_path = mocker.MagicMock(spec=Path)
        mock_domain_root_instance = mocker.MagicMock(path=mock_path)
        mock_domain_root_cls = mocker.patch(
            f'{_MODULE}.DomainRootExplorer',
            return_value=mock_domain_root_instance,
            )
        mock_main_root = mocker.MagicMock()
        domain_bookie = DomainBookkeeper(mock_main_root,
                                         cwd=mock_path)
        # pylint: disable=protected-access                # OK in tests
        domain_bookie._state_info['logger_prepared'] = log = mocker.MagicMock()
        domain_bookie._clean_state()

        # Check unchanged attributes
        assert domain_bookie._main_root is mock_main_root
        assert domain_bookie._state_info['logger_prepared'] is log
        # DomainRootExplorer should have been called twice, once
        # at creation of the domain_bookie instance above, once
        # as part of the _clean_state call.
        n_domain_root_calls = 2
        assert mock_domain_root_cls.call_count == n_domain_root_calls

    def test_init(self, mocker):
        """Test initialization."""
        mock_domain_root_instance = mocker.MagicMock()
        mock_domain_root_cls = mocker.patch(
            f'{_MODULE}.DomainRootExplorer',
            return_value=mock_domain_root_instance,
            )
        mock_main_root = mocker.MagicMock()
        mock_super_root_instance = mocker.MagicMock()
        mock_super_root_cls = mocker.patch(
            f'{_MODULE}.RootExplorer',
            return_value=mock_super_root_instance,
            )
        mock_path = mocker.MagicMock(spec=Path)
        domain_bookie = DomainBookkeeper(mock_main_root,
                                         cwd=mock_path)

        # Check delegation to super().__init__
        mock_super_root_cls.assert_called_once_with(
            path=mock_path,
            bookkeeper=domain_bookie,
            )
        # As well as replacement of the _root attribute with
        # a DomainRootExplorer instead of a RootExplorer.
        mock_domain_root_cls.assert_called_once_with(
            mock_super_root_instance.path,
            domain_bookie,
            mock_main_root,
            )
        # pylint: disable=protected-access                # OK in tests
        assert domain_bookie._root is mock_domain_root_instance
        assert domain_bookie._main_root is mock_main_root

    def test_run_in_subdomain(self, tmp_path, mocker):
        """Check delegation to private methods of run_in_subdomain."""
        mock_main_root = mocker.MagicMock()
        domain_bookie = DomainBookkeeper(mock_main_root,
                                         cwd=tmp_path)
        mock_return = mocker.MagicMock()
        run_one = mocker.patch.object(domain_bookie,
                                      '_run_one_domain',
                                      return_value=mock_return)
        mock_args = [mocker.MagicMock() for _ in range(10)]
        mock_kwargs = {f'kwarg_{k}': mocker.MagicMock() for k in range(10)}
        result = domain_bookie.run_in_subdomain(*mock_args, **mock_kwargs)

        assert result is mock_return
        run_one.assert_called_once_with(*mock_args, **mock_kwargs)
