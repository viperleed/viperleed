"""Tests for module domain_finder of viperleed.calc.bookkeeper."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2025-04-10'
__license__ = 'GPLv3+'

from pathlib import Path
import re

import pytest
from pytest_cases import fixture
from pytest_cases import parametrize

from viperleed.calc.bookkeeper.domain_finder import DomainFinder
from viperleed.calc.bookkeeper.domain_finder import MainPathNotFoundError
from viperleed.calc.bookkeeper.history.meta import DomainInfo
from viperleed.calc.constants import DEFAULT_DELTAS
from viperleed.calc.constants import DEFAULT_OUT
from viperleed.calc.constants import DEFAULT_SUPP
from viperleed.calc.constants import DEFAULT_TENSORS
from viperleed.calc.constants import LOG_PREFIX

from ...helpers import filesystem_from_dict


_MODULE = 'viperleed.calc.bookkeeper.domain_finder'


@fixture(name='bookkeeper')
def fixture_bookkeeper(mocker):
    """Return a fake bookkeeper."""
    return mocker.MagicMock()


@fixture(name='finder')
def fixture_finder(bookkeeper):
    """Return a DomainFinder with a fake bookkeeper."""
    finder = DomainFinder(bookkeeper)
    # The net one is just a convenience attribute so we don't need to
    # repeat over and over both the finder and the bookkeeper fixtures
    # in the tests below.
    finder.bookkeeper = bookkeeper
    return finder


@fixture(name='mock_finder_info')
def fixture_mock_finder_info(finder):
    """Set _domain_info to the given dict."""
    def _mock(info):
        # pylint: disable-next=protected-access           # OK in tests
        finder._domain_info = info
        return finder
    return _mock


class TestDomainFinder:
    """Tests for the DomainFinder class."""

    def test_collect_info_last_folder(self, finder, mocker):
        """Check the collection of domain information with a history folder."""
        mock_info = {'main': mocker.MagicMock()}
        bookkeeper = finder.bookkeeper
        bookkeeper.history.last_folder.metadata = mocker.MagicMock(
            domains=mock_info,
            )
        finder.collect_info()
        bookkeeper.update_from_cwd.assert_called_once_with(silent=True)
        assert finder.domain_info is mock_info['main']

    def test_collect_info_no_history(self, finder):
        """Check the collection of domain information without history."""
        bookkeeper = finder.bookkeeper
        bookkeeper.history.last_folder = None
        finder.collect_info()
        bookkeeper.update_from_cwd.assert_called_once_with(silent=True)
        assert finder.domain_info is None

    _find_domains = {  # _domain_info, called, result
        'not a domain calc': ({}, (), ()),
        'subdomain': ({'main': ('main_path', 'main_hash')},
                      ('_warn_running_in_subdomain',),
                      ()),
        }

    @parametrize('info,called,expect',
                 _find_domains.values(),
                 ids=_find_domains)
    # pylint: disable-next=too-many-arguments  # 2/6 are fixtures
    def test_find_domains(self, info, called, expect,
                          mock_finder_info, mocker):
        """Check the subcalls in find_domains."""
        finder = mock_finder_info(info)
        mocks = [mocker.patch.object(finder, attr) for attr in called]
        result = finder.find_domains(mocker.MagicMock())
        assert result == expect
        for mock in mocks:
            mock.assert_called_once()

    def test_find_domains_calls_from_main(self, mock_finder_info, mocker):
        """Check the subcalls in find_domains when called in a domain main."""
        finder = mock_finder_info(
            {'domains': (('domain_path', 'domain_hash'),)}
            )
        mock_from_main = mocker.patch.object(finder, '_find_domains_from_main')
        result = finder.find_domains(mocker.MagicMock())
        assert result is mock_from_main.return_value
        mock_from_main.assert_called_once()

    def test_find_potential_domains(self, finder, tmp_path, caplog):
        """Test the find_potential_domains method."""
        caplog.set_level(0)  # All messages
        not_domains = {
            DEFAULT_DELTAS: {},   # We don't explicitly check for this
            DEFAULT_OUT: {},
            DEFAULT_SUPP: {},
            DEFAULT_TENSORS: {},  # We don't explicitly check for this
            'history': {},
            'workhistory': {},
            'another_folder': {},
            }
        domains = {
            'Domain_1': {},    # An automatically-labeled domain
            'Domain_two': {},  # A domain with a user label
            'domain_after_calc': {
                # A root domain with all the expected contents at the
                # end of a calc run, not yet processed by bookkeeper.
                # (Having only one of the contents would be enough.)
                DEFAULT_OUT: {},
                DEFAULT_SUPP: {},
                },
            'domain_with_workhistory': {
                # A root domain with a workhistory folder,
                # not yet processed by bookkeeper.
                'workhistory': {},
                },
            'domain_after_calc_with_bookkeeper': {
                # A root domain already processed by bookkeeper.
                # (Having only one of the contents would be enough.)
                'history': {},
                'history.info': '',
                },
            'domain_with_log_file': {
                # A root domain in which calc was run manually.
                f'{LOG_PREFIX}_timestamp.log': '',
                },
            }
        filesystem_from_dict({**not_domains, **domains}, tmp_path)
        finder.bookkeeper.cwd = tmp_path
        domains_found = finder.find_potential_domains()
        assert not caplog.text
        assert len(domains_found) == len(domains)
        assert {d.name for d in domains_found} == set(domains)

    _is_subdomain = {
        'not a domain calc': ({}, False),
        'main domain folder': (
            {'domains': (('sub_path_1', 'sub_hash_1'),
                         ('sub_path_2', 'sub_hash_2'))},
            False,
            ),
        'subdomain folder': ({'main': ('path', 'hash')}, True),
        }

    @parametrize('info,expect', _is_subdomain.values(), ids=_is_subdomain)
    def test_is_subdomain(self, info, expect, mock_finder_info):
        """Check the expected value of the is_subdomain property."""
        finder = mock_finder_info(info)
        assert finder.is_subdomain == expect


class TestDomainFinderFromMain:
    """Tests for the _find_domains_from_main method."""

    @fixture(name='find')
    def fixture_find(self, mock_finder_info, existing_domains):
        """Return the subdomains found."""
        def _find(info=None):
            if info is None:
                info = {'domains': existing_domains}
            finder = mock_finder_info(info)
            # pylint: disable-next=protected-access       # OK in tests
            return finder._find_domains_from_main()
        return _find

    @fixture(name='existing_domains')
    def fixture_existing_domains(self):
        """Return information of some existing domains."""
        return (
            DomainInfo('existing_domain_path_one', 'existing_hash_one'),
            DomainInfo('existing_domain_path_two', 'existing_hash_two'),
            )

    def test_domain_folder_not_found(self, find, caplog):
        """Check the expected result when a domain subfolder is not found."""
        domains = find({'domains': (
            ('non_exiting_path', 'hash_of_non_existing'),
            )})
        assert not domains
        assert caplog.text

    def test_main_path_mismatched(self, finder, find, caplog, mocker):
        """Check result when the main path stored in history is mismatched."""
        def _mock_get_history_folder(*_):
            """Return a fake history folder."""
            folder = mocker.MagicMock()
            folder.metadata.domains = {'main': ('not-the-same-main-path',
                                                'main_hash')}
            return folder

        mocker.patch.object(finder,
                            '_get_history_folder',
                            _mock_get_history_folder)
        mock_ask_user = mocker.patch(f'{_MODULE}.ask_user_confirmation',
                                     return_value=True)
        assert find()
        assert caplog.text
        # Ask user only once, even if there's 2 mismatched
        # folders (see fixture_existing_domains).
        mock_ask_user.assert_called_once()

    def test_main_path_mismatched_user_says_no(self, finder, find, mocker):
        """Check behavior when users don't want to proceed."""
        check_info = mocker.patch.object(
            finder,
            '_check_subdomain_info_consistent',
            side_effect=MainPathNotFoundError(None),
            )
        mocker.patch(f'{_MODULE}.ask_user_confirmation', return_value=False)
        with pytest.raises(MainPathNotFoundError):
            find()
        check_info.assert_called_once()  # Even if there's 2 folders

    def test_no_domain_history_folder(self, finder, find, caplog, mocker):
        """Check the expected result when the history folder is not found."""
        mock_get_history = mocker.patch.object(finder,
                                               '_get_history_folder',
                                               return_value=None)
        domains = find({'domains': (
            DomainInfo('existing_domain_path',
                       'hash_of_non_existing_history_folder'),
            )})
        mock_get_history.assert_called_once_with(
            Path('existing_domain_path'),
            'hash_of_non_existing_history_folder',
            )
        assert domains
        assert caplog.text

    # pylint: disable-next=too-many-arguments    # All fixtures
    def test_success(self, existing_domains, finder, find,
                     tmp_path, caplog, mocker):
        """Test the successful detection of subdomains."""
        def _mock_get_history_folder(*_):
            """Return a fake history folder."""
            folder = mocker.MagicMock()
            folder.metadata.domains = {'main': (str(tmp_path), 'main_hash')}
            return folder

        mocker.patch.object(finder,
                            '_get_history_folder',
                            _mock_get_history_folder)
        finder.bookkeeper.cwd = tmp_path
        domains = find()
        assert not caplog.text
        assert domains
        assert all(isinstance(d, Path) for d in domains)
        assert all(d.name == i.path for d, i in zip(domains, existing_domains))


class TestDomainFinderFromSubdomain:
    """Tests for the _find_domains_from_subdomain method."""

    @fixture(name='find')
    def fixture_find(self, mock_finder_info):
        """Call _find_domains_from_subdomain."""
        def _find(info):
            finder = mock_finder_info(info)
            # pylint: disable-next=protected-access       # OK in tests
            return finder._find_domains_from_subdomain()
        return _find

    def test_main_path_not_found(self, finder, find, mocker):
        """Check complaints when the main path does not exist."""
        main_path = 'some_path'
        mocker.patch.object(finder,
                            '_get_history_folder',
                            side_effect=FileNotFoundError)
        with pytest.raises(MainPathNotFoundError, match=main_path):
            find({'main': DomainInfo(main_path, 'main_hash')})

    def test_no_history_in_main_path(self, finder, find, mocker):
        """Check complaints when the main path does not exist."""
        main_path = 'some_path'
        mocker.patch.object(finder,
                            '_get_history_folder',
                            return_value=None)
        assert not find({'main': DomainInfo(main_path, 'main_hash')})

    def test_success(self, finder, find, mocker):
        """Check that the expected domains are found."""
        main_path = 'some_path'
        mock_main_folder = mocker.MagicMock()
        all_domains = (
            DomainInfo('domain_1', 'hash_1'),
            DomainInfo('domain_2', 'hash_2'),
            DomainInfo('domain_3', 'hash_3'),
            DomainInfo('domain_4', 'hash_4'),
            )
        mock_main_folder.metadata.domains = {'domains': all_domains}
        mocker.patch.object(finder,
                            '_get_history_folder',
                            return_value=mock_main_folder)
        domains = find({'main': DomainInfo(main_path, 'main_hash')})
        assert domains
        assert len(domains) == len(all_domains)
        assert all(isinstance(d, Path) for d in domains)
        assert (d.name == i.path for d, i in zip(domains, all_domains))


class TestDomainFinderGetHistoryFolder:
    """Tests for the _get_history_folder method."""

    @fixture(name='mock_explorer')
    def fixture_patch_explorer(self, mocker):
        """Replace HistoryExplorer with a mock."""
        explorer_instance = mocker.MagicMock()
        explorer_cls = mocker.patch(f'{_MODULE}.HistoryExplorer',
                                    return_value=explorer_instance)
        return explorer_cls, explorer_instance

    @fixture(name='get')
    def fixture_get(self, finder):
        """Call _get_history_folder."""
        def _call(*args):
            # pylint: disable-next=protected-access       # OK in tests
            return finder._get_history_folder(*args)
        return _call

    def test_root_not_found(self, get, mocker):
        """Check complaints when the root path does not exist."""
        with pytest.raises(FileNotFoundError):
            get('not_a_root_folder', mocker.MagicMock())

    def test_success(self, get, mock_explorer, mocker):
        """Check the successful fetching of a history folder."""
        explorer_cls, explorer = mock_explorer
        mock_root = mocker.MagicMock(spec=Path)
        mock_hash = mocker.MagicMock()
        result = get(mock_root, mock_hash)
        assert result is explorer.subfolder_from_hash.return_value
        explorer_cls.assert_called_once_with(Path(mock_root))
        explorer.collect_subfolders.assert_called_once_with()
        explorer.subfolder_from_hash.assert_called_once_with(mock_hash)


class TestDomainFinderWarnSubdomain:
    """Tests for the _warn_running_in_subdomain helper method."""

    @fixture(name='warn')
    def fixture_warn(self, finder, mocker):
        """Call _warn_running_in_subdomain."""
        def _call(main_path='some_path'):
            # pylint: disable-next=protected-access       # OK in tests
            finder._warn_running_in_subdomain(mocker.MagicMock(), main_path)
        return _call

    # pylint: disable-next=too-many-arguments  # All fixtures
    def test_main_path_exists(self, finder, warn, tmp_path, caplog, mocker):
        """Check complaints when the main path is mismatched."""
        all_domains = Path('this_domain'), Path('domain_1'), Path('domain_2')
        mocker.patch.object(finder,
                            '_find_domains_from_subdomain',
                            return_value=all_domains)
        finder.bookkeeper.cwd = tmp_path/all_domains[0]
        warn()
        expect_log = (
            re.compile('Running bookkeeper in mode.* in a domain subfolder'),
            re.compile(r'also execute \'bookkeeper.*\' in the main path'),
            re.compile(r'other domain subfolders \(domain_1 and domain_2\)'),
            )
        log = caplog.text
        assert all(expect.search(log) for expect in expect_log)

    def test_main_path_not_found(self, finder, warn, caplog, mocker):
        """Check complaints when the main path is mismatched."""
        mocker.patch.object(finder,
                            '_find_domains_from_subdomain',
                            side_effect=MainPathNotFoundError(None))
        warn()
        expect_log = (
            re.compile('Running bookkeeper in mode.* in a domain subfolder'),
            re.compile('path to the .*main directory could not be identified'),
            )
        log = caplog.text
        assert all(expect.search(log) for expect in expect_log)
