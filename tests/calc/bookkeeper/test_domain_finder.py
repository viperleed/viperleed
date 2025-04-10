"""Tests for module domain_finder of viperleed.calc.bookkeeper."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2025-04-10'
__license__ = 'GPLv3+'

from viperleed.calc.bookkeeper.domain_finder import DomainFinder
from viperleed.calc.constants import DEFAULT_DELTAS
from viperleed.calc.constants import DEFAULT_OUT
from viperleed.calc.constants import DEFAULT_SUPP
from viperleed.calc.constants import DEFAULT_TENSORS
from viperleed.calc.constants import LOG_PREFIX

from ...helpers import filesystem_from_dict


class TestDomainFinder:
    """Tests for the DomainFinder class."""

    def test_find_potential_domains(self, tmp_path, caplog, mocker):
        """Test the find_potential_domains method."""
        caplog.set_level(0)  # All messages
        not_domains = {
            DEFAULT_DELTAS: {},   # We don't explicitly check for this
            DEFAULT_OUT: {},
            DEFAULT_SUPP: {},
            DEFAULT_TENSORS: {},  # We don't explicitly check for this
            'history': {},
            'workhistory': {},
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
        mock_bookkeeper = mocker.MagicMock(cwd=tmp_path)
        finder = DomainFinder(mock_bookkeeper)
        domains_found = finder.find_potential_domains()
        assert not caplog.text
        assert len(domains_found) == len(domains)
        assert {d.name for d in domains_found} == set(domains)
