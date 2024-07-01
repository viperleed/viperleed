"""Tests for the viperleed.calc command-line interface."""

__authors__ = (
    'Alexander M. Imre (@amimre)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2023-08-02'
__license__ = 'GPLv3+'


import pytest
from pytest_cases import fixture, parametrize

from viperleed.calc.cli import ViPErLEEDCalcCLI


@fixture(name='calc_parser')
def fixture_calc_parser():
    """Return a CLI argument parser for viperleed.calc."""
    return ViPErLEEDCalcCLI().parser


class TestCalcParser:
    """Tests for parsing CLI arguments of viperleed.calc."""

    def test_viperleed_calc_parse_version(self, calc_parser):
        """Check that requesting the version exits afterwards."""
        with pytest.raises(SystemExit):
            calc_parser.parse_args(['--version'])

    def test_viperleed_calc_parse_work(self, calc_parser, tmp_path):
        """Check interpretation of -w flag."""
        parsed = calc_parser.parse_args(['-w', str(tmp_path)])
        assert parsed.work == str(tmp_path)

    @parametrize(v_flag=('-v', '--verbose'))
    def test_viperleed_calc_parse_verbose(self, calc_parser, v_flag):
        """Check interpretation of -v flag."""
        assert calc_parser.parse_args([v_flag,]).verbose

    @parametrize(v_flag=('-vv', '--very-verbose'))
    def test_viperleed_calc_parse_very_verbose(self, calc_parser, v_flag):
        """Check interpretation of -vv flag."""
        assert calc_parser.parse_args([v_flag,]).very_verbose
