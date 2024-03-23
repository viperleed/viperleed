"""Tests for the viperleed command-line interface."""

__authors__ = (
    'Alexander M. Imre (@amimre)',
    )
__created__ = '2023-08-02'

import argparse

import pytest
from pytest_cases import parametrize

from viperleed.calc.cli import add_calc_parser_arguments
from viperleed.calc.bookkeeper import bookkeeper_cli_options
from viperleed.calc.bookkeeper import main as bookkeeper_main


@pytest.fixture(name='calc_parser')
def fixture_calc_parser():
    """Return a CLI argument parser for viperleed.calc."""
    parser = argparse.ArgumentParser()
    add_calc_parser_arguments(parser)
    return parser


@pytest.fixture(name='bookkeeper_parser')
def fixture_bookkeeper_parser():
    """Return a CLI argument parser for viperleed.calc.bookkeeper."""
    parser = argparse.ArgumentParser()
    bookkeeper_cli_options(parser)
    return parser


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


class TestBookkeeperParser:
    """Tests for parsing CLI arguments of viperleed.calc.bookkeeper."""

    @parametrize(flag=('-j', '--job-name'))
    def test_bookkeeper_parser_job_name(self, bookkeeper_parser, flag):
        """Check interpretation of --job-name flag."""
        job = 'test'
        assert bookkeeper_parser.parse_args([flag, job]).job_name == job

    @parametrize(flag=('-c', '--cont'))
    def test_bookkeeper_parser_continue(self, bookkeeper_parser, flag):
        """Check interpretation of --cont mode."""
        assert bookkeeper_parser.parse_args([flag,]).cont

    @parametrize(flag=('-d', '--discard'))
    def test_bookkeeper_parser_discard(self, bookkeeper_parser, flag):
        """Check interpretation of --discard mode."""
        assert bookkeeper_parser.parse_args([flag,]).discard

    def test_bookkeeper_parser_exclusive_options(self, bookkeeper_parser):
        """Check complaints when conflicting modes are given."""
        args = bookkeeper_parser.parse_args(['-c', '-d'])
        with pytest.raises(RuntimeError):
            bookkeeper_main(args)
