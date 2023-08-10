"""Tests for the viperleed command line interface.
"""
import argparse

import pytest

from viperleed.calc.__main__ import add_calc_parser_arguments
from viperleed.calc.bookkeeper import bookkeeper_cli_options
from viperleed.calc.bookkeeper import main as bookkeeper_main

authors = ["Alexander M. Imre"]
__created__ = "2023-08-02"


@pytest.fixture
def calc_parser():
    parser = argparse.ArgumentParser()
    add_calc_parser_arguments(parser)
    return parser

@pytest.fixture
def bookkeeper_parser():
    parser = argparse.ArgumentParser()
    bookkeeper_cli_options(parser)
    return parser

class TestCalcParser:
    def test_viperleed_calc_parse_version(self, calc_parser):
        with pytest.raises(SystemExit):
            calc_parser.parse_args(["--version"])

    def test_viperleed_calc_parse_work(self, calc_parser, tmp_path):
        assert calc_parser.parse_args(["-w", str(tmp_path)]).work == str(tmp_path)

    @pytest.mark.parametrize("v_flag", ["-v", "--verbose"])
    def test_viperleed_calc_parse_verbose(self, calc_parser, v_flag):
        assert calc_parser.parse_args([v_flag,]).verbose


    @pytest.mark.parametrize("v_flag", ["-vv", "--very-verbose"])
    def test_viperleed_calc_parse_very_verbose(self, calc_parser, v_flag):
        assert calc_parser.parse_args([v_flag,]).very_verbose


class TestBookkeeperParser:
    @pytest.mark.parametrize("flag", ["-j", "--job-name"])
    def test_bookkeeper_parser_job_name(self, bookkeeper_parser, flag):
        assert bookkeeper_parser.parse_args([flag, "test"]).job_name == "test"

    @pytest.mark.parametrize("flag", ["-c", "--cont"])
    def test_bookkeeper_parser_continue(self, bookkeeper_parser, flag):
        assert bookkeeper_parser.parse_args([flag,]).cont

    @pytest.mark.parametrize("flag", ["-d", "--discard"])
    def test_bookkeeper_parser_discard(self, bookkeeper_parser, flag):
        assert bookkeeper_parser.parse_args([flag,]).discard

    def test_bookkeeper_parser_exclusive_options(self, bookkeeper_parser):
        args = bookkeeper_parser.parse_args(["-c", "-d"])
        with pytest.raises(RuntimeError):
            bookkeeper_main(args)