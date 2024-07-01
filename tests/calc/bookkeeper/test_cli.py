"""Tests for the viperleed bookkeeper command-line interface."""

import pytest
from pytest_cases import fixture, parametrize

from viperleed.calc.bookkeeper.cli import BookkeeperCLI
from viperleed.calc.bookkeeper.mode import BookkeeperMode


@fixture(name='bookkeeper_parser')
def fixture_bookkeeper_parser():
    """Return a CLI argument parser for viperleed.calc.bookkeeper."""
    return BookkeeperCLI().parser


class TestBookkeeperParser:
    """Tests for parsing CLI arguments of viperleed.calc.bookkeeper."""

    @parametrize(flag=('-j', '--job-name'))
    def test_parser_job_name(self, bookkeeper_parser, flag):
        """Check interpretation of --job-name flag."""
        job = 'test'
        assert bookkeeper_parser.parse_args([flag, job]).job_name == job

    @parametrize(flag=('-a', '--archive'))
    def test_parser_archive(self, bookkeeper_parser, flag):
        """Check interpretation of --archive mode."""
        parsed = bookkeeper_parser.parse_args([flag,])
        assert parsed.archive
        assert parsed.mode is BookkeeperMode.ARCHIVE

    @parametrize(flag=('-c', '--clear'))
    def test_parser_clear(self, bookkeeper_parser, flag):
        """Check interpretation of --clear mode."""
        parsed = bookkeeper_parser.parse_args([flag,])
        assert parsed.clear
        assert parsed.mode is BookkeeperMode.CLEAR

    @parametrize(flag=('-d', '--discard'))
    def test_parser_discard(self, bookkeeper_parser, flag):
        """Check interpretation of --discard mode."""
        parsed = bookkeeper_parser.parse_args([flag,])
        assert parsed.discard
        assert parsed.mode is BookkeeperMode.DISCARD

    @parametrize(flag=('-df', '--discard-full'))
    def test_parser_discard_full(self, bookkeeper_parser, flag):
        """Check interpretation of --discard-full mode."""
        parsed = bookkeeper_parser.parse_args([flag,])
        assert parsed.discard_full
        assert parsed.mode is BookkeeperMode.DISCARD_FULL

    _exclusive = {
        'archive & clear': ('-a', '-c'),
        'discard & clear': ('-d', '-c'),
        'discard-full & archive': ('-df', '-a'),
        }

    @parametrize(flags=_exclusive.values(), ids=_exclusive)
    def test_parser_exclusive_options(self, bookkeeper_parser, flags):
        """Check complaints when conflicting modes are given."""
        with pytest.raises(SystemExit) as exc:
            bookkeeper_parser.parse_args(flags)
        assert exc.value.code == 2  # argparse error code
