"""Tests for the viperleed bookkeeper command-line interface."""

import pytest
from pytest_cases import fixture, parametrize

from viperleed.calc.bookkeeper.bookkeeper import Bookkeeper
from viperleed.calc.bookkeeper.cli import BookkeeperCLI
from viperleed.calc.bookkeeper.cli import StoreBookkeeperMode
from viperleed.calc.bookkeeper.mode import BookkeeperMode


_ARGPARSE_EXIT_WITH_ERROR = 2


@fixture(name='bookkeeper_parser')
def fixture_bookkeeper_parser():
    """Return a CLI argument parser for viperleed.calc.bookkeeper."""
    return BookkeeperCLI().parser


class TestBookkeeperParser:
    """Tests for parsing CLI arguments of viperleed.calc.bookkeeper."""

    @parametrize(flag=('-a', '--archive'))
    def test_archive(self, bookkeeper_parser, flag):
        """Check interpretation of --archive mode."""
        parsed = bookkeeper_parser.parse_args([flag,])
        assert parsed.archive
        assert parsed.mode is BookkeeperMode.ARCHIVE

    @parametrize(flag=('-c', '--clear'))
    def test_clear(self, bookkeeper_parser, flag):
        """Check interpretation of --clear mode."""
        parsed = bookkeeper_parser.parse_args([flag,])
        assert parsed.clear
        assert parsed.mode is BookkeeperMode.CLEAR

    @parametrize(flag=('-d', '--discard'))
    def test_discard(self, bookkeeper_parser, flag):
        """Check interpretation of --discard mode."""
        parsed = bookkeeper_parser.parse_args([flag,])
        assert parsed.discard
        assert parsed.mode is BookkeeperMode.DISCARD

    @parametrize(flag=('-df', '--discard-full'))
    def test_discard_full(self, bookkeeper_parser, flag):
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
    def test_exclusive_options(self, bookkeeper_parser, flags):
        """Check complaints when conflicting modes are given."""
        with pytest.raises(SystemExit) as exc:
            bookkeeper_parser.parse_args(flags)
        assert exc.value.code == _ARGPARSE_EXIT_WITH_ERROR

    @parametrize('flag,expect', (('-df', False), ('-y', True)))
    def test_user_confirmation_required(self, bookkeeper_parser, flag, expect):
        """Check correct interpretation of the '-y' CLI argument."""
        parsed = bookkeeper_parser.parse_args([flag,])
        assert parsed.skip_confirmation == expect


def test_action_raises_unexpected_mode(bookkeeper_parser, mocker):
    """Check complaints when called with an unexpected mode."""
    action = StoreBookkeeperMode(option_strings=None, dest='unknown')
    with pytest.raises(SystemExit):
        action(bookkeeper_parser, mocker.MagicMock(), None)


class TestBookkeeperRun:
    """Tests for calling the BookkeeperCLI."""

    @fixture(name='mock_run')
    def fixture_mock_run(self, mocker):
        """Return a BookkeeperCLI with its .run method replaced."""
        return mocker.patch.object(Bookkeeper, 'run')

    _mode = {
        '': (BookkeeperMode.ARCHIVE,  # The default
             {'requires_user_confirmation': True}),
        ('--discard',): (BookkeeperMode.DISCARD,
                         {'requires_user_confirmation': True}),
        ('-df',): (BookkeeperMode.DISCARD_FULL,
                   {'requires_user_confirmation': True}),
        ('-df', '-y'): (BookkeeperMode.DISCARD_FULL,
                        {'requires_user_confirmation': False}),
        }

    @parametrize('args,expect', _mode.items())
    def test_bookkeeper_run(self, args, expect, mock_run):
        """Check correct call to bookkeeper.run."""
        bookie = BookkeeperCLI()
        bookie(args)
        mode, kwargs = expect
        mock_run.assert_called_with(mode, **kwargs)
