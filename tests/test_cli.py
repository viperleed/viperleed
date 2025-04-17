"""Tests for module cli of viperleed."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2025-03-06'
__license__ = 'GPLv3+'

import io

from pytest_cases import fixture
from pytest_cases import parametrize

from viperleed.cli import ViPErLEEDMain
from viperleed.calc.lib.context import execute_in_dir

from .helpers import CustomTestException


def _collect_cli_names(cli_cls):
    """Yield the names of all children CLIs of `cli_cls`, recursively."""
    cli = cli_cls()
    if not cli.children:
        yield cli.cli_name
        return
    children_names = (_collect_cli_names(child)
                      for child in cli.children.values())
    children_names = (f'{cli.cli_name} {name}'
                      for children in children_names
                      for name in children)
    yield cli.cli_name  # Without any child, does --help
    yield from children_names


def _collect_concrete_clis(cli_cls):
    """Yield all concrete children CLIs of `cli_cls`, recursively."""
    cli = cli_cls()
    if not cli.children:
        yield cli_cls
        return
    yield from (concrete_cls
                for child_cls in cli.children.values()
                for concrete_cls in _collect_concrete_clis(child_cls))


ALL_CLI_NAMES = tuple(_collect_cli_names(ViPErLEEDMain))
CONCRETE_CLIS = {f'{c.__module__}.{c.__name__}': c
                 for c in _collect_concrete_clis(ViPErLEEDMain)}


@fixture(name='mock_stdin', autouse=True)
def fixture_mock_stdin(tmp_path, mocker):
    """Replace standard input with mocks, execute in tmp_path."""
    # Some of our utilities require user input
    mocker.patch('builtins.input', side_effect=CustomTestException)
    mocker.patch('sys.stdin', io.StringIO('user input on stdin'))
    with execute_in_dir(tmp_path):
        yield


def call_cli(cli_cls, argv=None):
    """Execute `cli_cls` with command-line arguments."""
    argv = argv or []
    cli = cli_cls()
    try:
        cli(argv)
    except (CustomTestException, SystemExit):
        pass


@parametrize(util=ALL_CLI_NAMES)
def test_can_call_main_viperleed(util):
    """Check that calling a viperleed `util` is successful.

    This test simulates a command-line call like 'viperleed poscar'.

    Parameters
    ----------
    util : str
        Space-separated command-line invocation of a viperleed CLI,
        e.g., 'viperleed calc'.

    Returns
    -------
    None.
    """
    _, *argv = util.split()  # strip "viperleed"
    call_cli(ViPErLEEDMain, argv)


@parametrize(cli_cls=CONCRETE_CLIS.values(), ids=CONCRETE_CLIS)
def test_call_child_util_cls(cli_cls):
    """Check the successful execution of a `cli_cls`.

    This test simulates a command-line call like
    '[python -m] viperleed.calc'.

    Parameters
    ----------
    cli_cls : ViPErLEEDCLI
        The utility to be called.

    Returns
    -------
    None.
    """
    call_cli(cli_cls)
