"""Pytest plugin for the --rng-seed option of the test suite."""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-10-27'
__license__ = 'GPLv3+'

import random

import pytest

RNG_SEED_OPTION = '--rng-seed'
WORKER_SEED_KEY = 'rng_seed'  # In pytest-xdist's workerinput


def pytest_addoption(parser):
    """Add the --rng-seed option."""
    parser.addoption(
        RNG_SEED_OPTION,
        dest='rng_seed',
        metavar='INT|random',
        default=None,
        help=('Derive the seeds of random-input tests from this base '
              'integer, or from a random one if "random" (printed at '
              'the end of the run, to reproduce it). Known failures '
              '(xfail) of random-input tests then report their real '
              'outcome. For local exploration, not for CI.'),
        )


def pytest_configure(config):
    """Resolve --rng-seed once.

    The seed is resolved in the controller only. A pytest-xdist worker
    takes the base that the controller resolved and handed over
    (pytest_configure_node).
    """
    workerinput = getattr(config, 'workerinput', None)
    if workerinput is not None:
        config.option.rng_seed = workerinput[WORKER_SEED_KEY]
        return
    value = config.option.rng_seed
    if value is None:
        return  # No option
    if value == 'random':
        value = random.SystemRandom().randrange(2**32)
    try:
        config.option.rng_seed = int(value)
    except ValueError:
        raise pytest.UsageError(
            f'{RNG_SEED_OPTION} must be an integer or "random", '
            f'not {value!r}'
            ) from None


@pytest.hookimpl(optionalhook=True)  # Only exists with pytest-xdist
def pytest_configure_node(node):
    """Hand the resolved base seed over to a pytest-xdist worker."""
    node.workerinput[WORKER_SEED_KEY] = node.config.option.rng_seed


def _rng_seed_message(config):
    """Return a line with the base seed, or None if not used."""
    base = config.option.rng_seed
    if base is None:
        return None
    return f'rng-seed base: {base} (reproduce with --rng-seed={base})'


def pytest_report_header(config):
    """Report the base seed, if used."""
    return _rng_seed_message(config)


@pytest.hookimpl(tryfirst=True)  # Before sections of other plugins
def pytest_terminal_summary(terminalreporter, config):
    """Report the base seed, if used, also with -q."""
    message = _rng_seed_message(config)
    if message:
        terminalreporter.write_line(message)


def pytest_collection_modifyitems(config, items):
    """Remove the xfail marks of random-input tests with fresh seeds."""
    if config.option.rng_seed is None:
        return
    for item in items:
        if 'draws' in item.fixturenames:
            item.own_markers[:] = [
                mark for mark in item.own_markers if mark.name != 'xfail'
                ]
