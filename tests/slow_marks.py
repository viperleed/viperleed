"""Pytest plugin marking slow and TensErLEED tests of the suite."""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-10-27'
__license__ = 'GPLv3+'

import collections
import sys

import pytest

from . import tiers

_COLLECTED_IDS = pytest.StashKey()
_STALE_MESSAGE = pytest.StashKey()  # On the controller, with xdist
_NODEID_SEP = '::'
STALE_KEY = 'stale_slow_cases'  # In pytest-xdist's workeroutput


def _split_case_id(nodeid):
    """Return the test part of `nodeid` and its '-'-padded case id."""
    test, bracket, params = nodeid.partition('[')
    if not bracket:
        return test, ''
    # Pad, to match whole '-'-separated parts of the case id
    return test, f'-{params[:-1]}-'


def _is_slow(item):
    """Return whether `item` is one of tiers.SLOW_CASES."""
    test, params = _split_case_id(item.nodeid)
    slow_parts = tiers.SLOW_CASES.get(test, ())
    return any(f'-{part}-' in params for part in slow_parts)


def pytest_itemcollected(item):
    """Remember the node ids of all items, before any deselection."""
    item.config.stash.setdefault(_COLLECTED_IDS, []).append(item.nodeid)


def _find_stale_slow_cases(nodeids):
    """Return descriptions of SLOW_CASES entries matching no node id.

    Only entries of modules with at least one collected node id are
    checked, so that running a subset of the suite is not an error.

    Parameters
    ----------
    nodeids : Iterable
        The node ids of all collected items.

    Returns
    -------
    stale : list
        One string per entry (test) or case part that matches nothing.
    """
    collected = collections.defaultdict(list)  # test: case ids
    modules = set()
    for nodeid in nodeids:
        test, params = _split_case_id(nodeid)
        modules.add(test.partition(_NODEID_SEP)[0])
        collected[test].append(params)
    stale = []
    for test, parts in tiers.SLOW_CASES.items():
        if test.partition(_NODEID_SEP)[0] not in modules:
            continue
        if test not in collected:
            stale.append(f'{test} (no such test)')
            continue
        stale.extend(
            f'{test}[{part}] (no such case)'
            for part in parts
            if not any(f'-{part}-' in p for p in collected[test])
            )
    return stale


def _stale_slow_cases_message(config):
    """Return an error message about stale SLOW_CASES, or None."""
    if any(_NODEID_SEP in arg for arg in config.args):
        return None  # Tests selected by node id: only part of a module
    if config.getoption('lf', default=False):
        return None  # --lf drops the passed items of a module early
    stale = _find_stale_slow_cases(config.stash.get(_COLLECTED_IDS, ()))
    if not stale:
        return None
    entries = '\n  '.join(stale)
    return (
        f'Stale entries in SLOW_CASES ({tiers.__file__}), which lists '
        f'the cases that take {tiers.SLOW_THRESHOLD} s or more: the '
        f'tests or cases below do not exist (renamed or removed?). '
        f'Update the module.\n  {entries}'
        )


def _check_slow_cases(config, items):
    """Report stale tiers.SLOW_CASES as a usage error.

    A pytest-xdist worker cannot raise it: that ends the controller
    with an internal error. Instead, a worker drops all its items
    (so that no test runs) and hands the message over to the
    controller, which reports it and ends with the same exit status
    as a serial run (pytest_testnodedown, pytest_sessionfinish).
    """
    message = _stale_slow_cases_message(config)
    if not message:
        return
    workeroutput = getattr(config, 'workeroutput', None)
    if workeroutput is None:
        raise pytest.UsageError(message)
    workeroutput[STALE_KEY] = message
    items.clear()


@pytest.hookimpl(optionalhook=True)  # Only exists with pytest-xdist
def pytest_testnodedown(node):
    """Remember the stale-entry message of a pytest-xdist worker."""
    message = getattr(node, 'workeroutput', {}).get(STALE_KEY)
    if message:
        node.config.stash[_STALE_MESSAGE] = message


def pytest_sessionfinish(session):
    """End a pytest-xdist run with stale SLOW_CASES as a usage error."""
    message = session.config.stash.get(_STALE_MESSAGE, None)
    if message:
        sys.stderr.write(f'ERROR: {message}\n')
        session.exitstatus = pytest.ExitCode.USAGE_ERROR


def pytest_collection_modifyitems(config, items):
    """Mark TensErLEED and slow tests; check SLOW_CASES."""
    _check_slow_cases(config, items)
    for item in items:
        if tiers.TENSORLEED_FIXTURES.intersection(item.fixturenames):
            item.add_marker(pytest.mark.tensorleed)
        if _is_slow(item):
            item.add_marker(pytest.mark.slow)
