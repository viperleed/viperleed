"""Pytest plugin marking slow, release, and TensErLEED tests."""

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

_COLLECTED = pytest.StashKey()
_STALE_MESSAGE = pytest.StashKey()  # On the controller, with xdist
_NODEID_SEP = '::'
STALE_KEY = 'stale_tier_cases'  # In pytest-xdist's workeroutput


def _split_case_id(nodeid):
    """Return the test part of `nodeid` and its '-'-padded case id."""
    test, bracket, params = nodeid.partition('[')
    if not bracket:
        return test, ''
    # Pad, to match whole '-'-separated parts of the case id
    return test, f'-{params[:-1]}-'


def _has_case(params, part):
    """Return whether `part` is a part of the '-'-padded `params`."""
    return f'-{part}-' in params


def _is_listed(item, cases):
    """Return whether `item` is one of the `cases` of SLOW_CASES."""
    test, params = _split_case_id(item.nodeid)
    parts = cases.get(test, ())
    return any(_has_case(params, part) for part in parts)


def _is_release_case(item, cases):
    """Return whether `item` is a case of RELEASE_CASES."""
    _, params = _split_case_id(item.nodeid)
    return any(
        _has_case(params, part)
        for fixture, parts in cases.items()
        if fixture in item.fixturenames
        for part in parts
        )


def pytest_itemcollected(item):
    """Remember the node ids and fixtures of all items.

    This is done before any deselection.
    """
    collected = item.config.stash.setdefault(_COLLECTED, [])
    collected.append((item.nodeid, frozenset(item.fixturenames)))


def _find_stale_cases(collected_items, cases):
    """Return descriptions of `cases` entries matching no node id.

    Only entries of modules with at least one collected node id are
    checked, so that running a subset of the suite is not an error.

    Parameters
    ----------
    collected_items : Iterable
        Pairs of node id and fixture names of all collected items.
    cases : dict
        The table SLOW_CASES of tiers.

    Returns
    -------
    stale : list
        One string per entry (test) or case part that matches nothing.
    """
    collected = collections.defaultdict(list)  # test: case ids
    modules = set()
    for nodeid, _ in collected_items:
        test, params = _split_case_id(nodeid)
        modules.add(test.partition(_NODEID_SEP)[0])
        collected[test].append(params)
    stale = []
    for test, parts in cases.items():
        if test.partition(_NODEID_SEP)[0] not in modules:
            continue
        if test not in collected:
            stale.append(f'{test} (no such test)')
            continue
        stale.extend(
            f'{test}[{part}] (no such case)'
            for part in parts
            if not any(_has_case(p, part) for p in collected[test])
            )
    return stale


def _find_stale_release_cases(collected, cases):
    """Return descriptions of `cases` entries matching no case.

    Only fixtures that at least one collected item uses are checked,
    so that running a subset of the suite is not an error.

    Parameters
    ----------
    collected : Iterable
        Pairs of node id and fixture names of all collected items.
    cases : dict
        The table RELEASE_CASES of tiers.

    Returns
    -------
    stale : list
        One string per case part that matches no case of its fixture.
    """
    stale = []
    for fixture, parts in cases.items():
        params = [
            _split_case_id(nodeid)[1]
            for nodeid, fixtures in collected
            if fixture in fixtures
            ]
        stale.extend(
            f'{fixture}[{part}] (no such case)'
            for part in parts
            if params and not any(_has_case(p, part) for p in params)
            )
    return stale


# Name of the tables of tiers, what they list, and how to find stale
# entries in them
_TABLES = (
    (
        'SLOW_CASES',
        f'the cases that take {tiers.SLOW_THRESHOLD} s or more',
        _find_stale_cases,
        ),
    (
        'RELEASE_CASES',
        'the cases that only run with --release',
        _find_stale_release_cases,
        ),
    )


def _stale_cases_message(config):
    """Return an error message about stale tier entries, or None."""
    if any(_NODEID_SEP in arg for arg in config.args):
        return None  # Tests selected by node id: only part of a module
    if config.getoption('lf', default=False):
        return None  # --lf drops the passed items of a module early
    collected = config.stash.get(_COLLECTED, ())
    messages = []
    for name, what, find_stale in _TABLES:
        stale = find_stale(collected, getattr(tiers, name))
        if not stale:
            continue
        entries = '\n  '.join(stale)
        messages.append(
            f'Stale entries in {name} ({tiers.__file__}), which lists '
            f'{what}: the tests or cases below do not exist (renamed '
            f'or removed?). Update the module.\n  {entries}'
            )
    return '\n'.join(messages) or None


def _check_stale_cases(config, items):
    """Report stale tiers.SLOW_CASES/RELEASE_CASES as a usage error.

    A pytest-xdist worker cannot raise it: that ends the controller
    with an internal error. Instead, a worker drops all its items
    (so that no test runs) and hands the message over to the
    controller, which reports it and ends with the same exit status
    as a serial run (pytest_testnodedown, pytest_sessionfinish).
    """
    message = _stale_cases_message(config)
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
    """End a pytest-xdist run with stale tier entries as an error."""
    message = session.config.stash.get(_STALE_MESSAGE, None)
    if message:
        sys.stderr.write(f'ERROR: {message}\n')
        session.exitstatus = pytest.ExitCode.USAGE_ERROR


def pytest_collection_modifyitems(config, items):
    """Mark TensErLEED, slow, and release tests; check the tables."""
    _check_stale_cases(config, items)
    for item in items:
        if tiers.TENSORLEED_FIXTURES.intersection(item.fixturenames):
            item.add_marker(pytest.mark.tensorleed)
        if _is_listed(item, tiers.SLOW_CASES):
            item.add_marker(pytest.mark.slow)
        if _is_release_case(item, tiers.RELEASE_CASES):
            item.add_marker(pytest.mark.release)
