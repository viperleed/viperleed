"""Skip tests on missing optional dependencies, or fail on request.

Tests that need an optional dependency (one of KNOWN_DEPENDENCIES,
the registry of their names) skip, with a reason naming the
dependency, when it is missing. With the command-line option
--require-all, those skips become failures, so that a run on a
complete environment cannot silently skip tests.
The option --allow-missing exempts the dependencies it names (comma
separated, may be repeated) from --require-all: their tests skip as
without it, e.g., on a machine with no compilers, with
--require-all --allow-missing=mpifort,gfortran. A name that is not
one of KNOWN_DEPENDENCIES is a usage error.

Call the `missing_dependency` fixture where a test finds that a
dependency is missing:

    def test_something(missing_dependency):
        if not shutil.which('gfortran'):
            missing_dependency('gfortran')

The options, their check and the fixture below are used by
tests/conftest.py.
"""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-10-21'
__license__ = 'GPLv3+'

import pytest

ALLOW_MISSING_OPTION = '--allow-missing'
KNOWN_DEPENDENCIES = ('gfortran', 'mpifort')  # Add new names here
REQUIRE_ALL_OPTION = '--require-all'


def _skip_or_fail(config, name, detail=''):
    """Skip the current test, or fail if --require-all was given."""
    reason = f'optional dependency {name} not available'
    if detail:
        reason = f'{reason}: {detail}'
    if config.getoption('require_all'):
        if name not in _allowed(config):
            pytest.fail(f'{reason} (required by {REQUIRE_ALL_OPTION})',
                        pytrace=False)
        reason = f'{reason} (allowed by {ALLOW_MISSING_OPTION})'
    pytest.skip(reason)


def _allowed(config):
    """Return the names given to --allow-missing."""
    given = config.getoption('allow_missing') or ()
    return {name.strip()
            for names in given
            for name in names.split(',')
            if name.strip()}


def pytest_addoption(parser):
    """Add the options that make missing dependencies fail tests."""
    known = '/'.join(KNOWN_DEPENDENCIES)
    parser.addoption(
        REQUIRE_ALL_OPTION,
        dest='require_all',
        action='store_true',
        default=False,
        help=('Fail, rather than skip, tests that need a missing '
              f'optional dependency ({known}).'),
        )
    parser.addoption(
        ALLOW_MISSING_OPTION,
        dest='allow_missing',
        action='append',
        metavar='NAME[,NAME...]',
        help=('Skip, even with --require-all, tests that need one of '
              'these missing optional dependencies (comma separated, '
              'may be repeated). Example: '
              '--allow-missing=mpifort,gfortran.'),
        )


def pytest_configure(config):
    """Reject unknown dependency names given to --allow-missing."""
    unknown = _allowed(config).difference(KNOWN_DEPENDENCIES)
    if not unknown:
        return
    unknown_names = ', '.join(sorted(unknown))
    valid = ', '.join(KNOWN_DEPENDENCIES)
    msg = (f'{ALLOW_MISSING_OPTION}: unknown optional dependencies '
           f'{unknown_names}. Valid names are {valid}.')
    raise pytest.UsageError(msg)


@pytest.fixture(name='missing_dependency')
def fixture_missing_dependency(request):
    """Return a function that skips (or fails) the test for `name`.

    The test fails if `name` is not one of KNOWN_DEPENDENCIES.
    """
    def _missing(name, detail=''):
        if name not in KNOWN_DEPENDENCIES:
            pytest.fail(f'{name!r} is not one of the KNOWN_DEPENDENCIES '
                        f'{KNOWN_DEPENDENCIES} of {__name__}',
                        pytrace=False)
        _skip_or_fail(request.config, name, detail)
    return _missing
