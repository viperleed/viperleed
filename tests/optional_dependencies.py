"""Skip tests on missing optional dependencies, or fail on request.

Tests that need an optional dependency (mpifort/gfortran) skip, with
a reason naming the dependency, when it is missing. With the
command-line option --require-all, those skips become failures, so
that a run on a complete environment cannot silently skip tests.

Call the `missing_dependency` fixture where a test finds that a
dependency is missing:

    def test_something(missing_dependency):
        if not shutil.which('gfortran'):
            missing_dependency('gfortran')

The option and the fixture below are used by tests/conftest.py.
"""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2026-10-21'
__license__ = 'GPLv3+'

import pytest

REQUIRE_ALL_OPTION = '--require-all'


def _skip_or_fail(config, name, detail=''):
    """Skip the current test, or fail if --require-all was given."""
    reason = f'optional dependency {name} not available'
    if detail:
        reason = f'{reason}: {detail}'
    if config.getoption('require_all'):
        pytest.fail(f'{reason} (required by {REQUIRE_ALL_OPTION})',
                    pytrace=False)
    pytest.skip(reason)


def pytest_addoption(parser):
    """Add the option that makes missing dependencies fail tests."""
    parser.addoption(
        REQUIRE_ALL_OPTION,
        dest='require_all',
        action='store_true',
        default=False,
        help=('Fail, rather than skip, tests that need a missing '
              'optional dependency (mpifort/gfortran).'),
        )


@pytest.fixture(name='missing_dependency')
def fixture_missing_dependency(request):
    """Return a function that skips (or fails) the test for `name`."""
    def _missing(name, detail=''):
        _skip_or_fail(request.config, name, detail)
    return _missing
