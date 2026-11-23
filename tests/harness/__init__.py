"""Machinery of the ViPErLEED test suite.

The modules of this package define the command-line options of the
test suite, the tiers of tests and their markers, the comparison of
written files to stored reference files, and the handling of missing
optional dependencies. Most of them are pytest plugins, loaded via
`pytest_plugins` in tests/conftest.py, which also calls the hooks of
optional_dependencies. Their own tests are in the tests subpackage.

Modules
-------
check_slow
    The --check-slow option: list the tests whose 'slow' marking
    may need to change.
markexpr
    Helper for the plugins that exclude a marker from the selection.
optional_dependencies
    Skip tests on missing optional dependencies, or fail them with
    --require-all (except those exempted with --allow-missing).
quick
    The --quick option (quick tier): deselect the 'slow' tests.
reference_files
    The check_reference fixture and the --rewrite-references option.
release
    The --release option (release tier): also run the 'release'
    tests.
rng_seed
    The --rng-seed option for random-input tests.
tier_marks
    Mark the slow, release, and TensErLEED tests, and report stale
    entries of the tables in tiers.
tiers
    Data: which tests belong to the slow, release, and TensErLEED
    tiers.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-11-23'
__license__ = 'GPLv3+'
