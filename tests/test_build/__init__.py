"""Tests for the contents of the build folder of the repository.

The tests otherwise follow the layout of the code they test (e.g.,
tests/calc for src/viperleed/calc), so this package would be called
tests/build. It is not, because pytest skips folders named build (as
well as dist, venv, and others; see the norecursedirs option of
pytest) when looking for tests: its tests would then run only when the
folder is given explicitly on the command line, never with
`pytest tests`.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-11-23'
__license__ = 'GPLv3+'
