"""Known bugs of viperleed.gui that make GUI tests fail.

The tests that fail because of a bug that is not fixed yet are marked
as expected failures with the marks defined here, rather than with
inline xfail marks, so that these bugs are listed in one place. When
fixing a bug, remove its mark, here and where it is used: strict
marks make the tests fail once the bug is fixed (XPASS).

Some bugs show up only with some versions of Python or of other
dependencies: their marks have a condition.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-11-23'
__license__ = 'GPLv3+'

import pytest
