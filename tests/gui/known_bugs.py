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

import sys

import pytest

_PY_310 = sys.version_info >= (3, 10)

# HoverAnnot.sizeHint passes numpy floats to QSize, whose constructor
# takes ints. Python >= 3.10 no longer converts floats implicitly. The
# TypeError is raised in a Qt virtual method: pytest-qt fails the test.
FLOATS_TO_QSIZE = pytest.mark.xfail(
    _PY_310,
    reason='HoverAnnot.sizeHint passes floats to QSize (Python >= 3.10)',
    )

# MPLFigureCanvas.mpl_to_qt_coordinates uses _dpi_ratio of matplotlib
# canvases, which does not exist any more in matplotlib >= 3.5 (the
# oldest supported).
MPL_TO_QT_COORDINATES = pytest.mark.xfail(
    raises=AttributeError,
    reason='MPLFigureCanvas.mpl_to_qt_coordinates uses _dpi_ratio',
    )

# The keyPressEvent of the old NewFileDialog checks for keys in
# [Key_Enter or Key_Return], i.e., only for Key_Enter. Return clicks
# the default button of the dialog rather than the focused one.
RETURN_KEY_OLD_NEW_FILE_DIALOG = pytest.mark.xfail(
    reason='Return does not click the focused button of the old '
           'NewFileDialog',
    )

# TextBox() without arguments uses `contents` before assignment (a
# typo in its __init__).
TEXTBOX_WITHOUT_ARGUMENTS = pytest.mark.xfail(
    raises=UnboundLocalError,
    reason='TextBox() without arguments uses contents before assignment',
    )
