"""Fixtures for the tests of viperleed.gui.leedsim."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-11-24'
__license__ = 'GPLv3+'

from pytest_cases import fixture


@fixture(name='check_annotations_follow_move')
def fixture_check_annotations_follow_move(qtc, qtg):
    """Return a function that checks that annotations follow a window.

    The function places the annotations at non-integer positions, then
    calls the moveEvent of the window.
    """
    def _check(window, annotations):
        for annotation in annotations:
            annotation.center = (10.4, 20.6)
        old_pos = window.pos()
        new_pos = old_pos + qtc.QPoint(5, 7)
        window.moveEvent(qtg.QMoveEvent(new_pos, old_pos))
        assert all(tuple(a.center) == (15, 28) for a in annotations)
    return _check
