"""Tests for module dialogbulk3dsym of viperleed.gui.leedsim.dialogs."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-11-18'
__license__ = 'GPLv3+'

from pytest_cases import fixture

from viperleed.gui.classes.lattice2d import Lattice2D

_C6 = '6-fold screw/rotations'
_RELATED_TO_C6 = {'3-fold screw/rotations', '2-fold screw/rotation'}
_HEXAGONAL = [[1, 0], [-0.5, 3**0.5 / 2]]


class TestBulk3DSymDialog:
    """Tests for the Bulk3DSymDialog class."""

    @fixture(name='dialog')
    def fixture_dialog(self, import_with_qt, show):
        """Return a shown dialog for a hexagonal bulk with p1 group."""
        module = import_with_qt('viperleed.gui.leedsim.dialogs.'
                                'dialogbulk3dsym')
        dialog = module.Bulk3DSymDialog()
        dialog.update_operations(Lattice2D(_HEXAGONAL, group='p1'))
        show(dialog)
        return dialog

    @fixture(name='extra_ops')
    def fixture_extra_ops(self, dialog):
        """Return the view listing the extra operations of `dialog`."""
        # pylint: disable-next=protected-access
        return dialog._Bulk3DSymDialog__all_widgets['extra_ops']

    @staticmethod
    def get_states(extra_ops, qtc):
        """Return {text: (checked, enabled)} for each operation."""
        model = extra_ops.model()
        checked = qtc.Qt.CheckState.Checked
        return {item.text(): (item.checkState() == checked,
                              item.isEnabled())
                for item in (model.item(i) for i in range(model.rowCount()))}

    def test_toggle_related(self, extra_ops, qtc):
        """Check that related operations follow a toggled one.

        Checking the 6-fold rotation checks and disables the 3- and
        2-fold rotations, but no other operation. Unchecking it again
        unchecks and enables them.
        """
        states = self.get_states(extra_ops, qtc)
        assert all(s == (False, True) for s in states.values())
        texts = list(states)
        c6_index = extra_ops.model().index(texts.index(_C6), 0)

        extra_ops.clicked.emit(c6_index)
        states = self.get_states(extra_ops, qtc)
        for text, state in states.items():
            expect = ((True, False) if text in _RELATED_TO_C6
                      else (text == _C6, True))
            assert state == expect, text

        extra_ops.clicked.emit(c6_index)
        states = self.get_states(extra_ops, qtc)
        assert all(s == (False, True) for s in states.values())

    def test_update_operations_checks_bulk(self, dialog, extra_ops, qtc):
        """Check that only the screws/glides of the bulk are checked."""
        bulk = Lattice2D(_HEXAGONAL, group='p1')
        bulk.group.set_screws_glides('r(3), m([1, 0])', bulk.cell_shape)
        dialog.update_operations(bulk)
        states = self.get_states(extra_ops, qtc)
        checked = {text for text, (is_checked, _) in states.items()
                   if is_checked}
        assert checked == {'3-fold screw/rotations',
                           'mirror/glide across [1, 0]'}
