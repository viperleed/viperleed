"""
===============================================
      ViPErLEED Graphical User Interface
===============================================
 *** module guilib.leedsim.dialogbulk3dsym ***

Created: 2021-05-17
Author: Michele Riva

"""

import numpy as np

import PyQt5.QtCore as qtc
import PyQt5.QtGui as qtg
import PyQt5.QtWidgets as qtw

from viperleed import guilib as gl


TEXT_FOR_OP = {gl.PlaneGroup.C6: '6-fold screw/rotations',
               gl.PlaneGroup.C4: '4-fold screw/rotations',
               gl.PlaneGroup.C3: '3-fold screw/rotations',
               gl.PlaneGroup.C2: '2-fold screw/rotation',
               gl.PlaneGroup.Mx: 'mirror/glide across [1, 0]',
               gl.PlaneGroup.My: 'mirror/glide across [0, 1]',
               gl.PlaneGroup.M11: 'mirror/glide across [1, 1]',    # == M45
               gl.PlaneGroup.M1m1: 'mirror/glide across [1, -1]',  # == Mm45
               gl.PlaneGroup.M10: 'mirror/glide across [1, 0]',
               gl.PlaneGroup.M01: 'mirror/glide across [0, 1]',
               gl.PlaneGroup.M21: 'mirror/glide across [2, 1]',
               gl.PlaneGroup.M12: 'mirror/glide across [1, 2]'}


class Bulk3DSymDialog(qtw.QDialog):
    """Dialog to input extra symmetry operations for bulk.


    Attributes
    ----------
    None.

    Methods
    -------
    update_operations(bulk)
        Update widgets to include all admissible operations
        that may be selected by the user. Should be called
        before opening/showing/executing the dialog.
    """

    def __init__(self, parent=None):
        super().__init__(parent)

        self.setWindowModality(qtc.Qt.WindowModal)
        flags = self.windowFlags()
        flags &= ~qtc.Qt.WindowCloseButtonHint #disable close button
        self.setWindowFlags(flags)
        self.setWindowTitle('Add bulk symmetry operations')

        # sub-widgets of which we need to keep a reference
        self.__all_widgets = {'conserve_sym': qtw.QCheckBox(),
                              'extra_ops': qtw.QListView(),
                              'buttons': {'done': qtw.QPushButton(),
                                          'cancel': qtw.QPushButton()}}

        # The following are set update_operations():
        self.__bulk = None      # gl.Lattice
        self.__extra_ops = []   # extra operations that may be selected

        self.__compose()
        self.__connect()

    def __compose(self):
        """Place widgets in the right spots."""
        description = qtw.QLabel()
        description.setText("Pick extra symmetry operations for the bulk. "
                            "\nThese are useful to account for symmetries "
                            "that are not captured by the plane group, e.g., "
                            "screw axes perpendicular to the surface, and "
                            "glide operations where mirroring is across a "
                            "plane perpendicular to the surface, and "
                            "translation has a pure z component.\nThey are "
                            "useful in case the sample has steps, and the "
                            "reconstruction rotates (or is mirrored) from one "
                            "terrace to the next [e.g., Fe3O4(001) or Si(001)]."
                            "\n")
        description.setWordWrap(True)
        description.setFont(gl.AllGUIFonts().labelFont)

        widgets = [description]

        sym_conserve = self.__all_widgets['conserve_sym']
        sym_conserve.setText("Constrain symmetry to bulk cell")
        sym_conserve.setFont(gl.AllGUIFonts().labelFont)
        sym_conserve.setChecked(True)
        widgets.append(sym_conserve)

        # Until we figure out how to handle the input (see Issue #15)
        # leave it checked and disabled
        sym_conserve.setEnabled(False)

        # Actual list of available extra symmetry operations
        widgets.append(self.__all_widgets['extra_ops'])

        # Done and Cancel buttons
        buttons = self.__all_widgets['buttons']
        buttons['done'].setText('&Done')
        buttons['cancel'].setText('&Cancel')
        widgets.extend([but for but in buttons.values()])

        # Make sure every graphical aspect of widgets is up to date
        for widget in widgets:
            widget.ensurePolished()

        # Change resizing policy for some widgets:
        # - Done and Cancel buttons should not resize
        for but in buttons.values():
            but.setSizePolicy(qtw.QSizePolicy.Fixed, qtw.QSizePolicy.Fixed)

        # Place widgets into the layout
        layout = qtw.QVBoxLayout()
        for widget in widgets:
            if widget in buttons.values():
                continue
            layout.addWidget(widget)

        buttons_layout = qtw.QHBoxLayout()
        for button in buttons.values():
            buttons_layout.addWidget(button)
            buttons_layout.setAlignment(qtc.Qt.AlignRight)
        layout.addLayout(buttons_layout)

        self.setLayout(layout)

    def __connect(self):
        """Hook up signals."""
        self.__all_widgets['extra_ops'].clicked.connect(self.__item_clicked)

        buttons = self.__all_widgets['buttons']
        buttons['cancel'].clicked.connect(self.close)
        buttons['done'].clicked.connect(self.__done_pressed)

    def __get_symmetry_related_operations(self, operation):
        """Return operations symmetry-related to the given one.

        'Related' means it is an operation that results from
        applying (operation*another_operation) where
        another_operation is one of those selected, or one
        of the 2D plane group operations of the bulk. The list
        will not contain operation itself.

        Parameters
        ----------
        operation : tuple
            2x2 matrix of a symmetry operation.

        Returns
        -------
        list
        """
        # This is very hard to do properly, when allowing the user
        # to input arbitrary operations as one would have to check
        # all the combinations of operations.
        # TODO: do this later when we figure out how to handle this.

        #       For now, look through the gl.PlaneGroup(s) for the
        #       current bulk.cell_shape, and pick the smallest group
        #       that contains all the selected operations (+ those of
        #       the plane group itself). Then use the operations of
        #       this group to automatically select the ones that are
        #       symmetry-related.
        extras = self.__all_widgets['extra_ops'].model().item
        all_ops = (*self.__bulk.group.operations(),
                   *[op
                     for i, op in enumerate(self.__extra_ops)
                     if extras(i).checkState()])

        # Get the smallest group that contains all the operations
        for group in gl.PlaneGroup.groupsForShape[self.__bulk.cell_shape]:
            group_ops = gl.PlaneGroup(group).operations()
            if all(op in group_ops for op in all_ops):
                break
            group_ops = []

        # For now stick to the easy ones: C6 --> C3 & C2; C4 -> C2
        related = []
        for related_op in group_ops:
            if related_op == operation:
                continue  # skip the one passed
            if related_op not in self.__extra_ops:
                continue  # skip those not in list
            related.append(related_op)

        return related

        # if operation in (gl.PlaneGroup.C6, gl.PlaneGroup.C4):
            # for related_op in (gl.PlaneGroup.C2, gl.PlaneGroup.C3):
                # try:
                    # related.append(self.__extra_ops.index(related_op))
                # except ValueError:
                    # pass

    def __done_pressed(self):
        """Return with code == qtw.QDialog.Accepted."""
        super().done(qtw.QDialog.Accepted)
        # self.accept()
        # self.close()

    def __item_clicked(self, index):
        """Toggle check-box of clicked item and those related.

        This is a slot to be connected to QListView.clicked.

        Parameters
        ----------
        index : QtCore.QModelIndex
        """
        item = self.__all_widgets['extra_ops'].model().itemFromIndex(index)

        if not item.isEnabled():
            return

        # toggle check state
        if item.checkState() == qtc.Qt.Unchecked:
            checked = qtc.Qt.Checked
        else:
            checked = qtc.Qt.Unchecked
        item.setCheckState(checked)

        # If the user checks an operation, others related to it via
        # either the bulk operations or others that are already selected
        # must also get selected automatically (and disabled).
        related_ops = self.__get_symmetry_related_operations(
            self.__extra_ops[index.row()]
            )

        items = self.__all_widgets['extra_ops'].model().item
        for op in related_ops:
            idx = self.__extra_ops.index(op)
            items(idx).setCheckState(checked)
            items(idx).setEnabled(not bool(checked))

    @property
    def n_extra_ops(self):
        """Return the number of extra operations that may be selected."""
        model = self.__all_widgets['extra_ops'].model()
        if model.item(0):
            return model.rowCount()
        return 0

    def update_operations(self, bulk):
        """Place admissible extra operations in the extra_ops widget.

        This method should be called before .show(), .open(), or
        .exec() to have the list of symmetry operations up to date.

        Parameters
        ----------
        bulk : gl.Lattice

        Returns
        -------
        None.
        """
        self.__bulk = bulk
        bulk_ops = self.__bulk.group.operations()

        # Get the plane group with most operations given the
        # current cell shape (always the last in groupsForShape)
        largest_group = gl.PlaneGroup(
            bulk.group.groupsForShape[bulk.cell_shape][-1]
            )
        self.__extra_ops = []
        for op in largest_group.operations():
            # Skip those that are already in the bulk group
            if op in bulk_ops:
                continue
            # Skip those that are the inverse of one that is
            # already present. They will always be added, but
            # we don't want the user to bother selecting both.
            # They are: Cm6 (for C6), Cm3 (for C3) and Cm4 (for C4)
            if op in (gl.PlaneGroup.Cm6, gl.PlaneGroup.Cm3, gl.PlaneGroup.Cm4):
                continue
            self.__extra_ops.append(op)

        # Add also those screws_glides that we may not have
        for op in bulk.group.screws_glides:
            if op in self.__extra_ops:
                continue
            self.__extra_ops.append(op)

        n_extra = len(self.__extra_ops)

        model = qtg.QStandardItemModel(n_extra, 1)  # rows, col
        for r, extra_op in enumerate(self.__extra_ops):
            txt = TEXT_FOR_OP.get(extra_op, None)

            # Handle arbitrary input
            if txt is None:
                # Not one of the known operations.
                # For now place the matrix straight as is.
                # TODO: adapt once we figure out how to handle this
                txt = str(extra_op)

            item = qtg.QStandardItem(txt)
            item.setFlags(qtc.Qt.ItemIsEnabled)

            if op in bulk.group.screws_glides:
                checked = qtc.Qt.Checked
            else:
                checked = qtc.Qt.Unchecked
            item.setData(checked, qtc.Qt.CheckStateRole)
            model.setItem(r, 0, item)

        self.__all_widgets['extra_ops'].setModel(model)

    def extra_operations(self):
        """Return the selected operations.

        This method should be called after the dialog
        returns with a qtw.QDialog.Accepted exit code.

        Returns
        -------
        str
            A string in the form expected by the property setter
            gl.PlaneGroup.screws_glides, i.e.,
                "r(#, #, ...), m([#, #], [#, #], ...)"
            NB: the setter should be called with a tuple having
            this as its first element, and bulk.cell_shape as
            the second element.
        """
        items = self.__all_widgets['extra_ops'].model().item

        glides = []
        screws = []

        for i in range(self.n_extra_ops):
            if items(i).checkState() == qtc.Qt.Unchecked:
                continue
            data = items(i).data(qtc.Qt.DisplayRole)
            if 'mirror' in data:       # take direction
                glides.append(data.split('mirror/glide across ')[1])
            elif 'rotation' in data:   # take order
                screws.append(data.split('-fold')[0])
            else:
                raise NotImplementedError(
                    f"Unknown operation {data}. Handling of arbitrary "
                    "symmetry operations is not yet sorted (Issue #15)."
                    )

        txt = []
        if screws:
            txt.append(f"r({', '.join(screws)})")
        if glides:
            txt.append(f"m({', '.join(glides)})")
        return ', '.join(txt)





