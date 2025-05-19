"""Module basewidgets of viperleed.guilib.widgets.

========================================
   ViPErLEED Graphical User Interface
========================================

Created: 2024-06-13
Author: Michele Riva
Author: Florian Dörr

This module defines basic, non-specific widgets.
"""

from PyQt5 import QtCore as qtc
from PyQt5 import QtWidgets as qtw

from viperleed.guilib.widgetslib import remove_spacing_and_margins


_ALIGN_CTR = qtc.Qt.AlignHCenter
_PIXEL_SPACING = 4


class ButtonWithLabel(qtw.QWidget):
    """QLabel and QPushButton."""

    def __init__(self, tight=True, **kwargs):
        """Initialise widget.

        Parameters
        ----------
        tight : bool, optional
            Whether the layout should not have
            spacing around it. Default is true.

        Returns
        -------
        None
        """
        super().__init__(**kwargs)
        self.label = qtw.QLabel()
        self.button = QNoDefaultPushButton()
        self._compose(tight)

    def _compose(self, tight):
        """Compose widget."""
        layout = qtw.QHBoxLayout()
        layout.addWidget(self.label)
        layout.addWidget(self.button)
        if tight:
            layout.setContentsMargins(0, 0, 0, 0)
        self.setLayout(layout)

    @qtc.pyqtSlot(str)
    def set_label_text(self, text):
        """Set the text of the label."""
        self.label.setText(text)

    @qtc.pyqtSlot(str)
    def set_button_text(self, text):
        """Set the text of the button."""
        self.button.setText(text)


class CollapsibleList(qtw.QScrollArea):
    """Base class for CollapsibleLists."""

    def __init__(self, parent=None):
        """Initialise widget.

        Parameters
        ----------
        parent : QObject
            The parent QObject of this widget.

        Returns
        -------
        None.
        """
        super().__init__(parent=parent)
        # The _views dict contains the displayed CollapsibleViews as
        # keys and the corresponding values are lists of the widgets
        # that are displayed right next to each CollapsibleView in the
        # CollapsibleList.
        self._views = {}
        self._widths = {}
        self._top_widget_types = []
        self._layout = None
        self.setWidgetResizable(True)
        self.clear()
        self.setFrameStyle(self.Panel | self.Sunken)

    @property
    def views(self):
        """Return views."""
        return self._views

    def clear(self):
        """Clear all views in the list and delete references to them."""
        self._views = {}
        self._layout = qtw.QVBoxLayout()
        self._layout.setSpacing(0)
        self._make_scroll_area()

    def insert_view(self, view):
        """Insert new view at the bottom of the list.

        Parameters
        ----------
        view : CollapsibleDeviceView
            The view that is to be inserted into the
            CollapsibleList.

        Returns
        -------
        None.
        """
        self._add_top_widgets_to_view(view)
        self._layout.insertWidget(self._layout.count()-1, view)

    def _add_top_widgets_to_view(self, view):
        """Add the top widget types to the CollapsibleView.

        Parameters
        ----------
        view : CollapsibleView
            The CollapsibleView to which the widgets will be attached.

        Returns
        -------
        None.
        """
        self.views[view] = []
        for widget_type in self._top_widget_types:
            widget = widget_type()
            view.add_top_widget(widget)
            self.views[view].append(widget)

    def _add_top_widget_types(self, *widg_types):
        """Extend the list of widgets to add.

        Must be called in the '__init__' of subclasses.

        Parameters
        ----------
        *widg_types : type(QWidget)
            A widget type of which an instance should be added
            next to the button of each CollapsibleView.

        Returns
        -------
        None.
        """
        self._top_widget_types.extend(widg_types)

    def _make_scroll_area(self):
        """Compose QScrollArea."""
        widget = qtw.QWidget(parent=self)
        self._layout.addStretch(1)
        widget.setLayout(self._layout)
        self.setWidget(widget)


class CollapsibleView(qtw.QWidget):
    """A widget that can be expanded/collapsed by the press of a button."""

    def __init__(self, parent=None):
        """Initialise widget.

        parent : QObject
            The parent QObject of this widget.

        Returns
        -------
        None.
        """
        super().__init__(parent=parent)
        self._button = QNoDefaultIconButton()
        self._frame = qtw.QFrame(parent=self)
        self._bottom_space = qtw.QSpacerItem(1, 1)
        self._compose_and_connect()

    @property
    def button(self):
        """Return the main button."""
        return self._button

    def add_collapsible_item(self, item):
        """Add widget to the widgets in the inner collapsible layout.

        Parameters
        ----------
        item : QWidget or QLayout
            Any widget or layout that should be
            added to the collapsible frame.

        Raises
        -------
        TypeError
            If the given item to add to the layout
            is neither a layout nor a widget.
        """
        if isinstance(item, qtw.QWidget):
            self._frame.layout().addWidget(item)
        elif isinstance(item, qtw.QLayout):
            self._frame.layout().addLayout(item)
        else:
            raise TypeError('Cannot add items to the collapsible view '
                            'that are neither a layout nor a widget.')

    def add_top_widget(self, widget, width=None, align=_ALIGN_CTR):
        """Add widget to the widgets in the outer top layout.

        Parameters
        ----------
        widget : QWidget
            Any widget that should be added next
            to the constantly displayed button.
        width : int or None, optional
            The pixel width the widget should be displayed with.
            Default is None, which means the widget itself
            determines the pixel width.
        align : Qt.AlignmentFlag, optional
            How the widget should be aligned horizontally.
            Default is central alignment.

        Returns
        -------
        None.
        """
        layout = qtw.QVBoxLayout()
        height = (self.button.sizeHint().height() -
                  widget.sizeHint().height())//2
        if height > 0:
            layout.addSpacing(height)
        # We use a surrounding widget in order to be
        # able to set a sizePolicy later on.
        surrounding_widget = qtw.QWidget(parent=self)
        widget_layout = qtw.QHBoxLayout()
        remove_spacing_and_margins(widget_layout)
        widget_layout.addWidget(widget)
        layout.addWidget(surrounding_widget)
        surrounding_widget.setLayout(widget_layout)
        layout.addStretch(1)
        self.layout().addLayout(layout)
        self.layout().addSpacing(_PIXEL_SPACING)
        self.set_top_widget_geometry(widget, width=width, align=align)

    def is_enabled(self):
        """Return whether the view is enabled."""
        return self.button.isEnabled()

    @qtc.pyqtSlot(int)
    @qtc.pyqtSlot(bool)
    def set_expanded_state(self, expanded):
        """Expand/collapse the collapsible part of this widget.

        Parameters
        ----------
        expanded : bool or int
            Converted to bool. Decides whether the collapsible
            widgets are visible/expanded. True means they become/stay
            visible.

        Returns
        -------
        None.
        """
        expanded = bool(expanded)
        self.button.setEnabled(expanded)
        self._toggle_expanded_state(expanded)

    def set_top_widget_geometry(self, widget, width=None, align=_ALIGN_CTR):
        """Adjust width and alignment of a widget in the outer top layout.

        Parameters
        ----------
        widget : QWidget
            The widget in the outer top layout whose geometry should be
            changed. (One of the widgets that have been added via
            add_top_widget.)
        width : int or None, optional
            The pixel width the widget should be displayed with.
            Default is None, which means the widget itself
            determines the pixel width.
        align : Qt.AlignmentFlag, optional
            How the widget should be aligned horizontally.
            Default is central alignment.

        Returns
        -------
        None.
        """
        for lay_i in range(1, self.layout().count()):
            layout = self.layout().itemAt(lay_i)
            try:
                lay_widget = layout.itemAt(1).widget()
            except AttributeError:
                # This means the item in question is a QSpacerItem.
                continue
            inner_layout = lay_widget.layout()
            inner_widget = inner_layout.itemAt(0).widget()
            if inner_widget is not widget:
                continue
            policy = lay_widget.sizePolicy()
            inner_layout.setAlignment(widget, align)
            if width:
                lay_widget.setMinimumWidth(width)
                lay_widget.setMaximumWidth(width)
            policy.setHorizontalPolicy(policy.Fixed)
            lay_widget.setSizePolicy(policy)
            return

    def _adjust_bottom_space(self):
        """Update spacing below the QFrame depending on frame visibility."""
        spacer_height = (1 if not self._frame.isVisible()
                         else self._button.sizeHint().height()/2)
        self._bottom_space.changeSize(1, round(spacer_height))
        self._bottom_space.invalidate()

    def _adjust_button_icon(self, button_up):
        """Change in which direction the icon in the top button is pointing.

        Parameters
        ----------
        button_up : bool
            True if the button is meant to point upwards.

        Returns
        -------
        None.
        """
        icon = (qtw.QStyle.SP_TitleBarShadeButton if button_up
                else qtw.QStyle.SP_TitleBarUnshadeButton)
        self.button.setIcon(self.style().standardIcon(icon))

    def _compose_and_connect(self):
        """Assemble children widgets and connect their signals."""
        # layout is the layout which contains the
        # button and the collapsible frame.
        layout = qtw.QVBoxLayout()
        remove_spacing_and_margins(layout)
        layout.addWidget(self.button)
        policy = self.button.sizePolicy()
        policy.setHorizontalPolicy(policy.Expanding)
        self.button.setSizePolicy(policy)
        self.button.setIcon(
            self.style().standardIcon(qtw.QStyle.SP_TitleBarUnshadeButton)
            )

        # inner_layout is the layout that is inside of the frame.
        # It is necessary to add widgets to the frame.
        inner_layout = qtw.QVBoxLayout()
        remove_spacing_and_margins(inner_layout)
        self._frame.setFrameStyle(self._frame.StyledPanel | self._frame.Plain)
        self._frame.setLayout(inner_layout)
        self._frame.setVisible(False)

        # frame_layout is the layout in which the frame is nested in.
        # It is necessary for spacing.
        frame_layout = qtw.QHBoxLayout()
        frame_layout.addSpacing(1)
        frame_layout.addWidget(self._frame)
        frame_layout.addSpacing(1)
        layout.addLayout(frame_layout)

        layout.addItem(self._bottom_space)
        # outer_layout is the main layout which holds button,
        # frame, and the top widgets.
        outer_layout = qtw.QHBoxLayout()
        remove_spacing_and_margins(outer_layout)
        outer_layout.addLayout(layout)
        self.setLayout(outer_layout)

        self.button.clicked.connect(self._on_button_clicked)

    @qtc.pyqtSlot()
    def _on_button_clicked(self):
        """Expand/collapse the CollapsibleView."""
        self._toggle_expanded_state(not self._frame.isVisible())

    def _toggle_expanded_state(self, expanded):
        """Expand/collapse the CollapsibleView."""
        self._frame.setVisible(expanded)
        self._adjust_bottom_space()
        self._adjust_button_icon(expanded)


class QCheckBoxInvertedSignal(qtw.QCheckBox):
    """QCheckBox with extra unchecked signal."""

    unchecked = qtc.pyqtSignal(bool)

    def __init__(self, **kwargs):
        """Initialise widget."""
        super().__init__(**kwargs)
        self.stateChanged.connect(self._emit_inverted_signal)

    @qtc.pyqtSlot(int)
    def _emit_inverted_signal(self, value):
        """Emit unchecked signal."""
        self.unchecked.emit(not bool(value))


class QNoDefaultDialogButtonBox(qtw.QDialogButtonBox):
    """QDialogButtonBox without default buttons."""

    def event(self, event):
        """Overwrite event to skip setting default."""
        if event.type() == qtc.QEvent.Show:
            self._unset_default_buttons()
            return qtw.QWidget().event(event)
        return super().event(event)

    def _unset_default_buttons(self):
        """Ensure no button is a default."""
        for button in self.buttons():
            button.setAutoDefault(False)
            button.setDefault(False)


class QNoDefaultPushButton(qtw.QPushButton):
    """QPushbutton that is not a default button."""

    def __init__(self, *args, **kwargs):
        """Initialise button."""
        super().__init__(*args, **kwargs)
        self.setAutoDefault(False)
        self.setDefault(False)


class QNoDefaultIconButton(QNoDefaultPushButton):
    """QNoDefaultPushButton with right-aligned icons."""

    def __init__(self, *args, **kwargs):
        """Initialise button."""
        super().__init__(*args, **kwargs)
        self.setLayout(qtw.QHBoxLayout())
        self._label = qtw.QLabel()
        self._label.setAlignment(qtc.Qt.AlignRight | qtc.Qt.AlignVCenter)
        self._label.setAttribute(qtc.Qt.WA_TransparentForMouseEvents, True)
        self.layout().addWidget(self._label)

    def setIcon(self, icon):        # pylint: disable=invalid-name
        """Set the desired icon on the button.

        Parameters
        ----------
        icon : QIcon
            The icon to be added to the button.
        """
        icon_size_real = icon.actualSize(self.sizeHint()*0.5)
        self._label.setPixmap(icon.pixmap(icon_size_real))


class QUncheckableButtonGroup(qtw.QButtonGroup):
    """QButtonGroup that can be unchecked."""

    def __init__(self, *args, **kwargs):
        """Initialise button group."""
        super().__init__(*args, **kwargs)
        self._pressed_button_was_checked = False
        self.buttonPressed.connect(self._remember_checked_state_of_button)

    def uncheck_buttons(self):
        """Uncheck all buttons."""
        self.setExclusive(False)
        self.checkedButton().setChecked(False)
        self.setExclusive(True)

    @qtc.pyqtSlot(qtw.QAbstractButton)
    def uncheck_if_clicked(self, button):
        """Uncheck if checked button was clicked."""
        was_checked = self._pressed_button_was_checked
        self._pressed_button_was_checked = False
        if was_checked and button == self.checkedButton():
            self.uncheck_buttons()

    @qtc.pyqtSlot(qtw.QAbstractButton)
    def _remember_checked_state_of_button(self, button):
        if button == self.checkedButton():
            self._pressed_button_was_checked = True
