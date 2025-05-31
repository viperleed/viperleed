"""Module widgetdecorators of viperleed.gui.

This module provides decorators that can be used to decorate QWidget
subclasses.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2020-01-29'
__license__ = 'GPLv3+'

import wrapt
import PyQt5.QtCore as qtc
import PyQt5.QtGui as qtg
import PyQt5.QtWidgets as qtw

from viperleed.gui.widgetslib import get_all_children_widgets
from viperleed.gui.decorators import ensure_decorates_class

mouse_event_types = (
    qtc.QEvent.MouseButtonPress,
    qtc.QEvent.MouseButtonRelease,
    qtc.QEvent.MouseButtonDblClick,
    qtc.QEvent.MouseMove,
    qtc.QEvent.Wheel
    )

@ensure_decorates_class(qtw.QWidget)
@wrapt.decorator
def receive_mouse_broadcast(cls, _, args, kwargs):
    """
    Decorator for QWidget classes used for declaring that the class wants to
    receive mouse events even if they are not originally targeted to the
    instance. Instances of the decorated class will receive mouse events
    triggered on the ancestor or on any descendant (child, grandchild,
    grand-grandchild, etc.) of the ancestor, provided the ancestor does
    @broadcast_mouse events.
    """
    setattr(cls, '_wants_broadcast', True)

    return cls(*args, **kwargs)

@ensure_decorates_class(qtw.QWidget)
@wrapt.decorator
def broadcast_mouse(cls, _, args, kwargs):
    """
    This is a class decorator intended to be used on a top-level window class.
    The decorated class will broadcast any mouse event to all children (as well
    as grand-children, grand-grand-children, etc..) widgets that are instances
    of a QWidget subclass decorated with @receive_mouse_broadcast.
    (A specific mouse event will not be sent twice to the target of the event.)

    After broadcasting is complete, the mouse events will be handled as
    implemented in the class definition.
    """

    # Keep references to a few original method names that will be redefined,
    # so that they can be called at appropriate places later on
    __orig_init = cls.__init__
    __orig_childEvent = cls.childEvent
    __orig_eventFilter = cls.eventFilter

    def __init__(self, *args, **kwargs):
        self.__broadcast_mouse_to = set()
        self.__last_broadcasted = dict()
        self.__polished_children = set()
        # All events generated within the instance will pass through its
        # eventFilter method defined below
        __orig_init(self, *args, **kwargs)

        self.installEventFilter(self)

    def childEvent(self, event):
        """
        Every time a child of type QWidget is added or removed, install/remove
        the eventFilter of this class on the child and all grand-children
        that want to @receive_mouse_broadcast.

        The original implementation of childEvent() of the decorated class
        is also called at the end.
        """
        # Notice that event.added() is fired at the very instant a widget is
        # added, but event.child() is not necessarily the right pointer to the
        # child as the widget might not be polished yet
        if isinstance(event.child(), qtw.QWidget):
            children = get_all_children_widgets(event.child())
            if (event.polished()
                    and event.child() not in self.__polished_children):
                self.__polished_children.update(children)
                self.__broadcast_mouse_to.update(child
                                                 for child in children
                                                 if hasattr(child,
                                                            '_wants_broadcast'))
                for child in children:
                    child.installEventFilter(self)
                    child.destroyed.connect(self.on_child_destroyed)
            elif event.removed() and event.child() in self.__polished_children:
                # Notice that event.child() will be in self.__polished_children
                # unless it is destroyed before the childEvent is processed
                # In that case the on_child_destroyed method takes care of
                # releasing the appropriate references
                self.__polished_children -= set(children)
                self.__broadcast_mouse_to -= set(children)
                for child in children:
                    child.removeEventFilter(self)
        # Create also a dictionary of "child: time-stamp of last event sent out"
        # that allows to not send multiple times the same event to the same
        # child, and prevents infinite recursion in case more than one child
        # wants to @receive_mouse_broadcast. Notice that the use of the sets
        # give a unique list of children here
        self.__last_broadcasted = {child: 0
                                   for child
                                   in self.__broadcast_mouse_to}
        # then call whatever implementation was there on the undecorated class
        # This should be a 'pass' except for QMainWindow and QSplitter
        __orig_childEvent(self, event)

    def on_child_destroyed(self, child):
        """
        This is the slot associated with the destroyed() signal of children
        widgets for which references are kept. The method releases these
        references so that the widgets cannot be accessed any longer after
        being destroyed.
        """
        self.__polished_children.discard(child)
        self.__broadcast_mouse_to.discard(child)

    def eventFilter(self, target, event):
        """
        Reimplementation of eventFilter that broadcasts all mouse events to
        children widgets that want to @receive_mouse_broadcast.

        The original implementation of the decorated class is also always
        called.
        """
        # make sure that all child events that need to be addressed by the
        # childEvent method are correctly handled. This part is needed for
        # QMainWindow: if the centralWidget was already set (and consequently
        # already handled by childEvent above) when a new widget is added
        # as its child, no child event is emitted without the next lines
        if event.type() in (qtc.QEvent.ChildPolished, qtc.QEvent.ChildRemoved):
            self.childEvent(event)

        if event.type() in mouse_event_types:
            for child, timestamp in self.__last_broadcasted.items():
                if child != target and timestamp != event.timestamp():
                    # Notice the filtering is on the basis of event.timestamp()
                    # that has a 1 ms resolution: This should be good enough
                    # even for mouseMoveEvent, as the maximum polling
                    # rate of mice is 1kHz to date (Jan. 2020)

                    if not child.underMouse():
                        edited_event = self.__edit_mouse_event(child, event)
                    else:
                        edited_event = event
                    if not self.__skip_event_broadcast(child, edited_event):
                        qtw.qApp.sendEvent(child, edited_event)
                self.__last_broadcasted[child] = event.timestamp()
        return __orig_eventFilter(self, target, event)

    def __edit_mouse_event(self, target, event):
        """
        Process the event by changing its location, so that widgets that are not
        under the mouse do not get animated (e.g., subclasses of QAbstractButton
        like QPushButton animate parts of the widget)
        """
        top_left = target.rect().topLeft()
        new_local_pos = top_left - qtc.QPoint(5, 5)
        window_pos = target.mapTo(target.window(), new_local_pos)
        screen_pos = target.mapToGlobal(new_local_pos)
        new_source = qtc.Qt.MouseEventSynthesizedByApplication
        if event.type() == qtc.QEvent.Wheel:
            new_event = qtg.QWheelEvent(new_local_pos, screen_pos,
                                        event.pixelDelta(), event.angleDelta(),
                                        event.buttons(), event.modifiers(),
                                        event.phase(), event.inverted(),
                                        new_source)
        else:
            new_event = qtg.QMouseEvent(event.type(), new_local_pos, window_pos,
                                        screen_pos, event.button(),
                                        event.buttons(), event.modifiers(),
                                        new_source)
        new_event.setTimestamp(event.timestamp())
        new_event.setAccepted(event.isAccepted())
        return new_event

    def __skip_event_broadcast(self, target, event):
        """
        This function is used to treat specific events to specific widgets in
        a different way.

        Now implemented:
           * QLineEdit subclasses that do not haveFocus() should not get
             mouseMove nor mouseDoubleClick events as the Qt implementation does
             not check whether the widget hasFocus() nor whether it is
             underMouse() in the event handlers, so that QLineEdit respond
             anyway (in both cases a portion of the text gets highlighted as if
             it was selected)

        For future extension to other cases that I have not thought about,
        return a logic OR of the conditions
        """
        if (not target.isEnabled()
            or target.testAttribute(qtc.Qt.WA_TransparentForMouseEvents)):
            return True

        skip_line_edit = (isinstance(target, qtw.QLineEdit)
                          and not target.hasFocus()
                          and event.type() in (qtc.QEvent.MouseButtonDblClick,
                                               qtc.QEvent.MouseMove))
        return skip_line_edit


    setattr(cls, '__init__', __init__)
    setattr(cls, 'childEvent', childEvent)
    setattr(cls, 'eventFilter', eventFilter)
    setattr(cls, 'on_child_destroyed', on_child_destroyed)
    setattr(cls, '__edit_mouse_event', __edit_mouse_event)
    setattr(cls, '__skip_event_broadcast', __skip_event_broadcast)

    return cls(*args, **kwargs)
