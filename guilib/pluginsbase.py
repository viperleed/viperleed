"""Module modulebase of viperleed.guilib.

======================================
  ViPErLEED Graphical User Interface
======================================

Defines the ViPErLEEDPluginBase class from which all ViPErLEED
modules should inherit. Any concrete implementation of a module
should use super().closeEvent() if it wants to accept a closeEvent()
rather than just .accept()ing the event. The AboutViPErLEED class is
also defined, a common 'About' dialog that is present on all modules
that are subclasses of ViPErLEEDPluginBase.

Author: Michele Riva
Created: 2021-06-29
"""

from PyQt5 import (QtWidgets as qtw,
                   QtCore as qtc,
                   QtGui as qtg)

from viperleed.gui import resources_path
from viperleed import guilib as gl


LOGO = resources_path('guilib/icons/viperleed_logo_circled_48x48.png')


def logo_one_line():
    """Return a QLabel with the one-liner ViPErLEED logo."""
    logo_pixmap = qtg.QPixmap(resources_path(
            'guilib/icons/viperleed_logo_oneline_200x38.png'
            ))
    logo = qtw.QLabel()
    logo.setPixmap(logo_pixmap)
    return logo


# TODO: find a nice way to set up a default menu with
# icons that are common to all ViPErLEEDPluginBase windows.
# The main difficulty is having instance-dependent tool-tips.

# Probably best to rather load settings from a settings.ini
# file (or similarly named) that must be contained in the
# __path__ of the __module__ that defines the plug-in class.

class ViPErLEEDPluginBase(qtw.QMainWindow):
    """Base class for the main window of a ViPErLEED plug-in."""

    module_closed = qtc.pyqtSignal(object)  # The class being destroyed
    
    # screen_changed tracks changes of screen. This signal carries
    # the old and the new screen objects. This  signal is emitted
    # only after the window is dropped onto the new screen. It can
    # be connected in subclasses, e.g., to adapt the size of the
    # window when the screen is changed.
    screen_changed = qtc.pyqtSignal(qtg.QScreen, qtg.QScreen)

    def __init__(self, parent=None, name='', description=''):
        """Initialize module.

        Parameters
        ----------
        parent : QWidget, optional
            The parent widget of this plug-in window. Typically
            not given as plug-ins are usually independent, top
            level windows. Default is None.
        name : str. optional
            Name of the plug-in. Currently unused. Default is
            an empty string.
        description : str, optional
            Verbose description of plug-in. Currently unused.
            Default is an empty string
        """
        super().__init__(parent)
        self.name = name
        self.about_action = qtw.QAction("About")
        self.__about = AboutViPErLEED(module_name=name,
                                      module_description=description)

        self.setAttribute(qtc.Qt.WA_DeleteOnClose)
        self.setWindowIcon(qtg.QIcon(LOGO))

        # Native windows are a little slower to draw, but offer
        # the advantage of giving access to windowHandler() (a
        # QWindow), which has a screenChanged signal.
        self.setAttribute(qtc.Qt.WA_NativeWindow)

        # __screens keeps track of the old
        # and new screens of the plug-in
        self.__screens = {'old': self.screen(),
                          'new': self.screen()}

        self.__compose()
        self.__connect()
        self.installEventFilter(self)

    @property
    def old_screen(self):
        """Return the last screen of the plug-in."""
        return self.__screens['old']
    
    @old_screen.setter
    def old_screen(self, old_screen):
        """Set the last screen of the plug-in."""
        self.__screens['old'] = old_screen

    @property
    def new_screen(self):
        """Return the new screen of the plug-in."""
        return self.__screens['new']

    @new_screen.setter
    def new_screen(self, new_screen):
        """Set the new screen of the plug-in.

        This is also the slot connected to the screenChanged
        signal of the native QWindow holding this QMainWindow

        Parameters
        ----------
        new_screen : QScreen
            The new screen object

        Emits
        -----
        screen_changed
            If the screen actually changed
        """
        self.__screens['new'] = new_screen
        if self.new_screen != self.old_screen:
            self.screen_changed.emit(self.old_screen, self.new_screen)
            self.old_screen = self.new_screen

    # TODO: check that it works nicely with multiple
    # monitors and different DPIs & scaling!
    def center_on_screen(self):
        """Center the window in the current screen."""
        window_rect = self.frameGeometry()
        screen_center = qtw.qApp.desktop().availableGeometry(self).center()
        window_rect.moveCenter(screen_center)
        self.move(window_rect.topLeft())

    def eventFilter(self,                # pylint: disable=invalid-name
                    watched_object, event):
        """Extend eventFilter to filter events.

        The following events are currently processed:
            atc.QEvent.NonClientAreaMouseButtonRelease
                Triggered when the window has been moved and the
                mouse button is release. The window is resized
                to a dimension that fits the current screen. The
                event is then processed normally.

        Returns
        -------
        reject_event : bool
            True if the event should not be processed.
        """
        if event.type() == qtc.QEvent.NonClientAreaMouseButtonPress:
            self.old_screen = self.screen()
        elif event.type() == qtc.QEvent.NonClientAreaMouseButtonRelease:
            self.new_screen = self.screen()
        return super().eventFilter(watched_object, event)

    def closeEvent(self, event):       # pylint: disable=invalid-name
        """Reimplement closeEvent to emit a module_closed."""
        self.module_closed.emit(self)
        super().closeEvent(event)

    def keyPressEvent(self, event):    # pylint: disable=invalid-name
        """Extend keyPressEvent to close 'About' on 'Esc'."""
        if event.key() == qtc.Qt.Key_Escape:
            self.__about.hide()
        super().keyPressEvent(event)

    def mousePressEvent(self, event):  # pylint: disable=invalid-name
        """Reimplement QMainWindow.mousePressEvent.

        The reimplementation closes the 'About' dialog if visible
        whenever a mouse click occurs on the window.

        This works so-so, as events may be filtered by children
        and not propagated. This is the case especially for
        disabled children widgets.

        Parameters
        ----------
        event : QMouseEvent
            The mouse-press event itself.

        Returns
        -------
        None.
        """
        if self.__about.isVisible():
            self.__about.close()
        super().mousePressEvent(event)

    def __compose(self):
        """Prepare a basic menu bar that all modules will have."""
        menu = self.menuBar()
        menu.addAction(self.about_action)

    def __connect(self):
        """Connect relevant signals to slots."""
        self.about_action.triggered.connect(self.__about.show)



class AboutViPErLEED(qtw.QWidget):
    """Window showing information about ViPErLEED."""

    def __init__(self, parent=None, module_name='', module_description=''):
        """Initialize module."""
        self.module_name = module_name
        self.module_description = module_description
        self.__close_btn = qtw.QPushButton('Close')

        super().__init__(parent)
        self.setWindowFlags(qtc.Qt.Tool
                            | qtc.Qt.WindowTitleHint
                            | qtc.Qt.CustomizeWindowHint
                            | qtc.Qt.WindowTitleHint)
        self.setWindowTitle('About ViPErLEED')
        # self.setWindowFlags(qtc.Qt.Window | qtc.Qt.FramelessWindowHint)

        #                    | (qtc.Qt.CustomizeWindowHint
        #                       & ~qtc.Qt.WindowMinMaxButtonsHint
        #                       & ~qtc.Qt.WindowCloseButtonHint)

        self.__compose()
        self.__connect()

    def __compose(self):
        """Set up children widgets."""
        layout = qtw.QVBoxLayout()
        layout.setSpacing(20)
        max_width = 600
        papers = (
            'X. Y, X. Y, and X. Y, ViPErLEED I, <I>Journal</I> <B>Vol</B>,'
            'pages (2021).',
            'X. Y, X. Y, and X. Y, ViPErLEED II, <I>Journal</I> <B>Vol</B>,'
            'pages (2021).',
            'X. Y, X. Y, and X. Y, ViPErLEED III, <I>Journal</I> <B>Vol</B>,'
            'pages (2021).',
            )
        contrib = ('Michele Riva', 'Florian Kraushofer', 'Michael Schmid',
                   'Lutz Hammer', 'Tilman Ki\u00dflinger', 'Florian D\u00f6rr',
                   'Bernhard Mayr',)

        txt = qtw.QLabel(
            'Copyright (2019\u22122021) ViPErLEED Team<p>'
            '<p>ViPErLEED (i.e., the Vienna Package for Erlangen LEED) is an '
            'open-source project that aims at making LEED-I(V) accessible to '
            'the broad scientific community for solving surface structures.<p>'
            'Should you find ViPErLEED useful, we are happy to receive '
            'citation to our papers:<p>'
            f'{"<br>".join(papers)}<p>'
            'ViPErLEED includes tools for measuring and analyzing LEED-I(V) '
            'data. It can calculate theoretical I(V) curves from a structural '
            'model, and optimize this model to best match experimental data. '
            'You can find more information by visiting our documentation '
            'page, or using the <U>H</U>elp.<p>'
            'The most recent code is available on our GitHub repository '
            '<a href="https://github.com/viperleed">github.com/viperleed</a>, '
            'and is released under the GNU General Public License '
            '<a href="https://www.gnu.org/licenses/gpl-3.0">version 3</a>'
            ' or later.<p>Bugs can be reported using the GitHub '
            '<a href="https://github.com/viperleed/viperleed/issues">Issues'
            '</a> or via email (<a href="mailto:riva@iap.tuwien.ac.at>'
            'riva@iap.tuwien.ac.at</a>).<p>ViPErLEED is developed as a '
            'collaboration between the Surface Physics group at the '
            '<a href="https://www.iap.tuwien.ac.at/www/surface/index">'
            'Institute of Applied Physics</a> of the TU Wien, and the '
            '<a href="https://www.fkp.physik.nat.fau.eu/">Chair of Solid '
            'State Physics</a> of the FAU Erlangen-N\u0252rnberg.<p>'
            f'Contributors: {", ".join(contrib)}'
            )
        txt_font = gl.AllGUIFonts().labelFont
        txt_font.setPointSize(9)
        txt.setFont(txt_font)
        txt.setWordWrap(True)
        txt.setMaximumWidth(max_width)
        txt.setOpenExternalLinks(True)
        txt.ensurePolished()

        layout.addWidget(logo_one_line(), alignment=qtc.Qt.AlignHCenter)
        layout.addWidget(txt)
        layout.addWidget(self.__close_btn)

        self.setLayout(layout)
        self.setMaximumWidth(max_width)
        palette = self.palette()
        palette.setColor(palette.Background, qtc.Qt.white)
        self.setAutoFillBackground(True)
        self.setPalette(palette)
        self.adjustSize()
        self.setFixedSize(self.size())

    def __connect(self):
        """Connect signals."""
        self.__close_btn.clicked.connect(self.close)
