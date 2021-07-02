"""Module modulebase of viperleed.guilib.

======================================
  ViPErLEED Graphical User Interface
======================================

Defines the ViPErLEEDModuleBase class from which all ViPErLEED
modules should inherit. Any concrete implementation of a module
should use super().closeEvent() if it wants to accept a closeEvent()
rather than just .accept()ing the event. The AboutViPErLEED class is
also defined, a common 'About' dialog that is present on all modules
that are subclasses of ViPErLEEDModuleBase.

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


class ViPErLEEDModuleBase(qtw.QMainWindow):
    """Base class for a ViPErLEED module."""

    module_closed = qtc.pyqtSignal(object)  # The class being destroyed

    def __init__(self, parent=None, name=''):
        """Initialize module."""
        super().__init__(parent)
        self.name = name
        self.about_action = qtw.QAction("About")
        self.__about = AboutViPErLEED(module_name=name,
                                      module_description='')

        self.setAttribute(qtc.Qt.WA_DeleteOnClose)
        self.setWindowIcon(qtg.QIcon(LOGO))

        self.__compose()

    def __compose(self):
        """Prepare a basic menu bar that all modules will have."""
        menu = self.menuBar()
        menu.addAction(self.about_action)

        self.about_action.triggered.connect(self.__about.show)

    def closeEvent(self, event):  # pylint: disable=invalid-name
        """Reimplement closeEvent to emit a module_closed."""
        self.module_closed.emit(self)
        super().closeEvent(event)

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
        print(layout.spacing())
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
            'and it is released under the GNU General Public License '
            '<a href="https://www.gnu.org/licenses/gpl-3.0">version 3</a>'
            ' or later.<p>ViPErLEED is developed as a '
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
