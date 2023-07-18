"""
=======================================
   ViPErLEED Graphical User Interface
=======================================
 *** module guilib.leedsim.widgets ***

Created: 2020-01-12
Author: Michele Riva

Blah blah TODO
"""

# TESTING
from time import perf_counter
# TESTING

import numpy as np
import matplotlib as mpl
#OR
#from matplotlib.patches import Circle as MplCircle
from matplotlib import colors as mpl_colors
import PyQt5.QtCore as qtc
import PyQt5.QtGui as qtg
import PyQt5.QtWidgets as qtw

# import guilib as gl
from viperleed import guilib as gl

TEST = False

class RealCanvas(gl.MPLFigureCanvas):
    
    def plotLattices(self):
        self.ax.cla()
        
        _win = self.window()  # this should be the LEED_GUI instance
        # if not isinstance(_win, gl.LEED_GUI):  # REMOVED 2020-01-25
            # raise RuntimeError("I expected a LEED_GUI instance, "
                               # f"but I found {_win}")
        
        rs = _win.real
        
        rotation = _win.get_angle()
        
        #### BULK #### 
        rotBulkB = rs.bulk.get_rotated_basis(rotation)
        
        # Plot bulk as lines. A line is p+x*v, where p is any point along 
        # the line, and v a vector parallel to the line
        alpha = 0.15
        lineColor = tuple(alpha*fg + (1-alpha)*bg
                          for (fg, bg) in zip(mpl_colors.to_rgb('gray'),
                                              mpl_colors.to_rgb('white')))
        
        x = np.array([-rs.fov, rs.fov])
        K = max(rs.bulk.hk.ravel())
        # the following is a fast code to plot lines, adapted from 
        # http://exnumerus.blogspot.com/2011/02/
        #                              /how-to-quickly-plot-multiple-line.html
        xlist = []
        ylist = []
        for i in [0, 1]:
            slope = [x*v for v in rotBulkB[i]]
            for k in range(-K, K):
                xlist.extend(slope[0] + k*rotBulkB[(i+1)%2, 0])
                xlist.append(None)
                ylist.extend(slope[1] + k*rotBulkB[(i+1)%2, 1])
                ylist.append(None)
        self.ax.plot(xlist, ylist, color=lineColor, zorder=-1)
        for unitVec in rotBulkB:
            self.ax.annotate("", xy=tuple(unitVec),
                        xytext=(0, 0),
                        arrowprops=dict(arrowstyle="->", color='gray')
                        )
        
        #### SURFACE ####
        windowScaleFactor = (self.rightSize().height()/434)**2
                            # this is used for scaling the size of the spots
        rotSurfB = rs.surf.get_rotated_basis(rotation)
        rotSurfL = rs.surf.get_rotated_lattice(rotation)
        
        # Plot surface as scatter
        self.ax.scatter(rotSurfL[:,0], rotSurfL[:,1],
                        s=8*windowScaleFactor, c='k', zorder=2)
        for unitVec in rotSurfB:
            self.ax.annotate("", xy=tuple(unitVec),
                             xytext=(0, 0),
                             arrowprops=dict(arrowstyle="->"))
        
        self.setAxLimits(rs.fov)
        self.ax.figure.canvas.draw_idle()


class LEEDCanvas(gl.MPLFigureCanvas):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        
        # Get rid of the spines completely in the LEED pattern
        self.setSpinesOn(False)
        
        # Define annotations that will appear when the mouse hovers above
        # LEED spots (including symmetry-equivalent). Notice that at most 12
        # symmetry-equivalent spots can exist at the same time
        self.annots = [gl.HoverAnnot(self.parentWidget()) for _ in range(12)]
        
        if TEST:
            self.mpl_connect('motion_notify_event', self.on_hover)
        
        # ## TESTING
        # self.mpl_connect('figure_enter_event', self.enter_figure)
        # self.mpl_connect('figure_leave_event', self.leave_figure)
        # self.mpl_connect('button_press_event', self.on_mousepress)
        
        tmpoffs = 0
        deltas = [np.array([np.cos(x + tmpoffs), np.sin(x + tmpoffs)])
                  for x in np.linspace(0, 2*np.pi,
                                       num=len(self.annots), endpoint=False)]
        for (annot, delta) in zip(self.annots, deltas):
            annot.head_pos = delta*10 + annot.center
        # ## END TESTING
        
    def initLEEDScreen(self, lineThick=0.1):
        # initialize LEED screen to a radius of 1
        if not isinstance(lineThick, (int, float)):
            raise ValueError("initLEEDScreen requires a non-negative"
                             " floating point input")
        elif lineThick < 0:
            raise ValueError("line thickness cannot be negative")
        else:
            screen = mpl.patches.Circle((0, 0), 1,  facecolor='white',
                                        edgecolor='k', lw=lineThick,
                                        zorder=-1)

        self.ax.add_patch(screen)
        self.setAxLimits(1)
        self.screen = screen
        
    def resetLEEDCanvas(self):
        self.ax.cla()  # clear the axes
        self.setSpinesOn(False)  # turn off again the spines
        self.ax.add_patch(self.screen)  # and add once more the screen circle
        
        _win = self.window()

        energy = _win.get_energy()  # read energy from input
        self.rotation = _win.get_angle()
        
        self.leed = _win.leed
        
        screenR = self.leed.screen_radius(energy)
        self.screen.set_radius(screenR)
        self.setAxLimits(screenR)
        
        windowScaleFactor = (self.rightSize().height()/434)**2
        
        # scaling factor for the markers, also depends on the window size
        self.markScale = np.sqrt(self.leed.max_energy/energy)*windowScaleFactor
        
        # Big dot at (0,0), on top of everything
        self.ax.scatter(0, 0, s=100*self.markScale, c='k', zorder=1e10)
        
        # plot the bulk spots as hollow circles
        bulkLEED = self.leed.rotate(self.rotation, 'bulk')
        self.bulkLEED = self.ax.scatter(bulkLEED[:,0], bulkLEED[:,1],
                                        s=25*self.markScale,
                                        facecolors='none',
                                        edgecolors='k')
        
        # Now clip the pattern with the circle
        self.bulkLEED.set_clip_on(True)
        self.bulkLEED.set_clip_path(self.screen)
    
    def plotLEED(self, hideDoms):
        self.resetLEEDCanvas()
        self.leed.rotate(self.rotation, 'surf')
        
        if hideDoms or self.leed.n_domains == 1:
            toPlot = self.leed.firstLEED
        else:
            toPlot = self.leed.domsLEED
        
        self.leedPat = []  # list of scatter plots
        for patt in toPlot:
            self.leedPat.append(
                            self.ax.scatter(patt.rotBeams[:,0],
                                            patt.rotBeams[:,1],
                                            s=5*self.markScale*patt.sizeScale,
                                            c=patt.color, marker=patt.marker
                                            )
                               )
        for patt in self.leedPat:
            patt.set_clip_on(True)
            patt.set_clip_path(self.screen)
            # Set pick radius for hovering annotations (scale with marker size).
            # Whenever the mouse is closer to a spot than pickradius points, the
            # mouse event will be considered 'contained' in the pattern to which
            # the spot belongs
            patt.set_pickradius(3*np.sqrt(self.markScale))
        self.ax.figure.canvas.draw_idle()
    
    def showEvent(self, event):
        super().showEvent(event)
        
        # Now set the center of the LEED screen as the 'center' of annotations
        # (i.e., annotations will always be radially symmetric around the screen
        # center)
        # NB: this needs to be UPDATED when I can handle off-normal incidence,
        #     as the center should then be the (0, 0) spot!
        #     It's most likely better to move/replicate this part in the
        #     resetLEEDCanvas as it will need to know where the (0, 0) spot is
        local_center = self.geometry().center()
        global_center = self.parentWidget().mapToGlobal(local_center)
        for annot in self.annots:
            annot.center = global_center
    
    def on_hover(self, event):
        if not self.is_on_screen(event):
            [annot.hide() for annot in self.annots]
            return
        show_annots = False
        for patt in self.leedPat:
            contained, ind = patt.contains(event)
            if contained:
                show_annots = True
                mpl_pos = self.closest_point(
                    self.scatter_coords(patt, ind['ind']),
                    (event.x, event.y))
                qt_pos = self.mpl_to_qt_coordinates(*mpl_pos)
                
                # TESTING
                self.annots[0].head_pos = self.mapToGlobal(qt_pos)
                # end TESTING
        
        for annot in self.annots:  # needs to select only those that have to be
                                   # shown
            annot.update_position()
            annot.setVisible(show_annots)
            annot.update()
    
    def is_on_screen(self, event):
        """
        Checks if the mouse event occurred inside the LEED screen.
        """
        # in principle it looks like one could use the matplotlib patch
        # .contains method to do this, but for reasons I don't understand it
        # seems to return True as soon as the mouse enters the bounding box of
        # the screen rather than its interior
        if event.inaxes != self.ax:
            return False
        event_r = np.sqrt(event.xdata**2 + event.ydata**2)
        return event_r <= self.screen.radius
    
    def scatter_coords(self, scatt_plot, indices):
        """
        Given a scatter plot and a list of indices it returns the (x, y)
        coordinates of the scatter points at the corresponding indices. The
        coordinates are given in the reference system of the canvas.
        """
        to_restore = scatt_plot.get_offset_position()
        scatt_plot.set_offset_position('data')
        mpl_pos = self.ax.transData.transform(
                    scatt_plot.get_offsets()[indices])
        scatt_plot.set_offset_position(to_restore)
        
        return mpl_pos
    
    def closest_point(self, points, close_to):  # move to base
        """
        Out of the array-like 'points' returns the one that is closest to
        'close_to'
        """
        dist = np.linalg.norm(np.subtract(points, close_to), axis=1)
        closest = np.where(abs(dist - dist.min()) <= 1e-8)
        
        return np.asarray(points)[closest][0]

    # TESTING
    if TEST:
        def enter_figure(self, event):
            for annot in self.annots:
                annot.show()
        
        def leave_figure(self, event):
            for annot in self.annots:
                annot.hide()
        
        def on_mousepress(self, event):
            x_qt, y_qt = self.mpl_to_qt_coordinates(event.x, event.y)
            print(f'(x, y) = MPL: ({event.x}, {event.y}), Qt: ({x_qt}, {y_qt})')
    # end TESTING
    



class RotationBlock(gl.TextBoxWithButtons):
    def __init__(self, parent=None):
        params = {
            'labelText': '&Rotation (\u00b0)',
            'textBoxText': '\u2014',
            'botButText': '\u21bb', #CW
            'topButText': '\u21ba', #CCW
            'parent': parent,
            'textBoxTip': 'Change rotation of lattices and pattern. '\
                          'Positive angles are counterclockwise',
            'topButTip': 'Rotate lattices and pattern 10° clockwise'\
                         'Hold Ctrl down for finer control',
            'botButTip': 'Rotate lattices and pattern 10° counterclockwise'\
                         'Hold Ctrl down for finer control',
            'textBoxWidth': 85
        }
        
        super().__init__(**params)
        self.setTips(**params)
        self.text.setLimits(-180, 180, 'cyclic')
        self.text.setStep(10, 'add')
        
        #set aliases for easier identification
        self.ccw = self.topBut
        self.cw = self.botBut
        
        #and build the bottom part
        self.makeBottomWidget()
    
    def makeBottomWidget(self):
        # bottomWidget looks like this:
        #
        #    'Align bulk'
        #     [1 0] hor ver
        #     [0 1] hor ver
        #
        # will be composed with a QGridLayout, that is then assigned to self.bottomWidget
        
        bWLay = qtw.QGridLayout()
        bWLay.setSpacing(0)
        bWLay.setContentsMargins(0, 0, 0, 0)
        
        # prepare the sub-widgets
        rotBulkLab = qtw.QLabel('Align bulk:')
        oneZeroLab = qtw.QLabel("[1 0]")
        zeroOneLab = qtw.QLabel("[0 1]")
        self.h10 = qtw.QPushButton("Hor.")
        self.v10 = qtw.QPushButton("Ver.")
        self.h01 = qtw.QPushButton("Hor.")
        self.v01 = qtw.QPushButton("Ver.")
        horVerButs = ((self.h10, self.v10), (self.h01, self.v01))
        
        # add the relevant ones to the list of sub-widgets
        self.subWidgs.extend([self.h10, self.v10, self.h01, self.v01,
                             *np.ravel(horVerButs)])
        
        # set fonts
        [lab.setFont(gl.AllGUIFonts().smallTextFont)
         for lab in [rotBulkLab,oneZeroLab,zeroOneLab]]
        [but.setFont(gl.AllGUIFonts().smallButtonFont)
         for but in np.ravel(horVerButs)]
        
        # set size policies and adjust the sizes to account for the new font
        for lab in [rotBulkLab, oneZeroLab, zeroOneLab]: 
            lab.setSizePolicy(qtw.QSizePolicy.Fixed, qtw.QSizePolicy.Preferred)
            lab.adjustSize()
        for but in np.ravel(horVerButs):
            but.setSizePolicy(qtw.QSizePolicy.Preferred, qtw.QSizePolicy.Fixed)
            but.setMaximumHeight(qtw.QPushButton().sizeHint().height()*.84)
        
        # set tooltips
        tips = ['Rotate lattices and pattern to bring the [1 0] '
                    + 'bulk vector horizontal',
                'Rotate lattices and pattern to bring the [1 0] '
                    + 'bulk vector vertical',
                'Rotate lattices and pattern to bring the [0 1] '
                    + 'bulk vector horizontal',
                'Rotate lattices and pattern to bring the [0 1] '
                    + 'bulk vector vertical']
        for (but, tip) in zip(np.ravel(horVerButs), tips):
            but.setStatusTip(tip)
            but.setToolTip(tip)
        
        # now add the widgets to the layout
        bWLay.addWidget(rotBulkLab, 0, 0, 1, 3, qtc.Qt.AlignLeft)
        bWLay.addWidget(oneZeroLab, 1, 0, 1, 1,
                        qtc.Qt.AlignHCenter | qtc.Qt.AlignLeft)
        bWLay.addWidget(zeroOneLab, 2, 0, 1, 1,
                        qtc.Qt.AlignHCenter | qtc.Qt.AlignLeft)
        [bWLay.addWidget(horVerButs[i][j], i+1, j+1, 1, 1, qtc.Qt.AlignCenter)
         for i in range(2) for j in range(2)]
        
        # set the minimum heights of rows and columns
        bWLay.setColumnMinimumWidth(0, oneZeroLab.width()*1.2)
        
        # Add the layout to the bottomWidget
        self.bottomWidget.setLayout(bWLay)
        self.bottomWidget.layout().activate()
    
    def on_botWidgPressed(self, pressed):
        direction = 0
        vert = False
        if pressed in (self.h01, self.v01):
            direction = 1
        if pressed in (self.v10, self.v01):
            vert = True
        angle = self.window().real.angle_for_horizontal_bulk(direction)
        if vert:
            angle += 90
        if qtw.qApp.keyboardModifiers() == qtc.Qt.ControlModifier:
            angle += 180
        self.text.updateText(angle)


class EnergyBlock(gl.TextBoxWithButtons):
    def __init__(self, parent=None):
        params = {'labelText': '&Energy (eV)',
                  'textBoxText': '\u2014',
                  'topButText': '\u25b2',  # energy up
                  'botButText': '\u25bc',  # energy down
                  'parent': parent,
                  'textBoxTip': 'Change primary energy of the electron beam',
                  'topButTip': 'Increase energy by 20%. ' \
                               'Hold Ctrl down for a 5% increase',
                  'botButTip': 'Decrease energy by 20%. ' \
                               'Hold Ctrl down for a 5% decrease',
                  'textBoxWidth': 85
                  }

        super().__init__(**params)
        self.setTips(**params)
        self.text.setStep(1.2, 'scale')
        
        #set aliases for easier identification
        self.enUp = self.topBut
        self.enDown = self.botBut
        
        # for reasons unknown the down arrow is typeset weirdly
        # (larger than it should) -> use smaller size
        self.enDown.setFont(gl.AllGUIFonts().smallButtonFont)
        
        self.makeBottomWidget()
    
    def makeBottomWidget(self):
        # The bottomWidget in this case contains one QlineEdit with
        # two text lines that state the minimum and maximum energy
        #
        self.limits = qtw.QLabel('Min = 10 eV\nMax = \u2014')
        self.limits.setFont(gl.AllGUIFonts().smallTextFont)
        self.limits.adjustSize()
        
        bwLay = qtw.QHBoxLayout()
        bwLay.setSpacing(0)
        bwLay.setContentsMargins(0,0,0,0)
        bwLay.addWidget(self.limits, qtc.Qt.AlignLeft | qtc.Qt.AlignTop)
        
        self.bottomWidget.setLayout(bwLay)
        self.subWidgs.extend([self.limits])
        
        self.layout().setRowMinimumHeight(5, self.limits.height())


class DomsBlock(qtw.QWidget):
    toggleText = ['Hide &Domains', 'Show &Domains']
    # These are the two possible strings that appear on the toggle button.
    # Pressing alt+d fires the toggle button
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self.hide = False
        
        # The next line is useful only if I ever want to place this in a 
        # layout
        self.setSizePolicy(qtw.QSizePolicy.Minimum, qtw.QSizePolicy.Minimum)
        self.text = qtw.QLabel('\u2014 inequivalent domain(s)')
        self.toggle = ToggleButton(self.toggleText[self.hide])
        
        self.subWidgs = (self.text, self.toggle)
        
        self.setTips(text='Open or drag-drop a LEED pattern file to'
                          'determine how many domains!',
                     but='Toggle visibility of symmetry-equivalent domains')
        self.compose()
    
    def setTips(self, **kwargs):
        if 'text' in kwargs.keys():
            self.text.setStatusTip(kwargs['text'])
        if 'but' in kwargs.keys():
            tip = kwargs['but']
            self.toggle.setToolTip(tip)
            self.toggle.setStatusTip(tip)
    
    def compose(self):
        #set fonts
        self.text.setFont(gl.AllGUIFonts().labelFont)
        self.toggle.setFont(gl.AllGUIFonts().buttonFont)
        
        #set sizes
        self.toggle.setSizePolicy(qtw.QSizePolicy.Minimum,
                                  qtw.QSizePolicy.Minimum)
        self.toggle.adjustSize()
        
        self.text.setAlignment(qtc.Qt.AlignRight | qtc.Qt.AlignBaseline)
        self.text.setSizePolicy(qtw.QSizePolicy.Expanding,
                                qtw.QSizePolicy.Minimum)
        self.text.adjustSize()
        self.text.setMinimumSize(self.text.size())
        self.toggle.setMinimumSize(self.toggle.size())
        
        #----- Set up the layout -----#
        lay = qtw.QVBoxLayout()
        lay.setSpacing(2)
        lay.setContentsMargins(0, 0, 0, 0)
        
        lay.addStretch(1)  # the stretches keep things together
        lay.addWidget(self.text) 
        lay.addWidget(self.toggle)
        lay.addStretch(1)
        lay.setAlignment(self.text, self.text.alignment())
        lay.setAlignment(self.toggle, qtc.Qt.AlignRight | qtc.Qt.AlignVCenter)
        
        self.setLayout(lay)
        
        self.adjustSize()
    
    def width(self):
        self.ensurePolished()
        return max(sub.width() for sub in self.subWidgs)
    
    def height(self):
        self.ensurePolished()
        return (sum(sub.height() for sub in self.subWidgs) 
                + self.layout().spacing())
    
    def minimumSizeHint(self):
        return self.sizeHint()
    
    def sizeHint(self):
        return qtc.QSize(self.width(), self.height())
    
    def togglePressed(self):
        self.hide = not self.hide
        self.toggle.setText(self.toggleText[self.hide])
        self.matricesPopup.updateMatrices(self.hide)
    
    def updateText(self, txt):
        if not isinstance(txt, str):
            raise
        
        self.text.setText(txt)
        self.text.setMinimumSize(self.text.sizeHint())
        
        trOld = self.geometry().topRight()
        self.adjustSize()
        newPos = (self.geometry().topLeft()
                  - self.geometry().topRight()
                  + trOld)
        self.move(newPos)

    def initPopup(self):
        leed = self.window().leed
        self.matricesPopup = MatricesPopup(leed.superlattices,
                                           leed.domColors,
                                           parent=self.window())
    
    #----------------------- Re-implement mouseEvents -----------------------#
    def mouseReleaseEvent(self, event):
        if self.underMouse():
            child = self.childAt(event.pos())
            if child == self.text:
                # Position and show matricesPopup
                if not self.matricesPopup.shown:
                    # not shown before --> place it at the standard position
                    tR = self.mapToGlobal(self.text.geometry().topRight())
                    sizePopup = self.matricesPopup.frameSize()
                    newPos = qtc.QPoint(tR.x() - 0.7*sizePopup.width(),
                                    tR.y() - sizePopup.height() - 25)
                    self.matricesPopup.move(newPos)
                    self.matricesPopup.dragPosition = newPos
                    self.matricesPopup.shown = True
                if self.matricesPopup.isVisible():
                    self.matricesPopup.hide()
                else:
                    self.matricesPopup.show()
        super().mouseReleaseEvent(event)


class ToggleButton(qtw.QPushButton):
    def sizeHint(self):
        self.ensurePolished()
        refbutton = qtw.QPushButton(self.text())
        refbutton.setFont(self.font())
        return refbutton.sizeHint()*1.2
    
    def minimumSizeHint(self):
        return self.sizeHint()


class MatricesPopup(qtw.QWidget):
    def __init__(self, matrices, colors,
                 fs=gl.AllGUIFonts().mathFont.pointSize(),
                 parent=None):
        self.parent = parent
        super().__init__()
        
        if colors is not None:  # there's more than one domain
            self.matrices = [gl.PainterMatrix(matrix, color=color, fs=fs)
                             for (matrix, color) in zip(matrices, colors)]
            # add one matrix to be displayed in black when the user toggles 
            # the visibility of the domains
            self.firstDomMatrix = gl.PainterMatrix(matrices[0], fs=fs)
        elif len(matrices) == 1:
            self.matrices = [gl.PainterMatrix(matrices[0], fs=fs)]
        else:  # some input error: colors are missing
            raise
        
        self.shown = False  # used for choosing whether is should pop up at 
                            # the standard position or not
        self.initPopup()
    
    def initPopup(self):
        # The next line makes the widget not steal focus from the main 
        # window when it is shown
        self.setAttribute(qtc.Qt.WA_ShowWithoutActivating)
        self.setParent(self.parent)  # This should come before the next 
                                     # lines to have the window stay on top 
                                     # of the main, but not on top of other 
                                     # applications
        
        #-------- Edit the appearance of the window --------#
        # Qt.Tool: thinner frame
        # Qt.CustomizeWindowHint: Remove titlebar, icon, and 
        #                         maximize/minimize/close buttons
        # Qt.WindowStaysOnTopHint: always keep it on top - do not use. 
        #                          Rather parent it to the MainWindow so 
        #                          that it stays only on top of it (and not 
        #                          on top of other apps)
        # Qt.WindowDoesNotAcceptFocus: prevent switching focus to the popup
        self.setWindowFlags(qtc.Qt.Tool
                            | qtc.Qt.CustomizeWindowHint
                            | qtc.Qt.WindowDoesNotAcceptFocus)
        gl.editStyleSheet(self, 'background-color: white;')
        
        # The next attribute is used in the mousePressEvent and 
        # mouseMoveEvent for dragging the window around
        self.dragPosition = self.frameGeometry().topLeft()
        
        #--------- Prepare layout for the matrices ---------#
        nMatrices = len(self.matrices)  # this is 1, 2, 3, 4, 6, 8, or 12
                                        # TRUE?
        if nMatrices in range(1, 4):  # 1 -- 3
            nRows = 1
        elif nMatrices in range(4, 10):  # 4 -- 8
            nRows = 2
        else:
            nRows = 3
        nCols = int(nMatrices/nRows)
        
        lay = qtw.QGridLayout()
        [lay.addWidget(matrix, index//nCols, index % nCols, qtc.Qt.AlignCenter)
         for (index, matrix) in enumerate(self.matrices)]
        
        if hasattr(self, 'firstDomMatrix'):
            # More than one domain. Add anyway the matrix of the first domain
            # for changing its color when domains are toggled. Initially 
            # hide (toggling the button will show it).
            lay.addWidget(self.firstDomMatrix, 0, 0, qtc.Qt.AlignCenter) 
            self.firstDomMatrix.hide()
        
        sizes = [(matrix.sizeHint().width(), matrix.sizeHint().height())
                 for matrix in self.matrices]
        (widths, heights) = zip(*sizes)
        
        [lay.setColumnMinimumWidth(col, max(widths))
         for col in range(lay.columnCount())]
        [lay.setRowMinimumHeight(row, max(heights))
         for row in range(lay.rowCount())]
        lay.setSizeConstraint(qtw.QLayout.SetFixedSize)  # This makes the window 
                                                         # not re-sizable
        self.setLayout(lay)
        self.layout().activate()
    
    def mousePressEvent(self, event):
        if event.button() == qtc.Qt.LeftButton:
            self.dragPosition = (event.globalPos()
                                - self.frameGeometry().topLeft())
    
    def mouseReleaseEvent(self, event):
        event.accept()
    
    def mouseDoubleClickEvent(self, event):
        event.accept()
    
    def mouseMoveEvent(self, event):
        if event.buttons() == qtc.Qt.LeftButton:
            self.move(event.globalPos() - self.dragPosition)
    
    def updateMatrices(self, domsHidden):
        if domsHidden:
            toHide = self.matrices[0]
            toShow = self.firstDomMatrix
        else:
            toHide = self.firstDomMatrix
            toShow = self.matrices[0]
        toHide.hide()
        toShow.show()


class HoverAnnot(qtw.QWidget):
    # Arrow properties:
    headW = 6  # total width of arrowhead
    headL = 1.2*headW  # total length of arrowhead
    arrowL = 40  # total length of arrow segment
    
    # quantities used for painting the frame around the text
    framePad = (5, 3)  # amount of space added around the text
    
    # Painting styles and colors:
    bgColor = qtg.QColor(qtc.Qt.white)
    bgColor.setAlphaF(0.8)
    lineColor = qtg.QColor.fromHsvF(0, 0, 0.15)  # somewhat gray
    lineWidth = 1.5
    gap = 3 + lineWidth
    
    paintPen = qtg.QPen()
    paintPen.setWidthF(lineWidth)
    paintPen.setBrush(lineColor)
    paintPen.setCapStyle(qtc.Qt.RoundCap)
    paintPen.setJoinStyle(qtc.Qt.RoundJoin)
    
    def __init__(self, parent=None, text='Test some longer text', center=(500, 500)):
        super().__init__()
        
        self.setAttribute(qtc.Qt.WA_ShowWithoutActivating)
        self.setAttribute(qtc.Qt.WA_TranslucentBackground)
        # With the next attribute, even if the widget is too large mouse
        # events will be passed on to the widgets directly below it (spatially)
        self.setAttribute(qtc.Qt.WA_TransparentForMouseEvents)
        self.setParent(parent)
        
        self.setWindowFlags(qtc.Qt.Tool
                            | qtc.Qt.CustomizeWindowHint
                            | qtc.Qt.FramelessWindowHint
                            | qtc.Qt.WindowDoesNotAcceptFocus)
        
        self.setFont(gl.AllGUIFonts().labelFont)

        # Global position of the point that, together with the arrow head_pos
        # defines the direction of the arrow (from head_pos to center)
        self.__center = center
        
        self.__head_pos = None  # Global position of the arrowhead
        self.__tip_to_center = None  # distance from tip to center of text box
        
        self.__text = ''
        self.frame = qtg.QPainterPath()  # Frame that will contain the text,
        self.framePos = qtc.QPoint()     # its local position (top-left),
        self.frameRadius = 0             # the curvature radius of its corners
        self.textBbox = qtc.QRectF()     # and the boundingRect of text
        
        self.text = text
    
    @property
    def center(self):
        """
        Global position of the point towards which the arrow points
        """
        return self.__center
    
    @center.setter
    def center(self, ctr):
        """
        Set the global position of the point towards which the arrow points to
        ctr
        """
        if (not gl.check_type(ctr, 'arraylike')
            and not isinstance(ctr, (qtc.QPoint, qtc.QPointF))):
            raise TypeError("head_pos must be array-like, QPoint or QPointF."
                            f"Found {type(ctr)} instead.")
        if isinstance(ctr, (qtc.QPoint, qtc.QPointF)):
            ctr = (ctr.x(), ctr.y())
        
        try:  # this accounts for the case in which head_pos is not defined yet
            delta_head = self.head_pos - self.center
        except TypeError:
            delta_head = None
        self.__center = np.asarray(ctr)
        
        if delta_head is not None:
            self.head_pos = self.center + delta_head
    
    @property
    def text(self):
        return self.__text
    
    @text.setter
    def text(self, text):
        if not isinstance(text, str):
            raise
        self.__text = text
        self.updateFrameOutline()
        self.updateFramePosition()
    
    @property
    def head_pos(self):
        """
        Global position of the arrowhead
        """
        return self.__head_pos
    
    @head_pos.setter
    def head_pos(self, pos):
        """
        Set the global position of the arrowhead to pos
        """
        if (not gl.check_type(pos, 'arraylike')
            and not isinstance(pos, (qtc.QPoint, qtc.QPointF))):
            raise TypeError("head_pos must be array-like, QPoint or QPointF."
                            f"Found {type(pos)} instead.")
        if isinstance(pos, (qtc.QPoint, qtc.QPointF)):
            pos = (pos.x(), pos.y())
        
        if not len(pos) == 2:
            raise ValueError("Invalid position. An array-like with length == "
                             f"2 is needed. Found {len(pos)} instead.")

        self.__head_pos = np.asarray(pos)

        self.updateFramePosition()
    
    @property
    def tip_pos(self):
        """
        Returns the correct local position of the tip of the arrow such that
        the whole annotation can be painted within the area of the widget.
        Painting of the arrow will start at this position.
        """
        theta = np.radians(self.angle)
        sin_theta = np.sin(theta)
        cos_theta = np.cos(theta)
        frame_h = self.frame.boundingRect().height() + self.lineWidth
        frame_w = self.frame.boundingRect().width() + self.lineWidth

        tip_to_center = self.__tip_to_center + self.lineWidth
        
        delta_h = frame_h/2 + tip_to_center*sin_theta
        delta_w = frame_w/2 + tip_to_center*cos_theta
        
        # set the minimum delta to 2 to take care of the arrowhead triangle
        # +1 px for anti-alias, and one more pixel on h to avoid issues when the
        # arrow is in the bottom-left quadrant and exiting from the short side
        delta_w = max(2, delta_w) + 1
        delta_h = max(2, delta_h) + 2
        
        # notice that I return rounded values so that, later on, the widget can
        # be moved correctly to this position: only int values are acceptable
        # for widget movements
        return qtc.QPointF(np.round(delta_w), np.round(delta_h))
    
    @property
    def direction(self):
        try:
            dir_vec = self.center - self.head_pos
        except TypeError:
            # one of the two above is undefined
            return None
        return dir_vec/np.linalg.norm(dir_vec)
    
    @property
    def angle(self):
        """
        Returns the angle in degrees
        """
        dir_vec = self.direction
        if dir_vec is None:
            return None
        return np.degrees(np.arctan2(dir_vec[1], dir_vec[0]))
    
    def sizeHint(self):
        self.ensurePolished()
        rect_size = self.frame.boundingRect().size()
        delta = self.__tip_to_center + self.lineWidth
        w = np.round(rect_size.width())
        h = np.round(rect_size.height())
        # notice that the return value is slightly over-sized even in the worst
        # case scenario, as one would strictly need to add only w/2 and h/2
        return qtc.QSize(w + delta, h + delta)
    
    def minimumSizeHint(self):
        return self.sizeHint()
    
    def updateFrameOutline(self):
        self.frame = qtg.QPainterPath()
        self.textBbox = qtg.QFontMetricsF(self.font()).boundingRect(self.text)
        
        bbox = qtc.QRectF(self.textBbox)
        bbox.adjust(-self.framePad[0], -self.framePad[1],
                     self.framePad[0], self.framePad[1])
        bbox.moveTo(self.framePad[0]/2, self.framePad[1]/2)
        self.frameRadius = 0.2*bbox.height()
        self.frame.addRoundedRect(bbox, self.frameRadius, self.frameRadius)
        
        self.textBbox.adjust(-0.5, -0.5, 0.5, 0.5)
        self.textBbox.moveTo(self.framePad[0] + 1, self.framePad[1] + 1)
    
    def updateFramePosition(self):
        # Get the right position for self.frame such that the center of the 
        # text box is coincident with the prolongation of the arrow
        
        if self.angle is not None:
            bbox = self.frame.boundingRect()
            w = bbox.width()
            h = bbox.height()
            
            if np.abs(np.sin(np.radians(self.angle))) > np.radians(1):
                # away from the horizontal direction by more than 1 degree
                dy = -h/2 * np.sign(self.direction[1])
            else:
                dy = 0
            
            if dy == 0:
                dx = -w/2 * np.sign(self.direction[0])
            elif np.abs(np.cos(np.radians(self.angle))) > np.radians(1):
                # away from the vertical direction by more than 1 degree
                dx = dy / np.tan(np.radians(self.angle))
            else:
                dx = 0
            
            if np.abs(dx) > w/2:
                # Shallow horizontal angle, the arrow should come out from the
                # vertical edge rather than from the horizontal one
                dx = -w/2 * np.sign(self.direction[0])
                dy = dx * np.tan(np.radians(self.angle))
            
            # Set up the gap between the arrow and the outer edge of the box, 
            # including a correction for the rounded corner
            gap = self.gap  # takes into account line width
            tt = np.abs(np.tan(np.radians(self.angle)))
            ct = np.abs(np.cos(np.radians(self.angle)))
            st = np.abs(np.sin(np.radians(self.angle)))
            r = self.frameRadius
            if tt > max((h - 2*r)/w, h/(w - 2*r)):
                (ddx, ddy) = (0.5*gap*ct, 0.5*gap*(1 + st))
            elif tt < min((h - 2*r)/w, h/(w - 2*r)):
                (ddx, ddy) = (0.5*gap*(1 + ct), 0.5*gap*st)
            else:  # arrow intersects at rounded corner -> different gap
                (ddx, ddy) = (0.8*gap*ct, 0.8*gap*st)
            
            dx -= ddx*np.sign(self.direction[0])
            dy -= ddy*np.sign(self.direction[1])
            
            dc = -(self.arrowL)*self.direction + np.array([dx, dy])
            dc -= np.array([2, 1.5])
            self.__tip_to_center = np.linalg.norm(dc)
            
            deltaCenter = self.tip_pos + qtc.QPointF(*dc)
            
            self.framePos = deltaCenter - qtc.QPointF(w/2, h/2)

    def paintEvent(self, event):
        painter = qtg.QPainter()
        
        painter.begin(self)
        
        painter.setRenderHint(qtg.QPainter.Antialiasing)
        painter.setFont(self.font())
        painter.setPen(self.paintPen)
        
        self.paintArrow(painter)
        self.paintTextBox(painter)
        
        painter.end()

    def paintArrow(self, painter):
        
        painter.save()
        painter.translate(self.tip_pos)
        painter.rotate(self.angle)
        
        painter.drawPolyline(qtc.QPointF(-self.headL, -self.headW/2),
                             qtc.QPointF(0, 0),
                             qtc.QPointF(-self.headL, self.headW/2))
        painter.drawLine(qtc.QPointF(0, 0), qtc.QPointF(-self.arrowL, 0))
        painter.restore()
    
    def paintTextBox(self, painter):
        painter.translate(self.framePos)
        
        painter.fillPath(self.frame, self.bgColor)
        painter.drawPath(self.frame)
        painter.drawText(self.textBbox, self.text)
    
    def update_position(self):
        """
        Moves the annotation such that self.tip_pos is at the same global
        position as self.head_pos. This function should be called before showing
        the annotation to make sure its position on the screen is correct
        """
        tip = self.tip_pos  # local position
        tip_global = self.mapToGlobal(qtc.QPoint(tip.x(), tip.y()))
        new_pos = self.pos()
        delta = qtc.QPoint(*self.head_pos) - tip_global
        new_pos += delta
        self.move(new_pos)

