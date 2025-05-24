"""
=======================================
   ViPErLEED Graphical User Interface
=======================================
 *** module guilib.leedsim.widgets ***

Created: 2020-01-12
Author: Michele Riva

Blah blah TODO
"""

import numpy as np
import matplotlib as mpl
#OR
#from matplotlib.patches import Circle as MplCircle

from viperleed.guilib.basewidgets import MPLFigureCanvas
from viperleed.guilib.leedsim.widgets.hoverannot import HoverAnnot

TEST = False


class LEEDCanvas(MPLFigureCanvas):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        # Get rid of the spines completely in the LEED pattern
        self.setSpinesOn(False)

        # Define annotations that will appear when the mouse hovers above
        # LEED spots (including symmetry-equivalent). Notice that at most 12
        # symmetry-equivalent spots can exist at the same time
        self.annots = [HoverAnnot(self.parentWidget()) for _ in range(12)]

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

        # TODO: there must be a better way. The PathCollection
        # that underlies the scatter plot has a set_sizes(array)
        # method that one could use to rescale the size of the
        # points without every time throwing away stuff and rebuilding
        # from scratch
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
