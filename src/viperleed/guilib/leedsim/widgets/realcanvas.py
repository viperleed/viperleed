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
from matplotlib import colors as mpl_colors

from viperleed.guilib.basewidgets import MPLFigureCanvas


class RealCanvas(MPLFigureCanvas):

    def plotLattices(self):
        self.ax.cla()

        _win = self.window()  # this should be the LEEDPatternSimulator instance

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
        rotSurfL = rs.surf.get_rotated_lattice_points(rotation)

        # Plot surface as scatter
        self.ax.scatter(rotSurfL[:,0], rotSurfL[:,1],
                        s=8*windowScaleFactor, c='k', zorder=2)
        for unitVec in rotSurfB:
            self.ax.annotate("", xy=tuple(unitVec),
                             xytext=(0, 0),
                             arrowprops=dict(arrowstyle="->"))

        self.setAxLimits(rs.fov)
        self.ax.figure.canvas.draw_idle()
