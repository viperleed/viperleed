"""Module measurement_plot of viperleed.gui.measure.widgets.

Defines the MeasurementPlot widget for displaying numerical data
acquired during measurement.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    'Florian Dörr (@FlorianDoerr)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2021-12-10'
__license__ = 'GPLv3+'

from collections import defaultdict

from matplotlib import colormaps
from matplotlib.lines import Line2D
import numpy as np
from PyQt5 import QtCore as qtc
from PyQt5 import QtWidgets as qtw

from viperleed.gui.basewidgets import MeasurementFigureCanvas as Canvas
from viperleed.gui.measure.classes.datapoints import DataPoints
from viperleed.gui.measure.classes.datapoints import QuantityInfo
from viperleed.gui.measure.widgets.checkcombobox import CheckComboBox
from viperleed.gui.widgetslib import AllGUIFonts



# TODO: temporarily one can adjust here the structure in which a
# measurement will be plotted (flat or step-wise). Will then become
# a combo box.
SEPARATE_STEPS = False
COLOR_FRACTION = 0.7
MARKERSIZE = 4


_MARKERS = (
    ('o', 'full'), ('o', 'none'),
    ('v', 'full'), ('v', 'none'),
    ('*', 'full'), ('*', 'none'),
    )

_COLORS = (colormaps['Greys'], colormaps['Blues'],
           colormaps['Oranges'], colormaps['Greens'],
           colormaps['Purples'], colormaps['Reds'],)

def _marker_style(marker, fill, color):
    """Return kwargs for Line2D given a marker and color."""
    empty = fill == 'none'
    return {'marker': marker, 'linestyle': '', 'color': color,
            'markerfacecolor': 'white' if empty else color,
            'markeredgecolor': color, 'markersize': MARKERSIZE}


class MeasurementPlot(qtw.QWidget):
    """A class that allows plotting of data."""

    def __init__(self, parent=None):
        """Initialise class."""
        super().__init__(parent=parent)
        self._ctrls = {'quantities': PlotComboBox(),
                       'no_data': qtw.QLabel("NO DATA")}
        self._glob = {'plot_lines': defaultdict(dict),}
        self.__markers = _MARKERS
        self.__ctrl_color = {}
        self.__data_points = None
        self.__canvas = Canvas()

        self.setWindowTitle("Measurement data plot")
        self.__compose()
        self._ctrls['quantities'].check_changed.connect(self.plot_all_data)

    @property
    def data_points(self):
        """Return the data points plotted."""
        return self.__data_points

    @data_points.setter
    def data_points(self, data_points):
        """Set the data points to be plotted.

        Setting also clears all axes, leaving the selection
        of plotted quantities unchanged. Call self.clear()
        to untick all explicitly.

        Parameters
        ----------
        data_points : DataPoints
            The data points containing measurement data to be plotted.

        Raises
        ------
        TypeError
            if data_points is not a DataPoints instance.
        """
        if not isinstance(data_points, DataPoints):
            raise TypeError(
                f"{self.__class__.__name__}: invalid type "
                f"{type(data_points).__name__!r} for data_points. "
                "Expected a DataPoints object."
                )
        self.__data_points = data_points
        self.clear(uncheck_all=False)
        self.plot_all_data()

    @property
    def lines(self):
        """Return a dict of plotted lines."""
        return self._glob['plot_lines']

    @property
    def plotted_quantities(self):
        """Return plotted quantities."""
        return self._ctrls['quantities'].selected_quantities

    def clear(self, uncheck_all=True):
        """Clear the plot."""
        self.__canvas.ax.clear()
        self._glob['plot_lines'] = defaultdict(dict)
        self.__canvas.figure.tight_layout()
        self.__canvas.draw_idle()
        if uncheck_all:
            self._ctrls['quantities'].uncheck_all()

    def plot_all_data(self):
        """Plot data on the canvas for the first time.

        This function is called if there is no plotted
        data yet or if the quantity that has to be plotted
        is changed. This method is much slower than calling
        plot_new_data()
        """
        if not self.data_points or not self.data_points.has_data:
            self._ctrls['no_data'].show()
            return
        self.__canvas.ax.clear()
        self._glob['plot_lines'] = defaultdict(dict)

        if self.plotted_quantities:
            if self.data_points.is_time_resolved:
                self.__plot_all_time_resolved_data()
            else:
                self.__plot_all_energy_resolved_data()

        self.__update_labels_and_legend()
        self.__canvas.draw_idle()

    def plot_new_data(self):
        """Plot new data to already-plotted data.

        This function is called if there is already
        plotted data and new data has to be added to
        the plot.
        """
        if (not self.data_points
            or not self.data_points.has_data
                or not self.plotted_quantities):
            return

        if self.data_points.is_time_resolved:
            self.__plot_new_time_resolved_data()
        else:
            self.__plot_new_energy_resolved_data()

        self.__update_labels_and_legend()
        self.__canvas.draw_idle()

    def __compose(self):
        """Prepare widget."""
        layout = qtw.QGridLayout()
        self._ctrls['quantities'].addItems(QuantityInfo.get_axis_labels('y'))
        self._ctrls['quantities'].setFont(AllGUIFonts().buttonFont)
        self._ctrls['quantities'].ensurePolished()

        self._ctrls['no_data'].setFont(AllGUIFonts().labelFont)
        self._ctrls['no_data'].ensurePolished()

        layout.addWidget(self._ctrls['quantities'], 0, 0)
        layout.addWidget(self.__canvas, 1, 0)
        layout.addWidget(self._ctrls['no_data'], 1, 0,
                         alignment=qtc.Qt.AlignCenter)
        self.setLayout(layout)
        self.__canvas.figure.tight_layout()

    def __make_legend(self):
        """Create the legend for the displayed plot."""
        controllers = []
        legend_elements = []
        for ctrl, color in self.__ctrl_color.items():
            try:
                ctrl_name = ctrl.name
            except AttributeError:
                ctrl_name = ctrl
            controllers.append(ctrl_name)
            legend_elements.append(
                # Line2D([], [], linestyle='None', label=ctrl_name,
                       # marker='o',
                       # markerfacecolor=color(COLOR_FRACTION),
                       # markersize=MARKERSIZE)
                Line2D([], [], label=ctrl_name, color=color(COLOR_FRACTION),
                       linewidth=4)
                )
        if len(self.plotted_quantities) == 1:
            return legend_elements

        # Add a spacer
        legend_elements.append(
            Line2D([], [], linestyle='', label="")
            )

        for marker, quantity in zip(self.__markers, self.plotted_quantities):
            style = _marker_style(*marker, 'black')
            legend_elements.append(
                Line2D([], [], label=quantity.label, **style)
                )
        return legend_elements

    def __plot_all_energy_resolved_data(self):
        axes = self.__canvas.ax
        data, nominal_energies = (
            self.data_points.get_energy_resolved_data(*self.plotted_quantities)
            )
        has_data = bool(data)
        self._ctrls['no_data'].setVisible(not has_data)
        if not has_data:
            return

        self.__ctrl_color = dict(zip(data, _COLORS))

        for marker, quantity in zip(self.__markers, self.plotted_quantities):
            for ctrl, measurements in data.items():
                if quantity not in measurements:
                    continue
                ctrl_data = measurements[quantity]
                if len(nominal_energies) == len(ctrl_data):
                    energies = nominal_energies
                else:
                    energies = nominal_energies[:-1]

                color = self.__ctrl_color[ctrl](COLOR_FRACTION)
                style = _marker_style(*marker, color)
                self.lines[quantity][ctrl], = axes.plot(energies, ctrl_data,
                                                        **style)

    def __plot_all_time_resolved_data(self):
        axes = self.__canvas.ax
        data, _ = self.data_points.get_time_resolved_data(
            *self.plotted_quantities,
            separate_steps=SEPARATE_STEPS
            )
        has_data = bool(data)
        self._ctrls['no_data'].setVisible(not has_data)
        if not has_data:
            return

        self.__ctrl_color = dict(zip(data, _COLORS))

        for marker, quantity in zip(self.__markers, self.plotted_quantities):
            for ctrl, measurements in data.items():
                if quantity not in measurements:
                    continue
                ctrl_data = measurements[quantity]
                ctrl_times = measurements[QuantityInfo.TIMES]

                if not SEPARATE_STEPS:
                    color = self.__ctrl_color[ctrl](COLOR_FRACTION)
                    style = _marker_style(*marker, color)
                    self.lines[quantity][ctrl], = axes.plot(ctrl_times,
                                                            ctrl_data, **style)
                    continue

                # SEPARATE_STEPS
                colors = self.__ctrl_color[ctrl](
                    np.linspace(0.2, 0.8, len(ctrl_data))
                    )
                # Cannot construct a big array of arrays and plot in
                # parallel bacause each step may (and usually will)
                # have a variable number of data.
                for color, times, values in zip(colors, ctrl_times, ctrl_data):
                    style = _marker_style(*marker, color)
                    axes.plot(times, values, **style)

    def __plot_new_energy_resolved_data(self):
        axes = self.__canvas.ax
        data, nominal_energies = (
            self.data_points.get_energy_resolved_data(*self.plotted_quantities)
            )
        has_data = bool(data)
        self._ctrls['no_data'].setVisible(not has_data)
        if not has_data:
            return

        self.__ctrl_color = dict(zip(data, _COLORS))                      # TODO: do we need this?

        for marker, quantity in zip(self.__markers, self.plotted_quantities):
            for ctrl, measurements in data.items():
                if quantity not in measurements:
                    continue
                ctrl_data = measurements[quantity]
                color = self.__ctrl_color[ctrl](COLOR_FRACTION)
                style = _marker_style(*marker, color)
                lines = self.lines[quantity]
                if ctrl not in lines:
                    lines[ctrl], = axes.plot(nominal_energies[-1],
                                             ctrl_data[-1], **style)
                lines[ctrl].set_data(nominal_energies, ctrl_data)

        axes.relim()
        axes.autoscale_view()

    def __plot_new_time_resolved_data(self):
        axes = self.__canvas.ax
        data, _ = self.data_points.get_time_resolved_data(
            *self.plotted_quantities,
            separate_steps=SEPARATE_STEPS
            )
        has_data = bool(data)
        self._ctrls['no_data'].setVisible(not has_data)
        if not has_data:
            return

        self.__ctrl_color = dict(zip(data, _COLORS))                      # TODO: do we need this?

        for marker, quantity in zip(self.__markers, self.plotted_quantities):
            for ctrl, measurements in data.items():
                if quantity not in measurements:
                    continue
                ctrl_data = measurements[quantity]
                ctrl_times = measurements[QuantityInfo.TIMES]
                if not SEPARATE_STEPS:
                    color = self.__ctrl_color[ctrl](COLOR_FRACTION)
                    style = _marker_style(*marker, color)
                    lines = self.lines[quantity]
                    if ctrl not in lines:
                        lines[ctrl], = axes.plot(ctrl_times, ctrl_data,
                                                 **style)
                        continue
                    lines[ctrl].set_data(ctrl_times, ctrl_data)
                    continue

                # SEPARATE_STEPS
                colors = self.__ctrl_color[ctrl](
                    np.linspace(0.2, 0.8, self.data_points.nr_steps_total)
                    )
                color_idx = (len(ctrl_data) - 1) // len(colors)
                style = _marker_style(*marker, colors[color_idx])
                # Enough to plot the last step for each
                axes.plot(ctrl_times[-1], ctrl_data[-1], **style)

        if not SEPARATE_STEPS:
            axes.relim()
            axes.autoscale_view()

    def __update_labels_and_legend(self):
        """Update axes labels and legend, if necessary."""
        axes = self.__canvas.ax
        if axes.legend_ is not None:
            # Axes and labels already up to date
            return

        if self.data_points.is_time_resolved:
            xlabel = f"Time ({QuantityInfo.TIMES.units})"
        else:
            xlabel = (f"{QuantityInfo.ENERGY.label} "
                      f"({QuantityInfo.ENERGY.units})")
        axes.set_xlabel(xlabel)

        if not self.plotted_quantities:
            axes.legend_ = None
            axes.set_ylabel("")
            self.__canvas.figure.tight_layout()
            return

        if len(self.plotted_quantities) == 1:
            yname = self.plotted_quantities[0].label
        else:
            yname = self.plotted_quantities[0].common_label
        axes.set_ylabel(f"{yname} ({self.plotted_quantities[0].units})")
        axes.legend(handles=self.__make_legend(), loc='upper left')
        self.__canvas.figure.tight_layout()


class PlotComboBox(CheckComboBox):
    """A CheckComboBox which can QuantityInfo objects."""

    @property
    def selected_quantities(self):
        """Return checked QuantityInfo objects."""
        checked = self.selected_items
        for index, item in enumerate(checked):
            checked[index] = QuantityInfo.from_label(item)
        return checked

    def toggle_checked(self, name):
        """Toggle the checked state of an entry with given name."""
        item = self.model().findItems(name)[0]
        super().toggle_checked(name)
        quantity = QuantityInfo.from_label(name)
        if self.is_item_checked(item):
            for i in range(self.count()):
                label = self.model().item(i, 0).text()
                other_quantity = QuantityInfo.from_label(label)

                if quantity.units != other_quantity.units:
                    disable = self.model().findItems(label)[0]
                    disable.setEnabled(False)
        elif not self.selected_items:
            for i in range(self.count()):
                enable = self.model().item(i, 0)
                enable.setEnabled(True)
