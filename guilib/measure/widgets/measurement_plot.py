"""Module measurement_plot of viperleed.guilib.measure.widgets.

===============================================
      ViPErLEED Graphical User Interface
===============================================

Created: 2021-12-10
Author: Michele Riva
Author: Florian Doerr

Defines the MeasurementPlot class.
"""

import PyQt5.QtCore as qtc
from PyQt5 import QtWidgets as qtw
import numpy as np
from matplotlib.cm import get_cmap
from matplotlib.lines import Line2D

from viperleed import guilib as gl
from viperleed.guilib.basewidgets import MeasurementFigureCanvas as Canvas
from viperleed.guilib.measure.datapoints import DataPoints, QuantityInfo
from viperleed.guilib.measure.widgets.checkcombobox import CheckComboBox

from matplotlib.markers import MarkerStyle


# TODO: temporarily one can adjust here the structure in which a
# measurement will be plotted (flat or step-wise). Will then become
# a combo box.
SEPARATE_STEPS = False
COLOR_FRACTION = 0.7


class MeasurementPlot(qtw.QWidget):
    """A class that allows plotting of data."""

    def __init__(self, parent=None):
        """Initialise class."""
        super().__init__(parent=parent)
        self._ctrls = {'quantities': PlotComboBox(),}
        self._glob = {'plot_lines': {}}
        self.__markers = ('o', 'P', '*', 'v', 'h', 's')
        self.__colors = [get_cmap('Greys'), get_cmap('Blues'),
                         get_cmap('Oranges'),  get_cmap('Greens'),
                         get_cmap('Purples'), get_cmap('Reds'),]
        self.__ctrl_color = {}
        self.__data_points = None
        self.__canvas = Canvas()

        self.__compose()
        self._ctrls['quantities'].check_changed.connect(self.plot_all_data)

    @property
    def data_points(self):
        """Return the data points plotted."""
        return self.__data_points

    @data_points.setter
    def data_points(self, data_points):
        """Set the data points to be plotted.

        Setting also clears all axes.

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
        self.__canvas.ax.clear()
        self._glob['plot_lines'] = {}
        self.__canvas.draw_idle()

    @property
    def plotted_quantities(self):
        """Return plotted quantities."""
        return self._ctrls['quantities'].selected_quantities

    def plot_all_data(self):
        """Plot data on the canvas for the first time.

        This function is called if there is no plotted
        data yet or if the quantity that has to be plotted
        is changed. This method is much slower than calling
        plot_new_data()
        """
        if not self.data_points or not self.data_points.has_data():
            return
        self.__canvas.ax.clear()
        self._glob['plot_lines'] = {}

        if self.data_points.is_time_resolved():
            self.__plot_all_time_resolved_data()
            label = f"Time ({QuantityInfo.TIMES.units})"
        else:
            self.__plot_all_energy_resolved_data()
            label = (f"{QuantityInfo.ENERGY.label} "
                     f"({QuantityInfo.ENERGY.units})")
        self.__canvas.ax.set_xlabel(label)

        if self.plotted_quantities:
            if len(self.plotted_quantities) == 1:
                self.__canvas.ax.set_ylabel(
                    f"{self.plotted_quantities[0].label} "
                    f"({self.plotted_quantities[0].units})")
            else:
                self.__canvas.ax.set_ylabel(
                    f"{self.plotted_quantities[0].common_label} "
                    f"({self.plotted_quantities[0].units})")
        self.__canvas.ax.legend(handles=self.__make_legend(), loc='center left', bbox_to_anchor=(1, 0.5))
        self.__canvas.draw_idle()

    def plot_new_data(self):
        """Plot new data to already-plotted data.

        This function is called if there is already
        plotted data and new data has to be added to
        the plot.
        """
        if not self.data_points or not self.data_points.has_data():
            return

        if self.data_points.is_time_resolved():
            self.__plot_new_time_resolved_data()
            label = f"Time ({QuantityInfo.TIMES.units})"
        else:
            self.__plot_new_energy_resolved_data()
            label = (f"{QuantityInfo.ENERGY.label} "
                     f"({QuantityInfo.ENERGY.units})")
        self.__canvas.ax.set_xlabel(label)

        if self.plotted_quantities:
            if len(self.plotted_quantities) == 1:
                self.__canvas.ax.set_ylabel(
                    f"{self.plotted_quantities[0].label} "
                    f"({self.plotted_quantities[0].units})")
            else:
                self.__canvas.ax.set_ylabel(
                    f"{self.plotted_quantities[0].common_label} "
                    f"({self.plotted_quantities[0].units})")
        self.__canvas.ax.legend(handles=self.__make_legend(), loc='center left', bbox_to_anchor=(1, 0.5))
        self.__canvas.draw_idle()

    def __compose(self):
        """Prepare widget."""
        layout = qtw.QVBoxLayout()
        self._ctrls['quantities'].addItems(QuantityInfo.get_axis_labels('y'))
        self._ctrls['quantities'].setFont(gl.AllGUIFonts().buttonFont)
        self._ctrls['quantities'].ensurePolished()

        layout.addWidget(self._ctrls['quantities'])
        layout.addWidget(self.__canvas)
        self.setLayout(layout)

    def __plot_all_energy_resolved_data(self):
        for i, measured_quantity in enumerate(self.plotted_quantities):
            marker = self.__markers[i]
            self._glob['plot_lines'][measured_quantity] = []
            try:
                data, nominal_energies = (
                    self.data_points.get_energy_resolved_data(
                        measured_quantity, include_energies=True
                        )
                    )
            except ValueError:
                self.__canvas.ax.text(0,0, 'No data')
                continue

            for j, ctrl_data in enumerate(data):
                color = self.__colors[j](COLOR_FRACTION)
                self._glob['plot_lines'][measured_quantity].extend(
                    self.__canvas.ax.plot(nominal_energies, ctrl_data,
                                          marker=marker, linestyle='',
                                          color=color)
                    )

    def __plot_all_time_resolved_data(self):
        fig = self.__canvas
        for i, measured_quantity in enumerate(self.plotted_quantities):
            marker = self.__markers[i]
            self._glob['plot_lines'][measured_quantity] = []
            try:
                data, times = (
                    self.data_points.get_time_resolved_data(
                        measured_quantity, separate_steps=SEPARATE_STEPS,
                        absolute_times=not SEPARATE_STEPS,
                        include_energies=False
                        )
                    )
            except ValueError:
                self.__canvas.ax.text(0, 0, 'No data')
                continue

            for j, (ctrl_times, ctrl_data) in enumerate(zip(times, data)):
                if SEPARATE_STEPS:
                    color = self.__colors[j](np.linspace(0.2, 0.8,
                                                         len(ctrl_data)))
                    for k, (step_time, step_data) in enumerate(zip(ctrl_times, ctrl_data)):
                        # Cannot construct a big array of arrays and
                        # plot in parallel bacause each step may (and
                        # usually will) have a variable number of data
                        fig.ax.plot(step_time, step_data, marker=marker,
                                    linestyle='', color=color[k])
                else:
                    color = self.__colors[j](COLOR_FRACTION)
                    self._glob['plot_lines'][measured_quantity].extend(
                        fig.ax.plot(ctrl_times, ctrl_data, marker=marker,
                                    linestyle='', color=color)
                        )

    def __plot_new_energy_resolved_data(self):
        fig = self.__canvas
        for i, measured_quantity in enumerate(self.plotted_quantities):
            marker = self.__markers[i]
            if measured_quantity not in self._glob['plot_lines']:
                self._glob['plot_lines'][measured_quantity] = []
            try:
                data, nominal_energies = (
                    self.data_points.get_energy_resolved_data(
                        measured_quantity, include_energies=True
                        )
                    )
            except ValueError:
                self.__canvas.ax.text(0,0, 'No data')
                continue

            n_controllers = len(data)
            # Loop over controllers
            for j, ctrl_data in enumerate(data):
                color = self.__colors[j](COLOR_FRACTION)
                if len(self._glob['plot_lines'][measured_quantity]) < n_controllers:
                    self._glob['plot_lines'][measured_quantity].extend(
                        fig.ax.plot(nominal_energies[-1], ctrl_data[-1],
                                    marker=marker, linestyle='', color=color)
                        )
                    continue
                line = self._glob['plot_lines'][measured_quantity][j]
                line.set_data(nominal_energies, ctrl_data)
        fig.ax.relim()
        fig.ax.autoscale_view()

    def __plot_new_time_resolved_data(self):
        fig = self.__canvas
        for i, measured_quantity in enumerate(self.plotted_quantities):
            marker = self.__markers[i]
            if measured_quantity not in self._glob['plot_lines']:
                self._glob['plot_lines'][measured_quantity] = []
            try:
                data, times = (
                    self.data_points.get_time_resolved_data(
                        measured_quantity, separate_steps=SEPARATE_STEPS,
                        absolute_times=not SEPARATE_STEPS,
                        include_energies=False
                        )
                    )
            except ValueError:
                self.__canvas.ax.text(0,0, 'No data')
                continue

            n_controllers = len(data)
            # Loop over controllers
            for j, (ctrl_times, ctrl_data) in enumerate(zip(times, data)):
                if SEPARATE_STEPS:
                    color = self.__colors[j](np.linspace(
                        0.2, 0.8, self.data_points.num_measurements
                        ))
                    # Enough to plot the last step for each
                    fig.ax.plot(ctrl_times[-1], ctrl_data[-1], marker=marker,
                                linestyle='', color=color[len(ctrl_data)-1])
                else:
                    color = self.__colors[j](COLOR_FRACTION)
                    if len(self._glob['plot_lines'][measured_quantity]) < n_controllers:
                        self._glob['plot_lines'][measured_quantity].extend(
                            fig.ax.plot(ctrl_times, ctrl_data, marker=marker,
                                        linestyle='', color=color)
                            )
                        continue
                    line = self._glob['plot_lines'][measured_quantity][j]
                    line.set_data(ctrl_times, ctrl_data)
        if not SEPARATE_STEPS:
            fig.ax.relim()
            fig.ax.autoscale_view()

    def __make_legend(self):
        """Create the legend for the displayed plot."""
        controllers = []
        quantities = {}
        # legend = [[], []]
        legend_elements = []
        j = 0
        for i, measured_quantity in enumerate(self.plotted_quantities):
            marker = self.__markers[i]
            # TODO: throws an error if key does not exist in data_points fix!!!!!
            for controller in self.data_points[0][measured_quantity]:
                if controller.serial.port_name not in controllers:
                    controllers.append(controller.serial.port_name)
                    legend_elements.append(Line2D([0], [0], marker='o', color='w', label=controller.serial.port_name, markerfacecolor=self.__colors[j](COLOR_FRACTION), markersize=7))
                    j += 1
            if measured_quantity not in quantities:
                quantities[measured_quantity] = marker
        for quantity, marker in quantities.items():
            legend_elements.append(Line2D([0], [0], marker=marker, color='w', label=quantity.label, markerfacecolor='#000000', markersize=7))
        return legend_elements


class PlotComboBox(CheckComboBox):
    """A CheckComboBox which can QuantityInfo objeczts."""

    check_changed = qtc.pyqtSignal()

    def __init__(self, parent=None):
        """ Initialize instance.

        Parameters
        ----------
        parent : QWidget
            The parent widget of this combo box.

        Returns
        -------
        None.
        """
        super().__init__(parent=parent)

    @property
    def selected_quantities(self):
        """Return checked QuantityInfo objects."""
        checked = self.selected_items
        for index, item in enumerate(checked):
            checked[index] = QuantityInfo.from_label(item)
        return checked

    def toggle_checked(self, name):     # pylint: disable=invalid-name
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
        else:
            if not self.selected_items:
                for i in range(self.count()):
                    enable = self.model().item(i, 0)
                    enable.setEnabled(True)
        self.check_changed.emit()

