"""Module measurement_plot of viperleed.guilib.measure.widgets.

===============================================
      ViPErLEED Graphical User Interface
===============================================

Created: 2021-12-10
Author: Michele Riva
Author: Florian Doerr

Defines the MeasurementPlot class.
"""

from PyQt5 import QtWidgets as qtw

from viperleed.guilib.basewidgets import MeasurementFigureCanvas as Canvas
from viperleed.guilib.measure.datapoints import DataPoints


# TODO: temporarily one can adjust here the structure in which a
# measurement will be plotted (flat or step-wise). Will then become
# a combo box.
SEPARATE_STEPS = False


class MeasurementPlot(qtw.QWidget):
    """A class that allows plotting of data."""

    def __init__(self, parent=None):
        """Initialise class."""
        super().__init__(parent=parent)
        self._contrls = {}
        self._glob = {'plot_lines': []}
        self.__markers = ('.', 'o', 'x', '+', '*')
        self.__data_points = None
        self.__canvas = Canvas()

        self.__compose()

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
        self._glob['plot_lines'] = []
        self.__canvas.draw_idle()

    def plot_all_data(self, measured_quantity):
        """Plot data on the canvas for the first time.

        This function is called if there is no plotted
        data yet or if the quantity that has to be plotted
        is changed. This method is much slower than calling
        plot_new_data()
        """
        if not self.data_points or not self.data_points.has_data():
            return
        self.__canvas.ax.clear()
        self._glob['plot_lines'] = []

        if self.data_points.is_time_resolved():
            self.__plot_all_time_resolved_data(measured_quantity)
        else:
            self.__plot_all_energy_resolved_data(measured_quantity)
        self.__canvas.draw_idle()

    # TODO: temporarily we will still use the measured_quantity
    # variable for plotting with the new class.
    # it is passed on from the Measure class
    def plot_new_data(self, measured_quantity):
        """Plot new data to already-plotted data.

        This function is called if there is already
        plotted data and new data has to be added to
        the plot.
        """
        if not self.data_points or not self.data_points.has_data():
            return

        if self.data_points.is_time_resolved():
            self.__plot_new_time_resolved_data(measured_quantity)
        else:
            self.__plot_new_energy_resolved_data(measured_quantity)
        self.__canvas.draw_idle()

    def __compose(self):
        """Prepare widget."""
        layout = qtw.QVBoxLayout()
        layout.addWidget(self.__canvas)

        self.setLayout(layout)

    def __plot_all_energy_resolved_data(self, measured_quantity):
        marker = self.__markers[0]
        data, nominal_energies = (
            self.data_points.get_energy_resolved_data(
                measured_quantity, include_energies=True
                )
            )
        for ctrl_data in data:
            self._glob['plot_lines'].extend(
                self.__canvas.ax.plot(nominal_energies, ctrl_data, marker)
                )

    def __plot_all_time_resolved_data(self, measured_quantity):
        fig = self.__canvas
        marker = self.__markers[0]
        data, times = (
            self.data_points.get_time_resolved_data(
                measured_quantity, separate_steps=SEPARATE_STEPS,
                absolute_times=not SEPARATE_STEPS, include_energies=False
                )
            )

        for ctrl_times, ctrl_data in zip(times, data):
            if SEPARATE_STEPS:
                for step_time, step_data in zip(ctrl_times, ctrl_data):
                    fig.ax.plot(step_time, step_data, marker)
            else:
                self._glob['plot_lines'].extend(
                    fig.ax.plot(ctrl_times, ctrl_data, marker)
                    )

    def __plot_new_energy_resolved_data(self, measured_quantity):
        fig = self.__canvas
        marker = self.__markers[0]  # TODO: will have to loop through different quantities
        data, nominal_energies = (
            self.data_points.get_energy_resolved_data(
                measured_quantity, include_energies=True
                )
            )
        n_controllers = len(data)
        # Loop over controllers
        for i, ctrl_data in enumerate(data):
            if len(self._glob['plot_lines']) < n_controllers:
                self._glob['plot_lines'].extend(
                    fig.ax.plot(nominal_energies[-1], ctrl_data[-1], marker)
                    )
                continue
            line = self._glob['plot_lines'][i]
            line.set_data(nominal_energies, ctrl_data)
        fig.ax.relim()
        fig.ax.autoscale_view()

    def __plot_new_time_resolved_data(self, measured_quantity):
        fig = self.__canvas
        marker = self.__markers[0]  # TODO: will have to loop through different quantities
        data, times = (
            self.data_points.get_time_resolved_data(
                measured_quantity, separate_steps=SEPARATE_STEPS,
                absolute_times=not SEPARATE_STEPS, include_energies=False
                )
            )
        n_controllers = len(data)
        # Loop over controllers
        for i, (ctrl_times, ctrl_data) in enumerate(zip(times, data)):
            if SEPARATE_STEPS:
                # Enough to plot the last step for each
                fig.ax.plot(ctrl_times[-1], ctrl_data[-1], marker)
            else:
                if len(self._glob['plot_lines']) < n_controllers:
                    self._glob['plot_lines'].extend(
                        fig.ax.plot(ctrl_times, ctrl_data, marker)
                        )
                    continue
                line = self._glob['plot_lines'][i]
                line.set_data(ctrl_times, ctrl_data)
        if not SEPARATE_STEPS:
            fig.ax.relim()
            fig.ax.autoscale_view()
