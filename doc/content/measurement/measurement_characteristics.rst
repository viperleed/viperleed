.. include:: /substitutions.rst

.. _measurement_characteristics:

===========================
Measurement characteristics
===========================

A *measurement*, as described in this section, is the entity that manages the
actual measurement process executed by the ViPErLEED measurement package. The
measurement object knows all devices involved, issues commands to them, and
receives data from them. Before a measurement can begin, it must be able to
communicate with the required hardware. All communication during a measurement
is performed through dedicated handlers for each device. These handlers allow
the measurement to issue generic commands to the software representations of
cameras and controllers. Each generic command is then translated into
device-specific instructions, fully abstracting the communication layer. This
design enables the use of any camera or controller, provided that its handler
is implemented on the software side.

Functional blocks
=================

A single measurement step consists of setting an energy, acquiring data, and
processing it. Processed numerical data is displayed as a plot. The latest
processed image from each camera is shown in the camera view. After processing,
the measurement decides if it has to collect more data by returning to setting
the energy, or if it can wrap things up and save the collected data. At any
point an error may occur, which will force the measurement to stop and attempt
to save the data it has collected so far. Detected errors are reported in the
GUI.

The functionality of the different measurement types can be broken down into
functional blocks, as shown in the figure below.

.. _fig_measurement_flow:
.. figure:: /_static/gui/Measurement_flow.svg
    :width: 90%
    :align: center

    Functional flow chart illustrating the blocks that make up a measurement.
    Progression goes from top to bottom.

Initialization
--------------

The first functional block is the initialization phase of the measurement. In
this phase settings are read from the configuration files and communication
with the hardware is established.

Preparation
-----------

Next is the preparation phase in which devices are instructed to perform their
own preparation and calibration. The preparation phase is complete once all
devices have completed their setup and are no longer busy. During preparation,
the beam energy is first set to the desired starting value.

Measurement loop
----------------

When the preparation has been completed, the measurement loop is engaged. At
the beginning of each iteration of this loop, the next beam energy of the
measurement is set. This can be a no-op if the required energy has already been
set, as is the case right after the end of the preparation. After setting an
energy, the measurement will wait for the respective settle times to elapse
before acquiring data. Setting the energy can be as simple as doing a single
abrupt energy step, or as complicated as setting multiple intermediate energies
with different times to wait afterwards.

After the target energy has been reached and the settle time has elapsed,
controllers and cameras begin acquiring data. Depending on the measurement type
this can be in triggered mode or continuous mode. In triggered mode the devices
will acquire a set number of measurements and return data once per request. In
continuous mode they will keep returning data until stopped.

After data acquisition, data is plotted in real time at the end of each
measurement step. The desired quantity to be plotted can be selected via the
data plot while the measurement is running. Cameras update their view with the
most recent image as well. If the measurement settings require further data
acquisition, then the measurement proceeds to setting the next beam energy. If
the desired end point has been reached, then the measurement class will enter
the finalization phase.

Finalization
------------

During finalization, all collected data is saved for later use. The devices are
disconnected, and the ViPErLEED package is returned to a state ready to start a
new measurement. The beam energy is always set to zero at the end of a
measurement. The full set of configuration files is stored alongside the
collected data to allow exact reproduction of the measurement.

Error handling
--------------

Errors may happen at any point during a measurement. A descriptive error report
to aid during troubleshooting has been added for many errors. Whenever able to,
ViPErLEED will attempt to salvage as much data as possible despite an error and
save it. The system is returned to a safe state by setting the electron beam
energy to zero. Many possible errors have been accounted for and will lead to
an abortion of the measurement. You are also able to abort the measurement at
any time, if required, which will also save all collected data and set the beam
energy to zero.

Energy-resolved and time-resolved measurements
==============================================

Measurements can either be energy-resolved or time-resolved. An energy-resolved
measurement will collect one measurement per controller for each requested
quantity at a certain energy. A time-resolved measurement in contrast will
remain at each energy for the desired amount of time and continue collecting
data during this time span.


For further information, visit :ref:`Best Practice <best_practice>`.
