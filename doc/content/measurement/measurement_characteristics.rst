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
    :width: 70%
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

Next is the preparation phase, in which devices are instructed to perform their
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

Measurements can be either energy-resolved or time-resolved. An energy-resolved
measurement will collect one measurement per controller for each requested
quantity at a certain energy. A time-resolved measurement, in contrast, will
remain at each energy for the desired amount of time and continue collecting
data during this time span.

Energy-resolved measurements
----------------------------

There are two types of energy-resolved measurements: the energy calibration and
the regular |LEED-IV| measurement. Data acquired in an energy-resolved
measurement is plotted against the beam energy.

Time-resolved measurements
--------------------------

Time-resolved measurements are divided into two subcategories: continuous and
triggered. Data is plotted against the time elapsed since the measurement was
initiated.

Step profile
============

When transitioning from one energy to another, the measurement waits for settle
times to elapse before triggering data acquisition in the controlled devices.
Settle times represent an unavoidable slowdown due to transients following a
change in beam energy.

By controlling how the beam energy is stepped, the overall time until the beam
energy becomes stable can be reduced for some electronics. To support this,
ViPErLEED provides several step profiles:

* **Abrupt** — A direct jump from one energy to another. This is the default.
* **Linear** — A transition using equidistant intermediate steps. Useful when
  lower intermediate energies with very short intermediate settle times reduce
  the total settle time.
* **Custom** — Allows specifying arbitrary intermediate energies and settle
  times. This can, for example, be used to deliberately overshoot in order to
  create destructive interference between transients.

.. _fig_profile:
.. figure:: /_static/gui/profile.svg
    :width: 70%
    :align: center

    The step profile can take a more complex form than a simple abrupt jump.

In total, the time until the next measurement is the duration of the step
profile plus the settle time. For cameras, the settle time is always the HV
settle time. For controllers, the applicable settle time depends on the
measurement type.

Device handling
===============

ViPErLEED is designed to control more than just a single controller and camera.
Devices are managed asynchronously, which means devices may return data at any
time and ViPErLEED will process data accordingly. At the same time, ViPErLEED
ensures that commands are given at appropriate times by tracking whether
devices are busy or not. Despite their asynchronous operation, device data
acquisition is automatically synchronized whenever required.

While cameras only acquire images and therefore do not differentiate in
behaviour, controllers may interact with the LEED control unit. It is therefore
necessary to differentiate between the **primary controller** and **secondary
controllers**. The **primary controller** is responsible for setting the beam
energy during a measurement and it may optionally acquire data as well. Other
devices are only triggered to perform data acquisition after the primary
controller reports that the target beam energy has been set. **Secondary
controllers** perform data acquisition only.

.. _fig_handling:
.. figure:: /_static/gui/hardware_handling.svg
    :width: 100%
    :align: center

    Device handling in ViPErLEED. Physical hardware is shown on the right.
    Collected data is passed to the GUI on the left via the software
    implementations of the devices in the center.

Each hardware device has a corresponding software abstraction. This software
representation defines the tasks required to control the device and is managed
by the measurement. Communication between the software representation and the
physical hardware is handled by a serial layer or device driver, which
translates these tasks into hardware-specific instructions.

Images collected by cameras are shown in the camera views, while numerical data
from controllers is presented in the data plot. All collected data is compiled
by the measurement and stored on disk once the measurement has finished.

----

For further information, visit :ref:`Best Practice <best_practice>`.
