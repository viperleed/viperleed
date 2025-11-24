.. include:: /substitutions.rst

.. _time_resolved:

==========================
Time-resolved measurements
==========================

To observe the evolution of a measurement over time, it is necessary to collect
time-resolved data. The core difference between an energy-resolved measurement
and a time-resolved measurement is that the time-resolved measurement acquires
multiple data points at the same energy. While an energy-resolved measurement
collects one data point and then moves on, a time-resolved measurement remains
at that energy and continues collecting data at fixed intervals for a defined
duration.

Time-resolved measurements are useful when trying to determine how long the
LEED optics need to relax after setting a new energy, or when performing a
long-term measurement to observe the energy drift of the optics as they heat up
during extended use. Another typical use case is the measurement of *I*\(*t*\)
curves, for example during a phase transition induced by the electron beam or
by temperature. The measurement package can perform these kinds of
measurements with the time-resolved measurement type.

A time-resolved measurement can be performed in two ways: as a **continuous**
measurement or as a **triggered** (non-continuous) measurement.

Continuous time-resolved measurement
====================================

The continuous time-resolved mode enables very fast sampling of quantities and
is ideal for observing rapid signal changes in real time. In this mode, the
controllers continuously return data without interruption. A typical use case
is determining the settle times of the beam energy and beam current. After
setting a new beam energy, the electronics may exhibit transient behavior; the
duration and shape of these transients can be characterized using this mode.

This mode acquires measurements as fast as possible and creates a continuous
data stream to the PC. The continuous time-resolved measurement does not
currently support camera use.

The behavior of controllers differs between primary and secondary controllers.
Primary controllers restart their acquisition sequence at each energy step
(their behavior is identical at all energies). Secondary controllers, however,
continue measuring while the primary controller sets the beam energy. This
provides additional information about transients.

.. _fig_continuous:
.. figure:: /_static/gui/StepTimingContinuous.svg
    :width: 90%
    :align: center

    Step timing for secondary controllers in continuous measurements. A primary
    controller would repeat the behaviour of the first step at all other energy
    steps.

Triggered time-resolved measurement
===================================

The triggered time-resolved mode acquires data at regular intervals for each
energy, potentially over long durations. Data acquisition is initiated by
periodic triggers rather than running continuously. Typical use cases include
observing the thermal drift of LEED electronics over time or monitoring a phase
transition while controlling the sample temperature.

Unlike the continuous mode, the triggered mode supports camera use. Each frame
or data point is acquired according to the defined timing interval, ensuring
synchronization between devices and stable acquisition over extended periods.

.. _fig_triggered:
.. figure:: /_static/gui/StepTimingTriggered.svg
    :width: 90%
    :align: center

    Step timing for devices in a triggered measurement.

----

For further information, visit :ref:`Best Practice <best_practice>`.
