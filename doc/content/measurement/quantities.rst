.. include:: /substitutions.rst

.. _quantities:

==========
Quantities
==========

Currents
========

|I0b| **/ Beam current**: The current emitted by the electron gun.

**Sample current**: The beam current arriving at the sample. The sample must be
biased to prevent secondary electron emission.

Energies
========

**Beam energy**: The kinetic energy of the electron beam directed at the
sample, determined by the acceleration voltage. For |LEED-IV| measurements, it
is important that the nominal (expected) energy matches the measured (real)
energy as closely as possible. The energy calibration measurement should be
performed if possible to improve the accuracy of the beam energy.

**Start energy**: The beam energy at which the measurement begins. For
|LEED-IV| measurements, select the first energy at which diffraction spots
begin to appear. If you experience charging effects at lower energies, select
an energy at which these effects no longer occur.

**Delta energy / ΔE**: The energy difference between two consecutively acquired
energies. In other words, the step height between energy points. Smaller delta
energies increase the information density and are therefore recommended for
|LEED-IV| measurements (e.g., 0.5 eV).

**End energy**: The final beam energy of the measurement, after which the
acquisition either concludes or restarts.

Times
=====

**HV settle time / Energy settle time**: The duration of the transient response
of the acceleration voltage after setting a new energy. This is the minimum
time cameras in |LEED-IV| measurements, and controllers during the energy
calibration, will wait before acquiring data.

|I0b| **settle time / Current settle time**: The duration of the transient
response of the beam current after setting a new energy. This is the minimum
time controllers will wait during |LEED-IV| measurements before acquiring data.
