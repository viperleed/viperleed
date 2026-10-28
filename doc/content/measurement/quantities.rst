.. include:: /substitutions.rst

.. _quantities:

==========
Quantities
==========

Currents
========

**Beam current** |I0b|: The current emitted by the electron gun. Note that this
differs from the value displayed by LEED electronics as "emission current":
Emission current refers to the electrons emitted from the filament (and mostly
collected at the anode inside the electron source). Thus, the emission current
(mA) is much higher than the beam current :math:`I_0` (µA or nA).

**Sample current**: The beam current arriving at the sample. In order to
measure the beam current as the sample current, the sample must be positively
biased to prevent secondary electron emission.

Energies
========

**Beam energy**: The kinetic energy of the electron beam directed at the
sample, determined by the filament voltage, except for a constant offset due to
the work functions, which is accounted for in |LEED-IV| calculations. For
|LEED-IV| measurements, the beam energy must be accurate within a few tenths of
a percent. Since most LEED electronics are not sufficiently accurate, an
energy calibration is required.

**Start energy**: The beam energy at which the measurement begins. For
|LEED-IV| measurements, select the first energy at which diffraction spots
begin to appear.

**Delta energy / ΔE**: The energy difference between two consecutively acquired
energies. In other words, the step height between energy points. Smaller delta
energies increase the information density and are therefore recommended for
|LEED-IV| measurements (e.g., 0.5 eV). For |LEED-IV| studies, a step of 0.5 eV
is recommended. If steps are larger, i.e., to reduce the measurement time for
very sensitive samples, the experimental |IV| curves should be interpolated to
0.5 eV steps or finer for comparison with theory
:cite:p:`imre2026improvedrfactor`.

**End energy**: The final beam energy of the measurement, after which the
acquisition either concludes or restarts. This energy should be chosen with
great care. A wider energy range adds more data, but damage to sensitive
samples by the electron beam (i.e., desorption or cracking of adsorbates, and
electron-stimulated desorption of oxygen from oxide surfaces) tends to increase
with higher beam energies.

Times
=====

**HV settling time / Energy settling time**: The time it takes for the actual
electron energy to reach a stable value after setting a new energy. This is the
minimum time cameras in |LEED-IV| measurements, and controllers during the
energy calibration and triggered time-resolved measurements, will wait before
acquiring data.

|I0b| **settling time / Current settling time**: The duration of the transient
response of the beam current after setting a new energy. This is the minimum
time controllers will wait during |LEED-IV| measurements before acquiring data.
