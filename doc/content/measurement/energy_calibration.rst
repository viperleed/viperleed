.. include:: /substitutions.rst

.. _energy_calibration:

==================
Energy calibration
==================

LEED electronics often exhibit a slight mismatch between the beam energy that
is nominally set and the acceleration voltage that is effectively applied to
the emitted electrons. To address this inherent characteristic of most LEED
optics, ViPErLEED includes an energy calibration function that automatically
determines the nominal energies required to achieve the desired target
energies.

The energy calibration is one of the available measurement types when starting
a new measurement. It ensures that the energy scale of |LEED-IV| data is
accurate and should be performed before any |LEED-IV| data acquisition.

Before starting an energy calibration, make sure that the LEED optics are
switched on and that your setup can measure the acceleration voltage of the
electron gun. Allow the electron gun filament sufficient time to heat up before
proceeding.

To perform the calibration for the correct energy range, select the desired
start and end beam energies for your |LEED-IV| measurement. Once an energy
calibration has been completed successfully, the resulting calibration data are
stored automatically and applied in all subsequent measurements.

If a calibration attempt yields implausible results, a warning will appear,
indicating that the energy calibration failed. In that case, the previously
stored calibration will continue to be used.

----

For further information, visit :ref:`Best Practice <best_practice>`.
