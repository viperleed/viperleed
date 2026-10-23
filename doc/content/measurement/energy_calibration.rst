.. include:: /substitutions.rst

.. _energy_calibration:

==================
Energy calibration
==================

For LEED I(V) studies, it is essential that the electron energy is correct. An
error of 0.5% on the energy scale will already have a clear impact on the
results. LEED electronics are usually not sufficiently accurate without
additional calibration. The ViPErLEED electronics therefore contain a dedicated
high-voltage input for energy calibration. For calibration (but **not** during
the actual measurement), this input should be connected to the filament voltage
output of the LEED electronics (for Specs ErLEED controllers that do not have
such an output, see :ref:`erleed_modification`). ViPErLEED includes an energy
calibration function that automatically determines the nominal energies
required to achieve the desired target energies.

The energy calibration is one of the available measurement types when starting
a new measurement. It ensures that the energy scale of |LEED-IV| data is
accurate and should be performed before any |LEED-IV| data acquisition.

Before starting an energy calibration, make sure that the LEED controller is
switched on for a sufficient duration to reach thermal equilibrium (typically
one hour; see **TODO link section where we explain this here** for how to
estimate this) and that your setup can measure the acceleration voltage of the
electron gun.

To perform the calibration for the correct energy range, select the desired
start and end beam energies for your |LEED-IV| measurement. Once an energy
calibration has been completed successfully, the resulting calibration data are
stored automatically and applied in all subsequent measurements.

If a calibration attempt yields implausible results, a warning will appear,
indicating that the energy calibration failed. In that case, the previously
stored calibration will continue to be used.

----

For further information, visit :ref:`Best Practice <best_practice>`.
