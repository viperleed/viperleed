.. include:: /substitutions.rst

.. _leed_iv:

======================
|LEED-IV| measurements
======================

|LEED-IV| measurements are the core functionality of the ViPErLEED measurement
package. Before performing such a measurement, make sure you have followed all
steps described in :ref:`Best Practice <best_practice>` to ensure high data
quality and to reduce acquisition time.

A |LEED-IV| measurement requires a camera for operation. A controller capable
of measuring the beam current is strongly recommended, as it allows the
acquired data to be normalized.

If measuring the beam current is not possible and you are using the ViPErLEED
hardware controller, you can measure the sample current instead to normalize
the data. In this case, first acquire the |LEED-IV| measurement with the
sample grounded and without measuring the sample current. Immediately
afterwards, repeat the exact same measurement with the sample connected to the
sample-current input instead of ground. Be aware that this will apply a +33 V
bias to your sample, meaning that images acquired with the sample-current input
connected are taken at an incorrect beam energy.

For further information, visit :ref:`Best Practice <best_practice>`.
