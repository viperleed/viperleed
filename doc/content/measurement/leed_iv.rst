.. include:: /substitutions.rst

.. _leed_iv:

======================
|LEED-IV| measurements
======================

|LEED-IV| measurements are the core functionality of the ViPErLEED measurement
package. Before performing such a measurement, make sure you have followed all
steps described in :ref:`Best Practice <best_practice>` to ensure high data
quality and to reduce acquisition time.

A |LEED-IV| measurement requires a camera to operate. A controller capable of
measuring the beam current :math:`I_0` is strongly recommended, as it allows
the acquired image intensities to be normalized by :math:`I_0`. (This is done
by the ViPErLEED Spot Tracker :cite:t:`schmidViPErLEEDPackageII2025`.)
