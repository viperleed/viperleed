.. include:: /substitutions.rst

.. _best_practice:

=============
Best Practice
=============

The energy calibration should be performed regularly, at least once per day
right before starting the data acquisition process. To eliminate the thermal
drift, the energy calibration can be done in between measurements as well.

To determine how long the thermal drift remains significant after powering up
the LEED optics, you can perform a triggered time-resolved measurement to
monitor the beam energy drift over several hours.

The strongest non-linearities in LEED optics typically occur at the lowest beam
energies. Whenever possible, avoid these highly non-linear regions when setting
the calibration range.

Use longer settle times and larger energy steps (ΔE) to keep the energy
calibration fast. Energy stability is more important than dense sampling;
additional data points do not meaningfully improve the calibration quality.

The energy calibration can be performed without a sample in front of the
electron gun. This is especially useful if your sample degrades quickly under
electron beam exposure. For a quick setup, prepare and align your sample as
usual, then move it out of the beam path along one axis and perform the
calibration. Once complete, you can immediately reposition the sample and begin
your |LEED-IV| measurement. An alternative is to set a very negative Wehnelt
voltage to repel emitted electrons.

Note that some electronics may require a short time for capacitors to discharge
and dipoles to relax after operating at high energies. In such cases, it may be
preferable to allow a brief downtime between measurements.
