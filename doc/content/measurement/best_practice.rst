.. include:: /substitutions.rst

.. _best_practice:

=============
Best practice
=============

The energy calibration should be performed regularly, at least once per day
right before starting the data acquisition process. To eliminate the thermal
drift, the energy calibration can be done in between measurements as well.

To determine how long the thermal drift remains significant after powering up
the LEED optics, you can perform a triggered time-resolved measurement to
monitor the beam energy drift over several hours.

The strongest non-linearities in LEED optics typically occur at the lowest beam
energies. These are also the energies with the fewest diffraction spots and
therefore the least data. Whenever possible, avoid these highly non-linear
regions when setting the calibration range.

Use longer settle times and larger energy steps (ΔE) to keep the energy
calibration fast. Energy stability is more important than dense sampling;
additional data points do not meaningfully improve the calibration quality.

The energy calibration can be performed without a sample in front of the
electron gun. This is especially useful if your sample degrades quickly under
electron-beam exposure. For a quick setup, prepare and align your sample as
usual, then move it out of the beam path along one axis and perform the
calibration. Once complete, you can immediately reposition the sample and begin
your |LEED-IV| measurement. An alternative is to set a very negative Wehnelt
voltage to repel emitted electrons.

Note that some electronics may require a short time for capacitors to discharge
and dipoles to relax after operating at high energies. In such cases, it may be
preferable to allow a brief downtime between measurements.

Make sure that you have recently (within the last few months) performed a bad
pixel calibration, and ensure that this file is applied to your images.

For |LEED-IV| measurements, use the shortest possible settle time that still
guarantees the transient is over, to keep the measurement fast. Furthermore, if
beam damage under extended exposure is an issue, you can reduce controller
sampling and increase camera gain (to reduce exposure time) to make the
measurement even faster. Note that reduced sampling and increased gain raise
the noise-to-signal ratio.

Before acquiring data with a camera, activate a camera live view and manually
go through the entire energy range you plan to measure at the desired gain.
Make sure the aperture is fully opened. If the camera detects overexposure at
any energy, decrease the exposure time. You want the camera to use the full
dynamic range without reaching overexposure.

For beam current measurements with the ViPErLEED hardware controller never
connect the beam energy port to your LEED optics simultaneously. The voltage
divider would draw current and therefore falsify the measurement.
