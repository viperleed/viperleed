.. include:: /substitutions.rst

.. _best_practice:

=============
Best practice
=============

Electronics
===========

To determine how long the thermal drift remains significant after powering up
the LEED optics, you can perform a triggered time-resolved measurement to
monitor the beam-energy drift over several hours.

The strongest non-linearities in LEED optics typically occur at the lowest beam
energies. These are also the energies with the fewest diffraction spots and
therefore the least data. Whenever possible, avoid these highly non-linear
regions when setting the calibration range.

Note that some electronics may require a short time for capacitors to discharge
and dipoles to relax after operating at high energies. In such cases, it may be
preferable to allow a brief downtime between measurements.

Make sure that you have recently (within the last few months) performed a bad
pixel calibration and ensure that this file is applied to your images.

Before acquiring data with a camera, activate a camera live view and manually
go through the entire energy range you plan to measure at the desired gain.
Make sure the aperture is fully opened. If the camera detects overexposure at
any energy, decrease the exposure time. The goal is to use the full dynamic
range of the sensor without reaching pixel saturation. Lower gain and higher
exposure times are generally preferable.

Never connect the beam-current and beam-energy ports to your LEED optics
simultaneously. The voltage divider on the beam-energy port draws current and
will falsify the beam-current measurement.

Energy calibration
==================

The energy calibration should be performed regularly, at least once per day
right before starting the data acquisition process. To eliminate the thermal
drift, the energy calibration can be done in between measurements as well.

Use longer settle times and larger energy steps (ΔE) during the energy
calibration to keep it fast. Energy stability is more important than dense
sampling; additional data points do not meaningfully improve the calibration
quality.

The energy calibration can be performed without a sample in front of the
electron gun. This is especially useful if your sample degrades quickly under
electron-beam exposure. For a quick workflow, prepare and align the sample as
usual, then move it out of the beam path along one axis and perform the
calibration. Once complete, you can immediately reposition the sample and begin
your |LEED-IV| measurement. An alternative is to set a very negative Wehnelt
voltage to repel emitted electrons.

|LEED-IV| measurements
======================

For |LEED-IV| measurements, use the shortest possible settle time that still
guarantees the transient is over, to keep the measurement fast. Furthermore, if
beam damage under prolonged exposure is a concern, you can reduce controller
sampling and increase camera gain (to reduce exposure time) to further
accelerate the measurement. Note that reducing sampling and increasing gain
both decrease the signal-to-noise ratio.

Measure the beam current to be able to normalize your data. If measuring the
beam current is not possible and you are using the ViPErLEED hardware
controller, you can measure the sample current instead. In this case, first
acquire the |LEED-IV| measurement with the sample grounded and without
measuring the sample current. Immediately afterwards, repeat the exact same
measurement with the sample connected to the sample-current input instead of
ground. Be aware that this will apply a +33 V bias to the sample, meaning that
images acquired with the sample-current input connected are taken at an
incorrect beam energy.

Use the |LEED-IV| measurement mode to measure the flat field. The flat field is
required to correct the measurement data later. To measure it, direct the
electron beam at the sample plate to produce uniform illumination. Ensure that
the distance between the sample plate and the LEED screen when measuring the
flat field is equal to the distance between the sample surface and the LEED
screen during the corresponding |LEED-IV| measurement. Measure the sample
thickness beforehand and move the sample plate by that distance closer to the
electron gun.

You can either use the |LEED-IV| measurement mode or a single snapshot to
measure stray light (dark frames) emitted by the electron gun and other light
sources. To acquire dark frames, increase the Wehnelt voltage such that no
electrons leave the electron gun. Usually, a single snapshot is sufficient for
dark-frame subtraction.

You can perform three |LEED-IV| measurements at slightly different distances
from the electron gun to average out the influence of the grid. A distance
variation of 5 mm is generally sufficient. A separate flat field should be
acquired for each of these |LEED-IV| measurements.

Measurement flow
================

To guarantee optimal data quality, it is best to follow the workflow outlined
below.

Check electronics and setup
---------------------------
Ensure everything is turned on and connected properly. Give the electronics
enough time to warm up. Check that the bad pixel file is up to date. Ensure
that you know the necessary settle times. Make sure as little stray light as
possible enters the chamber.

Adjust exposure time
--------------------
Turn the LEED optics on and manually ramp through your desired energy range
with the sample in place. Set the exposure time such that no pixels saturate at
any energy.

Energy calibration
------------------
Perform the energy calibration for your desired energy range.

|LEED-IV| measurement
---------------------
Acquire the |LEED-IV| measurement. If you can afford the time, perform three
|LEED-IV| measurements at different electron-gun-to-sample distances.

Flat field
----------
Acquire the flat field for each of the |LEED-IV| measurements. Do not forget to
correct for the sample thickness.

Dark frame
----------
Take a snapshot of the dark frame or acquire dark frames for each beam energy.
