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

If you are using the ViPErLEED hardware controller, make sure that it is
grounded at all times during operation. Do not connect the beam-energy
(=filament voltage) output of the LEED electronics to the ViPErLEED hardware if
you want to measure the beam current :math:`I_0`. The high-voltage input of the
hardware draws current and will falsify the beam-current measurement.

Camera
======

Make sure that you have recently (within the last few months) performed a bad
pixel calibration and ensure that this file is applied to your images.

Ensure that the lens aperture is fully open.

Before acquiring data with a camera, activate a camera live view and manually
go through the entire energy range you plan to measure at the desired gain. If
the camera detects overexposure (indicated by marking overexposed pixels red
and a red frame around the image) at any energy, decrease the gain (or, at the
minimum gain, the exposure time). If the maximum intensity is well below the
saturation limit, increase the exposure time unless the increase in measurement
time is an issue (e.g. because of electron-beam damage). Otherwise, increase
the gain. The goal is to use the full dynamic range of the sensor without
reaching pixel saturation.

Energy calibration
==================

If your energy calibration remains stable over multiple days (after giving the
electronics sufficient time to warm up on each day), you can skip redoing the
energy calibration; otherwise, it is usually sufficient to perform the
calibration once before acquiring |LEED-IV| data. To eliminate thermal drift,
the energy calibration can be done between measurements as well.

Use longer settling times (≈ 2 s) and larger energy steps (ΔE ≈ 10 eV) during
the energy calibration to keep it fast. Energy stability is more important than
dense sampling; additional data points do not substantially improve the
calibration quality.

The energy calibration can be performed without a sample in front of the
electron gun. This is especially useful if your sample degrades quickly under
electron-beam exposure. For a quick workflow, prepare and align the sample as
usual, then move it out of the beam path along one axis and perform the
calibration. Once complete, you can immediately reposition the sample and begin
your |LEED-IV| measurement. An alternative is to set a very negative Wehnelt
voltage to repel emitted electrons.

|LEED-IV| measurements
======================

For |LEED-IV| measurements, use the shortest possible settling time that still
guarantees the transient is over, to keep the measurement fast. Furthermore, if
beam damage under prolonged exposure is a concern, you can reduce controller
ADC averaging for :math:`I_0` and camera exposure to further accelerate the
measurement. If the maximum intensity is well below the saturation limit of the
camera, increase camera gain. Take care to stay within the range of the analog
gain of your camera (24 dB for the The Imaging Source cameras we have evaluated
in the hardware paper **TODO add reference**). Note that reducing averaging and
exposure time both decrease the signal-to-noise ratio.

Measure the beam current to be able to normalize your data. If measuring the
beam current is not possible and you are using the ViPErLEED hardware
controller, you can measure the sample current instead. In this case, first
acquire the |LEED-IV| measurement with the sample grounded and without
measuring the sample current. Immediately afterwards, repeat the exact same
measurement with the sample connected to the sample-current input instead of
ground. Be aware that this will apply a +33 V bias to the sample, meaning that
images acquired with the sample-current input connected are taken at an
incorrect beam energy.

Use the |LEED-IV| measurement mode to measure the flat fields. Flat-field
images are required to correct the measurement data for inhomogeneities of the
grids and other factors leading to spatially dependent sensitivity (e.g., dust
on the camera sensor). To measure flat fields, direct the electron beam at the
polycrystalline sample holder to produce a featureless LEED image. Ensure that
the distance between the sample plate and the LEED screen when measuring the
flat field is equal to the distance between the sample surface and the LEED
screen during the corresponding |LEED-IV| measurement. Measure the height
difference between the sample and holder beforehand and move the sample holder
by that difference closer to the electron gun.

You can either use the |LEED-IV| measurement mode or the average of a few
single snapshots to obtain dark frames. These images will be used to correct
for pixel-dependent dark current of the camera as well as stray light emitted
by the electron gun and other light sources. To acquire dark frames, increase
the Wehnelt voltage such that no electrons leave the electron gun. Usually, a
single snapshot is sufficient for dark-frame subtraction. If your LEED
electronics output for the beam current :math:`I_0` has an energy-dependent
offset, acquiring a movie of dark frames also provides the opportunity to
measure this offset as a function of energy (named :math:`I_{00}`).

You can perform two or three |LEED-IV| measurements at slightly different
distances from the electron gun to average out the influence of the grid on the
*I*\(*V*\) curves (and noise). A distance variation of 2 mm between the
measurements is usually sufficient. A separate flat field movie should be
acquired for each of these |LEED-IV| measurements.

Data acquisition flow
=====================

To guarantee optimal data quality, it is best to follow the workflow outlined
below.

Check electronics and setup
---------------------------
Ensure everything is turned on and connected properly. Give the electronics
enough time to warm up, then perform an energy calibration. Check that the bad
pixel file is up to date. Ensure that you know the necessary settling times.
Make sure as little stray light as possible enters the chamber. If possible,
the filaments of hot-cathode vacuum gauges should be switched off.

Adjust exposure time
--------------------
Turn the LEED optics on and manually ramp through your desired energy range
with the sample in place. Set the exposure time such that no or only very few
isolated pixels are saturated at any energy.

Energy calibration
------------------
If the energy stability of the LEED electronics is very poor, perform the
energy calibration shortly before the actual measurement. In this case, turn
the emitted beam off with a negative Wehnelt voltage to avoid any influence on
the sample by electrons or electron-desorbed molecules.

|LEED-IV| measurement
---------------------
Acquire the |LEED-IV| measurement. If you can afford the time, perform two or
three |LEED-IV| measurements at different electron-gun-to-sample distances.

Flat field
----------
Acquire the flat field for each of the |LEED-IV| measurements. Do not forget to
correct the position for height differences.

Dark frame
----------
Set a sufficiently negative Wehnelt voltage to suppress the beam and acquire a
few dark frames (as a short movie in time-dependent mode) or acquire dark
frames for each beam energy in |LEED-IV| mode. Exposure time and camera gain
must be the same as for the |LEED-IV| measurements. If you use higher gain or
exposure time for the flat-field movie, take an additional dark-frame movie
with these settings as well.
