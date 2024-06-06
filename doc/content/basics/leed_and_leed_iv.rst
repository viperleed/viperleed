.. include:: /substitutions.rst

.. _leed_and_leed_iv:

|LEED IV|
=========

|LEED IV| [also sometimes referred to as LEED :math:`I(E)`\ ] is the
quantitative analysis of low-energy electron diffraction (LEED) patterns.
For a detailed introduction on the theory and applications of LEED, and
|LEED IV| in particular, we suggest, for example, Chapter 4 in
:cite:t:`fausterSurfacePhysicsFundamentals2020`
or the overview by :cite:t:`heinzElectronBasedMethods2013`.

In essence, the goal of any |LEED-IV| calculation is the calculation of
energy-dependent electron-scattering amplitudes and intensities of diffracted
electron beams.
These intensity curves, taken as a function of the primary beam energy,
[often referred to as |IV| curves or spectra] are very sensitive to the
precise position and vibrational amplitudes of each atom in the surface
unit cell.

In ViPErLEED, these calculations are performed by the |calc| Python package,
which acts as a comprehensive wrapper and major feature extension to the
established :term:`TensErLEED` program package. Based on the input structure
and the desired calculation parameters, ViPErLEED can calculate the |LEED-IV|
spectra for a given surface structure (see also
:ref:`reference calculation<ref-calc>`), compare these to experimental data
(see also |R-factor| :ref:`calculation<r-factor_calculation>`), and perform
a structure optimization (see :ref:`search<sec_search>`) using the
:ref:`tensor-LEED approach<tensor_leed>`.

For computational details please have a look at the ViPErLEED paper (**TODO**)
and the original work describing TensErLEED by Blum and Heinz
:cite:p:`blumFastLEEDIntensity2001a`.
