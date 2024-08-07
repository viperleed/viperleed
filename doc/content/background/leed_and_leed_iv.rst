.. include:: /substitutions.rst

.. _leed_and_leed_iv:

|LEED IV|
=========

|LEED IV| — also sometimes referred to as LEED *I*\ (\ *E*\ ) — is the
quantitative analysis of energy-dependent low-energy electron diffraction
(LEED) patterns. It involves the :ref:`acquisition<hardware>` of LEED
patterns as a function of the energy of the electron beam — usually
referred to as "LEED videos" or "LEED movies". The intensity of each "spot"
is then :ref:`extracted<imagej_plugins>` to obtain so-called |IV| curves —
sometimes also referred to as spectra. Such |IV| curves are exceptionally
sensitive to the precise position, vibrational amplitude, and chemical element
of each atom in the surface unit cell.
|IV| curves can also be :ref:`calculated<viperleed_calc>` based on a test
structural model of the surface. Comparing calculated and experimental |IV|
curves allows one to accept or discard structural models for surfaces with
picometre-level accuracy.

For a detailed introduction on the theory and applications of LEED and
|LEED IV| we suggest, for example, Chapter 4 in
:cite:t:`fausterSurfacePhysicsFundamentals2020`
or the overview by :cite:t:`heinzElectronBasedMethods2013`.

In essence, the goal of a |LEED-IV| calculation is the determination of
energy-dependent electron-scattering amplitudes and of the corresponding
intensities of diffracted electron beams.
In ViPErLEED, |LEED-IV| calculations are performed by the |calc| Python
package, which acts as a comprehensive wrapper and major feature extension
to the established :term:`TensErLEED` program package. Based on the input
structure and the desired calculation parameters, ViPErLEED can calculate
the |LEED-IV| spectra for a given surface structure (see also :ref:`ref-calc`),
compare these to experimental data (see also |R-factor|
:ref:`calculation<r-factor_calculation>`), and perform
a structure optimization (see :ref:`sec_search`)
using the :ref:`tensor-LEED approach<tensor_leed>`.

For computational details, have a look at the relevant ViPErLEED paper
:cite:p:`viperleedCalc` and at the original work describing TensErLEED
by :cite:t:`blumFastLEEDIntensity2001a`.
