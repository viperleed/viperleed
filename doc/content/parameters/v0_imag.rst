.. include:: /substitutions.rst

.. _v0_imag:

V0_IMAG
=======

V0_IMAG defines the imaginary part of the inner potential (in electronvolts),
which is combined with the real part :ref:`V0_REAL<MUFTIN>`. (The original
equivalent Fortran parameters are VPI, VPIS, VPIO.)

**Default**: V0_IMAG = 4.5

**Syntax**:

::

   V0_IMAG = 4.0

**Accepted input**: Positive real value.

.. note::
   -  Typical values are between 4 and 5 eV.
   -  The value and energy dependence of ``V0_IMAG`` can be *estimated* as
      the half width at half maximum (FWHM) of non-overlapping peaks in the
      experimental |IV| curves :cite:p:`heinzElectronBasedMethods2013`.
