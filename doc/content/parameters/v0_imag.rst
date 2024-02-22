.. _v0_imag:

V0_IMAG
=======

V0_IMAG defines the imaginary part of the inner potential (in eV), which is combined with the real part :ref:`V0_REAL<MUFTIN>`. (The original equivalent Fortran parameters are VPI, VPIS, VPIO.)

**Default**: V0_IMAG = 4.5

**Syntax**:

::

   V0_IMAG = 4.0

**Accepted input**: Positive real value.

.. note::
   -  Typical values are between 4 and 5 eV.
   -  The value and energy dependence of ``V0_IMAG`` can be *estimated* as the half width at half maximum (FWHM) of non-overlapping peaks in the experimental IV curves :cite:p:`heinzElectronBasedMethods2013`.
      **TODO Alex, Michele, Michael**: perhaps this would be the way to get a decent guess for V0i? (Hmmm, if I try it on Cu(111), I get values between 6–8.5 eV, should be 4–4.5. On Pt111+Te I see few peaks that as narrow as the rule tells. At best, we might say that **sometimes**, it can be estimated ... of the **narrowest** non-overlapping peaks. -ms)
      **TODO Alex, Michele, Michael** we should not forget to look at the energy dependence at some point. Taking literature values, *V*\ <sub>0i<sub> for Pt changes from 3.1 to 5.7 eV between 50 and 250 eV. But that's something to do after publication!
