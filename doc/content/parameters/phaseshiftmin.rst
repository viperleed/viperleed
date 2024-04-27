.. _phaseshiftmin:

PHASESHIFT_EPS
==============

PHASESHIFT_EPS defines a cutoff criterion to determine the maximum angular
momentum number (LMAX=MLMAX in Fortran) to be used in the spherical harmonics
expansion, used for calculating the scattering matrices within each layer.

**Default**: 0.02

**Syntax**:

::

   PHASESHIFT_EPS = 0.03
   PHASESHIFT_EPS = r

**Acceptable values**: 0<``PHASESHIFT_EPS``\ <1, float, **OR** flags
**R**\ ough (0.1) / **N**\ ormal (0.05) / **D**\ efault (0.02) / **F**\ ine
(0.01) (only first letter of flags is read, case-insensitive).
Typical 0.001 to 0.1.

LMAX may further be limited by the :ref:`LMAX<LMAX>` parameter, which can be
used to define upper and lower bounds for acceptable LMAX, or also to fix LMAX
to a specific value (in the latter case, PHASESHIFT_EPS is ignored).
By default, LMAX is allowed to vary between 6 and 18.


Changelog
---------

.. versionchanged:: 0.12.0
      The default value of PHASESHIFT_EPS has been changed from 0.01 to 0.02.
