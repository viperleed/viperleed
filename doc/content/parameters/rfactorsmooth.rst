.. include:: /substitutions.rst

.. _rfactorsmooth:

R_FACTOR_SMOOTH
===============

R_FACTOR_SMOOTH defines the number of times that experimental |IV| curves
will be smoothed by 3-point smoothing during R-factor calculations, both
after the reference calculation and during the search.

**Default:** 0

**Allowed values:** 0 or positive integer

**Syntax:**

::

   R_FACTOR_SMOOTH = 4

It is recommended to smooth the experimental curves beforehand, such
that smoothing during TensErLEED operation is not necessary. See the
:ref:`R-factor<r-factor_calculation>` page for details.

If you do use smoothing through the R_FACTOR_SMOOTH parameter, check the
experimental curves in the :ref:`Rfactor_analysis.pdf<rfactoranalysis>` file.
This plots the curves after smoothing, as well as the Y-function used for
calculating the R factor.
