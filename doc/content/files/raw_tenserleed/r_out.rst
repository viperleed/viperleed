.. _r_out:

.. include:: /substitutions.rst

R_OUT
=====

``R_OUT`` is a file created during the 
:ref:`R-factor calculation<r-factor_calculation>` in :term:`TensErLEED`. 
It contains the R factors between the input :math:`I(V)` data 
(:ref:`file EXPBEAMS<EXPBEAMS>`) and theory calculation for 
each beam and each requested value of the  real part of the 
inner potential :math:`V_{0r}`.

:math:`V_{0r}` effectively shifts experiment and theory :math:`I(V)` curves
against each other on the energy axis and needs to be optimized along with
structural parameters to find a best fit structure. In TensErLEED this is
generally performed by brute force over a grid, i.e. by calculating the 
R-factors for each step on the requested :math:`V_{0r}` grid. See the parameter
:ref:`IV_SHIFT_RANGE<IVSHIFTRANGE>` for details.

The contents of ``R_OUT`` are read by |calc| to extract the R factors.
