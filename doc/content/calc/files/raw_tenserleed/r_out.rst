.. include:: /substitutions.rst

.. _r_out:

R_OUT
=====

``R_OUT`` is a file created during the |R-factor|
:ref:`calculation<r-factor_calculation>` in :term:`TensErLEED`.
It contains the |R factor|\ s between the input |IV| data
(file :ref:`EXPBEAMS`) and theory calculation for
each beam and each requested value of the  real part of the
inner potential |V0r|.

|V0r| effectively shifts experiment and theory |IV| curves against each other
on the energy axis and needs to be optimized along with structural parameters
to find a best fit structure. In TensErLEED this is performed by brute force
over a grid, i.e., by calculating the |R factor|\ s for each step on the
requested |V0r| grid. See the :ref:`IVSHIFTRANGE` parameter for details.

The contents of ``R_OUT`` are read by |calc| to extract the |R factor|\ s.
