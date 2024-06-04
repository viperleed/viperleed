.. include:: /substitutions.rst

.. _bulkdoubleiter:

BULKDOUBLING_MAX
================

BULKDOUBLING_MAX (corresponds to Fortran LITER) defines the maximum number
of layer-doubling iterations used to compute the scattering amplitudes of
the bulk. See :ref:`BULKDOUBLING_EPS<BULKDOUBLEEPS>`  for more information
on layer doubling.

**Default**: BULKDOUBLING_MAX = 10

**Syntax**:

::

   BULKDOUBLING_MAX = 8

**Acceptable values**: Positive integers

.. note::
    If ``BULKDOUBLING_MAX = n``, the layer-doubling scheme will end up stacking
    a maximum of :math:`2^n` layers for the bulk while seeking convergence
    of the scattering matrices. If no convergence is achieved with the default
    value (\ :math:`2^{10} = 1024` layers), it's likely that something else is
    wrong (e.g., geometry input/:ref:`POSCAR<POSCAR>` file). Large values for
    BULKDOUBLING_MAX will generally make the layer doubling converge to
    something, but this might be wrong (it might show up as very large noise
    on some |IV| curves, or some |IV| curves being a few orders of magnitude
    too much/too little intense)! In normal conditions, the convergence
    criterion on :ref:`BULKDOUBLING_EPS<BULKDOUBLEEPS>`  should be satisfied
    much earlier than after BULKDOUBLING_MAX iterations.
