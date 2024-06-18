.. include:: /substitutions.rst

.. _fdoptimization:

Full-dynamic optimization
=========================

|LEED-IV| calculations require setting some parameters that cannot be (exactly)
determined in experiment, but that are not accessible to variation within the
tensor-LEED approximation. These include :ref:`BEAMINCIDENCE`, :ref:`V0_IMAG`,
and the exact dimensions of the unit cell. In order to fit these parameters,
multiple full-dynamic (i.e., :ref:`"reference"<ref-calc>`) calculations must
be performed, instead.

.. note::
    For performance reasons, :file:`Tensors` are never output during
    a full-dynamic optimization.

To set up a full-dynamic optimization in |calc|, set

::

   RUN = 6

in the :ref:`PARAMETERS` file, and specify which quantity to vary with
the :ref:`OPTIMIZE` parameter. Options to control the behaviour of the
optimization are also available via :ref:`OPTIMIZE`.

The algorithm for full-dynamic optimization calculates |R factor|\ s for at
least three values of the quantity under variation, then fits a parabola. If
the minimum of the parabola is outside of the current scope (i.e., the points
already known), then more points are added to expand the scope such that at
least one point at either side of the minimum is present. Otherwise, the next
calculation is performed at the parabola minimum. The parabola fit is repeated
after obtaining any new data point. Convergence is reached when the new
predicted minimum is within a given distance (defined via :ref:`OPTIMIZE`)
of a point that was already calculated.

.. note::
    The full-dynamic optimization can be slow, as an entire :ref:`ref-calc`
    is executed for each of the values of the quantity under variation.

|R-factor| values corresponding to the calculated values of the quantity under
variation are output to :ref:`fdoptimizationdata`. :ref:`fdoptimizationbeams`
plots the corresponding |IV| curves together with the experimental data.
