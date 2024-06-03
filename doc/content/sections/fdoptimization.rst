.. _fdoptimization:

.. include:: /substitutions.rst

Full-dynamic optimization
=========================

|LEED-IV| calculations require setting some parameters that cannot be
(exactly) determined in experiment, but that also are not accessible
to variation within the Tensor LEED approximation. These include
:ref:`BEAM_INCIDENCE<BEAMINCIDENCE>`, :ref:`V0_IMAG<v0_imag>`
and the exact dimensions of the unit cell. In order to fit these
parameters, multiple full-dynamic (i.e. "reference") calculations have
to be performed. For performance reasons, Tensor output is always disabled
during full-dynamic (FD) optimization.

To set up a full-dynamic optimization, set

::

   RUN = 6

in the :ref:`PARAMETERS<PARAMETERS>`  file, and specify which parameter
to vary with the :ref:`OPTIMIZE<OPTIMIZE>`  parameter. More options to
control the behaviour of the optimization can also be set using
:ref:`OPTIMIZE<OPTIMIZE>`.

The algorithm for full-dynamic optimization is essentially to calculate
R-factors for at least three values of the parameter under variation,
then fit a parabola. If the minimum of the parabola is outside of the
current scope (i.e. the already known points), then more points are
added to expand the scope such that at least one point at either side
of the minimum is present. Otherwise, the next calculation will be performed
at the parabola minimum. The parabola fit is repeated after obtaining any new
data point. Convergence is reached when the new predicted minimum is within a
given distance (defined by :ref:`OPTIMIZE<OPTIMIZE>`) of a point that was
already calculated.

R-factor values corresponding to the calculated values of the parameter
under variation are output to
:ref:`FD_Optimization.csv and FD_Optimization.pdf<fdoptimizationdata>`.
I(V) curves are output to :ref:`FD_Optimization_beams.pdf<fdoptimizationbeams>`.
