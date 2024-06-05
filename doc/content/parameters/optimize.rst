.. include:: /substitutions.rst

.. _optimize:

OPTIMIZE
========

The OPTIMIZE parameter defines the behaviour of
:ref:`full-dynamic optimization<Fdoptimization>` runs.

**Syntax**:

::

   OPTIMIZE V0i = 0.5                                ! optimize the imaginary part of the inner potential, with an initial step of 0.5 eV.
   OPTIMIZE V0i = step 0.5                           ! same as above
   OPTIMIZE V0i = step -0.5                          ! same as above, but take the initial step towards smaller V0i rather than larger
   OPTIMIZE theta = step 1, convergence 0.05         ! optimize beam incidence polar angle theta with initial step of 1 degree, define convergence at 0.05 degrees
   OPTIMIZE phi = step 5, minpoints 5, maxpoints 8   ! optimize beam incidence azimuthal angle phi with intial step of 5 degrees, calculating at least 5 data points and at most 8 data points.
   OPTIMIZE ab = 0.02                                ! optimize the unit cell dimensions ab (scaling together), with an initial step of 2% scaling

**Acceptable values**: see below for each flag.

The left-hand side defines which parameter to optimize. Allowed values (not
case sensitive) are ``v0i`` for :ref:`V0_IMAG<v0_imag>`, ``theta`` and ``phi``
for the :ref:`BEAM_INCIDENCE<BEAMINCIDENCE>`  angles, ``a``, ``b`` and ``c``
for the unit cell vectors, as well as ``ab`` and ``abc`` for simultaneous,
proportional scaling of multiple unit cell vectors.

Input on the right-hand side is expected in the form ``flag value``, where
the different flags and allowed values are described below. Values for the
different flags can be set in one line with comma separation, as in the
examples above. Alternatively, you can set different flags by having multiple
lines ``OPTIMIZE x = flag value``. If only one value is given, this is assigned
as ``step``.

The minimum form is to define only what parameter to optimize and the initial
step width.

Values are always given relative to the input configuration, using

-  ``V0i``: eV
-  ``theta``, ``phi``: degrees
-  unit cell: values are scaling factors. For example, a ``step`` of 0.02
   corresponds to an initial data point scaling the unit vectors by a factor
   1.02, i.e. +2%.

step
----

The initial step width, and direction of the first step.

**Acceptable values:** non-zero float

The first data point is always calculated at the input configuration.
The second data point is defined entirely by the value of ``step``,
including the sign. For example, if the parameter under variation is
``V0i`` and the step is -0.5, then the second data point will be
calculated at ``V0i_initial - 0.5 eV``. After the initial step, the
sign is not used any more; the third data point will be calculated
at the same step width from either the first or second data point,
going in the direction of the minimum. Afterwards, a parabola is fit;
if the minimum is inside the range of known data point at that point,
then ``step`` is not used further. If it is outside, ``step`` will still
be used in the choice of further data points closer to the minimum.

convergence
-----------

Cutoff for stopping optimization.

**Acceptable values:** float

**Default:** 0.1 \* ``step``

During optimization, a new parabola is fit after every new data point.
If the new minimum is close to a point that was already calculated,
optimization will stop. The cutoff for this is given by the ``convergence``
value.

Note that ``convergence`` is given in units corresponding to the parameter,
e.g., a convergence 0.1 for V0i means convergence is reached when the minimum
is within 0.1 eV of a known point.

minpoints, maxpoints
--------------------

The minimum and maximum number of data points to calculate.

**Default:** 4 and 10

**Acceptable values:** positive integer

``minpoints`` defines a minimum number of data points. If less data points
have been calculated, the ``convergence`` criterion is ignored. This can be
used to ensure that the calculation does not stop with very few data points,
which may make the parabola fit unreliable. However, if the new parabola
minimum would fall within ``convergence`` of an already known data point,
a new point other than the minimum is chosen; this is meant to ensure that
new data points actually contribute new data to the parabola fit, rather
than re-calculating known data.

``maxpoints`` defines a maximum number of data points, after which
optimization will stop regardless of convergence. This may happen
if the |R factor| as a function of the parameter under variation
is not parabolic, or if ``convergence`` is set too low. Can also
be used to prevent excess computations.

maxstep
-------

Defines a maximum step width.

**Default:** 3 \* ``step``

**Acceptable values:** positive float

If the initial three data points do not describe a convex parabola, further
data points will be added outside the current scope with increasing step width.
This is meant to ensure that a better parabola is found if the initial data
points are noisy due to ``step`` having been chosen too small. ``maxstep``
defines the maximum step width to which the ``step`` will be extended.
