.. _search_convergence:

SEARCH_CONVERGENCE
==================

SEARCH_CONVERGENCE defines a number of parameters controlling the behaviour of the Tensor LEED search. 
See Ref. :cite:t:`kottckeNewApproachAutomated1997` for an explanation of the search algorithm used by TensErLEED.

Most importantly, the SEARCH_CONVERGENCE flag ``gaussian`` (corresponds to RMUT parameter in TenErLEED code) controls the step width during the search.
SEARCH_CONVERGENCE also defines *partial convergence* criteria, at which the search will be interrupted to re-scale the width of the gaussian distribution, and the re-scaling factors applied at that point.

**Syntax:**

::

   ! general case:
   SEARCH_CONVERGENCE flag = value [value]   ! see specifics for other flags below
   ! example:
   SEARCH_CONVERGENCE gaussian = 0.3 0.2
   SEARCH_CONVERGENCE dgen = 1000 1.5

   ! to turn automatic convergence off:
   SEARCH_CONVERGENCE = off
   ! if automatic convergence is turned off, the gaussian flag should be set manually
   SEARCH_CONVERGENCE gaussian = 0.01
   ! OR:
   SEARCH_CONVERGENCE gaussian = 0.01 1
   ! setting scaling factor 1 turns automatic convergence off implicitly

The specific flags are:

.. _rmut:

gaussian
--------

The flag ``gaussian`` is a scaling factor applied to the standard deviation of the
probability distribution (a normal distribution) of the step width during the search (see the TensErLEED paper by :cite:t:`blumFastLEEDIntensity2001a`).
It corresponds to the ``RMUT`` parameter in the TensErLEED code.

A larger value corresponds to a higher probability of taking a 
large step.
However, the width of the distribution function is not only controlled 
by the ``gaussian`` flag, but also has contributions from the number of independent 
parameters and a dynamic parameter that is varied during the search.
Two real values are accepted: The value to be used at the beginning of the search, 
and a scaling factor applied whenever partial convergence is reached.

**Default:** 0.1 0.5

**Allowed values:** ]0,∞[, ]0,1[

**Syntax:**

::

   ! initialize with=0.1, scale factor=0.3 whenever partial convergence is reached
   SEARCH_CONVERGENCE gaussian = 0.1 0.3

**SEARCH_CONVERGENCE** ``gaussian`` **is one of the main parameters controlling convergence**, and appropriate values are highly dependent on the system in use. 
By default, ``gaussian`` is initialized with a very high value of 0.5, and downscaled by 50% every time partial convergence is reached (see ``dgen`` flag below).
For rough (wide scan) searches, a value of 0.1 is appropriate. The finer the fit, the smaller the steps should be, so ``gaussian`` may be set at 0.01 or even lower.

**TODO Lutz: Could use some more advice-style info**

dgen
----

The flag ``dgen`` defines how many generations should pass without any changes 
to either the best configuration, or all configurations in the population, 
before the search is considered converged for the current value of ``gaussian``.
Once this condition is met, both ``gaussian`` and the ``dgen`` values themselves 
will be scaled by their respective scaling factors, and the search will be continued.
**Full convergence** of the search is considered to be reached when **partial convergence** 
is reached twice in a row without any changes in the meantime, i.e. when lowering 
``gaussian`` does not lead to any further improvement.
Once full convergence is reached, the search will be stopped.
Note that an upper limit of generations, 
independent of convergence, is given by :ref:`SEARCH_MAX_GEN<SEARCHGENMAX>`.

.. note::
   ``dgen`` also affects the output interval for the raw TenErLEED file :ref:`SD.TL<sdtl>`.
   Very values for ``dgen`` (:math:`\lesssim` 10) may slow down the search due to I/O overhead.

**Default:** partial convergence when best 10% of structures do not change in
200 generations, scaling factor = 0.75/(scaling factor of ``gaussian``).

**Allowed values:** one or two values > 1

**Syntax:**

::

   ! best 10% of structures don't change in 1000 generations, then for lowered gausssian flag in 1500 generations, then 2250, etc.
   SEARCH_CONVERGENCE dgen dec = 1000 1.5

   ! same as previous, "dgen" defaults to "dgen dec"
   SEARCH_CONVERGENCE dgen = 1000 1.5

   ! same as previous, but keep default scaling factor
   SEARCH_CONVERGENCE dgen = 1000

   ! best structure doesn't change in 1000 generations, then for lowered gausssian flag in 1500 generations, then 2250, etc.
   SEARCH_CONVERGENCE dgen best = 1000 1.5

   ! the entire population doesn't change in 200 generations, then for lowered gausssian flag in 400 generations, etc.
   SEARCH_CONVERGENCE dgen all = 200 2

The additional flags "all","dec", and "best" can be used to specify whether all configurations in the population, the best 10%, or only the best configuration should be considered.
If no additional flag is given, ``SEARCH_CONVERGENCE dgen`` will default to checking the best 10% of the population.
The scaling factor will default to the inverse of the scaling factor used by ``gaussian`` for any of the three.

Defining values for more than one convergence criterion is allowed; in that case, partial convergence will be considered to have been reached once *either* condition is met, but full convergence is reached only once *all* conditions are met.


Changelog
---------

.. versionchanged:: 0.12.0
    Changed default values for ``SEARCH_CONVERGENCE gaussian`` from 0.5 0.5 to
    0.1 0.5.
    Change default values for ``SEARCH_CONVERGENCE dgen`` from 100
    1/(scaling factor of ``gaussian``) to 200 0.75/(scaling factor of
    ``gaussian``).
