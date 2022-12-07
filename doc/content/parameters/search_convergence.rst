.. _search_convergence:

SEARCH_CONVERGENCE
==================

SEARCH_CONVERGENCE defines a number of parameters controlling the behaviour of the Tensor LEED search. 
See Ref. :cite:t:`kottckeNewApproachAutomated1997` for an explanation of the search algorithm used by TensErLEED.

Most importantly, :ref:`GAUSSIAN_WIDTH<rmut>` (defined with flag ``gaussian``, corresponds to RMUT parameter in TenErLEED code) controls the step width during the search.
SEARCH_CONVERGENCE also defines *partial convergence* criteria, at which the search will be stopped to re-scale :ref:`GAUSSIAN_WIDTH<rmut>`, and the re-scaling factors applied at that point.

**Syntax:**

::

   ! general case:
   SEARCH_CONVERGENCE flag = value [value]   ! see specifics for the different flags below
   ! example:
   SEARCH_CONVERGENCE gaussian = 0.3 0.2
   SEARCH_CONVERGENCE dgen = 1000 1.5

   ! to turn automatic convergence off:
   SEARCH_CONVERGENCE = off
   SEARCH_CONVERGENCE gaussian = 0.01    ! if automatic convergence is turned off, the GAUSSIAN_WIDTH parameter should be set manually
   ! OR:
   SEARCH_CONVERGENCE gaussian = 0.01 1  ! setting scaling factor 1 turns automatic convergence off implicitly

The specific flags are:

gaussian
--------

The flag ``gaussian`` controls :ref:`GAUSSIAN_WIDTH<rmut>`, a control parameter for the 
probability distribution of step width during the search.
A larger GAUSSIAN_WIDTH corresponds to a higher probability of taking a 
large step.
However, the width of the distribution function is not only controlled 
by GAUSSIAN_WIDTH, but also has contributions from the number of independent 
parameters and a dynamic parameter that is varied during the search.
Two values are accepted: The value to be used at the beginning of the search, 
and a scaling factor applied whenever partial convergence is reached.

**Default:** 0.5 0.5

**Allowed values:** one positive real, one real in range ]0,1[

**Syntax:**

::

   SEARCH_CONVERGENCE gaussian = 0.1 0.3  ! initialize with 0.1, scale with factor 0.3 whenever partial convergence is reached

**GAUSSIAN_WIDTH is one of the main parameters controlling convergence**, and appropriate values are highly dependent on the system in use. 
By default, :ref:`GAUSSIAN_WIDTH<rmut>` is initialized with a very high value of 0.5, and downscaled by 50% every time partial convergence is reached (see ``dgen`` flag below).

dgen
----

The flag ``dgen`` defines how many generations should pass without any changes 
to either the best configuration, or all configurations in the population, 
before the search is considered converged *for the current value of GAUSSIAN_WIDTH*. 
Once this condition is met, both :ref:`GAUSSIAN_WIDTH<rmut>` and the ``dgen`` values themselves 
will be scaled by their respective scaling factors, and the search will be continued. 
**Full convergence** of the search is considered to be reached when **partial convergence** 
is reached twice in a row without any changes in the meantime, i.e. when lowering 
GAUSSIAN_WIDTH does not lead to any further improvement. Once full convergence 
is reached, the search will be stopped. Note that an upper limit of generations, 
independent of convergence, is given by :ref:`SEARCH_MAX_GEN<SEARCHGENMAX>`.

**Default:** partial convergence when best 10% of structures do not change in 
100 generations, scaling factor = 1/(scaling factor of GAUSSIAN_WIDTH)

**Allowed values:** one or two values > 1

**Syntax:**

::

   SEARCH_CONVERGENCE dgen dec = 1000 1.5   ! best 10% of structures don't change in 1000 generations, then for lowered GAUSSIAN_WIDTH in 1500 generations, then 2250, etc.
   SEARCH_CONVERGENCE dgen = 1000 1.5       ! same as previous, "dgen" defaults to "dgen dec"
   SEARCH_CONVERGENCE dgen = 1000           ! like the previuos one, but keep default scaling factor
   SEARCH_CONVERGENCE dgen best = 1000 1.5  ! best structure doesn't change in 1000 generations, then for lowered GAUSSIAN_WIDTH in 1500 generations, then 2250, etc.
   SEARCH_CONVERGENCE dgen all = 200 2      ! the entire population doesn't change in 200 generations, then for lowered GAUSSIAN_WIDTH in 400 generations, etc.

The additional flags "dec", "best" and "all" can be used to specify whether all configurations in the population, the best 10%, or only the best configuration should be considered. If no additional flag is given, ``SEARCH_CONVERGENCE dgen`` will default to checking the best 10% of the population. The scaling factor will default to the inverse of the scaling factor used by :ref:`GAUSSIAN_WIDTH<rmut>` for any of the three. Defining values for more than one convergence criterion is allowed; in that case, partial convergence will be considered to have been reached once *either* condition is met, but full convergence is reached only once *all* conditions are met.
