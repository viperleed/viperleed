.. _search_cull:

SEARCH_CULL
===========

SEARCH_CULL defines a portion of the :ref:`population<SEARCHPOP>`  that will be dropped every time :ref:`partial convergence<SEARCH_CONVERGENCE>`  is reached in the search. The given fraction of the population (or number of members) with the worst R-factors will be removed, and replaced with new configurations generated from the surviving population, or at random. The mode of generating the new configurations can also be defined.

See Ref. `[1 <SEARCH_CULL#ref1>`__] for an explanation of the search algorithm used by TensErLEED.

**Default:** 0.1 genetic

**Allowed values:** Positive real < 1 (percentage of population), or positive integer (absolute number to be removed); optional: flag ``clone`` / ``genetic`` / ``random`` to define how to create new configurations

**Syntax:**

::

   SEARCH_CULL = 0.2 genetic   ! cull the worst-performing 20% of the population, replace by offspring from two random survivors
   SEARCH_CULL = 0.2           ! equivalent to the above, genetic is default
   SEARCH_CULL = 0.2 clone     ! cull the worst-performing 20% of the population, replace by clones of survivors (picked at random)
   SEARCH_CULL = 3             ! cull exactly three worst-performing configurations
   SEARCH_CULL = 2 random      ! cull exactly two worst-performing configurations, replace them with new configurations re-initialized at random.

For most systems, part of the population will get stuck in local minima. Culling some of those is generally safe, and may improve performance by increasing the population density in better-performing parts of parameter space.

Defining SEARCH_CULL via a percentage of the overall population is recommended, as this will be more robust to later changes to the population size.

Note that if the offsprings are generated as clones, this automatically improves the average r-factor of the population. If the offsprings are re-initialized as random configurations or via the genetic algorithm, the average r-factor may actually get worse.

For the ``genetic`` mode, each new configuration is generated as follows: Two parents are picked at random from the surviving population (such that the parents are not identical, if possible). Then, for each independent parameter, the new configuration inherits the value from one of the parents at random. Dependent parameters are then filled up in keeping with the constraints of the system. Note that this does not exclude the possibility of the offspring being identical to one of the parents, especially if the parents are similar to begin with.

--------------

M. Kottcke and K. Heinz, *A New Approach to Automated Structure Optimization in LEED Intensity Analysis* `Surf. Sci. 376, 352 (1997) <http://dx.doi.org/10.1016/S0039-6028(96)01307-6>`__.
