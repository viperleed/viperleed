.. _searchpop:

SEARCH_POPULATION
=================

SEARCH_POPULATION defines the number of trial structures used in parallel
during the search. The trial structures from this population can be mixed
in the genetic algorithm.

See Ref. :cite:p:`kottckeNewApproachAutomated1997` for an explanation of the
search algorithm used by TensErLEED.

**Default:** min(48, 15 + number of independent search parameters), rounded up
to a multiple of the number of available cores.

**Allowed values:** Positive integer

**Syntax:**

::

   SEARCH_POPULATION = 24

For best performance, SEARCH_POPULATION should be an integer multiple of the
number of available cores :ref:`N_CORES<NCORES>`.

Since each trial structure evolves by randomly modifying its parameters,
there is always a risk of some structures getting stuck in a local minimum,
especially when :ref:`SEARCH_CONVERGENCE gaussian<rmut>` is small.
Large populations reduce this risk, and it can therefore be useful to scale
up the population to deal with large parameter spaces. Note, however, that
every independent parameter adds an extra dimension to the parameter space,
and computational cost scales (roughly) linearly with the population size.
For a large number of independent parameters, scaling the population up
in direct proportion to the size of the parameter space is therefore not
realistic, so unreasonably large parameter spaces should be avoided.
