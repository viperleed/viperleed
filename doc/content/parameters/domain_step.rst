.. _domain_step:

DOMAIN_STEP
===========

DOMAIN_STEP defines the area fraction step (in %) for the domain search,
i.e. the grid on which the relative weight of domains is being varied.
See also the page about :ref:`Domain calculations<domain_calculation>`.

**Default:** 1

**Allowed values:** Integer between 1 and 100, must divide 100

**Syntax:**

::

   DOMAIN_STEP = 5

When fitting multiple domains, one weight parameter per domain is added to the
search. These parameters are varied freely, but renormalized to sum to 100%,
resulting in an effective N-1 free parameters, with N the number of domains.
DOMAIN_STEP defines the step width for these parameters, so increasing the
value may speed up the search.

DOMAIN_STEP should always be a divisor of 100 because the variation is
performed with integer-valued parameters, with the lowest value (1) always
corresponding to 0% and the highest value corresponding to 100%.
