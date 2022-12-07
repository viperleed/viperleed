.. _rmut:

GAUSSIAN_WIDTH
==============

GAUSSIAN_WIDTH (corresponds to RMUT parameter in fortran) is a control parameter for the probability distribution of step width during the search. A larger GAUSSIAN_WIDTH corresponds to a higher probability of taking a large step. However, the width of the distribution function is not only controlled by GAUSSIAN_WIDTH, but also has contributions from the number of independent parameters and a dynamic parameter that is varied during the search.

**Default:** 0.1

**Allowed values:** positive real

**Syntax:**

::

   GAUSSIAN_WIDTH = 0.05

**GAUSSIAN_WIDTH is one of the main parameters controlling convergence**, and should therefore not be left at the default, but rather adapted as needed. For rough (wide scan) searches, a value of 0.1 is appropriate. The finer the fit, the smaller the steps should be, so RMUT should then be set at 0.01 or even lower.

**TODO: Could use some more advice-style info from Lutz**
