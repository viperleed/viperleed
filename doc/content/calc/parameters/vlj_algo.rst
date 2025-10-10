.. _vlj_algo:

VLJ_ALGO
========

VLJ_ALGO sets the algorithm(s) used during the :ref:`tensor-LEED<tensor_leed>`
optimization using the viperleed-jax plugin. When an algorithm is specified as
a flag, that parameters settings can specified instead.

**Default:** ``CMAES, SLSQP``


**Syntax:**

::

   ! setting the algorithm to use
   VLJ_ALGO = CMAES          ! use only the CMA-ES algorithm
   VLJ_ALGO = CMAES, BFGS    ! use a two-step optimization with CMA-ES and BFGS
   VLJ_ALGO = SLSQP          ! use only the SLSQP algorithm

   ! settings for the algorithms
   VLJ_ALGO CMAES = pop 50, gens 100, ftol 1e-4
   VLJ_ALGO SLSQP = grad False


When using VLJ_ALGO without a flag, one or two algorithms can be specified on
the right-hand side, separated by commas. If two algorithms are specified,
a two-phase optimization will be performed where the second algorithm starts
from the best result of the first phase.
This can be useful to combine the global search capabilities of CMA-ES with the
local refinement capabilities of gradient-based algorithms such as
BFGS or SLSQP.

On a separate line, any of the algorithms can be specified as a flag to specify
settings to be used for the optimization using that algorithm. See below for
the allowed flags and their values for each algorithm.

:ref:`VLJ_CONFIG` `precondition`

CMAES
-----

The CMA-ES (Covariance Matrix Adaptation Evolution Strategy) algorithm is a
global optimization algorithm that is well-suited for high-dimensional and
noisy optimization problems. viperleed-jax uses the implementation of CMA-ES
from the :cite:t:`Clinamen2` library.

**Default:** ``pop 30, gens 200, ftol 1e-3``

**Syntax:**

::

   ! settings for the algorithms
   VLJ_ALGO CMAES = pop 50, gens 100, ftol 1e-4
   VLJ_ALGO CMAES = gens 400, ftol 1e-5


Available setting for CMAES are:

-  ``pop``: Population size, i.e. number of samples per generation.

-  ``gens``: Maximum number of generations to run the algorithm. If the
   algorithm converges before reaching this number, it will stop early.

-  ``ftol``: The convergence criterion on the R-factor value. The algorithm
   will stop if the R-factor does not improve by more than this value over
   the last generation.


SLSQP
-----

The SLSQP (Sequential Least Squares Programming) algorithm is a gradient-based
optimization algorithm that is well-suited for constrained optimization
problems. viperleed-jax uses the implementation of SLSQP from the
`SciPy <https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.minimize.html#scipy.optimize.minimize>`__
library.

**Default:** ``grad True, grad_damping 0.1``

Available setting for ``SLSQP`` are:

-  ``grad``: Toggle calculation of gradients using automatic differentiation.

-  ``grad_damping``: A factor applied to the computed gradient which can help to
   avoid the algorithm jumping outside of the basin around the global minimum.


BFGS
----

The BFGS (Broyden-Fletcher-Goldfarb-Shanno) algorithm is a different gradient-based
optimization algorithm well suited to high-dimensional problems. Unlike `SLSQP`
it always evaluates the gradient at every step. viperleed-jax uses the
implementation of BFGS from the
`SciPy <https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.minimize.html#scipy.optimize.minimize>`__
library.

**Default:** ``grad True``

Available setting for ``BFGS`` are:

-  ``grad``: Toggle calculation of gradients using automatic differentiation.


.. note::
    The ``grad`` switch for the ``SLSQP`` and ``BFGS`` algorithms can have
    significant impact on the performance of the optimization.
    If set to ``False``, the algorithm will use finite differences to
    approximate gradients. If set to ``True``, the algorithm will use
    automatic differentiation to compute gradients, which can be more
    efficient, but may significantly increase memory usage due to the need to
    store intermediate results for gradient calculation.
