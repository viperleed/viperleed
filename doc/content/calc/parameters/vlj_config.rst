.. _vlj_config:

VLJ_CONFIG
==========

VLJ_CONFIG allows to set several configuration parameters for the
:ref:`tensor-LEED<tensor_leed>` optimization using the viperleed-jax plugin.
These parameters control global settings around the calculation of the amplitude
changes and the optimization process.

**Syntax:**

::

   VLJ_CONFIG = recalc_ref_t_matrices True  ! use re-calculated t-matrices
   VLJ_CONFIG = precondition False          ! do not use preconditioning


Multiple parameters can be set on the right-hand side, separated by commas.

.. todo:
   Cite and refer to the SI of the viperleed-jax paper for more details on the
   available settings and their effects on the optimization process once
   published.

Available setting are:

-  ``precondition``: When using a two-phase optimization strategy (e.g. CMA-ES
   followed by SLSQP, see :ref:`vlj_algo`), this flag controls whether the
   :math:`R`-factor surface for the second phase is preconditioned using the
   covariance matrix estimated by CMA-ES.
   This is based on the idea that the covariance matrix in CMA-ES in the limit
   of many generations approaches the inverse of the local Hessian matrix of the
   sampled function :cite:p:`shirCovarianceHessianRelationEvolution2020`.

   Preconditioning can help improve the convergence of the second phase
   optimization.
   **Default**: True.

-  ``recalc_ref_t_matrices``: Toggle whether to re-calculate the t-matrices for
   the reference structure or to use the t-matrices as read in from the tensor
   files (calculated in the :ref:`ref-calc`).
   The t-matrix calculation in viperleed-jax is uses different (more efficient
   and accurate) algorithms than the TensErLEED reference calculation. This
   leads to slightly different results for the t-matrices of the reference
   structure. Since the calculation of amplitude changes uses the differences
   between the t-matrices of the reference structure and the perturbed
   structure, this can lead to slightly different results for the amplitude and
   thus the :math:`R`-factor.
   **Default**: False.

-  ``t-leed-l_max``: Maximum value of the angular momentum quantum number
   :math:`\ell` to be used in the tensor-LEED calculation. 
   **Default**: Use the maximum value of :ref:`LMAX`.

-  ``use_symmetry`` (experimental): Toggle whether to use the advanced symmetry
   based calculation of t-matrices and propagtors. This can significantly
   improve performance.
   **Default**: True.
