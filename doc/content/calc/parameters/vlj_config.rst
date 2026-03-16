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
   **Default**: *True*.

-  ``recalc_ref_t_matrices``: If *True*, re-calculate the t-matrices for
   the reference structure, else use the t-matrices as read in from the tensor
   files (calculated in the :ref:`ref-calc`).
   The t-matrix calculation in viperleed-jax is uses different (more efficient
   and accurate) algorithms than the TensErLEED reference calculation. This
   leads to slightly different results for the t-matrices of the reference
   structure. Since the calculation of amplitude changes uses the differences
   between the t-matrices of the reference structure and the perturbed
   structure, this can lead to slightly different results for the amplitude and
   thus the :math:`R`-factor.
   **Default**: *False*.

-  ``t-leed-l_max``: Cutoff value of the angular momentum quantum number
   :math:`\ell` to be used in the tensor-LEED calculation.
   **Default**: Use the highest value used in the reference calculation (as
   determined by :ref:`LMAX` and :ref:`phaseshiftmin`).

-  ``use_symmetry`` (experimental): Toggle whether to use the advanced symmetry
   based calculation of t-matrices and propagtors. This can significantly
   improve performance.
   **Default**: *True*.

-  ``occ_norm``: Select the method used to normalize occupational parameters
   that ensures the total occupation of all elements on one site is
   :math:`\leq 1`. Available options are:

   - ``mirror``: If the total occupation is larger than 1, mirror the occupation
     parameters on the :math:`\sum_i c_i = 1` (hyper-)plane.
   - ``project``: If the total occupation is larger than 1, project the
     occupation parameters onto the :math:`\sum_i c_i = 1` (hyper-)plane.

   **Default**: ``mirror``.

-  ``preoptimize_v0r``: If *True*, the initial guess for the offset of the real
   part of the inner potential :math:`V_0r` is sampled and optimized on the
   main optimization.
   The :math:`R`-factor is very sensitive to :math:`V_0r`. Pre-optimizing this
   parameter can help find a better starting point for the proceeding
   optimization.

   If the search segment is executed without a reference calculation during the
   same run, and no pre-optimization is done, the :math:`R`-factor at of the
   initial guess may be artificially high, since no information about the
   correct shift of :math:`V_0r` is available.

   **Default**: *True*.
