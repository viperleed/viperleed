.. include:: /substitutions.rst

.. _parabola_fit:

PARABOLA_FIT
============

.. warning::
    .. deprecated:: 0.11.0
       This functionality was experimental.
       Parabola fits don't appear to be very realiable at the moment and may be
       removed or significantly reworked in the future. Use at your own risk.
       Perform an :ref:`error calculation<error_calculation>` for more reliable
       1D |R-factor| data.

PARABOLA_FIT allows fitting the |R-factor| data over the N-dimensional space
of fit parameter values with a paraboloid.

**Default**: PARABOLA_FIT = off

**Syntax**:

::

   PARABOLA_FIT = type linearregression       ! set the regression method to linear regression
   PARABOLA_FIT = type lasso, alpha 1e-5      ! set the regression method to Lasso, with weight 1e-5 for the penalty function
   PARABOLA_FIT = type linear, localize 0.25  ! set the regression method to linear regression, use only the data points within in 1/4 of the displacement ranges, near the best known configuration.
   PARABOLA_FIT = off                         ! turn the parabola fit off entirely (may improve search performance)

**Acceptable values**: see below for each flag.

Generally, input is expected in the form ``flag value``, where the different
flags and allowed values are described below. Values for the different flags
can be set in one line with comma separation, as in the examples above.
Alternatively, you can set different flags by having multiple lines
``PARABOLA_FIT = flag value``.

The paraboloid model fits weights for all the independent parameters ``x``,
the ``x^2``, and the mixed terms ``x_i*x_j``, as well as one bias term.

type
----

The regression method to be used.

**Acceptable values** (not case sensitive): ``linear`` / ``linearregression``,
``ridge``, ``lasso``, ``elasticnet``

Regression methods are imported from the
`sklearn.linear_model module <https://scikit-learn.org/stable/modules/classes.html#module-sklearn.linear_model>`__.
A brief description is given here, but you can find more detailed information
on the documentation pages:
`LinearRegression <https://scikit-learn.org/stable/modules/generated/sklearn.linear_model.LinearRegression.html#sklearn.linear_model.LinearRegression>`__,
`Ridge <https://scikit-learn.org/stable/modules/generated/sklearn.linear_model.Ridge.html#sklearn.linear_model.Ridge>`__,
`ElasticNet <https://scikit-learn.org/stable/modules/generated/sklearn.linear_model.ElasticNet.html#sklearn.linear_model.ElasticNet>`__,
`Lasso <https://scikit-learn.org/stable/modules/generated/sklearn.linear_model.Lasso.html#sklearn.linear_model.Lasso>`__,
as well as this `overview page <https://scikit-learn.org/stable/modules/linear_model.html>`__

-  **Linear regression:** Unbiased least-squares fit.
-  **Ridge:** Least-squares fit with a penalty term on the sum over the squares
   of weights ``w1_i`` and ``w2_i``. This helps prevent overfitting, but does
   not artificially sparsen the number of used parameters.
-  **Lasso:** Least-squares fit with a penalty term on the sum over the
   absolute values of the weights ``w1_i`` and ``w2_i``. This tends to
   select a solution with as few non-zero weights as possible.
-  **Elastic Net:** Linear combination of the Lasso and Ridge methods.

Note that for all methods except ordinary linear regression, the penalty term
affects not only the curvature of the parabola, but also the position of the
minimum, as this is given by the ``w1_i * x_i`` term. Fitting is performed with
the parameter values centered around the current best configuration (i.e., the
combination yielding the lowest |R-factor| value in the search so far), so 
Ridge, Lasso and Elastic Net will all favour solutions close to the best 
known configuration.

alpha
-----

**Acceptable values**: positive float.

The prefactor for the penalty term of the ridge, lasso, and elastic net methods
(see above). If linear regression is used, alpha will remain unused. Note that
all methods listed above become ordinary linear regression for ``alpha = 0``.


..
   This section is commented out for now, because the feature is unused.

   localize
   --------

   **CURRENTLY NOT ACTIVE - best way to do something like this needs to**
   **be discussed.** Currently, the RR value ``RR = 8 * V0i / enrange`` is
   calculated, where V0i is the imaginary part of the inner potential and
   enrange the total energy range of all beams. Points farther than 3*RR
   from the best known |R factor| are discarded. Maybe a reasonable 'localize'
   parameter would be to re-define this prefactor to RR, i.e., influence the
   |R-factor| cutoff.

   **Acceptable values**: (0, 1), or 0 to deactivate.

   Limits which points on the |R-factor| landscape should be used for the fit.
   An interval centered around the best known configuration is placed on each
   displacements range, and only configurations that fall into this interval
   *for every parameter* are used for the fit. The value defines the fraction
   of the displacement range that should be used. For example, if you set
   ``localize`` to 0.25, then one quarter of the parameter space for each
   parameter, or (1/4^N) of the parameter space (for N independent parameters),
   will be used.

   If the |R-factor| landscape is rough, the ``localize`` flag can prevent
   points far from the minimum to affect the fit by only using points close
   to the global minimum. However, this requires that the best known
   configuration (identified by the search) is already close to the
   true global minimum.
