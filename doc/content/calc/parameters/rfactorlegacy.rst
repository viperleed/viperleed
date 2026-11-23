.. include:: /substitutions.rst

.. _rfactorlegacy:

R_FACTOR_LEGACY
===============

R_FACTOR_LEGACY toggles between using the old (run-time compiled) TensErLEED
|R factor| and a new python-native |R-factor| calculation implemented in
ViPErLEED. The new implementation uses spline interpolation, which yields
more stable (second) derivatives. However, this may result in slightly different
|R factor| values even when using the same :ref:`rfactortype`.

By default, the TensErLEED |R factor| is used. The default will change
to the python-native |R factor| once ViPErLEED supports structure
searches that do not use TensErLEED. R_FACTOR_LEGACY only selects the
|R factor| reported after a :ref:`reference calculation<ref-calc>` or a
:ref:`superposition calculation<super_pos>`, the one used for
:ref:`error curves<error_calculation>`, and the one minimized in the
:ref:`full-dynamic optimization<fdoptimization>`. The
TensErLEED structure search always uses the TensErLEED |R factor|,
whatever the value of R_FACTOR_LEGACY.

The python-native implementation calculates :math:`R_2` with its
literature definition :cite:p:`vanhoveLowEnergyElectronDiffraction1986`
(see |R-factor| :ref:`calculation<r-factor_calculation>`). The
:math:`R_2` of TensErLEED includes an additional factor 1/2: for the
same curves, the TensErLEED :math:`R_2` is half of the python-native
one. This also applies to the :math:`R_2` reported by the TensErLEED
structure search.

With R_FACTOR_LEGACY = ``False``, the |R-factor| calculation differs
from the TensErLEED one in these respects:

-  The |IV| curves are interpolated with cubic splines (not-a-knot
   boundary conditions), and the derivatives needed by the
   :math:`Y` functions are the analytic derivatives of these splines.
   TensErLEED instead uses local polynomial interpolation and
   finite-difference derivatives.
-  Both implementations scan the same inner-potential shifts (the values
   of :ref:`IV_SHIFT_RANGE<ivshiftrange>`, at its step) and compare the
   curves on a grid with the same step.
-  :ref:`R_FACTOR_SMOOTHING<rfactorsmooth>` is ignored: experimental
   curves are not smoothed. A warning is issued if R_FACTOR_SMOOTHING is
   set together with R_FACTOR_LEGACY = ``False``.

The new implementation supports only cubic interpolation. Combining
R_FACTOR_LEGACY = ``False`` with an INTPOL_DEG other than 3 is rejected
when the parameters are read. Set R_FACTOR_LEGACY to ``True`` to use other
interpolation degrees.

.. versionchanged:: 0.15.0
   R_FACTOR_LEGACY = ``False`` uses the python-native |R factor| that is
   part of ViPErLEED, for all four :ref:`R-factor types<rfactortype>`.
   Earlier versions used the python-native |R factor| only if an optional
   compiled extension was installed, and the TensErLEED |R factor|
   otherwise. Its :math:`R_2` is twice the TensErLEED value. INTPOL_DEG
   other than 3 is rejected if R_FACTOR_LEGACY is ``False``.

**Default:** True

**Allowed values:** True (old TensErLEED), False (new)

**Syntax:**

::

   R_FACTOR_LEGACY = True
