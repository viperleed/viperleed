.. include:: /substitutions.rst

.. _r-factor_calculation:

==============
The |R factor|
==============

The reliability factor (|R factor|) is a measure for the deviation between
two |IV| curves or two sets of |IV| curves.  Structure optimization
(\ :ref:`search<sec_search>` or :ref:`fdoptimization` sections) minimizes
the |R factor| between the calculated and experimental |IV| curves.

As the comparison of two curves is not an unambiguous task, multiple |R factor|
implementations exist. The :ref:`RFACTORTYPE` parameter can be set in
the :ref:`PARAMETERS` file to choose between
the smooth R factor |RS| :cite:p:`imreRFactor2026`,
Pendry's |R factor| |RP| :cite:p:`pendryReliabilityFactorsLEED1980`,
:math:`R_2` :cite:p:`vanhoveLowEnergyElectronDiffraction1986,spornAccuracyQuantitativeLEED1998`,
and the Zanazzi-Jona |R factor| :math:`R_\mathrm{ZJ}`
:cite:p:`mooreZanazziJonaPendryReliability1982`:

-   Both |RS| and |RP| compare the logarithmic derivatives of the |IV| data, using a fix
    for the divergence of the logarithm when the intensity approaches zero.
-   :math:`R_2` is based on the mean square difference of the |IV| curves
    (after appropriate scaling).
-   :math:`R_\mathrm{ZJ}` compares the first and second derivatives of the
    |IV| curves (after appropriate scaling).

.. note::
    Using |RP| is the **default setting**.
    For structure optimization, we recommend using |RS|, since it is better suited
    as a target function for minimization than |RP| (|RP| is noisy and
    has many local minima in the vicinity of the global minimum).
    We recommend **not using** :math:`R_2`, since it leads to less accurate results than
    minimization of |RS| or |RP| :cite:p:`spornAccuracyQuantitativeLEED1998,imreRFactor2026`.

.. _pendry_r:

The smooth |R factor| |RS| and Pendry's |R factor| |RP|
-------------------------------------------------------

Both |RS| and |RP| can have values between 0 and 2 and
are defined as
:cite:p:`pendryReliabilityFactorsLEED1980,heinzElectronBasedMethods2013,imreRFactor2026`

.. math::
    :label: eq_RP_def

    R = \frac{\displaystyle\sum\nolimits_\mathbf{g}{\int{\mkern-5mu\left(Y^\mathrm{theo}_\mathbf{g}(E) - Y^\mathrm{exp}_\mathbf{g}(E) \right)^2 dE }}}{\displaystyle\sum\nolimits_\mathbf{g}{\int{\mkern-5mu\left(Y^\mathrm{theo}_\mathbf{g}(E)^2 + Y^\mathrm{exp}_\mathbf{g}(E)^2\right) dE}}},

where :math:`\mathbf{g}` indexes the beams for which the |R factor| is
calculated, and :math:`Y(E)` is the :math:`Y` function. The |R factor|
can be calculated for all beams together or for each beam individually.\ [1]_
The :math:`Y` functions in Eq. :eq:`eq_RP_def` are based on the logarithmic
derivative :math:`\frac{d}{dE}\big(\ln{I(E)}\big) = \frac{I'(E)}{I(E)}`.
This makes |RS| and |RP| insensitive to differences in the absolute
intensities of the |IV| curves. The largest contributions to these |R factors|
come from differences in the *positions* of extrema (minima and maxima).

At deep minima of :math:`I(E)`, the logarithm :math:`\ln I(E)` diverges
towards :math:`-\infty` with :math:`I \to 0`. The logarithmic derivative also
diverges. Thus, :math:`\frac{d}{dE}\big(\ln{I(E)}\big)` cannot be directly used
as a :math:`Y` function in Eq. :eq:`eq_RP_def`. The |R factors| |RS| and |RP|
differ in their strategy how to avoid this divergence (for details, see
:cite:p:`pendryReliabilityFactorsLEED1980` and :cite:p:`imreRFactor2026`);
otherwise they are very similar, and their numeric values are also similar.
Typically, the value of |RS| for a given data set is marginally lower than
that of |RP| (by about 0.01). |RS| was designed to keep all desirable properties
of |RP|, but avoid its shortcomings, especially the noisiness of |RP|, which is
a problem for optimization :cite:p:`imreRFactor2026`.

An |RS| or |RP| value of zero corresponds to perfect agreement between curves.
|RS| and |RP| are close to one for statistically uncorrelated data, while values
larger than one indicate anticorrelation. For close-packed surfaces, |RS| or
|RP| values larger than 0.2 indicate a problem, such as an incorrect structural
model or structural distortions affecting coordinates that were not varied in
the search. |R factors| for more open, corrugated surfaces, such as
missing-row–reconstructed Pt(110), may be around 0.2. |RS| or |RP| values larger
than 0.25–0.30 should be taken as an indication of poor correspondence between
calculated and experimental beams. The best values of |RS| or |RP| obtained by
the `Erlangen group <https://www.fkp.physik.nat.fau.eu/research-schneider/>`__
are below 0.05.\ [2]_

The exact values of |RS| and |RP| slightly depend on the implementation,
specifically on how the derivatives are evaluated numerically. Therefore, the
|R factors| calculated by different implementations in |calc| (see
:ref:`rfactorlegacy`), as well as by :ref:`imagej_plugins`, may differ by about
0.01.

Note that some smoothing algorithms applied to both experimental and
calculated beams, such as the one suggested by
:cite:t:`pendryReliabilityFactorsLEED1980`, artificially reduce the |R factor|,
because they effectively raise the minima of the |IV| curves. At minima, where
the intensities approach zero, |RP| is extremely sensitive to small
differences; artificially increasing the intensity at minima thus gives
smaller |RP| values. Some LEED programs apply such a smoothing; in those
cases, smaller |R factors| than those obtained with ViPErLEED will be
reported, but this does not indicate a better agreement between calculated
and experimental |IV| curves.\ [3]_
By default, ViPErLEED applies no additional smoothing when calculating the
|R factor| (e.g., during structure optimization). Thus, the :ref:`EXPBEAMS`
file should already contain smoothed data.

.. tip::
    We highly recommend to smooth experimental data beforehand using the
    |IV|-curve editor of the :ref:`imagej_plugins`.
    Using the :ref:`RFACTORSMOOTH` parameter for smoothing the experimental
    |IV| curves is discouraged, as the smoothing algorithm applied there is
    inferior to that used by the |IV|-curve editor.

.. _r2_zj_r:

The |R factors| :math:`R_2` and :math:`R_\mathrm{ZJ}`
-----------------------------------------------------

For each beam, :math:`R_2` is defined as
:cite:p:`vanhoveLowEnergyElectronDiffraction1986`

.. math::
    :label: eq_R2_def

    R_2 = \frac{\int{\mkern-5mu\left(I^\mathrm{exp}(E)
                - c\,I^\mathrm{theo}(E)\right)^2 dE}}
               {\int{I^\mathrm{exp}(E)^2\,dE}},
    \qquad
    c = \frac{\int{I^\mathrm{exp}(E)\,dE}}{\int{I^\mathrm{theo}(E)\,dE}},

where the integrals run over the energy range in which both curves are
defined. The scaling factor :math:`c` removes the difference in the
absolute intensities of the two curves. Note that the :math:`R_2` calculated
by the TensErLEED |R-factor| program (used with
:ref:`R_FACTOR_LEGACY<rfactorlegacy>` = ``True``) is half of the value of
Eq. :eq:`eq_R2_def`.

:math:`R_\mathrm{ZJ}` compares the first and second derivatives of the
scaled curves, which makes it sensitive mostly to the positions and shapes
of the extrema :cite:p:`mooreZanazziJonaPendryReliability1982`. ViPErLEED
reports the reduced :math:`R_\mathrm{ZJ}`, i.e., the value divided by
0.027, as proposed by Zanazzi and Jona.

For both :math:`R_2` and :math:`R_\mathrm{ZJ}`, the value for a set of
beams is the average of the values of the individual beams, each weighted
with the energy range in which the experimental and calculated curves of
that beam overlap.

.. [1] Notice that the Pendry or Smooth |R factor| between two sets of beams is
       not the average of the |R factors| between beam pairs, as sums over all
       beams enter both the numerator and the denominator in 
       Eq. :eq:`eq_RP_def`.
       A weighted average (with the energy span of each beam used as its weight)
       is usually close to the overall |R factor|. In contrast, the overall
       :math:`R_2` and :math:`R_\mathrm{ZJ}` (see :ref:`below<r2_zj_r>`) are
       exactly such weighted averages.
.. [2] Unpublished data by Lutz Hammer and coworkers.
.. [3] See also `this discussion <https://github.com/viperleed/viperleed/issues/218>`_\ .
