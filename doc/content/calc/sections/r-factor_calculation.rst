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
Pendry's |R factor| |RP| :cite:p:`pendryReliabilityFactorsLEED1980`
and :math:`R_2` :cite:p:`spornAccuracyQuantitativeLEED1998`:

-   |RP| compares the logarithmic derivatives of the |IV| data, using a fix
    for the divergence of the logarithm when the intensity approaches zero.
-   :math:`R_2` is based on the mean square difference of the |IV| curves
    (after appropriate scaling).

.. note::
    Using |RP| is the **default setting** and highly encouraged since tests
    have shown that it leads to better results than :math:`R_2`
    :cite:p:`spornAccuracyQuantitativeLEED1998`.

.. _pendry_r:

The Pendry |R factor|
---------------------

The Pendry |R factor|, |RP|, can have values between 0 and 2 and
is defined as
:cite:p:`pendryReliabilityFactorsLEED1980,heinzElectronBasedMethods2013`

.. math::
    :label: eq_RP_def

    R_{\mathrm{P}} = \frac{\displaystyle\sum\nolimits_\mathbf{g}{\int{\mkern-5mu\left(Y^\mathrm{theo}_\mathbf{g}(E) - Y^\mathrm{exp}_\mathbf{g}(E) \right)^2 dE }}}{\displaystyle\sum\nolimits_\mathbf{g}{\int{\mkern-5mu\left(Y^\mathrm{theo}_\mathbf{g}(E)^2 + Y^\mathrm{exp}_\mathbf{g}(E)^2\right) dE}}},

where :math:`\mathbf{g}` indexes the beams for which the |R factor| is
calculated, and :math:`Y(E)` is the Pendry :math:`Y` function. The |R factor|
can be calculated for all beams together or for each beam individually.\ [1]_
The :math:`Y` functions in Eq. :eq:`eq_RP_def` are computed from the beam
intensities :math:`I(E)`, their derivatives :math:`I'(E)=\frac{dI}{dE}`, and
the imaginary part of the inner potential |V0i| (see parameter :ref:`V0_imag`)
as

.. math::
    :label: eq_Y_def

    Y(E) = \frac{I(E)/I'(E)}{[I(E)/I'(E)]^2 + V_{0\mathrm{i}}^2}.

The beam intensities enter Eq. :eq:`eq_Y_def` via their logarithmic
derivative :math:`\frac{d}{dE}\big(\ln{I(E)}\big) = \frac{I'(E)}{I(E)}`.
This makes the Pendry |R factor| insensitive to differences in the absolute
intensities of the |IV| curves. The largest contributions to |RP| come from
differences in the *positions* of extrema, especially *minima*.

An |RP| value of zero corresponds to perfect agreement between curves. |RP|
equals one for statistically uncorrelated data, while values larger than one
indicate anticorrelation. For close-packed surfaces, |RP| values larger than
0.2 indicate a problem, such as an incorrect structural model. |R factor|\ s
for more open, corrugated surfaces, such as missing-row-reconstructed Pt(110),
may be around 0.2. |RP| values larger than 0.25–0.30 should be taken as an
indication of poor correspondence between calculated and experimental beams.
The best values of |RP| obtained by the
`Erlangen group <https://www.fkp.physik.nat.fau.eu/research-schneider/>`__
are below 0.05.\ [2]_

Note that some smoothing algorithms applied to both experimental and
calculated beams, such as the one suggested by
:cite:t:`pendryReliabilityFactorsLEED1980`, artificially reduce the |R factor|,
because they effectively raise the minima of the |IV| curves. At minima, where
the intensities approach zero, |RP| is especially sensitive to small
differences; artificially increasing the intensity at minima thus gives
smaller |RP| values. Some LEED programs apply such a smoothing; in those
cases smaller |R factor|\ s than those obtained with ViPErLEED will be
reported, but this does not indicate a better agreement between calculated
and experimental |IV| curves.

.. todo::
    Refer to issue where we have discussed this (viperleed-betatest #8, after
    moving it to main.

By default, ViPErLEED applies no additional smoothing when calculating the
|R factor| (e.g., during structure optimization). Thus, the :ref:`EXPBEAMS`
file should already contain smoothed data.

.. tip::
    We highly recommend to smooth experimental data beforehand using the
    |IV|-curve editor of the :ref:`imagej_plugins`.
    Using the :ref:`RFACTORSMOOTH` parameter for smoothing the experimental
    |IV| curves is discouraged, as the smoothing algorithm applied there is
    inferior to that used by the |IV|-curve editor.

.. [1] Notice that the Pendry |R factor| between two sets of beams is not the
       average of the |R factor|\ s between beam pairs, as sums over all beams
       enter both the numerator and the denominator in Eq. :eq:`eq_RP_def`.
.. [2] Unpublished data by Lutz Hammer and coworkers.
