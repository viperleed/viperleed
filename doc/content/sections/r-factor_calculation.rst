.. _r-factor_calculation:

============
The R factor
============

The reliability factor (R factor) is a measure for the deviation between
two :math:`I(V)` curves or two sets of :math:`I(V)` curves. 
Structure optimization (:ref:`search section<sec_search>`) minimizes the R
factor between the calculated and experimental :math:`I(V)` curves.

As the comparison of two curves is not an unambiguous task, multiple R factors
implementations exist.
The used R factor is chosen via the :ref:`R_FACTOR_TYPE<RFACTORTYPE>` parameter.
ViPErLEED supports Pendry' R factor :math:`R_P`
:cite:p:`pendryReliabilityFactorsLEED1980` and :math:`R_2`
:cite:p:`spornAccuracyQuantitativeLEED1998`.

-  :math:`R_P` compares the logarithmic derivatives of the :math:`I(V)` data,
   with a fix for the divergence of the logarithm when the intensity approaches
   zero.
-  :math:`R_2` is based on the mean-square difference of the (scaled)
   :math:`I(V)` curves.

.. note::
    Using :math:`R_P` is the **default setting** and highly encouraged since
    tests have shown that it leads to better results than :math:`R_2`
    :cite:p:`spornAccuracyQuantitativeLEED1998`.

.. _pendry_r:

The Pendry R factor
-------------------

The Pendry R factor :math:`R_P` can have values between 0 and 2 and is defined
as follows
:cite:p:`pendryReliabilityFactorsLEED1980,heinzElectronBasedMethods2013`:

.. math:: 
    R_{\mathrm{P}} = \frac{\sum_{h,k}\int (Y^\mathrm{theo}_{h,k}(E) - Y^\mathrm{exp}_{h,k}(E) )^2 dE }{\sum_{h,k}\int (Y^\mathrm{theo}_{h,k}(E)^2 + Y^\mathrm{exp}_{h,k}(E)^2) dE}

The R-factor can be calculated for all beams :math:`g` together or individually.
:math:`Y(E)` is the Pendry Y-function that contains the beam intensities
:math:`I(E)`, the derivative :math:`I'(E)=\frac{dI{E}}{dE}`, and the imaginary
part of the inner potential :math:`V_{0\text{i}}` (see parameter
:ref:`V0_imag`).

.. math:: 
    Y(E) = \frac{I(E)/I'(E)}{[I(E)/I'(E)]^2 + V_{0\text{i}}^2}

:math:`R_P` becomes 0 in the case of perfect agreement between curves.
For uncorrelated data, :math:`R_P` = 1, while values larger than 1 indicate
anti-correlation.
The best values of :math:`R_P` obtained by the 
`Erlangen group <https://www.fkp.physik.nat.fau.eu/research-schneider/>`__ are
below 0.05 [#]_.
For close-packed surfaces, values above 0.2 indicate a problem such as an 
ncorrect structure model.
R factors for more open surfaces, such as missing-row-reconstructed Pt(110), can
be higher than 0.2.

Note that some smoothing algorithms applied to both, the experimental and the
calculated data, such as the smoothing suggested by Pendry in
:cite:p:`pendryReliabilityFactorsLEED1980`, apparently lowers the R factor,
because it effectively raises the minima of the :math:`I(V)` curves.
At minima, where the intensities approach zero, :math:`R_P` is especially 
sensitive to small differences; artificially increasing the intensity there
leads to lower values of :math:`R_P`.
Some LEED programs apply such a smoothing; in this case lower R factors than
those obtained with ViPErLEED will be reported, but this does not indicate a
better agreement between calculated and experimental data.

By default, ViPErLEED applies no additional smoothing upon structure
optimization and calculating the R factor.
We highly recommend that experimental data should be smoothed beforehand using
the I(V) curve editor in the :ref:`ViPErLEED Spot Tracker<spot_tracker>`.
The :ref:`EXPBEAMS.csv<EXPBEAMS>` file should thus already contain the
smoothed data.


(Using the :ref:`R_FACTOR_SMOOTH<RFACTORSMOOTH>` parameter for smoothing the
experimental :math:`I(V)` curves is discouraged;
the smoothing algorithm applied there is inferior to that used by the
:math:`I(V)` curve editor .)

.. [#] Unpublished data by Lutz Hammer and coworkers.
