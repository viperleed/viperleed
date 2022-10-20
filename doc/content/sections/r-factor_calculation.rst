.. _r-factor_calculation:

The R factor
============

The reliability factor (R factor) is a measure for the deviation between two *I*\ (*V*) curves or two sets of *I*\ (*V*) curves. Structure optimization (:ref:`search<Search>`) minimizes the R factor between the calculated and experimental *I*\ (*V*) curves.

As the comparison of two curves is not an unambiguous task, multiple R factors implementations exist. The used R factor is chosen via the :ref:`R_FACTOR_TYPE<RFACTORTYPE>`  parameter. ViPErLEED supports Pendry' R factor *R*\ :sub:`P` [1] and *R*\ :sub:`2`.

-  *R*\ :sub:`P` compares the logarithmic derivatives of the *I*\ (*V*) data, with a fix for the divergence of the logarithm when the intensity approaches zero.
-  *R*\ :sub:`2` is based on the mean-square difference of the (scaled) *I*\ (*V*) curves.

Using *R*\ :sub:`P` is the default setting and highly encouraged since tests have shown that it leads to better results than *R*\ :sub:`2` [2].

Pendry' R factor *R*\ :sub:`P` can have values between 0 and 2; it becomes 0 in the case of perfect agreement between curves. For uncorrelated data, *R*\ :sub:`P` = 1. Values larger than 1 indicate anti-correlation. The best values of *R*\ :sub:`P` obtained by the Erlangen group are below 0.05. For close-packed surfaces, values above 0.2 indicate a problem such as an incorrect structure model. R factors for more open surfaces, such as missing-row-reconstructed Pt(110), can be higher than 0.2.

Note that some smoothing algorithms applied to both, the experimental and the calculated data, such as the smoothing suggested by Pendry in Ref. [1], apparently lowers the R factor, because it effectively raises the minima of the *I*\ (*V*) curves. At minima approaching zero, *R*\ :sub:`P` is especially sensitive to small differences; artificially increasing the intensity there leads to lower values of *R*\ :sub:`P`. Some LEED programs apply such a smoothing; in this case lower R factors than those obtained with ViPErLEED will be reported, but this does not indicate a better agreement between calculated and experimental data. The ViPErLEED package uses only smoothing of the experimental *I*\ (*V*) curves to suppress noise; the smoothing algorithm employed does not raise the minima. Smoothing of the experimental data should be done by the *I*\ (*V*) curve editor in the Spot Tracker package; its output should be used as :ref:`EXPBEAMS.csv<EXPBEAMS>`  file. By default, ViPErLEED applies no additional smoothing upon structure optimization and calculating the R factor. (Using the :ref:`R_FACTOR_SMOOTH<RFACTORSMOOTH>`  parameter for smoothing the experimental *I*\ (*V*) curves is discouraged; the smoothing algorithm applied there is inferior to that used by the *I*\ (*V*) curve editor.)

| [1] J. B. Pendry, *Reliability Factors for LEED Calculations* `J. Phys. C: Solid State Phys. 13, 937 (1980) <http://dx.doi.org/10.1088/0022-3719/13/5/024>`__.
| [2] M. Sporn, E. Platzgummer, S. Forsthuber, M. Schmid, W. Hofer, and P. Varga, *The Accuracy of Quantitative LEED in Determining Chemical Composition Profiles of Substitutionally Disordered Alloys: A Case Study* `Surf. Sci. 416, 423 (1998) <http://dx.doi.org/10.1016/S0039-6028(98)00596-2>`__.
