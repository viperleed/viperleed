.. _rfactorplots:

The Rfactor_plots files
=======================

The Rfactor_plots files (``Rfactor_plots_refcalc.pdf`` and ``Rfactor_plots_superpos.pdf``) plot the results of the r-factor calculations. Theoretical and experimental I(V) curves are drawn for every beam, and the R-factor per beam is given.

``Rfactor_plots_refcalc.pdf`` is generated after the :ref:`R-factor calculation<r-factor_calculation>` following the :ref:`refercence calculation<ref-calc>` (section :ref:`RUN = 11<run>`), while ``Rfactor_plots_superpos.pdf`` is generated after the :ref:`R-factor calculation<r-factor_calculation>` following the :ref:`superposition<super_pos>` (section :ref:`RUN = 12<run>`).

Note that the curves are the ones output by the TensErLEED r-factor program (theo.column and exp.column), not the ones found in :ref:`EXPBEAMS.csv<EXPBEAMS>`  / :ref:`THEOBEAMS.csv<THEOBEAMS>`  / :ref:`FITBEAMS.csv<FITBEAMS>`, and may differ from the original data, in that an inner potential shift (:ref:`IV_SHIFT_RANGE<IVSHIFTRANGE>`), polynomial interpolation of both theoretical and experimental curves, as well as smoothing of experimental curves (:ref:`R_FACTOR_SMOOTH<RFACTORSMOOTH>`) and normalization may have been applied.

The appearance of the R-factor plots can be modified with the :ref:`PLOT_IV<PLOT_COLORS_RFACTOR>` parameter.


**TODO Alex: include an example here, possibly from Hematite once ready**