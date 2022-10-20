.. _rfactoranalysis:

The Rfactor_analysis files
==========================

The Rfactor_analysis files (Rfactor_analysis_refcalc.pdf and Rfactor_analysis_superpos.pdf) plot the results of the r-factor calculations, as well as rudimentary analysis.

The topmost plot on each page is the same as given in the :ref:`Rfactor_plots.pdf<Rfactorplots>`  files, i.e. the theoretical and experimental I(V) curves for a given beam, as well as the R-factor.

The second plot shows the `Y-functions <https://doi.org/10.1088/0022-3719/13/5/024>`__ of the I(V) curves, and the difference between the two. Since the R-factor is determined by integrating the squared difference between the Y-functions, this plot gives a rough idea where the main contributions to the R-factor come from.

Finally, the third plot on each page makes this more explicit by plotting the integral over (Y₁ - Y₂)² / (Y₁² + Y₂²), which is the exact term used in the R-factor. This allows quick identification of the areas contributing most to the R-factor.
