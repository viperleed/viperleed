.. _rfactoranalysis:

Rfactor_analysis.pdf
====================

The Rfactor_analysis files (Rfactor_analysis_refcalc.pdf and Rfactor_analysis_superpos.pdf) 
plot the results of the r-factor calculations, as well as rudimentary analysis.

The topmost plot on each page is the same as given in the :ref:`Rfactor_plots.pdf<Rfactorplots>` 
files, i.e. the theoretical and experimental I(V) curves for a given beam, as well as the R-factor.

The second plot shows the Y-functions (see e.g. :cite:t:`pendryReliabilityFactorsLEED1980`) 
of the I(V) curves, and the difference between the two. Since the R-factor is determined 
by integrating the squared difference between the Y-functions, this plot gives a 
rough idea where the main contributions to the R-factor come from.

Finally, the third plot on each page makes this more explicit by 
plotting the cumulative sum of the squared difference between the Y-functions 
:math:`\frac{(Y_1 - Y_2)^2}{(Y_1^2 + Y_2^2)}`, 
which is the term used in the Pendry R-factor.
This allows quick identification of the areas contributing most to the 
R-factor.


.. figure:: /_static/output_examples/Rfactor_analysis_refcalc_example.svg
   :width: 60%
   :align: center

   A R-factor analysis plots for a single beam calculated in the example on the 
   Hematite :math:`(012)-(1 \times 1)` surface discussed :ref:`here<example_Fe2O3>`.
