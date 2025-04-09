.. include:: /substitutions.rst

.. _searchprogresspdf:

Search-progress.pdf
===================

The Search-progress.pdf file is meant to give an overview of the current status
of the search, and of the convergence, while the search is running. It contains
several different plots, listed in the following.
See also :ref:`Search-report.pdf<searchreportpdf>`.

In addition to Search-progress.pdf a Search-progress.csv file is also created
which contains the same information in a tabular format. This file can be
imported into other programs for custom analysis.

R-factor / Generations
~~~~~~~~~~~~~~~~~~~~~~

The total number of generations is a reasonably good measure of the elapsed
time, so plotting the |R factor|\ s over generations gives an impression on how
much progress is being made. The best |R factor| (i.e., corresponding to the
current best-fit structure) is drawn as a black line. The blue line tracks
the average |R factor| over all members of the population; the grey area plots
the range between the current best and current worst structure. Note that both
the range and the average may sometimes increase, for example if the population
is re-initialized (partly) at random, or if new population members are
generated due to :ref:`SEARCH_CULL`.

At regular intervals, the |R-factor| scatter for the entire population is
plotted as points. The best-fit population is black, the worst-fit population
is red, and the rest is shaded relative to their distance from the two. If
multiple population members share the same |R factor|, they are drawn as a
bigger point.

Whenever the search is stopped to reduce the
:ref:`SEARCH_CONVERGENCE gaussian<rmut>` parameter, a vertical line is drawn
and the new width of the new value is indicated at the top of the plot.

Generations since change
~~~~~~~~~~~~~~~~~~~~~~~~

The x-axis is shared with the |R-factor| plot above. On the y-axis is the
number of generations that have passed since the last improvement to any
structure. So, every time any structure improves, a point is added to the
scatter plot, marking how many generations had passed before this improvement.
Color indicates whether the improvement occured in the best structure, in one
of the best 10% of structures, or any other structure (reflecting the
:ref:`SEARCH_CONVERGENCE gaussian<rmut>` criteria). If the delta values
become very large, then this means that few improvements are being found,
so the search is either converged or the current settings are badly suited
for finding improvements. These situations may automatically be handled by
the :ref:`SEARCH_CONVERGENCE`  settings, i.e., by
decreasing the :ref:`SEARCH_CONVERGENCE gaussian<rmut>` parameter whenever
the generation delta becomes large.

Parameter scatter plots
~~~~~~~~~~~~~~~~~~~~~~~

For each parameter under variation, a scatter plot is drawn for the most recent
state of the entire population. The y-axis of these scatter plots corresponds
to the displacements range of the given parameter. Small labels at the top and
bottom indicate the limits of the range (as displacements relative to the
positions of the reference calculation). Color coding of the points is the
same as for the right-most scatter in the |R-factor|/Generation plot, i.e.,
the current best configuration is drawn as black, the worst as red, and the
rest are shaded relative to their distance in R. In addition, the current
best-fit population is indicated by black arrows in order to identify it
unambiguously.

Parameter scatters are drawn as half-transparent if a parameter is coupled to
another parameter. See the :ref:`searchpars.info<searchparsinfo>`  file if you
are unsure how exactly parameters are coupled.

Point size represents the number of population members which have the same
value for the parameter. If more than one structure has the same parameter
value, the color (i.e., |R factor|) is indicated for the best one; for example,
if the best and the worst structure both have a given parameter at the center
of the displacement range, then the point at the center of the displacement
range is nevertheless drawn as black.

..
   Commented out for now, as N-dimensional parabola fits were deprecated in 0.11
   If :ref:`N-dimensional parabola fits<parabola_fit>` to the |R-factor| data
   were performed, then the parabola minimum and error bars are also indicated.
   The minimum is indicated by a red or green diamond (red if neither of the
   error bars is small enough to fit within the displacement range). Error bars
   are only drawn if at least one edge falls within the displacement range.
   The left error bar corresponds to the uncorrelated parameter error, the
   right to the correlated parameter error.


Example
~~~~~~~

Below are examples of a Search-progress.pdf file that was generated during the
example calculation discussed :ref:`here<example_Fe2O3>`.

.. figure:: /_static/output_examples/Search-progress.svg
   :width: 60%
   :align: center

   Example of the first page of a Search-progress.pdf file. The upper plot
   shows the |R factor| during the search, the lower plot shows the generations
   since the last change.
