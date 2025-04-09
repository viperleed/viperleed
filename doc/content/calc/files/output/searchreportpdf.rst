.. include:: /substitutions.rst

.. _searchreportpdf:

Search-report.pdf
=================

The Search-report.pdf file is meant to collate the overall progress of the
current optimization, stretching over multiple searches. If only one search
is performed, it is mostly redundant with the
:ref:`Search-progress.pdf<searchprogresspdf>` file.

In addition to Search-report.pdf a Search-report.csv file is also created
which contains the same information in a tabular format. This file can be
imported into other programs for custom analysis.

R-factor / Generations
~~~~~~~~~~~~~~~~~~~~~~

The |R factor| is plotted over the total number of generations, similar as
in the :ref:`Search-progress.pdf<searchprogresspdf>` file. Vertical lines
indicate the transition from one search to the next, and which search is
being performed is labelled at the top of the plot. As in
:ref:`Search-progress.pdf<searchprogresspdf>`, the black line represents
the |R factor| of the best-fit structure, the grey area is the range between
the best and worst |R factor|, and the blue line is the average |R factor| of
the entire population. Both the range and the average can and do often
increase when the next search begins.

Parameter scatter / Generations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The x-axis is shared with the |R-factor| plot above. This plot is supposed
to give a measure of how well a given search has converged. The 'parameter
scatter' being drawn is defined as the standard deviation of parameter
positions in the parameter scatter plots below, where each displacement
range is normalized to the interval [0, 1]. The black line indicates the
standard deviation for the *most* scattered parameter, i.e., the least
well-defined or well-converged parameter. The blue line tracks the mean
of those standard deviations. If the black line reaches zero, then this
means that the entire population has converged to the same values for
every parameter.

Parameter scatter plots
~~~~~~~~~~~~~~~~~~~~~~~

The parameter scatter plots in the Search-report.pdf file are direct copies
of the latest such plots in the :ref:`Search-progress.pdf<searchprogresspdf>`
files. If a search is performed in a loop, only the latest iteration of the
loop is found in the Search-report.pdf file.

Example
~~~~~~~

Below are is an example of a Search-report.pdf file that was generated a
structural optimization discussed in the example calculation discussed
:ref:`here<example_Fe2O3>`.


.. figure:: /_static/output_examples/Search-report_page_1.svg
   :width: 60%
   :align: center

   Example of the first page of a Search-report.pdf file. The upper plot
   shows the |R factor| over generations, the lower plot shows the number
   of generations since the last improvement.


.. figure:: /_static/output_examples/Search-report_page_2.svg
   :width: 60%
   :align: center

   Example of the second page of a Search-report.pdf file showing the
   parameter scatter plots. The upper three lines show z-displacements,
   the lowest line shows changes in vibration amplitude.
