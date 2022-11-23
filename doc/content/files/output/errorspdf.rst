.. _errorspdf:

Errors.pdf and Errors.csv
=========================

The ``Errors.csv`` file contains the data resulting from an :ref:`Error calculation<error_calculation>`.
Only one independent parameter is varied at the time, and the resulting R-factors of varying that parameter are given in the right-most column.
If multiple parameters were varied, the results are completely independent, and are printed back-to-back.
If multiple atoms are linked (e.g. by symmetry), they will be varied together, and the entire group will be listed in the left-most column "Atoms".
Note the for geometrical variations, the actual direction of the displacements may not be the same for the entire group.
The direction given in the ``Errors.csv`` file is the displacement direction for the *first* atom listed in the "Atoms" column.

The same data contained in ``Errors.csv`` is plotted in the ``Errors.pdf`` file.
Parameter are split by whether they concern geometry, vibrational amplitudes or site occupation.
However, all parameters of the same type (e.g. all geometrical displacements) are grouped together in the ``Errors.pdf`` file, so if you calculate very different displacements in the same error calculation, these would nevertheless be plotted together.

For each parameter type, the Errors.pdf file contains *one* plot in which the results for *all* parameters of that type are shown together, as well as separate single plots for each parameter.

When using Pendry's R factor, the variance provides an indication of the uncertainty of a given parameter.
:math:`\textrm{var}(R)` is defined as:

.. math::

  \frac{\textrm{var}(R)}{R} = \sqrt{ \frac{8 * V_{0i} }{ \delta E} }

where :math:`R` is the R-factor of the best-fit structurewhere, :math:`V_{0i}` is the imaginary part of the inner potential, and :math:`\delta E` is the total energy range of the beamset.
For more information, see :cite:t:`pendryReliabilityFactorsLEED1980`.
When using the Pendry R-factor, the value of :math:`R + \textrm{var}(R)` will also be drawn as a horizontal line in Errors.pdf.
