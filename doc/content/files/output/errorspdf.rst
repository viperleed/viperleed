.. include:: /substitutions.rst

.. _errorspdf:

Errors_summary.csv, Errors.zip, and Errors.pdf
==============================================

The files ``Errors_summary.csv``, ``Errors.zip``, and ``Errors.pdf``
contains the results of the :ref:`Error calculation<error_calculation>`.
Only one independent parameter is varied at a time, and the resulting
|R factor|\ s are collected. If multiple atoms are linked (e.g., by symmetry),
they will be varied together, and the entire group will be listed in labels.
Note that for geometrical variations, the actual direction of the displacements
may not be the same for the entire group.

.. note::
  When using the :ref:`Pendry <pendry_r>` |R factor|, the variance provides an
  indication of the uncertainty of a given parameter. The point of intersection
  between the error curve and :math:`R_{\mathrm{P,min}} + \mathrm{var}(R)` can
  be taken as a measure of the statistical error for each parameter
  :cite:p:`heinzElectronBasedMethods2013`.
  :math:`\mathrm{var}(R)` is defined as

  .. math::

    \frac{\mathrm{var}(R)}{R} = \sqrt{ \frac{8 * V_{0i} }{ \delta E} }

  where :math:`R` is the |R factor| of the best-fit structure, :math:`V_{0i}`
  is the imaginary part of the inner potential, and :math:`\delta E` is the
  total energy range of the beam set. For more information, see
  :cite:t:`pendryReliabilityFactorsLEED1980`.

  If the Pendry |R factor| is used, the value of :math:`R + \textrm{var}(R)`
  will be calculated for each error type and drawn as a horizontal line in
  ``Errors.pdf``. If ViPErLEED finds an intersection between the error curve
  for any parameter :math:`p` and the line :math:`R + \textrm{var}(R)`, the
  statistical errors for :math:`p` will be estimated and written to
  ``Errors_summary.csv`` and plotted in ``Errors.pdf``.

.. attention::
  The |R factor| values obtained in the error calculation contain
  errors from the :ref:`tensor-LEED approximation<tensor_leed>`.

Errors_summary.csv
==================

``Errors_summary.csv`` contains a summary of the error calculation results.
Every line lists one independently varied parameter with corresponding
information. The information given for each parameter is:

- atom numbers (as in :ref:`POSCAR<poscar>`),
- displacement mode (``geo``, ``vib``, ``occ``)
- displacement direction,
- minimum |R factor| of the error curve,
- variance of the |R factor| (see above; only
  applicable for Pendry |R factor|),
- parameter value at minimum |R factor| (labeled ``p_min``),
- statistical error estimates for the parameter fit (only if applicable,
  see above; labeled ``-Δp`` and ``+Δp``).

For geometrical displacements the column ``Direction`` lists the direction
requested in :ref:`DISPLACEMENTS<displacements>`.

The contents of ``Errors_summary.csv`` may look something like this:

.. code-block::

   Atoms, Mode, Direction,  R_min, var(R),  p_min,    -Δp,    +Δp
       1,  geo,         z, 0.0870, 0.0187, 0.0000, 0.0110, 0.0064
       2,  geo,         z, 0.0870, 0.0187, 0.0000, 0.0090, 0.0083
       1,  vib,       N/A, 0.0870, 0.0187, 0.0000,    N/A, 0.0161
       1,  occ,       N/A, 0.0883, 0.0190, 0.9500,    N/A,    N/A


Errors.zip
==========

Results for each individual parameter varied during the error calculation
are collected and stored in the ZIP archive ``Errors.zip``. Files are named
``Errors_{mode}_atoms#{ids}.csv``, where ``{mode}`` is one of ``geo``, ``vib``,
or ``occ`` and ``{ids}`` is the atom numbers (as in :ref:`POSCAR<poscar>`).
For geometrical displacements and vibrational amplitude changes, each file has
two columns: the first column lists the displacement and the right column the
corresponding |R factor|. For occupational errors, additional columns list the
occupations by chemical element. Geometrical displacements and vibrational
amplitude changes are given in units of Å, occupations in %.


Errors.pdf
==========

The same data contained in ``Errors_summary.csv`` and ``Errors.zip`` is plotted
in the ``Errors.pdf`` file (all placed in ``OUT``). Parameter are split by
geometry, vibrational amplitudes and site occupation. However, all parameters
of the same type (e.g., all geometrical displacements) are grouped together in
the ``Errors.pdf`` file, so if you calculate very different displacements in
the same error calculation, these would nevertheless be plotted together.

For each parameter type, the ``Errors.pdf`` file contains *one* plot in which
the results for *all* parameters of that type are shown together, as well as
separate single plots for each parameter.