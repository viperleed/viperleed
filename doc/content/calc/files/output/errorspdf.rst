.. include:: /substitutions.rst

.. _errorspdf_header:

Errors_summary.csv, Errors.zip, and Errors.pdf
==============================================

The files ``Errors_summary.csv``, ``Errors.zip``, and ``Errors.pdf``
contains the results of the :ref:`Error calculation<error_calculation>`.
Only one independent parameter is varied at a time, and the resulting
|R factor|\ s are collected. If multiple atoms are linked (e.g., by symmetry),
they will be varied together, and the entire group will be listed in labels.
Note that for geometric variations, the actual direction of the displacements
may not be the same for the entire group.


.. _errorscsv:

Errors_summary.csv
==================

``Errors_summary.csv`` contains a summary of the error calculation results.
Every line lists one independently varied parameter with corresponding
information. The information given for each parameter is:

- atom numbers (as in :ref:`POSCAR`),
- displacement mode (``geo``, ``vib``, ``occ``)
- displacement direction,
- minimum |R factor| of the error curve,
- standard error of the |R factor| (labeled ``var(R)``. See
  :ref:`error_calculation`; only applicable for Pendry's |R factor|),
- parameter value at minimum |R factor| (labeled ``p_min``),
- statistical error estimates for the parameter fit (labeled ``-Δp``
  and ``+Δp``. Only for Pendry's |R factor| and only if the error
  curve reaches ``R_min + var(R)`` in the displacement range considered.
  See also :ref:`error_calculation`.).

For geometric displacements the column ``Direction`` lists the direction
requested in :ref:`DISPLACEMENTS`.

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
or ``occ`` and ``{ids}`` is the atom numbers (as in :ref:`POSCAR`).
For geometric displacements and vibration-amplitude changes, each file has
two columns: the first column lists the displacement and the right column the
corresponding |R factor|. For occupational errors, additional columns list the
occupations by chemical element. Geometric displacements and
vibration-amplitude changes are given in units of Å, occupations in %.


.. _errorspdf:

Errors.pdf
==========

The same data contained in ``Errors_summary.csv`` and ``Errors.zip`` is plotted
in the ``Errors.pdf`` file (all placed in ``OUT``). Parameter are split by
geometry, vibration amplitudes, and site occupation. However, all parameters
of the same type (e.g., all geometric displacements) are grouped together in
the ``Errors.pdf`` file, so if you calculate very different displacements in
the same error calculation, these would nevertheless be plotted together.

For each parameter type, the ``Errors.pdf`` file contains *one* plot in which
the results for *all* parameters of that type are shown together, as well as
separate single plots for each parameter.