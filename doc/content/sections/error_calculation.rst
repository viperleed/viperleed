.. include:: /substitutions.rst

.. _error_calculation:

Error calculations
==================

Once a best-fit structure has been determined, it is useful to see how
strongly small changes to specific parameters affect the |R factor|.
In the error calculation, displacements (given by the
:ref:`DISPLACEMENTS file<DISPLACEMENTS>`) are applied to one parameter
at a time, and the |R factor| is output for each step of the variation 
range. If multiple parameters are linked (e.g., by symmetry), these 
parameters are treated as one, and varied together.

To run the error calculation, set :ref:`RUN<run>` as ``RUN = 5`` in the
:ref:`PARAMETERS<PARAMETERS>` file. It is recommended to first run a
reference calculation in the same run (e.g. ``RUN = 1 5``),
as the error curves may not be centered otherwise.

Required input files are the same as for running delta calculations and a
search, i.e., the structural input files, :ref:`experimental beams<EXPBEAMS>`,
a set of :ref:`Tensors<Tensorszip>` from a reference calculation, and a
:ref:`DISPLACEMENTS file<DISPLACEMENTS>` defining what parameters should
be varied. Note that defining multiple search sections in the
:ref:`DISPLACEMENTS file<DISPLACEMENTS>`, as is possible for the search,
is not allowed here: Instead, only the first search section of the
:ref:`DISPLACEMENTS file<DISPLACEMENTS>` will be read (or the last,
if the error calculation is run following a search). Defining geometrical,
vibrational and occupation variations all in the same
:ref:`DISPLACEMENTS file<DISPLACEMENTS>` is allowed, but the different
variations will be split up, so the result is the same as executing
multiple error calculations. This means, you cannot have error calculations
for multiple geometrical displacement directions (e.g., :math:`x` *and*
:math:`z`) at the same time since this would require multiple consecutive
search sections in the DISPLACEMENTS file.

The error calculation does *not* require a set of
:ref:`Delta files<Deltaszip>`, since the normal delta calculation routines
mix geometrical and vibrational displacements. Instead, the error calculation
will run the required delta calculations automatically, splitting the
geometrical and vibrational variations into separate delta files to
reduce computational cost.

The results of the error calculation will be output into the files
:ref:`Errors.pdf and Errors.csv<errorspdf>`.

.. note::

    -  If you find that the |R factor| is very insensitive to the displacement
       of a given atom (much less sensitivity than for other atoms with a
       similar depth and similar scattering properties), this is an indication
       that the respective atom is either absent or its position is **far**
       from the true position. In this case, you may want to consider
       increasing the displacement range for this atom. Note also the
       what is "far away" depends on how strongly the atom scatters
       (i.e., chemical species and depth), and in some cases may be
       as small as 0.1Å, e.g. for a z variation of a surface atom.

    -  Hydrogen is a very weak scatterer; the |R factor| depends only weakly
       on its position.
    -  If a site can be occupied by different chemical elements, the site
       occupation (i.e., the element concentrations) and vibration amplitude
       can be strongly correlated :cite:p:`blumSegregationOrderingFe12001`.
       In such a case, the increase of the |R factor| when changing one of
       these parameters is not a good indication for the error of that
       parameter.
    -  Simultaneous geometrical/vibrational variation of multiple chemical
       elements occupying the same site is possible, and the displacement
       ranges for the different elements may differ. However, all displacement
       ranges must have the same length.
