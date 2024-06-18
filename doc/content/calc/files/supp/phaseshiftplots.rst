.. include:: /substitutions.rst

.. _phaseshiftplots:

Phaseshifts-plots.pdf
=====================

The Phaseshifts-plots.pdf file is generated during the initialization, and
serves purely for information. It displays the phase shifts that occur when
an electron of a given energy and angular momentum number *L* scatters at
an atom in a given site. One plot is generated for each site and element
that can occupy that site. The data is displayed exactly as found in the
:ref:`PHASESHIFTS` file, without any processing.

The order of the plots is the same as in the :ref:`PHASESHIFTS` files, and
reflects the order that |calc| expects to find for the elements and/or sites.
If you supplied :ref:`PHASESHIFTS` manually, instead of generating them
automatically, make sure that the phase shifts have been correctly assigned
to their given elements/sites.
