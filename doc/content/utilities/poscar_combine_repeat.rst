.. _poscar_combine_repeat:

=========================================
Combine-POSCAR and POSCAR_get_bulk_repeat
=========================================

``Combine-POSCAR.py`` and ``POSCAR_get_bulk_repeat.py`` are two simple utility scripts for semi-manual modification of :ref:`POSCAR files<poscar>` included with ViPErLEED.


Combine-POSCAR Usage
====================

Combine-POSCAR takes a slab POSCAR and adds a bulk POSCAR on the bottom, rescaling the unit cell. 

To use the Combine-POSCAR utility, simply call the :term:`Python` script ``Combine-POSCAR.py`` from the command line (e.g. ``python Combine-POSCAR.py``).
The script will promt the user for the filenames of the surface slab and bulk :ref:`POSCAR files<poscar>`.
If the two structures are incompatible (if e.g. the :math:`\vec{a}` and :math:`\vec{b}` vectors mismatch), the utility will raise an error.
Otherwise it will write the combined POSCAR.

POSCAR_get_bulk_repeat Usage
============================

POSCAR_get_bulk_repeat reads a POSCAR file, asks at what c value the bulk starts, then automatically reduces the size of the POSCAR to non-redundant bulk layers only, and outputs the appropriate :ref:`N_BULK_LAYERS<blay>` and :ref:`BULK_REPEAT` values.

To use the POSCAR_get_bulk_repeat utility, simply call the :term:`Python` script ``POSCAR_get_bulk_repeat.py`` from the command line (e.g. ``python POSCAR_get_bulk_repeat.py``).
The script will instruct the user and promt for a POSCAR file name.
It will then read the structure and then promt for a cutoff c value and a desired tolerance for symmetry search.
If the bulk can be determined, the determined values for the paramters :ref:`BULK_REPEAT`, :ref:`N_BULK_LAYERS<blay>` and :ref:`LAYER_CUTS<ctrunc>` are output.
Additionally, the files ``POSCAR_bulk`` containing the bulk unit-cell and a file ``POSCAR_min`` containing the minimal surface slab will be written.
