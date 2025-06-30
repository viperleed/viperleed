.. _sitedef:

SITE_DEF
========

SITE_DEF defines the different types of sites for a given element.
Different sites can have different occupations and vibration amplitudes.

**Default:** all atoms have site name ``def``.

**Syntax:**

.. code-block:: none

   SITE_DEF Fe = tet 1 2, oct 4-10 12, surf top(2)
   SITE_DEF O = surf top(2)

The first argument is the element name in the :ref:`POSCAR` file for
sites that are being defined.
After the "=" sign, comma-delimited groups should contain names of the sites
and the atom numbers to be assigned (number of atom in the POSCAR list /
same as progressive numbers in :term:`VESTA`'s
Edit>Data>Structure parameters, **not** the :term:`VESTA` number per type,
i.e., not the number in 'Fe2', 'O16' etc.). The top(*N*) function selects
the *N* topmost atoms of the given POSCAR element. ``4-10`` or ``4:10`` will
select atoms 4 to 10, including 4 and 10. In the example, atoms 1 and 2 will
get site type 'tet', atoms 4 through 10 and atom 12 will get site type 'oct',
and the two topmost Fe atoms will get site type 'surf'. The two topmost O
atoms will get site type 'surf'.
For atoms identified by number, the element label must match the element label
specified in the POSCAR. Otherwise an error will be raised asking the user to
check the input.

If site types for multiple elements are being assigned, have multiple lines,
each starting with SITE_DEF and the respective elements left of the '=' sign,
as in the example.

Only atoms that need special treatment need to be listed. All others will get
site type ``'def``', which can be accessed normally. For most simple cases, a
simple assignment of the topmost (undercoordinated) atoms as special sites
should be sufficient.
The following example is for Fe2O3(012), where the topmost layer has two Fe
and two O atoms per unit cell:

::

   SITE_DEF Fe = atTop top(2)
   SITE_DEF O = atTop top(2)
