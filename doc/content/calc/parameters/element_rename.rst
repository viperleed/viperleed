.. _element_rename:

ELEMENT_RENAME
==============

ELEMENT_RENAME defines what chemical element should be used for an element
in the :ref:`POSCAR` file, if the name given in the POSCAR and the
chemical element symbol differ.

**Default:** Assume that each of the element names in the :ref:`POSCAR`
file matches the symbol of a chemical element (case-insensitive).

**Syntax:**

::

   ELEMENT_RENAME A = Fe
   ELEMENT_RENAME B = O

**Acceptable values**: On the left, :ref:`POSCAR`  names (exact match, apart 
from case). On the right, chemical elements from the periodic table (not
case sensitive).

In the first line of the example, A is the element name present in the
:ref:`POSCAR` file (line after the definition of the unit cell
vectors) to define a specific site, while Fe is the chemical element to
be used for these sites.

**Notes**:

-  The definition need not be made if the POSCAR element is named the same as
   the chemical element (i.e. ``ELEMENT_RENAME Fe = Fe`` is unnecessary). See
   :ref:`element name collision<ElementNameCollision>`  for more on the
   distinction between POSCAR elements and chemical elements.
-  If chemical elements are being defined for multiple POSCAR elements, have
   multiple lines starting with ELEMENT_RENAME and the respective POSCAR
   elements left of the '=' sign.
-  ELEMENT_RENAME is complementary to :ref:`ELEMENT_MIX`, so no POSCAR element 
   for which ELEMENT_RENAME is defined should appear as an ELEMENT_RENAME 
   parameter, and vice versa.
-  The VESTA program for viewing structures reads only the first two
   characters of the POSCAR element name
   :cite:p:`mommaVESTAThreedimensionalVisualization2011`. If you use VESTA,
   it is not a good idea to have "Atom1" and "Atom2" for different elements
   in the POSCAR, as both will be considered the same ("At").
