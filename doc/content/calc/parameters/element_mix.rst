.. _element_mix:

ELEMENT_MIX
===========

ELEMENT_MIX defines whether sites defined as one element in the
:ref:`POSCAR` file are occupied by multiple chemical elements.
See also :ref:`this page<occdelta>` for instructions on how to vary
the occupation of a site during structure optimization.

**Default:** None of the elements given in the :ref:`POSCAR` file are occupied
by multiple elements.

**Syntax:**

::

   ELEMENT_MIX A = La Sr

In the example, the element A should be an element defined in the :ref:`POSCAR`
file as the element for a list of sites, while La and Sr are the elements that
actually occupying these sites. In the :ref:`POSCAR` file, an atom in a given
site can only belong to one element. ELEMENT_MIX allows assigning a site to
multiple elements, which can have different properties (e.g., vibration
amplitudes), and read different phase-shift files during the LEED calculation.
The elements on the right-hand side should be actual chemical elements, with
the two-letter abbreviation as it is found in the periodic table. ELEMENT_MIX
is complementary to :ref:`ELEMENT_RENAME`, so no POSCAR element for which
ELEMENT_MIX is defined should appear as an ELEMENT_RENAME parameter, and
vice versa.

**Acceptable values**: Due to current limits in the code, a maximum of **five**
different species can occupy one site. If the occupations of the elements
defined in the in the :ref:`VIBROCC` file sum to less than one,
one of the five species will be a vacancy, so only **four** different elements
can occupy one site if their occupations do not sum to one.

If element splits for multiple elements of the POSCAR file are being defined,
have multiple lines starting with ELEMENT_MIX and the respective POSCAR
elements left of the '=' sign.

It is recommended to avoid overlap in the element names, e.g.,
``ELEMENT_MIX La = La Sr``. For more on this, see the page on
:ref:`element name collision<ElementNameCollision>`.

The actual occupations of the split sites are defined in the :ref:`VIBROCC`
file, and can be fitting parameters.
