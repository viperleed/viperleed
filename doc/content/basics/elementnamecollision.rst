.. include:: /substitutions.rst

.. _elementnamecollision:

Element names
=============

The way TensErLEED works leads to a distinction in what an "element" is.
Two (sometimes overlapping) categories exist:

-  Elements as defined in the :ref:`POSCAR file<POSCAR>` (here often
   referred to as "POSCAR elements"): These do *not* need to be chemical
   elements, but can in principle have arbitrary names.\ [#]_
-  Chemical elements, as they appear in the periodic table: Each of the POSCAR
   elements that does not already have names of chemical elements must be
   specified by either the :ref:`ELEMENT_MIX<ELSPLIT>` or the
   :ref:`ELEMENT_RENAME<ELDEF>` parameters.

The main reason why this distinction is necessary is mixed occupation of sites.
It is possible for specific atomic positions (as specified in the POSCAR file)
to contain mixtures of two or more elements. This split should be specified by
:ref:`ELEMENT_MIX<ELSPLIT>`. In these cases, it is strongly recommended to
rename the elements in the POSCAR file such that the POSCAR elements and
chemical elements don't overlap. For example, if  your POSCAR contains
atoms labeled ``La`` but your material contains  some mixture of ``La``
and ``Sr`` in those sites, you could rename ``La`` to ``A`` in the POSCAR,
then use ``ELEMENT_MIX A = La Sr`` to specify that these atoms are actually
either ``La`` or ``Sr``. If you do not do this, but simply use
``ELEMENT_MIX La = La Sr``, you produce an **element-name collision**.
This means that in your later input (e.g., in the :ref:`VIBROCC<vibrocc>`
and :ref:`DISPLACEMENTS<DISPLACEMENTS>` files) |calc| has no clear way to
distinguish whether you wanted to refer to the POSCAR element (i.e., all
chemical species in that position) or to one specific chemical element.

Note that element names are not case sensitive. This means that, for example,
``Pb``, ``PB`` and ``pb`` are considered equal. Thus, having a POSCAR element
``PB`` defined as ``ELEMENT_MIX PB = P B`` and the chemical element ``Pb`` in
the same calculation will lead to a name collision.

As mentioned above, the recommended approach is to **avoid** name collisions by
renaming the POSCAR elements where necessary. However, name collisions do not
normally cause the program to stop, although they will produce a warning.
If you insist to work with a name collision, you must make sure that your
later input is interpreted correctly. Generally, the program's approach is
to interpret element names as POSCAR elements, that is, to assign your input
to all atoms in the given site, wherever this makes sense. Chemical elements
only take precedence where interpreting as POSCAR elements would be
nonsensical, that is on the right-hand side of input in the
:ref:`VIBROCC<vibrocc>` file, and in :ref:`OCC_DELTA<OCCDELTA>`
blocks in the :ref:`DISPLACEMENTS<DISPLACEMENTS>`  file.

.. rubric:: Footnotes

.. [#] Mind the fact that the
       `POSCAR-file specification <https://www.vasp.at/wiki/index.php/POSCAR>`__
       only considers the first two characters of each element name. This is
       also reflected in `VESTA <https://jp-minerals.org/vesta/en/>`__, which
       reads only the first two characters. This means that ``Atom1`` and
       ``Atom2`` will both appear as species ``At`` both in VESTA
       :cite:p:`mommaVESTAThreedimensionalVisualization2011`
       and in a DFT calculation with VASP :cite:p:`vasp`. ViPErLEED will
       instead read the full name.
