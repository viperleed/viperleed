.. _elementnamecollision:

Element names
=============

The way TensErLEED works leads to a distinction in what an "element" is.
Two (sometimes overlapping) categories exist:

-  Elements as defined in the :ref:`POSCAR file<POSCAR>` (here often 
   referred to as "POSCAR elements"): These do *not* need to be chemical
   elements, but can in principle have arbitrary names. [#]_
-  Chemical elements, as they appear in the periodic table:
   If the POSCAR elements do not already have names of chemical 
   elements, these need to be specified by the
   :ref:`ELEMENT_MIX<ELSPLIT>`  and/or :ref:`ELEMENT_RENAME<ELDEF>` 
   parameters.

The main reason why this distinction is necessary is mixed site 
occupation.
When a specific atomic position (as specified in the POSCAR file) is not
chemically pure, but contains some mixture of two (or more) elements, 
this split needs to be specified by :ref:`ELEMENT_MIX<ELSPLIT>`. It is 
strongly recommended to rename the elements in the POSCAR file such that
the POSCAR elements and chemical elements don't overlap. For example, if
your POSCAR contains atoms labelled ``La``, but in reality, your 
material contains some mixture of ``La`` and ``Sr`` in those sites, you 
could rename ``La`` to ``A`` in the POSCAR, then use 
``ELEMENT_MIX A = La Sr`` to specify that these atoms are actually either
``La`` or ``Sr``. If you do not do this, but simply use 
``ELEMENT_MIX La = La Sr``, you produce an **element name collision**, 
meaning that in your later input (e.g. in the :ref:`VIBROCC<vibrocc>` 
and :ref:`DISPLACEMENTS<DISPLACEMENTS>`  files), there will be no clear 
way to distinguish whether you wanted to refer to the POSCAR element 
(i.e. all atoms in that position), or one specific chemical element.

Note that element names are not case-sensitive, i.e., 
``Pb``, ``PB`` and ``pb`` are considered equal. Thus, having a POSCAR 
element ``PB`` defined as ``ELEMENT_MIX PB = P B`` and the chemical 
element ``Pb`` in the same calculation will lead to a name collision.

As mentioned above, the recommended approach is to **avoid** name 
collisions (by renaming the POSCAR elements where necessary).
However, name collisions are generally not program-stopping, although 
they will produce a warning. If you insist to work with a name collision,
you must make sure that your later input is interpreted correctly. 
Generally, the program's approach is to interpret element names as 
POSCAR elements, i.e. assign your input to all atoms in the given site, 
wherever this makes sense.
Chemical elements will only take precedence 
where interpreting as POSCAR elements would be nonsensical, that is on 
the right-hand side of input in the :ref:`VIBROCC<vibrocc>` file, and 
in :ref:`OCC_DELTA<OCCDELTA>`  blocks in the 
:ref:`DISPLACEMENTS<DISPLACEMENTS>`  file.

.. rubric:: Footnotes

.. [#] Mind the fact that `VESTA <https://jp-minerals.org/vesta/en/>`__ 
       will only read in the first 2 characters of the names in a POSCAR
       file, so ``Atom1`` and ``Atom2`` will both appear as species 
       ``At`` in VESTA :cite:p:`mommaVESTAThreedimensionalVisualization2011`.
