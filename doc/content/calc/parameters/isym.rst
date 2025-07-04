.. _isym:

SYMMETRY_FIX
============

SYMMETRY_FIX allows you to constrain the symmetry of the displacements, 
vibration amplitudes, and concentrations during LEED optimization.

**Default**: SYMMETRY_FIX = True: Use the symmetry determined 
automatically from the :ref:`POSCAR` file.

**Syntax:**

::

   SYMMETRY_FIX = true
   SYMMETRY_FIX = false
   SYMMETRY_FIX = <group>

**Accepted values**: ``true``, ``True``, ``false``, ``False`` or any 
combination of small/capital letters or any abbreviation of the words 
``true`` and ``false`` (only the first character matters).
Alternatively, ``<group>`` can be the Hermann–Mauguin symbol for the 
planar group you want to constrain your symmetry to. Acceptable values 
are listed in the table below.

.. note:: 
  The program will evaluate the symmetry group of your input 
  :ref:`POSCAR`, and will allow only symmetry *reduction*, 
  i.e., the group you select must be a subgroup of the group of your 
  input :ref:`POSCAR`.
  Take a look :ref:`here<planegroups>` for (i) a graphical representation of the 
  operations in the group of your slab, (ii) a graphical representation 
  of the diplacements allowed for each plane group, and (iii) a list of 
  subgroups for each plane group.
  The symmetry group of your input slab is detected from the atomic 
  coordinates in your input 
  :ref:`POSCAR` (see also the parameter 
  :ref:`sym_eps`), and written to the POSCAR header after 
  initialization.

+-----------+------------------------------------------------------------------+
| ``group`` | Special cases: reduction from some groups requires \             |
|           | specification of which subgroup (see :ref:`here<planegroups>`)   |
+===========+==================================================================+
| ``p1``    | --                                                               |
+-----------+------------------------------------------------------------------+
| ``p2``    | --                                                               |
+-----------+------------------------------------------------------------------+
| ``pm``    | from ``pmm``, ``p4m`` or ``rcmm`` specify \                      |
|           | either ``pm[1 0]`` or ``pm[0 1]``                                |
+-----------+------------------------------------------------------------------+
| ``pg``    | from ``pgg``, ``p4g`` or ``rcmm`` specify \                      |
|           | either ``pg[1 0]`` or ``pg[0 1]``                                |
+-----------+------------------------------------------------------------------+
| ``cm``    | from ``p4m``, ``p4g`` or ``rcmm`` specify \                      |
|           | either ``cm[1 1]`` or ``cm[1 -1]``                               |
+-----------+------------------------------------------------------------------+
|           | from ``p3m1`` specify ``cm[1 -1]``, ``cm[2 1]``, or ``cm[1 2]``  |
+-----------+------------------------------------------------------------------+
|           | from ``p31m`` specify ``cm[1 0]``, ``cm[0 1]``, or ``cm[1 1]``   |
+-----------+------------------------------------------------------------------+
|           | from ``p6m`` specify ``cm[1 0]``, ``cm[0 1]``, \                 |
|           | ``cm[1 1]``, ``cm[1 -1]``, ``cm[1 2]``, or ``cm[2 1]``           |
+-----------+------------------------------------------------------------------+
| ``rcm``   | Note: not a real plane group, see :ref:`here<planegroups>` \     |
|           | (``cm`` with non-primitive cell. Calculation times will be \     |
|           | longer than with the primitive cell!)                            |
+-----------+------------------------------------------------------------------+
|           | from ``rcmm`` specify either ``rcm[1 0]`` or ``rcm[0 1]``        |
+-----------+------------------------------------------------------------------+
| ``pmm``   | --                                                               |
+-----------+------------------------------------------------------------------+
| ``pmg``   | from ``rcmm`` specify either ``pmg[1 0]`` or ``pmg[0 1]``        |
+-----------+------------------------------------------------------------------+
| ``pgg``   | --                                                               |
+-----------+------------------------------------------------------------------+
| ``cmm``   | from ``p6m`` specify ``cmm[1 2]`` (or ``cmm[1 0]``), \           |
|           | ``cmm[2 1]`` (or ``cmm[0 1]``), \                                |
|           | ``cmm[1 -1]`` (or ``cmm[1 1]``). Alternatives in \               |
|           | parentheses are equivalent                                       |
+-----------+------------------------------------------------------------------+
| ``rcmm``  | Note: not a real plane group, see :ref:`here<planegroups>` \     |
|           | (``cmm`` with non-primitive cell. Calculation times will be \    |
|           | longer than with the primitive cell!)                            |
+-----------+------------------------------------------------------------------+
| ``p4``    | --                                                               |
+-----------+------------------------------------------------------------------+
| ``p4m``   | --                                                               |
+-----------+------------------------------------------------------------------+
| ``p4g``   | --                                                               |
+-----------+------------------------------------------------------------------+
| ``p3``    | --                                                               |
+-----------+------------------------------------------------------------------+
| ``p3m1``  | --                                                               |
+-----------+------------------------------------------------------------------+
| ``p31m``  | --                                                               |
+-----------+------------------------------------------------------------------+
| ``p6``    | --                                                               |
+-----------+------------------------------------------------------------------+
| ``p6m``   | --                                                               |
+-----------+------------------------------------------------------------------+
