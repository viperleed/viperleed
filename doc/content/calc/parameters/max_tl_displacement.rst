.. _max_tl_displacement:

MAX_TL_DISPLACEMENT
===================

MAX_TL_DISPLACEMENT defines a maximum displacement which is acceptable for
consecutive searches without performing a new reference calculation. The total
displacements with respect to the last reference calculation are evaluated
after each block in :ref:`DISPLACEMENTS` for each individual atom, and an action
is triggered if at least one atom is displaced by more than the threshold. This
allows to e.g. break out of looped searches when atoms move very far from their
original positions, and to either stop at that point or perform a new reference
calculation and continue.

By default, the same threshold is used for geometric and vibration amplitude
displacements. A distinct value for vibration amplitudes can be given if
required.

Set ``MAX_TL_DISPLACEMENT action = ignore`` to disable this parameter. This may
be useful for testing, but is otherwise **not** recommended.

**Default**: MAX_TL_DISPLACEMENT = 0.15 (i.e., 0.15 Angstrom)

The default action is to perform a new reference calculation and continue if the
last reference calculation took less than 30 minutes, and to stop otherwise.

**Syntax:**

::

   ! only define different threshold values:
   MAX_TL_DISPLACEMENT = 0.1       ! sets the threshold to 0.1 Angstrom for both geo and vib
   MAX_TL_DISPLACEMENT = 0.1 0.05  ! 0.1 Angstrom for geometric displacements, 0.05 Angstrom for vibrations
   
   ! define thresholds using flags:
   MAX_TL_DISPLACEMENT geo = 0.1   ! same as above, 'geo' is the default
   MAX_TL_DISPLACEMENT vib = 0.01  ! re-define the vibration threshold, leaving geo unchanged
   
   ! define the action:
   MAX_TL_DISPLACEMENT action = stop           ! Always stop when a total displacement becomes too large
   MAX_TL_DISPLACEMENT action = ignore         ! Ignore large displacements, only print warnings
   MAX_TL_DISPLACEMENT action = refcalc 30m    ! The default. Perform a reference calculation and continue, but
                                               !   only if the last reference calculation took less than 30 minutes
   MAX_TL_DISPLACEMENT action = refcalc < 30m  ! Same as above, more explicit
   MAX_TL_DISPLACEMENT action = refcalc 2h     ! Perform new reference calculations if they take less than 2 hours.
   MAX_TL_DISPLACEMENT action = refcalc        ! Perform new reference calculations and continue, no matter how long it takes

**Acceptable values:** One or two positive floats if no flags are given.
One positive float for flags ``geo`` and ``vib``. 

Acceptable values for the ``action`` flag are:

- ``stop``: Stop the run and discard all remaining blocks in :ref:`DISPLACEMENTS`. If there are more entries in :ref:`RUN` after the search (e.g. another reference calculation), these will still be executed.
- ``ignore``: Proceed without a new reference calculation. A warning will be printed.
- ``refcalc``: Perform a reference calculation, then continue with the next block in :ref:`DISPLACEMENTS`. If a second (positive float) value is passed, this is interpreted as a time and compared to the time taken by the previous reference calculation. If the reference calculation took longer than this value, the action will be to ``stop`` instead. Values are interpreted as seconds by default, or as hours or minutes if ``h`` or ``m`` is appended to indicate the unit.

.. note:: 
    If ``MAX_TL_DISPLACEMENT`` is triggered, the ``action`` is ``refcalc`` and
    no reference calculation has so far been performed in the current run, then
    the maximum time will be ignored, and a reference calculation will be
    performed as if no time was set. The time taken for this reference
    calculation is then used if ``MAX_TL_DISPLACEMENT`` is triggered again.

.. versionadded:: 0.14.0
    This parameter is not available in versions <0.14.0. All earlier versions
    effectively behave like `MAX_TL_DISPLACEMENT action = ignore`.
