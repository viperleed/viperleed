.. _halting:

HALTING
=======

HALTING defines under which conditions execution should stop after a segment
(e.g., after initialization or refcalc). A low value means execution will be
stopped more readily.

**Default**: HALTING = 2 (stop for user confirmation of intermediate results)

**Syntax**:

::

   HALTING = 1

**Acceptable values**:

-  1 ... stop for minor warnings
-  2 ... stop for confirmation of some intermediate results, e.g., when an
   AUXPOSCAR file is first generated, or major hiccups that suggest that
   something has gone seriously wrong
-  3 ... only stop if program crashes, or required input or data is missing

The default value of 2 should be sensible for most users. A value of 1 might
make sense if the user is sceptical about the input. A value of 3 is **not**
recommended, but can be applied by experienced users who are very sure of
their input being correct, or who don't care about wasting computational
resources if it isn't.
