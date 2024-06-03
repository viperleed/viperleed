.. _stop:

.. include:: /substitutions.rst

STOP
====

The PARAMETERS file in the work folder will be regularly checked for updates. 
If STOP is set to True, the program will cancel ongoing delta calculations or 
searches, and stop after the current section.

**Default**: STOP = False

**Syntax**:

::

   STOP = True

or simply

::

   STOP

**Acceptable values**: no value required; simply starting a line with STOP will 
set it to True (unless the line is ``STOP = False``).

STOP provides a way to cleanly terminate the program when a keyboard interrupt 
is not possible. Note that if you are using a job script to execute |calc| in 
a work folder, you need to modify the PARAMETERS file in that work folder, not 
in the home folder, for STOP to take effect.

If STOP is already set when the program starts, the PARAMETERS file will be 
modified to comment out STOP, and a warning will be printed.

Even if STOP is set, the superpos and R-factor calculation sections will still 
run if they are next in line.
