.. include:: /substitutions.rst

.. _suppress_exe:

SUPPRESS_EXECUTION
==================

If SUPPRESS_EXECUTION is on, input files for TensErLEED will be written, but 
|calc| will not attempt to compile the TensErLEED code or execute any part of 
TensErLEED.

**Default**: SUPPRESS_EXECUTION = False (Compile and execute TensErLEED 
automatically)

**Syntax**:

::

   SUPPRESS_EXECUTION = True

**Acceptable values**: True, False, T, F (not case sensitive)

Note that SUPPRESS_EXECUTION will force the program to stop after each segment 
(e.g., after starting the reference calculation). This also means that output 
will be present in the form that TensErLEED outputs, and not automatically 
processed for easy analysis.
