.. _keep_refcalc_dirs:

KEEP_REFCALC_DIRS
=================

KEEP_REFCALC_DIRS is a boolean parameter that toggles whether the
subfolders in which the :ref:`reference calculation<ref-calc>`
for each energy step is executed are kept after the section finishes.
These directories usually contain little useful additional information
and are discarded by default.
The :ref:`reference calculation<ref-calc>` logs are instead collected in the
:ref:`refercence calculation log file<log_files_refcalc>`. However, the user
may want to keep the directories for debugging purposes.

**Default**: KEEP_REFCALC_DIRS = False

**Acceptable values**: ``True``, ``False``