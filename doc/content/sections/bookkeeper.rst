.. _bookkeeper:

Bookkeeper
==========

The bookkeeper is a small application made to run in between executions 
of ViPErLEED calculations. Its purpose is to collect input and 
output files from the previous run and store them in a "history" 
subfolder for future reference.

The bookkeeper can safely be run multiple times. If no new output 
is detected, it will simply exit without doing anything.

The bookkeeper may be called with the following arguments:

-  ``-c`` / ``--cont``: Continuation job. After storing input and output in the history, the bookkeeper will overwrite the POSCAR and VIBROCC files in the main folder with the latest POSCAR_OUT and VIBROCC_OUT files.
-  ``-h`` / ``--help``: Displays usage information.
-  ``-n`` / ``--name``: Defines a string to be appended to the directory name in the history folder being created. The name will also be stored in the history.info file.

.. _history_info:

history.info and notes.txt
--------------------------

The bookkeeper also creates or updates a ``history.info`` file, which stores the following information each time 
the bookkeeper runs:

-  **Notes**: Either added manually, or inserted automatically. When the bookkeeper runs and a ``notes`` or ``notes.txt`` file is present in the folder, the contents of that file will automatically be inserted here, and the notes file will be cleared. This allows taking notes concerning a currently running job in advance, without having to wait for it to finish and the bookkeeper to run.
-  **Tensors**: Which Tensors file was used in the tleedm job being stored (may be multiple if reference calculations repeat)
-  **Job ID**: A number given to the job being stored. Increments by 1 every time the bookkeeper runs. Resets to 1 when a new Tensors file is created, so a specific job is identified by the combined TensorNumber.JobID (see below)
-  **Run**: List of tleedm work segments that were run in the job. Taken from the end of the main tleedm log.
-  **R-factors** (if applicable): R-factors of the final reference calculation and/or superpos calculation, if any were run.
-  **Time**: Starting time of the tleedm job being moved. This serves as a reliable unique identifier, but is not as intuitive as the Tensors and job numbers.
-  **Folder**: The (main) folder created for this job in the "history" directory, i.e. where the bookkeeper is moving files to.


History organization
--------------------

History directory names are formatted as t\ ``XXX``.r\ ``YYY``\ \_\ ``timestamp``, where ``XXX`` is the (highest) Tensor number of the job, ``YYY`` is the Job ID, and the timestamp corresponds to the starting time of tleedm execution.

A similar system as the bookkeeper history is used by tleedm during 
execution when the order of execution loops backwards, i.e. by executing 
another reference calculation after a search, or running multiple 
searches consecutively.
In these cases, tleedm creates a "workhistory" directory and moves a 
snapshot of all input and output files that *may* be relevant and *may* 
get overwritten into a subfolder there. Note that this process may 
capture input and output files that are not from the latest segments; 
the purpose of this is to preserve information, but it does not reliably 
sort files by when they were created.

When the bookkeeper runs and a "workhistory" directory is present, the 
workhistory contents will be incorporated into the history.
Directories in the workhistory will be moved to history and renamed to 
t\ ``XXX``.r\ ``YYY``.\ ``ZZZ``\ \_\ ``segments``\ \_\ ``timestamp``, 
following the same basic formatting of the job directories in the 
history.
``XXX``, ``YYY`` and ``timestamp`` correspond to the Tensor number, Job 
ID and timestamp of the full job; final output will be in the directory 
t\ ``XXX``.r\ ``YYY``\ \_\ ``timestamp``. The ``ZZZ``\ \_\ ``segments`` 
directories store the snapshots at given times of executions, i.e. 
before the job returned to a segment that was already executed, or that 
may overwrite input.
``ZZZ`` simply numbers the different workhistory directories, while 
``segments`` gives a quick overview about what type of calculations were 
performed before this snapshot was created.

For example, if a first job on a system is run with ``RUN = 1-3``, and 
the DISPLACEMENTS file defines multiple consecutive searches, then the 
final history directory will be named t001.r001\_\ ``time``.
When the first search is finished, and the job loops around to execute 
further delta-amplitudes calculations and searches, output may be 
overwritten. 
Therefore, at this point, a snapshot of all input and output files is 
created and stored in workhistory; in the final history folder, this 
snapshot will appear as t001.r001.001_RDS\_\ ``time``. The ".001" marks 
that it is the first such snapshot taken for this job, and the "RDS"
indicates that a Reference calculation, Delta calculation, 
and Search were executed before the snapshot was taken.
If, after the second delta/search execution, the job loops back again, 
then a second snapshot will be created and end up in t001.r001.002_DS\_\ 
``time``. Again, if you are only interested in the final output, you 
can safely ignore all history directories following the t\ ``XXX``.r\ 
``YYY``.\ ``ZZZ``\ \_\ ``segments``\ \_\ ``timestamp`` formatting, and 
only check the t\ ``XXX``.r\ ``YYY``\ \_\ ``timestamp`` directory.
