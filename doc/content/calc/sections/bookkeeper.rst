.. include:: /substitutions.rst

.. |history|        replace:: :file:`history`
.. |info|           replace:: :file:`history.info`
.. |workhistory|    replace:: :file:`workhistory`
.. |state_files|    replace:: :ref:`PARAMETERS`, :ref:`POSCAR`, and :ref:`VIBROCC`
.. |ori_files|      replace:: :file:`PARAMETERS_ori`, :file:`POSCAR_ori`, and :file:`VIBROCC_ori`
.. |SUPP/ori|       replace:: :file:`SUPP/original_inputs/`

.. _bookkeeper:

Bookkeeper
==========

The |bookkeeper| is a built-in ViPErLEED utility that automatically runs before
and after each ViPErLEED calculation.
It moves input and output files to a |history| directory, keeps the main
folder organized, and tracks all relevant changes.

See :ref:`cli_bookkeeper` for details on how to run the |bookkeeper| manually.

Most importantly, the |bookkeeper| automatically moves and updates files
|state_files|.
By default, these files are overwritten by the outputs of the previous
calculation, so that a new calculation can be started without the need
of manually copying the output files.

The |bookkeeper| has four archiving-related modes:

**Archive**
    Stores the results of the previous calculation into the |history|
    directory, and overwrites the input files |state_files| with the
    results of the previous calculation stored in the |OUT| directory.
    The previous input files are renamed to |ori_files|, respectively.
    They may be used for comparison with the non-suffixed files.
    See :ref:`history_dir` for more details on how results are stored
    in |history|.

    Runs automatically **at the end of every calculation**.

**Clear**
    Removes files belonging to a previous run.\ [#1]_ This includes all
    :file:`*_ori` files, :file:`*.log` files, as well as the |OUT| and
    |SUPP| directories.

    Runs automatically **at the start of every calculation**.

**Discard**
    Discards the results of the previous run and restores the input files to
    their original state: |ori_files| are renamed to |state_files|,
    respectively. This is useful if the previous run was not successful
    and you want to start over.
    Input and output files of the previous run will still be stored in the
    |history| directory for reference. A ``DISCARDED`` line is added to the
    ``Notes`` section of the corresponding |info| entry.

    Needs to be run manually with ``viperleed bookkeeper --discard``.

**Discard full**
    Same as **Discard**, but also removes the input and output files of the
    previous run from the |history| directory as if it had never happened.
    It also deletes any :file:`Tensors/Tensors_<index>.zip` and
    :file:`Deltas/Deltas_<index>.zip` files that were created during the
    last run, unless they had been previously used for other calculations.
    Users are required to explicitly confirm deletion of the files and folders,
    unless the ``-y`` command-line argument is specified.

    Needs to be run manually with ``viperleed bookkeeper --discard-full``.

    .. warning::
        Running |bookkeeper| in **Discard full** mode in a folder that has
        been previously **Clear**\ ed will *not restore* the input files of
        the last non-removed execution. The files can be manually copied from
        the |SUPP/ori| directory of the corresponding "main" |history| folder.


All the modes that store results to |history| (i.e., **Archive**, **Clear**,
and **Discard**) check that the |state_files| files have not been edited
after |calc| was started. If any such file is found, it is renamed by adding
an ``_edited`` suffix. |bookkeeper| will warn if any such file is found. It
is up to you to decide whether to migrate these edits to the new input files,
or whether to delete the ``*_edited`` files.

See :ref:`other_bookie_modes` for the description of non-archiving-related
|bookkeeper| modes.

When |bookkeeper| is executed in the root directory of a multi-domain
calculation, all the domain subfolders are processed in the same mode,
except for those subfolders in which |bookkeeper| was manually invoked
explicitly.

.. versionchanged:: 0.13.0
    The |bookkeeper| behavior was overhauled and the names of the modes were
    changed. See :ref:`old_bookkeeper` for details concerning differences
    with respect to earlier versions.

.. versionchanged:: 0.13.0
    |bookkeeper| automatically propagates to domain subfolders when executed
    in the root of a multi-domain calculation.

.. [#1] Before removal, results are archived to |history| if
        no |history| directory exists for the run to be cleared.


.. _history_info:

|info| and ``notes.txt``
------------------------

In addition to storing results to |history|, the |bookkeeper| also creates
or updates a plain-text |info| file. The |info| file contains one entry per
each |calc| execution that was archived in |history| (see
:numref:`hist_info_example` for an example). Entries are structured
as follows:

``TENSORS``
    Which :file:`Tensors/Tensors_<index>.zip` file was created or used in the
    |calc| run (may be multiple if reference calculations repeat). ``None`` if
    the run did not involve any :file:`Tensors` file.

    This value may not be accurate for a :ref:`fdoptimization`. In fact, no
    :file:`Tensors` file is used or produced during a :ref:`fdoptimization`.
    However, the corresponding |info| entry will use the index of the most
    recent :file:`Tensors/Tensors_<index>.zip` present in the root directory
    at the time |bookkeeper| runs.

``JOB ID``
    A progressive identifier for each run. Increments by one every time the
    |bookkeeper| archives a new calculation to |history|. Resets to 1 when a
    new :file:`Tensors` file is created, so a specific job is identified by
    the combined ``(TENSORS, JOB ID)`` pair (see also :ref:`history_dir`).
    One ``JOB ID`` is present for each of the ``TENSORS``.

``RUN``
    List of |calc| work segments that were executed in the run.

``TIME``
    Starting time of the |calc| execution that was archived to |history|.
    This serves as a reliable unique identifier, but is not as intuitive
    as the ``TENSORS`` and ``JOB ID`` numbers.

    .. versionchanged:: 0.13.0
        Changed the default format of this field to ``yyyy-mm-dd HH:mm::ss``.
        Earlier versions used ``dd.mm.yy HH:mm:ss``. |bookkeeper| will add
        new |info| entries keeping the format consistent, but will warn about
        the format change. |info| files created with earlier |bookkeeper|
        versions can be updated automatically to the new format by running
        |bookkeeper| in ``--fix`` :ref:`mode<other_bookie_modes>`.

``R REF``/``R SUPER`` (if applicable)
   |R factor| of the final reference calculation and superpos calculation,
   if any were run. Includes the total |R factor| (i.e., considering all
   beams), as well as those for integer- and fractional-order beams, if
   a distinction was specified via the :ref:`searchbeams` parameter.

``FOLDER``
    The (main) folder created for this job in the |history| directory,
    i.e., where the |bookkeeper| moved the primary output files. See
    :ref:`history_dir` for more details concerning how results are
    stored in |history|.

``Notes``
    Added manually or inserted automatically. When the |bookkeeper| runs and
    a ``notes`` or ``notes.txt`` file is present in the folder, the contents
    of that file will be automatically inserted here, and the notes file will
    be cleared. This allows taking notes concerning a currently running job
    in advance, without having to wait for it to finish and the |bookkeeper|
    to run.

    A run that has been ``--discard``\ ed has one ``DISCARDED`` line in the
    ``Notes`` field.

.. versionchanged:: 0.13.0
    Earlier |bookkeeper| versions could also create a ``JOB NAME`` line
    (after ``JOB ID``) with the string given by the user as the ``--name``
    (before v0.12.0) or ``--job-name`` (v0.12.0–v0.12.2) command-line
    argument. The command-line argument has been removed in v0.13.0, but
    |bookkeeper| will still recognize one such line.


.. _history_dir:

History organization
--------------------

Whenever |bookkeeper| archives the results of a |calc| execution to |history|,
it creates one "main" folder with name format ``tTTT.rRRR_<timestamp>``, and
may produce additional folders named ``tTTT.rRRR.SSS_<segments>_<timestamp>``
containing :ref:`workhistory`.
Each folder is labeled by its ``TENSORS`` number (``TTT``), ``JOB ID``
(``RRR``), and the starting time of the |calc| execution (``<timestamp>``).
The latter is the same for both ``tTTT.rRRR_<timestamp>``- and
``tTTT.rRRR.SSS_<segments>_<timestamp>``-named folders.

An example of the |history| folder is shown in :numref:`hist_folder_example`,
corresponding to a first |calc| execution with ``RUN = 1-3 1``.
:numref:`hist_info_example` shows the corresponding entry added
by |bookkeeper| to the |info| file.

The "main" |history| folder (``t002.r001_250321-101512`` in
:numref:`hist_folder_example`) collects both the input files (in the root of
the folder) as well as the final results of a |calc| execution (i.e., the main
:file:`*.log` file and |SUPP|/|OUT| directories).

.. _hist_folder_example:
.. code-block:: console
    :caption:
        Example of the contents of |history| following a |calc| execution
        with ``RUN = 1-3 1``.

    history/
    ├── t001.r001.001_RDS_250321-101512/  <-- after first search iteration
    │   ├── OUT/                          <-- intermediate calc results
    │   │   └── ...
    │   └── SUPP/                         <-- intermediate calc results
    │       └── ...
    │
    ├── t001.r001.002_DS_250321-101512/  <-- after second search iteration
    │   ├── OUT/                         <-- intermediate calc results
    │   │   └── ...
    │   └── SUPP/                        <-- intermediate calc results
    │       └── ...
    │
    └── t002.r001_250321-101512/          <-- "main" history folder
        ├── SUPP/                         <-- calc results at end
        │   ├── original_inputs/
        │   │   ├── EXPBEAMS.csv
        │   │   ├── PARAMETERS
        │   │   ├── PHASESHIFTS
        │   │   └── POSCAR
        │   └── ...
        ├── OUT/                          <-- calc results at end
        │   └── ...
        ├── DISPLACEMENTS_from_root       <-- not found in original_inputs
        ├── EXPBEAMS.csv                  <-- original input file
        ├── PARAMETERS                    <-- original input file
        ├── PHASESHIFTS                   <-- original input file
        ├── POSCAR                        <-- original input file
        └── viperleed-calc_250321-101512.log


The input files stored in a "main" |history| folder are preferentially
collected from the |SUPP/ori| directory. If an input file is not found
in |SUPP/ori|, it is copied from the directory in which |bookkeeper|
executes and renamed by adding a ``_from_root`` suffix. This is meant
to prevent obscuring potential user edits to the root files since |calc|
was started.

.. note::
    |bookkeeper| will archive in the main |history| directory all the
    input files, irrespective of whether they were actually used during
    the corresponding calculation.

.. _hist_info_example:
.. code-block:: text
    :caption:
        Contents of the |info| file corresponding to the |history| of
        :numref:`hist_folder_example`.

    # TENSORS    1, 2
    # JOB ID     1, 1
    # RUN        0 1 11 2 3 31 12 2 3 31 12 1 11
    # TIME       2025-03-21 10:15:12
    # R REF      <R factor after last refcalc>
    # R SUPER    <R factor after last superpos>
    # FOLDER     t002.r001_250321-101512
    Notes:  ...user notes...

    ###########


.. _workhistory:

Intermediate |calc| results
***************************

A similar system as the |bookkeeper|'s history is used by |calc| during
execution. When the order of execution loops "backwards" — i.e., when
executing another reference calculation after a search, or running
multiple searches consecutively — |calc| creates a |workhistory|
directory and moves there a snapshot of all relevant outputs.

When the |bookkeeper| runs and a |workhistory| directory is present, the
|workhistory| contents will be incorporated into the |history| folder.
Directories in the |workhistory| will be moved to |history| and renamed to
``tTTT.rRRR.SSS_<segments>_<timestamp>``, following the same basic formatting
of the "main" directories in |history|. ``SSS`` is an incremental number for
|workhistory| directories generated during the same run, while ``<segments>``
gives a quick overview about what segments were performed before this snapshot
was created: ``R``\ eference calculation, ``D``\ elta amplitudes,
or ``S``\ earch.

:numref:`hist_folder_example` shows an example of the |history| folder
after executing |calc| for the first time with ``RUN = 1-3 1`` and with
a :ref:`DISPLACEMENTS` file that defines two consecutive search blocks.
When the first search is finished, and the job loops around to execute
further delta-amplitudes calculations and searches, |calc| collects in a
first |workhistory| folder the intermediate output files. This snapshot
appears as folder :file:`t001.r001.001_RDS_250321-101512` in |history|:
``.001`` marks that it is the first such snapshot taken for job ``.r001``,
and ``RDS`` indicates that a reference calculation, a delta-amplitudes
calculation, and a search were executed before the snapshot was taken.
The reference calculation generated a :file:`Tensors/Tensors_001.zip`
file (``t001``). When |calc| concludes the second search block and "loops
back" to run the last reference calculation, a second snapshot is taken,
yielding folder :file:`t001.r001.002_DS_250321-101512` in |history|: this
second (``.002``) intermediate snapshot used the same :file:`Tensors_001.zip`
file as the first one (``t001``) but only executed delta-amplitudes and
search segments (``DS``).

The last reference calculation produces a new a :file:`Tensors/Tensors_002.zip`
file. Its results constitute the final output of this |calc| execution, and
are stored in the "main" history folder, named :file:`t002.r001_250321-101512`.

.. versionchanged:: 0.13.0
    Earlier versions would assign the incorrect ``JOB ID`` to intermediate
    results of a |calc| execution with ``RUN = 2-3 1`` (and similar). See
    `this <https://github.com/viperleed/viperleed/pull/198#issuecomment-2716696302>`__
    comment for details.



.. _other_bookie_modes:

Other |bookkeeper| modes
------------------------

This section describes non-archiving-related modes for the |bookkeeper|
utility. They require manual invocation of the |bookkeeper| (see also
:ref:`cli_bookkeeper`).

**Fix**
    Edits the |history| folder and the |info| file to conform to the most
    recent format implemented in the |bookkeeper|, where possible. A list
    of the automatic edits that mode **Fix** can apply to the |info| file
    as of v0.13.0 can be found
    `here <https://github.com/viperleed/viperleed/pull/198#issuecomment-2403901360>`__.
    The only modification currently implemented for |history| folders consists
    in the addition of metadata (file :file:`bookkeeper_meta`) used internally
    by |bookkeeper| for cross-referencing |history| folders.


.. _old_bookkeeper:

Previous versions of |bookkeeper|
---------------------------------

This section outlines differences that have been introduced in |bookkeeper|
as of v0.13.0 compared to previous versions, and documents related changes
in |calc|. |bookkeeper| has been kept backward compatible as much as possible.
However, some breaking changes have been introduced in |calc| that make running
|bookkeeper| in a folder that was created with earlier versions of |calc| (or
|bookkeeper|) not entirely consistent with the results of executing
|bookkeeper| in a folder created with v0.13.0 and later versions.

The most important changes introduced in v0.13.0 concern the storage
of user-given input files by |calc|, the behavior of ``--discard-full``
mode, and the modification of the default interaction between |calc|
and |bookkeeper|.

.. _ori_inputs_v_0_13:

Storage of user-given input files
*********************************
While the user-given input files have been stored in |SUPP/ori| since
|calc| v0.7.0, storage would happen too late. As a result, the files
saved to |SUPP/ori| had already been processed by |calc|. See
`this <https://github.com/viperleed/viperleed/pull/198#issuecomment-2623894755>`__
and the following comments for more details.

Since |bookkeeper| v0.13.0 (and later versions) relies on the "originality"
of the files in |SUPP/ori|, executing it in a tree created by earlier |calc|
versions may potentially (i) archive in |history| input files that were already
processed by |calc|, and (ii) restore the incorrect |state_files| files when
executed in **Discard**/**Discard Full** modes. Especially in the latter case,
we suggest to carefully validate that the contents of the |state_files| input
files are as intended.

.. _discard_full_backwards_compatibility:

Behavior of ``--discard-full`` mode
***********************************
|bookkeeper| v0.13.0 introduced a mechanism to cross-reference the "main"
|history| folder to intermediate results (see :ref:`history_dir` for details).
Before running |bookkeeper| v0.13.0 (and later) in ``--discard-full`` mode in
a tree created by earlier |bookkeeper| versions, make sure to manually invoke
|bookkeeper| in ``--fix`` :ref:`mode <other_bookie_modes>`. This will add the
necessary cross-reference information where possible. Failing to execute
``bookkeeper --fix`` before ``--discard-full`` will only delete from |history|
the "main" folder of the most recent run, and remove from |info| the
corresponding entry. This will leave intermediate-results folders untouched.
Additionally, these folders will not have a corresponding |info| entry.

Interaction with |calc| v0.1.0 – v0.11.0
****************************************
Up to v0.11.0, |bookkeeper| would automatically execute only before a
|calc| run. It would, at that point, archive to |history|/|info| any
results of earlier |calc| executions. Then, it would remove any outputs
of the previous run from the root folder.
However, at that moment input files could have been modified relative to
the ones originally given when the previous |calc| execution started.
Such edited input files would then be archived to |history|.

The intended workflow for v0.1.0 – v0.11.0 that would give the correct
storage of files to |history| was:

  1. Execute |calc|.
  2. Manually execute |bookkeeper| in either "default" (i.e., the current
     **Clear** mode), ``--discard`` (i.e., the current **Discard Full** mode),
     or ``--cont`` modes (i.e., the current **Archive** mode). This manual
     execution would need to be done *before any edit to input files*.
  3. Start a new |calc| session.

Failing to manually execute |bookkeeper| after |calc| would cause the following
|calc| execution to always start from the same inputs as the previous one.
Manually executing |bookkeeper| after having introduced any edit to the input
files would cause storage of the incorrect input files to |history|.
|bookkeeper| v0.13.0 cannot detect the latter scenario.

Interaction with |calc| v0.12.0 – v0.12.2
*****************************************
In v0.12, |bookkeeper| would still automatically execute before a |calc| run
as in previous versions. Additionally, it would also run in ``--cont`` mode
**after** a |calc| run.

The rationale behind this choice is that it was easy to forget to manually
execute |bookkeeper| after |calc| and *before editing the input files* in
the root directory. As discussed in the previous sections, this would then
likely store the incorrect files to |history|.

This choice, however, turned out to be fairly poor for a number of reasons:

  1. The results of a calculation would be "hidden away" inside the |history|
     subfolder, rather than being readily available in the folder where |calc|
     was started.
  2. The |state_files| files in the root directory after execution would
     be the ones *at the end*, potentially after a structure optimization.
     This would make comparison of the progress of a calculation not
     straightforward.
  3. Manually running |bookkeeper| would have no effect after executing
     |calc| with default arguments, irrespective of the |bookkeeper| mode.

Executing |bookkeeper| v0.13.0 (and later versions) in a tree created by
versions 0.12.0–0.12.2 will not yield the expected results in modes
``--discard`` and ``--discard-full``: the input files will not be
restored to those at the beginning of the discarded run. These files
can, however, be manually copied from the corresponding |history| folder.
