.. _notes_for_developers:

.. |repo| replace:: ViPErLEED GitHub `repository <https://github.com//viperleed/viperleed>`__


====================
Notes for developers
====================

Make sure to work in a suitable virtual environment. See :ref:`use_venv`
for more information.


Running tests
-------------

ViPErLEED comes with a suite of tests and uses
`pytest <https://docs.pytest.org/en/>`__ for testing.
Tests that run :term:`TensErLEED` (marked ``tensorleed``) are not part of
a default run. To run them, a :term:`TensErLEED` tree should be available on
the testing machine and the ``VIPERLEED_TENSORLEED``
:term:`environment variable` should point to it. The tests need several
versions of :term:`TensErLEED` and a compiled :program:`eeasisss` (see also
:ref:`install_tensorleed` and :ref:`eeasisss_compile`). Assemble the tree
with the script :file:`build/scripts/assemble_tensorleed_tree.py`, from the
root of your local copy of the |repo|:

.. code-block:: bash

    python build/scripts/assemble_tensorleed_tree.py <path/to/tree>

The script downloads release tags of the ``viperleed-tensorleed`` repository,
copies the :term:`TensErLEED` versions used by the tests into the given
directory, compiles :file:`random_.o` for :term:`TensErLEED` v1.6.1 (from
the cross-platform :download:`random_.c</_static/tenserleed_windows/random_.c>`
of :ref:`install_tensorleed`, on all platforms), and compiles
:program:`eeasisss`. It needs Internet access, :program:`gcc`,
:program:`gfortran`, and :program:`make` (on Windows, the
:file:`make.bat` file shipped with :program:`eeasisss` is used instead;
the script has been tested on Linux only). For the Intel compilers, add
``--compilers intel`` (:program:`icx` and the ``intel`` target of
:program:`eeasisss`); ``--cc`` and ``--fc`` choose other C and Fortran
compilers (``--fc`` not on Windows). See
``python build/scripts/assemble_tensorleed_tree.py --help``. Only
:program:`gcc` and :program:`gfortran` are tested in continuous
integration. The directory must not exist or
be empty: the script does not overwrite an existing tree. If a step fails,
the script removes what it wrote, so that it can be run again with the same
directory. To update the tree after a change of the script, assemble it
into a new directory. At the end, the script warns if
``viperleed-tensorleed`` has release tags newer than those it uses (it
does not fail, also not if the tags cannot be listed): the script may then
need an update.

Tests are *not* distributed with the ``viperleed`` package on
`pypi <https://pypi.org/project/viperleed/>`__. To run the tests, clone the
|repo|, then install the test-related dependencies via

.. code-block:: bash

    python -m pip install --upgrade -r requirements/tests.txt

To execute tests, run the following from the root of your local copy
of the |repo|:

.. code-block:: bash

    python -m pytest tests

.. important::
    Notice that this will execute the tests against the version of ``viperleed``
    that you have currently **installed** on your system. If you are working on
    developing the code, and would like to test the most recent changes, make
    sure to first install ``viperleed`` in *editable* mode in your environment.
    See `editable_install`.

The command above runs the *default* tier: all tests except those that need
:term:`TensErLEED`. Tests are grouped in tiers via markers and the ``-m``
and ``--quick`` options:

.. code-block:: bash

    python -m pytest tests --quick        # quick tier
    python -m pytest tests                # default tier
    python -m pytest tests -m tensorleed  # TensErLEED tier
    python -m pytest tests -m "" --release  # everything

The quick tier leaves out tests that take about half a second or more (marked
``slow``). The ``--quick`` option adds ``and not slow`` to the ``-m``
expression in effect: it is equivalent to
``-m "not slow and not tensorleed and not release"``, and
``-m tensorleed --quick`` runs the TensErLEED tests that are not slow.

The TensErLEED tier requires ``VIPERLEED_TENSORLEED``: its tests fail (rather
than being skipped) when :term:`TensErLEED` or one of the versions they need
is missing. Optional dependencies that are not installed make some tests skip.
To turn those skips into failures, add the option ``--require-all``.

A few TensErLEED tests are in a fourth tier, the *release* tier (marker
``release``). They are deselected by default, also under ``-m tensorleed``,
and only run with the option ``--release``:

.. code-block:: bash

    python -m pytest tests -m tensorleed --release  # TensErLEED + release

Tests go in the release tier if they take very long and add little over the
other TensErLEED tests (currently, the calculation of the phase shifts of the
:math:`\mathrm{Al_2O_3}/\mathrm{NiAl}(111)` slab with a core hole). The cases
are listed in :file:`tests/tiers.py` (``RELEASE_CASES``) by fixture
(``run_phaseshift``) and case id, so that every test that uses the fixture
for that slab is a release test. Entries that match no case of a collected
test are reported as stale, like those of ``SLOW_CASES``. Naming a release
case by node id or with ``-m release`` also needs ``--release``; otherwise
nothing runs (a serial run says why). The option
combines with ``-n``, ``--quick`` and ``-m``: ``--release`` only removes the
deselection of release tests, so ``-m tensorleed`` without it remains the
command for the TensErLEED tier.

The TensErLEED tier also contains a full chain for :term:`TensErLEED` v1.6.1
and for the most recent version
(:file:`tests/calc/sections/test_chain.py`). A single fixture per version runs
initialization, reference calculation, delta amplitudes, and search for
Ag(100) in one go (``RUN = 0 1 2 3``), without any stored tensors or delta
files, so that the code paths for versions older than 1.7.0, and those for
newer versions, are exercised in all sections:

.. code-block:: bash

    python -m pytest tests/calc/sections/test_chain.py -m tensorleed --no-cov

The default tier also contains *stored-input* tests
(:file:`tests/calc/stored_inputs`). They need neither :term:`TensErLEED` nor
a compiler: starting from stored :file:`POSCAR`, :file:`PARAMETERS` and
related files, they run only the step of ViPErLEED that writes the input of
a :term:`TensErLEED` program (without executing it), and compare it with
reference files in :file:`tests/_test_data/stored_inputs/<step>/<case>/`.
The references cover the input of :program:`eeasisss`, the reference
calculation (:file:`refcalc-FIN`, :file:`refcalc-PARAM`, :file:`AUX*`,
:file:`BEAMLIST`), the delta amplitudes (:file:`delta-input`,
:file:`Delta_Compile_*-PARAM`), and the search (:file:`PARAM`,
:file:`restrict.f`, :file:`search.steu`, :file:`searchpars.info`).
Except for :program:`eeasisss`, each case has one subdirectory for the files
written for :term:`TensErLEED` 1.6.0--1.6.1 and one for 1.7.3 and later
(``1.6.0-1.6.1``, ``1.7.3-``). The :file:`PHASESHIFTS` block of
:file:`refcalc-FIN` is replaced by a single placeholder line. The tests of the
``1.7.3-`` group run at both ends of it (:term:`TensErLEED` 1.7.4 and 2.1.0)
against the same references. The tests pin timestamps and paths (the
timestamp, and the relative path of the :file:`seSernelius` file), and the
files are compared verbatim, so the references hold them as written.

Rewrite the references only if what ViPErLEED writes was intentionally
changed, or when adding a group for a new :term:`TensErLEED` version. From the
root of your local copy of the |repo|, run

.. code-block:: bash

    python -m pytest tests/calc/stored_inputs --rewrite-references

Tests whose reference differs, or is missing, rewrite it and are reported as
skipped (``reference rewritten: <path>``), and the run exits with a
non-zero status, so that a rewrite is never mistaken for a pass. Review the
change with ``git diff tests/_test_data/stored_inputs``, and then run the tests
again without the option. Do not rewrite references to make an unexplained
failure pass.

Coverage (including branches) is collected by default, via the
:program:`pytest-cov` plugin, which is among the test requirements. An HTML
report is written to :file:`htmlcov/`; add ``--cov-report=term`` to also get
a summary with the missing lines in the terminal. Coverage makes tests about
2--3 times slower. If you do not need it, turn it off with ``--no-cov``, which
is especially worthwhile for quick runs:

.. code-block:: bash

    python -m pytest tests --quick --no-cov

Slow tests are kept up to date with the option ``--check-slow``. The slow
tests that are individual cases of a parametrized test (as opposed to whole
tests or classes marked ``slow``) are listed in :file:`tests/tiers.py`
(``SLOW_CASES``): it maps each test to the parts of the case ids that make it
slow. At collection, an entry that matches no test or case (e.g., after
renaming a test or changing its parameters) is an error naming it; update the
module. Entries of modules that are not collected are not checked, nor are
runs that select tests by node id or with ``--lf``, so a subset of the tests
can still be run. To find tests that should be marked, or
that got fast and should not be, add ``--check-slow`` to a run (the more
quiet the machine, the better: timings are noisy). Coverage distorts the
durations, so ``--check-slow`` needs ``--no-cov``:

.. code-block:: bash

    python -m pytest tests --check-slow --no-cov -n 4
    python -m pytest tests/calc/symmetry --check-slow --no-cov

At the end of the run, the option lists the tests that are not marked
``slow`` but whose call took 1.2 times the threshold (``SLOW_THRESHOLD`` in
:file:`tests/tiers.py`) or more, and those marked ``slow`` whose call took
less than 0.8 times the threshold. The band around the threshold (the
``SLOW_BAND`` constant in the same module) is there because the time of a
test depends on the cases that ran before it; tests in the band keep their
marking. It only reports: it changes neither the selection nor the exit
status. It also works with ``-n`` (:program:`pytest-xdist`).


Continuous integration
''''''''''''''''''''''

The tests also run on GitHub Actions:

- The workflow ``Tests`` (:file:`.github/workflows/unit_tests.yml`) runs the
  default tier with ``-n auto`` and ``--require-all`` on every pull request
  and on every push to the ``master`` branch, on two Python versions (see the
  workflow file). Its result is the check ``Tests``. The coverage of the run
  with the newer Python is uploaded to
  `Codecov <https://codecov.io/gh/viperleed/viperleed>`__. After every push
  to ``master``, the job ``Tests (serial)`` also runs the default tier
  without :program:`pytest-xdist`, to catch tests that fail only when run
  serially. It is not part of the check ``Tests``: a pull request that
  breaks the serial run is noticed only after it is merged.
- The *TensErLEED tests* (:file:`.github/workflows/tensorleed_tests.yml`)
  run the TensErLEED tier on every pull request into ``master``, against a
  tree assembled with :file:`build/scripts/assemble_tensorleed_tree.py`.
  On pull requests, they include the release tests (``--release``, see
  above).
  Their result is the check ``TensErLEED tests``. The workflow does not run
  on pull requests into other branches. It starts when a pull request into
  ``master`` is opened or receives new commits: if the base of a pull
  request is changed to ``master`` later, push a commit or start the
  workflow by hand (see below). The tree is cached, but a cache made for a
  pull request is visible only to that pull request: the tree is assembled
  on the first run of each pull request, and later runs reuse it.

Whether a pull request can be merged with a failing check depends on the
rulesets of the |repo|, which the administrators of the repository set up.

To run the TensErLEED tests for another branch, e.g., a version branch
before opening the pull request into ``master``, start the workflow by hand:
on the *Actions* page of the |repo|, select ``TensErLEED tests``, click
:guilabel:`Run workflow`, and choose the branch. With the
`GitHub CLI <https://cli.github.com>`__, run

.. code-block:: bash

    gh workflow run tensorleed_tests.yml --ref <branch>

Runs started by hand skip the release tests unless the input ``release`` is
checked (with the GitHub CLI, add ``-f release=true``).


.. _editable_install:

Using editable installations
----------------------------

ViPErLEED is installed from source using the code in the |repo|. It
uses the :file:`pyproject.toml` file. To install an
`editable version <https://setuptools.pypa.io/en/latest/userguide/development_mode.html>`__,
navigate to your copy of the :file:`viperleed` repository, and use

.. code-block:: bash

    python -m pip install -e .[<options>]

Editable installations from :file:`pyproject.toml` files
require ``pip>=21.3``. Update your ``pip`` with

.. code-block:: bash

    python -m pip install --upgrade pip


Building this documentation
---------------------------

The source for this documentation is *not* distributed with the ``viperleed``
package on `pypi <https://pypi.org/project/viperleed/>`__. The sources of
the documentation are available in the |repo|. To build the documentation,
you will need to install additional dependencies. You can install them by
running

.. code-block:: bash

    python -m pip install --upgrade -r requirements/tests.txt

in your local copy of the |repo|.

Note that the documentation can only be built with :program:`Python` ≥3.9
because of dependency-resolution issues.

Navigate to the :file:`doc` subfolder of :file:`viperleed`, then

.. code-block:: bash

    make html

or

.. code-block:: bash

    make latexpdf

Producing the PDF documentation requires a working LaTeX
installation on your system.


Building ``viperleed`` for distribution
---------------------------------------

Install the distribution dependencies by running

.. code-block:: bash

    python -m pip install --upgrade -r requirements/dist.txt

in your local copy of the |repo|.

PyPi
''''

Follow the instructions under
`<https://packaging.python.org/en/latest/tutorials/packaging-projects/>`__.

.. todo:: be more specific


.. todo:: @amimre add instructions for pyinstaller


Installing all development dependencies
---------------------------------------

Development of ViPErLEED requires a few more dependencies than those in the
distribution version of the package. You can install all of them at once by
running

.. code-block:: bash

    python -m pip install --upgrade -r requirements/dev.txt

in your local copy of the |repo|.


Keeping copyright years up to date
----------------------------------

The |repo| comes with a `pre-commit <https://pre-commit.com>`__ hook that
sets the end year of the ``__copyright__`` module dunder of every Python file
you commit to the current year, if it is out of date. It only updates files
that are staged for the commit, so ``pre-commit run --all-files`` leaves
files you did not touch alone. Using it is optional. Like pre-commit
itself, the hook needs ``git`` to be available on your ``PATH``, which
may not be the case by default on Windows. To opt in, run

.. code-block:: bash

    python -m pip install pre-commit
    pre-commit install

once in your local copy of the |repo|. When the hook updates a file, the
commit is aborted and the updated files are listed. Stage them again (e.g.,
with ``git add``) and repeat the commit. If pre-commit instead reports that
it is rolling back fixes, a file you staged only partially has unstaged
changes next to its ``__copyright__`` dunder: stage or stash those changes
first, then repeat the commit.

The same set-up also tidies every file you commit: it strips trailing
whitespace, adds missing final newlines and unifies mixed line endings
(data files, e.g., for tests, are left alone), aborting the commit like the
copyright hook does. It also runs `ruff <https://docs.astral.sh/ruff/>`__
to reject syntax errors (including syntax newer than Python 3.7) and
undefined names.

Independently of whether you use the hook, pull requests get an advisory
check that flags out-of-date end years in the files they change. This
check does not block merging.
