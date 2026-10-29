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
a default run. To run them, the :term:`TensErLEED` source code (repository
``viperleed-tensorleed``) should be available on the testing machine and the
``VIPERLEED_TENSORLEED`` :term:`environment variable` should be set. See also
:ref:`install_tensorleed`, which also describes where the compiled
:ref:`eeasisss_compile` source code is expected.

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
``slow`` but whose call took the threshold (``SLOW_THRESHOLD`` in
:file:`tests/tiers.py`) or more, and those marked ``slow`` whose call took
less. It only reports: it changes neither the selection nor the exit status.
It also works with ``-n`` (:program:`pytest-xdist`).


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
