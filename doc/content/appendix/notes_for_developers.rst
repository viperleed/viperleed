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

To execute tests, navigate to the :file:`tests` folder of your local copy
of the |repo|, then run

.. code-block:: bash

    python -m pytest .

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

    python -m pytest . --quick        # quick tier
    python -m pytest .                # default tier
    python -m pytest . -m tensorleed  # TensErLEED tier
    python -m pytest . -m ""          # everything

The quick tier leaves out tests that take about half a second or more (marked
``slow``). The ``--quick`` option adds ``and not slow`` to the ``-m``
expression in effect: it is equivalent to ``-m "not slow and not tensorleed"``,
and ``-m tensorleed --quick`` runs the TensErLEED tests that are not slow.

The TensErLEED tier requires ``VIPERLEED_TENSORLEED``: its tests fail (rather
than being skipped) when :term:`TensErLEED` or one of the versions they need
is missing. Optional dependencies that are not installed make some tests skip.
To turn those skips into failures, add the option ``--require-all``.

Coverage is not collected by default. To produce a coverage report, run

.. code-block:: bash

    python -m pytest . --cov --cov-config=../pyproject.toml

(the :program:`pytest-cov` plugin is among the test requirements). What is
measured (including branches) is configured in :file:`pyproject.toml`; the
option ``--cov-config`` is not needed when running from the repository root
(``python -m pytest tests --cov``). A summary with the missing lines is
printed, and an HTML report is written to :file:`htmlcov/`.

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
quiet the machine, the better: timings are noisy):

.. code-block:: bash

    python -m pytest . --check-slow -n 4
    python -m pytest calc/symmetry --check-slow

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
