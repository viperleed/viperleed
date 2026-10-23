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
The :term:`TensErLEED` source code should be available on the testing machine
(repository ``viperleed-tensorleed``) and the ``VIPERLEED_TENSORLEED``
:term:`environment variable` should be set. See also :ref:`install_tensorleed`.

Also make sure to have compiled the :ref:`eeasisss_compile` source code.

Tests are *not* distributed with the ``viperleed`` package on
`pypi <https://pypi.org/project/viperleed/>`__. To run the tests, clone the
|repo|, then install the test-related dependencies via

.. code-block:: bash

    python -m pip install --upgrade -r requirements/tests.txt

To execute tests, navigate to your local copy of the |repo|,
then run

.. code-block:: bash

    python -m pytest .

.. important::
    Notice that this will execute the tests against the version of ``viperleed``
    that you have currently **installed** on your system. If you are working on
    developing the code, and would like to test the most recent changes, make
    sure to first install ``viperleed`` in *editable* mode in your environment.
    See `editable_install`.

If you also would like to produce coverage reports, install :program:`coverage`
via

.. code-block:: bash

    python -m pip install coverage

then, in your local copy of the ``viperleed`` repository, run

.. code-block:: bash

    python -m coverage run -m pytest .
    python -m coverage html


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

Independently of whether you use the hook, pull requests get an advisory
check that flags out-of-date end years in the files they change. This
check does not block merging.
