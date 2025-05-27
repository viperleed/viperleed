.. include:: /substitutions.rst

.. _installation:

============
Installation
============

This section provides a short guide on how to install required
components of ViPErLEED.

Installing the ViPErLEED ImageJ plugins
=======================================

The installation instructions for ImageJ and the ViPErLEED plugins can
be found in the :ref:`imagej_plugins` section.


Installing Python
=================

With the exception of the :ref:`imagej_plugins`, ViPErLEED requires a working
Python installation. You can install Python using a package manager or as
binary from the `Python website <https://www.python.org/downloads/>`__.

.. note::

    ViPErLEED requires Python 3.7 or newer, though we recommend using the
    latest available Python version.

We suggest Windows users to follow the notes in the
`official Python documentation <https://docs.python.org/3/using/windows.html#using-python-on-windows>`__:
use the full installer, and include the
`Python Launcher for Windows <https://docs.python.org/3/using/windows.html#python-launcher-for-windows>`__.


Installing the ViPErLEED Python package
=======================================

ViPErLEED is available on the :term:`Python Package Index (PyPI)<PyPI>`
and can be installed using the Python package manager ``pip``. We recommend
installing ViPErLEED in a
`virtual environment <https://docs.python.org/3/library/venv.html>`__
for easier dependency management. See also :ref:`use_venv`.

To install ViPErLEED from :term:`PyPI` using ``pip``, run

.. code-block:: console

    python -m pip install viperleed

in your terminal.\ [1]_

This will install the latest version of ViPErLEED and all required
dependencies. It will also automatically install the
:ref:`ViPErLEED command-line tools<command_line_tools>` which can be called
from the terminal using the ``viperleed`` command. The |calc| package for
|LEED-IV| calculation has also some non-Python dependencies that can be
installed as described in :ref:`install_tensor_leed_deps`.

If you want to also install the dependencies for running the ViPErLEED
graphical user interface, which also allows you to :ref:`measure<hardware>`
|LEED-IV| data and to simulate LEED patterns, run

.. code-block:: console

    python -m pip install viperleed[GUI]

in your terminal.

.. note::

    If you are using the :program:`zsh` shell (default on macOS and some Linux
    distributions), you need to quote or escape the square brackets as they are
    otherwise interpreted as a glob pattern:

    .. code-block:: bash

        python -m pip install "viperleed[GUI]"
        # or
        python -m pip install viperleed\[GUI\]

.. _wsl:

Windows Subsystem for Linux
===========================

To run |calc| under Windows, we recommend using the
:term:`Windows Subsystem for Linux<WSL>` (WSL), available starting from
Windows 10. Running |calc| natively on Windows is possible, but experimental
and *not recommended*.

.. note::
    The ViPErLEED graphical user interface, used, e.g., for acquiring |LEED-IV|
    :ref:`measurements<hardware>`, runs natively on Windows, and does not
    require WSL.

Installing WSL may require enabling some "developer" features (especially
on Windows 10). To this end, open the :guilabel:`Run` dialog (\ :kbd:`Win+R`),
type ``optionalfeatures``, then press :kbd:`Enter`. Make sure that
:guilabel:`Virtual Machine Platform` and
:guilabel:`Windows Subsystem for Linux` are selected in the
:guilabel:`Windows Features` dialog. You may need to restart your computer
after selecting them. Then, on the Microsoft Store, download your preferred
Linux distribution (e.g., Ubuntu). Once your distribution is downloaded,
you can install it by clicking on it. Follow the instructions in the terminal
to create a new user (note: user names should be lowercase only). You can then
start WSL by typing

.. code-block:: bat

   wsl

in your terminal. You can find more information about WSL on the official
`instructions by Microsoft <https://learn.microsoft.com/en-us/windows/wsl/install>`__.


.. _use_venv:

Using a Python virtual environment
==================================

A Python virtual environment is a confined workspace in which Python
dependencies can be installed without affecting other components of
the system. Several packages exist to create and manage virtual environments.
We list short instructions for a few of them below.

``venv``
--------

``venv`` is not the most popular virtual-environment package, but has the
advantage of being distributed as part of the Python standard library and
thus requires no external dependencies.

.. tip::
    It is a good idea to collect multiple virtual environments into one
    root directory of your choice. This makes it easier to switch between
    virtual environments later on.

To create a new virtual environment, navigate in your terminal to the
path where you would like the virtual-environment to be saved:

.. tab-set::

    .. tab-item:: Linux, macOS, Windows Subsystem for Linux
        :sync: unix

        .. code-block:: bash

            cd <path/to/folder/where/virtual/environment/should/be/saved>

    .. tab-item:: Windows
        :sync: win

        .. code-block:: bat

            cd "<path\to\folder\where\virtual\environment\should\be\saved>"

Then, create a fresh virtual environment with

.. tab-set::

    .. tab-item:: Linux, macOS, Windows Subsystem for Linux
        :sync: unix

        .. code-block:: bash

            pythonX.Y -m venv <virtual_env_name>

    .. tab-item:: Windows
        :sync: win

        .. code-block:: bat

            py -X.Y -m venv <virtual_env_name>

where ``X.Y`` is the Python version of your choice. Note that this is the
version of Python that will be available (via ``python``, or, on Windows,
also ``py``) in the virtual environment just created.

Virtual environments must be **activated** before they can be used:

.. tab-set::

    .. tab-item:: Linux, macOS, Windows Subsystem for Linux
        :sync: unix

        .. code-block:: bash

            source "path/to/<virtual_env_name>/bin/activate"

    .. tab-item:: Windows
        :sync: win

        .. code-block:: bat

            call "path\to\<virtual_env_name>\Scripts\activate.bat"

To disable a virtual environment, call ``deactivate`` instead.


.. [1]  If you are completely new to using a terminal, take a look at
        one of the many introductions available online. For example,
        `this one <https://www.techtarget.com/searchwindowsserver/definition/command-line-interface-CLI>`__.
