.. include:: /substitutions.rst

.. |oneAPI| replace:: Intel oneAPI®

.. _installation:

============
Installation
============

This section provides a short guide on how to install all required
components of ViPErLEED.


Installing Python
=================

With the exception of the :ref:`spot_tracker` and the related ImageJ plugins,
ViPErLEED requires a working Python installation. For more information on
installing the ViPErLEED ImageJ plugins see
:ref:`the corresponding section<spot_tracker>`. You can install Python using
a package manager or as binary from the
`Python website <https://www.python.org/downloads/>`__.

.. note::

    ViPErLEED requires Python 3.7 or newer, though we recommend using the
    latest available Python version.


Installing the ViPErLEED Python package
=======================================

ViPErLEED is available on the :term:`Python Package Index (PyPI)<PyPI>`
and can be installed using the Python package manager ``pip``. We recommend
installing ViPErLEED in a
`virtual environment <https://docs.python.org/3/library/venv.html>`__
for easier dependency management.

To install ViPErLEED from :term:`PyPI` using ``pip``, run

.. code-block:: console

    python -m pip install viperleed

in your terminal.\ [1]_

This will install the latest version of ViPErLEED and all required
dependencies. It will also automatically install the
:ref:`ViPErLEED command-line tools<command_line_tools>` which can
be called from the terminal using the ``viperleed`` command.

If you want to also install the dependencies for running the ViPErLEED
graphical user interface, which also allows you to :ref:`measure<hardware>`
|LEED-IV| data, run

.. code-block:: console

    python -m pip install viperleed[GUI]

in your terminal.


.. _install_tensorleed:

Tensor-LEED source code
=======================
The |calc| package acts as a wrapper and feature extension to the
:term:`TensErLEED` package. It requires TensErLEED source files to be present
on your machine. You can obtain all the necessary source code from the
`viperleed-tensorleed <https://github.com/viperleed/viperleed-tensorleed>`__
GitHub repository. We suggest downloading the most recent
`release <https://github.com/viperleed/viperleed-tensorleed/releases/latest>`__
into a local directory of your choice.
Note that you may need to uncompress the TensErLEED source files after
downloading.

If you want to use older version of TensErLEED, you can also download older
releases from the
`releases tab <https://github.com/viperleed/viperleed-tensorleed/releases/latest>`__
of the viperleed-tensorleed repository. You can have multiple versions of
TensErLEED on your system at the same time. |calc| will use the most recent
available version by default. To select a different version you can use the
:ref:`TL_VERSION<tl_version>` parameter.

|calc| will need to know where the tensor-LEED code is located on your machine.
You can either specify this each time you run |calc| via
:ref:`command-line arguments<cli_calc>` or you can define
the ``VIPERLEED_TENSORLEED`` environment variable:

.. tab-set::

  .. tab-item:: Linux, MacOS, Windows Subsystem for Linux

        You can set an environment variable via

        .. code-block:: bash

            export VIPERLEED_TENSORLEED="<path/to/your/local/copy/of/viperleed-tensorleed>"

        The above will only set the ``VIPERLEED_TENSORLEED`` environment
        variable for the current session of the shell. To make the changes
        permanent, you can add the ``export`` statement, e.g., to
        your ``~/.bashrc``, via

        .. code-block:: bash

            echo '# Set environment variable for viperleed-tensorleed' > ~/.bashrc
            echo 'export VIPERLEED_TENSORLEED="<path/to/your/local/copy/of/viperleed-tensorleed>"' > ~/.bashrc

  .. tab-item:: Windows, via CMD

        Execute

        .. code-block:: bat

            setx VIPERLEED_TENSORLEED "<path/to/your/local/copy/of/viperleed-tensorleed"

        in a CMD terminal. You then may need to reboot your system. Notice
        the use of ``setx`` rather than ``set``: the latter only sets the
        ``VIPERLEED_TENSORLEED`` environment variable for the active session.

  .. tab-item:: Windows, via System Properties

        Open up the Windows "Run" dialog (shortcut :kbd:`Win+R`), type
        ``sysdm.cpl``, then press :kbd:`Enter` to open the "System Properties"
        editor. Navigate to the "Advanced" tab, and click the
        "Environment Variables..." button. In the "User variables for
        <user name>" section, click "New...". Use ``VIPERLEED_TENSORLEED`` as
        the "Variable name". For "Variable value" use the full path to the
        folder containing the archives you have downloaded from the
        ``viperleed-tensorleed`` GitHub page. You may need to reboot
        your system.

.. _wsl:

Windows Subsystem for Linux
===========================

To run |calc| under Windows, we recommend using the
:term:`Windows Subsystem for Linux<WSL>` (WSL), available starting from
Windows 10. Running |calc| natively on Windows is possible, but experimental
and *not recommended*.

.. note::
    The ViPErLEED graphical user interface, e.g., used for acquiring |LEED-IV|
    :ref:`measurements<hardware>`, runs natively on Windows, and does not
    require WSL.

Installing WSL may require enabling some "developer" features (especially
on Windows 10). To this end, open the "Run" dialog (\ :kbd:`Win+R`), type
``optionalfeatures``, then press :kbd:`Enter`. Make sure that "Virtual
Machine Platform" and "Windows Subsystem for Linux" are selected in the
"Windows Features" dialog. You may need to restart your computer after
selecting them. Then, on the Microsoft Store, download your preferred
Linux distribution (e.g., Ubuntu). Once your distribution is downloaded,
you can install it by clicking on it. Follow the instructions in the terminal
to create a new user (note: user names should be lowercase only). You can then
start WSL by typing

.. code-block:: bat

   wsl

in your terminal. You can find more information about WSL on the official
`instructions by Microsoft <https://learn.microsoft.com/en-us/windows/wsl/install>`__.


.. _install_fortran_comp:

Fortran compilers
=================

|calc| needs to compile the TensErLEED source code *at run-time*. This requires
the presence of a suitable :term:`Fortran` 77 & 90 compiler on the system.
Unlike the original version of TensErLEED by
:cite:t:`blumFastLEEDIntensity2001a`, all TensErLEED versions supported by
|calc| (TensErLEED ≥ 1.6) also require :term:`BLAS` and :term:`LAPACK`
libraries to be available. |calc| supports :term:`gfortran` from the GNU
Compiler Collection (:term:`gcc`) and the Intel Fortran compiler :term:`ifort`
without additional configuration.
You can set the :ref:`FORTRAN_COMP<fortran_comp>` parameter to use any other
Fortran compiler installed on your system. |calc| will default to using
:term:`ifort` if available. Use the parameter :ref:`FORTRAN_COMP<fortran_comp>`
to adjust this behavior.

The :ref:`structure-search section<sec_search>`, which is usually the
computationally most expensive part of ViPErLEED and TensErLEED, supports
parallelized compilation and execution with :term:`MPI`. To use the :term:`MPI`
version of TensErLEED, you need to install an :term:`MPI` implementation
as well as the :term:`MPI` compiler corresponding to your Fortran compiler.
We recommend using Open MPI on Linux and MacOS.
The MPI compiler for :term:`gfortran` is :term:`mpifort`, for :term:`ifort`
it is :term:`mpiifort`.

If you are running |calc| on an :term:`HPC` system, appropriate Fortran
compilers and a preferred :term:`MPI` implementation are likely already
installed. For details regarding their usage, consult the documentation
for your HPC system as well as its administrators.

|calc| can run on Linux, MacOS and Microsoft Windows, but the installation
of the compilers differs significantly for each system.

.. note::

    -  If you are running on an system based on an Intel processor, we
       recommend using ``ifort``. It is known from experience to give
       better performance for TensErLEED.
    -  Using the :term:`MPI` version of TensErLEED is not strictly required,
       but **highly** recommended.
       Execution of the :ref:`structure search<sec_search>` may take
       significantly longer without :term:`MPI`. A working MPI
       implementation is necessary to make use of multi-processing
       in the :ref:`structure-search section<sec_search>`, even if
       you are working on a single node.



``ifort`` and ``mpiifort``
--------------------------

This section provides a guide on how to install the Intel Fortran compiler
:term:`ifort`, an MPI implementation, and the ``ifort`` MPI compiler
:term:`mpiifort`. All the necessary components are packaged as part of the
|oneAPI| toolkits.

ViPErLEED requires the |oneAPI| Base Toolkit and the |oneAPI| HPC Toolkit.

.. note::
    The toolkits are multiple gigabytes in size and
    will take a while to download and install.

The :term:`BLAS` and :term:`LAPACK` libraries are packaged in the |oneAPI|
Math Kernel Library (MKL), which is part of the |oneAPI| Base Toolkit.
The |oneAPI| HPC Toolkit contains the Fortran compilers as well as
an :term:`MPI` implementation.

The full documentation of the |oneAPI| is available from the
`Intel website <https://www.intel.com/content/www/us/en/develop/documentation/get-started-with-intel-oneapi-hpc-linux/top/before-you-begin.html#before-you-begin>`__.

.. _ifort_linux:

.. tab-set::

  .. tab-item:: Linux

    Installation of the Intel compilers and of the :term:`MPI` implementation
    for Linux can be performed using a few shell commands. In this manual, we
    use ``apt``, the standard package manager for Debian-based distributions.
    For installation instructions with other package managers see the
    `guides by Intel <https://www.intel.com/content/www/us/en/develop/documentation/installation-guide-for-intel-oneapi-toolkits-linux/top.html>`__.

    As a first step, update the package index:

    .. code-block:: bash

        sudo apt update && sudo apt upgrade

    After adding the |oneAPI| repository following the
    `instructions by Intel <https://www.intel.com/content/www/us/en/develop/documentation/installation-guide-for-intel-oneapi-toolkits-linux/top/installation/install-using-package-managers/apt.html#apt>`__,
    install the required packages:

    .. code-block:: bash

        sudo apt install intel-basekit -y
        sudo apt install intel-hpckit -y

    Once installation completes, we need to configure the system and add the
    compilers to our system path. First, we need to make sure the required
    build tools (e.g., ``cmake``) are present:

    .. code-block:: bash

        sudo apt install cmake pkg-config build-essential -y

    Then, we configure the |oneAPI| installation such that it is discovered
    by our environment. For this, we need to ``source`` the file
    ``/opt/intel/oneapi/setvars.sh`` that sets the required :term:`CLI`
    arguments. We recommend you do this by adding the following line to
    the end of the startup script of your shell (usually ``~/.bashrc``):

    .. code-block:: bash

        source <install-dir>/<toolkit-version>/oneapi-vars.sh

    Afterwards, the required compilers should be available for use.
    See also
    `this page by Intel <https://www.intel.com/content/www/us/en/docs/oneapi/programming-guide/2024-1/use-the-setvars-and-oneapi-vars-scripts-with-linux.html>`__
    for more details.

    You can check whether :term:`ifort` is present by using

    .. code-block:: bash

        which ifort

    If the result is a path, it means that the shell knows the compiler exists.
    You can do the same check with ``mpirun`` and ``mpiifort`` to check that
    they are properly configured as well.

  .. tab-item:: MacOS

    .. warning::
        Newer Macs using "Apple Silicon" ARM-based chips are incompatible
        with the Intel compilers (since they don't use Intel chips).
        Use :ref:`install_gcc` instead.

    Follow the
    `Intel guide <https://www.intel.com/content/www/us/en/develop/documentation/installation-guide-for-intel-oneapi-toolkits-macos/top.html>`__
    to install the |oneAPI| toolkits under MacOS. As for
    :ref:`Linux<ifort_linux>`, you will need to install the |oneAPI| Base
    Toolkit and the |oneAPI| HPC Toolkit.

  .. tab-item:: Windows Subsystem for Linux

    For information on installing WSL see :ref:`this section<wsl>`. Start WSL
    by typing

    .. code-block::

        wsl

    in your terminal, then follow the same instructions as in
    :ref:`the Linux section<ifort_linux>`.

  .. tab-item:: Windows, native

    .. warning::
        To run |calc| and TensErLEED under Windows, we recommend using the
        :term:`Windows Subsystem for Linux<WSL>` (WSL, available starting
        from Windows 10). Running natively on Windows is possible, but
        experimental and *not recommended*.

    .. warning::
        The :ref:`structure-optimization section<sec_search>` currently
        contains Python code that is incompatible with Windows. Therefore,
        a full |LEED-IV| calculation cannot be performed under Windows.

    Follow the Intel
    `guide <https://www.intel.com/content/www/us/en/develop/documentation/installation-guide-for-intel-oneapi-toolkits-windows/top.html>`__
    to install the |oneAPI| toolkits under Windows. As for Linux, you will
    need to install the |oneAPI| Base Toolkit and the |oneAPI| HPC Toolkit.

    The |oneAPI| toolkits require specific environment variables to be set
    before compilers can be used. The toolkits come with dedicated ``.bat``
    scripts that must be executed on each session of the terminal. For more
    information concerning which script to use, see
    `this guide <https://www.intel.com/content/www/us/en/docs/oneapi/programming-guide/2024-1/use-the-setvars-script-with-windows.html>`__
    from Intel.

    Notice that ``cmd.exe`` can use a mechanism similar to the ``.barchrc``
    startup script for Linux's ``bash``. This means that dedicated commands
    can be executed upon startup of each ``cmd`` session. To set this up:

    -  Open the Windows registry editor via the "Run" dialog (\ :kbd:`Win+R`):
       type ``regedit``, then press :kbd:`Enter`. Confirm the administrator
       permissions.
    -  Navigate to ``HKEY_CURRENT_USER\Software\Microsoft\Command Processor``.
       If ``Command Processor`` does not exist, right click on ``Microsoft``,
       select "New" → "Key", and name the new entry ``Command Processor``.
    -  Find the ``AutoRun`` entry in the right panel. If no ``AutoRun`` exists,
       create a "New" → "String value" by right clicking. Name it ``AutoRun``.
    -  Skip this passage if the ``AutoRun`` entry already has a value.
       "Modify..." the empty value of the ``AutoRun`` entry to
       ``%USERPROFILE%\cmd-autorun.bat``. ``%USERPROFILE%`` is the path
       to your user directory. You can find its value by typing

       .. code-block:: bat

          echo %USERPROFILE%

       in a terminal. On Windows 10 and later, you can also directly navigate
       to the location in Explorer by typing ``%USERPROFILE%`` in the Windows
       Start menu or in the address bar of an Explorer window.
    -  Navigate to the file path set as a value for the ``AutoRun`` entry.
       You can create a new file with the correct name if it does not exist.
       Then append the following lines to the end:

       .. code-block:: bat

         @echo off
         call "<full/path/to/correct/intel/oneapi/bat/file>"
         @echo on

    -  Close the registry editor.
    -  Open a new ``cmd.exe`` session, and test that the |oneAPI| compilers
       are now visible via

       .. code-block:: bat

         ifort --version
         mpiifort --version


.. _install_gcc:

``gfortran`` and ``mpifort``
----------------------------

This section provides a simple guide on how to install the GNU Fortran compiler
:term:`gfortran`, the Open MPI implementation, and the ``gfortran`` MPI wrapper
:term:`mpifort`.  See also the
`guide on the Fortran-language reference page <https://fortran-lang.org/en/learn/os_setup/install_gfortran/>`__
for how to install ``gfortran`` on various operating systems.


.. _gnu_linux:

.. tab-set::

  .. tab-item:: Linux

    First, using your package manager, update the package list and install the
    newest version of ``gfortran``. In this manual, we use ``apt``, the
    standard package manager for Debian-based distributions.\ [#]_

    .. code-block:: bash

        sudo apt update
        sudo apt install gfortran -y

    The compiler can be invoked with the ``gfortran`` command. You can show the
    version and check whether ``gfortran`` was installed properly using

    .. code-block:: bash

        gfortran --version

    In addition to ``gfortran``, you also need to install the :term:`BLAS`
    and :term:`LAPACK` libraries. Use

    .. code-block:: bash

        sudo apt install libblas-dev liblapack-dev

    Next, install Open MPI (or another MPI implementation of your
    choice) to make ``mpirun`` available:

    .. code-block:: bash

        sudo apt install openmpi-bin

    Finally, install the  ``gfortran`` MPI wrapper ``mpifort``:

    .. code-block:: bash

        sudo apt install libopenmpi-dev


  .. tab-item:: MacOS

    For running under MacOS, it is recommended to first install a package
    manager such as `brew <https://brew.sh>`__. This will also install the
    XCode command-line tools which are required for installing most other
    components.

    Using the ``brew`` command, you can then easily install ``gfortran`` and
    the Open MPI implementation (automatically including ``mpifort``).

    .. code-block:: bash

        brew install gfortran
        brew install open-mpi

    There is no need to install :term:`BLAS` and :term:`LAPACK`, as MacOS
    already ships with these libraries preinstalled.

    .. warning::
        If the XCode command-line tools are not installed before you install
        ``gfortran``, you will get an error stating that the ``-lSystem``
        library is not available. If this happens, make sure to first install
        the XCode command-line tools via

        .. code-block:: bash

            xcode-select --install

        and then reinstall ``gfortran`` with

        .. code-block:: bash

            brew reinstall gfortran


  .. tab-item:: Windows Subsystem for Linux

    Install WSL as described in :ref:`this section <wsl>`. Start WSL by typing

    .. code-block::

        wsl

    in your terminal, then follow the same instructions as for
    :ref:`Linux<gnu_linux>`.


  .. tab-item:: Windows (64 bit), native

    .. warning::
        To run |calc| and TensErLEED under Windows, we recommend using the
        :term:`Windows Subsystem for Linux<WSL>` (WSL, available starting
        from Windows 10). Running natively on Windows is possible, but
        experimental and *not recommended*.

    .. warning::
        The :ref:`structure-optimization section<sec_search>` currently
        contains Python code that is incompatible with Windows. Therefore,
        a full |LEED-IV| calculation cannot be performed under Windows.
        For this reason, this section does not describe a way to obtain
        an MPI implementation nor an the ``mpirun`` Fortran compiler.

    .. warning::
        The notes below are for installing the very basic, unoptimized
        LAPACK/BLAS versions. Hence, execution of the code will be
        rather **slow**.

    -  Install `MSys2 <https://www.msys2.org/>`__, which then installs MinGW,
       then open the MSys2 shell.
    -  Update MSys2 running

       .. code-block:: bash

            pacman -Syu

    -  Install ``gfortran`` and other useful packages via

       .. code-block:: bash

            pacman -S mingw-w64-x86_64-toolchain

    -  Add the ``"<path_to_mingw_installation>/mingw64/bin"`` entry to your
       ``Path`` environment variable. This way, calling ``gfortran`` from the
       terminal will find the one just installed with no need to explicitly
       passing the whole path.
       See the notes in :ref:`this section<install_tensorleed>` for how to
       access the environment-variable settings on Windows.
       Edit the ``Path`` environment variable by appending
       ``"<path_to_mingw_installation>/mingw64/bin"``. On Windows 7, use a
       semicolon as separator.

    -  Install developer tools, ``cmake`` and ``git`` with

       .. code-block:: bash

           pacman -S base-devel
           pacman -S mingw-w64-x86_64-cmake
           pacman -S git

    -  Clone the LAPACK git repository with

       .. code-block:: bash

            git clone https://github.com/msys2/MINGW-packages.git

       This is the basic, unoptimized version.
       There are ways to also build better versions . See
       `here <https://icl.cs.utk.edu/lapack-for-windows/lapack/>`__.
    -  Move to the LAPACK directory with

       .. code-block:: bash

            cd MINGW-packages/mingw-w64-lapack
    -  Build LAPACK and BLAS packages with

       .. code-block:: bash

            makepkg-mingw

       Should ``curl`` complain about some certificates, you can also
       `download <http://www.netlib.org/lapack/>`__ the LAPACK/BLAS source
       code as a ``.tar.gz`` archive.
       Take the version that ``curl`` complains about, and place the archive
       in the package folder, which you can find in
       ``<path_to_mingw_installation>/home/<user_name>/MINGW-packages/mingw-w64-lapack``.
       This build will take quite a while.

    -  Install LAPACK/BLAS packages with

       .. code-block:: bash

            pacman -U mingw-w64-x86_64-lapack-<REPLACE_WITH_VERSION>.pkg.tar.zst

       Note that the archive may have a different suffix. Run ``ls`` in the
       same folder to check the correct name.

    You can then test the LAPACK installation with

    .. code-block:: bash

       cd ~
       wget http://www.math.ucla.edu/~wotaoyin/software/lapack_test.cpp  # download
       g++ lapack_test.cpp -llapack -o lapack_test     # build
       ./lapack_test                                   # run

    For actually running, set the :ref:`FORTRAN_COMP<fortran_comp>` parameter
    in the :ref:`PARAMETERS file<parameters>` as follows:

    .. todo:: Michele: is -std=legacy required on native Windows?

    ::

       # -std=legacy makes it work for Fortran77
       FORTRAN_COMP = 'gfortran -O2 -std=legacy'
       # NOTE: order of LAPACK and BLAS is important!
       FORTRAN_COMP post = '-llapack -lblas -lpthread'


  .. tab-item:: Windows (32 bit), native

    There is currently no solution for installing ``gcc`` compilers with
    LAPACK/BLAS on 32-bit Windows.



.. _static_compile:

Compiling static files
======================

In addition to the TensErLEED source code, which is compiled *at run-time*,
ViPErLEED needs a few auxiliary programs that need compiling before a
calculation can be started.

EEASiSSS
--------

.. todo::
    Compile eeasisss at run time, instead of using a pre-compiled version?

This is the "Elastic Electron–Atom Scattering in Solids and Surface Slabs"
(EEASiSSS) program by :cite:t:`rundgrenElasticElectronatomScattering2007`.
Its source code is distributed by the ViPErLEED developers in the
``viperleed-tensorleed`` GitHub
`repository <https://github.com/viperleed/viperleed-tensorleed>`__
with permission from the author.

EEASiSSS is used by ViPErLEED during the :ref:`initialization<initialization>`
to generate the :ref:`PHASESHIFTS<phaseshifts>` file.

The EEASiSSS source code is included, together with TensErLEED, in the
``viperleed-tensorleed`` GitHub
`repository <https://github.com/viperleed/viperleed-tensorleed>`__.
See :ref:`this section<install_tensorleed>` for more information.

Install the Fortran compiler of your choice following
:ref:`these instructions<install_fortran_comp>`. Then proceed to
compilation from source as described in the following.

.. tab-set::

    .. tab-item:: Linux, MacOS, Windows Subsystem for Linux

        EEASiSSS can be compiled automatically using the provided ``Makefile``.

        Navigate to your local version of the ``viperleed-tensorleed``
        repository using

        .. code-block:: bash

            cd path/to/viperleed-tensorleed

        From there, call either ``make intel`` or ``make gcc`` to compile
        using the Intel or GCC Fortran compilers, respectively.

    .. tab-item:: Native Windows

        EEASiSSS can be compiled automatically using the provided ``make.bat``.

        Navigate to your local version of the ``viperleed-tensorleed``
        repository using

        .. code-block:: bat

            cd "path\to\viperleed-tensorleed"

        From there, call either ``make.bat intel`` or ``make.bat gcc`` to
        compile using the Intel or GCC Fortran compilers, respectively.

        An alternative to compilation from source is using the precompiled
        ``eeasisss_windows.exe`` executable available in the
        ``viperleed-tensorleed``
        `repository <https://github.com/viperleed/viperleed-tensorleed>`__.
        It should be renamed to ``eeasisss.exe`` after download.

        .. warning::
            The precompiled ``eeasisss_windows.exe`` executable was built on
            a machine that is likely different from yours. It may not work
            at all on your machine, or may produce unexpected results. It
            is safer to compile from source.


.. _rfactor_exentsion:

|R-factor| extension for ASE
----------------------------

For :ref:`using the atomic simulation environment<aseapi>` (ASE) with
ViPErLEED, you may need to compile the |R-factor| extension for |calc|.
This extension is used to calculate the |R factor| in conjunction with the ASE
package. It relies on `F2PY <https://numpy.org/doc/stable/f2py/index.html>`__,
which is installed by default with NumPy.

.. tab-set::

    .. tab-item:: Linux, MacOS, Windows Subsystem for Linux

        To build the |R-factor| extension module, navigate to your local
        copy of  the ``viperleed`` package in the terminal and call
        ``make`` in the ``extensions`` directory.

    .. tab-item:: Windows, native

        There are no automatic means to build the |R-factor| extension module
        on Windows. Native-Windows users must manually build the extension
        using F2PY. See the UNIX ``Makefile``  for build flags. You can find
        the ``Makefile`` in the ``extensions`` directory of your local copy of
        the ``viperleed`` package.


.. _mpirandom:

Randomizer library for TensErLEED < v1.7.4
------------------------------------------

.. note::
    Users wishing to run natively on Windows can skip this step,
    as the randomizers are used only in the structure-optimization
    :ref:`section<sec_search>`, which is currently incompatible with
    native execution on Windows due to limitations on the Python code.

TensErLEED versions up to v1.7.3 need the :term:`C`-object files called
``random_.o`` and ``MPIrandom_.o``. These files must be precompiled with C and
C MPI compilers, respectively. A ``Makefile`` is also provided for them. If you
followed the :ref:`instructions<install_fortran_comp>` for obtaining the
Fortran compilers, you should already have the necessary C compilers installed
(from either GCC or Intel).

To compile the randomizer library for TensErLEED version ``x.y.z``, go into
the respective `directory containing the TensErLEED source files and call either
``make intel`` or ``make gcc`` to compile them using the Intel or GCC :term:`C`
compilers, respectively.



.. [1]  If you are completely new to using a terminal, take a look at
        one of the many introductions available online. For example,
        `this one <https://www.techtarget.com/searchwindowsserver/definition/command-line-interface-CLI>`__.

.. [#]  For other distributions have a look, for example, at this tutorial
        `<https://fortran-lang.org/en/learn/os_setup/install_gfortran/>`__.
