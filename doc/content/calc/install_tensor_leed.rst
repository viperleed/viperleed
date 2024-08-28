.. include:: /substitutions.rst

.. |oneAPI| replace:: Intel oneAPI®


.. _install_tensor_leed_deps:

===================================
Installing tensor-LEED dependencies
===================================

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
Note that you may need to decompress the TensErLEED source files after
downloading. See also :numref:`list_tensorleed_folder` for the expected
folder structure.

If you want to use older version of TensErLEED, you can also download older
releases from the
`releases tab <https://github.com/viperleed/viperleed-tensorleed/releases>`__
of the ``viperleed-tensorleed`` repository. You can have multiple versions of
TensErLEED on your system at the same time. |calc| will use the most recent one
by default. To select a different version you can use the :ref:`TL_VERSION`
parameter.

When using multiple tensor-LEED versions, the folder containing the tensor-LEED
code is expected to have the structure in :numref:`list_tensorleed_folder`. The
name of the top-level folder is up to the user (\ :file:`my_tensorleed` in
:numref:`list_tensorleed_folder`).
However, folder names for the TensErLEED source code should be named
:file:`TensErLEED-v1.X.Y` (or :file:`TensErLEED-v1.XY`) for versions earlier
than v2.0.0 (see parameter :ref:`tl_version` for more details). The folder name
for later versions is only required to begin with :file:`TensErLEED`. The
top-level folder (\ :file:`my_tensorleed` in :numref:`list_tensorleed_folder`)
is expected to contain the compiled version of the :ref:`eeasisss_compile`
source code, obtained as described in :ref:`static_compile`.

.. _list_tensorleed_folder:
.. code-block::
    :caption:
        Expected structure of the directory
        containing the tensor-LEED source code.

    my_tensorleed/
    ├── eeasisss_src/       <-- EEASiSSS source code, compiled with make
    │   └── ...
    ├── TensErLEED-v2.0.0/
    │   ├── src/
    │   │    ├── ref-calc.f
    │   │    └── ...
    │   └── lib/
    │        ├── lib.tleed.f
    │        └── ...
    ├── TensErLEED-v1.7.3/  <-- optional older release
    │   ├── random_.c
    │   ├── Makefile        <-- For random_.c, versions < v.1.7.4 only
    │   └── ...
    ⋮
    └── Makefile            <-- For EEASiSSS


|calc| will need to know where the tensor-LEED code is located on your machine.
You can either specify this each time you run |calc| via
:ref:`command-line arguments<cli_calc>` or you can define
the :envvar:`VIPERLEED_TENSORLEED` environment variable.
See :ref:`set_envvar` for more details. With reference to
:numref:`list_tensorleed_folder`, the :envvar:`VIPERLEED_TENSORLEED`
environment variable (or the command-line argument) can equivalently point
to either ``my_tensorleed`` or to any of its :file:`TensErLEED*` subfolders.


.. _set_envvar:

Setting an environment variable
===============================

In this section we briefly describe how an :term:`environment variable`
is set on various platforms. The following examples are for the
:envvar:`VIPERLEED_TENSORLEED` environment variable that |calc|
can use to locate the :ref:`tensor-LEED code<install_tensorleed>`.

.. tab-set::

  .. tab-item:: Linux, macOS, Windows Subsystem for Linux

        You can set an environment variable via

        .. code-block:: bash

            export VIPERLEED_TENSORLEED="<path/to/your/local/copy/of/viperleed-tensorleed>"

        The above will only set the :envvar:`VIPERLEED_TENSORLEED` environment
        variable for the current session of the shell. To make the changes
        permanent, you can add the ``export`` statement, e.g., to
        your :file:`~/.bashrc`, via

        .. code-block:: bash

            echo '# Set environment variable for viperleed-tensorleed' >> ~/.bashrc
            echo 'export VIPERLEED_TENSORLEED="<path/to/your/local/copy/of/viperleed-tensorleed>"' >> ~/.bashrc

  .. tab-item:: Windows, via CMD

        Execute

        .. code-block:: bat

            setx VIPERLEED_TENSORLEED "<path/to/your/local/copy/of/viperleed-tensorleed"

        in a :program:`CMD` terminal. You then may need to reboot your
        system. Notice the use of ``setx`` rather than ``set``: the latter
        only sets the :envvar:`VIPERLEED_TENSORLEED` environment variable
        for the active session.

  .. tab-item:: Windows, via System Properties

        Open up the Windows :guilabel:`Run` dialog (shortcut :kbd:`Win+R`),
        type ``sysdm.cpl``, then press :kbd:`Enter` to open the
        :guilabel:`System Properties` editor. Navigate to the
        :guilabel:`Advanced` tab, and click the
        :guilabel:`Environment Variables...` button.
        In the :guilabel:`User variables for <user name>` section,
        click :guilabel:`New...`. Use ``VIPERLEED_TENSORLEED`` as
        the :guilabel:`Variable name`. For :guilabel:`Variable value`
        use the full path to the folder containing the archives you
        have downloaded from the ``viperleed-tensorleed`` GitHub page.
        You may need to reboot your system.


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
You can set the :ref:`FORTRAN_COMP` parameter to use any other Fortran compiler
installed on your system. |calc| will default to using :term:`ifort` if
available. Use the parameter :ref:`FORTRAN_COMP` to adjust this behavior.

The :ref:`structure-search section<sec_search>`, which is usually the
computationally most expensive part of ViPErLEED and TensErLEED, supports
parallelized compilation and execution with :term:`MPI`. To use the :term:`MPI`
version of TensErLEED, you need to install an :term:`MPI` implementation
as well as the :term:`MPI` compiler corresponding to your Fortran compiler.
We recommend using Open MPI on Linux and macOS.
The MPI compiler for :term:`gfortran` is :term:`mpifort`, for :term:`ifort`
it is :term:`mpiifort`.

If you are running |calc| on an :term:`HPC` system, appropriate Fortran
compilers and a preferred :term:`MPI` implementation are likely already
installed. For details regarding their usage, consult the documentation
for your HPC system as well as its administrators.

|calc| can run on Linux, macOS and Microsoft Windows, but the installation
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
    build tools (e.g., :program:`cmake`) are present:

    .. code-block:: bash

        sudo apt install cmake pkg-config build-essential -y

    Then, we configure the |oneAPI| installation such that it is discovered
    by our environment. For this, we need to :program:`source` the correct
    :file:`*vars.sh` file that sets the required :term:`CLI` arguments.
    We recommend you do this by adding the following line to the end of
    the startup script of your shell (usually :file:`~/.bashrc`):

    .. code-block:: bash

        source <full/path/to/intel/oneapi/shell/script>

    Follow
    `this guide by Intel <https://www.intel.com/content/www/us/en/docs/oneapi/programming-guide/2024-1/use-the-setvars-and-oneapi-vars-scripts-with-linux.html>`__
    to determine the correct shell script for your release version.

    Afterwards, the required compilers should be available for use.
    You can check whether :term:`ifort` is present by using

    .. code-block:: bash

        which ifort

    If the result is a path, it means that the shell knows the compiler exists.
    You can do the same check with ``mpirun`` and ``mpiifort`` to check that
    they are properly configured as well.

  .. tab-item:: macOS

    .. warning::
        Newer Macs using "Apple Silicon" ARM-based chips are incompatible
        with the Intel compilers (since they don't use Intel chips).
        Use :ref:`install_gcc` instead.

    Follow the
    `Intel guide <https://www.intel.com/content/www/us/en/develop/documentation/installation-guide-for-intel-oneapi-toolkits-macos/top.html>`__
    to install the |oneAPI| toolkits under macOS. As for
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

    Notice that :program:`cmd` can, in principle, use a mechanism similar
    to the :file:`.bashrc` startup script for Linux's :program:`bash`.
    However, there are limitations, and the |oneAPI| batch scripts
    **cannot be run automatically** together with :program:`cmd`.

    The easiest alternative is to create an alias batch file, and add its
    location to the :envvar:`Path` environment variable to limit the amount
    of typing. This way the alias batch file will be available in every
    :program:`cmd` session. Here the steps for an example setting:

    -   Navigate to your user-profile directory in :program:`Explorer`. To find
        the location of your user directory type

        .. code-block:: bat

            echo %USERPROFILE%

        in a terminal. On Windows 10 and later, you can also directly navigate
        to the location in :program:`Explorer` by typing ``%USERPROFILE%`` in
        the Windows Start menu or in the address bar of an :program:`Explorer`
        window.
    -   Create a folder :file:`.cmd_aliases`
    -   There, create a file :file:`activate_intel.bat` with the following
        contents

        .. code-block:: bat

            @echo off
            call "<full/path/to/correct/intel/oneapi/bat/file>" 1>nul
            echo Done
    -   Add the path to :file:`.cmd_aliases` to your :envvar:`Path`
        environment variable. See :ref:`set_envvar` for how to access the
        environment-variable settings on Windows. Edit the :envvar:`Path`
        environment variable by appending :file:`%USERPROFILE%\\.cmd_aliases`.
        On Windows 7, use a semicolon as separator.
    -   Open a new :program:`cmd` session. Activate |oneAPI| by typing

        .. code-block:: bat

            activate_intel.bat

        and wait for the "``Done``" reply. Then check that the |oneAPI|
        compilers are now visible via

        .. code-block:: bat

            ifort -D
            mpiifort -D
    
    .. important::
        You will have to
        
        .. code-block:: bat

            activate_intel.bat
        
        **every time** you open a **new** :program:`cmd` session to have
        the |oneAPI| compilers available.


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
    newest version of ``gfortran``. In this manual, we use :program:`apt`, the
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


  .. tab-item:: macOS

    For running under macOS, it is recommended to first install a package
    manager such as `brew <https://brew.sh>`__. This will also install the
    :program:`XCode` command-line tools which are required for installing
    most other components.

    Using the :program:`brew` command, you can then easily install ``gfortran``
    and the Open MPI implementation (automatically including ``mpifort``).

    .. code-block:: bash

        brew install gfortran
        brew install open-mpi

    There is no need to install :term:`BLAS` and :term:`LAPACK`, as macOS
    already ships with these libraries preinstalled.

    .. warning::
        If the :program:`XCode` command-line tools are not installed before
        you install ``gfortran``, you will get an error stating that the
        ``-lSystem`` library is not available. If this happens, make sure
        to first install the :program:`XCode` command-line tools via

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
        The notes below are for installing the very basic, non-optimized
        LAPACK/BLAS versions. Hence, execution of the code will be rather
        **slow**.

    -  Install `MSys2 <https://www.msys2.org/>`__, which then installs MinGW,
       then open the MSys2 shell.
    -  Update :program:`MSys2` running

       .. code-block:: bash

            pacman -Syu

    -  Install ``gfortran`` and other useful packages via

       .. code-block:: bash

            pacman -S mingw-w64-x86_64-toolchain

    -  Add the :file:`<path_to_mingw_installation>/mingw64/bin` entry to your
       :envvar:`Path` environment variable. This way, calling ``gfortran`` from the
       terminal will find the one just installed with no need to explicitly
       passing the whole path.
       See :ref:`set_envvar` for how to access the environment-variable
       settings on Windows. Edit the :envvar:`Path` environment variable
       by appending :file:`<path_to_mingw_installation>/mingw64/bin`.
       On Windows 7, use a semicolon as separator.

    -  Install developer tools, :program:`cmake`, and :program:`git` with

       .. code-block:: bash

           pacman -S base-devel
           pacman -S mingw-w64-x86_64-cmake
           pacman -S git

    -  Clone the LAPACK git repository with

       .. code-block:: bash

            git clone https://github.com/msys2/MINGW-packages.git

       This is the basic, non-optimized version. There are ways to also
       build better versions.
       See `here <https://icl.cs.utk.edu/lapack-for-windows/lapack/>`__.
    -  Move to the LAPACK directory with

       .. code-block:: bash

            cd MINGW-packages/mingw-w64-lapack
    -  Build LAPACK and BLAS packages with

       .. code-block:: bash

            makepkg-mingw

       Should :program:`curl` complain about some certificates, you can also
       `download <http://www.netlib.org/lapack/>`__ the LAPACK/BLAS source
       code as a ``.tar.gz`` archive.
       Take the version that :program:`curl` complains about, and place the
       archive in the package folder, which you can find in
       :file:`<path_to_mingw_installation>/home/<user_name>/MINGW-packages/mingw-w64-lapack`.
       This build will take quite a while.

    -  Install LAPACK/BLAS packages with

       .. code-block:: bash

            pacman -U mingw-w64-x86_64-lapack-<REPLACE_WITH_VERSION>.pkg.tar.zst

       Note that the archive may have a different suffix. Run ``ls`` in
       the same folder to check the correct name.

    You can then test the LAPACK installation with

    .. code-block:: bash

       cd ~
       wget http://www.math.ucla.edu/~wotaoyin/software/lapack_test.cpp  # download
       g++ lapack_test.cpp -llapack -o lapack_test     # build
       ./lapack_test                                   # run

    For actually running, set the :ref:`FORTRAN_COMP` parameter in the
    :ref:`PARAMETERS` file as follows:

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

.. _eeasisss_compile:

EEASiSSS
--------

.. todo::
    Compile eeasisss at run time, instead of using a pre-compiled version?

This is the "Elastic Electron–Atom Scattering in Solids and Surface Slabs"
(\ :program:`EEASiSSS`) program by
:cite:t:`rundgrenElasticElectronatomScattering2007`.
Its source code is distributed by the ViPErLEED developers
in the ``viperleed-tensorleed`` GitHub
`repository <https://github.com/viperleed/viperleed-tensorleed>`__
with permission from the author.
See also :ref:`this section<install_tensorleed>` for more information
on how to obtain the source code.

:program:`EEASiSSS` is used by ViPErLEED during the :ref:`initialization`
section to generate the :ref:`PHASESHIFTS` file.

Install the Fortran compiler of your choice following
:ref:`these instructions<install_fortran_comp>`. Then proceed to
compilation from source as described in the following.

.. tab-set::

    .. tab-item:: Linux, macOS, Windows Subsystem for Linux

        :program:`EEASiSSS` can be compiled automatically using
        the provided ``Makefile``. Navigate to your local version
        of the ``viperleed-tensorleed`` repository using

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

    .. tab-item:: Linux, macOS, Windows Subsystem for Linux

        To build the |R-factor| extension module, navigate to your local
        copy of  the ``viperleed`` package in the terminal and call
        ``make`` in the :file:`extensions` directory.

    .. tab-item:: Windows, native

        There are no automatic means to build the |R-factor| extension module
        on Windows. Native-Windows users must manually build the extension
        using F2PY. See the UNIX ``Makefile``  for build flags. You can find
        the ``Makefile`` in the :file:`extensions` directory of your local
        copy of the ``viperleed`` package.


.. _mpirandom:

Random numbers library for TensErLEED < v1.7.4
----------------------------------------------

.. note::
    Users wishing to run natively on Windows can skip this step,
    as the random numbers are used only in the structure-optimization
    :ref:`section<sec_search>`, which is currently incompatible with
    native execution on Windows due to limitations on the Python code.

TensErLEED versions up to v1.7.3 need the :term:`C`-object files called
``random_.o`` and ``MPIrandom_.o``. These files must be precompiled with C and
C MPI compilers, respectively. A ``Makefile`` is also provided for them. If you
followed the :ref:`instructions<install_fortran_comp>` for obtaining the
Fortran compilers, you should already have the necessary C compilers installed
(from either GCC or Intel).

To compile the random-number generation library for a certain TensErLEED
version, navigate in your terminal to the respective
:ref:`directory containing the TensErLEED source files<install_tensorleed>`
and call either ``make intel`` or ``make gcc`` to compile
them using the Intel or GCC :term:`C` compilers, respectively.


.. [#]  For other distributions have a look, for example, at this tutorial
        `<https://fortran-lang.org/en/learn/os_setup/install_gfortran/>`__.
