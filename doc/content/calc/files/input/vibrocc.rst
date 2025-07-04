.. include:: /substitutions.rst

.. _vibrocc:

=======
VIBROCC
=======

The |VIBROCC| file lists the starting guesses for vibration amplitudes
(in ångstrom) and site occupations. The minimum input is a vibration
amplitude for each element in the :ref:`POSCAR` file. If the :ref:`ELEMENT_MIX`
parameter is defined for an element in the :ref:`PARAMETERS` file, explicitly
assigning vibration amplitudes and occupations to all sub-elements is
recommended. See also :ref:`this page<occdelta>` for instructions on how
to vary the occupation of a site during structure optimization.

Additionally, the |VIBROCC| file can contain a block defining offsets in
vibration amplitudes, occupation, or position per element for specific
atoms.

A |VIBROCC| file containing only some starting guesses for vibration
amplitudes can be generated automatically using the :ref:`VIBR_AMP_SCALE`,
:ref:`T_EXPERIMENT` and :ref:`T_DEBYE` parameters in the
:ref:`PARAMETERS` file. See :ref:`below<vibrocc_auto>` for details.


Example
-------

**PARAMETERS file**

::

  ELEMENT_MIX M = Fe Ni
  SITE_DEF O = surf top(2)
  SITE_DEF M = surf top(2)

**VIBROCC file**

::

   = Vibration Amplitudes
   M_def = Fe 0.1, Ni 0.1
   M_surf = Fe 0.125, Ni 0.12          !some comment
   O_def = 0.19
   O_surf = 0.18

   = Occupations
   M* = Fe 0.8, Ni 0.2
   O_surf = 0.95

   = Search offsets
   POS 4 = Fe 0.0 0.0 0.01, Ni 0.0 0.0 -0.01     ! PFe_def
   OCC 4 = Fe 0.01, Ni -0.01     ! PFe_def

The two main blocks in the |VIBROCC| files are ``Vibration Amplitudes``
and ``Occupations``. Lines starting with ``=`` indicate the start of a
block.


Vibration amplitudes and occupations
------------------------------------

In each block, properties can be defined for each site type (left-hand side of
``=``).  The site types are labelled as ``El_sitename``, where ``El`` is an
element as found in the :ref:`POSCAR` file, and ``sitename`` is a site name
defined in the :ref:`PARAMETERS` file under :ref:`SITEDEF`. By default, an
asterisk (``*``) is interpreted as a wildcard character, so ``O*`` will access
both ``O_top`` and ``O_def``.

If required, the left-hand parameters can also be interpreted fully as regular
expressions (see also:
`python re syntax <https://docs.python.org/3.7/library/re.html>`__ and
`python re HOWTO <https://docs.python.org/3/howto/regex.html>`__). This feature
is turned off by default to avoid unintentional issues with e.g. full stops in
site names (not recommended!), but can be turned on by inserting a line
``= regex on`` at any point in the |VIBROCC| file, and disabled later by the
line ``= regex off``. Note that if regular expressions is on, the asterisk
``*`` will *not* be a wildcard character any more (the equivalent would be
``.*``)!

On the right-hand side of the ``=`` sign, you can either give only one value,
or give multiple values for different elements. Here, the elements are either
the ones found in the :ref:`POSCAR` file, or the ones defined in
:ref:`element_mix`. If element names in the |POSCAR| file and in ELEMENT_MIX
overlap, the assignment will nevertheless be made only for the chemical
element, see :ref:`element name collision<ElementNameCollision>`.
If only one value is given in the ``Vibration Amplitudes`` block, the
vibration amplitudes for all elements in this site will be set to this
value. If only one value is given in the ``Occupations`` block, this value
will be set for the main site element (e.g., O for the ``O_top`` site), or for
all main elements in a site affected by :ref:`element_mix`. The occupations for
all other elements will be set to zero for this site.

Total occupation in a site can be smaller than one, which will be interpreted
as the rest being vacancies. Defining an occupation greater than one will
produce a warning and may halt execution; if execution proceeds, the occupation
will be re-scaled to 1.

For simple systems, the ``Occupations`` block need not contain values for
elements with 100% site occupation, and can even be left out entirely. The
default value is 1.0 for the site's main element and 0.0 for all other
elements. If the site is affected by :ref:`element_mix`, the occupation will
be evenly split between the sub-elements defined in :ref:`element_mix`. A
simple example with 100% occupations and no :ref:`element_mix`  might therefore
look like this:

::

   = Vibration Amplitudes
   Fe_def = 0.10
   Fe_surf = 0.18
   O_def = 0.19
   O_surf = 0.18


.. _vibrocc_search_offsets:

Search offsets
--------------

Apart from starting values for vibration amplitudes and occupations, the
|VIBROCC| file can contain an additional block called ``Search offsets``.
This can be used to, *for a specific atom*, define position, vibration,
or occupation offsets from the site's values. This has two use cases:

-  If a parameter, e.g. the vibration amplitude, is varied independently for
   the different atoms sharing a site type, the search result will likely yield
   different values for these atoms. These values will be written to the
   :ref:`OUT/VIBROCC <vibrocc_out>` file to initialize a potential continuation
   job with the exact results from the previous search, instead of an average.
-  If there are multiple elements sharing a site via :ref:`element_mix`, the
   positions of the different chemical species may be different depending on
   the element. This cannot be mapped in the |POSCAR| file or the reference
   calculation of :term:`TensErLEED`, but can be mapped to the calculation
   via the ``Search offsets`` block by defining different values for different
   elements in the site.

**Example:**

::

   = Search offsets
   POS 4 = Fe 0.0 0.0 0.01, Ni 0.0 0.0 -0.01   ! for atom number 4, displace iron atoms by 0.01 A away from the bulk and Ni atoms 0.01 A towards the bulk.
   OCC 4 = Fe 0.01, Ni -0.01                   ! for atom number four, there is 1% more iron and 1% less nickel than defined for the site type

The syntax for this block differs somewhat from the vibration amplitudes and
occupations. On the left-hand side, each line is expected to contain:

-  A flag ``POS`` / ``VIB`` / ``OCC`` defining what type of parameter should
   be modified
-  An atom number (corresponding to the number in the |POSCAR| file)

On the right-hand side, the syntax is similar to the vibration amplitudes
and displacements blocks. For vibration amplitudes or occupations, one value
per element is expected, while for position offsets, three values per element
are expected. The three values for geometry are cartesian |x|, |y| and |z|
offsets, in ångströms, where positive |z| means away from the surface.


.. _vibrocc_out:

OUT/VIBROCC
-----------

After executing a search, the |VIBROCC| file found in the |OUT| folder contains
the vibration amplitudes and occupations of the best-fit structure found during
the (last) search (i.e., the one with the smallest |R factor|). If atoms in the
same site were allowed to vary independently, the vibrations and occupations
written for each site will be the average, and values for the single atoms
will be written as :ref:`vibrocc_search_offsets`.

When |VIBROCC| is :ref:`automatically generated<vibrocc_auto>` during
:ref:`initialization`, the resulting |VIBROCC| file is stored in the
|OUT| folder.

At the end of each |calc| execution, the file given as input for that run
is renamed to :file:`VIBROCC_ori`, while the (potentially) edited file
is copied to the root directory (from |OUT|) as a new |VIBROCC| file.
This ensures that further invocations of |calc| will automatically use
the output of previous executions as an input. You can manually call the
|bookkeeper| utility after a specific |calc| run if this behavior is not
desirable. See the :ref:`bookkeeper` page for more details.

.. note::
    A non-halted execution (i.e., one where :ref:`halting` was set to a
    value larger than the default) that includes a structure optimization
    will overwrite an auto-generated :file:`OUT/VIBROCC` file with the one
    found by the (last) optimization step. In this case, a copy of the
    auto-generated |VIBROCC| file can be found in |SUPP| (named
    :file:`VIBROCC_generated`).

.. versionchanged:: 0.13.0
    In earlier versions of |calc|, the automatically created |VIBROCC| file
    would only appear in the root directory after :ref:`initialization`,
    and only the |VIBROCC| file resulting from a structural optimization
    would be stored in |OUT|. This file used to be named :file:`VIBROCC_OUT`.


.. _vibrocc_auto:

Automatic generation of |VIBROCC|
---------------------------------

ViPErLEED can automatically generate a |VIBROCC| file containing starting
guesses for vibration amplitudes.
To do this, the experiment temperature :math:`T` (:ref:`T_EXPERIMENT`) and the
sample Debye temperature :math:`\Theta_\mathrm{D}` (:ref:`T_DEBYE`) must be
specified in :ref:`PARAMETERS`.
Additionally, :ref:`VIBR_AMP_SCALE` must be set if you are using non-default
sites (which is generally recommended).

Given these parameters and the atomic masses :math:`m`, the atomic vibration
amplitudes can be estimated as
:cite:p:`tongTheoryLowenergyElectron1975,vanhoveSurfaceCrystallographyLEED1979`

.. math::
    \langle u^2 \rangle _{T} \approx \frac{9 \hbar^2}{4 m k_\mathrm{B} \Theta_\mathrm{D}} \left[1+ 4\left(\frac{T}{\Theta_\mathrm{D}}\right)^2 \int_{0}^{\frac{\Theta_\mathrm{D}}{T}} \frac{x}{e^x - 1} dx\right].

Here :math:`\hbar` and :math:`k_\mathrm{B}` are the reduced Planck's and
Boltzmann's constants, respectively.

The integral can not be evaluated analytically, but a good approximation is
given by a combination of its low- and high-temperature limits

.. math::
    \langle u^2 \rangle _{T} \approx \sqrt{ (\langle u^2 \rangle _{T=0})^2 + (\langle u^2 \rangle _{T \rightarrow \infty})^2 }.

Evaluating these limits,

.. math ::
      \langle u^2 \rangle _{T=0} = \frac{9 \hbar^2}{4 m k_\mathrm{B} \Theta_\mathrm{D}},

      \langle u^2 \rangle _{T \rightarrow \infty} = \frac{9 \hbar^2 T}{m k_\mathrm{B} \Theta_\mathrm{D}^2},

gives

.. math ::
      \langle u^2 \rangle _{T} \approx \frac{9 \hbar^2}{4 m k_\mathrm{B} \Theta_\mathrm{D}} \sqrt{1+16\left(\frac{T}{\Theta_\mathrm{D}}\right)^2}.
