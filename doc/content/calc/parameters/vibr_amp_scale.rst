.. _vibr_amp_scale:

VIBR_AMP_SCALE
==============

VIBR_AMP_SCALE is used when
:ref:`vibrational amplitudes are automatically generated<vibrocc_auto>`
from the :ref:`T_EXPERIMENT` and :ref:`T_DEBYE` parameters.
E.g., vibrational amplitudes of the surface atoms are normally larger than in
the bulk, so one wants to scale the bulk vibration amplitudes.
VIBR_AMP_SCALE will **never** be used if the VIBROCC file exists and defines
a vibrational amplitude for the site in question.

**Default:** 1.0 for every site

**Allowed values:** one positive float value per site

**Syntax:**

::

   VIBR_AMP_SCALE = Fe_surf 1.3, O_surf 1.3
   ! OR
   VIBR_AMP_SCALE = *surf 1.3

Scaling factors can be defined on a single line (comma delimited pairs),
or on multiple lines.

The site types are labelled as ``El_sitename``, where ``El`` is an element
as found in the :ref:`POSCAR` file, and ``sitename`` is a site name
defined in the :ref:`PARAMETERS` file under :ref:`SITEDEF`.
Asterisks ``*`` are treated as wildcard characters, so ``*surf`` in the example
above will match both ``Fe_surf`` and ``O_surf``.
(The same convention as in the VIBROCC file.)

In the example above, if the vibrational amplitude for Fe has
been calculated to be 0.1 Å, the vibrational amplitude for the
Fe_surf atoms will be set to 0.13 Å.

.. note::

    The parameters :ref:`T_DEBYE`, :ref:`T_EXPERIMENT` and VIBR_AMP_SCALE
    will normally be used only once, to calculate an initial guess for
    vibrational amplitudes and generate a :ref:`VIBROCC` file. Afterwards, 
    all three parameters will automatically be commented out in the
    :ref:`PARAMETERS` file; the vibration amplitudes will be defined in the
    VIBROCC file instead. Even if the parameters were un-commented again, 
    they would never be used as long as a :ref:`VIBROCC` file is present.
