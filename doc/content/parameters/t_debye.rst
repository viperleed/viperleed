.. _t_debye:

T_DEBYE
=======

T_DEBYE is used when vibrational amplitudes are automatically generated with an empirical formula (see below).
T_DEBYE will **never** be used if the :ref:`VIBROCC file<viboccin>` exists and defines a vibrational amplitude for every site.

**Default:** No default. Execution cannot proceed if T_DEBYE is not 
defined and VIBROCC file is missing.

**Allowed values:** positive float, temperature in Kelvin

**Syntax:**

::

   T_DEBYE = 330

If no :ref:`VIBROCC file<VIBOCCIN>`  is found, or the VIBROCC file contains no value for the vibrational amplitude of a given element, the parameters T_DEBYE and `T_EXPERIMENT<t_experiment>` 
will be used to automatically generate the default vibrational amplitude for that element using:

.. math::
    u = \sqrt{(u^2)}

.. math::
    u^2 = \sqrt{\frac{(1 + 16(T/\theta)^2) * (9 \hbar)}{(4 m k_B \theta)}}


where :math:`\theta` is the Debye temperature and :math:`T` is the temperature of the experiment.

This yields one vibrational amplitude per chemical element.
However, vibrational amplitudes of e.g. surface atoms should normally be scaled by some factor with respect to atoms in the bulk. This scaling factor is defined by :ref:`VIBR_AMP_SCALE<vibr_amp_scale>`.

Note that the parameters T_DEBYE,
:ref:`T_EXPERIMENT<t_experiment>` and
:ref:`VIBR_AMP_SCALE<VIBR_AMP_SCALE>`
will normally be used only once, to calculate an initial guess for
vibrational amplitudes and generate a :ref:`VIBROCC file <viboccin>`. 
Afterwards, all three
parameters will automatically be commented out in the :ref:`PARAMETERS` file;
the vibration amplitudes will be defined in the VIBROCC file instead.
Even if the parameters were un-commented again, they would never be used
as long as a :ref:`VIBROCC file <viboccin>` is present.