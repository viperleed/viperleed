.. _t_experiment:

T_EXPERIMENT
============

T_EXPERIMENT is used when vibrational amplitudes are automatically generated with an empirical formula (see below). T_EXPERIMENT will **never** be used if the :ref:`VIBROCC file<VIBOCCIN>`  exists and defines a vibrational amplitude for every site.

**Default:** No default. 
Execution cannot proceed if T_EXPERIMENT is not defined *and* a 
VIBROCC file is missing.

**Allowed values:** positive float, temperature in Kelvin

**Syntax:**

::

   T_EXPERIMENT = 293

If no :ref:`VIBROCC file<VIBOCCIN>`  is found, or the VIBROCC file 
contains no value for the vibrational amplitude of a given element, the 
parameters T_EXPERIMENT and 
:ref:`T_DEBYE` will be used 
to automatically generate the default vibrational amplitude for that 
element using:

.. math::
    u = \sqrt{(u^2)}

.. math::
    u^2 = \sqrt{\frac{(1 + 16(T/\theta)^2) * (9 \hbar)}{(4 m k_B \theta)}}

where :math:`\theta` is the Debye temperature and :math:`T` is the 
temperature of the experiment.

This yields one vibrational amplitude per chemical element.
However, vibrational amplitudes of e.g. surface atoms should normally be
scaled by some factor with respect to atoms in the bulk. This scaling 
factor is defined by 
:ref:`VIBR_AMP_SCALE<VIBR_AMP_SCALE>`.

Note that the parameters T_EXPERIMENT, :ref:`T_DEBYE<T_DEBYE>` and 
:ref:`VIBR_AMP_SCALE <VIBR_AMP_SCALE>` will normally be used only once, 
to calculate an initial guess for vibrational amplitudes and generate a 
:ref:`VIBROCC file<VIBOCCIN>`. Afterwards, all three parameters will 
automatically be commented out in the :ref:`PARAMETERS file<parameters>`
; the vibration 
amplitudes will be defined in the :ref:`VIBROCC file<VIBOCCIN>` instead.
Even if the 
parameters were un-commented again, they would never be used as long as 
a :ref:`VIBROCC file<VIBOCCIN>` is present.
