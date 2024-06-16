.. _t_experiment:

T_EXPERIMENT
============

T_EXPERIMENT is used to
:ref:`automatically generate an initial guess for atomic vibrational amplitudes<vibrocc_auto>`
for the file :ref:`VIBROCC`. It is the temperature of the experiment in Kelvin.
T_EXPERIMENT will **never** be used if the :ref:`VIBROCC` file exists and 
defines a vibrational amplitude for every site.

**Default:** No default.
Execution cannot proceed if T_EXPERIMENT is not defined *and* a
VIBROCC file is missing.

**Allowed values:** positive float, temperature in Kelvin

**Syntax:**

::

   T_EXPERIMENT = 293

.. note::

    The parameters :ref:`T_DEBYE`, T_EXPERIMENT and :ref:`VIBR_AMP_SCALE`
    will normally be used only once, to calculate an initial guess for
    vibrational amplitudes and generate a :ref:`VIBROCC` file. Afterwards, 
    all three parameters will automatically be commented out in the 
    :ref:`PARAMETERS` file; the vibration amplitudes will be defined in the 
    VIBROCC file instead. Even if the parameters were un-commented again, 
    they would never be used as long as a :ref:`VIBROCC` file is present.
