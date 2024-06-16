.. _t_debye:

T_DEBYE
=======

T_DEBYE is used to
:ref:`automatically generate an initial guess for atomic vibrational amplitudes<vibrocc_auto>`
for the file :ref:`VIBROCC`. It is the Debye temperature, which is used to
calculate the vibrational amplitude of each atom in the system. T_DEBYE will
**never** be used if the :ref:`VIBROCC` file exists and defines a vibrational
amplitude for every site.

**Default:** No default. Execution cannot proceed if T_DEBYE is not
defined and VIBROCC file is missing.

**Allowed values:** positive float, temperature in Kelvin

**Syntax:**

::

   T_DEBYE = 330

.. note::

    The parameters T_DEBYE,
    :ref:`T_EXPERIMENT` and :ref:`VIBR_AMP_SCALE` will normally be used only
    once, to calculate an initial guess for vibrational amplitudes and generate
    a :ref:`VIBROCC` file. Afterwards, all three parameters will automatically
    be commented out in the :ref:`PARAMETERS` file; the vibration amplitudes
    will be defined in the VIBROCC file instead. Even if the parameters were
    un-commented again, they would never be used as long as a :ref:`VIBROCC`
    file is present.
