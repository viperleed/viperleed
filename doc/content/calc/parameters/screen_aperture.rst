.. _screen_aperture:

SCREEN_APERTURE
===============

SCREEN_APERTURE defines the aperture (in degrees) of the acceptance cone of
the LEED screen. This, in turn, determines the highest-order beams that will
be visible in the LEED experiment.

**Default**: SCREEN_APERTURE = 110 (i.e., the one of the ErLEED optics)

**Syntax**:

::

   SCREEN_APERTURE = 95

**Acceptable values**: integer/floating point number between 0 and 180

.. note::
   A beam with indices :math:`(h, k)` and corresponding in-plane reciprocal
   vector :math:`\mathbf{g} = h \mathbf{a}^* + k \mathbf{b}^*` is visible at
   energy :math:`E` (in electronvolts) if
   :math:`\theta <` ``SCREEN_APERTURE`` :math:`/ 2`,
   where :math:`\sin{\theta} = \frac{\hbar g}{\sqrt{2 m_e e E}}`.
   Here :math:`e` and :math:`m_e` are the elementary charge and electron
   mass respectively.
