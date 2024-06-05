.. _muftin:

=======
V0_REAL
=======

.. note::
   Some versions of TensErLEED currently do not support arbitrary
   definitions of V0_REAL, instead always using the "Rundgren form"
   (see below) with parameters from the PHASESHIFTS file.

**TODO**: Update documentation to match TensErLEED 2.0

V0_REAL is used to provide the real part of the inner potential of the solid.
(:ref:`V0_IMAG<v0_imag>`  is the imaginary part, instead)

**Default**: ``V0_REAL = RUNDGREN c0 c1 c2 c3``, where the ``c``\ \* values are
taken from the first line in PHASESHIFTS, as derived from the output of the
phase shifts calculation.

Syntax
------

..  code-block:: none

   V0_REAL = -1*max(-10.17, -0.08 - 74.19/sqrt(EE+19.18))
   V0_REAL = RUNDGREN -10.17 -0.08 -74.19 19.18

**Acceptable values**: The right-hand side should be any real-valued function
of the electron energy (in electronvolts). Use only ``EE``, ``ee``, ``eE``, or
``Ee`` to represent the electron energy. The expression will be interpreted by
Fortran, so follow Fortran syntax. Acceptable arithmetic/mathematic functions
are listed below. The special command ``RUNDGREN`` can be used to choose the
following functional form for the real part of the inner potential

V(EE) = :ref:`FILAMENT_WF<FILWF>`  - max(c0,c1+c2/sqrt(EE+c3)),

as per Eq. (A8) in Rundgren's paper,
Ref. :cite:p:`rundgrenOptimizedSurfaceslabExcitedstate2003`.

.. seealso::
    :cite:t:`rundgrenElasticElectronatomScattering2007,rundgrenLowenergyElectronDiffraction2021`

The same result can be obtained by the input

..  code-block:: none

   V0_REAL = -1*max(c0,c1+c2/sqrt(EE+c3))

Notice that, in this case, it's necessary that c0<0 and c1<0.

It is advisable to **stick to the Default** (i.e., do not define
V0_REAL), unless you have provided an externally generated
:ref:`PHASESHIFTS<PHASESHIFTS>` file. In this case, it is best to
define the parameter with the ``RUNDGREN`` command and copying the
c0-c3 constants from the first line of any of the PS.r.\* output
files of the phase-shift calculation tool (c0 is the second number,
c1 the third, and so on).

In all cases, the program will replace ``EE`` with
``E``\ +:ref:`FILAMENT_WF<FILWF>`, since the relevant
electron energy is the one in vacuum, with respect to Fermi.

**Acceptable math expressions**: all names are case insensitive, all angles
are in RADIANS, use parentheses '()' to indicate precedence of operation, as
well as for surrounding function arguments.

===================== ======= ===================================
Operation             Symbol  Syntax example
===================== ======= ===================================
Exponentiation        \*\*    a**b (a to the power of b)
Multiplication        \*      a*b (a times b)
Division              /       a/b (a divided by b)
Sum                   \+      a+b (a plus b)
Subtraction, Negation \-      a-b (a minus b), -a (negative of a)
Absolute value        abs()   abs(a) = \|a\|
Arc-cosine            acos()  acos(a), result in radians
Arc-sine              asin()  asin(a), result in radians
Arc-tangent           atan()  atan(a), result in radians
Complex conjugate     conjg() conjg(a+bi)=a-bi
Cosine                cos()   cos(a), a in radians
Hyperbolic cosine     cosh()  cosh(a)
Error function        erf()   erf(a)
Exponential           exp()   exp(a) = e to the power of a
Imaginary part        imag()  imag(a+bi)=a (real)
Natural logarithm     log()   log(a)
Maximum               max()   max(a,b,c,..)
Minimum               min()   min(a,b,c,..)
Sine                  sin()   sin(a), a in radians
Hyperbolic sine       sinh()  sinh(a)
Square root           sqrt()  sqrt(a)
Tangent               tan()   tan(a), a in radians
Hyperbolic tangent    tanh()  tanh(a)
===================== ======= ===================================
