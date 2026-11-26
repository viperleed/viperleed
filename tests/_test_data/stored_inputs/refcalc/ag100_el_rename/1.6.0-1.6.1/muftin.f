
C  Subroutine muftin contains explicit energy dependence of inner
C  potentials. The functional form should be left to the user entirely,
C  thus the subroutine is included in the input explicitly.

C  All quantities should be set in eV. Conversion to Hartree (atomic units)
C  is handled by ref-calc.f itself. Both the real and imaginary part should
C  be given as positive values (the program then handles the sign correctly).

      subroutine muftin(EEV,VO,VV,VPI,VPIS,VPIO)

C  global variable

C  EEV :  electron energy in the vacuum region
C  VO  :  difference between real part of inner potential in the bulk and in
C         topmost layer. Usually irrelevant -> set to zero.
C  VV  :  real part of the inner potential in the bulk.
C  VPI :  imaginary part of the inner potential.
C  VPIS:  in case different values for the imaginary part of the inner
C  VPIO:  potential for bulk and topmost layer were desired, VPIS would
C         correspond to the bulk value, VPIO would be the respective
C         value for the topmost layer. Usually, the effect is irrelevant,
C         i.e. both are set equal to VPI.

      real EEV,VO,VV,VPI,VPIS,VPIO

C  local variable

C  workfn: Work function of the LEED filament material. Theoretical predictions
C         for the inner potential usually have their energy zero fixed at the
C         Fermi level; the experimental energy scale (-> EEV) nominally
C         does the same, except for the fact that electrons need to overcome
C         the work function of the cathode first. Note that the this value is
C         thus formally determined by a LEED fit, yet don't trust its accuracy.

      real workfn

C  set work function of cathode
C  work function should be positive (added to exp. energy EEV)

      workfn = 2.65

C  set real part of inner potential

      VV = workfn-max(-10.85, ((0.00)+(-67.83)/sqrt(EEV+workfn+(13.72)))&
     &)

      write(6,*) workfn, EEV
      write(6,*) VV

c  set difference between bulk and overlayer inner potential

      VO = 0.

c  set imaginary part of inner potential - energy independent value used here

      VPI = 5.0

C  set substrate / overlayer imaginary part of inner potential

      VPIS = 5.0
      VPIO = 5.0

      return
      end