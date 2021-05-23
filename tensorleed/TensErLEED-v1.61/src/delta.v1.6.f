C  Tensor LEED delta amplitude calculation v1.2
C  VB, 13.04.00 
C  for use with lib.tleed.f v1.2, lib.delta.f v1.2 
C
C  as described in 
C
C  V. Blum, K. Heinz, submitted to Comp. Phys. Comm. (2000).
C
C  Please cite this reference when using this program.
C
C  Welcome to the Erlangen Tensor LEED program package. Permission is granted
C  to use these programs freely, and to modify them as you may find appropriate.
C  However, we cannot take any responsibility for possible bugs or other
C  limitations you may encounter. If you find a bug or worthwhile improvement,
C  please notify us at 
C
C                       tleed@fkp.physik.uni-erlangen.de
C
C  so we can make the updated code available to everyone. Likewise, in re-
C  distributing the code please refer back to us first to ensure the latest version
C  of the package is passed on.
C
C  What the program does:
C  * performs geometric, thermal and stoichiometric perturbations of the LEED
C    amplitudes at a single atomic site of a given reference structure, yielding
C    the appropriate amplitude differences

C  Based on Uli Loeffler's version of thermal Tensor LEED. Does not include 
C  anisotropic and anharmonic vibrations since both effects have too little
C  impact on spectra except in special cases. In those cases, please ask for an
C  appropriate version (note, however, that slight adaptions of that version may
C  be necessary).
C

      program delta

C  include necessary parameter for dimension statements

      include "PARAM"

C  set unchangeable constants

      include "GLOBAL"

C  MNDIM determines up to which LMAX variable dimensions for
C        Clebsch-Gordan coefficients are known to the program
C        must only be changed if data statement for NLMBS
C        is enlarged appropriately

      PARAMETER (MNDIM = 15)

C  set derived parameters

      PARAMETER (MLMAX1 = MLMAX+1)
      PARAMETER (MN = 2*MLMAX +1)
      PARAMETER (MNN = MN*MN)
      PARAMETER (MNN1 = MN, MNN2 = MLMAX1, MNN3 = MLMAX1)
      PARAMETER (MLMMAX = MLMAX1*MLMAX1)

      PARAMETER (MNFAC = 4*MLMAX + 3)
      PARAMETER (MNTOT = MNCSTEP*MNDEB)

C  variable and array declarations in order of appearance

C  variables derived from constants set before

C  EMACH: machine accuracy (variable derived from constant MEMACH)
C  NPSI : number of energy values for which phase shifts are tabulated in input
C  NEL  : number of elements for which phase shift values are tabulated in input
C  LMAX : maximum angular momentum to be used in calculation
C  L1,LMAX1,LMMAX,NFAC: derived from LMAX
C  NT0  : no. of TLEED output beams MNT0
C  NATOMS: Currently, NATOMS = 1 is the only possible choice. In principle, NATOMS is the
C         number of atoms in a reconstructed unit cell constructed from a non-reconstructed
C         layer of the reference calculation. This feature is currently not supported.
C  NCSTEP: number of geometric variations ('displacements') to be considered
C  NDEB : no. of thermal variations (different thermal vib. amps) to be considered
C  NTOT : total number of variations NDEB * NCSTEP

      REAL EMACH
      INTEGER NPSI,NEL
      INTEGER LMAX,L1,LMAX1,LMMAX,NFAC
      INTEGER NT0,NATOMS,NCSTEP,NDEB,NTOT

C  subroutine readbas

C  THETA,FI: incident angles of incoming electron beam onto surface
C  RAR1,RAR2: first and second reciprocal lattice vectors of substrate lattice, resp.
C  ES    : energy values for which phase shifts PHSS are tabulated in input
C  PHSS  : phase shift values for MNEL elements, MNPSI energy steps ES
C  PQFEX : Tensor LEED output beam list PQFEX (currently redundant since all TLEED beams
C          must be identical to fd. beams)
C  FORMIN:  0: Tensor file AMP is binary, 1: tensor file AMP is ascii
C  FORMOUT: 0: Delta amp. file DELWV is binary, 1: Delta amp. file DELWV is ascii 
C  CUNDISP: position of current atomic site in reference calculation - CDISP is
C  CDISP : displaced positions of current atomic site for variation
C  IEL   : element no. (in phase shifts supplied with input) that delta amplitudes
C          will be calculated for (not necessarily the same element as the one
C          used in the reference calculation!) - IEL = 0 means a vacancy will be assumed
C  DR0_A,DRPER_A,DRPAR_A: thermal vibration amplitude steps to be included in 
C          current variation - DR0_A is always 0., DRPER_A = DRPAR_A forced in readin
C  EI,EF: lower and upper limiting energies of energy range to be considered.
C  TV:    area of (overlayer) lateral unit cell - in case TLEED wrt smaller unit cell
C         is used, TVA from reference computation must be set.

      REAL THETA,FI
      REAL RAR1,RAR2
      DIMENSION RAR1(2),RAR2(2)
      REAL ES
      DIMENSION ES(MNPSI)
      REAL PHSS
      DIMENSION PHSS(MNPSI,MNEL,MLMAX1)
      REAL PQFEX
      DIMENSION PQFEX(2,MNT0)
      INTEGER FORMIN, FORMOUT
      REAL CUNDISP
      DIMENSION CUNDISP(MNATOMS,3)
      REAL CDISP
      DIMENSION CDISP(MNCSTEP,MNATOMS,3)
      INTEGER IEL
      REAL DR0_A,DRPER_A,DRPAR_A
      DIMENSION DR0_A(MNDEB),DRPER_A(MNDEB),DRPAR_A(MNDEB)
      REAL EI,EF
      REAL TV

C  subroutine OUTDELT_DWG

C  CDISP_DWG: aux. storage array for CDISP 
C  AID_DWG: aux. array for output of thermal vib. amps used in comp. for later
C           identification in delta-amp. file

      REAL CDISP_DWG,AID_DWG
      DIMENSION CDISP_DWG(MNTOT,MNATOMS,3),AID_DWG(MNTOT)

C  subroutines BELMG, CPPP

C  NLMBS : dimension of BELM as a function of LMAX
C  BELM  : Clebsh-Gordon coefficients for tmatrix(), generated by BELMG
C  NN1,NN2,NN3: variable dimensions for PPP
C  PPP   : Clebsh-Gordon coefficients for computation oftemperature-
C          dependent phase shifts
C  FAC   : auxiliary array for subroutines BELMG,CPPP

      INTEGER NLMBS
      DIMENSION NLMBS(MNDIM)
      REAL BELM
      DIMENSION BELM(MNLMB)
      INTEGER NN1,NN2,NN3
      REAL PPP
      DIMENSION PPP(MNN1,MNN2,MNN3)
      DOUBLE PRECISION FAC
      DIMENSION FAC(MNFAC)

C  variables for energy loop

C  subroutine INDATA

C  E    : computational energy inside crystal
C  CAF  : atomic t-matrix of current site as used in reference calculation
C  XIST : output amplitudes from reference calculation (only tleed beams)
C  VPIS,VPIO: imaginary part of the inner potential, substrate & overlayer, resp.
C  VO,VV: real part of the inner potential, overlayer & substrate, resp.
C  EEV :  current energy in eV

      REAL E
      COMPLEX CAF
      DIMENSION CAF(MLMAX1)
      COMPLEX XIST
      DIMENSION XIST(MNT0)
      REAL VPIS,VPIO,VO,VV
      REAL EEV

C  subroutine INAMP

C  PSQ  : lateral momentum of Tensor LEED beams relative to
C         incident beam (0,0)
C  AK2M,AK3M: (negative) absolute lateral momentum of Tensor LEED beams
C         (for use as incident beams in time-reversed LEED calculation)
C  ALM  : spherical wave amplitudes incident on current atomic site in reference calculation
C         (i.e., scattering path ends before scattering on that atom)
C  EXLM : spherical wave amplitudes incident from exit beam NEXIT in "time-reversed"
C         LEED experiment (or rather, all terms of Born series immediately after
C         scattering on current atom)
C  GTEMP: used here as a working space only - really designed for use as a working
C         space in subroutine tmatrix_dwg

      REAL PSQ
      DIMENSION PSQ(2,MNT0)
      REAL AK2M,AK3M
      DIMENSION AK2M(MNT0),AK3M(MNT0)
      COMPLEX ALM,EXLM
      DIMENSION ALM(MLMMAX),EXLM(MLMMAX,MNT0)
      COMPLEX GTEMP
      DIMENSION GTEMP(MLMMAX,MLMMAX)

C  subroutine TSCATF

C  VSITE : possible energy shift in phase shift computations - can be used to describe 
C          local variations of the muffin-tin-constant
C  T0    : must equal T if input vib. amplitudes are to be used properly - not 0. !!
C  T     : must equal T0 if input vib. amplitudes are to be used properly - not 0. !!
C          if T0 is set to the temperature that the vibs. were computed for, T could
C          in principle be used to simulate the temperature behaviour of 
C          a Debye-like phonon spectrum. Yet, this simply alters the vib. amplitude used 
C          for the DW factor, thus it only makes sense to either vary DRPER or T.
C  AF    : auxiliary array for TSCATF - here, tscatf stores the zero-temperature t-matrix
C          (no vibs at all) of current element IEL
C  NewCAF: working array in which current (displaced) atomic t-matrix is stored
C  PHS,DEL,CTAB,SUM,BJ: working spaces for TSCATF and PSTEMP - BJ also in MATEL_DWG!

      REAL VSITE,T0,T
      COMPLEX AF,NewCAF
      DIMENSION AF(MLMAX1),NewCAF(MLMAX1)
      REAL PHS
      DIMENSION PHS(MLMAX1)
      COMPLEX DEL,CTAB,SUM
      DIMENSION DEL(MLMAX1),CTAB(MNN3),SUM(MNN2)
      COMPLEX*16 BJ
      DIMENSION BJ(MNN1)

C  subroutine MATEL_DWG

C  DELWV : working space for computation and storage of amplitude differences
C  DELTAT: working space - storage of change in t-matrix through desired displacements
C  GTWOC,SYLM : other working spaces
C  NRATIO: originally, ratio between substrate and overlayer unit cell area. However,
C          currently all TLEED parts must be performed with overlayer symmetry only,
C          thus NRATIO = 1 and TV = TVB are the only safe choice
C  FAC1,FAC2,FAC3: auxiliary arrays for the computation of spherical harmonics in
C          subroutine sphrm4
C  LMAX21,LMMAX2: variable dimensions for FAC1,FAC2,FAC3

      COMPLEX DELWV
      DIMENSION DELWV(MNCSTEP,MNT0)
      COMPLEX DELTAT
      DIMENSION DELTAT(MLMMAX,MLMMAX)
      COMPLEX GTWOC,SYLM
      DIMENSION GTWOC(MLMMAX,MLMMAX),SYLM(MNN)
      INTEGER NRATIO
      REAL*8 FAC1,FAC2,FAC3
      DIMENSION FAC1(MN),FAC2(MNN),FAC3(MN)
      INTEGER LMAX21,LMMAX2

C  subroutine COPDELWV

C  DELWV_DWG: storage array for delta amplitudes, for all vibrational and geometric
C         variations (DELWV is designed to do this but holds only the geometric var.)

      COMPLEX DELWV_DWG(MNTOT,MNT0)

C  data statement

C  NLMBS is dimension of BELM as function of LMAX

      DATA NLMBS(1),NLMBS(2),NLMBS(3),NLMBS(4),NLMBS(5),NLMBS(6),
     *     NLMBS(7),NLMBS(8),NLMBS(9),NLMBS(10),NLMBS(11),NLMBS(12),
     *     NLMBS(13),NLMBS(14),NLMBS(15)
     *     /19,126,498,1463,3549,7534,14484,25821,43351,69322,106470,
     *      158067,227969,320664,441320/

C  common blocks

C  end header - begin actual computation

      EMACH  = MEMACH

C  set integers required for variable dimensions throughout code
C  (cf. PARAM)

      NPSI   = MNPSI
      NEL    = MNEL
      LMAX   = MLMAX
      L1     = LMAX+1
      LMAX1  = LMAX+1
      LMMAX  = MLMMAX
      NFAC   = MNFAC
      NT0    = MNT0
      NATOMS = MNATOMS
      NCSTEP = MNCSTEP
      NDEB   = MNDEB
      NTOT   = MNTOT
      LMAX21 = MN
      LMMAX2 = MNN

C  read information on reference calculation, and information on required perturbations
C  from stdin (unit 5)

      call ReadBas(EI,EF,RAR1,RAR2,TV,THETA,FI,FORMIN,PQFEX,NT0,
     +             ES,PHSS,NPSI,NEL,L1,LMAX,FORMOUT,IEL,CUNDISP,
     +             CDISP,NATOMS,NCSTEP,DR0_A,DRPER_A,DRPAR_A,
     +             NDEB)

C  read general information from tensor storage file AMP, unit 1

      CALL OPENIN(1,'AMP',FORMIN)

C  write header information for delta amplitude file DELWV, unit 10

      CALL OUTDELT_DWG(10,'DELWV',FORMOUT,NT0,PQFEX,THETA,FI,RAR1,
     1                     RAR2,NATOMS,CUNDISP,NCSTEP,CDISP,CDISP_DWG,
     2                     AID_DWG,NDEB,DR0_A,DRPER_A,DRPAR_A)

C  prepare calculational quantities for upcoming energy loop

      NLMB = NLMBS(LMAX)

      CALL BELMG(BELM,NLMB,LMAX,NFAC,FAC)

C  Calculate Clebsch-Gordan-Coeff. PPP for PSTEMP for use in TSCATF

      NN3=LMAX+1
      NN2=LMAX+1
      NN1=NN2+NN3-1

      CALL  CPPP (PPP, NN1, NN2, NN3, NFAC, FAC)

C  Default settings for energy loop (not all possible options of subroutines
C  are currently exploited)

C  in tscatf

      VSITE = 0.
      T0 = 100.
      T  = 100.

C  in matel_dwg

      NRATIO = 1      

C  note TV = TVB in readbas

C  read first energy, inner potential, momentum, reference t-matrix,
C  reference amplitudes from Tensor file AMP

      CALL INDATA(1,FORMIN,E,L1,CAF,NT0,XIST,VPIS,VPIO,VO,VV)

      IF (E.lt.0.) THEN
        write(6,*) "First energy in Tensor file AMP was zero. Stopped."
        STOP
      END IF

C  loop 1000: while E .gt. 0 perform loop over energies

 1000 CONTINUE

      EEV = (E - VV) * HARTREE

C  read Tensor wavefields for current energy from file AMP; GTEMP is
C  used as a working space with improper, but sufficient dimensions

          CALL INAMP(1,FORMIN,NT0,PSQ,AK2M,AK3M,ALM,EXLM,LMMAX,
     1               GTEMP,LMMAX)


      IF ((EEV.gt.(EI-1.E-3)).and.(EEV.lt.(EF+1.E-3))) THEN

        write(6,*) "Current energy ",EEV," eV within given ",
     +  "energy range - performing calculation."

C  loop over requested thermal vibration amplitudes

        DO IDEB = 1,MNDEB

C  create modified NewCAF from phase shifts via DW-correction
C  use only vibrational amplitude no. IDEB of element no. IEL
C  if IEL = 0, assume a vacancy on current site!

          IF (IEL.ne.0) THEN

            CALL TSCATF(IEL,L1,ES,PHSS,NPSI,E,VSITE,PPP,NN1,NN2,NN3,
     +      DR0_A(IDEB),DRPER_A(IDEB),DRPAR_A(IDEB),T0,T,AF,NewCAF,NEL,
     +      LMAX1,PHS,DEL,CTAB,SUM,BJ)

          ELSE

            DO IL = 1,MLMAX1
              NewCAF(IL) = 0.
            ENDDO

          END IF

C  now calculate ordinary geometries with NewCAF, just produced

         CALL MATEL_DWG(DELWV,NCSTEP,CAF,NewCAF,BJ,BELM,NLMB,
     1                DELTAT,GTWOC,SYLM,E,VV,VPIS,LMAX,L1,LMMAX,NT0,
     2                EXLM,ALM,AK2M,AK3M,NRATIO,TV,LMAX,L1,LMMAX,
     3                NATOMS,CDISP,CUNDISP,PSQ,GTEMP,FAC1,FAC2,FAC3,
     4                LMAX21,LMMAX2)

C  store the result of the geometric variation, DELWV, in a proper storage array
C  that includes the thermal vibration loop - DELWV_DWG

        CALL COPDELWV(DELWV_DWG,DELWV,IDEB,NCSTEP,NDEB,NT0)

        ENDDO

C  write results to output file DELWV

        CALL OUTRINT(10,FORMOUT,E,VV,VO,VPIS,NT0,NTOT,XIST,DELWV_DWG)

      ELSE

        write(6,*) "Current energy ",EEV," eV outside given ",
     +  "energy range - skipping calculation."

      END IF

C  read next energy, inner potential, momentum, reference t-matrix,
C  reference amplitudes from Tensor file AMP

      CALL INDATA(1,FORMIN,E,L1,CAF,NT0,XIST,VPIS,VPIO,VO,VV)

C  terminate only if negative E value is read! Otherwise, re-enter energy loop.

      IF (E.ge.0.) GOTO 1000

C  when arrived here, tensor file is finished.

      write(6,*) "Correct termination assumed."

      END
