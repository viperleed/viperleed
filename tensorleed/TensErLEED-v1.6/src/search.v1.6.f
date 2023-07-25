C  Tensor LEED optimization algorithm 
C  v1.2 (search version v106), VB 13.04.00
C  for use with lib.search.f v1.2, random_.c
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


**************************************************************************

C  original Author M. Kottcke

C  current version v106

C  Version R. Backofen, V. Blum, 06.09.95
C  A. Seubert v90 (Including capability of incoherent domain averaging; 5/97)
C  V. Blum v91: now a truly random sampling search algorithm; a constant
C  offset is added to the gaussian probability distribution that
C  determines the next parameter value. The slight distortion of the
C  gaussian distribution ensures that every new parameter value can be
C  reached from a certain starting point, despite the discrete nature
C  of the random values used. 15.5.97
C  V. Blum v92: handling of area fractions corrected so that now analogous 
C  to other parameters. C source code for random routine now uses random()
C  instead of rand(). Initialisation of random done by randominit (uses 
C  srandom(). (cf man-pages). R-Faktor subroutine corrected so that 
C  beam weights WB from WEXPEL are used correctly now. 21.10.97  
C  V. Blum v100: entire run now constructed by shell script, including
C  parameter restrictions (subroutine restrict) and r-factor control
C  file. WEXPEL no longer in use - now up-to-date-file rf.info. Note
C  that the obsolete quantity NSS is still dragged through the code. Anyone
C  wishing to become a hero should change that. RMUT now in search.steu.
C  16.12.97
C  V. Blum v101: enhanced portability (now works for Linux g77)
C  V. Blum v102: corrected computation of domain fractions
C  V. Blum, A. Seubert v103: Minor adjustments (e.g., corrected WB declaration)
C  V. Blum v104a: Minor corrections (type declarations adjusted)
C  V. Blum v105: energy dependent inner potential now taken over from LEED
C                calculation
C  V. Blum v106: minor adjustments for TensErLEED; subroutines now in lib.search.f

**************************************************************************

C  Welcome stranger!

C  This program searches for a global minimum of the r-factor of a given 
C  number of trial stuctures. The only algorithm implemented so far performs
C  random sampling with a downstep constraint. The program was written by
C  M. Kottcke.
C
C  It allows for the use of R-Pe or R2 during this search using different
C  beam groups, and comparison with theoretical spectra in the van Hove
C  format.
C
C  An energy-dependent inner potential of any shape can be used in the
C  search.
C
C  Chemistry can be varied for each atom individually. Note, however, that
C  the sum of concentrations for each place must then equal 1. In the case
C  of vacancies, a delta amplitude file must therefore also be generated
C  for the vacancy.
C
C  The treatment of the surface parameter is slightly different from usual
C  TLEED. Also, it will only work properly if the use of CUNDISP is taken
C  up in its original form - this is currently not the case. 
C
C  Much of this program has grown from LEED and TLEED codes developped by
C  generations of LEED scientists. At least the main program is now fully
C  documented, and all variables have been declared. This may not seem
C  necessary in many cases, however on the whole it improves the 
C  understanding of the purpose of each part greatly. Therefore please
C  comment each change, remove old comments immediately and declare 
C  and document each new variable in the program header immediately.
C
C  The generations of physicists to come in the future, boldy going where
C  no LEED analysis has gone before, will be indebted to you eternally.

*************************************************************************

C 

      PROGRAM SEARCH

C  Include global parameters for dimension statements etc.

      INCLUDE "PARAM"
      INCLUDE "GLOBAL"

*******************************************************************

C  Variables in order of appearance - please introduce ALL variables
C  here except for loop counters!!

********************************************************************

C  Variables for later use in dimension statements in subroutines

      INTEGER NSTEP,NFILES,NPRMK,NPS,NT0,NBED,NDATA
      INTEGER NCONCS,NNATOMS,NDATT,NBTD,NBMD, NDOM
      integer NPLACES(MNDOM)

C --- NPRAS(NDOM) is maximum number of parameters for each domain (cf. NPRMK)

      integer NPRAS(MNDOM)

C  Global variables used directly in search procedure

C  INIT is possible initialization for random number generator; INIT = 0 -> use system time
C       INIT > 0 -> init is used as initialising value for random, useful when testing

      INTEGER INIT

C  PARIND contains parameter value for each parameter (incl. concentration) and each
C         individual in current generation
C  PAROLD is similar for last generation

      INTEGER PARIND,PAROLD
      DIMENSION PARIND(MNPRMK,MPS),PAROLD(MNPRMK,MPS)


C  Variables for r-factor determination obtained from input (WEXPEL)

C  EMIN, EMAX are boundaries of energy range to be considered
C  EINCR is energy increment of interpolation grid for experimental as well
C        as theoretical data
C  IPR controls size of output during r-factor calculation such as Y-functions etc.
C  VI is optical potential value in eV of calculation (must be identical with 
C     VPI read from delta amp files which is in hartrees.)
C  V0RR is value of inner potential used in calculation (should be identical
C       with VV value read from delta amp files)
C  V01,V02 are lower and upper boundary of inner potential variation relative to
C          V0RR during r-factor calculation
C  VINCR is step width of this variation.
C  ISMOTH is number of smoothing steps of experimental data
C  EOT determines whether the data to compare search results to are experimental
C      or theoretical data
C  IBP contains the number of the theoretical beam IBP(IBE) that is compared to
C      experimental beam IBE
C  WB is weight of each exp beams in calculation of beam-averaged r-factors
C  KAV is averaging scheme for theoretical beams
C  NSS is 1 for all purposes here but passed on to too many other subroutines
C  BENAME is name of some beams and dummy array
C  MITTEL contains arrangement of exp beams into 2 beam groups, e.g. integer
C         and half order beams
C  OVLG is total energy range (sum of all used experimental beams)

      REAL EMIN,EMAX,OVLG
      REAL EINCR
      INTEGER IPR
      REAL VI,V0RR,V01,V02,VINCR
      INTEGER ISMOTH
      INTEGER EOT
      INTEGER IBP,KAV
      DIMENSION IBP(MNBED),KAV(MNBTD)
      REAL WB
      DIMENSION WB(MNBED)
      INTEGER NSS
      DIMENSION BENAME(5,MNBED)
      INTEGER MITTEL
      DIMENSION MITTEL(MNBED)
      

C  Arrays for data to compare search results to from WEXPEL via READE or READT

C  AE is array for measured experimental intensity data (or theor. comparison data
C     if EOT=1.) 
C  EE is array with measurement energy grid
C  NEE is number of experimental energies in each beam
C  NBEA contains information for averaging of experimental beams

      REAL AE
      DIMENSION AE(MNBED,MNDATA)
      REAL EE
      DIMENSION EE(MNBED,MNDATA)
      INTEGER NEE
      DIMENSION NEE(MNBED)
      INTEGER NBEA
      DIMENSION NBEA(MNBED)


C  Variables only used in READT, not in READE

C  MAXI is maximum intensity of each beam
C  IOFF is number of first non-zero intensity value for each beam

      REAL MAXI
      DIMENSION MAXI(MNBTD)
      INTEGER IOFF(MNBTD)


C  Variables calculated from input 'experimental' data in PREEXP

C  YE are Pendry Y-functions of averaged and smoothed experimental I(E) data
C     for each beam
C  TSE is integral value over exp data in R2 calculation for each beam
C  TSE2 is integral value over (squared) exp data in R2 calculation for each beam
C  TSEY2 is integral value over exp data in RPe calculation for each beam
C  XPL, YPL are arrays for later use in INTPOL, EXPAV, VARSUM that
C           must be stated globally simply to allow for variable dimensions 
C           in subroutines   (cheers, f77)
C  AEP contains 1st derivative of experimental data
C  NNN is array for later use in EXPAV that needs variable dimensions
C  NBE is number of experimental beams after averaging

      REAL YE
      DIMENSION YE(MNBED,MNDATA)
      REAL TSE,TSE2
      DIMENSION TSE(MNBED),TSE2(MNBED)
      REAL TSEY2
      DIMENSION TSEY2(MNBED)
      REAL XPL,YPL
      DIMENSION XPL(MNDATA),YPL(1,MNDATA)
      REAL AEP
      DIMENSION AEP(MNBED,MNDATA)
      INTEGER NNN
      DIMENSION NNN(MNDATA)
      INTEGER NBE


C  Variables for search procedure itself, obtained from file steu
C  via subroutine READSC

C  NPAR contains the number of independent parameters (as opposed to NPRMK)
C  INFILE contains delta amp file names to be used for calculation
C  NSURF determines whether current file influences the uppermost surface
C        of the crystal considered, i.e. if height of surface-vacuum interface
C        is changed due to buckling in this file
C  IFORM states whether or not a delta amp file input is formatted
C  PNUM is total number of parameters (including conc steps), must be equal MNPRMK
C  STAFLA determines whether random or given start configuration is used
C  WHICHG is flag to optimize for integer, half-order or total R-factor
C  WHICHR decides whether R2 or RPe is used for optimization
C  VARST is array containing number of grid points for each parameter
C  PARTYP is array containing the number of different parameters in each file
C  OUTINT forces output of current search data to SEADOC after OUTINT generations 
C         regardless of improvements made
C  FILREL is atom number that forces two atoms to be varied together,
C         i.e. with identical parameter choice only
C  NFIL(IDOM,IPLACE) is number of files for atomic site IPLACE in domain IDOM
C  CONC is array of concentrations to be used
C  DMISCH Schrittweite fuer Bildung von Mischungsverhaeltnissen (in Prozentpunkten)
C  MAXGEN: max. number of generations to be performed
C  SEANAME: name of search document file

      INTEGER NPAR
      CHARACTER*15 INFILE(MNDOM, MNPLACES, MNFILES)
      INTEGER NSURF,IFORM
      DIMENSION NSURF(MNDOM, MNPLACES), IFORM(MNDOM, MNPLACES, MNFILES)
      INTEGER PNUM,STAFLA,WHICHG,WHICHR
      INTEGER VARST,PARTYP
      DIMENSION VARST(MNPRMK),PARTYP(MNDOM, MNPLACES,MNFILES)
      INTEGER OUTINT,FILREL
      DIMENSION FILREL(MNDOM, MNPLACES)
      INTEGER NFIL
      DIMENSION NFIL(MNDOM, MNPLACES)
      REAL CONC
      DIMENSION CONC(MNDOM, MNCONCS,MNPLACES,MNFILES)
      REAL DMISCH
      integer MAXGEN
      character*10 SEANAME

C  Variables used to read delta amplitude files in ReadFile

C  CNTFIL is file number during readin of delta amplitudes
C  THETA, FI are determine angle of incident beam w.r.t. surface 
C            under consideration, in radians (not degrees!)
C  RAR1,RAR2 are reciprocal lattice vectors of surface
C  PQFEX is beam list
C  CUNDISP is undisplaced position of atom under consideration (probably unused)
C  CDISP are displaced positions of that atom in variation
C  AID is array that has a bright past in the history of TLEED and currently no future.
C  EMK,EMK0 contain energy steps in hartrees for computation of intensities
C           (EMK0 for test reasons only).
C  ESMK contains energy steps in eV for use in r-factor calculation and can be 
C       corrected for shifts of inner potential.
C  VPI is optical potential in hartrees
C  VV is inner potential value
C  VO is possible correction for overlayer inner potential
C  XISTMK is array to store amplitudes of beams for the undisplaced surface
C  DELMK contains delta amplitudes for all displacements for all beams, files, 
C        and places.
C  DATTNO is number of theoretical energies in interval [EMIN,EMAX] from WEXPEL
C         (only these are used in calculation)
C  DELWV,XIST are no longer used, they are simply dummy arrays for unformatted 
C             read-in (variable dimensions!)
C  CUNDVB,CDVB are for unformatted read-in (variable dimensions!)
C  quot is used for energy-dependent inner potential V0R=const+E/quot

      INTEGER CNTFIL
      REAL THETA,FI
      REAL RAR1,RAR2
      DIMENSION RAR1(2),RAR2(2)
      REAL PQFEX
      DIMENSION PQFEX(2,MNBTD)
      REAL CUNDISP,CDISP,AID
      DIMENSION CUNDISP(MNATOMS,3,MNDOM,MNPLACES,MNFILES)
      DIMENSION CDISP(MNCSTEP,MNATOMS,3,MNDOM,MNPLACES,MNFILES)
      DIMENSION AID(MNCSTEP,MNDOM,MNPLACES,MNFILES)
      REAL EMK0,EMK,ESMK
      DIMENSION EMK0(MNDATT),EMK(MNDATT),ESMK(MNDATT)
      REAL VPI,VO,VV
      DIMENSION VPI(MNDATT),VO(MNDATT),VV(MNDATT)
      COMPLEX XISTMK(MNDOM,MNDATT,MNBTD)
      COMPLEX DELMK(MNDATT,MNCSTEP,MNBTD,MNDOM,MNPLACES,MNFILES)
      INTEGER DATTNO
      COMPLEX DELWV(MNCSTEP,MNBTD)
      COMPLEX XIST(MNBTD)
      REAL CUNDVB,CDVB
      DIMENSION CUNDVB(MNATOMS,3),CDVB(MNCSTEP,MNATOMS,3)

      REAL quot


C  Variables for search procedure in main program

C  RMUT is constant that normalizes width of gaussian distribution 
C       to a useful value
C  AVERNEW is average search r-factor of current generation
C  RPEIND are current r-factors that are used for optimization, for each individual
C  WIDT is parameter for width determination of gaussian distribution determined
C        from current and previous generation, for each parameter
C  BV0 is optimized inner potential value
C  BRGES is total rfactor of current population 
C  BRINS is the integer rfactor
C  BRHAS is the half-order rfactor
C  ending 'OLD' means the appropriate parameter for the previous calculation
C  TCHANGE checks whether or not r-factor has improved since last generation
C  NDSL1, Norm,PARSUM for PMOLD initialisation

      REAL RMUT,AVERNEW
      REAL RPEIND,RPEOLD
      DIMENSION RPEIND(MPS),RPEOLD(MPS)
      REAL WIDT
      DIMENSION WIDT(MNPRMK)
      REAL BV0,BV0OLD
      DIMENSION BV0(MPS),BV0OLD(MPS)
      REAL BRGES,BRINS,BRHAS
      DIMENSION BRGES(MPS),BRINS(MPS),BRHAS(MPS)
      REAL BRGESOLD,BRINSOLD,BRHASOLD
      DIMENSION BRGESOLD(MPS),BRINSOLD(MPS),BRHASOLD(MPS)
      LOGICAL TCHANGE
      INTEGER PARSUM,NDSL1
      REAL Norm


C  Variables for randomizing of parameter values

C  WSK is gaussian probability distribution plus a small constant offset
C      in sea_rcd, only used there, but needs variable dimensions 

      REAL WSK
      DIMENSION WSK(MNCSTEP)

C  Variables used in processing parameter values for use in GetInt
C  and only defined for current individual (subroutine GetGrid, Mischung)

C  IFNUM is delta amplitude number in current delta amp file
C  NPARC stores away current concentration step number for each place
C        for easier access
C  PMISCH are parameter step numbers for domain fractions, computed from
C         PARIND so that normalisation to 1 is fulfilled
C  PMOLD  is storage array for last PMISCH values

      INTEGER IFNUM
      DIMENSION IFNUM(MNDOM,MNPLACES,MNFILES)
      INTEGER NPARC
      DIMENSION NPARC(MNDOM,MNPLACES)
      INTEGER PMISCH,PMOLD
      DIMENSION PMISCH(MNDOM,MPS),PMOLD(MNDOM,MPS)

C  Variables for computation of intensities from delta amplitudes in GetInt

C  ATSMK are theoretical output intensities for all beams
C  ATSAS has the same purpose for each domain

      REAL ATSMK, ATSAS
      DIMENSION ATSMK(MNBTD,MNDATT), ATSAS(MNDOM, MNBTD, MNDATT)


C  Variables for calculation of r-factor (subroutine Rfaktor)
C  All these are local variables with variable dimensions.

C  NETI is number of theoretical energies, historical reasons
C  PQ used to be theoretical beam identification - currently unused
C     but kept for various output routines
C  SYM used to store symmetry properties of beams but is useless in tensor LEED
C  NET is no of points in each theoretical beam after interpolation 
C  AT are theoretical intensities after interpolation to experimental grid
C  ET is working space containing part of energy array ESMK
C  ATP is first derivative of theoretical intensities for each beam
C  YT are theoretical Pendry y-functions for each beam
C  TSTY2 is an integral over theoretical y-function
C  ARM,ARPEM,ERANGM,RAZZM,RANNM,RAVPM,XRPEM are various r-factors
C        and energy ranges averaged over beam groups defined by MITTEL  
C  NST1, NST2 are number of data points outside theory-exp overlap
C        for each beam
C  R2 is obviously R2.
C  RPE is Pendry r-factor for each beam
C  EET is overlap between experiment and theory for each beam
C  TST is an integral over theoretical intensities

      INTEGER NETI
      REAL PQ
      DIMENSION PQ(2,MNBTD)
      INTEGER SYM
      DIMENSION SYM(MNBTD)
      INTEGER NET
      DIMENSION NET(MNBTD)
      REAL AT,ET,ATP
      DIMENSION AT(MNBTD,MNDATA),ET(MNBTD,MNDATA),ATP(MNBTD,MNDATA)
      REAL YT
      DIMENSION YT(MNBTD,MNDATA)
      REAL TSTY2
      DIMENSION TSTY2(MNBTD)
      REAL ARM(2),ARPEM(2),ERANGM(2),RAZZM(2),RANNM(2)
      REAL RAVPM(2),XRPEM(2)
      INTEGER NST1,NST2
      DIMENSION NST1(MNBMD),NST2(MNBMD)
      REAL R2,RPE,EET
      DIMENSION R2(MNBED),RPE(MNBED),EET(MNBED)
      REAL TST
      DIMENSION TST(MNBTD)

C  Variables for adoption of steepness factor in GetWid

C  NormInd is normalisation that unfortunately would like a variable
C          dimension

      INTEGER NormInd(MPS)

C Some auxiliary variables used for search document

      INTEGER   PARHELP
      DIMENSION PARHELP(MNPRMK,MPS)
      REAL      RPEHELP
      DIMENSION RPEHELP(MPS)
      REAL      BRGESHELP,BRINSHELP,BRHASHELP,BV0HELP
      DIMENSION BRGESHELP(MPS),BRINSHELP(MPS),BRHASHELP(MPS),
     +          BV0HELP(MPS)

***************************************************************************

C  end variable declarations and do something useful now ...

***************************************************************************

C  first thing: open control file

      OPEN(8,file='control.chem')

C  set dimension statements for subroutines to desired values
C  generally, none of these should be set too large in PARAM

C  NDATA,NDATT are possibly difficult dimension statements
C  both will be repeated later

      NSTEP = MNCSTEP
      NFILES = MNFILES
      NPRMK = MNPRMK
      NPS = MPS
      NCONCS = MNCONCS
      NBED=MNBED
      NNATOMS=MNATOMS
      NDATA=MNDATA
      NDATT=MNDATT
      NDOM = MNDOM

      NT0=MNBTD
      NBTD=MNBTD

      NBMD = MNBMD

C  initialize population here (may change in readsc)

      DO 1849 IPOP=1,MPS
      DO 1849 IPARAM=1,MNPRMK
         PAROLD(IPARAM,IPOP)=1
 1849    PARIND(IPARAM,IPOP)=1

C  initialize random function (is done in C, using system time)
C  only use randominit if random() is used

      call randominit(INIT)

C Modul 1: READIN INFORMATION FOR rfactor determination from WEXPEL,
C          such as beam grouping, energy ranges etc.

      CALL READRF(EMIN,EMAX,EINCR,
     +            IPR,VI,V0RR,V01,V02,VINCR,ISMOTH,EOT,
     +            IBP,WB,NBTD,NBED,KAV,NSS,MITTEL)


C Modul 2: READIN EXPERIMENTAL

C read data to be compared, either experimental or theoretical

      IF (EOT.eq.1.) THEN

        CALL READT(AE,EE,NBED,NEE,NBEA,NDATA,BENAME,IPR,MAXI,IOFF)

      ELSE

        CALL READE(AE,EE,NBED,NEE,NBEA,BENAME,IPR)

      END IF

C  end readin of exp. or theor. reference data


C Modul 3: PREPARE EXPERIMENTAL DATA FOR LATER USAGE

      CALL PREEXP(AE,EE,NBED,NEE,BENAME,NBEA,IPR,ISMOTH,
     +            EINCR,VI,YE,NDATA,TSE,TSE2,TSEY2,XPL,YPL,AEP,NNN,NBE)


C Modul 4: Readin control information for search algorithm

      CALL READSC(NDOM,NPLACES,NFILES,INFILE,NSURF,IFORM,
     +            PNUM,VARST,NPRMK,NPRAS,PARTYP,NPS,PARIND,STAFLA,
     +            OUTINT,FILREL,WHICHG,WHICHR,NFIL,NCONCS,CONC,DMISCH,
     +            MAXGEN,SEANAME,NPAR,RMUT,INIT)




C  Modul 5: Read in delta amplitudes from files

C  Now read in delta amplitudes from files, skipping those outside
C  energy range [EMIN,EMAX]

C  make sure dimensions are correct

      NDATT = MNDATT

C  CNTFIL is current file number

      CNTFIL = 10

      DO 2221 IDOM = 1, NDOM

        DO 2222 IPLACE = 1,NPLACES(IDOM),1

          DO 2223 IFILE = 1,NFIL(IDOM,IPLACE),1

            CNTFIL = CNTFIL + 1

            CALL ReadFile(CNTFIL,NDOM,IDOM,IPLACE,IFILE,NPLACES,
     +                    NFILES,NT0,NNATOMS,NSTEP,NDATT,IFORM,EMIN,
     +                    EMAX,INFILE,THETA,FI,RAR1,RAR2,PQFEX,
     +                    CUNDISP,CDISP,AID,EMK,EMK0,ESMK,VPI,VO,VV,
     +                    XISTMK,DELMK,DATTNO,XIST,DELWV,CUNDVB,CDVB)

 2223     CONTINUE     

 2222   CONTINUE

 2221 CONTINUE


C  Now check consistency of some of the data read from input files

      CALL CheckVal(NBTD,NBED,PNUM,NDOM,NPLACES,VI,VPI,V0RR,VV,VARST,
     +              NDATT,DATTNO)

cvb begin 
c    The following lines apply only in the energy dependent inner potential
c    was NOT included in the original LEED calculation. Otherwise, don't use
c    this part!!

C  do energy-dependent inner potential shift here
C  V0R = const. - E/quot.  'const' is VV from input, quot can be adjusted here.
C  EMK must not be corrected!

c      quot = 150.

c      DO 3335 IDATT = 1,DATTNO

c        ESMK(IDATT) = ESMK(IDATT) + ESMK(IDATT)/quot

c 3335 CONTINUE

cvb end

C  Open output file and write header (trivial but long if included here)

      OPEN(4,FILE = SEANAME,STATUS='UNKNOWN')

      call HeadDoc(WHICHR,WHICHG)


C  Modul 6: INITIALIZE SOME VALUES

C  IGEN is number of current generation

      IGEN=0

C   Throughout the main program, IPOP will always be counter for
C   current individual number; IPARAM will mean current parameter

      DO 1619 IPOP=1,MPS
      RPEIND(IPOP)=5.
 1619 RPEOLD(IPOP)=5.

      DO 7618 IPARAM=1,MNPRMK
 7618 WIDT(IPARAM)=1.
      
      AVERNEW=4.


C  determine parameters for first generation

C  If certain starting position is wanted (STAFLA = 1)
C  skip randomizing of parameters

      IF (STAFLA.eq.0) THEN

C  Determine parameters of new population

        CALL SEA_RCD(NDOM,NPS,NPRMK,NSTEP,PNUM,VARST,PARIND,RPEIND,
     +               WSK,WIDT,RMUT,NPAR)

      END IF

C  initialize PMOLD values prior to first generation so that restart runs can be 
C  successful

      DO IPOP=1,MPS

        PARSUM = 0
        DO IDOM=1,MNDOM,1
          PMOLD(IDOM,IPOP)= PARIND(NPRMK-NDOM+IDOM,IPOP) - 1
          PARSUM=PARSUM+PMOLD(IDOM,IPOP)
        ENDDO

        if (PARSUM.gt.0) then
          Norm = REAL(VARST(NPRMK)-1)/REAL(PARSUM)
        else
          Norm = 1.
        end if

        NDSL1 = 0
        DO IDOM = 1,NDOM-1
          PMOLD(IDOM,IPOP)=INT(REAL(PMOLD(IDOM,IPOP))*Norm)
          NDSL1=NDSL1+PMOLD(IDOM,IPOP)
        ENDDO

        PMOLD(NDOM,IPOP) = VARST(NPRMK)-1-NDSL1

      ENDDO

*************************************************************************

C START OPTIMIZATION LOOP

*************************************************************************

C START NEW GENERATION

 1492 CONTINUE

C  increment generation number

      IGEN=IGEN+1 

C  close and reopen control file for new generation 

      close(8,STATUS='keep')

      open (8,FILE='control.chem')

C  write current parameter values to control file

      write(8,*)

      write(8,'(a,i6,a)') 'Parameters of generation No.', IGEN, ':'

      do 2100 IPOP = 1, NPS

        IF (IGEN.eq.1) THEN
        write(8,'(500i3)') (PARIND(IPARAM, IPOP), IPARAM = 1, NPRMK)
        ELSE
        write(8,'(500i3)') (PARHELP(IPARAM, IPOP), IPARAM = 1, NPRMK)
        END IF

 2100 continue

C START FOR POPULATION LOOP - STAFLA=0 from now on

      DO 1848 IPOP=1,MPS

C  possibly restrict individual parameters in subroutine restrict

      CALL restrict(NPRMK,NPS,PARIND,IPOP)

C  compute TLEED intensities for different domains

      do 3000 IDOM = 1, NDOM

C  generate delta amp numbers in each file, IFNUM, and concentration step
C  numbers for each place, NPARC, in subroutine GetGrid, for later use in GetInt

        CALL GetGrid(NDOM,NPLACES,NFILES,NPRMK,NPRAS,NPS,IDOM,IPOP,
     .               NFIL,IFNUM,PARTYP,PARIND,VARST,NPARC,FILREL)


C  Compute intensities for current parameter values from delta amplitudes

        CALL GetInt(NDOM,IDOM,NPLACES,NFILES,NSTEP,NCONCS,NNATOMS,
     .              NDATT,NT0,NSURF,IFNUM,NFIL,CONC,CDISP,NPARC,
     .              DATTNO,EMK,VV,VO,VPI,THETA,FI,PQFEX,RAR1,RAR2,
     .              XISTMK,DELMK,ATSAS)

 3000 continue

C  average TLEED intensities incoherently according to weight of respective domains

      call Mischung(NDOM, NDATT, NT0, NPRMK, NPS, PARIND,
     .              PMISCH,DMISCH,ATSAS, ATSMK, DATTNO, IPOP,
     .              VARST)

C  This output possibility may be useful when in doubt about
C  computed r-factors.

C  Output current IFNUM

c      do IDOM = 1,NDOM
c        do IPLACE = 1,NPLACES(IDOM)
c          do IFILE = 1,NFIL(IDOM,IPLACE)

c            write(8,*) IDOM,IPLACE,IFILE,IFNUM(IDOM,IPLACE,IFILE)

c          enddo
c        enddo
c      enddo

C  Output some intensities to control file

c      do 1001 IESMK = 1, 5
c      write(8,'(F7.2,I6,4E13.6,/,5((5E13.6),/))') ESMK(IESMK),IPOP,
c     +              (ATSMK(IBEAM,IESMK),IBEAM=1,NT0)
c 1001 continue

C  end control output


C  Compute r-factor for current parameter set - produce RPEIND(IPOP)
C  Note RPEIND(IPOP) takes the value of the r-factor to be searched for, BARAV,
C  from subroutine Rfaktor

C  NETI is number of theoretical energies in subroutine Rfaktor

      NETI=DATTNO

C  renew dimension sizes just to make sure

      NDATA=MNDATA
      NDATT=MNDATT

      CALL RFAKTOR(ATSMK,NSS,NBTD,NBMD,NETI,PQ,KAV,SYM,ESMK,IPR,AE,
     +             NDATT,NET,AT,ET,EINCR,XPL,YPL,NDATA,ATP,YT,VI,NBED,
     +             TSTY2,V02,V01,VINCR,ARM,ARPEM,ERANGM,RAZZM,
     +             RANNM,RAVPM,XRPEM,NST1,NST2,MITTEL,IBP,R2,RPE,
     +             EE,NEE,EET,YE,WB,BRGES,BRINS,BRHAS,BV0,NPS,V0RR,
     +             IPOP,TST,TSE,TSE2,TSEY2,BENAME,NBE,WHICHG,
     +             WHICHR,RPEIND(IPOP),OVLG)

C  test statement
C      write(8,'("Rfaktor ",I5," ",I2," ok.")') IGEN,IPOP

C  This applies especially if half-order beams are considered and
C  the current structure has no superstructure but would be deadly in
C  T-T-comparison.

      IF (EOT.eq.0.) THEN

          IF (RPEIND(IPOP).EQ.0.0 ) THEN
            RPEIND(IPOP)=1.0
          ENDIF

      END IF


C  continue with next structure

 1848 CONTINUE

C  Possibly write parameter values to control file for test reasons

C        DO 3801 LPOP=1,MPS
C         WRITE(8,'(30I3,A3,F6.4)') (PARIND(IPARAM,LPOP),IPARAM=1,MNPRMK)
C     +                           ,'   ',BRGES(LPOP)
C         WRITE(8,'(30I3,A3,F6.4)') (PAROLD(IPARAM,LPOP),IPARAM=1,MNPRMK)
C     +                           ,'   ',BRGESOLD(LPOP)
C         WRITE(8,*)
C 3801   CONTINUE

C  compute width for PB distribution from current results

       CALL GetWid(WIDT,PARIND,PAROLD,RPEIND,RPEOLD,NormInd,
     +            NPS,NPRMK)

C  perform mutation?

      TCHANGE = .false.

      DO 2850 IPOP=1,MPS

        IF (RPEIND(IPOP).GE.RPEOLD(IPOP)) THEN

C  No improvement -> reset to old values

          RPEIND(IPOP)=RPEOLD(IPOP)
          BV0(IPOP)=BV0OLD(IPOP)
          BRGES(IPOP)=BRGESOLD(IPOP)
          BRINS(IPOP)=BRINSOLD(IPOP)
          BRHAS(IPOP)=BRHASOLD(IPOP)

          DO 2849 IPARAM=1,MNPRMK
            PARIND(IPARAM,IPOP)=PAROLD(IPARAM,IPOP)
 2849     CONTINUE

          DO 2848 IDOM = 1,NDOM
            PMISCH(IDOM,IPOP) = PMOLD(IDOM,IPOP)
 2848     CONTINUE

        ELSE

          TCHANGE = .true.

C  improvement -> store this generation

          RPEOLD(IPOP)=RPEIND(IPOP)
          BV0OLD(IPOP)=BV0(IPOP)
          BRGESOLD(IPOP)=BRGES(IPOP)
          BRINSOLD(IPOP)=BRINS(IPOP)
          BRHASOLD(IPOP)=BRHAS(IPOP)

          DO 2851 IPARAM=1,PNUM
            PAROLD(IPARAM,IPOP)=PARIND(IPARAM,IPOP)
 2851     CONTINUE

          DO 2852 IDOM = 1,NDOM
            PMOLD(IDOM,IPOP) = PMISCH(IDOM,IPOP)
 2852    CONTINUE

C  end storage of old values

        ENDIF

C  next population

 2850 CONTINUE

      IF (TCHANGE) then

C  compute current average R-factor AVERNEW

        AVERNEW = 0.

        DO 1782 IPOP=1,MPS
          AVERNEW=AVERNEW+RPEIND(IPOP)
 1782   CONTINUE

        AVERNEW=AVERNEW/FLOAT(MPS)

C order list of structures as function of (increasing) R-factor

        CALL ORDER(RPEIND,PARIND,BRGES,BRINS,BRHAS,BV0,
     +  RPEHELP,BRGESHELP,BRINSHELP,BRHASHELP,BV0HELP,PARHELP,
     +  PMISCH)

      END IF

C  write rfactors to file if improvement was made or output forced

      IF ( TCHANGE .OR. (MOD(IGEN,OUTINT).EQ.0) ) THEN

C  store this as last improved generation

       CALL OUTRF(PARHELP,NDOM,NPLACES,NFILES,NPRMK,NPRAS,NPS,NFIL,
     +  PARTYP,BRGESHELP,BRINSHELP,BRHASHELP,BV0HELP,AVERNEW,IGEN,
     +             WHICHG,WHICHR,PMISCH,DMISCH)

      ENDIF

C  Here ends a (nearly) endless loop - next generation

      if ((MAXGEN .eq. 0) .or. (IGEN .lt. MAXGEN)) THEN

C  Determine parameters of next population

        CALL SEA_RCD(NDOM,NPS,NPRMK,NSTEP,PNUM,VARST,PARIND,RPEIND,
     +               WSK,WIDT,RMUT,NPAR)

        goto 1492

      end if

      END

