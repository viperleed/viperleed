C  VAN HOVE R-FACTOR PROGRAM  /  MODIFIED BY G. BESOLD (JAN. '82)
C                                ADDITION       (JAN.  '87)              260187
C  MODIFIED BY RD: AVERAGE-R-FACTORS EXCHANGED                           141290
C  Allocation error in subroutine RINTAV removed 
C  via introduction of a second field ATAV.       L.Hammer Oct. 2018
C  OUTPUT expanded for error curve calculations   L.Hammer Mar. 2021
C  INPUT formats expanded to allow for 9999 exp. and theo. beams  L.Hammer Mar. 2021
C
C##############################################################################
C  R-FACTOR PROGRAM FOR COMPARING LEED EXPERIMENT AND THEORY
C##############################################################################

C  Three different R-factors (the R2-factor, the Zanazzi-Jona R-factor and the
C  Pendry R-factor), together with their average can be produced. Experimental
C  IV-curves from different experiments can be averaged together (the same 
C  energy values, but no energy limits are required). The input energies may be
C  in any order.
C  It is possible to specify gaps in the experimental IV-curves. These energy
C  ranges will be ignored in the R-factor calculation. The gaps will also be
C  taken into account in the processing of the experimenta data.
C  The experiment can be smoothed before interpolation onto a regular energy
C  grid. Individual experimental beams can be skipped (see NBEA).
C  The theoretical IV-curves (read in with the format used in the van Hove-Tong
C  LEED programs) for related beams can be averaged together to account for
C  different beams. Individual theoretical beams can be skipped (see IBP).
C  The averaged theoretical beams can be labelled as integer order or fractio-
C  nal order beam. In this case the average for all integer order and/or all
C  fractional order beams is calculated and written to output (see MITTEL).
C  Individual theoretical geometries can also be skipped (see NSSK).

C******************************************************************************
C  Input formats (see description of individual quantities further down)
C******************************************************************************

C  - &NL1 EMIN=...,...,&END (some compilers expect a symbol different from &)
C  - KAV (400(25I4,/))
C  - IBP (400(25I4,/))
C  - MITTEL (400(25I4,/))
C  - $NL2 NSSK=...,...,WB=...,..., $END
C  - TEXT (describing experiment) (CHARACTER*80)
C  - NBEA (400(25I4,/))
C  - FORMAT (for exp energies and intensities) (19A4)
C  - BENAME(IB) (5A4)
C  - NEE(IB),FAC (1I4,1E13.4)
C  - EE(IB,IE),AE(IB,IE) (FORMAT)
C  - repeat last step until NEE(IB) pairs of values have been read in
C  - read BENAME(IB+1),NEE(IB+1),... until NBED beams have been read in
C  - TEXT (describing theory) (CHARACTER*80)
C  - NBTD (1I3)
C  - LAB(IB),PQ(IB),SYM(IB) (1I3,2F10.5,1I3) (LAB not used)
C  - repeat last step until NBTD beams have been read in
C  - read energy, geometry and theoretical intensities for all beams 
C    (van Hove-output format)
C  - repeat last step until NS geometries have been read in
C  - repeat last two steps until an END-OF-FILE is encountered
C  - if END-OF-FILE is encountered read next set of theoretical data (skipping
C    TEXT and beam information) in van Hoove-output format until LIMFIL sets of
C    theoretical data have been read in

C******************************************************************************
C  Limitations
C******************************************************************************

C  Energies to be input in eV, intensities in any unit (may be beam-dependent),
C  cf purpose of quantity FAC in subroutine READE.


      PROGRAM RFACTOR

C  include problem dependent parameter statements for dimensions

      INCLUDE "PARAM"

C  variables in namelist NL1

C  EMIN,EMAX:(in eV) energies at sides of desired theoretical energy range 
C          (intensities outside this interval will be ignored)
C  EINCR : (in eV) energy grid step, a value of 1. eV or lower is recommended
C          (a 10-fold denser grid is used for the Zanazzi-Jona integrals)
C  LIMFIL: no. of files from which consecutive theoretical data will be read
C  IPR   : controls amount of printout on unit 6
C          IPR = 0: output of overview of input parameters and of all R-factors
C          IPR = 1: additional output of raw input exp energies and intensities
C          IPR = 2: additional printout of exp and theor data after averaging,
C                   smoothing and interpolation as well as Pendry Y-function
C  VI    : imaginary part of inner potential used in Pendry Y-function
C  V0RR  : real part of inner potential as given in LEED programs
C  V01,V02:limits for inner potential variation (V0 = -V0RR+V01,...,-V0RR+V02
C          V01.le.V02, V01 and V02 should be multiples of EINCR)
C  VINCR : step of inner potential variation (must be.gt.0)
C  ISMOTH: no of times, exp IV-curves will be smoothed by 3-point smoothing
C  EOT   : decides whether experimental (READE) or pseudoexperimental ie theo-
C          retical (READT) spectra are compared with theoretical spectra
C  PLOT  : for PLOT = 1 plots of the experimental data and the theoretical
C          best-fit spectra are drawn
C  GAP   : GAP = 1 indicates gaps in the experimental spectra

      REAL EMIN, EMAX, EINCR, VI, V0RR, V01, V02
      INTEGER LIMFIL, IPR, ISMOTH, EOT, GAP

C  quantities that are not supplied in namelists or input formats for 
C  experimental and theoretical spectra

C  NBED  : number of experimental beams (including those to be averaged over)
C  NBTD  : number of theoretical beams (including those to be averaged over)
C  KAV   : beam groupings for averaging of theo beams, KAV = 1 for first group,
C          KAV = I for Ith group to be averaged, KAV = 0 to ignore beams
C          (KAV.le.NBTD must hold)
C  IBP   : relationship between theor and exp beams, IBP(I) = J indicates that 
C          the Ith exp beam (after averaging over different experiments) corre-
C          sponds to the Jth theor beam (counted after domain-averaging), when
C          J = 0, the (exp) beam I is skipped
C  MITTEL: labelling of theoretical beams (counted after averaging), MITTEL = 1
C          for integer beams, MITTEL = 2 for fractional beams
C  NINT,NFRAC:number of integer and fractional beams respectively
C  TEXT  : title that gives a short description of either exp or theor data

      INTEGER NBED, NBTD, KAV(MNBTD), IBP(MNBED), MITTEL(MNBTD)
      INTEGER NINT,NFRAC
      CHARACTER*80 TEXT

C  variables in namelist NL2

C  NSSK  : list of geometries to be skipped, eg NSSK = 4,7,8 skips 4th,7th and
C          8th geometry
C  WB    : weights for different beams (counting by numbers NBEA), used to
C          obtain R-factors averaged over beams (should not include weights for
C          differing energy ranges, since these are already programmed in)
C  WR    : WR(1), WR(2) and WR(3) are the weights of the R2-, the Zanazzi-Jona-
C          and the Pendry-R-factor respectively used in obtaining the average
C          R-factor

      INTEGER NSSK (MNS)
      REAL WB(MNBED), WR(3)

C  variables in namelist NL3

C  NORM  : for NORM = 1 the experimental and theoretical spectra are plotted 
C          for each beam with the same integral intensity, for NORM = 2 they
C          have the same maximum intensity
C  INTMAX: maximum intensity used for plotting the spectra
C  PLSIZE: determines the aspect ratio of the plots
C  XTICS : interval for which tics are given on the plots' energy axis (eV)

      INTEGER NORM, XTICS
      REAL INTMAX,PLSIZE(2)

C  variables for experimental input (subroutines READE/READT, GAPIN)

C  AE    : intensities of experimental beams
C  EE    : energies for which intensities of experimental beams are given
C  NEE   : number of energy grid points of experimental beams (before and after
C          averaging and interpolating)
C  NBEA  : averaging scheme for experimental beams, NBEA(I) = NBEA(J) means
C          beams I and J will be averaged together (for I.gt.J the relation
C          NBEA(I).gt.NBEA(J) must hold), NBEA(I) = 0 is used to skip beam I
C  BENAME: name of experimental beams (before and after averaging)
C  NGAP  : total number of gaps in the experimental spectra
C  GBM   : beam in which gap IGAP is found (gaps must be input in order of the
C          beams in which they occur and in order of energy in the same beam)
C  EG1,EG2:energy range of gap IGAP (before and after averaging)
C  GCOUNT: number of gaps in experimental beam IB (before and after averaging)

      REAL AE(MNBED,MNGP), EE(MNBED,MNGP), BENAME(5,MNBED)
      INTEGER NEE(MNBED), NBEA(MNBED)
      INTEGER NGAP, GBM(MNGAP), GCOUNT(MNBED)
      REAL EG1(MNGAP), EG2(MNGAP)

C  variables for averaging, smoothing and interpolating experimental data

C  NBE   : number of experimental beams after averaging
C  NG1,NG2:energy indices of the gap boundaries (after averaging)
C  NNN   : working spaces
C  XPL,YPL,WORX: working spaces
C  A,E   : working space for intensities and energies of current beam

      INTEGER NBE, NG1(MNGAP), NG2(MNGAP)
      INTEGER NNN(MNGP)
      REAL XPL(MNGP), YPL(MNGP), WORX(MNGP), A(MNGP), E(MNGP)

C  variables used for functions of experimental data

C  AP,APP: working spaces for first and second derivative of current beam
C  AEP   : first derivative of experimental beams
C  AEPP  : second derivative of experimental beams
C  YE    : Pendry Y-function of intensities for experimental beams
C  TSE   : integral over intensities for experimental beams
C  TSE2  : integral over squared intensities for experimental beams
C  TSEY2 : integral over squared Pendry Y-function for experimental beams

      REAL AP(MNGP), APP(MNGP), AEP(MNBED,MNGP), AEPP(MNBED,MNGP)
      REAL YE(MNBED,MNGP), TSE2(MNBED), TSE(MNBED), TSEY2(MNBED)

C  variables for reading, averaging and interpolating theoretical intensities

C  ES    : energy grid to be read in
C  ATS   : intensities of theoretical beams on energy grid ES (for all
C          geometries)
C  ATSAV : averaged intensities of theoretical beams on energy grid ES 
C          for all geometries)
C  NSS   : no of geometries remaining after skipping those specified in NSSK
C  NETI  : number of energy grid points at time of reading in
C  VALUES: code specifying the various geometries (ie interlayer spacings)
C  PQ    : lateral momentum of beams in units of reciprocal lattice vectors
C  PQAV  : lateral momentum of averaged beams in units of reciprocal lattice vectors
C  SYM   : beam symmetry codes ! not used; has to be scrapped off
C  NBT   : number of theoretical beams after averaging
C  AM    : maximal intensity in theoretical beams
C  AT    : intensities of theoretical beams for current geometry
C  ET    : energy grid for theoretical beams of current geometry

      REAL ES(MNGP), ATS(MNS,MNBTD,MNET), ATSAV(MNS,MNBTD,MNET)
      INTEGER NSS, NETI
      REAL PQ(2,MNBTD),PQAV(2,MNBTD),VALUES(MNS)
      INTEGER SYM(MNBTD), NBT
      REAL AM, AT(MNBTD,MNGP), ET(MNBTD,MNGP)

C  variables used for functions of experimental data

C  NET   : no of energy grid points of theoretical beams after interpolating
C  ATP   : first derivative of theoretical beams
C  ATPP  : second derivative of theoretical beams
C  YT    : Pendry Y-function of intensities for theoretical beams
C  TST   : integral over intensities for theoretical beams
C  TSTY2 : integral over squared Pendry Y-function for theoretical beams

      INTEGER NET(MNBTD)
      REAL ATP(MNBTD,MNGP), ATPP(MNBTD,MNGP)
      REAL YT(MNBTD,MNGP), TST(MNBTD), TSTY2(MNBTD)

C  loop variables (geometry loop, inner potential loop, beam loop)

C  IS    : current geometry index
C  D12   : current geometry code
C  NV0   : number of variations of the real part of the inner potential ! rigid shift of energy theo vs. exp
C  IV    : current inner potential real part variation index
C  V0    : current inner potential real part variation
C  V0R   : current inner potential real part
C  V0GR  : inner potential of best R-factor within a geometry
C  IBE   : current experimental beam (counted after averaging)
C  IBT   : current theoretical beam (counted after averaging)
C  KMIT  : current beam labelling (integer or fractional beam)

      INTEGER IS
      REAL D12
      INTEGER NV0, IV
      REAL V0, V0R, V0GR
      INTEGER IBE, IBT, KMIT

C  variables used to specify common energy range of exp and theor beams

C  NE1,NE2:energy indices of common energy range for experimental beam
C  NT1,NT2:energy indices of common energy range for theoretical beam
C  NEET  : length of common energy interval in times of EINCR
C  EET   : length of common energy interval for all beams in eV
C  NGT1,NGT2:energy indices of gap boundaries for theoretical beam (for gaps in
C          the associated experimental beam)

      INTEGER NE1, NE2, NT1, NT2, NEET
      REAL EET(MNBED)
      INTEGER NGT1(MNGAP), NGT2(MNGAP)

C  variables used for integrals over experimental and theoretical intensities

C  EPSE  : maximum of first derivative of experimental intensity
C  SS,SU : integral of experimental or theoretical intensity in the region 
C          below NE1 or NT1 (SS) or above NE2 or NT2 (SU)
C  SS2,SU2:same as SS and SU2 for integral over squared intensities
C  SSY2,SUY2:same as SS and SU2 for integral over squared Pendry Y-function
C  SG,SGY2:integral over theoretical intensitiy (SG) or theoretical Pendry 
C          Y-function (SGY2) in the range of a gap
C  SE,ST : integral over experimental or theoretical intensity in common energy
C          range
C  SE2,ST2:same as SE and ST for integral over squared intensities
C  SEY2,STY2:same as SE and ST for integral over squared Pendry Y-function
C  C     : normalization factor experiment/theory
C  SEM,STM:cumulative integral over experimental or theoretical intensities in
C          common energy range, differentiated for integer and fractional beams
C  S2,SG2: integral over (exp-C*theo)**2 in common energy range (S2) or gap
C          region (SG2)
C  SRZJ,SGZJ:integral for Zanazzi-Jona R-factor in common energy range (SRZJ)
C          or gap region (SGZJ)
C  SY2   : integral over (expY -C*theoY)**2 in common energy range

      REAL EPSE
      REAL SS, SU, SS2, SU2, SSY2, SUY2, SG, SGY2
      REAL SE, ST, SE2, ST2, SEY2, STY2
      REAL C, SEM(2), STM(2), S2, SG2, SRZJ, SGZJ, SY2

C  variables for R-factors

C  R2    : R2 R-factor for all beams
C  RRZJ  : Zanazzi-Jona R-factor for all beams
C  RPE   : Pendry R-factor for all beams
C  RAZ,RAN:dividend and divisor of Pendry R-factor for current geometry
C  RAVB  : weighted average over all R-factors for all beams
C  ERANG : total common energy range for current geometry
C  AR    : R-factors (R2, Zanazzi-Jona, Pendry) averaged over all beams
C  RAV   : average R-factor for current geometry
C  RAZM,RANM,ARM,ERANGM,RAVM:same as RAZ, RAN, AR, ERANG, RAV differentiated
C          for all integer and fractional beams
C  CIFE,CIFT:relation between fractional and integer beam intensities for
C          experimental and theoretical data, respectively
C  BRAV  : best average R-factor
C  BGRAV : best average R-factor of a geometry
C  BRAVB : weighted average over all R-factors for all beams of best average
C          R-factor geometry
C  BV0   : inner potential variation of best average R-factor geometry
C  BCIFE : relation of fractional to integer experimental beam intensities for
C          best average R-factor geometry 
C  BCIFT : same as BCIFE for theoretical beam intensities
C  NSB   : geometry index of best average R-factor geometry

      REAL R2(MNBED), RRZJ(MNBED), RPE(MNBED), RAVB(MNBED)
      REAL RAZ, RAN, RAZM(2), RANM(2)
      REAL ERANG, ERANGM(2), AR(3), ARM(3,2), RAV, RAVM(2)
      REAL CIFE, CIFT , BRAV, BGRAV, BRAVB(MNBED), BV0, BCIFE, BCIFT
      INTEGER NSB

C  variables used for plotting the spectra (subroutines COLUMN and SCRIPT)

C  NTMIN,NTMAX: indices of the lower and upper energy of the common energy
C          interval as they are stored in ET, ie ET(IBT,NTMIN(IBT)) denotes the
C          lowest energy for which exp and theo intensities are available
C  NTMINB,NTMAXB: NTMIN and NTMAX for the geometry with best average R-factor
C  NEMIN,NEMAX,NEMINB,NEMAXB: same as NTMIN,... for experimental storage arrays
C  BET,BAT: storage array of theoretical energies and intensities for the
C          best-fit geometry
C  STST  : integral over theoretical intensity for each beam
C  TSTB  : integral over theo intensity for each beam for the best-fit geometry
C  SEST,TSEB: same as STST and TSTB for experimental intensities

      INTEGER NTMIN(MNBTD), NTMAX(MNBTD), NTMINB(MNBTD), NTMAXB(MNBTD)
      INTEGER NEMIN(MNBED), NEMAX(MNBED), NEMINB(MNBED), NEMAXB(MNBED)
      REAL BAT(MNBTD,MNGP), BET(MNBTD,MNGP), STST(MNBTD), TSTB(MNBTD)
      REAL SEST(MNBED), TSEB(MNBED)

      INTEGER NG1B(MNGAP),NG2B(MNGAP),NGT1B(MNGAP),NGT2B(MNGAP)

      NAMELIST /NL1/ EMIN,EMAX,EINCR,LIMFIL,IPR,VI,V0RR,V01,V02,VINCR,
     *               ISMOTH,EOT,PLOT,GAP
      NAMELIST /NL2/ NSSK,WB,WR
      NAMELIST /NL3/ NORM,INTMAX,PLSIZE,XTICS

      NS   = MNS
      NBED = MNBED
      NBTD = MNBTD

C  default values in namelist NL1

      EMIN   = 20.
      EMAX   = 220.
      EINCR  = 1.
      LIMFIL = 1
      IPR    = 1
      VI     = 4.
      V01    = -5.
      V02    = 5.
      VINCR  = 2.
      ISMOTH = 2
      PLOT   = 0
      GAP    = 0

C  default values in namelist NL2

      DO IB = 1,MNBED
        WB(IB)  = 1.
        BENAME(1,IB) = 0.
        BENAME(2,IB) = 0.
      ENDDO
      DO IS = 1,MNS
        NSSK(IS) = 0
      ENDDO
      WR(1) = 0.                                                         201103
      WR(2) = 0.                                                         201103
      WR(3) = 1.                                                         201103

C  default values in namelist NL3

      NORM = 1
      INTMAX = 999.99
      PLSIZE(1) = 1.
      PLSIZE(2) = 1.
      XTICS = 50

      NINT =  0
      NFRAC = 0

      OPEN(7,FILE='ROUT')
      OPEN(9,FILE='ROUTSHORT')
      OPEN(8,FILE='WEXPEL',STATUS='OLD')
      REWIND(8)

C  read in any non-default values of the NL1 variables (as supplied by user)

      READ(8,NL1)
      WRITE(6,NL1)

      WRITE(6,3)NBTD,NBED                                                040481
 3    FORMAT(13H*NBTD,NBED = ,2I4)                                       040481

C  read theoretical beam averaging information

      READ(8,4)(KAV(J), J = 1,NBTD)
      WRITE(6,5)(KAV(J), J = 1,NBTD)
 4    FORMAT(400(25I4,/))
 5    FORMAT(32H*AVERAGING SCHEME OF THEOR BEAMS,/,400(25I4,/))

C  read experimental beam averaging information

      READ(8,4)(NBEA(I),I=1,NBED)
      WRITE(6,6)(NBEA(I),I=1,NBED)
 6    FORMAT(31H*AVERAGING SCHEME OF EXP. BEAMS,/,400(25I4,/))

C  read relationship of experimental to theoretical beams

      READ(8,4)(IBP(J), J = 1,NBED)
      WRITE(6,7)(IBP(J), J = 1,NBED)
 7    FORMAT(36H*RELATIONSHIP OF THEOR AND EXP BEAMS,/,400(25I4,/))

C  read beam labelling

      READ(8,4) (MITTEL(J), J = 1,NBTD)
      WRITE(6,8) (MITTEL(J), J = 1,NBTD)
 8    FORMAT(34H*LABELLING OF THEOR BEAM GROUPINGS,/,400(25I4,/))

      DO I = 1,NBTD
        IF (MITTEL(I).EQ.1) THEN
          NINT = NINT + 1
        ELSE IF (MITTEL(I).EQ.2) THEN
          NFRAC = NFRAC + 1
        ELSE IF (MITTEL(I).NE.0) THEN
          WRITE(6,9) I
          STOP
        ENDIF
      ENDDO

 9    FORMAT(8HMITTEL(,I4,26H) MUST BE EITHER 0, 1 OR 2,/,
     *       36H   *****   PLEASE CORRECT   *****   )

C  read list of geometries to be skipped and weights for the R-factor averages

      READ(8,NL2)
C      WRITE(6,NL2)

C  NSS will be number of geometries remaining after skipping

      NSS = 0
      DO IS = 1,NS
        IF (NSSK(IS).EQ.0) NSS = NSS + 1
      ENDDO

C  read data for the plots of the experimental and theoretical best-fit spectra

      READ(8,NL3)

C  read and print description of experiment

      READ(8,22)TEXT
      WRITE(6,23)TEXT
 22   FORMAT(A80)
 23   FORMAT(1H*,A80)

C  read data from experiment or pseudoexperiment

      IF (EOT.LE.0.5) THEN                                               RTHTH
        CALL READE(AE,EE,NBED,MNGP,NEE,BENAME,IPR)                       RTHTH
      ELSE                                                               RTHTH
        CALL READT(AE,EE,NBED,MNGP,NEE,BENAME,IPR)                       RTHTH
      END IF                                                             RTHTH

C  read information about gaps in the experimental spectra 

      DO IB = 1,NBED
        GCOUNT(IB) = 0
      ENDDO

      IF (GAP.EQ.1) THEN
        CALL GAPIN(NGAP,GBM,EG1,EG2,GCOUNT,NBED)
      ENDIF

C  average data from different experiments and order by increasing energy

      CALL EXPAV(AE,EE,NBED,MNGP,NEE,BENAME,NBEA,NBE,IPR,XPL,NNN,GCOUNT,
     *           NGAP,EG1,EG2,NG1,NG2)

      IGAP = 0
      DO IB = 1,NBE

        IF (GCOUNT(IB).EQ.0) THEN

          NE = NEE(IB)
          DO IE = 1,NE
            A(IE) = AE(IB,IE)
            E(IE) = EE(IB,IE)
          ENDDO

C  smooth experimental data

          IF (ISMOTH.GT.0) THEN
            CALL SMOOTH(ISMOTH,A,E,NE,IPR,IB)
          ENDIF

C  interpolate experimental data on working grid (multiples of EINCR eV)

          CALL INTPOL(A,NE,E,MNGP,EINCR,IPR,XPL,YPL,IB)

C  produce 1st and 2nd derivatives of experimental data

          IF (NE.GT.0) THEN
            CALL DER(A,NE,AP,EINCR)
            CALL DER(AP,NE,APP,EINCR)
          ENDIF

          NEE(IB) = NE
          DO IE = 1,NE
            AE(IB,IE)   = A(IE)
            EE(IB,IE)   = E(IE)
            AEP(IB,IE)  = AP(IE)
            AEPP(IB,IE) = APP(IE)
          ENDDO

C  for experimental beams with gaps a more complicated scheme is necessary

        ELSE

          IF (ISMOTH.GT.0) THEN
            NZ = 0

C  smooth sections between the gaps

            DO IGC = IGAP+1,IGAP+GCOUNT(IB)

              NE = NG1(IGC) - (NZ+1)
              DO IE = 1,NE
                A(IE) = AE(IB,IE+NZ)
                E(IE) = EE(IB,IE+NZ)
              ENDDO

              CALL SMOOTH(ISMOTH,A,E,NE,IPR,IB)

              DO IE = 1,NE
                AE(IB,IE+NZ) = A(IE)
                EE(IB,IE+NZ) = E(IE)
              ENDDO
              NZ = NG2(IGC)

            ENDDO
            NE = NEE(IB) - (NZ+1)
            DO IE = 1,NE
              A(IE) = AE(IB,IE+NZ)
              E(IE) = EE(IB,IE+NZ)
            ENDDO

            CALL SMOOTH(ISMOTH,A,E,NE,IPR,IB)

            DO IE = 1,NE
              AE(IB,IE+NZ) = A(IE)
              EE(IB,IE+NZ) = E(IE)
            ENDDO

          ENDIF

C  produce energy working grid

          NE = NEE(IB)
          DO IE = 1,NE
            A(IE) = AE(IB,IE)
            E(IE) = EE(IB,IE)
          ENDDO
          CALL GRID(A,NE,E,NGP,EINCR,WORX,NPTS,IMIN)
          NEE(IB) = NE

C  interpolate intensities of sections between gaps onto working grid

          NZ = 0
          NEWNZ = 0
          ITIL = 0
          ITIH = 0
          DO IGC = IGAP+1,IGAP+GCOUNT(IB)
            IF (NG1(IGC).GT.IMIN) THEN

              NE1 = NG1(IGC) - IMIN - NZ
              DO IE = 1,NE1
                XPL(IE) = EE(IB,IE+NZ)
                YPL(IE) = AE(IB,IE+NZ)
              ENDDO
              NZ = NG2(IGC)

              NG1(IGC) = 0
              NG2(IGC) = NPTS
              DO IPTS = 1,NPTS
                IF (WORX(IPTS).LE.EG1(IGC)) THEN
                  NG1(IGC) = IPTS
                ENDIF
                IF (WORX(NPTS-IPTS+1).GE.EG2(IGC)) THEN
                  NG2(IGC) = NPTS - IPTS + 1
                ENDIF
              ENDDO
              EG1(IGC) = WORX(NG1(IGC))
              EG2(IGC) = WORX(NG2(IGC))

              LEN = NG1(IGC) - NEWNZ - 1
              DO ILEN = 1,LEN
                XVAL = WORX(ILEN+NEWNZ)
                A(ILEN) = YVAL(XVAL,YPL,XPL,NE1,ITIL,ITIH)
              ENDDO

C  produce 1st and 2nd derivative of intensities

              IF (LEN.GT.0) THEN
                CALL DER(A,LEN,AP,EINCR)
                CALL DER(AP,LEN,APP,EINCR)
              ENDIF

              DO ILEN = 1,LEN
                AE(IB,ILEN+NEWNZ)   = A(ILEN)
                AEP(IB,ILEN+NEWNZ)  = AP(ILEN)
                AEPP(IB,ILEN+NEWNZ) = APP(ILEN)
              ENDDO

C  in the gap region the intensities and 1st and 2nd derivatives are set zero

              DO IPTS = NG1(IGC),NG2(IGC)
                AE(IB,IPTS)   = 0.
                AEP(IB,IPTS)  = 0.
                AEPP(IB,IPTS) = 0.
              ENDDO
              NEWNZ = NG2(IGC)

C  if no non-zero intensity occurs until gap IGC, set gap data to zero

            ELSE IF (IMIN.GE.NG2(IGC)) THEN

              NG1(IGC) = 0
              NG2(IGC) = 0
              EG1(IGC) = 0.
              EG2(IGC) = 0.

            ENDIF
          ENDDO

C  interpolate intensities onto working grid and calculate 1st and 2nd
C  derivative for the section after the last gap

          NE1 = NE - IMIN - NZ
          DO IE = 1,NE1
            XPL(IE) = EE(IB,IE+NZ)
            YPL(IE) = AE(IB,IE+NZ)
          ENDDO

          LEN = NPTS - NEWNZ
          DO ILEN = 1,LEN
            XVAL = WORX(ILEN+NEWNZ)
            A(ILEN) = YVAL(XVAL,YPL,XPL,NE1,ITIL,ITIH)
          ENDDO

          IF (LEN.GT.0) THEN
            CALL DER(A,LEN,AP,EINCR)
            CALL DER(AP,LEN,APP,EINCR)
          ENDIF

          DO ILEN = 1,LEN
            AE(IB,ILEN+NEWNZ)   = A(ILEN)
            AEP(IB,ILEN+NEWNZ)  = AP(ILEN)
            AEPP(IB,ILEN+NEWNZ) = APP(ILEN)
          ENDDO
          NEE(IB) = NPTS

        ENDIF
        IGAP = IGAP + GCOUNT(IB)

      ENDDO

C  produce Pendry Y-function for experimental data

      CALL YPEND(AE,AEP,1,NBED,MNGP,1,NBE,NEE,EE,YE,VI,IPR)

C  produce some integrals over experimental data

      DO IB = 1,NBE

        IE2 = NEE(IB)
        IF (IE2.GT.0) THEN

          TSE(IB)  = 0.
          TSE2(IB) = 0.
          TSEY2(IB)= 0.

          CALL VARSUM(AE,AE,AE,AE,1,1,NBED,1,MNGP,1,1,IB,1,1,IE2,0,
     *                EINCR,0.,0.,1,TSE(IB),YPL)

          IF (WR(1).GT.1.E-6) THEN
            CALL VARSUM(AE,AE,AE,AE,1,1,NBED,1,MNGP,1,1,IB,1,1,IE2,0,
     *                  EINCR,0.,0.,2,TSE2(IB),YPL)
          ENDIF 

          IF (WR(3).GT.1.E-6) THEN
            CALL VARSUM(YE,AE,AE,AE,1,1,NBED,1,MNGP,1,1,IB,1,1,IE2,0,
     *                  EINCR,0.,0.,2,TSEY2(IB),YPL)
          ENDIF

        ENDIF

      ENDDO

      CLOSE(8)

C  read and write title of theoretical data (supplied by LEED program)

      READ(5,22) TEXT
      WRITE(7,22) TEXT
      WRITE(6,23) TEXT

C  NBTD is number of beams to be read in (supplied by the LEED program)

      READ(5,26) NBTD
 26   FORMAT(I5) 

C  read in calculated energies and intensities (in range EMIN to EMAX)

      CALL RDPL(ES,MNGP,ATS,NSS,NBTD,MNET,NETI,VALUES,PQ,KAV,SYM,        281279
     *          EMIN,EMAX,LIMFIL,NSSK,NS,MNBTD)                          281279

C  perform domain averaging and check for too high theoretical intensities

      CALL RINTAV(ATS,NSS,NBTD,NETI,PQ,PQAV,KAV,SYM,NBT,ES,IPR,ATSAV)

      DO IB = 1,NBT                                                      121280

        CALL MAXINT(ATSAV,NSS,NBTD,IB,NETI,AM)                             121280

        IF (AM.GT.1.) THEN                                               121280
          WRITE(6,28)IB,AM                                              121280
 28       FORMAT(1H ,//,29H MAX INTENSITY IN THEOR BEAM ,1I3,            121280
     *           12H IS SUSPECT(,1E13.5,25H)- **** STOP PROGRAM ****)    121280
          STOP
        ENDIF

      ENDDO

      WRITE(7,29)
 29   FORMAT(1X,/,6X,'BEAM',12X,'IBE',3X,'D12',5X,'V0R',4X,'EMIN',4X,
     *       'EMAX',5X,'OVL',5X,'RAV')                                   201103
C  start loop over geometries

      BRAV = 1.E2
      BGRAV = 1.E2

      DO 300 IS = 1,NSS

      WRITE(6,30)IS,VALUES(IS)
 30   FORMAT(14H*STRUCTURE NO.,1I5,23H CHARACTERIZED BY VALUE,1F7.4)

      D12=VALUES(IS)

C  copy energies and intensities for geometry IS

      CALL STRIP2(ATSAV,NSS,NBTD,NETI,IS,AT,ES,ET,MNGP)

      DO IB = 1,NBT
        DO IE = 1,NETI
          A(IE) = AT(IB,IE)
          E(IE) = ET(IB,IE)
        ENDDO

C  interpolate theoretical data on working grid

        NE = NETI

        CALL INTPOL(A,NE,E,MNGP,EINCR,IPR,XPL,YPL,IB)

C  produce 1st and 2nd derivative of theoretical spectra

        IF (NE.GT.0) THEN
          CALL DER(A,NE,AP,EINCR)
          CALL DER(AP,NE,APP,EINCR)
        ENDIF

        NET(IB) = NE
        DO IE = 1,NE
          AT(IB,IE)   = A(IE)
          ET(IB,IE)   = E(IE)
          ATP(IB,IE)  = AP(IE)
          ATPP(IB,IE) = APP(IE)
        ENDDO
      ENDDO

C  produce Pendry Y-function of theoretical data

      CALL YPEND(AT,ATP,1,NBTD,MNGP,1,NBT,NET,ET,YT,VI,IPR)

C  produce some integrals over theoretical data

      DO IB = 1,NBT

        IE2=NET(IB)
        IF (IE2.GT.0) THEN

          CALL VARSUM(AT,AT,AT,AT,1,1,NBTD,1,MNGP,1,1,IB,1,1,IE2,0,
     *                EINCR,0.,0.,1,TST(IB),YPL)

          CALL VARSUM(YT,AT,AT,AT,1,1,NBTD,1,MNGP,1,1,IB,1,1,IE2,0,
     *                EINCR,0.,0.,2,TSTY2(IB),YPL)

        ENDIF
      ENDDO

C  start loop over inner potential values

      NV0 = INT((V02 - V01)/VINCR + 0.0001) + 1
      V0  = V01

      DO 200 IV = 1,NV0

      WRITE(6,60)V0
60    FORMAT(34H*THEOR. INNER POTENTIAL SHIFTED BY,1F7.2,3H EV)

      V0R = -V0RR + V0

C  initialize some variables

      DO IR = 1,3
        AR(IR) = 0.
        DO IMIT = 1,2                                                    260187
          ARM(IR,IMIT) = 0.                                              260187
          ARM(IR,IMIT) = 0.                                              260187
        ENDDO
      ENDDO
      RAZ = 0.
      RAN = 0.
      ERANG  = 0.
      DO IMIT = 1,2                                                      260187
        RAZM(IMIT) = 0.                                                  260187
        RANM(IMIT) = 0.                                                  260187
        ERANGM(IMIT) = 0.                                                260187
        RAVM(IMIT) = 0.                                                  260187
        SEM(IMIT) = 0.
        STM(IMIT) = 0.
      ENDDO

C  start loop over (experimental) beams (or theoretical beams)

      IGAP = 0

      DO 100 IBE = 1,NBE

      R2(IBE)   = 0.
      RRZJ(IBE) = 0.
      RPE(IBE)  = 0.
      RAVB(IBE) = 0.
      NEMIN(IBE) = 1
      NEMAX(IBE) = MNGP
      SEST(IBE) = 1.

C  ascertain corresponding theoretical beam and skip beams if necessary

      IBT = IBP(IBE)

      IF (IBT.EQ.0) THEN
        WRITE(6,66) IBE,(BENAME(I,IBE),I=1,5)
 66     FORMAT(12H*EXP BEAM NO,1I3,10H SKIPPED (,5A4,1H))
        GO TO 100
      ENDIF

      WRITE(6,70) IBE,(BENAME(I,IBE),I=1,5),IBT,(PQAV(I,IBT),I=1,2)
 70   FORMAT(9H*BEAM NO.,1I3,10H IN EXP. (,5A4,20H), WHICH IS BEAM NO.,
     *       1I3,12H IN THEORY (,2F6.3,1H))

C  another initialisation of variables

      KMIT = MITTEL(IBT)                                                 201103

C  determine energy interval of length EET common to experiment and theory,
C  bounded by the grid points (NE1,NE2) and (NT1,NT2) respectively

      CALL COMNEI(EE,NBED,NEE,ET,NBTD,MNGP,NET,IBE,IBT,V0,EINCR,
     *            NE1,NE2,NT1,NT2)

      IF (NE2.GT.NE1) THEN                                               040280
        NEET = NT1 - NE1
        EET(IBE) = EE(IBE,NE2) - EE(IBE,NE1)

C  determine grid points of gaps in experimental spectra for theoretical grid
C  and reduce EET by the length of gaps occuring in the common energy interval


        IF (GCOUNT(IBE).NE.0) THEN

          DO IGC = IGAP+1,IGAP+GCOUNT(IBE)

            IF (NG1(IGC).LE.NE1) THEN
              IF (NG2(IGC).LT.NE2) THEN
                NGT1(IGC) = 0
                NGT2(IGC) = 0
                NE1 = MAX(NE1,NG1(IGC))
                NT1 = NE1 + NEET
                EET(IBE) = EE(IBE,NE2) - EE(IBE,NE1)
              ELSE
                EET(IBE) = 0.
              ENDIF
            ELSE IF (NG1(IGC).LE.NE2) THEN
              IF (NG2(IGC).LE.NE2) THEN
                NGT1(IGC) = NG1(IGC) + NEET
                NGT2(IGC) = NG2(IGC) + NEET
                EET(IBE) = EET(IBE) - (NG2(IGC)-NG1(IGC)) * EINCR
              ELSE
                NGT1(IGC) = 0
                NGT2(IGC) = 0
                NE2 = MIN(NE2,NG2(IGC))
                NT2 = NE2 + NEET
                EET(IBE) = EE(IBE,NE2) - EE(IBE,NE1)
              ENDIF
            ELSE
              NGT1(IGC) = 0
              NGT2(IGC) = 0
            ENDIF

          ENDDO

        ENDIF

      ELSE

        EET(IBE) = 0.

      ENDIF

      NTMIN(IBT) = NT1
      NTMAX(IBT) = NT2
      NEMIN(IBE) = NE1
      NEMAX(IBE) = NE2

C  skip beam if no overlap between experiment and theory occurs

      IF (EET(IBE).EQ.0.) THEN
        WRITE(6,75) IBE
        WRITE(7,80) (BENAME(I,IBE),I=1,5),IBE,D12,V0R,EMIN,EMAX,0.,0.
 75     FORMAT(48H  **** NO OVERLAP THEOR./EXP. IN (EXP.) BEAM NO.,1I5)
 80     FORMAT(5A4,I5,1F8.4,3F8.2,1F10.2,1F8.4)
        GO TO 100
      ENDIF

C  find MAX(ABS(derivative of exp intensities)) for Zanazzi-Jona R-factor

      CALL EPSZJ(AEP,NBED,MNGP,IBE,NE1,NE2,EPSE)

C  if a IV-curve is truncated, above integrals should be reduced accordingly

      NE = NEE(IBE)
      NT = NET(IBT)

      SS   = 0.
      SS2  = 0.
      SSY2 = 0.
      SU   = 0.
      SU2  = 0.
      SUY2 = 0.

      IF (NE1.GT.1) THEN

        CALL VARSUM(AE,AE,AE,AE,1,1,NBED,1,MNGP,1,1,IBE,1,1,NE1,0,
     *              EINCR,0.,0.,1,SS,YPL)

        IF (WR(1).GT.1.E-6) THEN
          CALL VARSUM(AE,AE,AE,AE,1,1,NBED,1,MNGP,1,1,IBE,1,1,NE1,0,
     *                EINCR,0.,0.,2,SS2,YPL)
        ENDIF

        IF (WR(3).GT.1.E-6) THEN
        CALL VARSUM(YE,AE,AE,AE,1,1,NBED,1,MNGP,1,1,IBE,1,1,NE1,0,
     *              EINCR,0.,0.,2,SSY2,YPL)
        ENDIF

      ENDIF

      IF (NE2.LT.NE) THEN

        CALL VARSUM(AE,AE,AE,AE,1,1,NBED,1,MNGP,1,1,IBE,1,NE2,NE,0,
     *              EINCR,0.,0.,1,SU,YPL)

        IF (WR(1).GT.1.E-6) THEN
        CALL VARSUM(AE,AE,AE,AE,1,1,NBED,1,MNGP,1,1,IBE,1,NE2,NE,0,
     *              EINCR,0.,0.,2,SU2,YPL)
        ENDIF

        IF (WR(3).GT.1.E-6) THEN
          CALL VARSUM(YE,AE,AE,AE,1,1,NBED,1,MNGP,1,1,IBE,1,NE2,NE,0,
     *                EINCR,0.,0.,2,SUY2,YPL)
        ENDIF

      ENDIF

      SE   = TSE(IBE) - SS - SU
      SE2  = TSE2(IBE) - SS2 - SU2
      SEY2 = TSEY2(IBE) - SSY2 - SUY2
      SEST(IBE) = SE

      SS   = 0.
      SSY2 = 0.
      SG   = 0.
      SGY2 = 0.
      SU   = 0.
      SUY2 = 0.

      IF (NT1.NE.1) THEN

        CALL VARSUM(AT,AT,AT,AT,1,1,NBTD,1,MNGP,1,1,IBT,1,1,NT1,0,
     *              EINCR,0.,0.,1,SS,YPL)

        IF (WR(3).GT.1.E-6) THEN
          CALL VARSUM(YT,AT,AT,AT,1,1,NBTD,1,MNGP,1,1,IBT,1,1,NT1,0,
     *                EINCR,0.,0.,2,SSY2,YPL)
        ENDIF

      ENDIF

      IF (GCOUNT(IBE).NE.0) THEN

        DO IGC = IGAP+1,IGAP+GCOUNT(IBE)

          CALL VARSUM(AT,AT,AT,AT,1,1,NBTD,1,MNGP,1,1,IBT,1,NGT1(IGC),
     *                NGT2(IGC),0,EINCR,0.,0.,1,SG,YPL)

          IF (WR(3).GT.1.E-6) THEN
            CALL VARSUM(YT,AT,AT,AT,1,1,NBTD,1,MNGP,1,1,IBT,1,
     *                  NGT1(IGC),NGT2(IGC),0,EINCR,0.,0.,2,SGY2,YPL)
          ENDIF

        ENDDO

      ENDIF

      IF (NT2.NE.NT) THEN

        CALL VARSUM(AT,AT,AT,AT,1,1,NBTD,1,MNGP,1,1,IBT,1,NT2,NT,0,
     *              EINCR,0.,0.,1,SU,YPL)

        IF (WR(3).GT.1.E-6) THEN
          CALL VARSUM(YT,AT,AT,AT,1,1,NBTD,1,MNGP,1,1,IBT,1,NT2,NT,0,
     *                EINCR,0.,0.,2,SUY2,YPL)
        ENDIF

      ENDIF

      ST   = TST(IBT) - SS - SG - SU
      STY2 = TSTY2(IBT) - SSY2 - SGY2 - SUY2
      STST(IBT) = ST 

C  calculate normalization factor experiment/theory

      C = SE / ST

      WRITE(6,112) C
 112  FORMAT(28H NORMALIZATION C EXP/THEOR =,1E13.5)

      IF (KMIT.NE.0) THEN
        SEM(KMIT) = SEM(KMIT) + SE
        STM(KMIT) = STM(KMIT) + ST
      ENDIF

C  produce integrals involving both, experiment and theory

      SG2  = 0.
      SGZJ = 0.
      SGY2 = 0.

      IF (WR(1).GT.1.E-6) THEN

        CALL VARSUM(AE,AT,AE,AE,1,1,NBED,NBTD,MNGP,1,1,IBE,IBT,NE1,NE2,
     *              NEET,EINCR,0.,C,3,S2,YPL)

        IF (GCOUNT(IBE).NE.0) THEN

          DO IGC = IGAP+1,IGAP+GCOUNT(IBE)
            CALL VARSUM(AE,AT,AE,AE,1,1,NBED,NBTD,MNGP,1,1,IBE,IBT,
     *                  NG1(IGC),NG2(IGC),NEET,EINCR,0.,C,3,SG2,YPL)
          ENDDO

        ENDIF

      ENDIF

      IF (WR(2).GT.1.E-6) THEN

        CALL VARSUM(AEP,ATP,AEPP,ATPP,1,1,NBED,NBTD,MNGP,1,1,IBE,IBT,
     *              NE1,NE2,NEET,EINCR,EPSE,C,4,SRZJ,YPL)

        IF (GCOUNT(IBE).NE.0) THEN

          DO IGC = IGAP+1,IGAP+GCOUNT(IBE)
            CALL VARSUM(AEP,ATP,AEPP,ATPP,1,1,NBED,NBTD,MNGP,1,1,IBE,
     *                  IBT,NG1(IGC),NG2(IGC),NEET,EINCR,EPSE,C,4,SGZJ,
     *                  YPL)
          ENDDO

        ENDIF

      ENDIF

      IF (WR(3).GT.1.E-6) THEN

        CALL VARSUM(YE,YT,YE,YE,1,1,NBED,NBTD,MNGP,1,1,IBE,IBT,NE1,NE2,
     *              NEET,EINCR,0.,1.,3,SY2,YPL)

        IF (GCOUNT(IBE).NE.0) THEN

          DO IGC = IGAP+1,IGAP+GCOUNT(IBE)
            CALL VARSUM(YE,YT,YE,YE,1,1,NBED,NBTD,MNGP,1,1,IBE,IBT,
     *                  NG1(IGC),NG2(IGC),NEET,EINCR,0.,1.,3,SGY2,YPL)
          ENDDO

        ENDIF

      ENDIF

      S2   = S2 - SG2
      SRZJ = SRZJ - SGZJ
      SY2  = SY2 - SGY2

C  produce R-factors (all are normalized to about 1 for anticorrelated curves,
C  ie for (SIN(E))**2 compared with (COS(E))**2 over one period)

      EET(IBE) = WB(IBE) * EET(IBE)

C  R-factor based on integral of (exp - C*theo)**2

      IF (WR(1).GT.1.E-6) THEN

        R2(IBE) = 0.5 * S2 / SE2
        AR(1) = AR(1) + EET(IBE)*R2(IBE)
        IF (KMIT.NE.0) THEN                                              201103
          ARM(1,KMIT) = ARM(1,KMIT) + EET(IBE)*R2(IBE)                   201103
        ENDIF                                                            201103

      ENDIF

C  reduced R-factor according to Zanazzi-Jona (mult. by 0.5)

      IF (WR(2).GT.1.E-6) THEN

        RRZJ(IBE) = SRZJ / (0.027*SE)
        AR(2) = AR(2) + EET(IBE)*RRZJ(IBE)
        IF (KMIT.NE.0) THEN                                              201103
          ARM(2,KMIT) = ARM(2,KMIT) + EET(IBE)*RRZJ(IBE)                 201103
        ENDIF                                                            201103

      ENDIF

C  R-factor according to Pendry (mult. by 0.5)

      IF (WR(3).GT.1.E-6) THEN

        RPE(IBE) = SY2 / (SEY2 + STY2)
        RAZ = RAZ + WB(IBE) * SY2
        RAN = RAN + WB(IBE) * (SEY2 + STY2)
        IF (KMIT.NE.0) THEN                                              201103
          RAZM(KMIT) = RAZM(KMIT) + WB(IBE) * SY2                        260187
          RANM(KMIT) = RANM(KMIT) + WB(IBE) * (SEY2 + STY2)              260187
        ENDIF                                                            260187

      ENDIF

C  average over above R-factors for current beam

      WS=0.
      DO I=1,3
        WS=WS+WR(I)      
      ENDDO
      RAVB(IBE) = (WR(1)*R2(IBE) + WR(2)*RRZJ(IBE) + WR(3)*RPE(IBE))/WS
      ERANG = ERANG + EET(IBE)
      IF (KMIT.NE.0) THEN                                                201103
         ERANGM(KMIT) = ERANGM(KMIT) + EET(IBE)                          260187
      ENDIF                                                              260187

      EET(IBE) = EET(IBE) / WB(IBE)

C  write R-factor averages for current beam to output

      WRITE(6,115)                                                      040280
      WRITE(6,120) V0,EET(IBE),R2(IBE),RRZJ(IBE),RPE(IBE),RAVB(IBE)
      WRITE(7,80) (BENAME(I,IBE),I=1,5),IBE,D12,V0R,EMIN,EMAX,
     *               EET(IBE),RAVB(IBE)

115   FORMAT(45H   V0        EET     R2    RRZJ    RPE    RAV)
120   FORMAT(1H ,1F7.2,1F10.4,4F7.4)

      IGAP = IGAP + GCOUNT(IBE)

C  end of loop over beams

 100  CONTINUE

C  calculate average R-factor RAV for current geometry

      WRITE(6,140)
140   FORMAT(23H*AVERAGE OVER ALL BEAMS)

      DO I = 1,2
        AR(I) = AR(I)/ERANG
        DO IMIT = 1,2
          IF (ERANGM(IMIT).NE.0) THEN
            ARM(I,IMIT) = ARM(I,IMIT)/ERANGM(IMIT)                       260187
          ENDIF
        ENDDO
      ENDDO

      IF (WR(3).GT.1E-6) THEN
        AR(3) = RAZ/RAN
        DO IMIT = 1,2
          IF (RANM(IMIT).NE.0) THEN
            ARM(3,IMIT) = RAZM(IMIT)/RANM(IMIT)                          260187
          ENDIF
        ENDDO
      ENDIF

      RAV = (WR(1)*AR(1) + WR(2)*AR(2) + WR(3)*AR(3)) / WS
      DO IMIT = 1,2
        RAVM(IMIT) = (WR(1)*ARM(1,IMIT) + WR(2)*ARM(2,IMIT) +
     *                WR(3)*ARM(3,IMIT)) / WS
      ENDDO

C  calculate relation fractional beams/integer beams for experiment and theory

      IF (NINT.NE.0.AND.NFRAC.NE.0) THEN
         CIFE = (SEM(2) * FLOAT(NINT)) / (SEM(1) * FLOAT(NFRAC))
         CIFT = (STM(2) * FLOAT(NINT)) / (STM(1) * FLOAT(NFRAC))
      ENDIF

C  write R-factors for current geometry to output

      IF (NINT.NE.0) THEN
        WRITE(6,120) V0,ERANGM(1),(ARM(I,1),I=1,3),RAVM(1)
        WRITE(7,125) -1,D12,V0R,EMIN,EMAX,ERANGM(1),RAVM(1)              141290
      ENDIF

      IF (NFRAC.NE.0) THEN
        WRITE(6,120) V0,ERANGM(2),(ARM(I,2),I=1,3),RAVM(2)
        WRITE(7,126) -1,D12,V0R,EMIN,EMAX,ERANGM(2),RAVM(2)              141290
      ENDIF
  
      WRITE(6,120) V0,ERANG    ,(AR(I),I=1,3), RAV
      WRITE(7,127)  0,D12,V0R,EMIN,EMAX,ERANG,RAV                        141290

      IF (NFRAC.NE.0.AND.NINT.NE.0) THEN
        WRITE(6,128) CIFE,CIFT
      ENDIF

 125  FORMAT('AV.-INT             ',I5,F8.4,3F8.2,1F10.2,F8.4,'  <---')             260187
 126  FORMAT('AV.-FRAC ORD        ',I5,F8.4,3F8.2,1F10.2,F8.4,'  <---')             260187
 127  FORMAT('AVERAGE             ',I5,F8.4,3F8.2,1F10.2,F8.4,'  <---')
 128  FORMAT('RELATION FRAC TO INT BEAM INTENSITIES:',F6.3,
     *       ' FOR EXP DATA,',F6.3,' FOR THEOR DATA')
 129  FORMAT(' ')
 130  FORMAT(F8.4)


C  check if RAV is a better average R-factor than current best average BRAV

      IF (RAV.LT.BGRAV) THEN
        BGRAV = RAV
        V0GR = V0
      ENDIF

      IF (RAV.LT.BRAV) THEN
        BRAV = RAV
        NSB  = IS
        BV0  = V0

        IGAP = 0

        DO IBE=1,NBE
          IBT = IBP(IBE)

          IF (GCOUNT(IBE).NE.0) THEN
            DO IGC = IGAP+1,IGAP+GCOUNT(IBE)
              NG1B(IGC) = NG1(IGC)
              NG2B(IGC) = NG2(IGC)
              NGT1B(IGC) = NGT1(IGC)
              NGT2B(IGC) = NGT2(IGC)
            ENDDO
          ENDIF
          IGAP = IGAP + GCOUNT(IBE)

          NEMINB(IBE) = NEMIN(IBE)
          NEMAXB(IBE) = NEMAX(IBE)
          TSEB(IBE) = SEST(IBE)
          BRAVB(IBE) = RAVB(IBE)

          IF (IBT.NE.0) THEN 
             NTMINB(IBT) = NTMIN(IBT)
             NTMAXB(IBT) = NTMAX(IBT)
             DO IET=1,NET(IBT)
                BAT(IBT,IET) = AT(IBT,IET)
                BET(IBT,IET) = ET(IBT,IET)
             ENDDO
             TSTB(IBT) = STST(IBT)
          ENDIF
        ENDDO


        IF (NFRAC.NE.0.AND.NINT.NE.0) THEN
         BCIFE = CIFE
         BCIFT = CIFT
        ENDIF

      ENDIF

C  end of inner potential loop

 200  V0 = V0 + VINCR

      WRITE(7,129) 
      WRITE(9,130) BGRAV
      BGRAV = 1.E2

C  end of geometry loop

 300  CONTINUE

C  write output of experimental and best-fit theoretical data in column format

      IF (PLOT.EQ.1) THEN

        CALL COLUMN(NEMINB,NEMAXB,NTMINB,NTMAXB,EE,AE,BET,BAT,NBED,NBTD,
     *              MNGP,NBE,IBP,GCOUNT,NG1B,NG2B,NGT1B,NGT2B,MNGAP,
     *              TSEB,TSTB,BV0,INTMAX,NORM)

C  create gnuplot and latex scripts for plotting the spectra

        CALL SCRIPT(EMIN,EMAX,INTMAX,PLSIZE,XTICS,NBE,BENAME,BRAVB,NBED,
     *              TEXT,BRAV,IBP)

      ENDIF

C  write best average R-factor to output files

      WRITE(6,175) NSB,VALUES(NSB),BV0,BRAV
      WRITE(7,175) NSB,VALUES(NSB),BV0,BRAV
      IF (NFRAC.NE.0.AND.NINT.NE.0) THEN
        WRITE(7,128) BCIFE,BCIFT
      ENDIF

 175  FORMAT(30HBEST GEOMETRY ENCOUNTERED, NO.,1I5,
     *       25H, CHARACTERIZED BY VALUE ,1F7.4,/,
     *       27HWITH INNER POTENTIAL SHIFT ,1F7.2,
     *       24H EV, AVERAGE R-FACTOR = ,1F7.4)
 180  FORMAT(20H*CORRECT TERMINATION)

      WRITE(6,180)

      STOP

      END

