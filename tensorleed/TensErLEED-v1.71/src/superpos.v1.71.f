C  Tensor LEED intensity calculation v1.7
C  VB, 13.04.00
C  for use with lib.superpos.f v1.7
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
***********************************************************************************
C
C  The program works but is not well maintained since it has mostly been superseded
C  by the search algorithm included with TensErLEED. Code does, in places, not keep
C  conventions of structured programming, and may thus seem complicated although
C  performing a simple task in principle. This has historical reasons; the program
C  was designed to run on extremely low-memory machines (mid-eighties style PCs), 
C  and its current "structure" allows just that. Please refer to the comments below 
C  to understand its full functionality.
C
C******************* file "superpos.f"
C  RHRINT V1.0  12/03/87  W.OED
C  slightly modificated by U.Loeffler 29.04.91:
C  CHECK-Criterion set to 1.e-9 instead of 1.e-10, since deltaamplitudes
C  calculated by the cyber are slightly different from those calculated
C  by the cray.
C
C  Modified by U.Loeffler 08/05/91:
C  Modified and updated version 1.1 :
C  The calculation of different geometries is documented now a bit.
C  There is a file called DOC, which contains a bijection of number
C  of geometry and corresponding variation.
C  There are several comments to document the order of geometries.
C  In this version the variation of layer distances is also possible.
C  Modified by U.Loeffler 05/06/91: Even more layer distance variations are
C                                   possible
C  Modified 04/95, V.Blum: 'Concentrations' for delta amplitudes from
C  different atoms now included in the input file
C
C-----------------------------------------------------------------------
C             M A I N
C-----------------------------------------------------------------------
C
      PROGRAM TRINT
C
C  PROGRAM RHRINT READS IN ALL DELTA AMPS NEEDED FOR ONE TEST-STRUCTURE
C  AND PRODUCES ALL INTENSITIES IN VAN-HOVE PUNCH OUTPUT FORMAT.
C
C  EXAMPLE FOR CONTROLFILE (UNIT 5).
C>  1  1               FILEANZAHL,VarAll: JEDER MIT JEDEM ? (1/0)
C>  1                  No of concentration steps
C> 1.0000              concentration step #1
C>DELWVRH  2  1  0     FILENAME1,VARIATIONEN,OBERFLAECHE?,FORMATIERT?
C>  1  2               NUMMERN DER AUSGEW. GEOMETRIEN
C
C   If the number of variations is set to a negative number, the variations
C   of that file are identical to that (neg) number.
C>   e.g. :
C>FILE1    2  1  0     FILENAME1,VAR,OBERFLAECHE ?,FORMATIERT?
C>  1  2               NUMMERN DER AUSGEW. GEOMETRIEN
C>FILE2    3  1  0     FILENAME1,VAR,OBERFLAECHE ?,FORMATIERT?
C>  1  2  3            NUMMERN DER AUSGEW. GEOMETRIEN
C>FILE3   -2  1  0     FILENAME1,VAR,OBERFLAECHE ?,FORMATIERT?
C>  0
C>FILE4    1  1  0     FILENAME1,VAR,OBERFLAECHE ?,FORMATIERT?
C>  7
C means that the variations of file3 are the same as for file2,
C i.e. 1  2  3.
C NOTE: The referenced filenumber must be less than the current number !
C       The last file must not be varied !
C
C If flag VarAll is set, all combinations of geometries will be calculated.
C If not, just the order of geometries will be calculated.
C In the case above, following combinations will be performed:
C 1117,2227,1337, ie the max. number of variations is the number of
C Geometries.
C
C  This version includes the possibility to mix deltaamplitudes from
C  different files. This allows for statistical variations within
C  the average t-matrix approximation without use of too large 
C  amplitude files. The block is included after the FILEANZAHL/VarAll
C  line, containing the number of concentration steps to be used and
C  then the weight factor to be multiplied with the deltaamplitude of
C  each file.
C  

      INCLUDE "PARAM"
      INCLUDE "GLOBAL"

C
C  INFILE,OUTFILE CONTAIN ALL FILENAMES
C
      CHARACTER*15 INFILE(MNFILES)
C
C  IFORM CONTAINS INFORMATION, WHETHER FILE IS WRITTEN FORMATTED OR NOT
C  IFVAR CONTAINS THE NUMBER OF VARIATIONS TO BE COMBINED
C  IFNUM CONTAINS ALL SELECTED VARIATIONS
C
      INTEGER IFORM(MNFILES),IFVAR(MNFILES),IFNUM(MNCSTEP+1,MNFILES)
C
C  NCSTEPS CONTAIN THE NUMBER OF VARIATIONS ON ACTUAL FILE
C  NATOMS  CONTAIN THE NUMBER OF ATOMS ON ACTUAL FILE
C  NSURF CONTAIN FLAGS.EQ.1 IF SURFACE-LAYER
C
      INTEGER NCSTEP(MNFILES),NATOMS(MNFILES),NSURF(MNFILES)
C
C  PQFEX CONTAINS THE NUMBERS OF OUTPUTBEAMS
C  RAR1,RAR2 CONTAIN THE RECIPROCAL LATTICE VECTORS
C  TRAR1,TRAR2 WORKINGSPACES
C
      DIMENSION PQFEX(2,MNT0),RAR1(2),RAR2(2),TRAR1(2),TRAR2(2)
C
C  CONC are the concentrations to be used for the various delta amp
C  files
C  
      DIMENSION CONC(MNCONCS,MNFILES)
C
C  CUNDISP,CDISP,AID CONTAIN ALL GEOMETRY-PARAMETERS
C
      DIMENSION CUNDISP(MNATOMS,3,MNFILES)
      DIMENSION CDISP(MNCSTEP,MNATOMS,3,MNFILES)
      DIMENSION AID(MNCSTEP,MNFILES)
C
C  DELWV CONTAINS ALL DELTA AMP'S
C
      COMPLEX DELWV(MNCSTEP,MNT0,MNFILES)
C
C  XIST CONTAINS ALL AMPLITUDES OF UNRECONSTRUCTED SURFACE
C
      COMPLEX XIST(MNT0)
C
      COMPLEX BKZ,AKZ(MNT0),DELACT,PRE
      DIMENSION AK2(MNT0),AK3(MNT0),APERP(MNT0),ATS(MNT0)
      INTEGER ICNT(MNFILES),VarAll,Finish
C
   10 FORMAT (2I4)
   11 FORMAT (1000I5)
   20 FORMAT (A15,3I3)
C------------------------------------------------------------------------------
C
C  NFILES NUMBER OF FILES IN ACTUAL VARIATION
C
      OPEN(1,FILE='INPUT')
      OPEN(2,FILE='OUTPUT')
      OPEN(3,FILE='DOC')
C     OPEN(5,FILE='CONTRIN')
C     OPEN(6,FILE='INTENS')
      VarAll = 1
C     label 1000 is start of a new run.
C     EOF of the DELWV-files causes the program to jump to this place.
C     If now EOF of SP_IN is reached, the program terminates.
 1000 READ (5,10,ERR=2001,END=2000) NFILES,VarAll
      READ (5,10,ERR=2001,END=2000) NCONCS

      IF (NFILES.gt.MNFILES) THEN
        write(3,*) "Dimension error: Number of files required",
     +             " exceeds MNFILES!"
        STOP
      END IF

      IF (NCONCS.gt.MNCONCS) THEN
        write(3,*) "Dimension error: Number of concentration ",
     +             " steps required exceeds MNCONCS!"
        STOP
      END IF

      DO 98 ICONC=1,NCONCS
      READ (5,'(10F7.4)',ERR=2001,END=2000) (CONC(ICONC,I),I=1,NFILES)
 98   CONTINUE
C-------------------------------------------------------------------------------
C  LOOP OVER INPUTFILENAMES AND SELECTED GEOMETRIES
C
      MaxNum = 0
      Max = 0
      DO 100 IFIL=1,NFILES
           READ (5,20) INFILE(IFIL),IFVAR(IFIL),NSURF(IFIL),
     1                  IFORM(IFIL)
C           write(3,*) "Complete"
           READ (5,11) (IFNUM(I,IFIL),I=1,IFVAR(IFIL))

           DO I = 1,IFVAR(IFIL)
             IF (IFNUM(I,IFIL).gt.MNCSTEP) THEN
        write(3,*) "Dimension error: Required Delta amp. no. ",
     +             I," in file no. ",IFIL," exceeds MNCSTEP!"
               STOP
             END IF
           ENDDO

           IF (IFVAR(IFIL).gt.Max) THEN
              Max = IFVAR(IFIL)
              MaxNum = IFIL
           ENDIF
  100 CONTINUE
C
      Lay=0
      DO 105 IFIL=1,NFILES
C        if IFVAR(IFIL) .lt. 0 => means layer variation (special treatment)
         IF (IFVAR(IFIL).lt.0) THEN
C           preset IFNUM(i,IFIL) with IFNUM of corresponding file
            DO 102 I=1,IFVAR(-IFVAR(IFIL))
               IFNUM(i,IFIL) = IFNUM(i,-IFVAR(IFIL))
  102       CONTINUE
C           set flag Lay
            Lay=1
         ENDIF
  105 CONTINUE
C
C     check whether last file is varied
      IF (IFVAR(NFILES).lt.0) THEN
         WRITE(2,*) ' Last file must not be varied !'
         STOP ' Last file must not be varied !'
      ENDIF
C
C------------------------------------------------------------------------------
C  NOW OPEN ALL FILES
C
      DO 110 IFIL=1,NFILES
           IF (IFORM(IFIL).EQ.0) THEN
              OPEN (10+IFIL,FILE=INFILE(IFIL),FORM='UNFORMATTED',
     1              STATUS='OLD')
           ELSE
              OPEN (10+IFIL,FILE=INFILE(IFIL),FORM='FORMATTED',
     1              STATUS='OLD')
           ENDIF
           REWIND (10+IFIL)
  110 CONTINUE
C
C  NOW READ IN HEADERS OF ALL FILES
C
C  LOOP over files
      DO 120 IFIL=1,NFILES
           IF (IFORM(IFIL).EQ.0) THEN
              READ (10+IFIL) TTHETA,TFI,TRAR1,TRAR2,INT0,
     1                    NATI,NCST
           NATOMS(IFIL) = NATI
           NCSTEP(IFIL) = NCST
           ELSE
              READ (10+IFIL,112) TTHETA,TFI,(TRAR1(I),I=1,2),
     1                        (TRAR2(I),I=1,2)
  112         FORMAT (6E13.7)
              READ (10+IFIL,114) INT0,NATOMS(IFIL),NCSTEP(IFIL)
  114         FORMAT (3I5)
           ENDIF
C
C  NOW STORE DATA OR CHECK, WETHER FILES ARE OK
C
           IF (INT0.gt.MNT0) THEN
        write(3,*) "Dimension error: Number of beams in file no. ",
     +             IFIL," exceeds MNT0!"
             STOP
           END IF

           IF (NCSTEP(IFIL).gt.MNCSTEP) THEN
        write(3,*) "Dimension error: Number of variations in file ",
     +             IFIL," exceeds MNCSTEP!"
             STOP
           END IF

           IF (NATOMS(IFIL).gt.MNATOMS) THEN
        write(3,*) "Dimension error: Number of atoms in file no. ",
     +             IFIL," exceeds MNATOMS!"
        STOP
           END IF

           IF (IFIL.EQ.1) THEN
              THETA = TTHETA
              FI    = TFI
              NT0   = INT0
              DO 116 I=1,2
                   RAR1(I) = TRAR1(I)
                   RAR2(I) = TRAR2(I)
  116         CONTINUE
           ELSE
              CHECK=ABS(THETA-TTHETA)
              CHECK=CHECK+ABS(FI-TFI)
              DO 118 I=1,2
                   CHECK=CHECK+ABS(RAR1(I)-TRAR1(I))
                   CHECK=CHECK+ABS(RAR2(I)-TRAR2(I))
  118         CONTINUE
              IF ((CHECK.GT.1E-9).OR.(INT0.NE.NT0)) THEN                290491
                 WRITE(2,*) ' ERROR : IMPROPER FILE COMBINATION '       290491
                 STOP 'IMPROPER FILE-COMBINATION!'
              ENDIF
           ENDIF
C
C  AND READ IN REMAINING PART OF FILEHEADER USING SUBROUTINE INDELT
C
           CALL INDELT(10+IFIL,IFORM(IFIL),NT0,NATOMS(IFIL),
     1                 NCSTEP(IFIL),PQFEX,CUNDISP(1,1,IFIL),
     2                 CDISP(1,1,1,IFIL),AID(1,IFIL))
C
  120 CONTINUE
C     End LOOP over files
C
C-----------------------------------------------------------------------------
C  NOW WRITE HEADER FOR OUTPUTFILE
C
      IF (NFILES.GT.20) GOTO 124
      WRITE (6,122) (INFILE(I),I=1,NFILES)
  122 FORMAT (20A15)
      GOTO 130
  124 WRITE (6,126) (INFILE(I),I=1,20)
  126 FORMAT (20A15,'and further DEL files')
  130 WRITE (6,134) NT0
  134 FORMAT (I5)
      DO 138 I=1,NT0
           WRITE (6,136) I,PQFEX(1,I),PQFEX(2,I),1
  136      FORMAT (I5,2F10.5,I5)
  138 CONTINUE
C      WRITE (6,122)
C
C-------------------------------------------------------------------------------
C  START LOOP OVER ENERGIES. THE PROGRAM RETURNS ALWAYS TO THIS LINE
C  UNTIL EOF IS REACHED ON ANY INPUTFILE (There is E=-1 set in INRINT)
C
C     set flag first energy (FE) to 1
      FE = 1.
  500 CONTINUE
C
C  READ IN ALL DELTA AMP'S BELONGING TO ONE ENERGY
C
      DO 510 IFIL=1,NFILES
           CALL INRINT(10+IFIL,IFORM(IFIL),TE,VV,VPI,NT0,
     1                 NCSTEP(IFIL),XIST,DELWV(1,1,IFIL))
C
C  IF ACTUAL ENERGY .LT. 0 EOF IS REACHED, START NEW RUN
C  EOF => INRINT sets E=-1. then jump to a completely new  run
           IF (TE.LT.0) GOTO 1000
C
           IF (IFIL.EQ.1) THEN
              E = TE
           ELSE
              IF (TE.NE.E) THEN
                 WRITE(2,*) ' ILLEGAL ENERGY READ ON INPUTFILE!',IFIL
                 STOP 'ILLEGAL ENERGY READ ON INPUTFILE!'
              ENDIF
           ENDIF
  510 CONTINUE
      EEV = (E-VV) * HARTREE
C
C  CALCULATE INCIDENT WAVEVECTOR
C
      AK  = SQRT(AMAX1(2.0*E-2.0*VV,0.))
      C   = AK*COS(THETA)
      BK2 = AK*SIN(THETA)*COS(FI)
      BK3 = AK*SIN(THETA)*SIN(FI)
C
C  AND PREPARE CALCULATION OF OVERLAYER PROPAGATOR
C
      BKZ = CMPLX(2.0*E-BK2*BK2-BK3*BK3,-2.0*VPI)
      BKZ = CSQRT(BKZ)
C
C  CALCULATE WAVEVECTORS OF OUTPUTBEAMS, PREPARE CALCULATION OF
C  OVERLAYER PROPAGATOR
C
      DO 520 N=1,NT0
           AK2(N) = BK2 + PQFEX(1,N)*RAR1(1) + PQFEX(2,N)*RAR2(1)
           AK3(N) = BK3 + PQFEX(1,N)*RAR1(2) + PQFEX(2,N)*RAR2(2)
           AK = 2.0*E - AK2(N)*AK2(N) - AK3(N)*AK3(N)
           AKZ(N) = CMPLX(AK,-2.0*VPI)
           AKZ(N) = CSQRT(AKZ(N))
           APERP(N) =  AK-2.0*VV
  520 CONTINUE
C
C-----------------------------------------------------------------------------
C     NOW PRESET VARIATION COUNTER AND START LOOP OVER CONCENTRATIONS
C
      IOUT = 1

      DO 1900 ICONC=1,NCONCS

C     NOW START LOOP OVER ALL FILES AND VARIATIONS
C     FIRST PRESET FILE COUNTERS
C
      DO 530 I=1,NFILES
           ICNT(I) = 1
  530 CONTINUE
C
C----------------------------------------------------------------------------
C     Here starts the loop over different geometries.
C     (The order of geometries is produced at the end of this loop !)
 1500 CONTINUE
C-----------------------------------------------------------------------------
C  CALCULATE ACTUAL AMPLITUDE
C
      CXDISP = 1000.
      DO 550 N=1,NT0
           DELACT = XIST(N)
           DO 540 IFIL=1,NFILES
C               ICNT(IFIL) is the actual number of variation in file IFIL
C               e.g. the 3rd number of variations containing information
C               which variation to use.
                IACT = IFNUM(ICNT(IFIL),IFIL)
C               IFNUM is the actual variation in IFIL
C               e.g. the 4th variation in DEL13...
                IF (IACT.EQ.0) GOTO 540
                CALL ADDWV(DELACT,DELWV(1,1,IFIL),NT0,NCSTEP(IFIL),
     1                     N,IACT,CONC(ICONC,IFIL))
                IF ((N.EQ.1).AND.(NSURF(IFIL).EQ.1)) THEN
                   CALL XMIN(CXDISP,CDISP(1,1,1,IFIL),NCSTEP(IFIL),
     1                       NATOMS(IFIL),IACT)
                ENDIF
  540      CONTINUE
C
C  AND CALCULATE NOW THE RESULTING INTENSITY
C
C          if (AK - 2VV) >= 0. (Beam emerges )
           IF (APERP(N).GE.0) THEN
C             somehow fussy, but uses CXDISP=XMIN, if surface;
C             CXDISP=0. otherwise.
              IF (CXDISP.EQ.1000.) CXDISP = 0.
              A = SQRT(APERP(N))
              PRE = (BKZ + AKZ(N)) * CXDISP
              PRE = CEXP(CMPLX(0.0,-1.0)*PRE)
              ATEST = CABS(DELACT)
              IF (ATEST.GT.1.E+10) THEN
C                preset  maximum intensity ! caution,if occures !
                 ATS(N) = 1.E+20
              ELSE
                 RPRE = CABS(PRE)
C                actual intesity
                 ATS(N) = ATEST * ATEST * RPRE * RPRE * A / C
              ENDIF
           ELSE
C             no intesity
              ATS(N) = 0.0
           ENDIF
  550 CONTINUE
C-----------------------------------------------------------------------------
C  AND WRITE RESULT TO OUTPUTFILE
C
      WRITE (6,552) EEV,FLOAT(IOUT)/10000.,(ATS(N),N=1,NT0)
  552 FORMAT (1F7.2,1F7.4,4E14.5,/,2000(5E14.5,/))
C
C     write documentation file (if first energy .eq. 1)
      IF (FE.eq.1.) THEN
         WRITE(3,560) IOUT
  560    FORMAT(' Geometry number ',I5,':')
         WRITE(3,561) ICONC
 561     FORMAT(' Concentration no. ',I5)
         DO 570 I=1,NFILES
            WRITE(3,565) I,IFNUM(ICNT(I),I)
  565       FORMAT(' File ',I5,'; Variation no : ',I5)
  570    CONTINUE

         WRITE(3,*) ' '
      ENDIF
C
      IOUT = IOUT+1
C
C------------------------------------------------------------------------------
C  INCREMENT ALL COUNTERS
C  This is the part mentioned above, which produces the order of geometries.
C
      new = 1
      Finish = 0
      DO 600 IFIL=1,NFILES
C        if IFVAR(IFIL) .lt. 0 => means layer variation (special treatment)
C        do not create new combination (new=0)
         IF ((IFVAR(IFIL).lt.0).or.(new.eq.0)) THEN
C           if really layer var :
            IF (IFVAR(IFIL).lt.0) THEN
C              use -IFVAR(IFIL) as corresponding filenumber
               ICNT(IFIL) = ICNT(-IFVAR(IFIL))
C              note : IFNUM(count,IFIL) has been preset above
            ENDIF
C              otherwise do nothing !
         ELSE
C           increment counter
            ICNT(IFIL) = ICNT(IFIL) + 1
            IF (ICNT(IFIL).GT.IFVAR(IFIL)) THEN
C              if VarAll=0, endcriterion is, whether max. number is reached
C              means jump to next concentration
               IF((IFIL.EQ.MaxNum).and.(VarAll.eq.0)) THEN
                  Finish = 1
               ENDIF
C              if counter reaches its limit, reset this counter,
C              use next file and increment that counter (see above)
               ICNT(IFIL) = 1
C              (jump to next concentration, if last variation is done and set FE=0.)
               IF ((IFIL.EQ.NFILES).and.(VarAll.eq.1)) THEN
                  GOTO 1900
               ENDIF
C              if VarAll=0, endcriterion is, whether max. number is reached
C              means jump to next concentration
               IF((IFIL.EQ.NFILES).and.(Finish.eq.1)
     1         .and.(VarAll.eq.0)) THEN
                  GOTO 1900
               ENDIF
            ELSE
C              normal case: vary all.
               IF (VarAll.eq.1) THEN
C                 normaly if no layer variations, just leave loop
                  IF (Lay.eq.0) GOTO 610
C                 if variations, check for layer variations, but do not create
C                 new combination
                  New = 0
               ENDIF
            ENDIF
         ENDIF
  600 CONTINUE
  610 CONTINUE
C
C---- jump to next geometry --------------------------------------------------
      GOTO 1500
C-------------------------------------------------------------------------------
C
C     Continue concentration loop

 1900 CONTINUE

C  jump to next energy, which is no longer first energy
   
      FE=0.
      GOTO 500
C
C-----------------------------------------------------------------------------
C
 2000 WRITE(2,*) '.  CORRECT TERMINATION'
      STOP '.  CORRECT TERMINATION'
 2001 WRITE(2,*) ' SOME ERROR OCCURED WHILE READING SP_IN '
      STOP ' SOME ERROR OCCURED WHILE READING SP_IN '
      END
C
