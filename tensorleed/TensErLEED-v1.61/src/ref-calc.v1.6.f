C  Tensor LEED reference calculation 
C  v1.2, VB 13.04.00
C  for use with lib.tleed.f v1.2
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
C  * This TLEED Tensor calculation program works for any multi-component 
C    surface LEED problem I have encountered so far. Yet, as the number 
C    of applications increases, new problems may turn up, requiring some 
C    additional coding. Please, in that process, avoid obsolete code that
C    is simply dragged along.
C
C  What the program does not do:
C
C  * The f77 implicit convention for integer (IJKLMNO) and real variables
C    is not kept. Each variable is declared explicitly, regardless of its
C    name (except for some counters, "I..." ). Subroutines mostly use the
C    convention, sometimes don't. Message: Check each variable's type 
C    before you use it.
C
      program leed

C  include parameter statements for dimensions etc

      INCLUDE "PARAM"

C  set unchangeable constants

      INCLUDE "GLOBAL"

C  MNDIM determines up to which LMAX variable dimensions for
C        Clebsch-Gordan coefficients are known to the program
C        must only be changed if data statements for NCA,NLMS
C        are to be enlarged 

      PARAMETER (MNDIM = 15)

C  set derived parameters

      PARAMETER (MNFAC = 4*MLMAX + 3)
      PARAMETER (MN = 2*MLMAX+1, MNN = MN*MN)
      PARAMETER (MNN1=MN,MNN2=MLMAX1,MNN3=MLMAX1)
      PARAMETER (MKLM = MN*(MN+1)/2)
      PARAMETER (MLEV = MLMAX1*(MLMAX1+1)/2)
      PARAMETER (MLEV2 = 2 * MLEV)
      PARAMETER (MLOD = MLMMAX-MLEV)

C  variable and array declarations in order of appearance

C  storage for variable dimensions (see PARAM)

      INTEGER IDEG,NL1,NL2,NL
      INTEGER KNBS,KNT,NPUN,NT0,NPSI,NEL,LMAX1
      INTEGER NSITE,NLTYPE,NSTACK,NNSUB,NBRAV

      REAL EMACH

C  variables for subroutine READIN():

C  title: descriptive title for calculation to be performed
C  TVA  : area of (1x1) unit mesh
C  TVB  : area of superlattice unit mesh
C  ARA1,ARA2: 2D real-space lattice vectors of (1x1) unit mesh
C  RAR1,RAR2: 2D reciprocal lattice vectors of (1x1) unit mesh
C  SS1,SS2,SS3,SS4: 4 possible registry shifts for substrate layer
C         calculation (in units of ARA1,ARA2)
C  ARB1,ARB2: 2D real-space lattice vectors of superlattice unit mesh
C  RBR1,RBR2: 2D reciprocal lattice vectors of superlattice unit mesh
C  SO1,SO2,SO3: 3 possible registry shifts for superlattice layers
C         calculation (in units of ARA1,ARA2!)
C  ASB  : convergence criterion - should be set .le. the smallest interlayer
C         distance for which layer doubling is used (enters into TST)
C  FR   : obsolete
C  ASE  : distance between first layer coordinate and onset of inner potential
C  V,VL,JJS : variables for subroutine SLIND
C  KNB  : number of beams in each beam set
C  SPQF : lateral momentum of beams in units of reciprocal lattice vectors
C  SPQ  : lateral momentum of beams in atomic units
C  KSYM : symmetry factor for beam group in beam list
C  TSTS : convergence criterion 
C  TST  : modified convergence criterion TSTS for layer doubling
C  NPU  : list of beams for which full dynamical output is desired 
C         normally identical with Tensor LEED output beam list
C  THETA,FI : (polar) incidence angles of incoming beam
C  VV :   muffin tin constant
C  VO :   obsolete - should be zero.
C  EPS,LITER : convergence criteria for substrate layer doubling sequence
C  TI,TF,DT,T0 : obsolete
C  LMAX : highest angular momentum to be used in calculation - should be
C         equal to MLMAX
C  L1    : LMAX+1
C  ES    : energies for which phase shifts are tabulated
C  PHSS  : l-dependent phase shifts for up to three elements
C  VPI,VPIS,VPIO : imaginary part of inner potential (damping)
C  EI,EF,DE: initial, final energy and step width for computation
  
      CHARACTER*80 TITLE
      REAL TVA,TVB
      REAL ARA1,ARA2,RAR1,RAR2
      DIMENSION ARA1(2),ARA2(2),RAR1(2),RAR2(2)
      REAL SS1,SS2,SS3,SS4
      DIMENSION SS1(2),SS2(2),SS3(2),SS4(2)
      REAL ARB1,ARB2,RBR1,RBR2
      DIMENSION ARB1(2),ARB2(2),RBR1(2),RBR2(2)
      REAL SO1,SO2,SO3
      DIMENSION SO1(2),SO2(2),SO3(2)
      REAL ASB
      DIMENSION ASB(3)
      REAL FR,ASE
      REAL V
      DIMENSION  V(MNL,2)
      COMPLEX VL
      DIMENSION VL(MNL,2)
      INTEGER JJS
      DIMENSION JJS(MNL,MIDEG)
      INTEGER KNB
      DIMENSION KNB(MKNBS)
      REAL SPQF,SPQ
      DIMENSION SPQF(2,MKNT),SPQ(2,MKNT)
      INTEGER KSYM
      DIMENSION KSYM(2,MKNT)
      REAL TST,TSTS
      INTEGER NPU
      DIMENSION NPU(MNPUN)
      REAL THETA,FI
      REAL VO,VV
      REAL EPS
      INTEGER LITER
      INTEGER LMAX,L1
      REAL ES,PHSS
      DIMENSION ES(MNPSI),PHSS(MNPSI,MNEL,MLMAX1)
      REAL VPI,VPIS,VPIO
      REAL EI,EF,DE

C  variables for subroutine READTL

C  IFORM : Formatted (1) or unformatted (0) Tensor LEED output
C  PQFEX : List of Tensor LEED beams (should be identical with list NPU)

      INTEGER IFORM
      REAL PQFEX
      DIMENSION PQFEX(2,MNT0)

C  variables for subroutine READGEO

C  CONC  : concentration of element IEL on site type ISITE
C  VIB   : vibrational amplitude of element IEL on site type ISITE
C  NSUB  : number of (Bravais) subplanes in layer type ILTYPE
C  LBRAV : layer type number ILTYPE is Bravais layer no. LBRAV (0 if composite layer)
C  LCOMP : layer type number ILTYPE is composite layer no. LCOMP (0 if Bravais layer)
C  LATT  : lattice type of layer type ILTYPE ( 1 = overlayer, 2 = substrate)
C  STYPE : site type of layer type ILTYPE, subplane ISUB
C  SUBPOS: 3D position of layer type ILTYPE, subplane ISUB
C  TSLAB : 0: bulk calculation required, 1: no bulk calc. required
C  ASA   : substrate interlayer vector (used in subas)
C  TOPLAYB: layer type of top "bulk" layer (dblgas)
C  BOTLAYB: layer type of bottom "bulk" layer (dblgas)
C  ASBULK: Interlayer vector between top and bottom bulk layer. Reason: Bulk may
C          be stacked of two different layers alternatingly. The vector between
C          those two layers is ASBULK. The bulk is then composed of the resulting
C          entity, repeatedly stacked using ASA (subras).
C  LTYPE : layer type of stacked layer no. ISTACK
C          note reverse order: ISTACK = 1 is top layer, ISTACK = NSTACK ist
C          lowest layer
C  LDIST : interlayer vector between layer ISTACK and lower layer ISTACK+1
C  TENS  : tensor output for layer ISTACK desired?

      REAL CONC,VIB
      DIMENSION CONC(MNSITE,MNEL),VIB(MNSITE,MNEL)
      INTEGER NSUB,LATT,STYPE,LBRAV,LCOMP
      DIMENSION NSUB(MNLTYPE),LATT(MNLTYPE),STYPE(MNLTYPE,MNSUB)
      DIMENSION LBRAV(MNLTYPE),LCOMP(MNLTYPE)
      REAL SUBPOS
      DIMENSION SUBPOS(MNLTYPE,MNSUB,3)
      INTEGER TSLAB,TOPLAYB,BOTLAYB
      REAL ASA
      DIMENSION ASA(3)
      REAL ASBULK
      DIMENSION ASBULK(3)
      INTEGER LTYPE
      DIMENSION LTYPE(MNSTACK)
      REAL LDIST
      DIMENSION LDIST(MNSTACK,3)
      INTEGER TENS
      DIMENSION TENS(MNSTACK)
      CHARACTER*40 LAYFILE
      DIMENSION LAYFILE(MNSTACK,MNSUB)

C  file channel numbers

C  OUTNO: file channel numbers for tensor output files LAYFILE(ISTACK,ISUB)

      INTEGER OUTNO
      DIMENSION OUTNO(MNSTACK,MNSUB)

C  general variables for computation

C  T0,T : temperature - no longer used (vibs set directly instead) but
C         needed as default input in TSCATF
C  AIDENT: identification of current spectrum in output file, useless here
C         since only one spectrum is computed
C  LMMAX: NO. OF SPHERICAL HARMONICS TO BE CONSIDERED
C  KLM  : dimension for lattice sums
C  LEV  : number of even (lm) components
C  LOD  : number of odd  (lm) components
C  LEV2 : 2*LEV
C  NCA  : lmax-dependent dimension of CAA (cf data statement later)
C  NCAA : dimension for composite layer Clebsh-Gordans, from CAAA()
C  CAA  : Clebsh-Gordon-coefficients for composite layer
C  N,NN,NLM: dimensions for computation of Clebsh-Gordons for matrices
C         X and TAU, CELMG()
C  NLMS : lmax-dependent dimension of X and TAU Clebsh-Gordans CLM
C  CLM  : Clebsh-Gordan coefficients for matrices X and TAU (RSMF,RTINV)
C  YLMC : special ylm for celmg
C  FAC2,FAC1: auxiliary arrays in CELMG
C  NFAC : variable dimension for array FAC
C  FAC  : array that contains values n! - formerly dragged through code
C         via a common block
C  LX,LXI,LT,LXM: permutations of the (L,M) sequence used to compute
C         intra-layer scattering more efficiently 
C  NN1,NN2,NN3: variable dimensions for gaunt coefficients PPP (cf. TSCATF)
C  PPP  : gaunt coefficients for temperature-dependent phaseshifts

      REAL T0,T
      REAL AIDENT
      INTEGER LMMAX,KLM,LEV,LOD,LEV2,NCAA
      INTEGER NCA
      DIMENSION NCA(MNDIM)
      REAL CAA
      DIMENSION CAA(MNLMO)
      INTEGER N,NN,NLM
      INTEGER NLMS
      DIMENSION NLMS(MNDIM)
      REAL CLM,YLMC,FAC2,FAC1
      DIMENSION CLM(MNLM),YLMC(MNN),FAC2(MNN),FAC1(MN)
      INTEGER NFAC
      DOUBLE PRECISION FAC
      DIMENSION FAC(MNFAC)
      INTEGER LX,LXI,LXM,LT
      DIMENSION LX(MLMMAX),LXI(MLMMAX),LXM(MLMMAX),LT(MLMMAX)
      INTEGER NN1,NN2,NN3
      REAL PPP
      DIMENSION PPP(MNN1,MNN2,MNN3)

C  variables used in energy loop

C  E    : current energy
C  EEV  : vacuum energy in eV (for output purposes)
C  DCUTS,DCUTO: limiting radii for lattice sums in FMAT - 
C               DCUTS, DCUTO identical here.
C  KSQ : total squared momentum of beams outside the crystal
C  AK   : lateral momentum of incident beam
C  AK2,AK3: x and y momentum of incident beam
C  PSQ  : lateral momentum of Tensor LEED beams relative to
C         incident beam (0,0)
C  AK2M,AK3M: (negative) absolute lateral momentum of Tensor LEED beams
C         (for use as incident beams in time-reversed LEED calculation)

      REAL E,EEV
      REAL DCUTS,DCUTO
      REAL KSQ,AK,AK2,AK3
      REAL PSQ,AK2M,AK3M
      DIMENSION PSQ(2,MNT0),AK2M(MNT0),AK3M(MNT0)

C  variables for atomic t-matrix calculation

C  DR0  : zero point vibrational amplitude - must be zero here
C  TMAT : atomic t-matrix for each site ISITE
C  VSITE: site-dependent inner potential shift - set 0 here !
C         if ever used, please set in input!
C  DRPER,DRPAR: working vibrational amplitudes for TSCATF
C  AF,CAF : zero-temperature and temperature-corrected atomic
C         t-matrices for current element IEL in tscatf; used differently
C         in RSMF
C  PHS,DEL,CTAB,SUM,BJ: working spaces for TSCATF and PSTEMP
C  TMAT : storage array for site-specific atomic t-matrices including ATA

      REAL DR0
      REAL VSITE
      COMPLEX AF,CAF
      DIMENSION AF(MLMAX1),CAF(MLMAX1)
      REAL PHS
      DIMENSION PHS(MLMAX1)
      COMPLEX DEL,CTAB,SUM
      DIMENSION DEL(MLMAX1),CTAB(MNN3),SUM(MNN2)
      COMPLEX*16 BJ
      DIMENSION BJ(MNN1)
      COMPLEX TMAT
      DIMENSION TMAT(MNSITE,MLMAX1)

C  variables for NEXIT loop

C  subroutine decide

C  NEXIT: counter for all beams that a Tensor LEED calculation is
C         required for: the incident beam and all (reversed) Tensor
C         LEED output beams PQFEX. 
C  PSQ1,PSQ2: lateral momentum components of current beam NEXIT relative
C         to last calculated "incident beam"
C  AK21,AK31: momentum components the of incident beam for which
C         layer diffraction matrices were last calculated
C  EMERGE: set to 1 if current beam is non-evanescent outside crystal
C  KNOWN: set to 1 if diffraction properties up to second layer are
C         already known from previous computation NEXIT (this one concerns
C         only the case of off-normal incidence or the - currently 
C         unimplemented - construction of a superstructure not present in
C         the reference structure via Tensor LEED).

      INTEGER NEXIT      
      REAL PSQ1,PSQ2
      REAL AK21,AK31
      INTEGER EMERGE,KNOWN

C  subroutine beams

C  NB   : beam set size used in computation
C  PQ   : lateral momentum of actually used beams in atomic units
C  PQF  : lateral momentum of actually used beams in units of substrate reciprocal
C         lattice vectors
C  SYM  : symmetry weight for beam group in beam list (in Tensor LEED, this kind
C         of symmetry is currently not exploited - do not use any symmetry code but 1!)
C  NPUC : list of beams that need to be punched out at current energy
C  MPU  : no. of beams that need to be punched out at current E
C  NT   : total number of beams included in calculation at current energy
C  NP   : maximum number of beams in a single beamset at current E (used as a dimension later)

      INTEGER NB
      DIMENSION NB(MKNBS)
      REAL PQ,PQF,SYM
      DIMENSION PQ(2,MKNT),PQF(2,MKNT),SYM(2,MKNT)
      INTEGER NPUC
      DIMENSION NPUC(MNPUN)
      INTEGER MPU,NT,NP

C  subroutine fmat

C  NLS  : number of sublattices to be treated in FMAT (either 1 or NL, use NL!)
C  FLMS : output lattice sums F_s(l,m) from fmat 
C  SCC,SA: working spaces for FMAT - see there!

      INTEGER NLS
      COMPLEX FLMS,SCC,SA
      DIMENSION FLMS(MNL,MKLM),SCC(MIDEG,4*MLMAX+1),SA(2*MLMAX*MIDEG)

C  subroutine rsmf (and rtinv)

C  IT   : leave set to zero
C  ID   : compute scattering matrices for registries 1, ..., ID - currently
C         not enabled. Might be useful for substrate type layer in subref, hcp stacking
C  LAY  : substrate (1) or overlayer (2) type layer - LAY = LATT(ILTYPE)
C  NLL  : number of sublattices for rsmf (NL for substrate, 1 for overlayer)
C  NA   : offset in beam list PQ for subsequent rsmf (to treat beam sets separately)
C  NS   : offset in scattering matrices for subsequent rsmf (to treat beam sets separately)
C  NM   : internal layer diffraction matrix dimension for rsmf
C  NAA  : number of beams in subsequent rsmf calculation
C  ROP,TOP: reflection and transmission matrices for bravais layer by rsmf
C  AMULT,CYLM: complex working spaces 
C  FLM  : total lattice sum, added from individual FLMS
C  XEV,XOD: meant to be even and odd (lm)-components of X-matrix - note XEV has confusing
C         dimensions since also used in rtinv in a different fashion
C  YLM  : used as Y-function & aux. array
C  YLME,YLMO: working spaces (angular eigenfunctions etc)
C  XEVST,XODST: intermediate storage for X-matrix from rsmf
C  IPLE,IPLO: also "working spaces" 
C  REPSTO,TRAPSTO,REMSTO,TRAMSTO: storage arrays for layer diffraction matrices
C         (reflection, incidence from +x, transmission, inc. +x, 
C         reflection, inc -x, transmission,inc. -x)
C  XEVSTO,XODSTO: permanent storage array for msmf x-matrices for each layer
C         unlike TENSOV, these must not be used as pointers.

      INTEGER IT
      INTEGER ID,LAY,NLL,NA,NS,NM,NAA
      COMPLEX ROP,TOP
      DIMENSION ROP(MKNT,MKNT),TOP(MKNT,MKNT)
      COMPLEX AMULT,CYLM,FLM
      DIMENSION AMULT(MKNT),CYLM(MKNT,MLMMAX),FLM(MKLM)
      COMPLEX XEV,XOD
      DIMENSION XEV(MLEV,MLEV2),XOD(MLOD,MLOD)
      COMPLEX YLM,YLME,YLMO
      DIMENSION YLM(MLMMAX),YLME(MLEV),YLMO(MLOD)
      COMPLEX XEVST,XODST
      DIMENSION XEVST(MLEV,MLEV),XODST(MLOD,MLOD)
      INTEGER IPLE,IPLO
      DIMENSION IPLE(MLEV),IPLO(MLOD)
      COMPLEX REPSTO,TRAPSTO,REMSTO,TRAMSTO
      DIMENSION REPSTO(MNLTYPE,MKNT,MKNT),TRAPSTO(MNLTYPE,MKNT,MKNT)
      DIMENSION REMSTO(MNLTYPE,MKNT,MKNT),TRAMSTO(MNLTYPE,MKNT,MKNT)
      COMPLEX XEVSTO,XODSTO
      DIMENSION XEVSTO(MLEV,MLEV,MNBRAV2),XODSTO(MLOD,MLOD,MNBRAV2)

C  subroutine rtinv

C  NOPT : 1 should be default (calc matrix for all beams) - 2 would assume a symmetric
C         layer, 3 just 1 inc beam. 2 and 3 can not be used in Tensor LEED calculation.
C  NEW  : would save computing time for tau-matrices, but can be a dangerous bug unless
C         set to 1
C  NLAY : dimension - number of sublayers in current rtinv
C  NLAY2,LMG,L2M: variable dimensions for rtinv
C  NTAU : number of different site types (t-matrices) in next rtinv call
C  TAUNO: ID of site type ISITE in current rtinv (t-matrices are passed via LPS and TSF,
C         contracted versions of STYPE and TMAT, to avoid unnecessary calls to TAUMAT)
C  TSF  : temperature-corrected atomic t-matrices for use in RTINV 
C         - contains only t-matrices actually used in next rtinv
C  LPS  : contains ID of t-matrix TSF(LPS(ISUB),...) belonging to subplane ISUB
C  LPSS : later storage array for reordered sublayer chemistry (subroutine srtlay)
C  NORD : map between input sublayer number NORD(ILAY) and sublayer number ILAY 
C         after srtlay reordering
C  POS  : internal array for RTINV that contains positions of sublayers of current layer
C  POSS : same as POS, but ordered by rising z (sim. to LPSS)
C  LMT  : variable dimension for tau-matrix in rtinv (should be NTAU*LMMAX)
C  TV1  : unit cell area of composite layer
C  ROM,TOM: reflection and transmission matrix for composite layer, but incidence from
C           -x (cf. ROP TOP, declared for rsmf)
C  TAU  : contains tau-matrices for the sublayers of rtinv
C  TAUG,TAUGM: derivates of tau matrix used during rtinv
C  MGH  : contains identification of sublayer propagator positions in matrix GH
C  DRL,SDRL: inter-subplanar vectors after reordering by srtlay, and (superfluous) storage
C  TEST : cutoff for rec. lattice sum in rtinv 
C  GH   : matrix containing inter-subplanar propagators for composite layer
C  RG   : inter-subplanar plane-wave phase shifts
C  TS,TG,VT,TH: working spaces for rtinv
C  IPL  : working space for matrix inversions 
C  XH   : working space for subroutine CXMTXT
C  HGHD : working space for subroutine GHD
C  RG1,RG2: working spaces for subroutine MFOLT
C  TENSOV: storage for intra-layer wavefield components of composite layer type
C         LCOMP(ILTYPE). In subroutines, Tensov is used as a pointer to the current
C         Tensor location. While not conforming to structured programming, this 
C         method saves a considerable amount of time.

      INTEGER NOPT,NEW
      INTEGER NLAY,NLAY2,LMG,NTAU,L2M
      INTEGER TAUNO
      DIMENSION TAUNO(MNSITE)
      COMPLEX TSF
      DIMENSION TSF(MNSUB,MLMAX1)
      INTEGER LPS,LPSS,NORD
      DIMENSION LPS(MNSUB),LPSS(MNSUB),NORD(MNSUB)
      DIMENSION POS(MNSUB,3),POSS(MNSUB,3)
      INTEGER LMT
      REAL TV1
      COMPLEX ROM,TOM
      DIMENSION ROM(MKNT,MKNT),TOM(MKNT,MKNT)
      COMPLEX TAU,TAUG,TAUGM
      DIMENSION TAU(MLMT,MLEV),TAUG(MLMT),TAUGM(MLMT)
      INTEGER MGH
      DIMENSION MGH(MNSUB,MNSUB)
      REAL DRL,SDRL
      DIMENSION DRL(MNSUB2,3),SDRL(MNSUB2,3)
      INTEGER NUGH,NGEQ,NGOL
      DIMENSION NUGH(MNSUB2),NGEQ(MNSUB2),NGOL(MNSUB2)
      REAL TEST
      DIMENSION TEST(MNSUB2)
      COMPLEX GH,RG,TS,TG,VT,TH
      DIMENSION GH(MLMG,MLMMAX),RG(3,MNSUB,MKNT),TS(MLMN),TG(2,MLM2N)
      DIMENSION VT(MLM2N),TH(MLMNI,MLMNI)
      INTEGER IPL
      DIMENSION IPL(MLMNI)
      COMPLEX XH,HGHD,RG1,RG2
      DIMENSION XH(MLEV),HGHD(MN),RG1(MNSUB),RG2(MNSUB)
      COMPLEX TENSOV
      DIMENSION TENSOV(2,MLMN,MKNT,MNCOMP)

C  dimensions for layer stacking

C  subroutines dblgas & subras (for bulk)

C  RAMP,TAPP,RAPM,TAMM,RBMP,TBPP,RBPM,TBMM:
C  reflection and transmission matrices for bulk layers (RA later used
C  for DBG calls too) - nomenclature: R/T = reflection/transmission,
C  A/B = top/bottom bulk layer, P/M = exit to +/- x, P/M: incidence
C  from +/-x
C  S1,S2,S3,S4,PP,XS,JNT: used as working spaces (intermediate storage)
C       in subroutines DBLGAS,SUBRAS,DBGT and ADREF2T
C  REBELO: reflection matrix of crystal below current layer
C  REABOV: reflection matrix after current layer stacking step
C  ASV  : current interlayer vector for DBGT and ADREF2T
C  EFSCAM,EFSCAP: storage arrays for effective scattering matrices of each
C         stacked layer ISTACK - to be used later when retrieving the wave
C         field inside the crystal for various incident beams   
C         Used as pointer to the correct array in subroutines (saves time & space)
C  WK   : working space for ADREF2T
C  XI   : outgoing amplitudes from crystal, for current incident beam
C  APLUS,AMINUS,ADOWN: for current incident beam, total wavefield
C         impinging on current layer (towards +x), current layer (towards -x),
C         and impinging downwards, onto next layer (towards +x)

      COMPLEX RAMP,TAPP,RAPM,TAMM,RBPM,TBMM,RBMP,TBPP
      DIMENSION RAMP(MKNT,MKNT),TAPP(MKNT,MKNT)
      DIMENSION RAPM(MKNT,MKNT),TAMM(MKNT,MKNT)
      DIMENSION RBMP(MKNT,MKNT),TBPP(MKNT,MKNT)
      DIMENSION RBPM(MKNT,MKNT),TBMM(MKNT,MKNT)
      COMPLEX S1,S2,S3,S4
      DIMENSION S1(MKNT,MKNT),S2(MKNT,MKNT),S3(MKNT,MKNT),S4(MKNT,MKNT)
      COMPLEX PP,XS
      DIMENSION PP(2,MKNT),XS(MKNT)
      INTEGER JNT
      DIMENSION JNT(MKNT)
      COMPLEX REBELO,REABOV
      DIMENSION REBELO(MKNT,MKNT),REABOV(MKNT,MKNT)
      REAL ASV
      DIMENSION ASV(3)
      COMPLEX EFSCAM,EFSCAP
      DIMENSION EFSCAM(MKNT,MKNT,MNSTACK),EFSCAP(MKNT,MKNT,MNSTACK)
      COMPLEX WK
      DIMENSION WK(3,MKNT)
      COMPLEX XI,APLUS,AMINUS,ADOWN
      DIMENSION XI(MKNT),APLUS(MKNT),AMINUS(MKNT),ADOWN(MKNT)

C  full dyn. output: rint and outxist

C  AT   : outgoing intensities
C  ATP  : intensities of requested output beams only, written to file no. 7 (fd.out)
C  XIST : storage array for outgoing amplitudes, Tensor LEED
C         beams only

      REAL AT,ATP
      DIMENSION AT(MKNT),ATP(MNPUN)
      COMPLEX XIST
      DIMENSION XIST(MNT0)

C  Tensor component production, and output

C  A0LM  : "initially incident spherical wave amplitudes (Bravais layer atom)"
C  ALM   : "final incident spherical wave amplitudes (Bravais layer atom)"
C  AEV,AOD: working space of final1 (even and odd lm-components of ALM)

      COMPLEX A0LM,ALM
      DIMENSION A0LM(MLMMAX),ALM(MLMMAX)
      COMPLEX AEV,AOD
      DIMENSION AEV(MLEV),AOD(MLOD)

Ctest
C  for outmat:

c      INTEGER NOUT
c      DIMENSION NOUT(5)

C  data statements

C  NCA IS DIMENSION OF CAA AS FUNCTION OF LMAX

      DATA NCA(1),NCA(2),NCA(3),NCA(4),NCA(5),NCA(6),NCA(7),
     +     NCA(8),NCA(9),NCA(10),NCA(11),NCA(12),NCA(13),
     +     NCA(14),NCA(15)/
     +     1,70,264,759,1820,3836,7344,13053,21868,34914,53560,
     +     79443,114492,160952,221408/

C  NLMS IS DIMENSION OF CLM AS FUNCTION OF LMAX

      DATA NLMS(1),NLMS(2),NLMS(3),NLMS(4),NLMS(5),NLMS(6),NLMS(7),
     +     NLMS(8),NLMS(9),NLMS(10),NLMS(11),NLMS(12),NLMS(13),
     +     NLMS(14),NLMS(15)/
     +     0,76,284,809,1925,4032,7680,13593,22693,36124,55276,
     +     81809,117677,165152,226848/

C  common blocks

      COMMON E,AK21,AK31,VPI,TV,EMACH
      COMMON /SL/ARA1,ARA2,ARB1,ARB2,RBR1,RBR2,NL1,NL2
      COMMON /MS/LMAX,EPS,LITER,SO1,SO2,SO3,SS1,SS2,SS3,SS4
      COMMON /ADS/FR,ASE,VPIS,VPIO,VO,VV
      COMMON /BT/IT,T,T0,DRPER,DRPAR,DR0
      COMMON /MPT/ NA,NS,ID,LAY,L1,NTAU,TSTS,TV1,DCUTS,NOPT,NEW

C  end header - begin actual computation

      EMACH = MEMACH

C  set integers required for variable dimensions throughout code
C  (cf. PARAM)

      IDEG = MIDEG         
      NL1  = MNL1
      NL2  = MNL2
      NL   = MNL

      KNBS = MKNBS
      KNT  = MKNT
      NPSI = MNPSI

      NPUN = MNPUN
      NT0  = MNT0

      NEL  = MNEL
      LMAX1= MLMAX1

      NSITE= MNSITE
      NLTYPE=MNLTYPE
      NNSUB= MNSUB
      NSTACK=MNSTACK

      NBRAV= MNBRAV 

C*************************************************************
C read information required for entire calculation
C*************************************************************

C  subroutine READIN reads lateral lattice information, general
C  calculational features, beam list and phase shifts

      CALL READIN(TITLE,TVA,RAR1,RAR2,TVB,IDEG,
     +NL,V,VL,JJS,KNBS,KNB,KNT,SPQF,KSYM,SPQ,TST,TSTS,NPUN,NPU,
     +THETA,FI,NPSI,ES,PHSS,L1,NEL,LMAX1,VPI,EI,EF,DE)

C  subroutine READTL reads information required for generic 
C  TENSOR LEED, i.e. Tensor LEED beam list and output format
C  READTL subroutine is presently only retained because superstructures
C  not known to the fd calculation can be calculated using Tensor LEED.
C  Since this feature is currently not implemented, the TLEED list PQFEX
C  is superfluous. Either include 'TLEED superstructures' later or do not
C  read in PQFEX, instead use appropriate beams NPU from fd list SPQF.

      CALL READTL(IFORM,NT0,PQFEX)

C  subroutine READGEO reads structural information for the reference
C  calculation (i.e. the vertical atomic arrangement within the 
C  2D lateral unit mesh)

      CALL READGEO(NEL,NSITE,NLTYPE,NNSUB,NSTACK,NBRAV,
     +             CONC,VIB,NSUB,LBRAV,LCOMP,LATT,STYPE,SUBPOS,
     +             TSLAB,ASA,TOPLAYB,BOTLAYB,ASBULK,LTYPE,LDIST,TENS,
     +             LAYFILE,ASB,NL)

C*************************************************************
C open output files, write title and general information
C*************************************************************

      OPEN(7,FILE='fd.out')
      CALL HEAD(7,TITLE,NPUN,NPU,KNT,SPQF,KSYM)

      IOUT = 11
      DO ISTACK = 1,NSTACK,1
        IF (TENS(ISTACK).eq.1) THEN
          DO ISUB = 1,NSUB(LTYPE(ISTACK)),1

            OUTNO(ISTACK,ISUB) = IOUT
            CALL OPENOUT(OUTNO(ISTACK,ISUB),LAYFILE(ISTACK,ISUB),IFORM)
            IOUT=IOUT+1

          ENDDO
        END IF
      ENDDO

C*************************************************************
C prepare some general quantities for the calculation
C*************************************************************

C  correct beam selection criterion TST (formerly done in READIN)

      IF ((TSLAB.eq.1).and.(NSTACK.eq.1)) THEN
        TST = 0.
      ELSE
        TST=ALOG(TST)/(AMIN1(ASB(1),ASA(1)))
        TST=TST*TST
      END IF

C  temperatures no longer used but carried through the code, 
C  thus set to default 100 - note they must be equal and not 0
C  if vibrational amplitudes from input are to be meaningful

      T0 = 100.
      T  = 100.

C  default code for fd. spectrum (fd.out)

      AIDENT = 0.0001

C  counters, variable array dimensions etc.

      LMMAX=(LMAX+1)**2
      KLM=(2*LMAX+1)*(2*LMAX+2)/2
      LEV=(LMAX+1)*(LMAX+2)/2
      LOD=LMMAX-LEV
      LEV2=2*LEV

C  compute CLEBSCH-GORDON COEFFICIENTS FOR COMPOSITE LAYER, CAA

      NCAA=NCA(LMAX)
      CALL CAAA(CAA,NCAA,LMMAX)

C  compute CLEBSCH-GORDON COEFFICIENTS FOR MATRICES X AND TAU, CLM

      N=2*LMAX+1
      NN=N*N
      NLM=NLMS(LMAX)
      NFAC = MNFAC
      CALL CELMG(CLM,NLM,YLMC,FAC2,NN,FAC1,N,LMAX,NFAC,FAC)

C  compute PERMUTATIONS OF (L,M)-SEQUENCE LX,LXI,LT,LXM 

      CALL LXGENT(LX,LXI,LT,LXM,LMAX,LMMAX)

C  compute CLEBSCH-GORDON COEFFICIENTS FOR COMPUTATION OF TEMPERATURE-
C  DEPENDENT PHASE SHIFTS, CPPP

      NN3=LMAX+1
      NN2=LMAX+1
      NN1=NN2+NN3-1

      CALL  CPPP (PPP,NN1,NN2,NN3,NFAC,FAC)

C*************************************************************
C  Loop over energies - calculate from initial to final energy
C  with a step width de
C*************************************************************

      EEV = EI

 100  CONTINUE

c  energy dependent inner potential is set in subroutine muftin

      call muftin(EEV,VO,VV,VPI,VPIS,VPIO)     

      VO   =     VO/HARTREE
      VV   =     VV/HARTREE
      VPI  = -  VPI/HARTREE
      VPIS = - VPIS/HARTREE
      VPIO = - VPIO/HARTREE

c  shift scattering energy by inner potential value VV

      E = EEV/HARTREE + VV

      write(6,*) "Energy in crystal: ", E, " Hartrees."
      write(6,*) "That's ", EEV, " eV in the vacuum."

C  SET LIMITING RADII ON LATTICE SUMS, POSSIBLY DIFFERENT FOR SUBSTRATE
C  (DCUTS) AND OVERLAYER (DCUTO) (currently DCUTO = DCUTS)
C  used in FMAT ...

      DCUTS =  - 5.0 * SQRT(2.0*E)/(AMIN1(VPIS, - 0.05))
      DCUTO =  - 5.0 * SQRT(2.0*E)/(AMIN1(VPIO, - 0.05))

C  compute components of incident vector parallel to surface

      KSQ = 2.0*(E-VV) 
      AK  = SQRT(KSQ) * SIN(THETA)
      AK2 = AK * COS(FI)
      AK3 = AK * SIN(FI)

C  compute components of Tensor LEED output beams parallel to
C  surface

      DO IBEAM = 1,NT0

        DO I = 1,2
          PSQ(I,IBEAM) = PQFEX(1,IBEAM) * RAR1(I) +
     +                   PQFEX(2,IBEAM) * RAR2(I)
        ENDDO

        AK2M(IBEAM) = - (AK2 + PSQ(1,IBEAM))
        AK3M(IBEAM) = - (AK3 + PSQ(2,IBEAM))

      ENDDO

C*************************************************************
C  Calculate atomic t-matrices for each site type
C*************************************************************

      DR0 = 0.
      VSITE = 0.

      DO ISITE = 1,NSITE

        DO IL = 1,L1,1
          TMAT(ISITE,IL) = 0.
        ENDDO

        DO IEL = 1,NEL,1

          IF (CONC(ISITE,IEL).gt.0.) THEN
            DRPER = VIB(ISITE,IEL)
            DRPAR = VIB(ISITE,IEL)

            CALL TSCATF(IEL,L1,ES,PHSS,NPSI,E,VSITE,PPP,NN1,NN2,NN3,
     +                  DR0,DRPER,DRPAR,T0,T,AF,CAF,NEL,LMAX1,PHS,DEL,
     +                  CTAB,SUM,BJ)                                  

            DO IL = 1,L1
              TMAT(ISITE,IL) = TMAT(ISITE,IL)+CONC(ISITE,IEL)*CAF(IL)
            ENDDO

          END IF

        ENDDO

      ENDDO

CTest
c      write(6,*) "atomic t-matrices: "
c      DO ISITE=1,NSITE
c        DO IL=1,LMAX1
c          write(6,*) ISITE,IL,TMAT(ISITE,IL)
c        ENDDO
c      ENDDO

C*****************************************************************
C  begin loop over Tensor LEED output beams (incl. incident beam, NEXIT = 0)
C  NEXIT
C*****************************************************************

      DO NEXIT = 0, NT0,1

C*************************************************************
C  determine what part of the computation is necessary for
C  current beam and energy
C  mark 200: skip computation up to second layer
C  mark 300: skip entire computation for current beam

C*************************************************************

      CALL DECIDE(NEXIT,KNOWN,EMERGE,PSQ1,PSQ2,AK2,AK3,AK21,AK31,
     +            KSQ,AK2M,AK3M,NT0,RBR1,RBR2)

      IF (KNOWN.eq.1) GO TO 200
      IF (EMERGE.eq.0) GO TO 300

C*************************************************************
C  begin computation of scattering properties up to 
C  2nd layer
C*************************************************************

      write(6,*) "Compute layer diffraction matrices, NEXIT =", NEXIT

C  choose appropriate beams

      CALL BEAMS(KNBS,KNB,SPQ,SPQF,KSYM,KNT,NPU,NPUN,AK21,AK31,E,TST,
     1NB,PQ,PQF,SYM,NPUC,MPU,NT,NP)

C  compute lattice sums
C  (NLS = NL is necessary, cf. appendix B of Van Hove / Tong 1979)

      NLS = NL

      CALL  FMAT (FLMS,V,JJS,NL,NLS,DCUTS,IDEG,LMAX,KLM,SCC,SA)

Ctest 
C  for outmat:

c      NOUT(1) = 1
c      NOUT(2) = 2
c      NOUT(3) = 3
c      NOUT(4) = 4
c      NOUT(5) = 5

C*************************************************************
C  compute (plane-wave) scattering matrices for each layer type
C  and (angular momentum) x-matrices 
C*************************************************************

      DO ILTYPE = 1,NLTYPE,1

C  initialize working matrices

            DO I = 1,MKNT
              DO J = 1,MKNT
                ROP(I,J) = 0.
                TOP(I,J) = 0.
                ROM(I,J) = 0.
                TOM(I,J) = 0.
              ENDDO
            ENDDO

C  initialize storage

            DO I = 1,MKNT
              DO J = 1,MKNT
                REPSTO(ILTYPE,I,J)  = 0.
                TRAPSTO(ILTYPE,I,J) = 0.
                REMSTO(ILTYPE,I,J)  = 0.
                TRAMSTO(ILTYPE,I,J) = 0.
              ENDDO
            ENDDO

        IF (NSUB(ILTYPE).eq.1) THEN

C          use rsmf() to compute scattering matrices for a bravais layer

C  set default values for temp-treatment - obsolete 
C  DRPER,DRPAR,DR0 not used - leave IT = 0 !

          IT=0                                                                   
          DRPER=VIB(STYPE(ILTYPE,1),1)
          DRPAR=VIB(STYPE(ILTYPE,1),1)

          DO IL=1,LMAX1
            CAF(IL) = TMAT(STYPE(ILTYPE,1),IL)
            AF(IL) = TMAT(STYPE(ILTYPE,1),IL)
          ENDDO

          ID  = 1
          NM  = NT
          LAY = LATT(ILTYPE)

C  slightly different calls for substrate or overlayer types

          IF (LAY.eq.1) THEN

            TV  = TVB                                                               
            NLL = 1
            NA  = 0
            NS  = 0
            NAA = NT

            CALL RSMF(ROP,TOP,ROP,TOP,ROP,TOP,ROP,TOP,NAA,NM,NM,
     +           AMULT,CYLM,PQ,SYM,NT,FLMS,FLM,V,NL,NA,NS,ID,NLL,
     +           AF,CAF,L1,LX,LXI,LMMAX,KLM,XEV,XOD,LEV,LOD,YLM,
     +           YLME,YLMO,IPLE,IPLO,CLM,NLM,LAY,XEVST,XODST)

          ELSE IF (LAY.eq.2) THEN

            TV  = TVA                                                               
            NLL = NL
            NA  = 0
            NS  = 0

            DO IBS = 1,KNBS

              NAA = NB(IBS)

              CALL RSMF(ROP,TOP,ROP,TOP,ROP,TOP,ROP,TOP,NAA,NM,NM,
     +           AMULT,CYLM,PQ,SYM,NT,FLMS,FLM,V,NL,NA,NS,ID,NLL,
     +           AF,CAF,L1,LX,LXI,LMMAX,KLM,XEV,XOD,LEV,LOD,YLM,
     +           YLME,YLMO,IPLE,IPLO,CLM,NLM,LAY,XEVST,XODST)

              NA  = NA + NB(IBS)
              NS  = NS + NB(IBS)

            ENDDO

          END IF

CTest
c          write(6,*) "Diffraction matrices of layer type ",ILTYPE
c          call OUTMAT(ROP,TOP,ROP,TOP,NOUT,NT)

          DO I = 1,MKNT
            DO J = 1,MKNT
              REPSTO(ILTYPE,I,J)  = ROP(I,J)
              TRAPSTO(ILTYPE,I,J) = TOP(I,J)
              REMSTO(ILTYPE,I,J)  = ROP(I,J)
              TRAMSTO(ILTYPE,I,J) = TOP(I,J)
            ENDDO
          ENDDO

C  store X-matrix separately

          IF ((NL.eq.1).or.(LAY.eq.1)) THEN
            DO ILEV = 1,MLEV
              DO JLEV = 1,MLEV
                XEVSTO(ILEV,JLEV,LBRAV(ILTYPE)) = XEVST(ILEV,JLEV)
              ENDDO
            ENDDO
            DO ILOD = 1,MLOD
              DO JLOD = 1,MLOD
                XODSTO(ILOD,JLOD,LBRAV(ILTYPE)) = XODST(ILOD,JLOD)
              ENDDO
            ENDDO
          END IF

        ELSE 

C  composite layer: rtinv
C  set some default quantities 

          ID   = 1
          NOPT = 1
          NEW  = 1

          NM   = NT

          NLAY = NSUB(ILTYPE)
          NLAY2= NLAY * (NLAY-1)/2
          LMG  = NLAY2*LMMAX*2 
          LMN  = NLAY * LMMAX
          LM2N = 2 * LMN
          LMNI = NLAY*LMMAX
          L2M  = MN

C  set sublattice characterisation

          NTAU = 0

          DO ISITE = 1,NSITE,1
            TAUNO(ISITE) = 0
          ENDDO

          DO ISUB=1,NSUB(ILTYPE),1

            IF (TAUNO(STYPE(ILTYPE,ISUB)).eq.0) THEN

              NTAU = NTAU+1
              TAUNO(STYPE(ILTYPE,ISUB)) = NTAU

              DO IL=1,LMAX1
                TSF(NTAU,IL) = TMAT(STYPE(ILTYPE,ISUB),IL)
              ENDDO

            END IF

            LPS(ISUB) = TAUNO(STYPE(ILTYPE,ISUB))

          ENDDO

C  copy sublayer positions onto correct grid POS

          CALL COPYPOS(POS,SUBPOS,NLTYPE,NNSUB,NLAY,ILTYPE)

          LMT = NTAU*LMMAX
          LAY = LATT(ILTYPE)

          IF (LAY.eq.1) THEN

            TV1 = TVB
            NA  = 0
            NS  = 0
            NAA = NT

            CALL RTINV(ROP,TOP,ROM,TOM,ROP,TOP,ROM,TOM,ROP,TOP,
     +      ROM,TOM,NAA,NM,NM,AMULT,CYLM,PQ,SYM,NT,FLMS,FLM,NL,
     +      LX,LXI,LT,LXM,LMMAX,KLM,XEV,LEV,LEV2,TAU,LMT,TAUG,
     +      TAUGM,CLM,NLM,POS,POSS,MGH,NLAY,DRL,SDRL,NUGH,NGEQ,
     +      NGOL,NLAY2,TEST,GH,LMG,RG,TS,LMN,TG,LM2N,VT,CAA,NCAA,
     +      TH,LMNI,IPL,TSF,LPS,LPSS,NORD,NNSUB,LMAX1,XH,HGHD,
     +      L2M,RG1,RG2,TENSOV(1,1,1,LCOMP(ILTYPE)))

          ELSE IF (LAY.eq.2) THEN

            TV1 = TVA
            NA = 0
            NS = 0

            DO IBS = 1,KNBS

              NAA = NB(IBS)

              CALL RTINV(ROP,TOP,ROM,TOM,ROP,TOP,ROM,TOM,ROP,TOP,
     +        ROM,TOM,NAA,NM,NM,AMULT,CYLM,PQ,SYM,NT,FLMS,FLM,NL,
     +        LX,LXI,LT,LXM,LMMAX,KLM,XEV,LEV,LEV2,TAU,LMT,TAUG,
     +        TAUGM,CLM,NLM,POS,POSS,MGH,NLAY,DRL,SDRL,NUGH,NGEQ,
     +        NGOL,NLAY2,TEST,GH,LMG,RG,TS,LMN,TG,LM2N,VT,CAA,NCAA,
     +        TH,LMNI,IPL,TSF,LPS,LPSS,NORD,NNSUB,LMAX1,XH,HGHD,
     +        L2M,RG1,RG2,TENSOV(1,1,1,LCOMP(ILTYPE)))

              NA  = NA + NB(IBS)
              NS  = NS + NB(IBS)

            ENDDO

          END IF

CTest
c          write(6,*) "Diffraction matrices of layer type ",ILTYPE
c          call OUTMAT(ROP,TOP,ROM,TOM,NOUT,NT)

          DO I = 1,MKNT
            DO J = 1,MKNT
              REPSTO(ILTYPE,I,J)  = ROP(I,J)
              TRAPSTO(ILTYPE,I,J) = TOP(I,J)
              REMSTO(ILTYPE,I,J)  = ROM(I,J)
              TRAMSTO(ILTYPE,I,J) = TOM(I,J)
            ENDDO
          ENDDO

        END IF

      ENDDO

C*************************************************************
C  next, compute intra-crystal scattering matrices up to 
C  second layer in plane wave representation
C*************************************************************

      IF (TSLAB.eq.0) THEN

C  compute bulk reflection matrix using subras instead of subref for generality
C  Since an alternating stacking sequence of two different layer types must be
C  possible in bulk (e.g. hcp stacking sequence, ordered alloys), a call to
C  DBLGAS is necessary before!)

        DO IBEAM = 1,MKNT
          DO JBEAM = 1,MKNT

            RAMP(IBEAM,JBEAM) = REPSTO(TOPLAYB,IBEAM,JBEAM)
            TAPP(IBEAM,JBEAM) = TRAPSTO(TOPLAYB,IBEAM,JBEAM)
            RAPM(IBEAM,JBEAM) = REMSTO(TOPLAYB,IBEAM,JBEAM)
            TAMM(IBEAM,JBEAM) = TRAMSTO(TOPLAYB,IBEAM,JBEAM)

            RBMP(IBEAM,JBEAM) = REPSTO(BOTLAYB,IBEAM,JBEAM)
            TBPP(IBEAM,JBEAM) = TRAPSTO(BOTLAYB,IBEAM,JBEAM)
            RBPM(IBEAM,JBEAM) = REMSTO(BOTLAYB,IBEAM,JBEAM)
            TBMM(IBEAM,JBEAM) = TRAMSTO(BOTLAYB,IBEAM,JBEAM)

          ENDDO
        ENDDO

        NA = 0

        CALL DBLGAS(RAMP,TAPP,RAPM,TAMM,RBMP,TBPP,RBPM,TBMM,ASBULK,NT,
     +              NA,PQ,NT,S1,S2,S3,S4,XS,PP,JNT,NT,EMACH)

        CALL SUBRAS(RAMP,TAPP,RBPM,TBMM,NT,S1,S2,S3,S4,PP,XS,JNT,
     +              NT,NA,NT,PQ,ASA)

        DO IBEAM = 1,MKNT,1
          DO JBEAM = 1,MKNT,1
            REABOV(IBEAM,JBEAM) = RAMP(IBEAM,JBEAM)
          ENDDO
        ENDDO

      ELSE

        DO IBEAM = 1,MKNT,1
          DO JBEAM = 1,MKNT,1
C           JBEAM produced my headache ??
            REABOV(IBEAM,JBEAM) = CMPLX(0.,0.)
          ENDDO
        ENDDO

      END IF

C  stack remaining layers up to second layer, feeding REBELO with the
C  reflection matrix of the previous step, REABOV, and save the effective
C  scattering matrices of each layer

CTest
c          write(6,*) "Bulk reflection matrix"
c          call OUTMAT(REABOV,REABOV,REABOV,REABOV,NOUT,NT)

      DO ISTACK = NSTACK,2,-1

        DO IBEAM = 1,MKNT
          DO JBEAM = 1,MKNT

            RAMP(IBEAM,JBEAM) = REPSTO(LTYPE(ISTACK),IBEAM,JBEAM)
            TAPP(IBEAM,JBEAM) = TRAPSTO(LTYPE(ISTACK),IBEAM,JBEAM)
            RAPM(IBEAM,JBEAM) = REMSTO(LTYPE(ISTACK),IBEAM,JBEAM)
            TAMM(IBEAM,JBEAM) = TRAMSTO(LTYPE(ISTACK),IBEAM,JBEAM)

            REBELO(IBEAM,JBEAM) = REABOV(IBEAM,JBEAM)

          ENDDO
        ENDDO

        DO I = 1,3
          ASV(I) = LDIST(ISTACK,I)
        ENDDO

        NA = 0

        CALL DBGT(RAMP,TAPP,RAPM,TAMM,REBELO,REABOV,ASV,NT,
     +            NA,PQ,NT,EFSCAM(1,1,ISTACK),EFSCAP(1,1,ISTACK),
     +            XS,PP,JNT,NT,EMACH)

CTest
c          write(6,*) "reflection matrix with layer ",ISTACK
c          call OUTMAT(REABOV,REABOV,REABOV,REABOV,NOUT,NT)

      ENDDO

C*************************************************************
C  continue computation here if diffraction matrices up to 2nd
C  layer were already known
C*************************************************************

 200    CONTINUE

C*************************************************************
C  add top layer and compute outgoing amplitudes XI
C*************************************************************

        DO IBEAM = 1,MKNT
          DO JBEAM = 1,MKNT

            RAMP(IBEAM,JBEAM) = REPSTO(LTYPE(1),IBEAM,JBEAM)
            TAPP(IBEAM,JBEAM) = TRAPSTO(LTYPE(1),IBEAM,JBEAM)
            RAPM(IBEAM,JBEAM) = REMSTO(LTYPE(1),IBEAM,JBEAM)
            TAMM(IBEAM,JBEAM) = TRAMSTO(LTYPE(1),IBEAM,JBEAM)

          ENDDO
        ENDDO

        DO I = 1,3
          ASV(I) = LDIST(1,I)
        ENDDO

      CALL ADREF2T (RAMP,TAPP,RAPM,TAMM,REABOV,S1,WK,XS,JNT,XI,
     +PQ,NT,EMACH,ASV,APLUS,AMINUS,ADOWN,PSQ1,PSQ2)           

C*************************************************************
C  full dyn. part of calculation finished here - if NEXIT = 0,
C  write all available output before proceeding
C*************************************************************

      IF (NEXIT.eq.0) THEN

C  7 is output unit - file fd.out

        CALL RINT(NT,XI,AT,ATP,PQ,PQF,SYM,VV,THETA,FI,MPU,NPUC,EEV,
     +            AIDENT,1,7)

        DO ISTACK = 1,NSTACK,1
          IF (TENS(ISTACK).eq.1) THEN
            DO ISUB = 1,NSUB(LTYPE(ISTACK)),1

              DO IL=1,LMAX1
                CAF(IL) = TMAT(STYPE(LTYPE(ISTACK),ISUB),IL)
              ENDDO

              CALL OUTXIST(OUTNO(ISTACK,ISUB),IFORM,E,PQF,PQFEX,
     +                     NT0,NT,XI,XIST,L1,CAF)

            ENDDO
          END IF
        ENDDO

      END IF

C******************************************************************
C  next, produce and write Tensor components: output wavefield ALM
C  read Rous's PhD thesis for that, but beware a mix-up: his formula
C  is derived for an initial definition of delta t, ALM is correct for 
C  his later definition (and so is his computational chapter)
C******************************************************************

      DO ISTACK = 1,NSTACK,1
        
        IF (ISTACK.ge.2) THEN

C  compute (plane-wave) wave-field onto and from current layer using the
C  known wave-field from above

          DO IBEAM = 1,MKNT

            APLUS(IBEAM) = ADOWN(IBEAM)

          ENDDO

          CALL MULTAMP(APLUS,AMINUS,EFSCAM(1,1,ISTACK),NT)
          CALL MULTAMP(APLUS,ADOWN,EFSCAP(1,1,ISTACK),NT)

        END IF

C  produce ALM and output if desired

        IF (TENS(ISTACK).eq.1) THEN

          IF (NSUB(LTYPE(ISTACK)).eq.1) THEN

            ISUB = 1

            DO ILEV = 1,MLEV
              DO JLEV = 1,MLEV
                XEVST(ILEV,JLEV) = 
     +          XEVSTO(ILEV,JLEV,LBRAV(LTYPE(ISTACK)))
              ENDDO
            ENDDO
            DO ILOD = 1,MLOD
              DO JLOD = 1,MLOD
                XODST(ILOD,JLOD) = 
     +          XODSTO(ILOD,JLOD,LBRAV(LTYPE(ISTACK)))
              ENDDO
            ENDDO

            CALL FINAL1(ALM,A0LM,APLUS,AMINUS,CYLM,XODST,XEVST,
     +                  IPLO,IPLE,LMAX,LMMAX,LOD,LEV,NT,LX,AEV,
     +                  AOD,EMACH)

            IF (NEXIT.eq.0) THEN
              CALL OUTAMP(OUTNO(ISTACK,ISUB),IFORM,NEXIT,0.,
     +                  0.,ALM,LMMAX,AK2,AK3)
            ELSE
              CALL OUTAMP(OUTNO(ISTACK,ISUB),IFORM,NEXIT,PSQ(1,NEXIT),
     +                  PSQ(2,NEXIT),ALM,LMMAX,AK2M(NEXIT),AK3M(NEXIT))
            END IF

          ELSE

            LMN  = NSUB(LTYPE(ISTACK)) * LMMAX

            DO ISUB = 1,NSUB(LTYPE(ISTACK))

              DO IL=1,LMAX1
                CAF(IL) = TMAT(STYPE(LTYPE(ISTACK),ISUB),IL)
              ENDDO

              CALL FINALOV(TENSOV(1,1,1,LCOMP(LTYPE(ISTACK))),ALM,
     +             APLUS,AMINUS,NT,ISUB,LMN,LMMAX,CAF,LMAX,LXM,E,VPI)

              IF (NEXIT.eq.0) THEN
                CALL OUTAMP(OUTNO(ISTACK,ISUB),IFORM,NEXIT,0.,
     +                  0.,ALM,LMMAX,AK2,AK3)
              ELSE
                CALL OUTAMP(OUTNO(ISTACK,ISUB),IFORM,NEXIT,
     +                      PSQ(1,NEXIT),PSQ(2,NEXIT),ALM,LMMAX,
     +                      AK2M(NEXIT),AK3M(NEXIT))
              END IF

            ENDDO

          END IF

        END IF

      ENDDO

C*************************************************************
C  end NEXIT LOOP 
C*************************************************************

 300    CONTINUE

      ENDDO

C*************************************************************
C  delimit output files after each energy
C*************************************************************

      DO ISTACK = 1,NSTACK,1
        IF (TENS(ISTACK).eq.1) THEN
          DO ISUB = 1,NSUB(LTYPE(ISTACK)),1
            CALL DELIMIT(OUTNO(ISTACK,ISUB),IFORM)
          ENDDO
        END IF
      ENDDO

C*************************************************************
C  end energy loop - next energy, or terminate, if E.gt.EF
C*************************************************************

      EEV = EEV+DE

      IF (EEV.le.EF) THEN
        GO TO 100
      ELSE
        write(6,*) "Program ends with energy ", EEV-DE," eV."
      END IF       

C*************************************************************
C  that's all, folks!
C*************************************************************

      END
