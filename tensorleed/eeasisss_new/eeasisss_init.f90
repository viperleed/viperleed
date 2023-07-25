  !=======================================================================
  !This program is free software under the terms of the GNU General Public
  !License as published by the Free Software Foundation.
  !Author: John O. Rundgren, jru@KTH.se ,
  !        KTH Royal Institute of Technology, Stockholm, Sweden.
  !Version: 28 March 2021.
  !-----------------------------------------------------------------------
  ! Adapted for ViPErLEED by Alexander M. Imre, October 2021
  !-----------------------------------------------------------------------
 module param
 !double precision definition.
 integer,parameter :: dp=selected_real_kind(15,307)
 !
 !I/O items.
 character(len=127) :: logfil,outdir,atdir
 integer :: logunit,thread
 !
 !basic constants.
 real(dp),allocatable :: z(:)
 character(len=2),parameter :: Elem(92)=(/&
 'H ','He','Li','Be','B ','C ','N ','O ','F ','Ne','Na','Mg','Al','Si',&
 'P ','S ','Cl','Ar','K ','Ca','Sc','Ti','V ','Cr','Mn','Fe','Co','Ni',&
 'Cu','Zn','Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y ','Zr','Nb','Mo',&
 'Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I ','Xe','Cs','Ba',&
 'La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb',&
 'Lu','Hf','Ta','W ','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po',&
 'At','Rn','Fr','Ra','Ac','Th','Pa','U '/)
  real(dp),parameter :: bohr=0.5291772083d0,rydb=13.60569172d0

 !select program operations.
 character(len=1) :: bos,SpinPS,Chg,Pot,WF
 character(len=80) :: compound
 logical,parameter :: linux=.true., windows=.false. 

 !Energy input.
 !emesh = mesh of uniform energy grid
 integer  :: lmax,nthread
 real(dp) :: eev1,eev2,emesh,relerr,abserr

 !idZ(ir)     = nuclear number,
 !idA(ir)     = atom A's number in logfile print-out,
 !idZA(ir)    = atom identifier = idZ(ir)//'.'//idA(ir)
 !idElemA(ir) = Elem identifier = Elem(nint(idZ(ir)))//'.'//idA(ir),
 !            see subroutine(StructureAndMethods)
 integer,allocatable :: idZ(:),idA(:)
 character,allocatable :: idZA(:)*5,idElemA(:)*5
end module param

!======================================================================
module maths
use param,only : dp
contains
  !--------------------------------------------------------------------   
  subroutine trapez(rx,dx,nx,f,j1,j2,s)
  !integration by trapezoidal rule.
  !interval [j1,j2] is subinterval of [1,nx].
  implicit none
  integer :: nx,j1,j2,i
  real(dp) :: rx(nx),f(nx),s(nx),dx,dxh,a1,a2
  dxh=dx*0.5d0
  s(j1)=0; a1=f(j1)*rx(j1)
  a1 = 0.d0 ! was not initialized
  do i=j1+1,j2
    a2=f(i)*rx(i)
    s(i)=s(i-1)+a1+a2
    a1=a2
  enddo
  s(j1:j2)=s(j1:j2)*dxh
  return
  end subroutine trapez
  !--------------------------------------------------------------------   
  function sph_av(rs,R,rx,dx,nx)
  !spherical average by spherical trapezoidal rule.
  implicit none
  integer :: nx,i
  real(dp) :: sph_av,rs(nx),R,rx(nx),f(nx),s(nx),dx,dxh,a1,a2,p
  !position of R in exponential grid.
  p = (log(R/rx(1)))/dx + 1.d0
  i = int(p)
  !position of R relative to next lower grid point i.
  p = p-i
  dxh=dx*0.5d0
  !rx**2 for spherical average and **1 for exponential grid.
  do i=1,nx
    f(i) = rs(i)*rx(i)**3
  enddo
  s(1) = 0.d0
  a1 = f(1)
  do i=2,nx-1
    a2 = f(i)
    s(i) = s(i-1)+a1+a2
    a1 =a2
  enddo
  a2 = f(nx-1) + (f(nx)-f(nx-1)) * p/(rx(nx)-rx(nx-1))
  s(nx) = s(nx-1)+a1+a2
  s(1:nx) = s(1:nx)*dxh
  !4*pi in nominator and denominator is cancelled.
  sph_av = s(nx) * 3.d0/R**3
  return
  end function sph_av
  !------------------------------------------------------------------- 
  subroutine pois(rx,dx,nx,z,rho4pir2,V)
  !integrating Poissson's equation calculating V.
  implicit none
  integer :: nx,i
  real(dp) :: rx(nx),dx,z,rho4pir2(nx),V(nx)
  real(dp) :: a1,a2,b1,b2,c,Q(nx),s(nx)
  c=0.5d0*dx
  Q(1)=0.d0; a1=rho4pir2(1)*rx(1)
  s(1)=0.d0; b1=rho4pir2(1)
  do i=2,nx
    a2=rho4pir2(i)*rx(i)
    Q(i)=Q(i-1)+a1+a2
    a1=a2
    b2=rho4pir2(i)
    s(i)=s(i-1)+b1+b2
    b1=b2
  enddo
  Q=Q*c
  s=s*c
  V=2.d0*(Q/rx-s+s(nx)-z/rx) !in Ry.
  return
  end subroutine pois
  !-------------------------------------------------------------------
  function ipol(x,nx,r,d,f,messg)
  use param,only : dp
  implicit none
  integer  :: nx,i1  
  real(dp) :: ipol,x,r(nx),d,f(nx),p1,p
  character(len=*) :: messg
  !ipol is exponential-grid lagrangian 3-point interpolation.
  !Ref.: Abramowitz and Stegun, Sec. 25.2.11.
  !
  if(x-r(nx)>0.d0)then
    write(611,900)'egl3 stop: x-r(nx),r(nx)=',x-r(nx),r(nx)
    write(611,902)'error location: ',messg
    stop
  endif
  p1=(log(x/r(1)))/d+1.d0
  i1=nint(p1)
  if(i1==nx) i1=i1-1
  p=p1-i1
  ipol = &
  0.5d0*p*(p-1.d0)*f(i1-1)+(1.d0-p*p)*f(i1)+0.5d0*p*(p+1.d0)*f(i1+1)
  return
  900 format(a,es9.2,f7.4)
  902 format(2a)
  end function ipol
  !-----------------------------------------------------------------
  !finding subscript .ge. that of x in array a.
  function iofx(x,a,la,ha)
  implicit none
  integer :: iofx,i,il,ih,la,ha
  real(dp) :: x,a(la:ha)
  il=la; ih=ha
1 if(ih-il>1)then
    i=(ih+il)/2
    if(a(i)>x)then; ih=i; else; il=i; endif
    goto 1
  endif 
  iofx=ih    
  return
  end function iofx
end module maths

!==============================================================================
module space
use param,only : dp
implicit none
!nieq   = # inequivalent atoms in the structure,
!neq(i) = # equivalent atoms for i=1:nieq,
!nlat   = # lattice points in the unit cell,
!nshell = # neighbor atomic shells, >= 2,
integer,parameter :: nshell=8
integer :: nieq,nlat
integer,allocatable :: neq(:),ieq(:)
!rx(i,ir)  )= radial grid rx(i,ir)=r(1,ir)*exp(dx(ir)*(i-1)),
!             where ir=1:nieq and i=1:nx(ir),

!rmt(:)  = MT radius,
!rmtovl  = overlap exterior to NN distance ad(2,:),
!rmin(:) = min MT radius,
!rmax(:) = max MT radius,
!rx(:,:) = radial grid,
!dx(:)   = logarithmic increment in radial grid,
!nx(:)   = (input) # grid points in atomic calculation,
!        = (calc.) outermost grid point determined by superposition z,
!nxx     = max nx(:).
!real(dp) :: rmtovl
real(dp),allocatable :: rx(:,:),dx(:),rmin(:),rmax(:),rmt(:),rmtovl(:)
integer,allocatable  :: nx(:)
integer :: nxx
!unit cell,
!rcX(i,j)=the i'th coordinate of the j'th axis of the unit cell,
!rk(i,j) =the i'th coordinate of the j'th atom in the unit cell,
!volUC   =volume of unit cell (inclusive vacuum slab),
!volXC   =volume of bulk unit cell (filled with crystal electrons).
!neighbors in crystal,
!  ia(j) = atomic type in j'th shell,
!  na(j) = number of atoms in j'th shell,
!  ad(j) = distance to j'th shell,
!  ncon  = number of shells included,
integer,allocatable ::  na(:,:),ia(:,:),ncon(:)
real(dp),allocatable :: rk(:,:),ad(:,:)
real(dp) :: rc(3,3),volUC,volXC
contains
 !----------------------------------------------------------------------
 !NEIGHBOR_SHELLS is nearest neighbor data for atoms in a crystal
 !structure.
 !Ref.: A. Barbieri and M.A. Van Hove, Phase shift calculation package,
 !      http://www.icts.hkbu.edu.hk/vanhove/
 !The main part of the subroutine is a fortran90 editon of Ref.

 !input:
 !rc(i,j) =the i'th coordinate of the j'th axis of the unit cell,
 !rk(i,j) =the i'th coordinate of the j'th atom in the unit cell,
 !neq(ir) =the number of equivalent atoms of type ir.
 !output:
 !ia(j,ir)=the type of atoms in the j'th neighbor shell,
 !na(j,ir)=the number of atoms in the j'th shell,
 !ad(j,ir)=the radius of the j'th shell,
 !ncon(ir)=no. of ir type shells.
 subroutine NeighborShells
 use param,only : dp
 use param,only : z,Elem
 implicit none
 integer :: j,jr,jjr,k,kr,ic,iic,jjc,jc,ir,nc,kc,i,jx,jy,jz
 real(dp) :: rj(3),dr,r,a,b
 real(dp),parameter :: bohr=0.5291772083d0
 real(dp),parameter :: rlimit=15.d0
 !
 allocate(ad(nshell,nieq),na(nshell,nieq),ia(nshell,nieq),ncon(nieq))
 ia=0
 na=0
 ad=1.d+33
 do jx=-nshell,nshell
  do jy=-nshell,nshell
    do jz=-nshell,nshell
      do j=1,3
        rj(j)=jx*rc(j,1)+jy*rc(j,2)+jz*rc(j,3)
      enddo
      !rj is current unit cell origin. for each atom in this unit cell,
      !find displacement r from kr-type atom in basic unit cell.
      j=0
      do jr=1,nieq
        do jjr=1,neq(jr)
          j=j+1
          k=1
          do kr=1,nieq
            r=sqrt((rj(1)+rk(1,j)-rk(1,k))**2+&
                   (rj(2)+rk(2,j)-rk(2,k))**2+&
                   (rj(3)+rk(3,j)-rk(3,k))**2 )
            !fixed limit for r.
            if(r>rlimit)then
              k=k+neq(kr)
              cycle
            endif
            !compare r with nearest neighbor distances already found.
            ic=0
            do
              ic=ic+1
              if(ic>nshell)then
                k=k+neq(kr)
                exit
              endif
              dr=r-ad(ic,kr)
              if(dr >= 1.d-03)then
                cycle
              elseif(abs(dr) < 1.d-04)then
                if(ia(ic,kr)/=jr)then
                  cycle
                else
                  na(ic,kr)=na(ic,kr)+1
                  k=k+neq(kr)
                  exit
                endif
              elseif(dr <= -1.d-04)then
                if(ic<nshell)then
                  iic=ic+1
                  do jjc=iic,nshell
                    jc=nshell+iic-jjc
                    ia(jc,kr)=ia(jc-1,kr)
                    na(jc,kr)=na(jc-1,kr)
                    ad(jc,kr)=ad(jc-1,kr)
                  enddo
                endif
                ia(ic,kr)=jr
                na(ic,kr)=1
                ad(ic,kr)=r
                k=k+neq(kr)
                exit
              endif
            enddo !ic
          enddo !kr
        enddo !jjr
      enddo !jr
    enddo !jz
  enddo !jy
 enddo !jx
 do ir=1,nieq
   ncon(ir)=0
   do ic=1,nshell
     if(na(ic,ir)==0) cycle
     ncon(ir)=ncon(ir)+1
   enddo
 enddo
 !display.
 write(61,'(/a)')'NNshells from JR.'
 do ir=1,nieq
  j=0
  do i=2,ncon(ir)
    j=j+na(i,ir)  !j counts the neighbors about atom ir.
  enddo
  write(61,'(72("="))')
  write(61,200) ir,Elem(nint(z(ir))),ncon(ir)-1,j
  nc=ncon(ir)
  ic=(nc-1)/8+1
  kc=0
  do i=1,ic
    jc=kc+1
    kc=min(nc,kc+8)
    write(61,208) (ia(j,ir),Elem(nint(z(ia(j,ir)))),j=jc,kc)
    write(61,206) (na(j,ir),j=jc,kc)
    write(61,204) (ad(j,ir),j=jc,kc)
    write(61,202) (ad(j,ir)*bohr,j=jc,kc)
    write(61,*)
  enddo !i
  write(61,'(a,20i2)')'ia',(ia(j,ir),j=1,kc)
  write(61,'(a,20i2)')'na',(na(j,ir),j=1,kc)
  write(61,*)
 enddo !ir
 !
 write(61,'(a)')'NN distances (B).'
 write(61,212) (i,i=1,nieq)
 write(61,214) (ad(2,i),i=1,nieq)
 return
 !
 !calculating volXC.
 !dz of volUC.
 a=rc(3,3)-rc(3,1)
 !dz of volXC.
 b=maxval(rk(3,:))-minval(rk(3,:))+0.5d0*(ad(2,1)+ad(2,nieq))
 volXC=volUC*(b/a)
 write(61,'(/a,2(1x,f0.2))')'volUC,volXC(B**3) =',volUC,volXC
 return
 200 format('   atom   ',i2,'.',a,i6,' neighbor shells,',i3,' neighbors')
 202 format('dist.(A)',10f8.4:)
 204 format('dist.(B)',10f8.4:)
 206 format('  number ',10(3x,i2,3x:))
 208 format('      NN  ',10(i2,'.',a,3x:))
 212 format(10(2x,i2,4x:))
 214 format(10(f8.4:))
 end subroutine NeighborShells
!----------------------------------------------------------------------
      SUBROUTINE NBR !(IA,NA,AD,NCON,NRR,NR,RC,RK,N,RMAX,MC)
      use param,only : dp
      use param,only : z,Elem
      implicit none
      integer :: MC,LIMC,NR,NC,I,IIC,IC,IR,J,JJC,JX,JY,JZ,K,KC,KR, &
        ITC,JJR,JC,JNR,JR
      real(dp),parameter :: bohr=0.5291772083d0
      real(dp) :: RAD,A1,A2,A3,RCMIN,RMAX,AS,AX,AY,AZ,R,RJ(3),DR
     !DIMENSION IA(MC,NR),NA(MC,NR),AD(MC,NR),NCON(NR),NRR(NR)
     !dimension rj(3),RC(3,3),RK(3,nlat)
     !NR=nieq
     !NRR(NR)=neq(nieq)
     !MC=nshell = 20 !Barbieri-VanHove
  !allocate(ad(nshell,nieq),na(nshell,nieq),ia(nshell,nieq),ncon(nieq))
!
!  ROUTINE TO SUPPLY NEAREST NEIGHBOUR DATA FOR ATOMS IN
!  A CRYSTAL STRUCTURE, GIVEN?
!  RC(I,J)? THE I'TH COORDINATE OF THE J'TH AXIS OF THE UNIT CELL
!  RK(I,J)? THE I'TH COORDINATE OF THE J'TH ATOM IN THE UNIT CELL
!  NRR(IR)? THE NUMBER OF TYPE-IR ATOMS IN THE UNIT CELL
!  THE INFORMATION RETURNED, FOR A TYPE-IR ATOM, IS
!  NCON(IR)? THE NUMBER OF NEAREST NEIGHBOUR SHELLS OF A TYPE-IR
!  ATOM INCLUDED, OUT TO A DISTANCE OF RMAX, BUT .LE. MC
!  IA(J,IR)? THE TYPE OF ATOMS IN THE J'TH NEIGHBOURING SHELL
!  NA(J,IR)? THE NUMBER OF ATOMS IN THE J'TH SHELL
!  AD(J,IR)? THE RADIUS OF THE J'TH SHELL
!
!  INITIALISATION
      RAD(A1,A2,A3)=SQRT(A1*A1+A2*A2+A3*A3)
      MC=nshell
      NR=nieq
      rmax=10.d0
      RCMIN=1.0E6
      DO 1 I=1,3
1     RCMIN=min(RCMIN,RAD(RC(1,I),RC(2,I),RC(3,I)))
      DO 2 IR=1,NR
      DO 2 IC=1,MC
      IA(IC,IR)=0
      NA(IC,IR)=0
2     AD(IC,IR)=1.0d06
!  SEARCH OVER ADJACENT UNIT CELLS TO INCLUDE MC NEAREST NEIGHBOURS
      ITC=INT(RMAX/RCMIN)+1
      LIMC=ITC+ITC+1
      AS=-FLOAT(ITC+1)
      AX=AS
      DO 10 JX=1,LIMC
      AX=AX+1.0d0
      AY=AS
      DO 10 JY=1,LIMC
      AY=AY+1.0d0
      AZ=AS
      DO 10 JZ=1,LIMC
      AZ=AZ+1.0d0
      DO 3 J=1,3
3     RJ(J)=AX*RC(J,1)+AY*RC(J,2)+AZ*RC(J,3)
!  RJ IS CURRENT UNIT CELL ORIGIN.  FOR EACH ATOM IN THIS UNIT CELL
!  FIND DISPLACEMENT R FROM KR-TYPE ATOM IN BASIC UNIT CELL
      J=0
      DO 10 JR=1,NR
      JNR=neq(JR)  !NRR=neq
      DO 10 JJR=1,JNR
      J=J+1
      K=1
      DO 9 KR=1,NR
      R=RAD(RJ(1)+RK(1,J)-RK(1,K), &
            RJ(2)+RK(2,J)-RK(2,K), &
            RJ(3)+RK(3,J)-RK(3,K))
      IF(R.GT.RMAX)GOTO 9
!  COMPARE R WITH NEAREST NEIGHBOUR DISTANCES ALREADY FOUND
      IC=0
4     IC=IC+1
      IF(IC.GT.MC) GOTO 9
      DR=R-AD(IC,KR)
      IF(ABS(DR).LT.1.0d-04) DR=0.0d0
      IF(DR) 6,5,4
5     IF(IA(IC,KR).NE.JR) GOTO 4
      NA(IC,KR)=NA(IC,KR)+1
      GOTO 9
6     IF(IC.EQ.MC) GOTO 8
      IIC=IC+1
      DO 7 JJC=IIC,MC
      JC=MC+IIC-JJC
      IA(JC,KR)=IA(JC-1,KR)
      NA(JC,KR)=NA(JC-1,KR)
7     AD(JC,KR)=AD(JC-1,KR)
8     IA(IC,KR)=JR
      NA(IC,KR)=1
      AD(IC,KR)=R
9     K=K+neq(KR)  !NRR=neq
10     CONTINUE
      DO 12 IR=1,NR
      NCON(IR)=0
      DO 11 IC=1,MC
      IF(NA(IC,IR).EQ.0) GOTO 12
11    NCON(IR)=NCON(IR)+1
12    CONTINUE
 !display.
 write(61,'(/a)')'NBshells from BVH.'
 do ir=1,nieq
   write(61,'(80("-"))')
   j=0
   do i=2,ncon(ir)
     j=j+na(i,ir)  !j counts the neighbors about atom ir.
   enddo
   write(61,200) ir,Elem(nint(z(ir))),ncon(ir)-1,j
   nc=ncon(ir)
   ic=(nc-1)/8+1
   kc=0
   do i=1,min(ic,nshell)
     jc=kc+1
     kc=min(nc,kc+8)
     write(61,208) (ia(j,ir),Elem(nint(z(ia(j,ir)))),j=jc,kc)
     write(61,206) (na(j,ir),j=jc,kc)
     write(61,204) (ad(j,ir),j=jc,kc)
     write(61,202) (ad(j,ir)*bohr,j=jc,kc)
     write(61,*)
   enddo
   write(61,'(a,20i2)')'ia',(ia(j,ir),j=1,kc)
   write(61,'(a,20i2)')'na',(na(j,ir),j=1,kc)
   write(61,*)
 enddo
 200 format('   atom   ',i2,'.',a,i6,' neighbor shells,',i3,' neighbors')
 202 format('dist.(A)',10f8.4:)
 204 format('dist.(B)',10f8.4:)
 206 format('  number ',10(3x,i2,3x:))
 208 format('      NN  ',10(i2,'.',a,3x:))
      RETURN
      END subroutine NBR
end module space
!==============================================================================
module enrgy
use param,only : dp,z,lmax,eev1,eev2,emesh,linux,windows
use space,only : rx,dx,nx,rmt,ad,ia,neq,nlat
implicit none
!Vnist = initial atomic potential corresponding to rho,
!      = tends to zero at radius infty, Ref. Eric Shirley, NIST,
real(dp) :: Vcry0
real(dp),allocatable :: atrho4pir2(:,:),rho(:,:),rs(:,:),&
  QvsR(:,:),Vnist(:,:),Vcry(:,:),VcryR(:)
!
!Energy grid items.
!eev(1:ne) = energy grid (in eV),
!eev1,eev2 = energy interval for phase shift calculation,
!emv0      = Einc-V0(inner potential),
!m         = energy position of V0 dip,
!xcfac     = correction factor on Vxc.
integer  :: ne,m,kappa,iatom
real(dp) :: Einc
real(dp),allocatable :: eev(:),emv0(:),v0ev(:),xcfac(:)
!
!Phase shift quantum numbers etc.
!iatom =atom subscript,
!kappa =spin-orbit coupling quantum number (dirac eq.):
!      =-l-1, l=0,1,2,... and =l, l=1,2,...(dirac eq.),
!      =l, l=0,1,2,... (schroedinger eq.),
!
!anu = atomic mass of chemical elements (1:92).
real(dp),dimension(92),parameter :: anu=(/&
  1.d0,  4.d0,  7.d0,  9.d0, 11.d0, 12.d0, 14.d0, 16.d0, 19.d0, 20.d0, 23.d0, 24.d0, 27.d0, 28.d0,&
 31.d0, 32.d0, 35.d0, 40.d0, 39.d0, 40.d0, 45.d0, 48.d0, 51.d0, 52.d0, 55.d0, 56.d0, 59.d0, 59.d0,&
 64.d0, 65.d0, 70.d0, 73.d0, 75.d0, 79.d0, 80.d0, 84.d0, 85.d0, 88.d0, 89.d0, 91.d0, 93.d0, 96.d0,&
 98.d0,101.d0,103.d0,106.d0,108.d0,112.d0,115.d0,119.d0,122.d0,128.d0,127.d0,131.d0,133.d0,137.d0,&
139.d0,140.d0,141.d0,144.d0,145.d0,150.d0,152.d0,157.d0,159.d0,163.d0,165.d0,167.d0,169.d0,173.d0,&
175.d0,178.d0,181.d0,184.d0,186.d0,190.d0,192.d0,195.d0,197.d0,201.d0,204.d0,207.d0,209.d0,209.d0,&
210.d0,222.d0,223.d0,226.d0,227.d0,232.d0,231.d0,238.d0/)
contains
 !----------------------------------------------------------------------
 subroutine Free_Atom_Overlap
 !Free_Atom_Overlap uses overlapping scatterers whose radii rmt are 
 !variables fitted to experimental intensity spectra.
 !its scattering potential and charge density (P&CD) is a free-atom
 !P&CD.
 use param,only : atdir,outdir,Chg,Pot,idZA,Elem
 use maths,only : trapez,pois,iofx
 use space,only : rx,dx,nx,nxx,nieq
 implicit none
 character :: atm*2
 integer  :: i,ir
 integer,allocatable :: mx(:)
 real(dp) :: HFrmin,HFrmax,r1,a
 real(dp),allocatable :: HFz(:)
 real(dp),parameter :: thrd=1.d0/3.d0,rydb=13.60569172d0
 real(dp),parameter :: pi=acos(-1.d0),fpi=4.d0*pi
 !
 !CHARGE DENSITY data, first run.
 !hf and eeasisss share directory atdir for chgden files.
 write(61,*)
 do ir=1,nieq
   atm=Elem(nint(z(ir)))
   open(50,file=trim(atdir)//'chgden'//trim(atm),status='old')
   read(50,*)
   read(50,*) z(ir),nx(ir),HFrmin,HFrmax
   close(50)
   write(61,898)'z,nx,HFrmin,HFrmax',z(ir),nx(ir),HFrmin,Hfrmax
   r1=HFrmin*(HFrmax/HFrmin)**(1.d0/dble(nx(ir)))
   dx(ir)=(log(HFrmax/HFrmin))/nx(ir)
   write(61,899)'r1,dx ',r1,dx(ir)
 enddo
 !
 !DEFINITION OF ARRAY SPACE from program HF.
 !inside entry: HFz,atrho4pir2
 !goes to eeas: rx.
 nxx=maxval(nx)
 allocate(HFz(nxx),rx(nxx,nieq),atrho4pir2(nxx,nieq),rho(nxx,nieq), &
   rs(nxx,nieq),Vnist(nxx,nieq),mx(nieq))
 do ir=1,nieq
   do i=1,nx(ir)
     rx(i,ir)=r1*exp(dx(ir)*dble(i-1))
   enddo
 enddo
 !
 !CHARGE DENSITY data, second run.
 write(61,'(/a)')"isuper from Shirley's Hartree-Fock program"
 do ir=1,nieq
   atm=Elem(nint(z(ir)))
   open(50,file=trim(atdir)//'chgden'//trim(atm),status='old')
   read(50,*)
   read(50,*) z(ir),nx(ir),HFrmin,HFrmax
   do i=1,nx(ir)
     read(50,*) a,atrho4pir2(i,ir)
   enddo
   close(50)
   !check nuclear charge Z.
   call trapez(rx(1,ir),dx(ir),nx(ir),atrho4pir2(1,ir),1,nx(ir),HFz) ! AMI: Something goes wrong here... Integreation is off by a bit which caused below check to fail - thus removed (chgden files are static anyways)
   write(61,901) ir,z(ir),HFz(nx(ir))
   ! AMI: below checks if integrated charge matches atomic number. This fails for heavy atoms for some reason.
   ! While such a sanity check makes sense in principal, there are two issues:
   ! 1) the atomic density files are not changed in between runs anyways, so it is useless to do this every time
   ! 2) something goes wrong in the integration here (boundaries?). I wrote a short python script that integrates
   !    and checks the values in the charge density files. Differences between actual and integrated Z are less than
   !    0.001% in ALL cases. (see charge_density_check.py in atlib/)
   if(abs(z(ir)-HFz(nx(ir))) > 0.08)then ! Why 0.08 anyways?
   !  write(61,'(a, f9.4)')'Difference in nuclear charge: ', abs(z(ir)-HFz(nx(ir)))
     write(61,'(a)')'stop - Difference in nuclear charge bigger than 0.08' !;stop ! AMI: uncomment to cause stop
   endif
   !avoiding numerical noise of large-radius atrho4pir2.
   do i=nx(ir),nx(ir)/2,-1
     if(atrho4pir2(i,ir) > 1.d-12)then
       mx(ir)=i
       exit
     endif
   enddo

 enddo
 close(10)
 !
 !radial size of "noise-free" atrho4pir2; Ref.: E.Shirley, hf.f90 code. 
 i=minval(mx)
 nx=i
 write(61,'(/a)')'Noise-free atrho4pir2 after E.Shirley.'
 write(61,'(a,20(1x,i0))')'nx =',nx(1:nieq)
 !
 !check atrho4pir2.
 if(Chg=='x')then
   do ir=1,nieq
     open(10,file=trim(outdir)//'/'//'atrho4pir2.'//idZA(ir),status='unknown')
     write(10,902) (rx(i,ir),atrho4pir2(i,ir),i=1,nx(ir))
   enddo
 endif
 !
 !POTENTIAL Vnist.
 do ir=1,nieq
   call pois(rx(1,ir),dx(ir),nx(ir),z(ir),atrho4pir2(1,ir),Vnist(1,ir))
   !do normalization to zero of large-radius "noise-free" Vnist.
   a = Vnist(nx(ir),ir)
   Vnist(1:nx(ir),ir) = Vnist(1:nx(ir),ir) - a
   !enddo normalization.
 enddo
 !check Vnist.
 if(Pot=='y')then
   do ir=1,nieq
     open(10,file=trim(outdir)//'/'//'Vnist.'//idZA(ir),status='unknown')
     do i=nx(ir),1,-1
       a=Vnist(i,ir)*rydb; write(10,902) rx(i,ir),a; if(a<-10.d0) exit
     enddo
     close(10)
   enddo
 endif
 901 format('atom ',i2,' z=',f9.6,' HFz=',f13.10)
 902 format(2es14.6)
 !
 !EXCHANGE-CORRELATION: rs.
 !rs = radius of sphere accommodating a single electron, 
 !rs satisfying (4pi/3)*rs**3 = 1./rho.
 do ir=1,nieq
   do i=1,nx(ir)
     rho(i,ir)=atrho4pir2(i,ir)/(fpi*rx(i,ir)**2)
     rs(i,ir)=( (0.75d0/pi)/rho(i,ir) )**thrd
   enddo
 enddo
 close(50)
 898 format(a,f5.1,i5,2es22.14)
 899 format(a,2es22.14)
 return
 end subroutine Free_Atom_Overlap
end module enrgy
!=========================================================================================================================
!FORTRAN MODULE SDATA with PROGRAM RESHAPE SDATA
 module sdata
 use param,only : dp
! implicit none
!Self-energy data by B.E.Sernelius.
!Ref.: K.W.Shung, B.E.Sernelius, and G.D.Mahan, PRB 36,4499 (1987).
integer,parameter :: nsr=17,nsp=151
real(dp),dimension(nsr),parameter :: sr=(/&
0.00d0,0.01d0,0.02d0,0.05d0,0.10d0,0.20d0,0.30d0,0.40d0,0.50d0,0.70d0,1.00d0,1.50d0,2.00d0,3.00d0,&
4.00d0,5.00d0,6.00d0/)
real(dp),dimension(nsp),parameter :: sp=(/&
0.00d0,0.02d0,0.04d0,0.06d0,0.08d0,0.10d0,0.12d0,0.14d0,0.16d0,0.18d0,0.20d0,0.22d0,0.24d0,0.26d0,&
0.28d0,0.30d0,0.32d0,0.34d0,0.36d0,0.38d0,0.40d0,0.42d0,0.44d0,0.46d0,0.48d0,0.50d0,0.52d0,0.54d0,&
0.56d0,0.58d0,0.60d0,0.62d0,0.64d0,0.66d0,0.68d0,0.70d0,0.72d0,0.74d0,0.76d0,0.78d0,0.80d0,0.82d0,&
0.84d0,0.86d0,0.88d0,0.90d0,0.92d0,0.94d0,0.96d0,0.98d0,1.00d0,1.02d0,1.04d0,1.06d0,1.08d0,1.10d0,&
1.12d0,1.14d0,1.16d0,1.18d0,1.20d0,1.22d0,1.24d0,1.26d0,1.28d0,1.30d0,1.32d0,1.34d0,1.36d0,1.38d0,&
1.40d0,1.42d0,1.44d0,1.46d0,1.48d0,1.50d0,1.52d0,1.54d0,1.56d0,1.58d0,1.60d0,1.62d0,1.64d0,1.66d0,&
1.68d0,1.70d0,1.72d0,1.74d0,1.76d0,1.78d0,1.80d0,1.82d0,1.84d0,1.86d0,1.88d0,1.90d0,1.92d0,1.94d0,&
1.96d0,1.98d0,2.00d0,2.02d0,2.04d0,2.06d0,2.08d0,2.10d0,2.12d0,2.14d0,2.16d0,2.18d0,2.20d0,2.22d0,&
2.24d0,2.26d0,2.28d0,2.30d0,2.32d0,2.34d0,2.36d0,2.38d0,2.40d0,2.42d0,2.44d0,2.46d0,2.48d0,2.50d0,&
2.52d0,2.54d0,2.56d0,2.58d0,2.60d0,2.62d0,2.64d0,2.66d0,2.68d0,2.70d0,2.72d0,2.74d0,2.76d0,2.78d0,&
2.80d0,2.82d0,2.84d0,2.86d0,2.88d0,2.90d0,2.92d0,2.94d0,2.96d0,2.98d0,3.00d0/)
!
real(dp),dimension(nsp,nsr) :: sdat
real(dp),dimension(nsp*nsr),parameter :: sdat1=(/&
-0.66340d0,-0.66330d0,-0.66310d0,-0.66260d0,-0.66200d0,-0.66120d0,-0.66020d0,-0.65910d0,-0.65770d0,-0.65620d0,-0.65450d0,&
-0.65260d0,-0.65050d0,-0.64830d0,-0.64580d0,-0.64320d0,-0.64030d0,-0.63720d0,-0.63400d0,-0.63050d0,-0.62680d0,-0.62290d0,&
-0.61880d0,-0.61450d0,-0.60990d0,-0.60500d0,-0.60000d0,-0.59460d0,-0.58900d0,-0.58310d0,-0.57700d0,-0.57050d0,-0.56370d0,&
-0.55660d0,-0.54920d0,-0.54130d0,-0.53310d0,-0.52450d0,-0.51540d0,-0.50580d0,-0.49570d0,-0.48500d0,-0.47370d0,-0.46160d0,&
-0.44870d0,-0.43480d0,-0.41970d0,-0.40310d0,-0.38440d0,-0.36250d0,-0.33170d0,-0.30140d0,-0.28060d0,-0.26330d0,-0.24850d0,&
-0.23530d0,-0.22350d0,-0.21290d0,-0.20310d0,-0.19420d0,-0.18590d0,-0.17820d0,-0.17110d0,-0.16450d0,-0.15820d0,-0.15240d0,&
-0.14690d0,-0.14180d0,-0.13690d0,-0.13230d0,-0.12790d0,-0.12380d0,-0.11990d0,-0.11620d0,-0.11260d0,-0.10930d0,-0.10610d0,&
-0.10300d0,-0.10010d0,-0.09727d0,-0.09459d0,-0.09203d0,-0.08958d0,-0.08722d0,-0.08496d0,-0.08280d0,-0.08071d0,-0.07870d0,&
-0.07678d0,-0.07493d0,-0.07314d0,-0.07142d0,-0.06977d0,-0.06817d0,-0.06662d0,-0.06513d0,-0.06369d0,-0.06230d0,-0.06096d0,&
-0.05965d0,-0.05840d0,-0.05718d0,-0.05600d0,-0.05486d0,-0.05375d0,-0.05267d0,-0.05163d0,-0.05062d0,-0.04964d0,-0.04869d0,&
-0.04777d0,-0.04687d0,-0.04600d0,-0.04515d0,-0.04433d0,-0.04353d0,-0.04275d0,-0.04199d0,-0.04125d0,-0.04054d0,-0.03984d0,&
-0.03916d0,-0.03849d0,-0.03785d0,-0.03722d0,-0.03660d0,-0.03600d0,-0.03542d0,-0.03485d0,-0.03429d0,-0.03375d0,-0.03322d0,&
-0.03270d0,-0.03220d0,-0.03170d0,-0.03122d0,-0.03075d0,-0.03029d0,-0.02984d0,-0.02940d0,-0.02897d0,-0.02855d0,-0.02814d0,&
-0.02773d0,-0.02734d0,-0.02696d0,-0.02658d0,-0.02621d0,-0.02585d0,-0.02549d0,-0.02515d0,-0.61980d0,-0.61970d0,-0.61950d0,&
-0.61920d0,-0.61870d0,-0.61800d0,-0.61720d0,-0.61620d0,-0.61500d0,-0.61360d0,-0.61210d0,-0.61040d0,-0.60840d0,-0.60630d0,&
-0.60400d0,-0.60150d0,-0.59880d0,-0.59590d0,-0.59280d0,-0.58950d0,-0.58600d0,-0.58220d0,-0.57830d0,-0.57410d0,-0.56970d0,&
-0.56500d0,-0.56010d0,-0.55500d0,-0.54960d0,-0.54390d0,-0.53790d0,-0.53170d0,-0.52510d0,-0.51820d0,-0.51100d0,-0.50340d0,&
-0.49550d0,-0.48710d0,-0.47830d0,-0.46910d0,-0.45940d0,-0.44910d0,-0.43830d0,-0.42690d0,-0.41490d0,-0.40230d0,-0.38920d0,&
-0.37560d0,-0.36180d0,-0.34750d0,-0.33290d0,-0.31990d0,-0.30740d0,-0.29630d0,-0.28700d0,-0.27060d0,-0.25900d0,-0.24830d0,&
-0.23830d0,-0.22900d0,-0.22050d0,-0.21250d0,-0.20500d0,-0.19800d0,-0.19150d0,-0.18530d0,-0.17950d0,-0.17400d0,-0.16880d0,&
-0.16390d0,-0.15930d0,-0.15490d0,-0.15060d0,-0.14660d0,-0.14280d0,-0.13920d0,-0.13570d0,-0.13240d0,-0.12920d0,-0.12620d0,&
-0.12330d0,-0.12050d0,-0.11780d0,-0.11520d0,-0.11270d0,-0.11030d0,-0.10800d0,-0.10580d0,-0.10370d0,-0.10170d0,-0.09967d0,&
-0.09776d0,-0.09591d0,-0.09412d0,-0.09240d0,-0.09073d0,-0.08911d0,-0.08754d0,-0.08603d0,-0.08456d0,-0.08314d0,-0.08176d0,&
-0.08042d0,-0.07912d0,-0.07786d0,-0.07663d0,-0.07544d0,-0.07428d0,-0.07316d0,-0.07206d0,-0.07100d0,-0.06996d0,-0.06895d0,&
-0.06797d0,-0.06701d0,-0.06608d0,-0.06517d0,-0.06429d0,-0.06343d0,-0.06259d0,-0.06177d0,-0.06097d0,-0.06018d0,-0.05942d0,&
-0.05860d0,-0.05794d0,-0.05723d0,-0.05653d0,-0.05585d0,-0.05518d0,-0.05453d0,-0.05389d0,-0.05327d0,-0.05266d0,-0.05206d0,&
-0.05147d0,-0.05090d0,-0.05034d0,-0.04979d0,-0.04925d0,-0.04873d0,-0.04821d0,-0.04773d0,-0.04723d0,-0.04675d0,-0.04628d0,&
-0.04581d0,-0.04536d0,-0.04493d0,-0.04450d0,-0.04407d0,-0.60270d0,-0.60270d0,-0.60250d0,-0.60210d0,-0.60170d0,-0.60110d0,&
-0.60030d0,-0.59930d0,-0.59820d0,-0.59690d0,-0.59540d0,-0.59380d0,-0.59190d0,-0.58990d0,-0.58770d0,-0.58530d0,-0.58270d0,&
-0.57980d0,-0.57680d0,-0.57360d0,-0.57020d0,-0.56650d0,-0.56270d0,-0.55860d0,-0.55430d0,-0.54970d0,-0.54490d0,-0.53990d0,&
-0.53450d0,-0.52900d0,-0.52310d0,-0.51700d0,-0.51050d0,-0.50380d0,-0.49670d0,-0.48930d0,-0.48150d0,-0.47330d0,-0.46480d0,&
-0.45580d0,-0.44640d0,-0.43660d0,-0.42630d0,-0.41570d0,-0.40470d0,-0.39340d0,-0.38190d0,-0.37030d0,-0.35840d0,-0.34630d0,&
-0.33390d0,-0.32240d0,-0.31130d0,-0.30100d0,-0.29130d0,-0.28190d0,-0.27220d0,-0.26220d0,-0.25250d0,-0.24330d0,-0.23470d0,&
-0.22660d0,-0.21900d0,-0.21190d0,-0.20520d0,-0.19880d0,-0.19290d0,-0.18730d0,-0.18190d0,-0.17690d0,-0.17210d0,-0.16750d0,&
-0.16320d0,-0.15900d0,-0.15510d0,-0.15130d0,-0.14770d0,-0.14430d0,-0.14100d0,-0.13790d0,-0.13480d0,-0.13190d0,-0.12910d0,&
-0.12650d0,-0.12390d0,-0.12140d0,-0.11900d0,-0.11670d0,-0.11450d0,-0.11230d0,-0.11030d0,-0.10830d0,-0.10630d0,-0.10450d0,&
-0.10270d0,-0.10090d0,-0.09924d0,-0.09760d0,-0.09601d0,-0.09440d0,-0.09297d0,-0.09150d0,-0.09012d0,-0.08875d0,-0.08742d0,&
-0.08613d0,-0.08488d0,-0.08366d0,-0.08247d0,-0.08132d0,-0.08019d0,-0.07910d0,-0.07803d0,-0.07699d0,-0.07598d0,-0.07499d0,&
-0.07403d0,-0.07309d0,-0.07218d0,-0.07128d0,-0.07041d0,-0.06956d0,-0.06872d0,-0.06791d0,-0.06712d0,-0.06634d0,-0.06558d0,&
-0.06484d0,-0.06411d0,-0.06340d0,-0.06271d0,-0.06203d0,-0.06136d0,-0.06071d0,-0.06008d0,-0.05945d0,-0.05884d0,-0.05824d0,&
-0.05765d0,-0.05707d0,-0.05651d0,-0.05596d0,-0.05541d0,-0.05488d0,-0.05435d0,-0.05384d0,-0.05333d0,-0.05284d0,-0.05235d0,&
-0.05187d0,-0.05141d0,-0.57140d0,-0.57130d0,-0.57110d0,-0.57080d0,-0.57040d0,-0.56980d0,-0.56910d0,-0.56820d0,-0.56720d0,&
-0.56600d0,-0.56470d0,-0.56310d0,-0.56140d0,-0.55950d0,-0.55750d0,-0.55520d0,-0.55280d0,-0.55010d0,-0.54730d0,-0.54430d0,&
-0.54100d0,-0.53760d0,-0.53390d0,-0.53000d0,-0.52590d0,-0.52160d0,-0.51700d0,-0.51220d0,-0.50720d0,-0.50190d0,-0.49630d0,&
-0.49050d0,-0.48440d0,-0.47800d0,-0.47130d0,-0.46440d0,-0.45710d0,-0.44960d0,-0.44180d0,-0.43380d0,-0.42560d0,-0.41710d0,&
-0.40860d0,-0.39990d0,-0.39120d0,-0.38240d0,-0.37350d0,-0.36440d0,-0.35530d0,-0.34590d0,-0.33650d0,-0.32740d0,-0.31860d0,&
-0.31020d0,-0.30230d0,-0.29460d0,-0.28780d0,-0.28130d0,-0.27630d0,-0.26990d0,-0.26200d0,-0.25410d0,-0.24650d0,-0.23920d0,&
-0.23230d0,-0.22580d0,-0.21950d0,-0.21360d0,-0.20800d0,-0.20260d0,-0.19750d0,-0.19270d0,-0.18800d0,-0.18360d0,-0.17940d0,&
-0.17540d0,-0.17150d0,-0.16780d0,-0.16430d0,-0.16090d0,-0.15760d0,-0.15450d0,-0.15150d0,-0.14860d0,-0.14580d0,-0.14320d0,&
-0.14060d0,-0.13810d0,-0.13570d0,-0.13330d0,-0.13110d0,-0.12890d0,-0.12680d0,-0.12480d0,-0.12280d0,-0.12090d0,-0.11910d0,&
-0.11730d0,-0.11550d0,-0.11390d0,-0.11220d0,-0.11060d0,-0.10910d0,-0.10760d0,-0.10610d0,-0.10470d0,-0.10330d0,-0.10200d0,&
-0.10070d0,-0.09939d0,-0.09815d0,-0.09694d0,-0.09575d0,-0.09460d0,-0.09348d0,-0.09238d0,-0.09131d0,-0.09026d0,-0.08924d0,&
-0.08825d0,-0.08727d0,-0.08632d0,-0.08539d0,-0.08448d0,-0.08359d0,-0.08272d0,-0.08186d0,-0.08103d0,-0.08021d0,-0.07941d0,&
-0.07863d0,-0.07786d0,-0.07711d0,-0.07638d0,-0.07566d0,-0.07495d0,-0.07425d0,-0.07357d0,-0.07290d0,-0.07225d0,-0.07161d0,&
-0.07098d0,-0.07036d0,-0.06975d0,-0.06915d0,-0.06857d0,-0.06799d0,-0.06743d0,-0.06687d0,-0.06633d0,-0.06579d0,-0.54000d0,&
-0.54000d0,-0.53980d0,-0.53950d0,-0.53910d0,-0.53850d0,-0.53780d0,-0.53700d0,-0.53600d0,-0.53490d0,-0.53370d0,-0.53230d0,&
-0.53070d0,-0.52890d0,-0.52700d0,-0.52490d0,-0.52270d0,-0.52020d0,-0.51760d0,-0.51480d0,-0.51170d0,-0.50850d0,-0.50510d0,&
-0.50150d0,-0.49770d0,-0.49370d0,-0.48940d0,-0.48500d0,-0.48030d0,-0.47540d0,-0.47030d0,-0.46500d0,-0.45950d0,-0.45380d0,&
-0.44790d0,-0.44180d0,-0.43560d0,-0.42920d0,-0.42270d0,-0.41620d0,-0.40950d0,-0.40290d0,-0.39620d0,-0.38940d0,-0.38260d0,&
-0.37580d0,-0.36880d0,-0.36180d0,-0.35470d0,-0.34740d0,-0.34010d0,-0.33300d0,-0.32610d0,-0.31940d0,-0.31310d0,-0.30700d0,&
-0.30140d0,-0.29600d0,-0.29120d0,-0.28650d0,-0.28240d0,-0.28030d0,-0.27580d0,-0.26920d0,-0.26250d0,-0.25590d0,-0.24950d0,&
-0.24330d0,-0.23740d0,-0.23170d0,-0.22630d0,-0.22110d0,-0.21610d0,-0.21140d0,-0.20680d0,-0.20250d0,-0.19830d0,-0.19430d0,&
-0.19050d0,-0.18680d0,-0.18330d0,-0.17990d0,-0.17660d0,-0.17350d0,-0.17050d0,-0.16750d0,-0.16470d0,-0.16200d0,-0.15940d0,&
-0.15680d0,-0.15440d0,-0.15200d0,-0.14970d0,-0.14750d0,-0.14530d0,-0.14320d0,-0.14120d0,-0.13920d0,-0.13730d0,-0.13550d0,&
-0.13370d0,-0.13190d0,-0.13020d0,-0.12850d0,-0.12690d0,-0.12530d0,-0.12380d0,-0.12230d0,-0.12090d0,-0.11950d0,-0.11810d0,&
-0.11670d0,-0.11540d0,-0.11410d0,-0.11290d0,-0.11170d0,-0.11050d0,-0.10930d0,-0.10820d0,-0.10700d0,-0.10590d0,-0.10490d0,&
-0.10380d0,-0.10280d0,-0.10180d0,-0.10080d0,-0.09988d0,-0.09894d0,-0.09802d0,-0.09712d0,-0.09623d0,-0.09537d0,-0.09452d0,&
-0.09369d0,-0.09287d0,-0.09207d0,-0.09129d0,-0.09051d0,-0.08976d0,-0.08902d0,-0.08829d0,-0.08757d0,-0.08687d0,-0.08618d0,&
-0.08550d0,-0.08483d0,-0.08417d0,-0.08353d0,-0.08289d0,-0.08227d0,-0.08166d0,-0.50230d0,-0.50220d0,-0.50210d0,-0.50180d0,&
-0.50140d0,-0.50090d0,-0.50020d0,-0.49940d0,-0.49860d0,-0.49750d0,-0.49640d0,-0.49510d0,-0.49370d0,-0.49210d0,-0.49040d0,&
-0.48860d0,-0.48650d0,-0.48440d0,-0.48200d0,-0.47960d0,-0.47690d0,-0.47410d0,-0.47110d0,-0.46800d0,-0.46470d0,-0.46120d0,&
-0.45760d0,-0.45390d0,-0.44990d0,-0.44590d0,-0.44170d0,-0.43750d0,-0.43310d0,-0.42860d0,-0.42410d0,-0.41950d0,-0.41490d0,&
-0.41020d0,-0.40560d0,-0.40090d0,-0.39610d0,-0.39140d0,-0.38660d0,-0.38180d0,-0.37690d0,-0.37200d0,-0.36700d0,-0.36190d0,&
-0.35670d0,-0.35150d0,-0.34620d0,-0.34100d0,-0.33600d0,-0.33110d0,-0.32640d0,-0.32190d0,-0.31770d0,-0.31370d0,-0.30990d0,&
-0.30640d0,-0.30310d0,-0.30040d0,-0.29800d0,-0.29590d0,-0.29480d0,-0.29560d0,-0.28990d0,-0.28420d0,-0.27840d0,-0.27260d0,&
-0.26680d0,-0.26130d0,-0.25590d0,-0.25070d0,-0.24570d0,-0.24090d0,-0.23630d0,-0.23190d0,-0.22760d0,-0.22350d0,-0.21960d0,&
-0.21580d0,-0.21210d0,-0.20860d0,-0.20520d0,-0.20190d0,-0.19880d0,-0.19570d0,-0.19270d0,-0.18990d0,-0.18710d0,-0.18440d0,&
-0.18180d0,-0.17930d0,-0.17690d0,-0.17450d0,-0.17220d0,-0.17000d0,-0.16780d0,-0.16570d0,-0.16370d0,-0.16170d0,-0.15970d0,&
-0.15780d0,-0.15600d0,-0.15420d0,-0.15250d0,-0.15080d0,-0.14910d0,-0.14750d0,-0.14590d0,-0.14440d0,-0.14280d0,-0.14140d0,&
-0.13990d0,-0.13850d0,-0.13720d0,-0.13580d0,-0.13450d0,-0.13320d0,-0.13200d0,-0.13070d0,-0.12950d0,-0.12830d0,-0.12720d0,&
-0.12600d0,-0.12490d0,-0.12380d0,-0.12280d0,-0.12170d0,-0.12070d0,-0.11970d0,-0.11870d0,-0.11770d0,-0.11680d0,-0.11590d0,&
-0.11490d0,-0.11400d0,-0.11320d0,-0.11230d0,-0.11140d0,-0.11060d0,-0.10980d0,-0.10900d0,-0.10820d0,-0.10740d0,-0.10660d0,&
-0.10590d0,-0.10510d0,-0.10440d0,-0.10370d0,-0.47820d0,-0.47810d0,-0.47800d0,-0.47770d0,-0.47730d0,-0.47680d0,-0.47620d0,&
-0.47550d0,-0.47470d0,-0.47370d0,-0.47270d0,-0.47150d0,-0.47020d0,-0.46880d0,-0.46720d0,-0.46560d0,-0.46380d0,-0.46180d0,&
-0.45980d0,-0.45760d0,-0.45530d0,-0.45280d0,-0.45030d0,-0.44760d0,-0.44480d0,-0.44190d0,-0.43880d0,-0.43570d0,-0.43250d0,&
-0.42920d0,-0.42590d0,-0.42250d0,-0.41900d0,-0.41550d0,-0.41200d0,-0.40850d0,-0.40490d0,-0.40130d0,-0.39780d0,-0.39410d0,&
-0.39050d0,-0.38690d0,-0.38320d0,-0.37940d0,-0.37560d0,-0.37180d0,-0.36780d0,-0.36380d0,-0.35980d0,-0.35560d0,-0.35150d0,&
-0.34740d0,-0.34340d0,-0.33950d0,-0.33580d0,-0.33220d0,-0.32880d0,-0.32560d0,-0.32260d0,-0.31970d0,-0.31700d0,-0.31480d0,&
-0.31270d0,-0.31080d0,-0.30940d0,-0.30850d0,-0.30780d0,-0.30900d0,-0.30780d0,-0.30260d0,-0.29710d0,-0.29160d0,-0.28610d0,&
-0.28070d0,-0.27540d0,-0.27030d0,-0.26540d0,-0.26060d0,-0.25610d0,-0.25160d0,-0.24740d0,-0.24330d0,-0.23930d0,-0.23550d0,&
-0.23180d0,-0.22830d0,-0.22480d0,-0.22150d0,-0.21830d0,-0.21520d0,-0.21210d0,-0.20920d0,-0.20640d0,-0.20370d0,-0.20100d0,&
-0.19840d0,-0.19590d0,-0.19350d0,-0.19110d0,-0.18880d0,-0.18650d0,-0.18440d0,-0.18220d0,-0.18020d0,-0.17810d0,-0.17620d0,&
-0.17430d0,-0.17240d0,-0.17060d0,-0.16880d0,-0.16700d0,-0.16530d0,-0.16370d0,-0.16210d0,-0.16050d0,-0.15890d0,-0.15740d0,&
-0.15590d0,-0.15450d0,-0.15310d0,-0.15170d0,-0.15030d0,-0.14900d0,-0.14770d0,-0.14640d0,-0.14520d0,-0.14390d0,-0.14270d0,&
-0.14150d0,-0.14040d0,-0.13920d0,-0.13810d0,-0.13700d0,-0.13600d0,-0.13490d0,-0.13390d0,-0.13280d0,-0.13180d0,-0.13090d0,&
-0.12990d0,-0.12890d0,-0.12800d0,-0.12710d0,-0.12620d0,-0.12530d0,-0.12440d0,-0.12360d0,-0.12270d0,-0.12190d0,-0.12110d0,&
-0.12030d0,-0.46080d0,-0.46080d0,-0.46060d0,-0.46040d0,-0.46000d0,-0.45960d0,-0.45900d0,-0.45830d0,-0.45760d0,-0.45670d0,&
-0.45580d0,-0.45470d0,-0.45350d0,-0.45220d0,-0.45090d0,-0.44940d0,-0.44780d0,-0.44610d0,-0.44430d0,-0.44240d0,-0.44040d0,&
-0.43830d0,-0.43610d0,-0.43380d0,-0.43140d0,-0.42900d0,-0.42650d0,-0.42390d0,-0.42130d0,-0.41870d0,-0.41600d0,-0.41320d0,&
-0.41050d0,-0.40770d0,-0.40490d0,-0.40210d0,-0.39930d0,-0.39650d0,-0.39370d0,-0.39080d0,-0.38790d0,-0.38500d0,-0.38200d0,&
-0.37900d0,-0.37590d0,-0.37280d0,-0.36960d0,-0.36640d0,-0.36300d0,-0.35970d0,-0.35620d0,-0.35290d0,-0.34960d0,-0.34640d0,&
-0.34330d0,-0.34040d0,-0.33760d0,-0.33490d0,-0.33240d0,-0.33000d0,-0.32780d0,-0.32590d0,-0.32420d0,-0.32260d0,-0.32130d0,&
-0.32040d0,-0.31960d0,-0.31940d0,-0.31990d0,-0.32140d0,-0.32150d0,-0.31610d0,-0.31090d0,-0.30550d0,-0.30020d0,-0.29490d0,&
-0.28980d0,-0.28480d0,-0.28000d0,-0.27530d0,-0.27080d0,-0.26640d0,-0.26220d0,-0.25820d0,-0.25420d0,-0.25040d0,-0.24680d0,&
-0.24320d0,-0.23980d0,-0.23640d0,-0.23320d0,-0.23010d0,-0.22700d0,-0.22410d0,-0.22130d0,-0.21850d0,-0.21580d0,-0.21320d0,&
-0.21060d0,-0.20810d0,-0.20570d0,-0.20340d0,-0.20110d0,-0.19890d0,-0.19670d0,-0.19460d0,-0.19250d0,-0.19050d0,-0.18860d0,&
-0.18660d0,-0.18480d0,-0.18290d0,-0.18110d0,-0.17940d0,-0.17770d0,-0.17600d0,-0.17440d0,-0.17280d0,-0.17120d0,-0.16970d0,&
-0.16820d0,-0.16670d0,-0.16530d0,-0.16390d0,-0.16250d0,-0.16110d0,-0.15980d0,-0.15850d0,-0.15720d0,-0.15600d0,-0.15470d0,&
-0.15350d0,-0.15230d0,-0.15120d0,-0.15000d0,-0.14890d0,-0.14780d0,-0.14670d0,-0.14560d0,-0.14460d0,-0.14360d0,-0.14260d0,&
-0.14160d0,-0.14060d0,-0.13960d0,-0.13870d0,-0.13770d0,-0.13680d0,-0.13590d0,-0.13500d0,-0.13410d0,-0.44770d0,-0.44760d0,&
-0.44750d0,-0.44730d0,-0.44690d0,-0.44650d0,-0.44600d0,-0.44540d0,-0.44470d0,-0.44390d0,-0.44310d0,-0.44210d0,-0.44100d0,&
-0.43990d0,-0.43870d0,-0.43740d0,-0.43600d0,-0.43450d0,-0.43290d0,-0.43130d0,-0.42960d0,-0.42780d0,-0.42590d0,-0.42400d0,&
-0.42210d0,-0.42000d0,-0.41800d0,-0.41590d0,-0.41380d0,-0.41160d0,-0.40940d0,-0.40720d0,-0.40500d0,-0.40280d0,-0.40060d0,&
-0.39830d0,-0.39610d0,-0.39380d0,-0.39150d0,-0.38920d0,-0.38690d0,-0.38450d0,-0.38210d0,-0.37960d0,-0.37700d0,-0.37450d0,&
-0.37180d0,-0.36910d0,-0.36630d0,-0.36350d0,-0.36060d0,-0.35780d0,-0.35510d0,-0.35240d0,-0.34980d0,-0.34730d0,-0.34500d0,&
-0.34270d0,-0.34060d0,-0.33860d0,-0.33680d0,-0.33510d0,-0.33370d0,-0.33230d0,-0.33130d0,-0.33050d0,-0.32970d0,-0.32950d0,&
-0.32960d0,-0.32980d0,-0.33100d0,-0.33440d0,-0.33210d0,-0.32670d0,-0.32150d0,-0.31630d0,-0.31110d0,-0.30600d0,-0.30100d0,&
-0.29610d0,-0.29140d0,-0.28680d0,-0.28240d0,-0.27810d0,-0.27400d0,-0.27000d0,-0.26610d0,-0.26230d0,-0.25870d0,-0.25520d0,&
-0.25170d0,-0.24840d0,-0.24520d0,-0.24210d0,-0.23910d0,-0.23610d0,-0.23330d0,-0.23050d0,-0.22780d0,-0.22520d0,-0.22260d0,&
-0.22010d0,-0.21770d0,-0.21530d0,-0.21300d0,-0.21080d0,-0.20860d0,-0.20640d0,-0.20430d0,-0.20230d0,-0.20030d0,-0.19840d0,&
-0.19650d0,-0.19460d0,-0.19280d0,-0.19100d0,-0.18930d0,-0.18760d0,-0.18590d0,-0.18430d0,-0.18270d0,-0.18110d0,-0.17950d0,&
-0.17800d0,-0.17660d0,-0.17510d0,-0.17370d0,-0.17230d0,-0.17090d0,-0.16960d0,-0.16830d0,-0.16700d0,-0.16570d0,-0.16450d0,&
-0.16330d0,-0.16210d0,-0.16090d0,-0.15970d0,-0.15860d0,-0.15750d0,-0.15640d0,-0.15530d0,-0.15420d0,-0.15320d0,-0.15210d0,&
-0.15110d0,-0.15010d0,-0.14910d0,-0.14810d0,-0.14720d0,-0.14630d0,-0.42930d0,-0.42920d0,-0.42910d0,-0.42890d0,-0.42870d0,&
-0.42830d0,-0.42790d0,-0.42740d0,-0.42690d0,-0.42620d0,-0.42560d0,-0.42480d0,-0.42400d0,-0.42310d0,-0.42220d0,-0.42120d0,&
-0.42010d0,-0.41910d0,-0.41790d0,-0.41670d0,-0.41550d0,-0.41430d0,-0.41300d0,-0.41170d0,-0.41030d0,-0.40900d0,-0.40760d0,&
-0.40620d0,-0.40480d0,-0.40340d0,-0.40200d0,-0.40060d0,-0.39920d0,-0.39770d0,-0.39630d0,-0.39480d0,-0.39330d0,-0.39180d0,&
-0.39030d0,-0.38880d0,-0.38720d0,-0.38550d0,-0.38390d0,-0.38220d0,-0.38040d0,-0.37860d0,-0.37670d0,-0.37470d0,-0.37270d0,&
-0.37070d0,-0.36860d0,-0.36660d0,-0.36460d0,-0.36260d0,-0.36070d0,-0.35890d0,-0.35720d0,-0.35560d0,-0.35400d0,-0.35260d0,&
-0.35130d0,-0.35010d0,-0.34910d0,-0.34810d0,-0.34740d0,-0.34680d0,-0.34630d0,-0.34620d0,-0.34620d0,-0.34630d0,-0.34670d0,&
-0.34740d0,-0.34850d0,-0.35060d0,-0.35660d0,-0.35290d0,-0.34720d0,-0.34200d0,-0.33690d0,-0.33180d0,-0.32680d0,-0.32200d0,&
-0.31720d0,-0.31260d0,-0.30810d0,-0.30370d0,-0.29950d0,-0.29540d0,-0.29140d0,-0.28760d0,-0.28380d0,-0.28020d0,-0.27670d0,&
-0.27320d0,-0.26990d0,-0.26670d0,-0.26350d0,-0.26050d0,-0.25750d0,-0.25460d0,-0.25180d0,-0.24900d0,-0.24630d0,-0.24370d0,&
-0.24120d0,-0.23870d0,-0.23630d0,-0.23390d0,-0.23160d0,-0.22930d0,-0.22710d0,-0.22500d0,-0.22290d0,-0.22080d0,-0.21880d0,&
-0.21680d0,-0.21490d0,-0.21300d0,-0.21120d0,-0.20940d0,-0.20760d0,-0.20590d0,-0.20410d0,-0.20250d0,-0.20080d0,-0.19920d0,&
-0.19770d0,-0.19610d0,-0.19460d0,-0.19310d0,-0.19160d0,-0.19020d0,-0.18880d0,-0.18740d0,-0.18600d0,-0.18470d0,-0.18340d0,&
-0.18210d0,-0.18080d0,-0.17960d0,-0.17840d0,-0.17710d0,-0.17600d0,-0.17480d0,-0.17360d0,-0.17250d0,-0.17140d0,-0.17030d0,&
-0.16920d0,-0.16810d0,-0.16710d0,-0.41320d0,-0.41320d0,-0.41310d0,-0.41300d0,-0.41280d0,-0.41260d0,-0.41230d0,-0.41200d0,&
-0.41160d0,-0.41120d0,-0.41080d0,-0.41030d0,-0.40980d0,-0.40920d0,-0.40870d0,-0.40810d0,-0.40750d0,-0.40680d0,-0.40620d0,&
-0.40550d0,-0.40490d0,-0.40420d0,-0.40350d0,-0.40280d0,-0.40210d0,-0.40150d0,-0.40080d0,-0.40010d0,-0.39940d0,-0.39870d0,&
-0.39800d0,-0.39730d0,-0.39660d0,-0.39590d0,-0.39520d0,-0.39450d0,-0.39370d0,-0.39290d0,-0.39210d0,-0.39130d0,-0.39040d0,&
-0.38950d0,-0.38850d0,-0.38750d0,-0.38650d0,-0.38540d0,-0.38420d0,-0.38300d0,-0.38180d0,-0.38050d0,-0.37920d0,-0.37790d0,&
-0.37660d0,-0.37530d0,-0.37410d0,-0.37290d0,-0.37190d0,-0.37080d0,-0.36990d0,-0.36900d0,-0.36820d0,-0.36750d0,-0.36690d0,&
-0.36640d0,-0.36600d0,-0.36580d0,-0.36560d0,-0.36560d0,-0.36580d0,-0.36600d0,-0.36650d0,-0.36700d0,-0.36770d0,-0.36860d0,&
-0.36980d0,-0.37140d0,-0.37390d0,-0.37810d0,-0.38520d0,-0.37700d0,-0.37110d0,-0.36580d0,-0.36070d0,-0.35570d0,-0.35090d0,&
-0.34610d0,-0.34150d0,-0.33700d0,-0.33260d0,-0.32830d0,-0.32410d0,-0.32010d0,-0.31620d0,-0.31240d0,-0.30860d0,-0.30500d0,&
-0.30150d0,-0.29810d0,-0.29480d0,-0.29150d0,-0.28840d0,-0.28530d0,-0.28230d0,-0.27940d0,-0.27650d0,-0.27370d0,-0.27100d0,&
-0.26840d0,-0.26580d0,-0.26320d0,-0.26080d0,-0.25830d0,-0.25600d0,-0.25370d0,-0.25140d0,-0.24920d0,-0.24700d0,-0.24490d0,&
-0.24280d0,-0.24080d0,-0.23880d0,-0.23680d0,-0.23490d0,-0.23300d0,-0.23120d0,-0.22940d0,-0.22760d0,-0.22580d0,-0.22410d0,&
-0.22240d0,-0.22080d0,-0.21920d0,-0.21760d0,-0.21600d0,-0.21450d0,-0.21300d0,-0.21150d0,-0.21000d0,-0.20860d0,-0.20720d0,&
-0.20580d0,-0.20440d0,-0.20310d0,-0.20180d0,-0.20040d0,-0.19920d0,-0.19790d0,-0.19670d0,-0.19540d0,-0.19420d0,-0.19300d0,&
-0.40190d0,-0.40190d0,-0.40190d0,-0.40180d0,-0.40180d0,-0.40170d0,-0.40160d0,-0.40150d0,-0.40130d0,-0.40120d0,-0.40110d0,&
-0.40090d0,-0.40080d0,-0.40060d0,-0.40050d0,-0.40040d0,-0.40020d0,-0.40010d0,-0.40000d0,-0.39990d0,-0.39980d0,-0.39970d0,&
-0.39960d0,-0.39960d0,-0.39950d0,-0.39950d0,-0.39940d0,-0.39940d0,-0.39940d0,-0.39940d0,-0.39930d0,-0.39930d0,-0.39930d0,&
-0.39920d0,-0.39920d0,-0.39910d0,-0.39900d0,-0.39890d0,-0.39880d0,-0.39860d0,-0.39840d0,-0.39820d0,-0.39790d0,-0.39760d0,&
-0.39720d0,-0.39680d0,-0.39640d0,-0.39590d0,-0.39540d0,-0.39480d0,-0.39420d0,-0.39360d0,-0.39310d0,-0.39250d0,-0.39200d0,&
-0.39150d0,-0.39100d0,-0.39060d0,-0.39020d0,-0.38990d0,-0.38970d0,-0.38950d0,-0.38940d0,-0.38930d0,-0.38930d0,-0.38950d0,&
-0.38960d0,-0.39000d0,-0.39040d0,-0.39090d0,-0.39150d0,-0.39210d0,-0.39280d0,-0.39360d0,-0.39460d0,-0.39560d0,-0.39680d0,&
-0.39820d0,-0.40000d0,-0.40220d0,-0.40520d0,-0.40980d0,-0.42080d0,-0.41930d0,-0.41160d0,-0.40520d0,-0.39960d0,-0.39430d0,&
-0.38920d0,-0.38430d0,-0.37950d0,-0.37490d0,-0.37040d0,-0.36610d0,-0.36180d0,-0.35770d0,-0.35360d0,-0.34970d0,-0.34590d0,&
-0.34210d0,-0.33850d0,-0.33500d0,-0.33150d0,-0.32810d0,-0.32480d0,-0.32160d0,-0.31850d0,-0.31540d0,-0.31240d0,-0.30950d0,&
-0.30660d0,-0.30390d0,-0.30110d0,-0.29840d0,-0.29580d0,-0.29330d0,-0.29080d0,-0.28830d0,-0.28590d0,-0.28350d0,-0.28120d0,&
-0.27900d0,-0.27670d0,-0.27460d0,-0.27240d0,-0.27030d0,-0.26830d0,-0.26630d0,-0.26430d0,-0.26230d0,-0.26040d0,-0.25850d0,&
-0.25670d0,-0.25490d0,-0.25310d0,-0.25130d0,-0.24960d0,-0.24790d0,-0.24630d0,-0.24460d0,-0.24300d0,-0.24140d0,-0.23990d0,&
-0.23830d0,-0.23680d0,-0.23530d0,-0.23380d0,-0.23240d0,-0.23100d0,-0.22960d0,-0.22820d0,-0.39990d0,-0.39990d0,-0.39990d0,&
-0.39990d0,-0.39990d0,-0.39990d0,-0.39990d0,-0.40000d0,-0.40000d0,-0.40000d0,-0.40010d0,-0.40020d0,-0.40020d0,-0.40030d0,&
-0.40040d0,-0.40060d0,-0.40070d0,-0.40090d0,-0.40110d0,-0.40130d0,-0.40150d0,-0.40170d0,-0.40200d0,-0.40230d0,-0.40260d0,&
-0.40290d0,-0.40320d0,-0.40350d0,-0.40380d0,-0.40420d0,-0.40450d0,-0.40480d0,-0.40520d0,-0.40550d0,-0.40580d0,-0.40610d0,&
-0.40640d0,-0.40660d0,-0.40690d0,-0.40710d0,-0.40720d0,-0.40740d0,-0.40750d0,-0.40760d0,-0.40770d0,-0.40770d0,-0.40770d0,&
-0.40760d0,-0.40750d0,-0.40740d0,-0.40720d0,-0.40710d0,-0.40690d0,-0.40680d0,-0.40670d0,-0.40660d0,-0.40650d0,-0.40650d0,&
-0.40650d0,-0.40650d0,-0.40660d0,-0.40670d0,-0.40690d0,-0.40720d0,-0.40750d0,-0.40780d0,-0.40820d0,-0.40880d0,-0.40940d0,&
-0.41000d0,-0.41080d0,-0.41160d0,-0.41240d0,-0.41330d0,-0.41420d0,-0.41520d0,-0.41630d0,-0.41750d0,-0.41890d0,-0.42050d0,&
-0.42220d0,-0.42430d0,-0.42680d0,-0.43000d0,-0.43430d0,-0.44120d0,-0.45510d0,-0.44820d0,-0.44070d0,-0.43400d0,-0.42790d0,&
-0.42230d0,-0.41700d0,-0.41190d0,-0.40700d0,-0.40230d0,-0.39780d0,-0.39330d0,-0.38900d0,-0.38480d0,-0.38080d0,-0.37680d0,&
-0.37290d0,-0.36910d0,-0.36550d0,-0.36190d0,-0.35840d0,-0.35490d0,-0.35160d0,-0.34830d0,-0.34510d0,-0.34200d0,-0.33890d0,&
-0.33590d0,-0.33300d0,-0.33010d0,-0.32730d0,-0.32460d0,-0.32190d0,-0.31930d0,-0.31670d0,-0.31410d0,-0.31170d0,-0.30920d0,&
-0.30680d0,-0.30450d0,-0.30220d0,-0.29990d0,-0.29770d0,-0.29550d0,-0.29340d0,-0.29130d0,-0.28920d0,-0.28720d0,-0.28520d0,&
-0.28320d0,-0.28130d0,-0.27940d0,-0.27750d0,-0.27570d0,-0.27390d0,-0.27210d0,-0.27040d0,-0.26870d0,-0.26700d0,-0.26530d0,&
-0.26370d0,-0.26200d0,-0.26040d0,-0.25890d0,-0.25730d0,-0.40670d0,-0.40670d0,-0.40680d0,-0.40680d0,-0.40690d0,-0.40700d0,&
-0.40720d0,-0.40730d0,-0.40750d0,-0.40780d0,-0.40800d0,-0.40830d0,-0.40860d0,-0.40890d0,-0.40930d0,-0.40970d0,-0.41010d0,&
-0.41060d0,-0.41100d0,-0.41150d0,-0.41210d0,-0.41260d0,-0.41320d0,-0.41380d0,-0.41440d0,-0.41500d0,-0.41570d0,-0.41630d0,&
-0.41700d0,-0.41760d0,-0.41830d0,-0.41900d0,-0.41970d0,-0.42030d0,-0.42100d0,-0.42160d0,-0.42230d0,-0.42290d0,-0.42360d0,&
-0.42420d0,-0.42470d0,-0.42530d0,-0.42590d0,-0.42640d0,-0.42690d0,-0.42730d0,-0.42780d0,-0.42820d0,-0.42860d0,-0.42890d0,&
-0.42930d0,-0.42960d0,-0.43000d0,-0.43030d0,-0.43070d0,-0.43100d0,-0.43140d0,-0.43180d0,-0.43230d0,-0.43270d0,-0.43320d0,&
-0.43370d0,-0.43430d0,-0.43490d0,-0.43550d0,-0.43620d0,-0.43690d0,-0.43770d0,-0.43850d0,-0.43940d0,-0.44030d0,-0.44130d0,&
-0.44230d0,-0.44340d0,-0.44440d0,-0.44550d0,-0.44670d0,-0.44790d0,-0.44920d0,-0.45050d0,-0.45200d0,-0.45350d0,-0.45520d0,&
-0.45700d0,-0.45900d0,-0.46120d0,-0.46380d0,-0.46690d0,-0.47050d0,-0.47510d0,-0.48150d0,-0.49340d0,-0.50360d0,-0.49810d0,&
-0.49100d0,-0.48400d0,-0.47730d0,-0.47090d0,-0.46500d0,-0.45940d0,-0.45400d0,-0.44890d0,-0.44400d0,-0.43920d0,-0.43460d0,&
-0.43020d0,-0.42590d0,-0.42170d0,-0.41760d0,-0.41360d0,-0.40980d0,-0.40600d0,-0.40230d0,-0.39870d0,-0.39520d0,-0.39170d0,&
-0.38840d0,-0.38510d0,-0.38190d0,-0.37870d0,-0.37560d0,-0.37260d0,-0.36970d0,-0.36670d0,-0.36390d0,-0.36110d0,-0.35840d0,&
-0.35570d0,-0.35300d0,-0.35050d0,-0.34790d0,-0.34540d0,-0.34300d0,-0.34050d0,-0.33820d0,-0.33590d0,-0.33360d0,-0.33130d0,&
-0.32910d0,-0.32690d0,-0.32480d0,-0.32270d0,-0.32060d0,-0.31860d0,-0.31660d0,-0.31460d0,-0.31260d0,-0.31070d0,-0.30880d0,&
-0.30700d0,-0.30510d0,-0.41860d0,-0.41860d0,-0.41870d0,-0.41880d0,-0.41890d0,-0.41910d0,-0.41920d0,-0.41950d0,-0.41970d0,&
-0.42000d0,-0.42040d0,-0.42070d0,-0.42110d0,-0.42150d0,-0.42200d0,-0.42250d0,-0.42300d0,-0.42360d0,-0.42410d0,-0.42480d0,&
-0.42540d0,-0.42610d0,-0.42680d0,-0.42750d0,-0.42820d0,-0.42900d0,-0.42970d0,-0.43050d0,-0.43130d0,-0.43210d0,-0.43290d0,&
-0.43370d0,-0.43450d0,-0.43540d0,-0.43620d0,-0.43700d0,-0.43780d0,-0.43860d0,-0.43940d0,-0.44020d0,-0.44100d0,-0.44180d0,&
-0.44250d0,-0.44320d0,-0.44400d0,-0.44470d0,-0.44530d0,-0.44600d0,-0.44660d0,-0.44730d0,-0.44790d0,-0.44850d0,-0.44910d0,&
-0.44970d0,-0.45030d0,-0.45090d0,-0.45160d0,-0.45220d0,-0.45290d0,-0.45360d0,-0.45430d0,-0.45500d0,-0.45580d0,-0.45660d0,&
-0.45740d0,-0.45830d0,-0.45920d0,-0.46010d0,-0.46110d0,-0.46210d0,-0.46320d0,-0.46430d0,-0.46540d0,-0.46650d0,-0.46770d0,&
-0.46890d0,-0.47010d0,-0.47130d0,-0.47260d0,-0.47400d0,-0.47540d0,-0.47690d0,-0.47840d0,-0.48010d0,-0.48180d0,-0.48360d0,&
-0.48550d0,-0.48770d0,-0.49010d0,-0.49270d0,-0.49560d0,-0.49900d0,-0.50300d0,-0.50780d0,-0.51420d0,-0.52380d0,-0.54550d0,&
-0.54160d0,-0.53660d0,-0.53010d0,-0.52330d0,-0.51660d0,-0.51000d0,-0.50370d0,-0.49760d0,-0.49190d0,-0.48630d0,-0.48100d0,&
-0.47590d0,-0.47100d0,-0.46620d0,-0.46160d0,-0.45720d0,-0.45280d0,-0.44860d0,-0.44450d0,-0.44060d0,-0.43670d0,-0.43290d0,&
-0.42920d0,-0.42560d0,-0.42200d0,-0.41860d0,-0.41520d0,-0.41190d0,-0.40860d0,-0.40550d0,-0.40240d0,-0.39930d0,-0.39630d0,&
-0.39340d0,-0.39050d0,-0.38770d0,-0.38490d0,-0.38220d0,-0.37950d0,-0.37690d0,-0.37430d0,-0.37180d0,-0.36930d0,-0.36680d0,&
-0.36440d0,-0.36200d0,-0.35970d0,-0.35740d0,-0.35520d0,-0.35290d0,-0.35070d0,-0.34860d0,-0.34650d0,-0.34440d0,-0.43170d0,&
-0.43170d0,-0.43180d0,-0.43190d0,-0.43200d0,-0.43220d0,-0.43240d0,-0.43270d0,-0.43290d0,-0.43330d0,-0.43360d0,-0.43400d0,&
-0.43440d0,-0.43490d0,-0.43540d0,-0.43590d0,-0.43650d0,-0.43710d0,-0.43770d0,-0.43840d0,-0.43900d0,-0.43970d0,-0.44050d0,&
-0.44120d0,-0.44200d0,-0.44280d0,-0.44360d0,-0.44450d0,-0.44530d0,-0.44620d0,-0.44700d0,-0.44790d0,-0.44880d0,-0.44970d0,&
-0.45060d0,-0.45150d0,-0.45240d0,-0.45330d0,-0.45420d0,-0.45500d0,-0.45590d0,-0.45680d0,-0.45770d0,-0.45850d0,-0.45930d0,&
-0.46020d0,-0.46100d0,-0.46180d0,-0.46260d0,-0.46330d0,-0.46410d0,-0.46490d0,-0.46560d0,-0.46640d0,-0.46720d0,-0.46790d0,&
-0.46870d0,-0.46950d0,-0.47040d0,-0.47120d0,-0.47210d0,-0.47290d0,-0.47380d0,-0.47470d0,-0.47570d0,-0.47670d0,-0.47770d0,&
-0.47870d0,-0.47980d0,-0.48090d0,-0.48210d0,-0.48320d0,-0.48440d0,-0.48560d0,-0.48680d0,-0.48810d0,-0.48940d0,-0.49060d0,&
-0.49200d0,-0.49340d0,-0.49480d0,-0.49620d0,-0.49770d0,-0.49930d0,-0.50090d0,-0.50260d0,-0.50440d0,-0.50640d0,-0.50840d0,&
-0.51060d0,-0.51290d0,-0.51540d0,-0.51820d0,-0.52120d0,-0.52470d0,-0.52870d0,-0.53330d0,-0.53910d0,-0.54690d0,-0.55950d0,&
-0.57810d0,-0.57410d0,-0.57020d0,-0.56420d0,-0.55790d0,-0.55130d0,-0.54480d0,-0.53830d0,-0.53200d0,-0.52590d0,-0.52000d0,&
-0.51430d0,-0.50880d0,-0.50350d0,-0.49840d0,-0.49350d0,-0.48870d0,-0.48400d0,-0.47950d0,-0.47510d0,-0.47090d0,-0.46680d0,&
-0.46270d0,-0.45880d0,-0.45490d0,-0.45120d0,-0.44750d0,-0.44390d0,-0.44040d0,-0.43700d0,-0.43370d0,-0.43040d0,-0.42720d0,&
-0.42400d0,-0.42090d0,-0.41790d0,-0.41490d0,-0.41200d0,-0.40910d0,-0.40630d0,-0.40350d0,-0.40080d0,-0.39810d0,-0.39550d0,&
-0.39290d0,-0.39040d0,-0.38790d0,-0.38540d0,-0.38300d0,-0.38060d0,-0.37820d0,-0.44480d0,-0.44480d0,-0.44490d0,-0.44500d0,&
-0.44510d0,-0.44530d0,-0.44550d0,-0.44570d0,-0.44600d0,-0.44640d0,-0.44670d0,-0.44710d0,-0.44760d0,-0.44800d0,-0.44860d0,&
-0.44910d0,-0.44970d0,-0.45030d0,-0.45090d0,-0.45160d0,-0.45220d0,-0.45300d0,-0.45370d0,-0.45450d0,-0.45530d0,-0.45610d0,&
-0.45690d0,-0.45780d0,-0.45860d0,-0.45950d0,-0.46040d0,-0.46130d0,-0.46220d0,-0.46310d0,-0.46410d0,-0.46500d0,-0.46590d0,&
-0.46690d0,-0.46780d0,-0.46870d0,-0.46970d0,-0.47060d0,-0.47150d0,-0.47240d0,-0.47330d0,-0.47420d0,-0.47510d0,-0.47600d0,&
-0.47690d0,-0.47780d0,-0.47860d0,-0.47950d0,-0.48030d0,-0.48120d0,-0.48210d0,-0.48300d0,-0.48390d0,-0.48480d0,-0.48570d0,&
-0.48660d0,-0.48760d0,-0.48850d0,-0.48950d0,-0.49050d0,-0.49150d0,-0.49260d0,-0.49370d0,-0.49480d0,-0.49590d0,-0.49710d0,&
-0.49830d0,-0.49950d0,-0.50080d0,-0.50200d0,-0.50330d0,-0.50460d0,-0.50590d0,-0.50720d0,-0.50860d0,-0.51000d0,-0.51140d0,&
-0.51290d0,-0.51440d0,-0.51590d0,-0.51750d0,-0.51910d0,-0.52090d0,-0.52270d0,-0.52460d0,-0.52660d0,-0.52860d0,-0.53080d0,&
-0.53320d0,-0.53570d0,-0.53850d0,-0.54150d0,-0.54470d0,-0.54850d0,-0.55270d0,-0.55760d0,-0.56360d0,-0.57120d0,-0.58280d0,&
-0.60870d0,-0.60460d0,-0.60110d0,-0.59620d0,-0.59040d0,-0.58430d0,-0.57800d0,-0.57160d0,-0.56520d0,-0.55890d0,-0.55270d0,&
-0.54670d0,-0.54080d0,-0.53520d0,-0.52960d0,-0.52430d0,-0.51910d0,-0.51410d0,-0.50930d0,-0.50450d0,-0.49990d0,-0.49550d0,&
-0.49110d0,-0.48690d0,-0.48270d0,-0.47870d0,-0.47480d0,-0.47090d0,-0.46720d0,-0.46350d0,-0.45990d0,-0.45640d0,-0.45300d0,&
-0.44960d0,-0.44630d0,-0.44310d0,-0.43990d0,-0.43680d0,-0.43370d0,-0.43070d0,-0.42780d0,-0.42490d0,-0.42210d0,-0.41930d0,&
-0.41650d0,-0.41380d0,-0.41120d0,-0.40860d0/)
 end module sdata

! AMI: removoved program reshape_sdata – was unused and causing issues.

! AMI: got rid of module Eparallel
!==============================================================================
module DifferentialEvolution
  use param,only : dp
  implicit none

! subroutine DEM(errsub,XCmin,XCmax,bestmem_XC,Dim_XC)
! Differential Evolution for Optimal Control Problems.
!
! Written by Dr. Feng-Sheng Wang, Department of Chemical Engineering,
! National Chung Cheng University, Chia-Yi 621, Taiwan. 
! Ref.:
! R. STORN and K.V. PRICE, J. Global Optimization 11, 341-59 (1997). 
! http://www1.icsi.berkeley.edu/~storn/code.html
!
! intent(in) from external module:
!   Dim_XC : Dimension of decision parameters.
! intent(in):
!   NP : Population size.
!   itermax : Maximum number of iteration.
!   strategy=1,2,...,6 : Strategy of mutation operations, 2 is default.
!   method(1) = 0, Fixed mutation scaling factors F_XC,
!                  = 1, Random mutation scaling factors F_XC=[0, 1],
!                  = 2, Random mutation scaling factors F_XC=[-1, 1]. 
!   method(2) = 1, Random combination factor F_CR used for 
!                         strategy = 6 in mutation operation, 
!                   = 0, fixed combined factor provided by the user .
!   method(3) = 1, Saving results in external data file.
!                   = other, displaying results only.
!   errf : User provided subroutine errf(xc,fitness), where "xc" is 
!       decision parameter vector and "fitness" is fitness value.
!   XCmin(Dim_XC) : Lower bound of decision parameters.
!   XCmax(Dim_XC) : Upper bound of decision parameters.
!   VTR : Desired fitness value.
!   CR_XC : Crossover factor for decision parameters.
!   refresh : Intermediate output produced after "refresh" iterations,
!       with "refresh < 1" no intermediate output.
!
! intent(in) and modified (depending on method):
!   F_XC : Mutation scaling factor for decision parameters.
!   F_CR : Mutation scaling factor for crossover factor, with strategy=6.
!
! intent(out):
!   bestval : Best value of fitness in error function.
!   bestmen_XC(Dim_XC) : Best decision parameters.
!   nfeval : Number of function calls.
!   itval : Iteration for which the bestval occurred.

  integer :: NP,itermax,strategy,refresh,method(3),nfeval,itval
  real(dp) :: VTR,CR_XC,F_XC,F_CR,bestval
contains
  !--------------------------------------------------------------------
  subroutine DEM(errsub,XCmin,XCmax,bestmem_XC,Dim_XC)
  implicit none
  integer,intent(in) :: Dim_XC
  real(dp),intent(in)  :: XCmin(Dim_XC),XCmax(Dim_XC)  
  real(dp),intent(out) :: bestmem_XC(Dim_XC)
  !    
  integer :: i, ibest, iter
  integer, dimension(NP) :: rot, a1, a2, a3, a4, a5, rt
  integer, dimension(4) :: ind
  real(dp),parameter :: bohr=0.5291772083d0  
  real(dp) :: tempval
  real(dp), dimension(NP,Dim_XC) :: pop_XC,bm_XC,mui_XC,mpo_XC,&
    popold_XC, rand_XC, ui_XC
  real(dp), dimension(NP) :: val
  real(dp), dimension(Dim_XC) :: bestmemit_XC, rand_C1, tmp
  external errsub
  !
  ! Initialize a population. 
  pop_XC=0.0d0
  do i=1,NP
    call random_number(rand_C1)
    pop_XC(i,:)=XCmin+rand_C1*(XCmax-XCmin)
  end do

  ! Evaluate fitness functions and find the best member.
  val=0.0d0
  nfeval=0
  ibest=1
  tmp(:)=pop_XC(1,:)
  call errsub(tmp, val(1))
  bestval=val(1)
  nfeval=nfeval+1
  do i=2,NP
    tmp(:)=pop_XC(i,:)
    call errsub(tmp,val(i))
    nfeval=nfeval+1
    if (val(i) < bestval) then
      ibest=i
      bestval=val(i)
    end if
  end do
  bestmemit_XC=pop_XC(ibest,:)
  bestmem_XC(:)=bestmemit_XC(:)
  itval=1

  ! Perform evolutionary computation.
  bm_XC=0.0d0
  rot=(/(i,i=0,NP-1)/)
  iter=1
  do while (iter <= itermax)
    popold_XC=pop_XC
    ! Mutation operation.
    ind=randperm(4)
    a1=randperm(NP)
    rt=mod(rot+ind(1),NP)
    a2=a1(rt+1)
    rt=mod(rot+ind(2),NP)
    a3=a2(rt+1)
    rt=mod(rot+ind(3),NP)
    a4=a3(rt+1)
    rt=mod(rot+ind(4),NP)
    a5=a4(rt+1)
    bm_XC=spread(bestmemit_XC, DIM=1, NCOPIES=NP)

    ! Generating a random sacling factor.
    select case (method(1))
      case (1)
        call random_number(F_XC)
      case(2)
        call random_number(F_XC)
        F_XC=2.0d0*F_XC-1.0d0
    end select

    ! select a mutation strategy.
    select case (strategy)
    case (1)
      ui_XC=bm_XC + &
          F_XC*(popold_XC(a1,:)-popold_XC(a2,:))
    case default
      ui_XC=popold_XC(a3,:) + &
          F_XC*(popold_XC(a1,:)-popold_XC(a2,:))
    case (3)
      ui_XC=popold_XC + &
          F_XC*(bm_XC-popold_XC + popold_XC(a1,:)-popold_XC(a2,:))
    case (4)
      ui_XC=bm_XC + &
          F_XC*(popold_XC(a1,:)-popold_XC(a2,:) + &
          popold_XC(a3,:)-popold_XC(a4,:))
    case (5)
      ui_XC=popold_XC(a5,:) + &
          F_XC*(popold_XC(a1,:)-popold_XC(a2,:) + &
          popold_XC(a3,:)-popold_XC(a4,:))
    case (6) ! A linear crossover combination of bm_XC and popold_XC
      if (method(2) == 1) call random_number(F_CR) 
      ui_XC=popold_XC + &
          F_CR*(bm_XC - popold_XC) + &
          F_XC*(popold_XC(a1,:) - popold_XC(a2,:))
    end select

    ! Crossover operation.
    call random_number(rand_XC)
    mui_XC=0.0d0
    mpo_XC=0.0d0
    where (rand_XC < CR_XC)
      mui_XC=1.0d0
      ! mpo_XC=0.0d0
    elsewhere
      ! mui_XC=0.0d0
      mpo_XC=1.0d0
    end where
    ui_XC=popold_XC*mpo_XC+ui_XC*mui_XC

    ! Evaluate fitness functions and find the best member.
    do i=1,NP
      ! Confine each of feasible individuals in the lower-upper bound.
      ui_XC(i,:)=max(min(ui_XC(i,:),XCmax),XCmin)
      tmp(:)=ui_XC(i,:)
      call errsub(tmp,tempval)
      nfeval=nfeval+1
      if (tempval < val(i)) then
        pop_XC(i,:)=ui_XC(i,:)
        val(i)=tempval
        if (tempval < bestval) then
          bestval=tempval
          itval=iter
          bestmem_XC(:)=ui_XC(i,:)
        end if
      end if
    end do
    bestmemit_XC(:)=bestmem_XC(:)
    if( (refresh > 0) .and. (mod(iter,refresh)==0)) then
      if (method(3)==1) write(61,203) iter
      write(*,203) iter
      do i=1,Dim_XC
        if (method(3)==1) write(61,202) i, bestmem_XC(i)*bohr
        write(*,202) i,bestmem_XC(i)*bohr
      end do
      if (method(3)==1) write(61,201) bestval
        write(*,201) bestval
    end if
    iter=iter+1
    if(bestval <= VTR .and. refresh > 0) then
      exit
    endif
  end do
  return
  201 format(2x, 'bestval =', es9.1 /)
  202 format(5x, 'bestmem_XC(', I3, ') =', es10.3)
  203 format(2x, 'No. of iteration  =', I8)
  end subroutine DEM
  !-------------------------------------------------------------------
  function randperm(num)
  implicit none
  integer, intent(in) :: num
  integer :: number, i, j, k
  integer, dimension(num) :: randperm
  real(dp), dimension(num) :: rand2
  intrinsic random_number
  call random_number(rand2)
  do i=1,num
     number=1
     do j=1,num
        if (rand2(i) > rand2(j)) then
          number=number+1
        end if
     end do
     do k=1,i-1
        if (rand2(i) <= rand2(k) .and. rand2(i) >= rand2(k)) then
          number=number+1
        end if
     end do
     randperm(i)=number
  end do
  return
  end function randperm
end module DifferentialEvolution

!======================================================================
