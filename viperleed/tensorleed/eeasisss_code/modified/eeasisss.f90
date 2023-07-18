!VERSION: November 25, 2006.
!
!THIS IS A BETA VERSION. IT SHOULD NOT BE REDISTRIBUTED.
!IN CASE !ANYBODY IS INTERESTED IN THE PROGRAM, PLEASE CONTACT 
!THE AUTHOR, !jru@kth.se .
!
!program EEASiSSS,
!'Elastic Electron-Atom Scattering in Solids and Surface Slabs',
!author: John Rundgren.
!The author appreciates acknowledgement in publications by citing:
!J. Rundgren, Phys. Rev. B 68 125405 (2003).
!----------------------------------------------------------------------- 
module param
!basic constants.
!alfa=fine-structure constant, c=speed of light, cm2=1/c**2.
real*8 :: pi
real*8,parameter :: rydb=13.60569172d0,&
                    bohr=0.5291772083d0,bohr3=bohr**3
real*8,parameter :: alfa=7.297352533d-03,c_light=2.d0/alfa,&
                    cm2=alfa*alfa*0.25d0
real*8 :: cm2tmp

!select program operations.
character(len=1) :: BulkOrSlab,RelOrNonrel,Output,RhoPot,&
  SpinPhaseShift,Omega,Theta,Sherman,CrossSection,PSvsE,UvsR,&
  NMselect
!identify the atoms.
character(len=9),allocatable :: id_at(:)
 
!deltasup=discriminator on superpositon=low limit on charge density and
!         abs(potential); regulates the #neighbors in superposition.
!regulates the #neighbors in superposition
!radfac=factor in Mattheiss prescription (structure),
!deltarmt=allowed interatomic distance interval (mt_optimisation);
!         the value 0.05 AA gives a mean rmt for 12 nneighbors in Be.
real*8,parameter :: deltasup=1.d-04
real*8,parameter :: radfac=0.9d0
real*8,parameter :: deltarmt=0.05/bohr
real*8 :: NMlambda,NMeps,NMepsit
integer :: NMiter
 
!parameters of the integration routines ode.f and rkf45.f .
!node  =max. no. of integration steps.
integer,parameter :: node=4000
real*8 :: macheps
 
!wave_start=start of u1(kappa,r) integration,
!eta_cutoff=phase shift cut-off,
!relerr0   =temporary relative error for investigating the magnitude
!           of wave amplitude and phase,
!relerr1   =relative error of wave_integrate_ode
!relerr2   =relative error of phase_integrate_ode
real*8,parameter :: wave_start=1.d-12
real*8,parameter :: eta_cutoff=1.d-10
real*8,parameter :: relerr0=1.d-02,relerr1=1.d-08,relerr2=1.d-08

!cls =core level shift.
real*8,allocatable :: cls(:)

!self-energy data.
integer :: nsr,nsp
real*8,allocatable :: sp(:),sr(:),sdat(:,:)
 
!eev(1:ne)=energy grid,
![ep1,ep2]=energy interval for phase shift calculation,
![ie1,ie2]=subscript interval for phase shift calculation,
!dep      =energy step in phase shift table.
integer :: ne,ie1,ie2
real*8  :: ep1,ep2,dep
real*8,allocatable,dimension(:) :: eev
end module param
!----------------------------------------------------------------------
 module dimen
!nieq  =no. of inequivalent atoms in the structure,
!nocc  =no. of occupied atom types in a slab,
!neq(i)=no. of equivalent atoms for i=1:nieq,
!ncell =no. of atoms in the unit cell,
!nx    =no. of grid points in atomic calculation,
!nshell=no. of neighbor shells in subroutine neigbor_sum,
!lmax  =max orbital quantum number l,
integer :: nieq,nocc,ncell,nx,lmax,nshell
integer,allocatable :: neq(:)
end module dimen
!-----------------------------------------------------------------------
module radii
!rx(1:nx)   =radial grid, rx(i)=r(1)*exp(dx*(i-1)),
!dx    =logarithmic increment in radial grid,
!rmt(1:nieq)=muffin-tin radius,
!rxt(1:nieq)=external atomic radius,
!nxt(1:nieq)=grid point just outside rxt,
!qmt(1:nocc)=ElCharge within MT radius
real*8 :: dx
real*8 ,allocatable,dimension(:) :: rx,rmt,rxt,qmt
integer,allocatable,dimension(:) :: nxt
end module radii
!-----------------------------------------------------------------------
module scatt
!iatom =atom subscript,
!l     =orbital quantum number,
!kappa =spin-orbit coupling quantum number (dirac eq.):
!      =-l-1, l=0,1,2,... and =l, l=1,2,...(dirac eq.),
!      =l, l=0,1,2,... (schroedinger eq.),
!emv0  =Eprimary-V0(inner potential),
!q     =wave number corresponding to emv0.,
!a1    =aprime/a in effective Dirac potential veff.
integer :: iatom,l,kappa
real*8  :: emv0,q,a1
end module scatt
!-----------------------------------------------------------------------
module bessl
!bess=temporary bessel function.
real*8,allocatable :: bess(:)
end module bessl
!-----------------------------------------------------------------------
module poten
!mtpot ='MT potential' or 'MT potential * radius',
!ispot =interstitial potential,
!v0tab =inner potential.
character :: xc*1
real*8,target,allocatable,dimension(:,:) :: mtpot,ispot
real*8,allocatable :: z(:),v0tab(:)
end module poten
!-------------------------------------------------------------------------      
module radeq
!radial Dirac and Schroedinger waves.
!Ref.: M.E.Rose, Relativistic Electron Theory (Wiley, 1961).
use param,only : RelOrNonrel,cm2tmp
use radii,only : rx,dx,nxt
use scatt,only : iatom,l,kappa,emv0,q,a1
use bessl,only : bess
use poten,only : mtpot
contains

  subroutine waveq(x,y,dydx)
  implicit none
  integer :: i
  real*8  :: y(2),dydx(2),x,p,vmv0x,vmv0,qbx,emv
  real*8,parameter :: a6=1.d0/6.d0
  !lagrangian 4-point interpolation.
  p=log(x/rx(1))/dx+1.d0
  i=max(int(p),2); i=min(i,nxt(iatom)-2)
  p=p-i
  vmv0x=   -a6*p*(p-1.d0)*(p-2.d0)*mtpot(i-1,iatom)&
        +0.5d0*(p*p-1.d0)*(p-2.d0)*mtpot(i  ,iatom)&
        -0.5d0*p*(p+1.d0)*(p-2.d0)*mtpot(i+1,iatom)&
                  +a6*p*(p*p-1.d0)*mtpot(i+2,iatom)
  emv=emv0-vmv0x/x
  qbx=kappa/x
  dydx(1)=-y(1)*qbx+y(2)*(1.d0+cm2tmp*emv)
  dydx(2)= y(2)*qbx-y(1)*emv
  return
  end subroutine waveq
end module radeq
!-----------------------------------------------------------------------
module error
contains

  real*8 function mt_err(x)
  !MT_ERR is mean square MT core jump to a common MT floor.
  use dimen,only : neq,nocc
  use radii,only : rxt,nxt
  use poten,only : mtpot,ispot
  implicit none
  integer :: n,neqsum
  real*8 :: x(nocc),vrmt(nocc),stp(nocc),avstp,v0,vol,dum,integral_egt,&
            interpol_egl3
  !v0   =interstitial potential.
  !stp  =potential step at MT radius.
  !avstp=average potential step at MT radius.
  v0=0.d0; vol=0.d0; neqsum=0.d0
  do n=1,nocc
    if(x(n)>rxt(n)) x(n)=rxt(n) !simplex radius <= charge radius.
    vol=vol+neq(n)*(rxt(n)**3-x(n)**3)/3.d0
    v0=v0+neq(n)*&
       integral_egt(3,nxt(n),ispot(1,n),x(n),rxt(n),dum)
    vrmt(n)=interpol_egl3(x(n),mtpot(1,n))
    neqsum=neqsum+neq(n)
  enddo
  v0=v0/vol
  stp=vrmt-v0
  !step average.
  avstp=0.d0
  do n=1,nocc
    avstp=avstp+neq(n)*stp(n)
  enddo
  avstp=avstp/neqsum
  !mean square relative step.
  mt_err=0.d0
  do n=1,nocc
    mt_err=mt_err+neq(n)*(stp(n)-avstp)**2
  enddo
  mt_err=mt_err/neqsum
  return
  end function mt_err

  real*8 function v0_err(x)
  !ERR_V0 is mean square error of a four-parameter fit to the energy-
  !dependent inner potential.
  use param,only : rydb,eev,ne
  use poten,only : v0tab
  implicit none
  integer :: ie
  real*8  :: x(4)
  !v0_err expression modified against negative energy shifts x(3).
  v0_err=0.d0
  do ie=1,ne
    v0_err=v0_err+(v0tab(ie)*rydb-&
          max(x(1),x(2)+x(3)/sqrt(max(eev(ie)+x(4),1.d0))))**2
  enddo
  v0_err=v0_err/ne
  end function v0_err
end module error
!-----------------------------------------------------------------------
program EEASiSSS

!'Elastic Electron-Atom Scattering in Solids and Surface Slabs',
!author: John Rundgren.

use param,only : BulkOrSlab,RelOrNonrel,Output,RhoPot,SpinPhaseShift,&
  Omega,Theta,Sherman,CrossSection,PSvsE,UvsR,macheps,pi,rydb,bohr,&
  cm2,ne,eev,id_at,cls,nsr,nsp,sr,sp,sdat,ie1,ie2,ep1,ep2,dep,&
  NMselect,NMlambda,NMiter,NMeps,NMepsit
use dimen,only : nshell,nieq,neq,nocc,ncell,nx,lmax
use radii,only : rx,dx
use poten,only : z
implicit none
character :: title*34,xcinput*80
integer :: i,j,jj,ir,inr,iir,nx1,l,last_interval,nes,ne1,k
real*8  :: rc(3,3),spa,spaa,rcc1,rcc2,rcc3,voluc,baseuc,zatm,rmin,rmax,&
  rmin1,rmax1,dum,zsum,w1,w2,es1,des,q,ery,xi,fi,qstart,qstop,&
  dq,rmt79,dummy,vatm,rsfac
character,allocatable :: chgden(:)*80
integer,allocatable :: na(:,:),ia(:,:),ncon(:)
real*8,allocatable,dimension(:)   :: valieq,eev1
real*8,allocatable,dimension(:,:) :: rk,rk_tmp,ad,atrho4pir2,swrho4pi,&
  swpot,rs,myxc

pi=acos(-1.d0)
macheps=epsilon(dummy)
open(61,file='logfile',status='replace')
write(61,'(a)')'STRUCTURE:'
read(5,'(/a)') title; write(61,'(a)') title
!spaa   =lattice constant in Angstroms,
!spa    =lattice constant in a.u.,
!rc(i,j)=i'th coordinate of the j'th axis of unit cell (au),
!rk(i,j)=i'th coordinate of the j'th atom in unit cell (au),
read(5,*) BulkOrSlab,spaa,nshell
spa=spaa/bohr; write(61,200) spaa,spa,BulkOrSlab,nshell
do j=1,3
  read(5,*) (rc(i,j),i=1,3)
enddo
rc=spa*rc
!nieq      =number of ineqivalent atoms in unit cell,
!nocc      =number of occupied types in a slab,
!neq(ir)   =number of equivalent atoms of type ir,
!z(ir)     =atomic number of type ir,
!valieq(ir)=valence charge of type ir,
!ncell     =total number of atoms in unit cell.
read(5,*) nieq,nocc

allocate(z(nieq),neq(nieq),cls(nieq),valieq(nieq),id_at(nieq))
allocate(chgden(nieq),rk_tmp(3,10000))

ncell=0
do ir=1,nieq
  read(5,*) neq(ir),z(ir),valieq(ir),cls(ir),chgden(ir)
  do j=1,neq(ir)
    ncell=ncell+1
    read(5,*) rk_tmp(1:3,ncell)
    rk_tmp(1:3,ncell)=rk_tmp(1:3,ncell)*spa
  enddo
enddo
!atom identification.
call id_atom

!are AtomType z coordinates in nondecreasing order?
do i=2,ncell
  if(rk_tmp(3,i)-rk_tmp(3,i-1)<0.d0)then
    write(61,*)'stop: AtomType z coordinates must be nondecreasing'
    write(61,*)i,rk_tmp(3,i)-rk_tmp(3,i-1)
    stop
  endif
enddo

!reduction to nonredundant dimensions, that is, 10000 to ncell.
allocate(rk(3,ncell))
rk(1:3,1:ncell)=rk_tmp(1:3,1:ncell)
deallocate(rk_tmp)
!volume of unit cell, rcc1 and rcc2 are || to the surface;
rcc1=rc(2,1)*rc(3,2)-rc(3,1)*rc(2,2)
rcc2=rc(3,1)*rc(1,2)-rc(1,1)*rc(3,2)
rcc3=rc(1,1)*rc(2,2)-rc(2,1)*rc(1,2)
!voluc=volume of unit cell,
!baseuc=area of unit cell || surface, used in case('s').
select case(BulkOrSlab)
 case('b')
  voluc=abs(rc(1,3)*rcc1+rc(2,3)*rcc2+rc(3,3)*rcc3)
  baseuc=0.d0
 case('s')
  rc(3,3)=rk(3,ncell)-rk(3,1)
  voluc=abs(rc(1,3)*rcc1+rc(2,3)*rcc2+rc(3,3)*rcc3)
  baseuc=abs(rc(1,1)*rc(2,2)-rc(2,1)*rc(1,2))
end select
write(61,201) ((rc(i,j)*bohr,i=1,3),j=1,3)
write(61,202) voluc*bohr**3
write(61,203) ncell
jj=0
do ir=1,nieq
  write(61,204) ir,neq(ir),z(ir),valieq(ir),cls(ir)
  inr=neq(ir)
  do iir=1,inr
    jj=jj+1
    write(61,205) (rk(i,jj)*bohr,i=1,3)
  enddo
enddo
!cls, input in eV, calculation in rydbergs.
 cls=cls/rydb
  
  200 format ('lattice constant=',f8.4,' AA, ',f8.4,' au',&
    &', bulk/slab=',a,', nshell=',i3)
  201 format(/'axes of unit cell (AA)'/(4x,3f9.4))
  202 format('UnitCellVolume(AA**3)',f9.3)
  203 format('no. atoms in unit cell',i3)
  204 format(/'type',i2,', NoAtoms',i2,', AtomicNo',f5.1,&
            &', Valence(eV)',f5.1,', CoreLevelShift(eV)',f6.2/&
            &8x,'Coordinates(AA):')
  205 format(3f8.4)
  208 format(/'allowed interatomic distance interval (AA):',f6.2)
  206 format( 'atomic occupation number',8(f4.0,2x):/(24x,8(f4.0,2x):))
  207 format( 'core level shift',8x,8f6.2:/(24x,8f6.2:)) 
  209 format('unit-cell volume',f10.4,', occupied volume',f10.4)

!CHARGE DENSITY indata.
do ir=1,nieq
  j=len_trim(chgden(ir))
  open(50,file=chgden(ir)(1:j))
  read(50,*)
  !read exponential grid.
  read(50,*) zatm,nx,rmin,rmax
  if(abs(zatm-z(ir))>1.d-03)then
    write(61,*)&
    'charge_density stop: ','atom',ir,', zsup=',z(ir),' but zatm=',zatm
    stop
  endif
  if(ir==1)then
  
    allocate(rx(nx),atrho4pir2(nx,nieq),swrho4pi(nx,nieq))
    allocate(swpot(nx,nieq))
  
    rmin1=rmin; rmax1=rmax; nx1=nx
    dx=log(rmax/rmin)/(nx-1)
    do i=1,nx
      rx(i)=rmin*exp(dx*(i-1))
    enddo
  else
    if(nx1/=nx.or.abs(rmin-rmin1)+abs(rmax-rmax1)>1.d-06)then
      write(61,*)'charge_density stop: '&
           ,'all atoms must have the same rmin,rmax,nx.'; stop
    endif
  endif
  !read charge density.
  read(50,*) (dum,atrho4pir2(i,ir),i=1,nx)
  zsum=0.d0
  w1=0.d0
  do i=2,nx
    w2=atrho4pir2(i,ir)*rx(i)
    zsum=zsum+w1+w2
    w1=w2
  enddo
  zsum=zsum*0.5d0*dx
  w2=zatm-valieq(ir)
  write(61,'(a,i2,a,f7.4,a,f7.4)')&
    'charge_density: atom ',ir,', zatm-val=',w2,', zel=',zsum
  if(abs(w2-zsum)>1.d-03)then
    write(61,'(a)')'stop'; stop
  endif
  call rho_center(atrho4pir2,ir)
  close(50)
enddo

!SCATTERING indata.
read(5,*)
!exchange-correlation.
read(5,*) xcinput
j=len_trim(xcinput)
open(50,file=xcinput(1:j))
read(50,*) nsr,nsp
allocate(sr(nsr),sp(nsp),sdat(nsp,nsr))
read(50,*) sr
do i=1,nsp
  read(50,*) sp(i),sdat(i,1:nsr)
enddo
close(50)
!no. phase shifts. 
read(5,*) lmax
!modes of calculation.
read(5,*) RelOrNonrel
    write(61,'(2a)')'RelOrNonrel=',RelOrNonrel
read(5,*) Output
    write(61,'(2a)')'Output=',Output
read(5,*) SpinPhaseShift
    write(61,'(2a)')'SpinPhaseShift=',SpinPhaseShift
read(5,*) Omega
    write(61,'(2a)')'Omega=',Omega
read(5,*) Theta
    write(61,'(2a)')'Theta=',Theta
read(5,*) Sherman
    write(61,'(2a)')'Sherman=',Sherman
read(5,*) CrossSection
    write(61,'(2a)')'CrossSection=',CrossSection
read(5,*) RhoPot
    write(61,'(2a)')'RhoPot=',RhoPot
read(5,*) PSvsE
    write(61,'(2a)')'PSvsE=',PSvsE
read(5,*) UvsR
    write(61,'(2a)')'UvsR=',UvsR
read(5,*) ep1,ep2,dep
    write(61,'(a,3f6.0)')'ep1,ep2,dep=',ep1,ep2,dep
read(5,*) es1,des,nes
    write(61,'(a,2f6.0,i6)')'es1,des,nes=',es1,des,nes
read(5,*)
    write(61,'(a)')'POTENTIAL OPTIMIZATION (Nelder-Mead)'
read(5,*) NMselect
    write(61,'(2a)')'NMselect=',NMselect
read(5,*) NMlambda
    write(61,'(a,f4.2)')'NMlambda=',NMlambda
read(5,*) NMiter
    write(61,'(a,i2)')'NMiter=',NMiter
read(5,*) NMeps
    write(61,'(a,1pe6.0)')'NMeps=',NMeps
read(5,*) NMepsit
    write(61,'(a,f6.4)')'NMepsit=',NMepsit
select case(Output)
 case('p')
  allocate(eev1(202))
  i=1; eev1(1)=3
  do
    if(eev1(i)==45)exit
    i=i+1; eev1(i)=eev1(i-1)+3.d0
  enddo
  do
    if(eev1(i)==60.d0)exit
    i=i+1; eev1(i)=eev1(i-1)+5.d0
  enddo
  do
    if(eev1(i)==100.d0)exit
    i=i+1; eev1(i)=eev1(i-1)+10.d0
  enddo
  do
    if(eev1(i)==300.d0)exit
    i=i+1; eev1(i)=eev1(i-1)+20.d0
  enddo
  do
    if(eev1(i)==600.d0)exit
    i=i+1; eev1(i)=eev1(i-1)+30.d0
  enddo
  do
    if(eev1(i)==1000.d0)exit
    i=i+1; eev1(i)=eev1(i-1)+40.d0
  enddo
  do
    if(eev1(i)==2000.d0)exit
    i=i+1; eev1(i)=eev1(i-1)+50.d0
  enddo
  do
    if(eev1(i)==10000.d0)exit
    i=i+1; eev1(i)=eev1(i-1)+100.d0
  enddo
  do
    if(eev1(i)==20000.d0)exit
    i=i+1; eev1(i)=eev1(i-1)+200.d0; ne1=i
  enddo
  if(ep2>20000.d0) ep2=20000.d0
endselect
select case(Output)
 case('p')
  !least no. of energy values enclosing the intervals [ep1,ep2].
  ie1=10000; ie2=0
  call energy_interval(ep1,ep2,ie1,ie2)
  ne=ie2-ie1+1
  allocate(eev(ne))
  do i=ie1,ie2
    eev(i-ie1+1)=eev1(i)
  enddo
  deallocate(eev1)
 case('s','d')
  if(es1+des*(nes-1)>20000.d0) nes=1+int((20000.d0-es1)/des)
  !energy grid uniform in energy.
  ne=nes
  allocate(eev(ne))
  do i=1,ne
    eev(i)=es1+(i-1)*des
  enddo
end select

!EXECUTION.
allocate(ad(nshell,nieq),na(nshell,nieq),ia(nshell,nieq),ncon(nieq))
allocate(rs(nx,nieq),myxc(nx,nieq))
call structure_I&
  (atrho4pir2,swrho4pi,swpot,rs,rsfac,myxc,ad,na,ia,ncon,rc,rk,&
  valieq,voluc,baseuc)
call scattering&
  (swrho4pi,swpot,rs,rsfac,myxc,z,ad,na,ia,ncon,valieq,title)
stop
contains

  subroutine rho_center(rho4pir2,ir)
  !correction: rho(1:16) calculated by E.L.Shirley's atomic program
  !contains minor discontiuities. A spline calculation on rho(17:30) is
  !used for smoothing rho(1:16) by extrapolation.
  !The smoothed rho can be studied on plots of the files a4pirho,
  !a4pirhop, and a4pirhob, see subroutine structure_I.
  integer :: ir,i1,i2,n
  real*8  :: rho4pir2(nx,nieq),x(14),f(14),fp(14),fb(14),ft(14)
  i1=17; i2=30; n=i2-i1+1
  j=0
  do i=i1,i2
    j=j+1
    x(j)=i
    f(j)=rho4pir2(i,ir)/(rx(i)*rx(i))
  enddo
  call spcoef(n,x,f,1.d+33,1.d+33,fp,fb,ft,'main')
  last_interval=1
  do i=1,i1-1
    xi=i
    call svalue(n,x,f,fp,fb,ft,xi,fi,dum,last_interval,'extrapolation')
    rho4pir2(i,ir)=fi*rx(i)*rx(i)
  enddo
  return
  end subroutine rho_center
  
  subroutine energy_interval(e1,e2,i1,i2)
  !enclose interval [e1,e2] in energy grid interval [eev1(i1),eev1(i2)].
  integer :: i1,i2
  real*8  :: e1,e2
  if(e1<=eev1(1))then
    i1=1
  else
    do i=2,ne1
      if(e1<eev1(i))then
        i1=i-1
        exit
      elseif(e1==eev1(i))then
        i1=i
        exit
      endif
    enddo
  endif
  if(e2>=eev1(ne1))then
    i2=ne1
  else
    do i=ne1-1,1,-1
      if(e2>eev1(i))then
        i2=i+1
        exit
      elseif(e2==eev1(i))then
        i2=i
        exit
      endif
    enddo
  endif
  select case(Output)
  case('p')
    !spcoef requires i2-i1+1>=3 for end-point derivatives.
    if(i2-i1+1<4)then
      i1=max(i1-1,1)
      i2=min(i2+1,ne1)
    endif
    if(i1==ne1-2) i1=ne1-3
    if(i2==3) i2=4
  endselect
  return
  end subroutine energy_interval

  subroutine id_atom  !subscript of ineqvivalent atom in filename.
  do ir=1,nieq
    id_at(ir)='.00.     '
    i=z(ir)
    if(i<10)then
      write(id_at(ir)(3:3),'(i1)') i
    elseif(i<100)then
      write(id_at(ir)(2:3),'(i2)') i
    endif
    if(ir<10 )then
      write(id_at(ir)(5:5),'(i1)') ir
    elseif(ir<100)then
      write(id_at(ir)(5:6),'(i2)') ir
    elseif(ir<1000)then
      write(id_at(ir)(5:7),'(i3)') ir
    elseif(ir<10000)then
      write(id_at(ir)(5:8),'(i4)') ir
    elseif(ir<100000)then
      write(id_at(ir)(5:9),'(i5)') ir
    endif
  enddo
  return
  end subroutine id_atom
end program EEASiSSS
!-----------------------------------------------------------------------
subroutine structure_I&
  (atrho4pir2,swrho4pi,swpot,rs,rsfac,myxc,ad,na,ia,ncon,rc,rk,&
  valieq,voluc,baseuc)
!STRUCTURE calculates muffin-tin potentials from superposed atomic
!potentials.
!(1) Lattice: A 3D structure is generated from a given unit cell.
!    A surface slab of n layers is modeled by a unit cell whose depth
!    is equal to n-1 interlayer separations plus a "large" distance.
!    Hence, the atoms situated near the surface of the slab get
!    negligible overlap with the atoms next to the bulk cut-off.
!(2) Superposition: The charge density of an atom A in the unit cell
!    and the charge density contributions from all neighbours of A are
!    added. The same procedure for the Hartree potential.
!(3) Spherically symmetric approximation: A's superposition density and
!    superposition potential are considered as expanded in multipoles
!    around A. The zeroth multipole is kept.
!(4) The outer radius of an atom in the material is determined by the
!    requirement that the atom is neutral. This is the so-called Norman
!    atomic radius.

!Ref.: L.F. Mattheiss, Phys. Rev. 133, A1399 (1964).
!Ref.: T.L. Loucks, Augmented plane wave method (Benjamin, New York,
!      1967).
!Ref.: J.G. Norman, Mol. Phys. 31, 1191 (1975).
!Ref.: A. Barbieri and M.A. Van Hove, Phase shift package,
!      http://electron.lbl.gov/software.
!Ref.: DLPHASE v1.1, CCP3 Library code,
!      http://www.cse.clrc.ac.uk/Activity.
use param,only : Output,RhoPot,pi,rydb,bohr,bohr3,radfac,id_at,cls
use dimen,only : nieq,neq,nocc,ncell,nshell,nx
use radii,only : rx,dx,rxt,nxt,rmt,qmt
use poten,only : z,mtpot,ispot
implicit none
character :: xc*1
integer :: na(nshell,nieq),ia(nshell,nieq),ncon(nieq),ir,i,k,n
real*8  :: rc(3,3),rk(3,ncell),atrho4pir2(nx,nieq),swrho4pi(nx,nieq),&
  swpot(nx,nieq),ad(nshell,nieq),rs(nx,nieq),myxc(nx,nieq),&
  valieq(nieq),x,x1,d,fp1,fpn,voluc,baseuc,vdum,volmt,volcs,volis,&
  qxmt,pi2,fpi,fpi3,facmy,ufpi,rsfac,beta,myx,dum,&
  integral_egt,interstitial_pot
real*8,parameter :: thrd=1.d0/3.d0
integer,allocatable :: imax(:)
real*8,allocatable,dimension(:,:) :: atpot,atrho4pi
pi2=pi*pi; fpi=4.d0*pi; fpi3=pi/0.75d0; ufpi=1.d0/fpi
facmy=(18.d0/pi2)**thrd

allocate(rxt(nieq),nxt(nieq),rmt(nieq),atpot(nx,nieq),atrho4pi(nx,nieq))
allocate(imax(nieq),qmt(nocc))

imax(1:nieq)=nx
call Poisson(atrho4pir2,atpot,imax)
do ir=1,nieq 
  do i=1,nx
    atrho4pi(i,ir)=atrho4pir2(i,ir)/(rx(i)*rx(i))
  enddo
enddo
!superposition of atomic charge densities.
!rc(i,j) =i'th coordinate of the j'th axis,
!rk(i,j) =i'th coordinate of the j'th atom,
!ncon(ir)=number of ir type shells,
!ia(j,ir)=atomic type in j'th shell,
!na(j,ir)=number of atoms in j'th shell,
!ad(j,ir)=distance to j'th shell.
call neighbor_shells(ia,na,ad,ncon,rc,rk,atrho4pi,voluc,baseuc)
imax(1:nieq)=(log(ad(2,1:nieq)*radfac/rx(1)))/dx+1.d0
imax(1:nieq)=min(imax(1:nieq),nx)
call write_AtRho(atrho4pi,imax)
call write_AtPot(atpot,imax)
call neighbor_sum(ia,na,ad,ncon,imax,atrho4pi,swrho4pi)
!A sphere of radius rxt containing Z-valieq electrons.
!Potential inside sphere from swrho4pi using Poisson's equation.
swpot=0.d0
do ir=1,nieq
  call ChargeAndPoisson(imax(ir),swrho4pi(1,ir),z(ir),valieq(ir),&
    swpot(1,ir),nxt(ir),rxt(ir))
  n=nxt(ir)
  atrho4pir2(1:n,ir)=swrho4pi(1:n,ir)*rx(1:n)*rx(1:n)
  !core level shift.
  if(cls(ir)/=0.d0) swpot(1:n,ir)=swpot(1:n,ir)+cls(ir)
enddo
call write_SWpot(swpot,nxt)
call write_SWrho(swrho4pi,nxt)
write(61,*)
write(61,906)'            atom ',(ir,ir=1,nieq)
write(61,906)' NearestNeighbor ',(ia(2,ir),ir=1,nieq)
write(61,908)'        Znucleus ',z(1:nieq)
write(61,908)'    ChgElectrons ',swpot(1,1:nieq)*rx(1)*0.5d0+valieq(1:nieq)
write(61,910)'   ChgRadius(AA) ',rxt(1:nieq)*bohr
write(61,908)'CLS (up, >0)(eV) ',cls(1:nieq)*rydb
906   format(a,i6,  7i7:  /(16x,8i7:))
908   format(a,f6.2,7f7.2:/(16x,8f7.2:))
910   format(a,f6.3,7f7.3:/(16x,8f7.3:))

allocate(mtpot(nx,nieq),ispot(nx,nieq))
!Fast-electron energy level and zeroth approximation muffin-tin radius.
!COMMENT: mtpot,ispot,xc are used in mt_optimisation,mt_err.
mtpot=swpot
ispot=swpot
xc='n'
call mt_optimisation(xc,ad,na,ia,ncon,vdum,swrho4pi)
swpot=mtpot

!interstitial charge density.
!qxmt=charge exterior to MT spheres
!volmt=volume of MT spheres,
!volcs=volume of charge shells,
!volis=volume of interstice,

qxmt=0.d0; volmt=0.d0; volcs=0.d0
do n=1,nieq
  qxmt=qxmt+neq(n)*&
          integral_egt(3,nxt(n),swrho4pi(1,n),rmt(n),rxt(n),dum)
  volmt=volmt+neq(n)*rmt(n)**3
  volcs=volcs+neq(n)*(rxt(n)**3-rmt(n)**3)
enddo
volmt=fpi3*volmt
volcs=fpi3*volcs
volis=voluc-volmt
rsfac=(volis/volcs)**thrd
write(61,932)'     VolUCell(AA**3) ',voluc*bohr3
write(61,932)'VolMuffinTins(AA**3) ',volmt*bohr3
write(61,932)'VolInterstice(AA**3) ',volis*bohr3
write(61,932)'    VolShells(AA**3) ',volcs*bohr3
write(61,932)'        rsShells(au) ',(fpi3*qxmt/volcs)**(-thrd)
write(61,932)'    rsInterstice(au) ',(fpi3*qxmt/volis)**(-thrd)
  932 format(a,f8.2)
!preparing exchange-correlation calculation:
!The Hedin-Lundqvist local density approximation.
!rs            -- radius of sphere accommodating one electron,
!                 (4pi/3)*rs**3=1/charge density,
!X, BETA, MYXC -- see
!L.Hedin and B.I.Lundqvist, J.Phys.C4,2064(1971),
!J. Rundgren Phys.Rev.B 68, 125405 (2003).
do n=1,nieq
  do i=1,nxt(n)
    rs(i,n)=(3.d0/swrho4pi(i,n))**thrd
    myx=facmy/rs(i,n)
    x=rs(i,n)/21.d0; beta=1.d0+0.7734d0*x*log(1.d0+1.d0/x)
    myxc(i,n)=-beta*myx
  enddo
enddo
return
contains

  subroutine poisson_equation(atrho4pir2,atpot,imax)
  !POISSON_EQUATION calculates Hartree potential from charge density.
  implicit none
  integer :: j,k,ir,imax(nieq),mx
  real*8  :: atrho4pir2(nx,nieq),atpot(nx,nieq)
  real*8,allocatable :: wa(:),wb(:)
  allocate(wa(nx),wb(nx))
  do ir=1,nieq
    mx=imax(ir)
    do j=1,mx
      k=mx-j+1
      atpot(j,ir)=atrho4pir2(k,ir)
      wa(j)=rx(k)*atrho4pir2(k,ir)
    enddo
    call integral_ugb(wa,wb)
    call integral_ugb(atpot(1,ir),wa)
    do j=1,mx
      k=mx-j+1
      atpot(j,ir)=2.0d0*(wa(k)-wb(k)/rx(j))
    enddo
  enddo
  return
  end subroutine poisson_equation

  subroutine integral_ugb(y,z)
  !INTEGRAL_UGB is uniform-grid integration by Bode's rule.
  !Ref.: Abramowitz and Stegun, sec. 25.4.
  implicit none
  integer :: i,j
  real*8  :: y(nx),z(nx),dxh,w1,w2
  real*8,parameter :: a7=7.d0/22.5d0, a32=32.d0/22.5d0, a12=12.d0/22.5d0
  !start by trapezoidal integration.
  dxh=dx*0.5d0
  z(1)=0.d0
  w1=z(1)
  do i=2,4
    w2=dxh*y(i)
    z(i)=z(i-1)+w1+w2
    w1=w2
  enddo
  !integration by Bode's rule.
  do i=5,8
    do j=i,nx,4
      z(j)=z(j-4)+dx*(a7*(y(j-4)+y(j))+a32*(y(j-3)+y(j-1))+a12*y(j-2))
    enddo
  enddo
  return
  end subroutine integral_ugb

  subroutine Poisson(atrho4pir2,atpot,imax)
 !solves Poissons equation for the spherical well.
  implicit none
  integer :: imax(nieq),n,i
  real*8  :: atrho4pir2(nx,nieq),atpot(nx,nieq),q(nx),s(nx),dh,w1,w2
 !trapezoidal rule with exponential grid of logarithmic increment dx.
  do ir=1,nieq
   !charge integral in the solution to Poisson's equation.
    n=imax(ir)
    dh=0.5d0*dx    
    q(1)=atrho4pir2(1,ir)*rx(1)/3.d0
    w1=dh*atrho4pir2(1,ir)*rx(1)
    do i=2,n
      w2=dh*atrho4pir2(i,ir)*rx(i)
      q(i)=q(i-1)+w1+w2
      w1=w2
    enddo
   !shell integral in the solution to Poisson's equation.
    s(1)=atrho4pir2(1,ir)/2.d0
    w1=dh*atrho4pir2(1,ir)
    do i=2,n
      w2=dh*atrho4pir2(i,ir)
      s(i)=s(i-1)+w1+w2
      w1=w2
    enddo
   !solution to Poisson's equation.    
    w2=s(n)
    atpot(1:n,ir)=2.d0*( (q(1:n)-z(ir))/rx(1:n)+w2-s(1:n) )
  enddo
  return
  end subroutine Poisson

  subroutine ChargeAndPoisson(imax,rho4pi,z,val,pot,nxt,rxt)
!CHARGE_RADIUS calculates the radius of an ion containing charge z and
!solves Poissons equation for the spherical well.
  implicit none
  integer :: imax,i,j,nxt
  real*8  :: rho4pi(nx),pot(nx),q(imax),s(imax),z,val,dh,rxt,&
    q1,q2,s1,s2,p,zv
!trapezoidal rule with exponential grid of logarithmic increment dx.
  dh=0.5d0*dx  
  q(1)=rho4pi(1)*rx(1)**3/3.d0
  q1=dh*rho4pi(1)*rx(1)**3
  s(1)=rho4pi(1)*rx(1)**2/2.d0
  s1=dh*rho4pi(1)*rx(1)**2
  do i=2,imax
   !charge integrated.
    q2=dh*rho4pi(i)*rx(i)**3
    q(i)=q(i-1)+q1+q2
    q1=q2
   !shell integral in the solution to Poisson's equation.    
    s2=dh*rho4pi(i)*rx(i)**2
    s(i)=s(i-1)+s1+s2
    s1=s2
  enddo
  s2=s(imax)
  pot(1:imax)=2.d0*((q(1:imax)-z)/rx(1:imax)+s2-s(1:imax)) 
 !radius where internal charge is equal to Z-valence.
  zv=z-val
  nxt=imax
  do i=1,imax
    if(q(i)>=zv)then
      nxt=i
      exit
    endif
  enddo
 !linear interpolation of j between q(nxt-1) and q(nxt).
  p=(zv-q(nxt-1))/(q(nxt)-q(nxt-1))  
 !first-order Taylor series for radius.
  rxt=rx(nxt-1)+p*(rx(nxt)-rx(nxt-1))
 !potential is shifted on the energy scale to zero maximum.
  q2=maxval(pot(1:nxt))
  pot(1:imax)=pot(1:imax)-q2
  return
  end subroutine ChargeAndPoisson
  
  subroutine write_AtRho(atrho4pi,imax)
  implicit none
  integer :: imax(nieq),i,ir
  real*8 :: atrho4pi(nx,nieq)
  select case(RhoPot)
  case('y')
    do ir=1,nieq
      open(10,file='rho.A'//id_at(ir),status='replace')
      do i=1,imax(ir)
        if(atrho4pi(i,ir)>5.d0) cycle      
        write(10,900) rx(i)*bohr,atrho4pi(i,ir)*ufpi
      enddo
    enddo
  endselect
  close(10)
  return
  900 format(1p2e14.6)
  910 format(i4,1pe14.6,0p)
  end subroutine write_AtRho
 
  subroutine write_SWrho(swrho4pi,nxt)
  implicit none
  integer :: nxt(nieq),i,ir
  real*8 :: swrho4pi(nx,nieq)
  select case(RhoPot)
  case('y')
    do ir=1,nieq   
      open(10,file='rho.SW'//id_at(ir),status='replace')
      do i=1,nxt(ir)
        if(swrho4pi(i,ir)>5.d0) cycle      
        write(10,900) rx(i)*bohr,swrho4pi(i,ir)*ufpi 
      enddo
    enddo
    close(10)
  endselect
  return
  900 format(1p4e14.6)
  910 format(i4,1pe14.6,0p)
  end subroutine write_SWrho

  subroutine write_AtPot(atpot,imax)
  implicit none
  integer :: imax(nieq)
  real*8 :: atpot(nx,nieq) 
  if(RhoPot=='y')then
    do ir=1,nieq
      open(10,file='pot.A'//id_at(ir),status='replace')
      do i=1,imax(ir)
        if(atpot(i,ir)<-5.d0) cycle
        write(10,900) rx(i)*bohr,atpot(i,ir)
      enddo
    enddo
  endif
  close(10)
  return
  900 format(1p3e14.6)
  910 format(i4,1pe14.6,0p)
  end subroutine write_AtPot
    
  subroutine write_SWpot(swpot,nxt)
  implicit none
  integer :: nxt(nieq),i,ir
  real*8 :: swpot(nx,nieq)
  select case(RhoPot)
  case('y')
    do ir=1,nieq   
      open(10,file='pot.SW'//id_at(ir),status='replace')
      do i=1,nxt(ir)
        if(swpot(i,ir)<-5.d0) cycle      
        write(10,900) rx(i)*bohr,swpot(i,ir) 
      enddo
      open(10,file='rad.SW'//id_at(ir),status='replace')
      write(10,900) rxt(ir)*bohr,-1.d0
      write(10,900) rxt(ir)*bohr, 1.d0
      close(10)
    enddo
  endselect
  return
  900 format(1p3e14.6)
  910 format(i4,1pe14.6,0p)
  end subroutine write_SWpot
  
  real*8 function xxxrho4pi(i,pot,n)
  implicit none
  integer :: i,n
  real* 8 :: pot(n),didr,dpdi,d2pdi2
  didr=1.d0/(dx*rx(i))
  !three-point differentiation with respect to log. grid subscript.
  if(i<=n-1)then
    dpdi=0.5d0*(pot(i+1)-pot(i-1))
    d2pdi2=pot(i+1)-2.d0*pot(i)+pot(i-1)
  else
    dpdi=1.5d0*pot(i)-2.d0*pot(i-1)+0.5d0*pot(i-2)
    d2pdi2=pot(i)-2.d0*pot(i-1)+pot(i-2)
  endif
  !Poisson's differential operator, 0.5 is V(Ry) back to rho.
  xxxrho4pi=-0.5d0*(d2pdi2+dx*dpdi)*didr**2
  return
  end function xxxrho4pi
  
end subroutine structure_I
!-----------------------------------------------------------------------
subroutine neighbor_shells(ia,na,ad,ncon,rc,rk,atrho4pi,voluc,baseuc)
!NEIGHBOR_SHELLS is nearest neighbor data for atoms in a crystal
!structure.
!Ref.: A. Barbieri and M.A. Van Hove, Phase shift package,
!      http://electron.lbl.gov/software.
!Ref.: DLPHASE v1.1, CCP3 Library code,
!      http://www.cse.clrc.ac.uk/Activity.

!New in this code is ARMAX(1:nieq). These radii sort out the necessary
!and sufficient lattice points for getting neighbor contributions of
!a given accuracy to charge density and potential.
!The main part of the subroutine is a fortran90 editon of the Ref.

!input:
!rc(i,j) =the i'th coordinate of the j'th axis of the unit cell,
!rk(i,j) =the i'th coordinate of the j'th atom in the unit cell,
!neq(ir) =the number of equivalent atoms of type ir.
!output:
!ia(j,ir)=the type of atoms in the j'th neighbor shell,
!na(j,ir)=the number of atoms in the j'th shell,
!ad(j,ir)=the radius of the j'th shell,
!ncon(ir)=no. of ir type shells, limited by ARMAX(ir) and NSHELL.
use param,only : BulkOrSlab,bohr,deltasup
use dimen,only : nieq,neq,nshell,ncell,nx
use radii,only : rx
implicit none
integer :: ia(nshell,nieq),na(nshell,nieq),ncon(nieq),j,jr,&
           jjr,k,kr,ic,iic,jjc,jc,ir,nc,kc,i,it,itz,limc,jx,jy,jz
real*8  :: rc(3,3),rk(3,ncell),rj(3),ad(nshell,nieq),&
           atrho4pi(nx,nieq),rmax,rcmin,dr,r,voluc,baseuc,rchg
!charge density >DELTASUP determines rmax.
rmax=0.d0
do ir=1,nieq
  i=nx+1
  do
    i=i-1
    if(abs(atrho4pi(i,ir))>deltasup)exit
  enddo
  rmax=max(rmax,rx(i))
enddo
write(61,901) deltasup,rmax*bohr
901  format('atomic charge density >',1p,e6.0,0p,&
           &' implies rmax(AA)=',f5.2)
           
!find sufficiently large loop index for translation vectors RC.
select case(BulkOrSlab)
 case('b')
  rcmin=1.d+33
  do i=1,3 !three dimensions.
    rcmin=min(rcmin,abs_value(rc(1,i),rc(2,i),rc(3,i)))
  enddo
  it=(rmax+rmax)/rcmin+1.d0
  itz=it
 case('s')
  rcmin=1.d+33
  do i=1,2 !two dimensions.
    rcmin=min(rcmin,abs_value(rc(1,i),rc(2,i),rc(3,i)))
  enddo
  it=(rmax+rmax)/rcmin+1.d0
  itz=0
end select

!search over adjacent unit cells to include nshell nearest neighbors.
ia=0
na=0
ad=1.d+33
do jx=-it,it
  do jy=-it,it
    do jz=-itz,itz
      do j=1,3
        rj(j)=jx*rc(j,1)+jy*rc(j,2)+jz*rc(j,3)
      enddo
      !rj is current unit cell origin.  for each atom in this unit cell,
      !find displacement r from kr-type atom in basic unit cell.
      j=0
      do jr=1,nieq
        do jjr=1,neq(jr)
          j=j+1
          k=1
          do kr=1,nieq
            r=abs_value(rj(1)+rk(1,j)-rk(1,k),&
                        rj(2)+rk(2,j)-rk(2,k),&
                        rj(3)+rk(3,j)-rk(3,k))
            !fixed limit for r.
            if(r>rx(nx))then !16.d0)then
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
              dr=r-ad(ic,kr); if(abs(dr)<1.d-04) dr=0.0d0
              if(dr>0.d0)then
                cycle
              elseif(dr==0.d0)then
                if(ia(ic,kr)/=jr)then
                  cycle
                else
                  na(ic,kr)=na(ic,kr)+1
                  k=k+neq(kr)
                  exit
                endif
              else
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
    if(na(ic,ir)>0) ncon(ir)=ncon(ir)+1
  enddo
enddo
!display.
do ir=1,nieq
  j=0
  do i=2,ncon(ir)
    j=j+na(i,ir)  !j counts the neighbors about atom ir.
  enddo
  write(61,209) ir,ncon(ir)-1,j
  nc=ncon(ir)
  ic=(nc-1)/8+1
  kc=0
  do i=1,ic
    jc=kc+1
    kc=min(nc,kc+8)
    write(61,210) (ad(j,ir)*bohr,j=jc,kc)
    write(61,211) (na(j,ir),j=jc,kc)
    write(61,212) (ia(j,ir),j=jc,kc)
  enddo
enddo
!volume of unit cell is completed, in case 's' or baseuc/=0,
!by the spill-out of atomic volume at the first atomic plane.
if(baseuc/=0.d0)then
  dr=0.d0; nc=0
  do ir=1,nieq
    dr=dr+neq(ir)*ad(2,ir)
    nc=nc+neq(ir)
  enddo
  dr=dr/nc
  voluc=voluc+0.25d0*dr*baseuc
  !the spill-out is taken to be half of (one atomic radius=0.5*dr)
endif
return
209   format(/'type',i3,', ',i3,' distance vs type neighbor shells, '&
            &,i3' neighbors')
210   format('dist.(AA)',f7.4,10f8.4)
211   format('number   ',i4,3x,10(i5,3x))
212   format('type     ',i4,3x,10(i5,3x))
contains

  real*8 function abs_value(a1,a2,a3)
  implicit none
  real*8 :: a1,a2,a3
  abs_value=sqrt(a1*a1+a2*a2+a3*a3)
  end function abs_value
end subroutine neighbor_shells
!-----------------------------------------------------------------------
subroutine neighbor_sum(ia,na,ad,ncon,imax,elem,xtal)
!NEIGHBOR_SUM superposes charge-density (or potential) contributions
!from neighbor atoms. The present code updates the routine 'sumax' of
!the following program packages,
!Ref.: A. Barbieri and M.A. Van Hove, Phase shift package,
!      http://electron.lbl.gov/software.
!Ref.: DLPHASE v1.1, CCP3 Library code,
!      http://www.cse.clrc.ac.uk/activity.
!Ref.: T.L. Loucks, Augmented plane wave method (Benjamin, New York,
!      1967), eqs. 3.22,3.26,3.28.

!elem(i,ir) = atrho4pi or atpot for the element ir,
!xtal(i,ir) = swrho4pi or swpot for the element ir,
!     ncon  = number of shells included,
!     ia(j) = atomic type in j'th shell,
!     na(j) = number of atoms in j'th shell,
!     ad(j) = distance to j'th shell,
!     imax  = summation up to rx(imax).
use dimen,only : nieq,nshell,nx
use radii,only : rx,dx
implicit none
integer :: ia(nshell,nieq),na(nshell,nieq),ncon(nieq),ja,ic,i,ir,&
  imax(nieq),jmax,m,j1,j2,j
real*8  :: xtal(nx,nieq),elem(nx,nieq),ad(nshell,nieq),x1,x2,sum,dum,&
  integral_egt,p1,p2,elem1,elem2,w1,w2,hh
  
hh=0.5d0*dx
do ir=1,nieq
  jmax=imax(ir)
  xtal(1:jmax,ir)=elem(1:jmax,ir)
  do ja=2,ncon(ir)
    ic=ia(ja,ir)
    do i=1,jmax
      x1=min(abs(rx(i)-ad(ja,ir)),rx(nx)) !rx(nx) is correct.
      x2=min(rx(i)+ad(ja,ir),rx(nx))
      if(x1<x2)then
        sum=integral_egt(2,nx,elem(1,ic),x1,x2,dum)
        xtal(i,ir)=xtal(i,ir)+na(ja,ir)*sum*0.5d0/(ad(ja,ir)*rx(i))
      endif
    enddo
  enddo
enddo
return
end subroutine neighbor_sum
!-----------------------------------------------------------------------
subroutine mt_optimisation(xc,ad,na,ia,ncon,V0mVfe,swrho4pi)
!MT_OPTIMISATION calculates muffin-tin radii and interstitial potential.
!The MT floor of a particular atom is the average superposition
!potential in the shell between the muffin-tin radius and the charge
!radius. The interstitial potential is an average of the atomic
!MT floors, each atom weighted with its occupation in the structure.

!SELECTION 'Norman' or 'n':
!Consider a particular atom A and its nearest neighbors B. For each
!couple A,B touching spheres are designed with radii proportional to
!the Norman atomic radii of A and B. The set of A radii thus obtained
!are averaged and taken to be the MT radius of A.
!
!The potential steps at the MT radii are generally different.

!SELECTION 'intersection' or 'i':
!Consider a particular atom A and its nearest neighbors B. For each
!couple A,B touching spheres are designed with radii determined by
!the point where the two potentials intersect on the line between the
!nuclei. The set of A radii thus obtained are averaged and taken to be
!the MT radius of A.
!
!The potential steps at the MT radii are generally different.

!SELECTION 'optimum' or 'o':
!Selection 'optimum' starts from selection 'intersection'. The MT radii
!are optimised in such a way that each particular MT floor is equal to
!the interstitial potential.
!
!There is a common potential step at all MT radii.

!xc='no':
!The interstitial static potential, shifted to zero, is the fast-
!electron level of the surface slab.

!xc='yes':
!The interstitial total potential (static + exchange-correlation) is
!the inner potential V0mVfe, and the MT cores are shifted to zero step at
!the MT radius for all SELECTIONs above.
use param,only : pi,deltarmt,bohr,rydb,&
  NMselect,NMlambda,NMeps,NMepsit,NMiter
use dimen,only : nshell,nx,nieq,neq,nocc
use radii,only : rxt,nxt,rmt,rx,dx,qmt
use poten,only : mtpot,ispot
use error
implicit none
character :: xc*1,MTradius*1
integer :: na(nshell,nieq),ia(nshell,nieq),ncon(nieq),i1,i2,i,lo,hi,&
           iter,count,ns,nmt,neqsum
real*8  :: ad(nshell,nieq),r1,r2,r12,err1,err2,h,vol,v0,V0mVfe,lambda,&
  eps,norm,avvrmt,interpol_egl3,integral_egt,swrho4pi(nx,nieq),dum
real*8,allocatable :: v(:,:),y(:),vrmt(:),tmp(:),rmt0(:)
allocate(v(0:nocc,nocc),y(0:nocc),vrmt(nieq),tmp(nieq),rmt0(nieq))
MTradius='o' !'optimized'/'intersection'/'norman'

select case(MTradius)
 case('n','N')
  !Norman radius: 
  !atom i's MT radius is proportional to its charge radius together with
  !the charge radius(radii) of the nearest neighbor(s). With several
  !neighbors ia(ns,i) the MT radius is a weighted average.
  !COMMENT 1:nieq charge radii are used to calculate 1:nocc MT radii.
  do i=1,nieq
    rmt(i)=0.d0; nmt=0
    do ns=2,ncon(i)
      if(abs(ad(2,i)-ad(ns,i))>deltarmt+0.001d0)exit
      rmt(i)=rmt(i)+na(ns,i)*ad(ns,i)*rxt(i)/(rxt(i)+rxt(ia(ns,i)))
      nmt=nmt+na(ns,i)
    enddo
    rmt(i)=rmt(i)/nmt
  enddo
 case('i','o')
  !Intersection radius:
  !atom i's MT radius from the point(s) where its potential intersects
  !the potential(s) of the nearest neighbor(s). With several neighbors
  !ia(ns,i) the MT radius is a weighted average.
  !COMMENT 1:nieq potentials are used to calculate 1:nocc MT radii.
  do i=1,nieq
    rmt(i)=0.d0; nmt=0
    do ns=2,ncon(i)
      if(abs(ad(2,i)-ad(ns,i))>deltarmt+0.001d0)exit
      rmt(i)=rmt(i)+na(ns,i)*&
                    intersect_pot(mtpot,ia(1,i),ia(ns,i),ad(ns,i))
      nmt=nmt+na(ns,i)
    enddo
    rmt(i)=rmt(i)/nmt
  enddo
end select
write(61,912) rmt*bohr
!A single atomic type in the structure.
select case(MTradius)
 case('o')
  if(nieq==1)then
    MTradius='i'
    write(61,'(/a/)')"nieq=1 => SelectRadius='i'."
  endif
end select 

select case(MTradius)
 case('o')
  !optimisation of MT radii starting from intersection radii.
  !The downhill simplex method of Nelder and Mead.
  rmt0=rmt
  if(NMselect=='*')then 
    lambda=1.d0+NMlambda
    eps=NMeps
    do iter=1,NMiter
      do i=0,nocc
        v(i,1:nocc)=rmt(1:nocc)*0.5d0+rmt0(1:nocc)*0.5d0
      enddo
      do i=1,nocc
        v(i,i)=v(i,i)*lambda
      enddo
      do i=0,nocc
        tmp(1:nocc)=v(i,1:nocc)
        y(i)=mt_err(tmp)
      enddo
      call simplex(mt_err,v,y,nocc,eps,lo,hi,norm,count)
      rmt0=rmt
      rmt(1:nocc)=v(lo,1:nocc)
      write(61,902) rmt(1:nocc)*bohr
      write(61,904) count,norm*bohr
      eps=eps*NMepsit
      lambda=sqrt(lambda)
    enddo
  elseif(NMselect=='+')then
    lambda=NMlambda
    eps=NMeps
    do iter=1,NMiter
      do i=0,nocc
        v(i,1:nocc)=rmt(1:nocc)*0.5d0+rmt0(1:nocc)*0.5d0
      enddo
      do i=1,nocc
        v(i,i)=v(i,i)+lambda
      enddo
      do i=0,nocc
        tmp(1:nocc)=v(i,1:nocc)
        y(i)=mt_err(tmp)
      enddo
      call simplex(mt_err,v,y,nocc,eps,lo,hi,norm,count)
      rmt0=rmt
      rmt(1:nocc)=v(lo,1:nocc)      
      write(61,902) rmt(1:nocc)*bohr
      write(61,904) count,norm*bohr
      eps=eps*NMepsit
      lambda=lambda*0.5d0
    enddo
  else
    write(61,'(3a)')'NMselect=',NMselect,', in error, stop'; stop
  endif
end select
  do i=1,nocc
    qmt(i)=integral_egt(3,nx,swrho4pi(1,i),rx(1),rmt(i),dum)
  enddo
  write(61,903) -qmt(1:nocc)

select case(xc)
 case('n')
  !Interstitial potential is shifted to zero for i=1:nieq.
  !COMMENT mtpot(1:nx,1:nieq) needed for 1:nocc intersection radii.
  v0=interstitial_pot(mtpot)
  do i=1,nieq
    mtpot(1:nxt(i),i)=mtpot(1:nxt(i),i)-v0
  enddo
  v0=0.
 case('y')
  !'mtpot at r=rmt' is shifted to zero.
  v0=interstitial_pot(ispot)
  avvrmt=0.d0; neqsum=0
  do i=1,nocc
    vrmt(i)=interpol_egl3(rmt(i),mtpot(1,i))
    mtpot(1:nxt(i),i)=mtpot(1:nxt(i),i)-vrmt(i)
    avvrmt=avvrmt+neq(i)*vrmt(i)
    neqsum=neqsum+neq(i)
  enddo
  eps=sum((vrmt(1:nocc)-avvrmt)**2)
  eps=sqrt(eps/neqsum)
  write(61,914) v0*rydb
  write(61,915) (vrmt(1:nocc)-v0)*rydb
  write(61,916)
  if(nocc>1)then
    avvrmt=avvrmt/neqsum
    eps=sum((vrmt(1:nocc)-avvrmt)**2); eps=sqrt(eps/neqsum)
    write(61,917) eps*rydb
  endif
end select
V0mVfe=v0
return
902   format('  optimum rMT AA ',f6.3,7f7.3:/(16x,8f7.3:))
903   format('  ElCharge in MT ',f6.2,7f7.2:/(16x,8f7.2:))
904   format(' iter,rmt_rms AA ',i4,2x,1p2e7.0)
910   format('   Norman rMT AA ',f6.3,7f7.3:/(16x,8f7.3:))
912   format(' intersec rMT AA ',f6.3,7f7.3:/(16x,8f7.3:))
913   format('  optimum rMT AA ',f6.3,7f7.3:/(16x,8f7.3:))
914   format('       V0-Vfe eV ',f6.2) 
915   format('      Vrmt-V0 eV ',f6.2,7f7.2:/(16x,8f7.2:))
916   format('      Vrmt-V0->0 ')
917   format(' rms(Vrmt-V0) eV ',1pe6.0)
contains

  real*8 function interstitial_pot(pot)
  implicit none
  real*8  :: pot(nx,nieq),v0,vol,dum,integral_egt
  real*8,parameter :: thrd=1.d0/3.d0
  v0=0.d0; vol=0.d0
  do i=1,nocc
    v0=v0+neq(i)*&
       integral_egt(3,nxt(i),pot(1,i),rmt(i),rxt(i),dum)
    vol=vol+neq(i)*thrd*(rxt(i)**3-rmt(i)**3)
  enddo
  interstitial_pot=v0/vol
  return
  end function interstitial_pot
  
  real*8 function intersect_pot(pot,i1,i2,r12)
  implicit none
  integer :: i1,i2,ii,j
  real*8  :: pot(nx,nieq),r12,a,b
  !intersection by grid search, start from half interatomic distance.
  j=nint((log(r12*0.5d0/rx(1)))/dx+1.0d0)
  a=interpol_egl3(r12-rx(j),pot(1,i2))-pot(j,i1)
  if(abs(a)<1.d-04)then
    intersect_pot=rx(j)
    return
  elseif(a>0.d0)then
    do
      if(a<=0.d0)exit
      b=a
      j=j+1
      a=interpol_egl3(r12-rx(j),pot(1,i2))-pot(j,i1)
    enddo
    intersect_pot=rx(j)-(rx(j)-rx(j-1))*abs(a)/(abs(a)+abs(b))
    return
  elseif(a<0.d0)then
    do
      if(a>=0.d0)exit
      b=a
      j=j-1
      a=interpol_egl3(r12-rx(j),pot(1,i2))-pot(j,i1)
    enddo
    intersect_pot=rx(j)+(rx(j+1)-rx(j))*abs(a)/(abs(a)+abs(b))
    return
  endif
  end function intersect_pot
end subroutine mt_optimisation
!-----------------------------------------------------------------------
subroutine structure_II&
  (EmVfe,V0mVfe,swrho4pi,swpot,rs,rsfac,myxc,ad,na,ia,ncon)
use param,only : nsp,nsr,sp,sr,sdat,rydb
use dimen,only : nx,nieq,nshell
use radii,only : rx,dx,rmt,nxt
use poten,only : mtpot,ispot
implicit none
character :: xc*1,XCrenorm*1
integer :: na(nshell,nieq),ia(nshell,nieq),ncon(nieq),n,i,imt
real*8 :: swrho4pi(nx,nieq),swpot(nx,nieq),rs(nx,nieq),myxc(nx,nieq),&
  ad(nshell,nieq),rsfac,EmVfe,V0mVfe,p2,sep,se1
real*8,allocatable :: f(:),g(:)
allocate(f(nx),g(nx))
!muffin-tin potential MTPOT and interstitial potential ISPOT.
XCrenorm='y'
do n=1,nieq
  imt=(log(rmt(n)/rx(1)))/dx+1.d0
  do i=1,nxt(n)
    f(i)=xcpot(i,EmVfe,rs(i,n),myxc(i,n),p2,sep,se1)
!    if(i>=imt)then
!      write(89,'(i5,1p4e12.4)') i,f(i),myxc(i,n)/se1
!      write(91,*) i,p2,rs(i,n)
!    endif
  enddo
  call xc_smooth(f,g)
  do i=1,nxt(n)
    mtpot(i,n)=swpot(i,n)+f(i)
!    write(90,*) i,f(i)
    if(XCrenorm=='y')then
      ispot(i,n)=swpot(i,n)+f(i)/rsfac
    else
      ispot(i,n)=swpot(i,n)+f(i)
    endif
  enddo
enddo
xc='y'
call mt_optimisation(xc,ad,na,ia,ncon,V0mVfe,swrho4pi)
deallocate(f,g)
return
contains

  real*8 function xcpot(i,eryd,rs,myxc,p2,sep,se1)
  !XCPOT is excited-state exchange-correlation potential.
  !Ref.: L.Hedin and B.I.Lundqvist, J.Phys.C 4,2064 (1971).
  !Ref.: R.E.Watson, J.F.Herbst, L.Hodges, B.I.Lundqvist, and
  !      J.W.Wilkins, Phys.Rev.B13, 1463 (1976), Appendix B.
  !Ref.: J.Neve, J.Rundgren, and P.Westrin, J.Phys.C 15, 4391 (1982).
  !pF=fermi momentum,
  !p =electron momentum in units of pF, used in self_energy table,
  !se(p,rs)=self energy in units of pF**2,
  !myxc    =ground-state exchange-correlation potential,
  !myxc*(se(p,rs)/se(1,rs))=excited-state exch.-corr. potential.
  !Hedin-Lundqvist equation:
  !  (p*pF)**2+myxc*(se(p,rs)/se(1,rs)=eryd+pF**2+myxc
  !below rewritten in the form
  !  f=(p-1)*(p+1)+[myxc*(se(p,rs)/se(1,rs)-myxc-eryd]/pF**2
  implicit none
  integer :: ie,n,i
  real*8  :: se1,eryd,eeff,rs,myxc,pF,p1,p2,f1,f2,fp1,h,vxc,a,sep
  real*8,parameter :: pi9_4thrd=1.919158292677504d0
 !Newton-Raphson method using derivative.
  if(i==3000)write(*,*)'a'
  pF=pi9_4thrd/rs
  a=eryd+pF*pF
  p1=sqrt(a+1.d0)/pF
  p2=sqrt(a)/pF
  se1=self_energy(1.d0,rs)
  vxc=myxc*self_energy(p1,rs)/se1
  a=1.d0/(pF*pF)
  f1=(p1-1.d0)*(p1+1.d0)+(vxc-myxc-eryd)*a
  vxc=myxc*self_energy(p2,rs)/se1
  f2=(p2-1.d0)*(p2+1.d0)+(vxc-myxc-eryd)*a
  fp1=(f2-f1)/(p2-p1)
  ie=0
  if(i==3000)write(*,*)'b'
  do
    ie=ie+1
    if(i==3000)write(*,*)'ie=',ie
    h=-f1/fp1
    p2=p1+h
    vxc=myxc*self_energy(p2,rs)/se1
    f2=(p2-1.d0)*(p2+1.d0)+(vxc-myxc-eryd)*a
    if(abs(f2)<1.d-10)exit
    fp1=(f2-f1)/h
    p1=p2
    f1=f2
  enddo
  sep=self_energy(p2,rs)
  xcpot=vxc
  return
  end function xcpot

  real*8 function self_energy(ps,rs)
  !VXCDATA extracts data from the Vxc tables supplied by Hedin-Lundqvist
  !and Sernelius.
  implicit none
  real*8 :: ps,rs
  if(rs>sr(nsr))then
    write(61,*)'VXCDATA stop: rs>6., new code needed.'; stop
  endif
  if(ps<=sp(nsp))then
    self_energy=interpol_tl(ps,rs)
  else
    self_energy=extrapol_tl(ps,rs)
  endif
  return
  end function self_energy
  
  real*8 function extrapol_tl(ps,rs)
  !EXTRAPOL_TL=inverse-wavenumber extrapolation from given Vxc tables.
  implicit none
  real*8,parameter :: ps1=2.80d0,ps2=3.00d0
  real*8 :: ps,rs,v1,v2,g,b,a
  v1=interpol_tl(ps1,rs)
  v2=interpol_tl(ps2,rs)      !value at table margin.
  g=abs(v2*(ps2-ps1)/(v2-v1)) !gradient at table margin.
  b=ps2*(g-ps2)
  a=v2*sqrt(ps2*g)
  extrapol_tl=a/sqrt(ps*ps+b)
  return
  end function extrapol_tl
  
  real*8 function interpol_tl(ps,rs)
  !INTERPOL_TL=triangulated linear interpolation in given Vxc tables.
  !Ref.: L.F.Shampine, R.C.Allen, and S.Pruess, Fundamentals of
  !      numerical computing (Wiley,1997).
  implicit none
  integer :: k,k1,k2,l1,l2,l,i
  real*8  :: ps,rs,a1,a2,a3,p21,r21,s,ta,a,b
  !rectangle containing ps,rs.
  l1=1; l2=nsr
  do
    l=(l1+l2)/2; if(sr(l)>rs)then; l2=l; else; l1=l; endif
    if(l2-l1==1)exit
  enddo
  k1=1; k2=nsp
  do
    k=(k1+k2)/2; if(sp(k)>ps)then; k2=k; else; k1=k; endif
    if(k2-k1==1)exit
  enddo
  !triangles: rectangle with diagonal.
  p21=sp(k2)-sp(k1); r21=sr(l2)-sr(l1); s=r21/p21; ta=r21*p21 
  if(rs-sr(l1)>=s*(ps-sp(k1)))then
   !upper triangle.
    a1=phi(ps,rs,sp(k2),sr(l2),sp(k1),sr(l2),ta)
    a2=phi(ps,rs,sp(k1),sr(l2),sp(k1),sr(l1),ta)
    a3=phi(ps,rs,sp(k1),sr(l1),sp(k2),sr(l2),ta)
    a=a1*sdat(k1,l1)*sr(l1)+a2*sdat(k2,l2)*sr(l2)+a3*sdat(k1,l2)*sr(l2)
  else
   !lower triangle.
    a1=phi(ps,rs,sp(k2),sr(l1),sp(k2),sr(l2),ta)
    a2=phi(ps,rs,sp(k2),sr(l2),sp(k1),sr(l1),ta)
    a3=phi(ps,rs,sp(k1),sr(l1),sp(k2),sr(l1),ta)
    a=a1*sdat(k1,l1)*sr(l1)+a2*sdat(k2,l1)*sr(l1)+a3*sdat(k2,l2)*sr(l2)
  endif
  !triangles: rectangle with antidiagonal.
  s=-s
  if(rs-sr(l2)>=s*(ps-sp(k1)))then
   !upper triangle.
    a1=phi(ps,rs,sp(k2),sr(l1),sp(k2),sr(l2),ta)
    a2=phi(ps,rs,sp(k2),sr(l2),sp(k1),sr(l2),ta)
    a3=phi(ps,rs,sp(k1),sr(l2),sp(k2),sr(l1),ta)
    b=a1*sdat(k1,l2)*sr(l2)+a2*sdat(k2,l1)*sr(l1)+a3*sdat(k2,l2)*sr(l2)
  else
   !lower triangle.
    a1=phi(ps,rs,sp(k1),sr(l1),sp(k2),sr(l1),ta)
    a2=phi(ps,rs,sp(k2),sr(l1),sp(k1),sr(l2),ta)
    a3=phi(ps,rs,sp(k1),sr(l2),sp(k1),sr(l1),ta) 
    b=a1*sdat(k1,l2)*sr(l2)+a2*sdat(k1,l1)*sr(l1)+a3*sdat(k2,l1)*sr(l1)
  endif
  interpol_tl=0.5d0*(a+b)
  return
  end function interpol_tl
  
  real*8 function phi(x,y,xa,ya,xb,yb,ta)
  real*8 :: x,y,xa,ya,xb,yb,ta
  phi=(xa*yb-xb*ya+(ya-yb)*x+(xb-xa)*y)/ta
  end function phi
  
  subroutine xc_smooth(f,g)
 !binomial smoothing of second order, with coefficients 1,2,1;
 !when executed many times, it converges to the normal distribution.
  implicit none
  integer :: i,j
  real*8 :: f(nx),g(nx)
  real*8,parameter :: thrd=1.d0/3.d0
  do j=1,32
   !f goes to g.
    g(1)=(f(2)+f(1)+f(1))*thrd
    do i=2,nxt(n)-1
      g(i)=(f(i-1)+f(i+1)+f(i)+f(i))*0.25d0
    enddo
    i=nxt(n)
    g(i)=(f(i-1)+f(i)+f(i))*thrd
   !g goes back to f.
    f(1)=(g(2)+g(1)+g(1))*thrd
    do i=2,nxt(n)-1
      f(i)=(g(i-1)+g(i+1)+g(i)+g(i))*0.25d0
    enddo
    i=nxt(n)
    f(i)=(g(i-1)+g(i)+g(i))*thrd
  enddo
  return
  end subroutine xc_smooth
end subroutine structure_II
!-----------------------------------------------------------------------
subroutine scattering&
  (swrho4pi,swpot,rs,rsfac,myxc,z,ad,na,ia,ncon,valieq,title)
!SCATTERING calculates electron-atom scattering phase shifts on an
!energy scale referred to an energy-dependent inner potential, equal to
!Hedin-Lundqvist's local-density-approximation exchange-correlation
!potential.
!Ref.: J.Neve, J. Rundgren, and P. Westrin, J. Phys. C 15, 4391 (1982).

!    Vfe=fast-electron level
!       =zero point of the energy scale,
!      E=Eprimary-Vfe
!  V0(E)=Vinterstitial-Vfe
!       =energy-dependent inner potential,
!      q=sqrt(E-V0(E))
!       =electron wave number in the crystal,
!psr,ps =phaseshift versus q and orbital quantum number l,
!psrm,psrp=phaseshift versus q and relativistic quantum number kappa.
use param,only : Output,RelOrNonrel,SpinPhaseShift,PSvsE,pi,bohr,rydb,&
  cm2,cm2tmp,wave_start,eta_cutoff,ne,eev,id_at,relerr1,relerr2
use dimen,only : nshell,nieq,neq,nocc,nx,lmax
use radii,only : rx,rxt,nxt,rmt,dx,qmt
use scatt,only : iatom,emv0,q
use bessl,only : bess
use poten,only : mtpot,ispot,v0tab
implicit none
character :: xc*1,idl*3
character(len=*) :: title
logical :: ex
integer :: na(nshell,nieq),ia(nshell,nieq),ncon(nieq),i,n,ie,i1,i2,&
           imt,l,lmin,lmaxq,je,m,lx
integer,allocatable :: lxq(:,:)
real*8  :: z(nieq),valieq(nieq),ad(nshell,nieq),swrho4pi(nx,nieq),&
           swpot(nx,nieq),rs(nx,nieq),myxc(nx,nieq),v0r(4),EmVfe,&
           V0mVfe,h,qr,d,xc_pot,rsfac
real*8,allocatable,dimension(:)     :: bj,by
real*8,allocatable,dimension(:,:,:) :: psr,psrm,psrp,psn

allocate(v0tab(ne))

write(61,901) ne
write(61,902) wave_start,eta_cutoff
write(61,903) relerr1,relerr2
901   format(/'SCATTERING: ne=',i3)
902   format(12x,'wave_start=',1pe8.1,', eta_cutoff=',e8.1)
903   format(12x,'relerr1,2=',1p2e8.1)
if(Output=='p'.or.Output=='d')then
  open(10,file='RMTvsE',status='replace')
  open(12,file='QMTvsE',status='replace')  
endif
if(Output=='d')then
  open(20,file='LvsKRstart')
  open(21,file='LMAXvsKRmt',status='replace')
endif

!loop over energy grid:
!calculation of energy-dependent potential and phase shifts.
lx=0
do ie=1,ne
  write(61,'(18("----"))')
  write(61,'(a,i3,1x,f6.0)')'ie,E=',ie,eev(ie)
  write(*,'(i3,1x,f6.0)') ie,eev(ie)
  EmVfe=eev(ie)/rydb
  
  call structure_II&
    (EmVfe,V0mVfe,swrho4pi,swpot,rs,rsfac,myxc,ad,na,ia,ncon)
  
  !potential close to MT radius, for illustration in CPC paper.
  if(ie==-1.or.ie==-35)then
    if(ie== 1)open(99,file='mtpot_lo')
    if(ie==35)open(99,file='mtpot_hi')
    imt=(log(rmt(1)/rx(1)))/dx+1.d0
    h=rmt(1)*rmt(1)
    do i=imt-15,imt
      write(99,909) rx(i)-rmt(1),mtpot(i,1)*h
    enddo
    write(99,909) 0.d0,0.d0
  endif
  if(Output=='p'.or.Output=='d')then
    write(10,910) eev(ie),rmt(1:nocc)*bohr
    write(12,912) eev(ie),-qmt(1:nocc)    
  endif
  v0tab(ie)=V0mVfe
909   format(1p2e14.6)
910   format(f6.0,24f7.4)
912   format(f6.0,24f8.3)
  
  !allocate phase shift arrays corresponding to eta_cutoff.
  if(ie==1)then
    select case(RelOrNonrel)
      case('r'); cm2tmp=cm2
      case('n'); cm2tmp=0.d0
    endselect
    emv0=eev(ne)/rydb-V0mVfe
    qr=sqrt(emv0*(1.d0+emv0*cm2tmp))*maxval(rmt(1:nocc))
    l=min(10+int(1.07d0*qr),lmax) !lmax is input value.
    lmax=min(l,lmax)              !lmax can be put too large.
    allocate(lxq(ne,nocc))
    allocate(psr(ne,0:l,nocc),psn(ne,0:l,nocc))
    allocate(psrm(ne,0:l,nocc),psrp(ne,0:l,nocc))
    psr=0.d0; psn=0.d0; psrm=0.d0; psrp=0.d0
    write(61,'(12x,a,i3)')'lmax ',lmax
  endif
  
  !calculation of phase shifts, mtpot->total potential times radius.
  do n=1,nocc
    mtpot(1:nxt(n),n)=mtpot(1:nxt(n),n)*rx(1:nxt(n))
  enddo
  emv0=EmVfe-V0mVfe
  select case(RelOrNonrel)
  case('r')
    cm2tmp=cm2
    q=sqrt(emv0*(1.d0+emv0*cm2tmp))
    do iatom=1,nocc
      qr=q*rmt(iatom)
      lmaxq=min(10+int(1.07d0*qr),lmax)
      allocate(bj(0:lmaxq+1),by(0:lmaxq+1),bess(0:lmaxq+1))
      call sph_bessel_j(lmaxq+1,qr,bj)
      call sph_bessel_y(lmaxq+1,qr,by)
      call phaseshift('-',ie,lmaxq,lx,bj,by,psrm)
      call phaseshift('+',ie,lmaxq,lx,bj,by,psrp)
      lxq(ie,iatom)=lx
      deallocate(bj,by,bess)
    enddo
  case('n')
    cm2tmp=0.d0 
    q=sqrt(emv0)
    do iatom=1,nocc
      qr=q*rmt(iatom)
      lmaxq=min(10+int(1.07d0*qr),lmax)
      allocate(bj(0:lmaxq+1),by(0:lmaxq+1),bess(0:lmaxq+1))
      call sph_bessel_j(lmaxq+1,qr,bj)
      call sph_bessel_y(lmaxq+1,qr,by)
      call phaseshift('n',ie,lmaxq,lx,bj,by,psn)
      lxq(ie,iatom)=lx
      deallocate(bj,by,bess)
    enddo
  end select
enddo !ie
close(10); close(20); close(21)

deallocate(mtpot)

!Jumps of pi removed from the phaseshift versus energy curves.
select case(RelOrNonrel)
 case('r')
  call phaseshift_nojump(psrm)
  call phaseshift_nojump(psrp)
  !the spin averaged phase shift PSR is a weighted average of
  !PSRM and PSRP.
  do ie=1,ne
    do n=1,nocc
      psr(ie,0,n)=psrm(ie,0,n)
      do l=1,lmax
        i=nint((psrm(ie,l,n)-psrp(ie,l,n))/pi)
        psrm(ie,l,n)=psrm(ie,l,n)-i*pi
        psr(ie,l,n)=((l+1)*psrm(ie,l,n)+l*psrp(ie,l,n))/(l+l+1)
      enddo
    enddo
  enddo
 case('n')
  call phaseshift_nojump(psn)
end select

select case(Output)
 case('p')
  call inner_potential(v0r)
  select case(RelOrNonrel)
  case('r')
    call phaseshift_tab(title,psr,'r',v0r)
    select case(SpinPhaseShift)
    case('y')
      call phaseshift_tab(title,psrm,'-',v0r)
      call phaseshift_tab(title,psrp,'+',v0r)
    end select
  case('n')
    call phaseshift_tab(title,psn,'n',v0r)
  end select
  call write_psvse
endselect
select case(Output)
 case('p','s')
  call sigma(lmaxq,lxq,psrm,psrp,psn)
endselect
return
contains
  
  subroutine id_orbit(l)  !orbital quantum number L in filename.
  integer :: l
  idl='000'
  if(l<10)then
    write(idl(3:3),'(i1)') l
  elseif(l<100)then
    write(idl(2:3),'(i2)') l
  else
    write(idl(1:3),'(i3)') l
  endif
  return
  end subroutine id_orbit
  
  subroutine write_psvse
  implicit none
  select case(PSvsE)
  case('y')
    do n=1,nocc
      do l=0,lmaxq
        call id_orbit(l)
        open(11,file='PSvsE.L'//idl//id_at(n))
        select case(RelOrNonrel)
        case('r')
          do ie=1,ne
            emv0=eev(ie)-v0tab(ie)*rydb
            write(11,980) emv0,psr(ie,l,n),psrm(ie,l,n),psrp(ie,l,n)
        enddo
        case('n')
          do ie=1,ne
            emv0=eev(ie)-v0tab(ie)*rydb
            write(11,980) emv0,psn(ie,l,n)
          enddo
        endselect !RelOrNonrel
      enddo       !lmaxq
    enddo         !nocc
  endselect       !PSvsE
980   format(1p4e14.6,0p)
  end subroutine write_psvse
end subroutine scattering
!-----------------------------------------------------------------------
subroutine phaseshift(tag,ie,lmaxq,lx,bj,by,delta)
!PHASESHIFT calculates phase shifts and radial wave functions.
use param,only : RelOrNonrel,Output,UvsR,bohr,c_light,cm2tmp,node,&
  eta_cutoff,ne,id_at,macheps,relerr0,relerr1,relerr2
use dimen,only : nieq,nocc,lmax
use radii,only : rx,rmt
use scatt,only : iatom,l,kappa,emv0,q,a1
use bessl,only : bess
use radeq
implicit none
logical :: w_xec
character(len=*) :: tag
integer,parameter :: l_new_start=5
integer :: ie,l1,lmin,lmaxq,i,i_start,lx,n,ns1,ns2,ns3,ns4
real*8  :: delta(ne,0:lmax,nocc),bj(0:lmax+1),by(0:lmax+1),y(2),&
           y_start(2),t_end,qr,abserr,w_max(2),ps1,ps2,ps3,ps4
real*8,allocatable :: or(:),ou(:,:)
allocate(or(node),ou(2,node))

qr=q*rmt(iatom)
select case('Output')
 case('d')
  write(61,'("tag.atom=",a,", KRMT=",f6.2)') tag//id_at(iatom),qr
  write(61,'("  l    ps        no")')
endselect

lmin=0; if(tag=='+')lmin=1; i_start=1
do l=lmin,lmaxq
  kappa=l; if(tag=='-')kappa=-l-1
  call start_integrate(l,l_new_start,i_start,y_start)
  if(Output=='d'.and.tag/='+'.and.ie==ne.and.iatom==1)then
    write(20,904) q*rx(i_start),l  !LvsKRstart
    !LvsKRstart is virtually independent of the energy, hence ie==ne.
  endif
  
  t_end=rmt(iatom)
  w_xec=.true.; w_max=(/1.d0,1.d0/)
  call wave_integrate_ode&
    (i_start,y_start,t_end,y,relerr0,w_xec,w_max,ps1,ns1)
  w_xec=.false.
  call wave_integrate_ode&
    (i_start,y_start,t_end,y,relerr1,w_xec,w_max,ps1,ns1)
  call write_wave(y,ps1,ns1)
  delta(ie,l,iatom)=ps1
  if(Output=='P'.or.Output=='d')then
    write(61,902) l,ps1,ns1
  endif
  if(l>=6) then
    lx=lmaxq
    if(abs(ps1)<eta_cutoff)then
      lx=l-1
      if(Output=='d'.and.tag/='+'.and.iatom==1)then
        write(21,904) qr,lx !LMAXvsKRmt
      endif
      exit
    endif
  endif
enddo
deallocate(or,ou)
return
900   format(3x,1p4e12.4)
902   format(i5,1pe12.4,0p,i5)
904   format(e15.8,i4)
910   format(i3,4i5,1p4e12.4)
contains

  subroutine wave_integrate_ode&
    (i_start,y_start,t_end,y,relerr,w_xec,w_max,ps,step_count)
  implicit none
  logical :: w_xec
  integer,parameter :: neqn=2,nwork=100+21*neqn
  integer :: i_start,flag,iwork(5),step_count
  real*8  :: t,y(2),y_start(2),tout,t_end,work(nwork),relerr,abserr,&
             ps,a,b,h_start,w_max(2)
  t=rx(i_start); tout=t_end
  y=y_start
  h_start=0.d0
  abserr=10.d0*macheps*max(w_max(1),w_max(2))
  step_count=0; flag=1 
  call ode(waveq,neqn,y,t,tout,relerr,abserr,flag,nwork,work,iwork,&
           h_start,or,ou,step_count,node)
  if(flag>=3)then
    write(61,*)'wave_integrate_ode: l,flag=',l,flag; stop
  endif
  if(w_xec)then
    w_max=0.d0
    do i=1,step_count
      w_max(1)=max(w_max(1),abs(ou(1,i)))
      w_max(2)=max(w_max(2),abs(ou(2,i)))
    enddo
  endif
  a=-(kappa+l+1)*y(1)+rmt(iatom)*(1.d0+cm2tmp*emv0)*y(2)
  b=qr*y(1)
  ps=atan((a*bj(l)+b*bj(l+1))/(a*by(l)+b*by(l+1)))
  return
  end subroutine wave_integrate_ode
  
  subroutine write_wave(y,ps,step_count)
  implicit none
  character :: idl*3
  integer :: step_count
  real*8 :: y(2),ps,c,s,a,g,gp,gnrm,gnrmp,bjl,bjlp,byl,bylp
  select case(UvsR)
  case('y')
    if(ie==ne.and.tag/='+')then
      !normalisation of u1(l,r) to r*j(l,q*r) at muffin-tin radius.
      g=y(1)/rmt(iatom)
      gp=(-(kappa+1)*y(1)/rmt(iatom)+&
         (1.d0+cm2tmp*emv0)*y(2))/rmt(iatom)
      bjl=bj(l); bjlp=(l/qr)*bjl-bj(l+1)
      byl=by(l); bylp=(l/qr)*byl-by(l+1)
      c=cos(ps); s=sin(ps)
      gnrm=c*bjl-s*byl; gnrmp=q*(c*bjlp-s*bylp)
      if(l<l_new_start.and.abs(gnrmp)>=abs(q*gnrm))then
        a=gnrmp/gp
      else
        a=gnrm/g
      endif
      call id_orbit(l,idl)
      open(11,file='UvsR.L'//idl//id_at(iatom))
      do i=1,step_count
        call sph_bessel_j(l,q*or(i),bess)
        write(11,940) or(i)*bohr,ou(1,i)*a,or(i)*bess(l),ou(2,i)*a/c_light
      enddo
    endif
  endselect
  return
940   format(1p4e14.6)
  end subroutine write_wave
  
  subroutine id_orbit(l,idl)  !orbital QN in filename.
  implicit none
  character :: idl*3
  integer :: l
  idl='000'
  if(l<10)then
    write(idl(3:3),'(i1)') l
  elseif(l<100)then
    write(idl(2:3),'(i2)') l
  elseif(l<1000)then
    write(idl(1:3),'(i3)') l
  endif
  return
  end subroutine id_orbit
end subroutine phaseshift
!-----------------------------------------------------------------------
subroutine start_integrate(l,l_new_start,i_start,y_start)
!ODE_INIT initiates the wavefunctions at r=t1, output is y=y_start(1:2).
!the wavefunctions are y(1)=u1 and y(2)=c*u2 in the notation u1,u2
!of Rose; the large Dirac wavefunction component is G=u1/r=y(1)/r;
!in the limit of c=infinity, G is the Schroedinger wavefunction.
use param,only : pi,alfa,c_light,cm2tmp,wave_start,ne
use dimen,only : nieq,lmax
use radii,only : rx,dx,nxt
use scatt,only : iatom,kappa,emv0,q
use bessl,only : bess
use poten,only : z,mtpot
implicit none
integer :: l,l_new_start,i,i_start,j
real*8  :: t1,y_start(2),q1,alfaz,gamma,aa,z0,a0,b0,v0,gg,kg,ee

!potential rV(r) approximated by -2*z0+v0*r .
z0=-0.5d0*(rx(2)*mtpot(1,iatom)-rx(1)*mtpot(2,iatom))/(rx(2)-rx(1))
v0=(mtpot(2,iatom)-mtpot(1,iatom))/(rx(2)-rx(1))
!magnitude of j(l,kr) for small values of kr.
do i=i_start,nxt(iatom)
  call sph_bessel_j(l+1,q*rx(i),bess)
  if(abs(rx(i)*bess(l))>wave_start)then
    j=i
    exit
  endif
enddo
i_start=max(j,1)
if(l<l_new_start)then
  if(cm2tmp>0.d0)then
    !initial condition of Dirac equation.
    !Ref.: M.E.Rose, Relativistic Electron Theory (Wiley, 1961).
    alfaz=alfa*z0
    gamma=sqrt(kappa*kappa-alfaz*alfaz)
    if(kappa>0)then
      a0=alfaz
      b0=kappa+gamma
    else
      a0=kappa-gamma
      b0=alfaz
    endif
    gg=gamma+gamma+1.d0
    kg=kappa/gg
    gg=gamma/gg
    ee=1.d0+(emv0-v0)*cm2tmp
    y_start(1)=rx(1)*(a0+((1.d0-kg)*ee-gg)*c_light*b0*rx(1))
    y_start(2)=rx(1)*(b0+((1.d0+kg)*ee-gg)*c_light*a0*rx(1))*c_light
  else
    !initial condition of Schroedinger equation.
    !Ref.: J.B.Pendry, Low energy electron diffraction (Academic,1974).
    ee=z0*z0/(l+1)-0.5d0*(emv0-v0)
    t1=rx(1)
    y_start(1)=(1.d0 +(-z0/(l+1)+(ee/(l+l+3))*t1)*t1)*t1
    y_start(2)=l+l+1+(-2.d0*z0+ee*t1)*t1
  endif
else
  !initial condition when l>=lbessel:
  !u1(l,r) is similar to r*j(l,kr) and approximately zero over a radial
  !interval [0,t1]. The integration of u1 is stable for r>r1, because
  !the linearly independent companion wavefunction to u1 is similar to
  !r*y(l,kr) and quickly vanishing. Ref.: this work.
  t1=rx(i_start)
  y_start(1)=t1*bess(l)
  y_start(2)=( (kappa+l+1)*bess(l) - q*t1*bess(l+1) )/&
             (1.d0+cm2tmp*(emv0-mtpot(i_start,iatom)/t1))
                               !mtpot is tot.pot.*radius.
endif
return
end subroutine start_integrate
!-----------------------------------------------------------------------
subroutine phaseshift_nojump(delta)
!PHASESHIFT_NOJUMP removes jumps of pi.
!Ref.: W. Moritz,
!Ref.: A. Barbieri and M.A. Van Hove, Phase shift package, 
!      http://electron.lbl.gov/software.
!Ref.: DLPHASE v1.1, CCP3 Library code,
!      http://www.cse.clrc.ac.uk/Activity.
!Ref.: K. Heinz and V. Blum, Erlangen phase shift packaage,
!      private communication.
use param,only : pi,ne
use dimen,only : nocc,lmax
implicit none
integer :: n,l,npi,i
real*8  :: delta(ne,0:lmax,nocc),dif
real*8,allocatable :: del(:)

allocate(del(ne))

do n=1,nocc
  do l=0,lmax
    del(1:ne)=delta(1:ne,l,n)
    npi=0
    do i=2,ne
      dif=del(i)-del(i-1)
      if(dif>1.9d0)then
        npi=npi-1
      elseif(dif<-1.9d0)then
        npi=npi+1
      endif
      delta(i,l,n)=del(i)+pi*npi
    enddo
  enddo !l
enddo !n
return
end subroutine phaseshift_nojump
!-----------------------------------------------------------------------
subroutine inner_potential(v0r)
!INNER_POTENTIAL approximates the inner potential by the expression
!	V0mVfe=a+b/sqrt(e+c).
!Optimum parameters a,b,c are found by the downhill simplex method of
!Nelder and Mead.
use param,only : RelOrNonrel,rydb,ne,eev
use poten,only : v0tab
use error
implicit none
integer :: j,i,n1,n2,iter,ie,ie20,lo,hi,count
integer,parameter :: n=4
real*8 :: v0r(4),v(0:4,4),y(0:4),tmp(4),eps,fit,norm
real*8,allocatable :: v0ev(:)
allocate(v0ev(ne))

v0ev=v0tab*rydb
!constant inner potential.
open(60,file='V0vsE',status='replace')
if(abs(v0ev(1)-v0ev(ne))<0.05d0)then
  write(61,'(a,f6.2)')'constant inner potential, V0mVfe=',v0ev(1)
  v0r(1)=v0ev(1); v0r(2)=v0ev(1); v0r(3)=0.d0; v0r(4)=0.d0
  write(60,'(3e12.4)') eev(1),v0ev(1),v0ev(1)
  write(60,'(3e12.4)') eev(ne),v0ev(ne),v0ev(ne)
  close(60)
  return
endif

!V0(E) by the downhill simplex method of Nelder and Mead;
!the design of an initial simplex is discussed by W.H.Press et al.,
!Numerical Recipes.
do i=0,n
  v(i,1)=-12.d0
  v(i,2)=  0.d0
  v(i,3)=-60.d0
  v(i,4)=  8.d0
enddo
!initial simplex iter=1, eps=1.e-04 and restart iter=2, eps=1.e-06.
eps=1.d-04
do iter=1,2
  v(1,1)=v(1,1)+2.d0
  v(2,2)=v(2,2)-0.5d0
  v(3,3)=v(3,3)-5.d0
  v(4,4)=v(4,4)+2.d0
  do i=0,n
    tmp(1:n)=v(i,1:n)
    y(i)=v0_err(tmp)
  enddo
  call simplex(v0_err,v,y,n,eps,lo,hi,norm,count)
  write(61,'(a,i3,1pe7.0)')'simplex count,acc=',count,norm
  eps=eps*1.d-02
enddo
!do j=1,n
!  v0r(j)=sum(v(1:n,j))/n
!enddo
v0r(1:n)=v(lo,1:n)

!print.
eps=0.d0
do ie=1,ne
  if(eev(ie)+v0r(4)>=0.d0)then
    fit=v0r(2)+v0r(3)/sqrt(eev(ie)+v0r(4))
    if(fit>=v0r(1))then
      eps=max(eps,abs(v0ev(ie)-fit))
    endif
    write(60,'(1p3e12.4)') eev(ie),v0ev(ie),max(v0r(1),fit)
  endif
enddo
write(61,900) v0r
write(61,902) ne,v0ev(1)-v0ev(ne)
write(61,904) eps
900   format('V0mVfe=max(c0,c1+c2/sqrt(e+c3)), where ci=',4f8.2)
902   format('V0mVfe(1)-V0mVfe(',i2,')=',f6.2,' eV')
904   format('rms misfit= ',f4.2,' eV')
return
end subroutine inner_potential
!-----------------------------------------------------------------------
subroutine phaseshift_tab(title,delta,tag,v0r)
!PHASESHIFT_TAB transforms the set of phase shifts delta(l,E-V0) onto
!a uniform E-V0 grid using spline interpolation.
use param,only : rydb,ne,eev,id_at,ep1,ep2,dep,pi
use dimen,only : nieq,nocc,lmax
use poten,only : v0tab
implicit none
character(len=*) :: title,tag
integer :: n,i,i1,i2,l,lmin,idep,ie,j,last_interval,k
real*8  :: delta(ne,0:lmax,nocc),v0r(4),a1,prim,e1,e2,pih
real*8,allocatable,dimension(:)   :: emv0,db,dc,dd,ep
real*8,allocatable,dimension(:,:) :: ps
pih=pi*0.5d0
     
allocate(emv0(ne),db(ne),dc(ne),dd(ne))

emv0=eev-v0tab*rydb
e1=emv0(1)
idep=dep
i1=(int(e1)/idep)*idep; a1=i1; if(a1<e1)i1=i1+idep
i2=(int(emv0(ne)+0.001d0)/idep)*idep
j=(i2-i1/idep)+1
  
allocate(ep(1:j),ps(1:j,0:lmax))
!!!!!! A. Imre 17.09.21: lmin was uninitialized? Should be 0 I assume.
lmin = 0

do n=1,nocc  
  !spline interpolation.
  ps=0.d0
  do l=lmin,lmax
    call spcoef(ne,emv0,delta(1,l,n),1.d+33,1.d+33,db,dc,dd,&
               'phaseshift_tab')
    last_interval=1
    j=0
    do i=i1,i2,idep
      j=j+1
      ep(j)=i
      if(ep(j)>ep2)then
        j=j-1
        exit
      endif
      call svalue(ne,emv0,delta(1,l,n),db,dc,dd,ep(j),ps(j,l),prim,&
                  last_interval,'phaseshift_tab')
    enddo
  enddo
  !phase shift table.
  open(60,file='PS.'//tag//id_at(n),status='replace')
  write(60,900) lmax+1,v0r,'PS.'//tag//id_at(n),title
  do l=0,lmax
    if(ps(1,l)<=-pih) ps(:,l)=ps(:,l)+pi
    if(ps(1,l)>pih) ps(:,l)=ps(:,l)-pi
  enddo
  do i=1,j
    write(60,910) ep(i),(ps(i,l),l=0,lmax)
  enddo
enddo !n
return
900   format(i2,f7.2,3f8.2,1x,a,1x,a)
910   format(f8.1,8f8.4/(10f8.4))
end subroutine phaseshift_tab
!-----------------------------------------------------------------------
subroutine sigma(lmaxq,lxq,psrm,psrp,psn)
!SIGMA calculates differential cross section and total cross section in
!atomic units (Bohr radius squared).
use param,only : RelOrNonrel,Omega,Theta,Sherman,CrossSection,&
  pi,rydb,bohr,cm2,eta_cutoff,ne,eev,id_at
use dimen,only : nieq,nocc,lmax
use poten,only : v0tab
implicit none
character :: io*5,no*2
integer,parameter :: ntheta=360
integer :: ie,i,l,n,lx,k,itheta,lmaxq,lxq(ne,nocc)
real*8  :: psrm(ne,0:lmax,nocc),psrp(ne,0:lmax,nocc),&
           psn(ne,0:lmax,nocc),dtheta,fg,fgsum,d,emv0,tpi,frac
real*8,allocatable :: facr(:),facn(:),p(:,:,:),dsdo(:,:),shf(:,:)
complex*16 :: f,g,s
complex*16,allocatable :: ckn(:),ckp(:),cnr(:)
complex*16,parameter :: ci=(0.d0,1.d0)
allocate(facr(ne),facn(ne))
allocate(ckn(0:lmaxq),ckp(0:lmaxq),cnr(0:lmaxq))
allocate(dsdo(0:ntheta,ne),p(0:ntheta,0:lmax,0:1),shf(0:ntheta,ne))

!prefactor |1/ik|**2 accompanying exp(i*ps)*sin(ps).
do ie=1,ne
  emv0=eev(ie)/rydb-v0tab(ie)
  facn(ie)=1.d0/emv0
  facr(ie)=1.d0/(emv0*(1.d0+emv0*cm2))
enddo
dtheta=pi/ntheta; tpi=pi+pi
call legendre_polynomial
do n=1,nocc
  select case(RelOrNonrel)
  case('r')
    !relativistic calculation.
    do ie=1,ne
      do l=0,lxq(ie,n)
        ckn(l)=exp(ci*psrm(ie,l,n))*sin(psrm(ie,l,n))
        ckp(l)=exp(ci*psrp(ie,l,n))*sin(psrp(ie,l,n))
      enddo
      do itheta=0,ntheta
        f=ckn(0)
        g=(0.d0,0.d0)
        do l=1,lxq(ie,n)
          f=f+((l+1)*ckn(l)+l*ckp(l))*p(itheta,l,0)
          g=g+(ckn(l)-ckp(l))*p(itheta,l,1)
        enddo
        fg=abs(f)**2+abs(g)**2
        dsdo(itheta,ie)=fg*facr(ie)
        shf(itheta,ie)=ci*(f*conjg(g)-conjg(f)*g)/fg
      enddo !itheta
    enddo !ie
  case('n')
    !nonrelativistic calculation.
    do ie=1,ne
      do l=0,lxq(ie,n)
        cnr(l)=exp(ci*psn(ie,l,n))*sin(psn(ie,l,n))
      enddo
      do itheta=0,ntheta
        s=cnr(0)
        do l=1,lxq(ie,n)
          s=s+(l+l+1)*cnr(l)*p(itheta,l,0)
        enddo
        dsdo(itheta,ie)=abs(s)**2*facn(ie)
      enddo !itheta
    enddo !ie
  end select
  
  select case(CrossSection)
  case('y')
    open(62,file='sigma.'//RelOrNonrel//id_at(n),status='replace')
  end select
  frac=180.d0/ntheta
  do ie=1,ne
    call id_energy
    
    select case(Omega)
    case('y')
      open(63,file='dsdo'//io//'.'//RelOrNonrel//id_at(n),status='replace')
      do itheta=0,ntheta
        write(63,900) itheta*frac,log10(dsdo(itheta,ie))
      enddo
    end select
    
    select case(Theta)
    case('y')
      open(64,file='dsdt'//io//'.'//RelOrNonrel//id_at(n),status='replace')
      do itheta=0,ntheta
        write(64,900) itheta*frac,dsdo(itheta,ie)*sin(itheta*dtheta)*tpi
      enddo
    end select
    
    select case(RelOrNonrel)
    case('r')
      select case(Sherman)
      case('y')
        open(65,file='shf'//io//'.'//RelOrNonrel//id_at(n),status='replace')
        do itheta=0,ntheta
          write(65,900) itheta*frac,shf(itheta,ie)
        enddo
      end select
    end select
    
    select case(CrossSection)
    case('y')
      write(62,910) eev(ie),integral_ugt(dsdo(0,ie),ntheta+1,dtheta)
    end select
    
  enddo !ie
enddo !n
  900 format(f6.2,1p20e14.6)
  910 format(f6.0,1pe14.6)
return
contains

  subroutine legendre_polynomial
  !LEGENDRE_POLYNOMIAL with calculation in quartic precision and
  !output in double precision.
  
  !dtheta      =pi/ntheta,
  !p(theta,l,m)=legendre function of degree l, and order m.
  !qp(l,0)     =legendre polynomial of degree l for a given angle,
  !qp(l,1)     =associated legendre polynomial of degree l and order 1.

  !The recurrence relations of the legendre polynomials are used.
  !The stability the recurrence was tested by a calculation from 
  !l=0 to 200 and back to 0. The maximum disagreement for l=0:200 and
  !theta=1:179 degrees was found to be 4.E-30. One has reason to
  !believe that all digits are correct in the double precision output. 
  !Ref.: Abramowitz and Stegun, page XIII and Sec. 8.5.3.
  integer :: i,l,ll
  real*8 :: x,c,s,dtheta
  real*8,allocatable :: qp(:,:)
  dtheta=acos(-1.d0)/ntheta
  !theta equal to 0 and pi.
  p(0,0:lmax,0)=1.d0
  do l=0,lmax,2
    p(ntheta,l,0)=1.d0
  enddo
  do l=1,lmax,2
    p(ntheta,l,0)=-1.d0
  enddo
  p(0,0:lmax,1)=0.d0; p(ntheta,0:lmax,1)=0.d0
  !theta between 0 and pi.
  allocate(qp(0:lmax,0:1))
  do i=1,ntheta-1
    x=dtheta*i
    c=cos(x);    s=sin(x)
    qp(0,0)=1.d0; qp(0,1)=0.d0
    qp(1,0)=c;    qp(1,1)=s
    ll=1
    do l=2,lmax
      ll=ll+2
      qp(l,0)=(ll*c*qp(l-1,0)-(l-1)*qp(l-2,0))/l
      qp(l,1)=(ll*c*qp(l-1,1)- l   *qp(l-2,1))/(l-1)
    enddo
    p(i,0:lmax,0)=qp(0:lmax,0)
    p(i,0:lmax,1)=qp(0:lmax,1)
  enddo
  return
  end subroutine legendre_polynomial
  
  real*8 function integral_ugt(f,n,d)
  !INTEGRAL_UGT is uniform-grid trapezoidal integration.
  !Ref.: Abramowitz and Stegun, chap. 25.2.
  integer :: n,i
  real*8  :: f(n),d,sum,w1,w2
  sum=0.d0
  w1=f(1)
  do i=2,n
    w2=f(i)
    sum=sum+w1+w2
    w1=w2
  enddo
  integral_ugt=0.5d0*d*sum
  return
  end function integral_ugt
  
  subroutine id_energy  !energy value in filename.
  io='00000'
  i=eev(ie)
  if(i<10)then
    write(io(5:5),'(i1)') i
    elseif(i<100)then
      write(io(4:5),'(i2)') i
    elseif(i<1000)then
      write(io(3:5),'(i3)') i
    elseif(i<10000)then
      write(io(2:5),'(i4)') i
    elseif(i<100000)then
      write(io(1:5),'(i5)') i
    endif
    return
    end subroutine id_energy
end subroutine sigma
!-----------------------------------------------------------------------
real*8 function integral_egt(m,n,f,x1,x2,f1)
!INTEGRAL_EGT is exponential-grid trapezoidal integration;
!f(r)*(r**m) is integrated.
!Ref.: Abramowitz and Stegun, chap. 25.2.
use radii,only : rx,dx
implicit none
integer,intent(in) :: m
integer :: n,i1,i2,i
real*8  :: f(n),x1,x2,f1,f2,p1,p2,sum,w1,w2

!p is a subscript value on a continuous exponential-grid.
p1=log(x1/rx(1))/dx+1.d0
p2=log(x2/rx(1))/dx+1.d0
i1=max(int(p1),2); i1=min(int(p1),n-1)
i2=max(int(p2),2); i2=min(int(p2),n-1)
!trapezoidal integration on the interval [rx(i1),rx(i2)].
sum=0.d0
w1=f(i1)*rx(i1)**m
do i=i1+1,i2
  w2=f(i)*rx(i)**m
  sum=sum+w1+w2
  w1=w2
enddo
!end-point corrections:
!lagrangian 3-point interpolation at x1 and x2, followed by trapezoidal
!integration on [rx(i1),x1] (subtraction) and on [rx(i2),x2] (addition).
if(i1==1) i1=2
p1=p1-i1
p2=p2-i2
f1=0.5d0*p1*(p1-1.d0)*f(i1-1)+(1.d0-p1*p1)*f(i1)+0.5d0*p1*(p1+1.d0)*&
   f(i1+1)
f2=0.5d0*p2*(p2-1.d0)*f(i2-1)+(1.d0-p2*p2)*f(i2)+0.5d0*p2*(p2+1.d0)*&
   f(i2+1)
sum=sum-(f(i1)*rx(i1)**m+f1*x1**m)*p1+(f(i2)*rx(i2)**m+f2*x2**m)*p2
integral_egt=0.5d0*dx*sum
return
end function integral_egt
!-----------------------------------------------------------------------
real*8 function interpol_egl3(x1,f1)
!INTERPOL_EGL3 is 'exponential-grid lagrangian 3-point interpolation.
!Ref.: Abramowitz and Stegun, Sec. 25.2.11.
use dimen,only : nx
use radii,only : rx,dx
implicit none
integer :: i1
real*8  :: f1(nx),x1,p1,p,f1x
p1=log(x1/rx(1))/dx+1.d0
i1=p1
p=p1-i1
f1x=0.5d0*p*(p-1.d0)*f1(i1-1)+(1.d0-p*p)*f1(i1)+0.5d0*p*(p+1.d0)*&
    f1(i1+1)
interpol_egl3=f1x
return
end function interpol_egl3
!-----------------------------------------------------------------------
subroutine sph_bessel_j(lmax,x,bj)
!SPH_BESSEL_J is spherical bessel function j(l,x).
implicit none
integer :: l,ll,lx,lmax,llp1,llp3
real*8  :: bj(0:lmax),x,xi,xl,x2h,x2hh,cl,f,s
real*8,allocatable :: aux(:)
!small arguments: limit determined by mashine precision 2.22E-16 .
if(x<1.0188d-02)then
  x2h=x*x*0.5d0; x2hh=x2h*0.5d0
  bj(0)=1.d0-(x2h/3.d0)*(1.d0-x2hh/5.d0)
  xl=1.d0; cl=1.d0
  do l=1,lmax
    xl=xl*x; ll=l+l; cl=cl/(ll+1)
    bj(l)=cl*xl*(1.d0-(x2h/(ll+3))*(1.d0-x2hh/(ll+5)))
  enddo
  return
endif
!large arguments.
xi=1.d0/x
if(x>dfloat(lmax))then
  !upward recurrence.
  bj(0)=sin(x)*xi
  if(lmax==0) return
  bj(1)=(bj(0)-cos(x))*xi
  if(lmax==1) return
  ll=3
  do l=2,lmax
    bj(l)=(ll*xi)*bj(l-1)-bj(l-2)
    ll=ll+2
  enddo
else
  !downward recurrence.
  if(lmax<100)then
    lx=lmax+10+int(x*(1.d0-dfloat(lmax)*0.0075d0))
  else
    lx=lmax+10+int(x*0.25d0)
  endif
  allocate(aux(0:lx))
  aux(lx)=0.d0
  aux(lx-1)=1.d-30
  llp1=lx+lx-1; llp3=lx+lx+1
  s=llp1*aux(lx-1)**2
  do l=lx-2,0,-1
    llp1=llp1-2; llp3=llp3-2
    aux(l)=(llp3*xi)*aux(l+1)-aux(l+2)
    s=s+llp1*aux(l)*aux(l)
  enddo
  !normalisation by sum (2*l+1)*j(l,x)**2=1, l=0:lx.
  f=1.d0/sqrt(s)
  bj(0:lmax)=aux(0:lmax)*f
endif
return
end subroutine sph_bessel_j
!-----------------------------------------------------------------------
subroutine sph_bessel_y(lmax,x,by)
!SPH_BESSEL_Y is spherical bessel function y(l,x).
implicit none
integer :: l,ll,lmax
real*8  :: by(0:lmax),x,xi
xi=1.d0/x
by(0)=-cos(x)*xi
if(lmax==0) return
by(1)=(by(0)-sin(x))*xi
if(lmax==1) return
ll=3
do l=2,lmax
  by(l)=(ll*xi)*by(l-1)-by(l-2)
  ll=ll+2
enddo
return
end subroutine sph_bessel_y
