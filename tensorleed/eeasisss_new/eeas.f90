  !=======================================================================
  !This program is free software under the terms of the GNU General Public
  !License as published by the Free Software Foundation.
  !Author: John O. Rundgren, jru@KTH.se ,
  !        KTH Royal Institute of Technology, Stockholm, Sweden.
  !Version: 28 March 2021.
  !-----------------------------------------------------------------------
  program eeas

  !program EEAS
  !Elastic Electron-Atom Scattering in Solids and Solid Surfaces,
  !author: John Rundgren, jru@kth.se .
  !The author appreciates acknowledgement in publications by citation:
  !J. Rundgren, Phys.Rev.B68,125405(2003),
  !             Phys.Rev.B76,195441(2007).
  implicit none
  character(len=1) :: Pot,WF
  character(len=127) :: outdir
  character(len=3) :: cie
  character,allocatable :: idElemA(:)*5,idZA(:)*5

  integer,parameter :: dp=selected_real_kind(15,307)
  integer :: nxx,nieq,ne,nlat,nshell,nsp,nsr,nthread,lmax
  integer :: ie,ithread
  integer,allocatable :: neq(:),nx(:),nxR(:),ia(:,:)

  real(dp),parameter :: rydb=13.60569172d0,pi=acos(-1.d0)
  real(dp) :: Einc,Vcry0,Vxc0,gam,zet
  real(dp),allocatable :: ad(:,:),eev(:),v0ev(:),sdat(:,:),sp(:),sr(:),&
    rmt(:),rx(:,:),dx(:),z(:),rho(:,:),rs(:,:),fxc(:),&
    Vcry(:,:),VxcR(:),Vtot(:,:),Vxc(:,:)

  !PScal variables.
  integer :: iatom,kappa,noru,nory

  !alfa=fine-structure constant, c=speed of light, cm2=1/c**2.
  real(dp),parameter :: alfa=7.297352533d-03
  real(dp),parameter :: clight=2.d0/alfa, cm2=alfa*alfa*0.25d0

  !phase shifts.
  integer :: l,ir
  real(dp) :: q,emv0
  real(dp),allocatable :: psu(:,:),psd(:,:),psu1(:,:),psd1(:,:)
  real(dp) :: relerr,abserr

  !LINDHARD's const = Lconst.
  real(dp) :: pF,Esdat,Lsdat
  real(dp),allocatable :: EmVxc0(:),Lrs(:),Ltest(:),Lconst(:)

  !EEAS import from EEASiSSS.
  open(1732,file='../uinp1',form='unformatted',status='unknown')
  read(1732) nxx,nieq,ne,nlat,nshell,nsp,nsr,nthread,lmax
  close(1732)
  !
  allocate(neq(nieq),nx(nieq),ia(nshell,nieq),ad(nshell,nieq), &
    idElemA(nieq),idZA(nieq),eev(ne),v0ev(ne),sdat(nsp,nsr),sp(nsp),sr(nsr), &
    rx(nxx,nieq),dx(nieq),rmt(nieq),z(nieq),rho(nxx,nieq),rs(nxx,nieq), &
    Vtot(nxx,nieq),Vcry(nxx,nieq),Vxc(nxx,nieq),VxcR(nieq),nxR(nieq), &
    fxc(nieq),EmVxc0(nieq),Lconst(nieq),Lrs(nieq),Ltest(nieq), &
    psu(0:lmax,nieq),psd(0:lmax,nieq), &
    psu1(0:lmax,nieq),psd1(0:lmax,nieq) )
  !
  open(1732,file='../uinp2',form='unformatted',status='unknown')
  read(1732) outdir,Pot,WF,idElemA,idZA,neq,nx,dx,rx,ad,ia, &
    sp,sr,sdat,eev,Vcry,Vcry0,rmt,nxR,z,rho,rs,fxc,relerr,abserr
  close(1732)
  !
  !THREAD definition.
  open(20,file='threadno',status='unknown')
  read(20,*) ithread
  close(20)
  !
  !LOGFILE and RESULT of eeas(ithread,ie).
  psu=0.d0; psd=0.d0; psu1=0.d0; psd1=0.d0
  do ie=ithread,ne,nthread
    write(cie,'(i0)') ie
    open(611,file='../ulog'//trim(cie),status='unknown')
    Einc=eev(ie)/rydb
    call MTvsE(ie)
    call PSvsE(ie)
    open(20,file='../udat'//trim(cie),status='unknown',access='stream') 
    write(20) Vxc0, &
              ((psu(l,ir),l=0,lmax),ir=1,nieq), &
              ((psd(l,ir),l=0,lmax),ir=1,nieq), &
              ((psu1(l,ir),l=0,lmax),ir=1,nieq), &
              ((psd1(l,ir),l=0,lmax),ir=1,nieq)
    if(ie==ne)then
    write(611,'(/20("===="))')
    write(611,'(a)')'Lindhard correlation potential -Lconst/p at Einc = max'
    write(611,'(a)')'      for the atoms of the given crystal'
    write(611,'(20("----"))')
    write(611,'(a,20(2x,a):)')'atom   =',idElemA(:)
    write(611,'(a)')
    write(611,'(a)')'Upper energy boundary of Sernelius data base sdat(p/pF,rs),'
    write(611,'(a)')      'Einc-Vxc0(eV)'
    write(611,'(a,20(1x,f6.0):)')'      =',EmVxc0(:)*rydb
    write(611,'(a)')
    write(611,'(a)')'Electron-space radius rs at MT radius (rs in a0),'
    write(611,'(a,20(1x,f6.2):)')'MT rs =',Lrs(:)
    write(611,'(a)')
    write(611,'(a)')"Lindhard const extrapolated from sdat (MT pot. & Vxc0 step-free),"
    write(611,'(a,20(1x,f6.2):)')'const =',Ltest(:)
    write(611,'(a)')
    write(611,'(a)')'Lindhard const per atom: -sqrt((pi^3)*rho),'
    write(611,'(a,20(1x,f6.2):)')'const =',Lconst(:)
  endif
  enddo
  stop
contains

!=======================================================================
! MUFFIN-TIN SPHERES
!-----------------------------------------------------------------------
  subroutine MTvsE(ie)
  implicit none
  integer :: ie,ir,i
  real(dp) :: Einc
  !
  Einc = max(eev(ie)/rydb,1.d-04)
  write(611,*)
  write(611,900)' ie Einc(eV)',ie,eev(ie)
  !write(*  ,900)' ie Einc(eV)',ie,eev(ie)
  !
  !Vxc(i,ir) generated.        
  call calc_Vxc(Einc)
  !
  do ir=1,nieq
    VxcR(ir) = ipol(rmt(ir),nxR(ir),rx(1,ir),dx(ir),Vxc(1,ir),'Vxc')
  enddo
  Vxc0 = sum(neq*VxcR)/dble(nlat)
  !Normalization to MT-interstice continuity.
  do ir=1,nieq
    Vxc(1:nxR(ir),ir) = Vxc(1:nxR(ir),ir) - (VxcR(ir) - Vxc0)
  enddo
  ! 
  write(611,822)'VxcR(eV)',(VxcR(ir)*rydb,ir=1,nieq)
  write(611,823)'Vxc0(eV) = ',Vxc0*rydb
  if(Pot=='y'.and.ie==ne)then
    do ir=1,nieq
      open(10, &
      file='../'//trim(outdir)//'/'//'Vxc.'//idZA(ir),status='unknown')
      do i=1,nxR(ir)
        if(Vxc(i,ir)*rydb > -100.d0)then
          write(10,910) rx(i,ir),Vxc(i,ir)*rydb 
        endif
      enddo
    enddo
  endif
  !
  !Phase shift generated from XC modulated carrier potential.
  do ir=1,nieq
    Vtot(1:nxR(ir),ir) = &
    Vcry(1:nxR(ir),ir) + Vxc(1:nxR(ir),ir) * fxc(ir)
  enddo
  emv0=Einc-Vxc0
  q=sqrt(emv0 + cm2*emv0**2)
  !
  write(611,'(a)')' ie Einc(eV)     emv0     Vxc0'
  write(611,901) ie,eev(ie),emv0*rydb,Vxc0*rydb
  900 format(a/i3,f9.4)
  901 format(i3,3f9.4)
  822 format(a/(10f8.3:))
  823 format(a,f0.3,2(2x,f0.3))
  910 format(2es14.6)
  return
  end subroutine MTvsE
  !-------------------------------------------------------------------
  subroutine calc_Vxc(Einc)
  !Sernelius's XC potential;
  implicit none
  integer :: ir,i
  real(dp) :: Einc,f(nxx),rho(nieq)
  !
  do ir=1,nieq
    do i=1,nxR(ir)
      Vxc(i,ir) = sigma(Einc,rs(i,ir))
    enddo
    !
    !smooth Vxc(i,ir) by successive 3-point binomial convolution.
    do i=1,nxR(ir)
      f(i)=Vxc(i,ir)
    enddo
    call smooth(f,nxR(ir))
    do i=1,nxR(ir)
      Vxc(i,ir)=f(i)
    enddo
    !
    !Lindhard high-energy correlation pot. =-Lconst/p.
    if(ie==ne)then
      EmVxc0(ir) = Esdat-Vxc0
      Ltest(ir)  = Lsdat
      Lrs(ir)    = rs(nxR(ir),ir)
      rho(ir)    = (0.75d0/pi)/rs(nxR(ir),ir)**3
      Lconst(ir)  = -sqrt(pi**3 *rho(ir))
    endif
  enddo
  return
  end subroutine calc_Vxc
  !-------------------------------------------------------------------
  function sigma(Einc,rs)
  implicit none
  integer :: it
  real(dp) :: sigma,Einc,rs,EF,p_pF1,p_pF2,d,s0,s1,s21,s22,a1,a2
  real(dp),parameter :: pi=acos(-1.d0),pi9_4th=(pi*2.25d0)**(1.d0/3.d0)
  !
  !The energy balance is written
  !p/pF  = sqrt{ 1 + Einc/EF + [sdat(1,rs) - sdat(p/pF,rs)] * rs }
  !sdatw = work space corresponding to sdat.
  !
  pF=pi9_4th/rs
  EF=pF*pF
  s0=1.d0+Einc/EF
  s1=sdatw(1.d0,rs)*rs   !sigma in right-hand part of energy balance.
  !
  p_pF1=sqrt(s0)
  s21=sdatw(p_pF1,rs)*rs !sigma in left-hand part of energy balance.
  a1=s21/rs
  do it=1,4
    p_pF2=sqrt(s0 + s1 - s21)
    s22=sdatw(p_pF2,rs)*rs
      !
      !it termination; a1,a2 = sdat values -0.7xyz to -0.0xyz .
      a2=s22/rs
      d=a2-a1
      if(abs(d) < 0.0005d0)then
        !write(611,*)'it ',it
        exit
      endif
      a1=a2
      !
    p_pF1=p_pF2
    s21=s22
  enddo
  sigma = s22*EF
  !sdat termination.
  Esdat = (3.d0*pF)**2
  return
  end function sigma
  !-------------------------------------------------------------------
  function sdatw(ps,rs)  !ps,rs from sigma
  implicit none
  real(dp) :: sdatw,ps,rs,ps3,v3,g3,a,b
  !sdatw extending values inside sdat(sp,sr) to energies higher
  !than sp(nps).
  !
  if(ps<=sp(nsp))then
    sdatw=sdatb(ps,rs)
    return
  else
    !sdatw -> Hartree-Fock exchange potential \prop -1/p**2 and
    !         Lindhard correlation potential \prop -1/p,
    !         the latter remains at TEM energy .
    !model:
    ! sdatw = a/(ps+b), a<0, b>0,
    ! gradient g = -a/(ps+b)**2, g>0,
    !
    ! a = v3*(ps3+b) = -v3*v3/g3, check a<0,
    ! b = -ps3 - v3/g3, check ps-ps3 >= 0 and -v3/g3 > 0.
    !
    ps3 = 3.00d0
    v3 = sdatb(ps3,rs)
    g3 = 10.d0* &    
    (0.5d0*sdatb(2.80d0,rs)-2.d0*sdatb(2.90d0,rs)+1.5d0*sdatb(3.00d0,rs))
    a = -v3*v3/g3
    b = -ps3-v3/g3
    sdatw = a/(ps+b)
    !LINDHARD approximation.
    Lsdat = a*pF
  endif
  return
  end function sdatw
  !---------------------------------------------------------------------
  function sdatb(ps,rs)
  implicit none
  integer :: k,k1,h,h1
  real(dp) :: sdatb,ps,rs,q,p
  !bivariate interpolation from
  !NIST Handbook of Mathematical Functions, 25.2.66, p. 882.
  !
  !k close to ps and h close to rs. 
  k = minloc(abs(ps-sp(:)),1)
  h = minloc(abs(rs-sr(:)),1)
  !
  if(ps-sp(k) >= 0.d0)then
    if(k==nsp)then
      k1 = k; q = 1.d0
    else
      k1 = k+1; q = (ps-sp(k))/(sp(k1)-sp(k))
    endif
  else
    if(k==1)then
      k1 = k; q = 1.d0
    else
      k1 = k-1; q = (ps-sp(k))/(sp(k)-sp(k1))
    endif
  endif
  !
  if(rs-sr(h) >= 0.d0)then
    if(h==nsr)then
      h1 = h; p = 1.d0
    else
      h1 = h+1; p = (rs-sr(h))/(sr(h1)-sr(h))
    endif
  else
    if(h==1)then
      h1 = h; p = 1.d0
    else
      h1 = h-1; p = (rs-sr(h))/(sr(h)-sr(h1))
    endif
  endif
  sdatb = (1.d0-p)*(1.d0-q)*sdat(k ,h ) + &
                 p*(1.d0-q)*sdat(k ,h1) + &
                 q*(1.d0-p)*sdat(k1,h ) + &
                        p*q*sdat(k1,h1)
  return
  end function sdatb
  !---------------------------------------------------------------------
  subroutine smooth(f,nx)
  implicit none
  integer :: i,j,nx
  real(dp) :: f(nx),g(nx)
  real(dp),parameter :: thrd=1.d0/3.d0
  !Curve smoothing by binomial convolution of coefficients 1,2,1.
  !When executed many times, it converges to the normal distribution.
  !
  do j=1,64
    !f goes to g.
    g(1)=(f(2)+f(1)+f(1))*thrd
    do i=2,nx-1
      g(i)=(f(i-1)+f(i+1)+f(i)+f(i))*0.25d0
    enddo
    g(nx)=(f(nx-1)+f(nx)+f(nx))*thrd
    !g goes back to f.
    f(1)=(g(2)+g(1)+g(1))*thrd
    do i=2,nx-1
      f(i)=(g(i-1)+g(i+1)+g(i)+g(i))*0.25d0
    enddo
    f(nx)=(g(nx-1)+g(nx)+g(nx))*thrd
  enddo
  return
  end subroutine smooth
  !-------------------------------------------------------------------
  function ipol(x,nx,r,d,f,messg)
  implicit none
  integer  :: nx,i1  
  real(dp) :: ipol,x,r(nx),d,f(nx),p1,p
  character(len=*) :: messg
  !ipol is exponential-grid lagrangian 3-point interpolation.
  !Ref.: Abramowitz and Stegun, Sec. 25.2.11.
  !
  if(x-r(nx)>0.d0)then
    write(611,900)'ipol stop: x-r(nx),r(nx)=',x-r(nx),r(nx)
    write(611,902)'error location: ',messg
    stop
  endif
  p1=(log(x/r(1)))/d+1.d0
  i1=nint(p1)
  if(i1==0) i1=1; if(i1==nx) i1=i1-1
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
  !
  !=====================================================================
  ! PHASE SHIFTS
  !---------------------------------------------------------------------
  subroutine PSvsE(ie)
  implicit none
  integer ::  ie,ir,l
  real(dp) :: bj(0:lmax+1),by(0:lmax+1),bjp(0:lmax+1),byp(0:lmax+1),qr
  !
  noru=0; nory=0
  do ir=1,nieq
    qr=q*rmt(ir)
    call sph_bessel_j(lmax,qr,bj)
    call sph_bessel_y(lmax,qr,by)
    do l=0,lmax
      bjp(l)=dble(l)/qr*bj(l)-bj(l+1)
      byp(l)=dble(l)/qr*by(l)-by(l+1)
    enddo
    !
    Vtot(1:nxR(ir),ir) = Vtot(1:nxR(ir),ir) * rx(1:nxR(ir),ir)
    call PScalc(ie,ir,'su',bj,by,bjp,byp)
    call PScalc(ie,ir,'sd',bj,by,bjp,byp)
  enddo
  write(611,'(3(a,i0))') &
  'ODE STEPu = ',noru,', STEPy = ',nory,', RADIAL grid = ',maxval(nxR)
  return
  end subroutine PSvsE
  !-------------------------------------------------------------------
  subroutine PScalc(ie,ir,tag,bj,by,bjp,byp)
  implicit none
  character :: tag*2,idL*3,fil*127
  integer,parameter :: neqn=2,nwork=100+21*neqn,node=2000
  integer :: ie,ir,i,l,lmin,nor,nor1,flag,iwork(5)
  real(dp) :: bj(0:lmax+1),by(0:lmax+1),bjp(0:lmax),byp(0:lmax), &
    work(nwork),t,t1,or(node),ou(2,node),or1(node),ou1(2,node), &
    u(2),dudx(2),y(2),dydx(2),ps,ps1,A,cnst
  !
  iatom=ir
  zet=2.d0*z(ir)/clight
      if(tag=='su')then; lmin=0 
  elseif(tag=='sd')then; lmin=1 
  endif
  do l=lmin,lmax
        if(tag=='su')then; kappa=-l-1 
    elseif(tag=='sd')then; kappa=l
    endif
    gam=sqrt(dble(kappa**2)-zet**2)
    do i=3,nxR(ir)
      t1=rx(i,ir)
      if(t1**gam > rx(1,ir)) exit
    enddo
 !WAVE FUNCTION u = zero AT RADIUS r = t1 varying with orbit l.
    t=t1
    A=1.d0
    u(1)=A
    u(2)=A*(dble(kappa)+gam)/zet
    flag=1
    call ODE(waveq,neqn,u,t,rmt(ir),relerr,abserr,flag, &
             nwork,work,iwork,or,ou,nor,node)
    if(flag>=3)then
      write(611,'(a,i0,a,i0)')'ODE1: l=',l,', flag=',flag; stop
    endif
    !ODE step count.
    noru=max(noru,nor)
    !phase shift.
    call waveq(t,u,dudx)
    ps = atan((u(1)*(bj(l)+q*t*bjp(l)) - dudx(1)*t*bj(l))/ &
              (u(1)*(by(l)+q*t*byp(l)) - dudx(1)*t*by(l)))
        if(tag=='su')then; psu(l,ir)=ps 
    elseif(tag=='sd')then; psd(l,ir)=ps 
    endif
 !END WAVE FUNCTION u = zero.
 !
 !WAVE FUNCTION y = unity AT RADIUS r = t1 varying with orbit l.
    t=t1
    y=1.d0
    flag=1
    call ODE(waveq1,neqn,y,t,rmt(ir),relerr,abserr,flag, &
             nwork,work,iwork,or1,ou1,nor1,node)
    !ODE step count.
    nory=max(nory,nor1)
    !phase shift.
    call waveq1(t,y,dydx)
    u(1) = y(1)
    dudx(1) = (y(1)*gam/t + dydx(1))
    ps1 = atan((u(1)*(bj(l)+q*t*bjp(l)) - dudx(1)*t*bj(l))/ &
               (u(1)*(by(l)+q*t*byp(l)) - dudx(1)*t*by(l)))
        if(tag=='su')then; psu1(l,ir)=ps1 
    elseif(tag=='sd')then; psd1(l,ir)=ps1
    endif
 !END WAVE FUNCTION y = unity.
  !
  !wave function display.
    if(WF=='y'.and.ie==ne)then
      write(idL,'(i0)') l
      if(tag=='su') fil = '../'//trim(outdir)//'/u'//trim(idZA(ir))//'.'//idL
      if(tag=='sd') fil = '../'//trim(outdir)//'/v'//trim(idZA(ir))//'.'//idL 
      open(11,file=trim(fil),status='unknown')
      !ou normalized so that ou(1,nor) = interstitial wave function.
      cnst=t*(bj(l)*cos(ps)-by(l)*sin(ps))/ou(1,nor)
      do i=1,nor
        write(11,940) or(i),ou(1,i)*cnst,ou(2,i)*cnst
      enddo
      close(11)
    endif
  enddo !l
  return

  940 format(3es14.6)
  end subroutine PScalc
  !--------------------------------------------------------------------
  subroutine waveq(x,u,dudx)
  !WAVE FUNCTION u = zero AT RADIUS r = zero.
  implicit none
  integer  :: i
  real(dp) :: u(2),dudx(2),x,p,Vtotx,qbx,emvc
  real(dp),parameter :: a6=1.d0/6.d0
  !
  p=log(x/rx(1,iatom))/dx(iatom)+1.d0
  i=nint(p)
  i=max(i,2); i=min(i,nxR(iatom)-2)
  p=p-i
  !Vtotx = Vtot * x for convenient 4-point interpolation.
  Vtotx=   -a6*p*(p-1.d0)*(p-2.d0)*Vtot(i-1,iatom)&
        +0.5d0*(p*p-1.d0)*(p-2.d0)*Vtot(i  ,iatom)&
        -0.5d0*p*(p+1.d0)*(p-2.d0)*Vtot(i+1,iatom)&
                  +a6*p*(p*p-1.d0)*Vtot(i+2,iatom)
  !Vtotx/x = back-division by x.
  emvc=(Einc-Vtotx/x)/clight
  qbx=dble(kappa)/x
  dudx(1)=-u(1)*qbx+u(2)*(clight+emvc)
  dudx(2)= u(2)*qbx-u(1)*emvc
  return
  end subroutine waveq
  !--------------------------------------------------------------------
  subroutine waveq1(x,y,dydx)
  !WAVE FUNCTION y = unity AT RADIUS r = zero.
  implicit none
  integer  :: i
  real(dp) :: y(2),dydx(2),x,p,Vtotx,emvc
  real(dp),parameter :: a6=1.d0/6.d0
  !
  p=log(x/rx(1,iatom))/dx(iatom)+1.d0
  i=nint(p)
  i=max(i,2); i=min(i,nxR(iatom)-2)
  p=p-i
  !Vtotx = Vtot * x for convenient 4-point interpolation.
  Vtotx=   -a6*p*(p-1.d0)*(p-2.d0)*Vtot(i-1,iatom)&
        +0.5d0*(p*p-1.d0)*(p-2.d0)*Vtot(i  ,iatom)&
        -0.5d0*p*(p+1.d0)*(p-2.d0)*Vtot(i+1,iatom)&
                  +a6*p*(p*p-1.d0)*Vtot(i+2,iatom)
  !Vtotx/x = back-division by x.
  emvc = (Einc-Vtotx/x)/clight
  dydx(1) = (dble(kappa)+gam) * (-y(1)/x + y(2)*(clight+emvc)/zet)
  dydx(2) = (dble(kappa)-gam) * ( y(2)/x - y(1)*emvc/zet)
  return
  end subroutine waveq1
  !--------------------------------------------------------------------
  subroutine sph_bessel_j(lmax,x,bj)
  !SPH_BESSEL_J is spherical bessel function j(l,x).
  implicit none
  integer :: l,ll,lx,lmax,llp1,llp3
  real(dp)  :: bj(0:lmax+1),x,xi,xl,x2h,x2hh,cl,f,s,almax
  real(dp),allocatable :: aux(:)
  almax=lmax
  !small arguments: limit determined by machine precision 2.22E-16 .
  if(x<1.0188d-02)then
    x2h=x*x*0.5d0; x2hh=x2h*0.5d0
    bj(0)=1.d0-(x2h/3.d0)*(1.d0-x2hh/5.d0)
    xl=1.d0; cl=1.d0
    do l=1,lmax+1
      xl=xl*x; ll=l+l; cl=cl/(ll+1)
      bj(l)=cl*xl*(1.d0-(x2h/(ll+3))*(1.d0-x2hh/(ll+5)))
    enddo
    return
  endif
  !large arguments.
  xi=1.d0/x
  if(x>almax)then
    !upward recurrence.
    bj(0)=sin(x)*xi
    if(lmax==0) return
    bj(1)=(bj(0)-cos(x))*xi
    if(lmax==1) return
    ll=3
    do l=2,lmax+1
      bj(l)=(ll*xi)*bj(l-1)-bj(l-2)
      ll=ll+2
    enddo
  else
    !downward recurrence.
    if(lmax<100)then
      lx=lmax+10+int(x*(1.d0-almax*0.0075d0))
    else
      lx=lmax+10+int(x*0.25d0)
    endif
    allocate(aux(0:lx)) !lx variable.
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
    bj(0:lmax+1)=aux(0:lmax+1)*f
  endif
  return
  end subroutine sph_bessel_j
  !--------------------------------------------------------------------
  subroutine sph_bessel_y(lmax,x,by)
  !SPH_BESSEL_Y is spherical bessel function y(l,x).
  implicit none
  integer :: l,ll,lmax
  real(dp)  :: by(0:lmax+1),x,xi
  xi=1.d0/x
  by(0)=-cos(x)*xi
  if(lmax==0) return
  by(1)=(by(0)-sin(x))*xi
  if(lmax==1) return
  ll=3
  do l=2,lmax+1
    by(l)=(ll*xi)*by(l-1)-by(l-2)
    ll=ll+2
  enddo
  return
  end subroutine sph_bessel_y
!=======================================================================
!subroutine ODE
!=======================================================================
!integrates a system of neqn first order ordinary differential eqs.:
!  dy(i)/dt = f(t,y(1:neqn))
!  y(i) given at an initial point t.
!
!Ref.: L.F. Shampine and M.K. Gordon, Computer solution of ordinary 
!      differential equations: the initial value problem (Freeman,1975).
!Ref.: http//www.netlib.org, search for ode/ode.f.
!
!the subroutine integrates from  t  to  tout .  on return the
!parameters in the call list are set for continuing the integration.
!the user has only to define a unknown value  tout  and call  ode  again.
!
!the differential equations are actually solved by a suite of codes
!de ,  step , and  intrp .  ode  allocates virtual storage in the
!arrays  work  and  iwork  and calls  de .  de  is a supervisor which
!directs the solution.  it calls on the routines  step  and  intrp
!to advance the integration and to interpolate at output points.
!step  uses a modified divided difference form of the adams pece
!formulas and local extrapolation.  it adjusts the order and step
!size to control the local error per unit step in a generalized
!sense.  normally each call to  step  advances the solution one step
!in the direction of  tout .  for reasons of efficiency  de
!integrates beyond  tout  internally, though never beyond
!t+10*(tout-t), and calls  intrp  to interpolate the solution at
!tout .  an option is provided to stop the integration at  tout  but
!it should be used only if it is impossible to continue the
!integration beyond  tout .
!
!the parameters represent:
!   f -- double precision subroutine f(t,y,yp) to evaluate
!             derivatives yp(i)=dy(i)/dt
!   neqn (integer*4)-- number of equations to be integrated,
!   y(*) (real(dp))-- solution vector at t,                
!   t (real(dp))-- independent variable,                    
!   tout (real(dp))-- point at which solution is desired,   
!   relerr,abserr (real(dp))-- relative and absolute error tolerances for
!        local error test .  at each step the code requires
!          dabs(local error) .le. dabs(y)*relerr + abserr
!        for each component of the local error and solution vectors,
!   iflag (integer*4)-- indicates status of integration,
!   work(*) (real(dp))  -- arrays to hold information internal to
!   iwork(*) (integer*4)    which is necessary for subsequent calls,
!   OR(*) (real(dp)) -- integration points chosen by the program,
!   OU(*) (real(dp)) -- the corresponding values of y(1),
!   NOR (integer*4) -- the number of integration points,
!   NODE (integer*4) -- dimension of OR and OU.
!
!first call to ode:
!the user must provide storage in his calling program for the arrays
!in the call list,
!   y(neqn), work(100+21*neqn), iwork(5),
!declare  f  in an external statement, supply the double precision
!subroutine f(t,y,yp)  to evaluate
!   dy(i)/dt = yp(i) = f(t,y(1),y(2),...,y(neqn))
!and initialize the parameters:
!   neqn  -- number of equations to be integrated
!   y(*)  -- vector of initial conditions
!   t     -- starting point of integration
!   tout  -- point at which solution is desired
!   relerr,abserr -- relative and absolute local error tolerances
!   iflag -- +1,-1.  indicator to initialize the code.  normal input
!            is +1.  the user should set iflag=-1 only if it is
!            impossible to continue the integration beyond  tout .
!all parameters except  f ,  neqn  and  tout  may be altered by the
!code on output so must be variables in the calling program.
!
!output from ode:
!   neqn -- unchanged
!   y(*) -- solution at  t
!   t    -- last point reached in integration. normal return has t=tout.
!   tout -- unchanged
!   relerr,abserr -- normal return has tolerances unchanged.  iflag=3
!        signals tolerances increased
!   iflag = 2 -- normal return.  integration reached  tout
!         = 3 -- integration did not reach  tout  because error
!                tolerances too small.  relerr ,  abserr  increased
!                appropriately for continuing
!         = 4 -- integration did not reach  tout  because more than
!                NODE steps needed
!         = 5 -- integration did not reach  tout  because equations
!                appear to be stiff
!         = 6 -- invalid input parameters (fatal error)
!        the value of  iflag  is returned negative when the input
!        value is negative and the integration does not reach  tout ,
!        i.e., -3, -4, -5.
!   work(*),iwork(*) -- information generally of no interest to the
!                       user but necessary for subsequent calls.
!
!subsequent calls to ode:
!subroutine  ode  returns with all information needed to continue
!the integration.  if the integration reached  tout , the user need
!only define a unknown  tout  and call again.  if the integration did not
!reach  tout  and the user wants to continue, he just calls again.
!the output value of  iflag  is the appropriate input value for
!subsequent calls.  the only situation in which it should be altered
!is to stop the integration internally at the unknown  tout , i.e.,
!change output  iflag=2  to input  iflag=-2 .  error tolerances may
!be changed by the user before continuing.  all other parameters must
!remain unchanged.
subroutine ode(f,neqn,y,t,tout,relerr,abserr,iflag,nwork,work,&
               iwork,OR,OU,NOR,NODE)
      implicit none
      logical :: start,phase1,nornd      
      integer :: neqn,nwork,NODE,iyy,iwt,ip,iyp,iypout,iphi,NOR,&
        istart,iphase,iwork(5),ialpha,ibeta,isig,iv,iw,ig,ipsi,ix,ih,&
        ihold,itold,idelsn,iflag
      real(dp) :: y(neqn),t,tout,relerr,abserr,work(nwork),&
        OR(NODE),OU(2,NODE),twou,fouru
      external f
      data ialpha,ibeta,isig,iv,iw,ig,iphase,ipsi,ix,ih,ihold,istart,&
        itold,idelsn/1,13,25,38,50,62,75,76,88,89,90,91,92,93/
      twou=2.d0*epsilon(1.d0); fouru=2.d0*twou
      iyy = 100
      iwt = iyy + neqn
      ip = iwt + neqn
      iyp = ip + neqn
      iypout = iyp + neqn
      iphi = iypout + neqn
      NOR=0
      if(iabs(iflag) .eq. 1) go to 1
      start = work(istart) .gt. 0.0d0
      phase1 = work(iphase) .gt. 0.0d0
      nornd = iwork(2) .ne. -1
    1 call de(f,neqn,y,t,tout,relerr,abserr,&
        iflag,work(iyy),work(iwt),work(ip),work(iyp),work(iypout),&
        work(iphi),work(ialpha),work(ibeta),work(isig),work(iv),&
        work(iw),work(ig),phase1,work(ipsi),work(ix),work(ih),&
        work(ihold),start,work(itold),work(idelsn),iwork(1),nornd,&
        iwork(3),iwork(4),iwork(5),OR,OU,NOR,NODE,twou,fouru)
      work(istart) = -1.0d0
      if(start) work(istart) = 1.0d0
      work(iphase) = -1.0d0
      if(phase1) work(iphase) = 1.0d0
      iwork(2) = -1
      if(nornd) iwork(2) = 1
      return
  end subroutine ode
!-----------------------------------------------------------------------
subroutine de(f,neqn,y,t,tout,relerr,abserr,iflag,yy,wt,p,yp,&
              ypout,phi,alpha,beta,sig,v,w,g,phase1,psi,x,h,hold,start,&
              told,delsgn,ns,nornd,k,kold,isnold,OR,OU,NOR,NODE,&
              twou,fouru)

!ode  merely allocates storage for  de  to relieve the user of the
!inconvenience of a long call list.
!NODE  is the maximum number of steps allowed in one call to ode .
      implicit none
      logical :: stiff,crash,start,phase1,nornd      
      integer  :: neqn,NODE,iflag,isn,kle4,isnold,NOR,kold,k,ns
      real(dp) :: t,tout,relerr,abserr,eps,told,del,absdel,tend,&
        releps,abseps,delsgn,x,h,twou,fouru,hold
      real(dp) :: y(neqn),yy(neqn),wt(neqn),phi(neqn,16),p(neqn),&
        yp(neqn),ypout(neqn),psi(12),alpha(12),beta(12),sig(13),v(12),&
        w(12),g(13),OR(NODE),OU(2,NODE)
      external f
!test for improper parameters.
      if(neqn .lt. 1) go to 10
      !if(t .eq. tout) go to 10  !JR
      if(abs(t-tout)<1.d-14) goto 10
      if(relerr .lt. 0.0d0  .or.  abserr .lt. 0.0d0) go to 10
      eps = dmax1(relerr,abserr)
      if(eps .le. 0.0d0) go to 10
      if(iflag .eq. 0) go to 10
      isn = isign(1,iflag)
      iflag = iabs(iflag)
      if(iflag .eq. 1) go to 20
      !if(t .ne. told) go to 10  !JR
      if(abs(t-told)>1.d-14) goto 10
      if(iflag .ge. 2  .and.  iflag .le. 5) go to 20
   10 iflag = 6
      return
!on each call set interval of integration and counter for number of
!steps. adjust input error tolerances to define weight vector for
!subroutine  step.
   20 continue
      del = tout - t
      absdel = dabs(del)
      tend = t + 10.0d0*del
      if(isn .lt. 0) tend = tout
      kle4 = 0
      stiff = .false.
      releps = relerr/eps
      abseps = abserr/eps
      if(iflag .eq. 1) go to 30
      if(isnold .lt. 0) go to 30
      if(delsgn*del .gt. 0.0d0) goto 50
!on start and restart also set work variables x and yy(*), store the
!direction of integration and initialize the step size.
   30 continue
      start = .true.
      x = t
      yy(1:neqn)=y(1:neqn)
      delsgn = dsign(1.0d0,del)
      h = dsign(dmax1(dabs(tout-x),fouru*dabs(x)),tout-x)
      if(NOR.eq.0)then
        NOR=NOR+1; OR(NOR)=x; OU(1:neqn,NOR)=y(1:neqn)
      endif
!if already past output point, interpolate and return.
   50 continue
      if(dabs(x-t) .lt. absdel) goto 60
      call intrp(x,yy,tout,y,ypout,neqn,kold,phi,psi)
      !NOR is not incremented when OR and OU replace overshoot values.
      OR(NOR)=tout; OU(1:neqn,NOR)=y(1:neqn)
      iflag = 2
      t = tout
      told = t
      isnold = isn
      return
!if cannot go past output point and sufficiently close,
!extrapolate and return.
   60 continue
      if(isn .gt. 0  .or.  dabs(tout-x) .ge. fouru*dabs(x)) goto 80
      h = tout - x
      call f(x,yy,yp)
      y(1:neqn)=yy(1:neqn)+h*yp(1:neqn)
      iflag = 2
      t = tout
      told = t
      isnold = isn
      return
!test for too many steps.
   80 continue
      if(NOR .lt. NODE-2) goto 100
      iflag = isn*4
      if(stiff) iflag = isn*5
      write(611,'(a,i0,a,i0)')&
           'ode: more than ',NODE,' steps, iflag=',iflag
      y(1:neqn)=yy(1:neqn)
      t = x
      told = t
      isnold = 1
      return
!limit step size, set weight vector and take a step.
  100 continue
      h = dsign(dmin1(dabs(h),dabs(tend-x)),h)
      wt(1:neqn) = releps*dabs(yy(1:neqn)) + abseps
      call step(x,yy,f,neqn,h,eps,wt,start,hold,k,kold,crash,phi,p,yp,&
                psi,alpha,beta,sig,v,w,g,phase1,ns,nornd,twou,fouru)
      NOR=NOR+1; OR(NOR)=x; OU(1:neqn,NOR)=yy(1:neqn)
!test for tolerances too small.
      if(.not.crash) go to 130
      iflag = isn*3
      relerr = eps*releps
      abserr = eps*abseps
      y(1:neqn)=yy(1:neqn)
      t = x
      told = t
      isnold = 1
      return
!augment counter on number of steps and test for stiffness.
  130 continue
      kle4 = kle4 + 1
      if(kold .gt. 4) kle4 = 0
      if(kle4 .ge. 50) stiff = .true.
      goto 50
end subroutine de
!-----------------------------------------------------------------------
!SUBROUTINE STEP integrates a system of first order ordinary
!differential equations one step, normally from x to x+h, using a
!modified divided difference form of the adams pece formulas.  local
!extrapolation is used to improve absolute stability and accuracy.
!the code adjusts its order and step size to control the local error
!per unit step in a generalized sense.  special devices are included
!to control roundoff error and to detect when the user is requesting
!too much accuracy.
!
!the parameters represent:
!   x     -- independent variable (real(dp))
!   y(*)  -- solution vector at x (real(dp))
!   yp(*) -- derivative of solution vector at  x  after successful
!            step (real(dp))
!   neqn  -- number of equations to be integrated (integer*4)
!   h     -- appropriate step size for next step. normally determined by
!            code (real(dp))
!   eps   -- local error tolerance.  must be variable (real(dp))
!   wt(*) -- vector of weights for error criterion (real(dp))
!   start -- logical variable set .true. for first step,  .false.
!            otherwise (logical*4)
!   hold  -- step size used for last successful step (real(dp))
!   k     -- appropriate order for next step (determined by code)
!   kold  -- order used for last successful step
!   crash -- logical variable set .true. when no step can be taken,
!            .false. otherwise.
!the arrays  phi, psi  are required for the interpolation subroutine
!intrp.  the array p is internal to the code.  all are real(dp)
!
!input to  step:
!
!   first call:
!the user must provide storage in his driver program for all arrays
!in the call list, namely,
!   y(neqn),wt(neqn),phi(neqn,16),p(neqn),yp(neqn),psi(12)
!
!the user must also declare  start  and  crash  logical variables
!and  f  an external subroutine, supply the subroutine  f(x,y,yp)
!to evaluate
!   dy(i)/dx = yp(i) = f(x,y(1),y(2),...,y(neqn))
!and initialize only the following parameters:
!   x     -- initial value of the independent variable
!   y(*)  -- vector of initial values of dependent variables
!   neqn  -- number of equations to be integrated
!   h     -- nominal step size indicating direction of integration
!            and maximum size of step.  must be variable
!   eps   -- local error tolerance per step.  must be variable
!   wt(*) -- vector of non-zero weights for error criterion
!   start -- .true.
!
!step  requires the l2 norm of the vector with components
!local error(l)/wt(l)  be less than  eps  for a successful step.  the
!array  wt  allows the user to specify an error test appropriate
!for his problem.  for example,
!   wt(l) = 1.0  specifies absolute error,
!         = dabs(y(l))  error relative to the most recent value of
!              the l-th component of the solution,
!         = dabs(yp(l))  error relative to the most recent value of
!              the l-th component of the derivative,
!         = dmax1(wt(l),dabs(y(l)))  error relative to the largest
!              magnitude of l-th component obtained so far,
!         = dabs(y(l))*relerr/eps + abserr/eps  specifies a mixed
!              relative-absolute test where  relerr  is relative
!              error,  abserr  is absolute error and  eps =
!              dmax1(relerr,abserr) .
!
!   subsequent calls:
!subroutine  step  is designed so that all information needed to
!continue the integration, including the step size  h  and the order
!k , is returned with each step.  with the exception of the step
!size, the error tolerance, and the weights, none of the parameters
!should be altered.  the array  wt  must be updated after each step
!to maintain relative error tests like those above.  normally the
!integration is continued just beyond the desired endpoint and the
!solution interpolated there with subroutine  intrp .  if it is
!impossible to integrate beyond the endpoint, the step size may be
!reduced to hit the endpoint since the code will not take a step
!larger than the  h  input.  changing the direction of integration,
!i.e., the sign of  h , requires the user set  start = .true. before
!calling  step  again.  this is the only situation in which  start
!should be altered.
!
!output from  step:
!
!   successful step:
!the subroutine returns after each successful step with  start  and
!crash  set .false. .  x  represents the independent variable
!advanced one step of length  hold  from its value on input and  y
!the solution vector at the unknown value of  x .  all other parameters
!represent information corresponding to the unknown  x  needed to
!continue the integration.
!
!   unsuccessful step:
!when the error tolerance is too small for the machine precision,
!the subroutine returns without taking a step and  crash = .true. .
!an appropriate step size and error tolerance for continuing are
!estimated and all other information is restored as upon input
!before returning.  to continue with the larger tolerance, the user
!just calls the code again.  a restart is neither required nor
!desirable.
subroutine step(x,y,f,neqn,h,eps,wt,start,hold,k,kold,crash,phi,p,yp,&
                psi,alpha,beta,sig,v,w,g,phase1,ns,nornd,twou,fouru)
      implicit none
      logical  :: start,crash,phase1,nornd      
      integer  :: neqn,l,k,kold,ifail,kp1,kp2,km1,km2,ns,nsp1,i,im1,&
        iq,nsm2,j,limit1,nsp2,limit2,ip1,kunknown
      real(dp) :: y(neqn),wt(neqn),phi(neqn,16),p(neqn),yp(neqn),&
        psi(12),alpha(12),beta(12),sig(13),w(12),v(12),g(13),gstr(13),&
        two(13)
      real(dp) :: x,h,eps,hold,twou,fouru,p5eps,round,sum,absh,&
        realns,temp1,temp2,temp3,temp4,temp5,temp6,reali,tau,xold,erk,&
        erkm1,erkm2,erkp1,err,rho,hunknown,r
      external f
      data two/2.0d0,4.0d0,8.0d0,16.0d0,32.0d0,64.0d0,128.0d0,256.0d0,&
        512.0d0,1024.0d0,2048.0d0,4096.0d0,8192.0d0/
      data gstr/0.500d0,0.0833d0,0.0417d0,0.0264d0,0.0188d0,0.0143d0,&
        0.0114d0,0.00936d0,0.00789d0,0.00679d0,0.00592d0,0.00524d0,&
        0.00468d0/

!begin BLOCK 0
!check if step size/error tolerance is too small for machine precision.
!if first step, initialize phi array and estimate a starting step size.
!
!if step size is too small, determine an acceptable one.
      crash = .true.
      if(dabs(h) .ge. fouru*dabs(x)) go to 5
      h = dsign(fouru*dabs(x),h)
      return
    5 p5eps = 0.5d0*eps
!if error tolerance is too small, increase it to an acceptable value.
      round = 0.0d0
      do 10 l = 1,neqn
        round = round + (y(l)/wt(l))**2
   10 continue
      round = twou*dsqrt(round)
      if(p5eps .ge. round) go to 15
      eps = 2.0*round*(1.0d0 + fouru)
      return
   15 crash = .false.
      g(1)=1.0d0
      g(2)=0.5d0
      sig(1)=1.0d0
      if(.not.start) go to 99
!initialize. compute appropriate step size for first step.
      call f(x,y,yp)
      sum = 0.0d0
      do 20 l = 1,neqn
        phi(l,1) = yp(l)
        phi(l,2) = 0.0d0
        sum = sum + (yp(l)/wt(l))**2
   20 continue
      sum = dsqrt(sum)
      absh = dabs(h)
      if(eps .lt. 16.0d0*sum*h*h) absh = 0.25d0*dsqrt(eps/sum)
      h = dsign(dmax1(absh,fouru*dabs(x)),h)
      hold = 0.0d0
      k = 1
      kold = 0
      start = .false.
      phase1 = .true.
      nornd = .true.
      if(p5eps .gt. 100.0d0*round) go to 99
      nornd = .false.
      do 25 l = 1,neqn
        phi(l,15) = 0.0d0
   25 continue
   99 ifail = 0
!end BLOCK 0
!
!begin BLOCK 1
!compute coefficients of formulas for this step. avoid computing
!those quantities not changed when step size is not changed.
!                  ***
  100 kp1 = k+1
      kp2 = k+2
      km1 = k-1
      km2 = k-2
!ns is the number of steps taken with size h, including the current one.
!when k.lt.ns, no coefficients change.
      !if(h .ne. hold) ns = 0 !JR
      if(abs(h-hold)>1.d-14) ns=0
      if(ns.le.kold)   ns=ns+1
      nsp1 = ns+1
      if (k .lt. ns) go to 199
!compute those components of alpha(*),beta(*),psi(*),sig(*) which are
!changed.
      beta(ns) = 1.0d0
      realns = ns
      alpha(ns) = 1.0d0/realns
      temp1 = h*realns
      sig(nsp1) = 1.0d0
      if(k .lt. nsp1) go to 110
      do 105 i = nsp1,k
        im1 = i-1
        temp2 = psi(im1)
        psi(im1) = temp1
        beta(i) = beta(im1)*psi(im1)/temp2
        temp1 = temp2 + h
        alpha(i) = h/temp1
        reali = i
        sig(i+1) = reali*alpha(i)*sig(i)
  105 continue
  110 psi(k) = temp1
!compute coefficients g(*). initialize v(*) and set w(*).
!g(2) is set in data statement.
      if(ns .gt. 1) go to 120
      do iq = 1,k
        temp3 = iq*(iq+1)
        v(iq) = 1.0d0/temp3
        w(iq) = v(iq)
      enddo
      go to 140
!if order was raised, update diagonal part of v(*).
  120 if(k .le. kold) go to 130
      temp4 = k*kp1
      v(k) = 1.0d0/temp4
      nsm2 = ns-2
      if(nsm2 .lt. 1) go to 130
      do 125 j = 1,nsm2
        i = k-j
        v(i) = v(i) - alpha(j+1)*v(i+1)
  125 continue
!update v(*) and set w(*).
  130 limit1 = kp1 - ns
      temp5 = alpha(ns)
      do 135 iq = 1,limit1
        v(iq) = v(iq) - temp5*v(iq+1)
        w(iq) = v(iq)
  135 continue
      g(nsp1) = w(1)
!compute the g(*) in the work vector w(*).
  140 nsp2 = ns + 2
      if(kp1 .lt. nsp2) go to 199
      do 150 i = nsp2,kp1
        limit2 = kp2 - i
        temp6 = alpha(i-1)
        do 145 iq = 1,limit2
          w(iq) = w(iq) - temp6*w(iq+1)
  145   continue
        g(i) = w(1)
  150 continue
  199 continue
!end BLOCK 1
!
!begin BLOCK 2
!predict a solution p(*), evaluate derivatives using predicted solution,
!estimate local error at order k and errors at orders k, k-1, k-2
!as if constant step size were used.
!
!change phi to phi star.
      if(k .lt. nsp1) go to 215
      do 210 i = nsp1,k
        temp1 = beta(i)
        do 205 l = 1,neqn
          phi(l,i) = temp1*phi(l,i)
  205   continue
  210 continue
!predict solution and differences.
  215 do 220 l = 1,neqn
        phi(l,kp2) = phi(l,kp1)
        phi(l,kp1) = 0.0d0
        p(l) = 0.0d0
  220 continue
      do 230 j = 1,k
        i = kp1 - j
        ip1 = i+1
        temp2 = g(i)
        do 225 l = 1,neqn
          p(l) = p(l) + temp2*phi(l,i)
          phi(l,i) = phi(l,i) + phi(l,ip1)
  225   continue
  230 continue
      if(nornd) go to 240
      do 235 l = 1,neqn
        tau = h*p(l) - phi(l,15)
        p(l) = y(l) + tau
        phi(l,16) = (p(l) - y(l)) - tau
  235 continue
      go to 250
  240 do 245 l = 1,neqn
        p(l) = y(l) + h*p(l)
  245 continue
  250 xold = x
      x = x + h
      absh = dabs(h)
      call f(x,p,yp)
!estimate errors at orders k,k-1,k-2.
      erkm2 = 0.0d0
      erkm1 = 0.0d0
      erk = 0.0d0
      do 264 l = 1,neqn
        temp3 = 1.0d0/wt(l)
        temp4 = yp(l) - phi(l,1)
!!      if(km2)265,260,255
!!255   erkm2 = erkm2 + ((phi(l,km1)+temp4)*temp3)**2
!!260   erkm1 = erkm1 + ((phi(l,k)+temp4)*temp3)**2
!!265   erk = erk + (temp4*temp3)**2
        if(km2>0)then
          erkm2 = erkm2 + ((phi(l,km1)+temp4)*temp3)**2
          erkm1 = erkm1 + ((phi(l,k)+temp4)*temp3)**2
        elseif(km2==0)then
          erkm1 = erkm1 + ((phi(l,k)+temp4)*temp3)**2
        endif
        erk = erk + (temp4*temp3)**2
  264 continue
!!    if(km2)280,275,270
!!270 erkm2 = absh*sig(km1)*gstr(km2)*dsqrt(erkm2)
!!275 erkm1 = absh*sig(k)*gstr(km1)*dsqrt(erkm1)
!!280 temp5 = absh*dsqrt(erk)
      if(km2>0)then
        erkm2 = absh*sig(km1)*gstr(km2)*dsqrt(erkm2)
        erkm1 = absh*sig(k)*gstr(km1)*dsqrt(erkm1)
      elseif(km2==0)then
        erkm1 = absh*sig(k)*gstr(km1)*dsqrt(erkm1)
      endif
      temp5 = absh*dsqrt(erk)      
      err = temp5*(g(k)-g(kp1))
      erk = temp5*sig(kp1)*gstr(k)
      kunknown = k
!test if order should be lowered.
!!    if(km2)299,290,285
!!285 if(dmax1(erkm1,erkm2) .le. erk) kunknown = km1
!!    go to 299
!!290 if(erkm1 .le. 0.5d0*erk) kunknown = km1
!!299 if(err .le. eps) go to 400
      if(km2>0)then
        if(dmax1(erkm1,erkm2) .le. erk) kunknown = km1
      elseif(km2==0)then
        if(erkm1 .le. 0.5d0*erk) kunknown = km1
      endif
      if(err .le. eps) go to 400
!test if step successful.  
!end BLOCK 2
!
!begin BLOCK 3
!the step is unsuccessful.  restore  x, phi(*,*), psi(*) .
!if third consecutive failure, set order to one. if step fails more than
!three times, consider an optimal step size.  double error tolerance and
!return, if estimated step size is too small for machine precision.
!
!restore x, phi(*,*) and psi(*).
      phase1 = .false.
      x = xold
      do 310 i = 1,k
        temp1 = 1.0d0/beta(i)
        ip1 = i+1
        do 305 l = 1,neqn
          phi(l,i) = temp1*(phi(l,i) - phi(l,ip1))
  305   continue
  310 continue
      if(k .lt. 2) go to 320
      do 315 i = 2,k
        psi(i-1) = psi(i) - h
  315 continue
!on third failure, set order to one.  thereafter, use optimal step size.
  320 ifail = ifail + 1
      temp2 = 0.5d0
!!    if(ifail - 3) 335,330,325
!!325 if(p5eps .lt. 0.25d0*erk) temp2 = dsqrt(p5eps/erk)
!!330 kunknown = 1
!!335 h = temp2*h
      if(ifail-3>0)then
        if(p5eps .lt. 0.25d0*erk) temp2 = dsqrt(p5eps/erk)
        kunknown = 1
      elseif(ifail-3==0)then
        kunknown=1
      endif
      h = temp2*h
      k = kunknown
      if(dabs(h) .ge. fouru*dabs(x)) go to 340
      crash = .true.
      h = dsign(fouru*dabs(x),h)
      eps = eps + eps
      return
  340 go to 100
!end BLOCK 3
!
!begin BLOCK 4
!the step is successful.  correct the predicted solution, evaluate
!the derivatives using the corrected solution and update the
!differences. determine best order and step size for next step.
  400 kold = k
      hold = h
!correct and evaluate.
      temp1 = h*g(kp1)
      if(nornd) go to 410
      do l = 1,neqn
        rho = temp1*(yp(l) - phi(l,1)) - phi(l,16)
        y(l) = p(l) + rho
        phi(l,15) = (y(l) - p(l)) - rho
      enddo
      go to 420
  410 do l=1,neqn
        y(l)=p(l)+temp1*(yp(l)-phi(l,1))
      enddo
  420 continue
      call f(x,y,yp)
!update differences for next step.
      do l=1,neqn
        phi(l,kp1)=yp(l)-phi(l,1)
        phi(l,kp2)=phi(l,kp1)-phi(l,kp2)
      enddo
      do i=1,k
        do l=1,neqn
          phi(l,i)=phi(l,i)+phi(l,kp1)
        enddo
      enddo
!estimate error at order k+1 unless:
!in first phase when always raise order,
!already decided to lower order,
!step size not constant so estimate unreliable.
      erkp1 = 0.0d0
      if(kunknown .eq. km1  .or.  k .eq. 12) phase1 = .false.
      if(phase1) go to 450
      if(kunknown .eq. km1) go to 455
      if(kp1 .gt. ns) go to 460
      do 440 l = 1,neqn
        erkp1 = erkp1 + (phi(l,kp2)/wt(l))**2
  440 continue
      erkp1 = absh*gstr(kp1)*dsqrt(erkp1)
!using estimated error at order k+1, determine appropriate order for
!next step.
      if(k .gt. 1) go to 445
      if(erkp1 .ge. 0.5d0*erk) go to 460
      go to 450
  445 if(erkm1 .le. dmin1(erk,erkp1)) go to 455
      if(erkp1 .ge. erk  .or.  k .eq. 12) go to 460
!here erkp1 .lt. erk .lt. dmax1(erkm1,erkm2) else order would have been
!lowered in block 2. thus order is to be raised. raise order.
  450 k = kp1
      erk = erkp1
      go to 460
!lower order.
  455 k = km1
      erk = erkm1
!with unknown order determine appropriate step size for next step.
  460 hunknown = h + h
      if(phase1) go to 465
      if(p5eps .ge. erk*two(k+1)) go to 465
      hunknown = h
      if(p5eps .ge. erk) go to 465
      temp2 = k+1
      r = (p5eps/erk)**(1.0d0/temp2)
      hunknown = absh*dmax1(0.5d0,dmin1(0.9d0,r))
      hunknown = dsign(dmax1(hunknown,fouru*dabs(x)),h)
  465 h = hunknown
      return
!end BLOCK 4
end subroutine step
!-----------------------------------------------------------------------
subroutine intrp(x,y,xout,yout,ypout,neqn,kold,phi,psi)
!
!the methods in subroutine  step  approximate the solution near  x
!by a polynomial.  subroutine  intrp  approximates the solution at
!xout  by evaluating the polynomial there.
!
!input to intrp:
!   the user provides storage in the calling program for the arrays 
!       y(neqn),yout(neqn),ypout(neqn),phi(neqn,16),psi(12)
!   xout     -- point at which solution is desired.
!
!output from  intrp:
!   yout(*)  -- solution at  xout
!   ypout(*) -- derivative of solution at  xout
!the remaining parameters are returned unaltered from their input
!values.  integration with  step  may be continued.
      implicit none
      integer :: neqn,ki,kold,kip1,i,j,jm1,limit1,l
      real(dp) :: y(neqn),yout(neqn),ypout(neqn),phi(neqn,16),&
        psi(12),g(13),w(13),rho(13),x,xout,hi,term,temp1,temp2,temp3,&
        psijm1,gamma,eta
      data g(1)/1.0d0/,rho(1)/1.0d0/
      hi = xout - x
      ki = kold + 1
      kip1 = ki + 1
!initialize w(*) for computing g(*)
      do 5 i = 1,ki
        temp1 = i
        w(i) = 1.0d0/temp1
    5 continue
      term = 0.0d0
!compute g(*)
      do 15 j = 2,ki
        jm1 = j - 1
        psijm1 = psi(jm1)
        gamma = (hi + term)/psijm1
        eta = hi/psijm1
        limit1 = kip1 - j
        do 10 i = 1,limit1
          w(i) = gamma*w(i) - eta*w(i+1)
   10   continue
        g(j) = w(1)
        rho(j) = gamma*rho(jm1)
        term = psijm1
   15 continue
!interpolate
      do 20 l = 1,neqn
        ypout(l) = 0.0d0
        yout(l) = 0.0d0
   20 continue
      do 30 j = 1,ki
        i = kip1 - j
        temp2 = g(i)
        temp3 = rho(i)
        do 25 l = 1,neqn
          yout(l) = yout(l) + temp2*phi(l,i)
          ypout(l) = ypout(l) + temp3*phi(l,i)
   25   continue
   30 continue
      do 35 l = 1,neqn
        yout(l) = y(l) + hi*yout(l)
   35 continue
      return
end subroutine intrp

end program eeas
