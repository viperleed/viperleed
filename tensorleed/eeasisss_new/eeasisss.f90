!=========================================================================
  !This program is free software under the terms of the GNU General Public
  !License as published by the Free Software Foundation.
  !Author: John O. Rundgren, jru@KTH.se ,
  !        KTH Royal Institute of Technology, Stockholm, Sweden.
  !Version: 28 March 2021.
  !-----------------------------------------------------------------------
  subroutine EEASiSSS(input_file,log_file,output_dir,atom_dir)

  !Elastic Electron-Atom Scattering in Solids and Solid Surfaces,
  !author: John Rundgren, jru@kth.se .
  !The author appreciates acknowledgement in publications by citation:
  !J. Rundgren, Phys.Rev.B68,125405(2003),
  !             Phys.Rev.B76,195441(2007).
  !use ifport

  use param,only : dp,z,outdir,atdir,lmax,nthread,Chg,Pot,WF,SpinPS, &
    idZA,idElemA,relerr,abserr,linux,windows

  use maths,only : ipol,pois,iofx,sph_av

  use space,only : NeighborShells,NBR,nshell,nieq,neq,rx,dx,nx,nxx, &
    nlat,ad,ia,rmt,rmtS,rmin,rmax

  use enrgy,only : Free_Atom_Overlap,rho,rs,eev,ne,v0ev,emv0,fxc, &
    Vnist,Vcry,VcryR,Vcry0

  use sdata,only : nsp,nsr,sp,sr,sdat,sdat1

  !AMI: got rid of Eparallel
  !AMI instead use eeas from eeas.f90
  use eeas, only: read_eeas

  use DifferentialEvolution,only : DEM,NP,refresh,VTR,bestval,nfeval,itval

  implicit none
  character(len=*), intent(in) :: input_file,output_dir,atom_dir
  character(len=*) :: log_file
  character(len=2) :: ci
  character(len=4) :: cie
  character(len=127) :: txt
  integer :: iter,ie,ir,isys,i,j,l,NPfac
  integer,allocatable :: nxR(:)
  real(dp) :: v0coef(8),acc
  real(dp),allocatable :: Vxc0(:),XCmin(:),XCmax(:),bestmem_XC(:), &
    rmtSexp(:),psu(:,:,:),psd(:,:,:),psl(:,:,:),psu1(:,:,:),psd1(:,:,:), &
    eps(:,:),a(:)
  real(dp),parameter :: rydb=13.60569172d0
  external errMT
  !
  write(outdir,'(127(" "))')
  outdir = trim(output_dir)
  write(atdir,'(127(" "))')
  atdir = trim(atom_dir)
  !
  !OPEN LOG FILE.
  open(61,file=trim(outdir)//trim(log_file),status='unknown')
  write(61,'(2a)')'log_file: ',trim(outdir)//trim(log_file)
  write(61,'(a/80("-"))')'READ inputX:'
  !
  !INPUT.
  call StructureAndMethods(input_file)
  !
  !SCATTERING POTENTIALS.
  call NeighborShells
  !call NBR
  call Free_Atom_Overlap
  !
  !eev,emv0,v0ev allocated in StructureAndMethods.
  allocate(rmt(nieq),Vcry(nxx,nieq),VcryR(nieq), &
    Vxc0(ne),nxR(nieq),rmtSexp(nieq), &
    psu(ne,0:lmax,nieq),psd(ne,0:lmax,nieq),psl(ne,0:lmax,nieq), &
    psu1(ne,0:lmax,nieq),psd1(ne,0:lmax,nieq), &
    eps(0:lmax,nieq),a(ne) )
  !
  !COULOMBIC POTENTIAL INITIATED.
  write(61,'(/20("===="))')
  write(61,'(a)')'CRYSTAL POTENTIAL: MT radii and potential levels'
  write(61,'(20("----"))')
  do ir=1,nieq
    Vcry(1:nx(ir),ir) = Vnist(1:nx(ir),ir)
  enddo
  write(61,822)'rmax(B)',(rmax(ir),ir=1,nieq)
  write(61,822)'rmin(B)',(rmin(ir),ir=1,nieq)
  write(61,822)'rmtS',(rmtS(ir),ir=1,nieq)
  !DEM PARAMETRS.
  NPfac=max(10,nieq)
  NP=NPfac*nieq
  allocate(XCmin(nieq),XCmax(nieq),bestmem_XC(nieq))
  XCmax=rmax
  XCmin=rmin
  refresh=-1; VTR=0.0001d0
  !DO DEM LOOP.
  do iter=1,2
    write(61,'(/a,i0,a)')'eeasiss ITERATION ',iter
    call DEM(errMT,XCmin,XCmax,bestmem_XC,nieq)
    if(bestval>=1.d+33)then
      write(61,'(a)')'DEM failed, fitness = infty, stop.'; stop
    endif
    !NEW rmt.
    rmt=bestmem_XC
    do ir=1,nieq
      rmtSexp(ir)=(rmt(ir)+rmt(ia(2,ir)))/ad(2,ir)-1.d0
    enddo
    write(61,825)'fitness =',bestval*rydb
    write(61,810)'#calls,#improvements=',nfeval,itval
    write(61,820)'atom'  ,(idElemA(ir),ir=1,nieq)
    write(61,820)'NNatom',(idElemA(ia(2,ir)),ir=1,nieq)
    write(61,822)'rmt(B)',(rmt(ir),ir=1,nieq)
    write(61,822)'rmtNN(B)',(rmt(ia(2,ir)),ir=1,nieq)
    write(61,822)'rmtSexp',(rmtSexp(ir),ir=1,nieq)
    !Calculating VcryR and Vcry0.
    do ir=1,nieq
      VcryR(ir) = ipol(rmt(ir),nx(ir),rx(1,ir),dx(ir),Vcry(1,ir),'Vcry')
    enddo
    Vcry0 = sum(neq*VcryR)/dble(nlat)
    !Normalization to no MT-to-interstice potential steps.
    do ir=1,nieq
      Vcry(1:nx(ir),ir) = Vcry(1:nx(ir),ir) - (VcryR(ir) - Vcry0)
    enddo
  enddo !iter
  write(61,'(a)')'No MT-to-interstice potential steps.'
  !ENDDO DEM LOOP.
  !
  !Carrier potential Vcry (E invariant) normalized to Vcry0 = 0.
  do ir=1,nieq
    Vcry(1:nx(ir),ir) = Vcry(1:nx(ir),ir) - Vcry0
    VcryR(ir) = ipol(rmt(ir),nx(ir),rx(1,ir),dx(ir),Vcry(1,ir),'Vcry')
  enddo
  write(61,822)'VcryR(eV)',(VcryR(ir)*rydb,ir=1,nieq)
  Vcry0 = sum(neq*VcryR)/dble(nlat)
  write(61,823)'Vcry0(eV) = ',Vcry0*rydb
  write(61,*)
  810 format(a,2(1x,i0),a,es9.2,1x,es9.2)
  820 format(a/(10(3x,a)))
  822 format(a/(10f8.4:))
  922 format(a/(10i8:))
  823 format(a,f0.4,2(2x,f0.4))
  825 format(a,es7.0,a,es7.0)
  902 format(2es14.6)
  !COULOMBIC POTENTIAL FINISHED.
  !
  write(61,'(a)')'nxR(ir) enclosing rmt(ir).'
  do ir=1,nieq
    nxR(ir) = iofx(rmt(ir),rx(1,ir),1,nx(ir))
  enddo
  write(61,922)'nxR',(nxR(ir),ir=1,nieq)
  if(Pot=='y')then
    do ir=1,nieq
      open(10,file=trim(outdir)//'/'//'Vcry.'//idZA(ir),status='unknown')
      do i=1,nxR(ir)
        if(Vcry(i,ir)*rydb > -100.d0)then
          write(10,902) rx(i,ir),Vcry(i,ir)*rydb 
        endif
      enddo
    enddo
  endif
  if(Chg=='y')then
    do ir=1,nieq
    open(10,file=trim(outdir)//'/'//'rho.'//idZA(ir),status='unknown')
      do i=1,nxR(ir)
        write(10,902) rx(i,ir),rho(i,ir)
      enddo
    enddo
    do ir=1,nieq
    open(10,file=trim(outdir)//'/'//'rs.'//idZA(ir),status='unknown')
      do i=1,nxR(ir)
        write(10,902) rx(i,ir),rs(i,ir)
      enddo
    enddo
  endif
  !
  !SELF-ENERGY.
  sdat=reshape(sdat1,shape=(/nsp,nsr/))
goto 1001
  !plot of 'spx'//trim(ci) by means of xmgrace.
  do i=1,17
    write(ci,'(i0)') i
    open(62,file='spx'//trim(ci),status='unknown')    
    write(62,'(2f9.5)') (sp(j),sdat(j,i),j=1,151)
  enddo
  close(62)
1001 continue
  !
  write(61,'(/a)')'eeasisss: write files 1414 to eeas.'
  open(1414,file='uinp1',form='unformatted',status='unknown')
  write(1414) nxx,nieq,ne,nlat,nshell,nsp,nsr,nthread,lmax
  !
  open(1414,file='uinp2',form='unformatted',status='unknown')
  write(1414) outdir,Pot,WF,idElemA,idZA,neq,nx,dx,rx,ad,ia, &
    sp,sr,sdat,eev,Vcry,Vcry0,rmt,nxR,z,rho,rs,fxc,relerr,abserr 
  close(1414)
  write(61,'(a)')'eeasisss: files 1414 written.'
  !
  !ELASTIC ELECTRON-ATOM SCATTERING
  call read_eeas() ! Calls routines in eeas.f90
  !
  !EEASiSSS reads EEAS.
  write(61,'(/20("===="))')
  write(61,'(a)')'EEAS'
  write(61,'(20("----"))')
  !logfile.
  do ie=1,ne
    write(cie,'(i0)') ie
    !logfile.
    open(19,file='ulog'//trim(cie),status='unknown')
    do
      read(19,'(a)',end=1) txt
      write(61,'(a)') trim(txt)
    enddo
    1 continue
  enddo
  close(19)
  !data.
  do ie=1,ne
    write(cie,'(i0)') ie
    open(20,file='udat'//trim(cie),status='unknown',access='stream') 
    read(20) Vxc0(ie), &
             ((psu(ie,l,ir),l=0,lmax),ir=1,nieq), &
             ((psd(ie,l,ir),l=0,lmax),ir=1,nieq), &
             ((psu1(ie,l,ir),l=0,lmax),ir=1,nieq), &
             ((psd1(ie,l,ir),l=0,lmax),ir=1,nieq)
  enddo
  close(20)
  !
  write(61,'(/20("===="))')
  write(61,'(a)')'PS accuracy'
  write(61,'(20("----"))')
  write(61,'(a,es6.0,a,es6.0)') &
    'ODE tolerances, relerr=',relerr,' abserr=',abserr
  write(61,'(a,i2)')'#phase shifts=',lmax+1
  write(61,'(a)')'PS accuracy on subset ie,kappa:'
  write(61,'(a)')'ir'
!difference vs spin and ie.
  do ir=1,nieq
    do l=0,lmax
      do ie=1,ne
        a(ie)=0.d0
        a(ie)=max(a(ie),abs(psu(ie,l,ir))-abs(psu1(ie,l,ir)))
        a(ie)=max(a(ie),abs(psd(ie,l,ir))-abs(psd1(ie,l,ir)))
      enddo
      eps(l,ir)=maxval(a) 
    enddo
  enddo
  do ir=1,nieq
    write(61,'(/i2,16(1x,es6.0):)') ir,(eps(l,ir),l=0,lmax)
  enddo
  write(61,'(20("----"))')
!difference vs l and ir.
  acc=0.d0
  do ir=1,nieq
    do l=0,lmax
      acc=max(acc,eps(l,ir))
    enddo
  enddo
  write(61,'(a,1x,es6.0)')'PS accuracy on the set ir,l,ie,kappa:',acc
  write(61,'(20("----"))')
  !
  !ENERGY GRID in vacuum and in crystal.
  v0ev(:)=Vxc0(:)*rydb
  emv0(:)=eev(:)-v0ev(:)
  !
  !Vxc0Einc.
  open(63,file=trim(outdir)//'/Vxc0Einc',status='replace')
  write(63,930) (eev(i),v0ev(i),i=1,ne)
  write(61,'(/a)')'Vxc0Einc written.'
  write(61,*)
  !
  !Vxc0EincAprx.
  if(ne>=10)then
    call Vxc0EincAprx(Vxc0,v0coef)
    write(61,'(/a)')'Vxc0EincAprx written (Aprx for approximation).'
  else
    write(61,'(/a)')'ne<10, no Vxc0EincAprx calc.'
  endif
  930 format(2es14.6)
  !
  !PHASE SHIFTS.
  !Jumps of pi removed from the phaseshift versus energy curves.
  !psl = spin averaged = spinless phase shift.
  call PSnojump(psu,psd)
  do ir=1,nieq
    psl(1:ne,0,ir)=psu(1:ne,0,ir)
    do l=1,lmax
      do i=1,ne
        psl(i,l,ir)=((l+1)*psu(i,l,ir)+l*psd(i,l,ir))/(l+l+1)
      enddo
    enddo
  enddo
  !
  call PStab('sl',psl,v0coef)
  if(SpinPS=='y')then
    call PStab('su',psu,v0coef)
    call PStab('sd',psd,v0coef)
  endif
  write(61,'(/a)')'PStab written.'
  !
  !Cleaning up.
  if(windows)then ! AMI TODO: get rid of this mess...
    isys = SYSTEM('del uinp*')
    isys = SYSTEM('del ulog*')
    isys = SYSTEM('del udat*')
    isys = SYSTEM('del eeas.cmd')
    do i=1,nthread
      write(ci,'(i0)') i
      isys = SYSTEM('rmdir /S /Q '//trim('thread'//ci))
    enddo
  elseif(linux)then
    !isys = SYSTEM('rm uinp*') ! AMI: Clean up now accomplished in Python.
    !isys = SYSTEM('rm ulog*') ! AMI: Clean up now accomplished in Python.
    !isys = SYSTEM('rm udat*') ! AMI: Clean up now accomplished in Python.
    !isys = SYSTEM('rm eeas.sh') ! AMI: File no longer present
    do i=1,1
      write(ci,'(i0)') i
      !isys = SYSTEM('rm -r '//trim('thread'//ci))
    enddo
  endif
  !
  write(61,'(/a)')'eeasisss end stop.'
  stop
  end subroutine EEASiSSS
  !-------------------------------------------------------------------
  subroutine errMT(X,fitness)
  use param,only : dp
  use maths,only : ipol
  use space,only : nieq,neq,nlat,ad,ia,nx,rx,dx,rmtS
  use enrgy,only : Vcry,VcryR,Vcry0
  implicit none
  integer  :: ir
  real(dp) :: X(nieq),harvest,fitness
  intrinsic random_number
  !
  !rmtS in percent of ad(2,ir).
  do ir=1,nieq
    if(ad(2,ir)*(1.d0+rmtS(ir))-X(ir)-X(ia(2,ir)) < 0.d0)then
      call random_number(harvest)
      fitness=(1.d0+harvest)*1.d+33
      return
    endif
  enddo
  do ir=1,nieq
    VcryR(ir) = ipol(X(ir),nx(ir),rx(1,ir),dx(ir),Vcry(1,ir),'Vcry')
  enddo
  Vcry0 = sum(neq*VcryR)/dble(nlat)
  fitness = maxval(abs(VcryR-Vcry0))
  return
  end subroutine errMT
  !---------------------------------------------------------------------
  subroutine StructureAndMethods(input_file)
  use param,only : dp,Elem,z,compound,idZ,idA,idZA,idElemA,bos,SpinPS, &
    Chg,Pot,WF,lmax,nthread,eev1,eev2,emesh,relerr,abserr

  use space,only : rc,rk,dx,nx,nieq,neq,ieq,nlat,rmin,rmax,volUC,rmtS

  use enrgy,only : eev,ne,v0ev,emv0,fxc

  use DifferentialEvolution,only : method,itermax,strategy,F_XC,F_CR,CR_XC

  implicit none
  character(len=*), intent(in) :: input_file
  character(len=127) :: txt
  integer :: i,ir,j
  real(dp) :: UOL,rcc1,rcc2,rcc3,rcX(3,3)
  real(dp),allocatable :: rk_tmp(:,:)
  real(dp),parameter :: bohr=0.5291772083d0
  !  
 open(5,file=trim(input_file),status='unknown')
 read(5,*) txt(1:10)
   write(61,'(a)') txt(1:10)          !STRUCTURE:
 read(5,'(a)') compound
   compound=adjustl(trim(compound))   !CRYSTAL.
   write(61,'(a)') trim(compound)
 read(5,'(a)') txt
   write(61,'(a)') trim(txt)          !DIRECTORY.
 !
 !Bohr units of length inside crystal.
 !UOL = unit converter:
 !Bohr units = UOL * Researcher's units at input.
 read(5,*) UOL
   write(61,200) UOL
 200 format('UOL=',f10.6,', Bohr units of length inside crystal.')
 !
 !rc(i,j)=i'th coordinate of the j'th unit cell vector (B).
 write(61,'(a)')"UC vectors (input units)"
 do j=1,3
   read(5,*) rcX(1:3,j)
   write(61,'(3f12.7)') (rcX(i,j),i=1,3)
 enddo
 rc=rcX*UOL
 write(61,'(a)')'UC vectors (B)'
 write(61,'(3(3f12.7/))') ((rc(i,j),i=1,3),j=1,3)
 !
 !nieq    =number of ineqivalent atoms in unit cell,
 !neq(ir) =number of equivalent atoms of type ir,
 !nlat    =number of lattice points in unit cell,
 !z(ir)   =atomic number of type ir.
 read(5,*) nieq !
 !
 allocate(z(nieq),neq(nieq),dx(nieq),idZA(nieq),idElemA(nieq),nx(nieq), &
   rmin(nieq),rmax(nieq),ieq(nieq),rmtS(nieq),fxc(nieq), &
   rk_tmp(3,1000),idZ(nieq),idA(nieq))
 nlat=0
 !rmin, rmax in Bohr.
 do ir=1,nieq
   read(5,*) neq(ir),idZ(ir),idA(ir)
   z(ir)=dble(idZ(ir))
   !
   !atom identifiers.
   write(idZA(ir), '(i0,a,i0)') idZ(ir),'.',idA(ir)
   write(idElemA(ir),'( a,a,i0)') trim(Elem(idZ(ir))),'.',idA(ir)
   !
   write(61,203) ir,neq(ir)
   203 format(i0, 2x,i0, 2x,a)
   !
   !inequivalent lattice points.
   do j=1,neq(ir)
     nlat=nlat+1
     if(nlat>1000)then
       write(61,'(a)')'nlat>dimension 1000, stop'; stop
     endif
     read(5,*) rk_tmp(1:3,nlat)
     rk_tmp(1:3,nlat)=rk_tmp(1:3,nlat)*UOL
     write(61,205) rk_tmp(1:3,nlat)*bohr,rk_tmp(1:3,nlat)
   enddo
 enddo
 205 format(3f9.4,' Ang',2x,3f9.4,' Bohr')
 !
 do ir=1,nieq
   read(5,*)     idA(ir),rmin(ir),rmax(ir),rmtS(ir),fxc(ir)
   write(61,210) idA(ir),rmin(ir),rmax(ir),rmtS(ir),fxc(ir),&
                 trim(Elem(idZ(ir)))
 enddo
 210 format(i2,4f8.4,2x,a)
!
 !OPTIONS.
  read(5,*) txt
    write(61,'(a)') trim(txt)
  read(5,*) bos
    write(61,'(2a)')'bos==',bos
    if(.not.(bos=='b'.or.bos=='s'))goto 3
  read(5,*) SpinPS
    write(61,'(2a)')'SpinPS==',SpinPS
    if(.not.(SpinPS=='y'.or.SpinPS=='n'))goto 3    
  read(5,*) Chg
    write(61,'(2a)')'Chg==',Chg
    if(.not.(Chg=='y'.or.Chg=='n'))goto 3    
  read(5,*) Pot
    write(61,'(2a)')'Pot==',Pot
    if(.not.(Pot=='y'.or.Pot=='n'))goto 3
  read(5,*) WF
    write(61,'(2a)')'WF==',WF
    if(.not.(WF=='y'.or.WF=='n'))goto 3
  goto 4
  3 write(61,'(a)')'stop, incorrect input character'; stop
  4 continue
  !
  read(5,*) eev1,eev2,emesh
    write(61,'(3(a,f0.0))')'eev1=',eev1,', eev2=',eev2,', emesh=',emesh
  read(5,*) nthread,lmax
    write(61,'(2(a,i0))')'nthread=',nthread,', lmax=',lmax
  read(5,*) relerr,abserr
    write(61,'(2(a,es6.0))')'relerr=',relerr,', abserr=',abserr
  read(5,*)
    write(61,'(a)')'DIFFERENTIAL EVOLUTION'
  read(5,*) F_XC,CR_XC
    write(61,'(a,2f6.2)')'F_XC, CR_XC=',F_XC,CR_XC
  read(5,*) method
    write(61,'(a,3i2)')'method=',method
  read(5,*) strategy,F_CR
    write(61,'(a,i3,f6.2)')'strategy, F_CR=',strategy,F_CR
  read(5,*) itermax
    write(61,'(a,1x,i0)')'itermax=',itermax
  close(5)
 !
 write(61,'(/a,i3)')'#atoms in UC=',nlat
 allocate(rk(1:3,nlat))
 do j=1,nlat
   rk(1:3,j)=rk_tmp(1:3,j)
 enddo
 deallocate(rk_tmp)
 !volUC = volume of unit cell, where "unit cell" signifies:
 !for bulk, crystallographic unit cell;
 !for slab, crystallographic unit cell inclusive vacuum slab;
 !rcc1 and rcc2 || to the surface.
 rcc1=rc(2,1)*rc(3,2)-rc(3,1)*rc(2,2)
 rcc2=rc(3,1)*rc(1,2)-rc(1,1)*rc(3,2)
 rcc3=rc(1,1)*rc(2,2)-rc(2,1)*rc(1,2)
 volUC=abs(rc(1,3)*rcc1+rc(2,3)*rcc2+rc(3,3)*rcc3)
 write(61,'(a,f8.2)')'volUC (B**3)=',volUC
  !
  !UNIFORM ENERGY GRID.
  ne=1+nint((eev2-eev1)/emesh)
  allocate(eev(ne),v0ev(ne),emv0(ne))
  do i=1,ne
    eev(i)=eev1+(i-1)*emesh
  enddo
  return
 end subroutine StructureAndMethods
 !----------------------------------------------------------------------
  subroutine PSnojump(psu,psd)
  !Author: John Rundgren (2016).
  use param,only : dp,lmax
  use space,only : nieq
  use enrgy,only : ne
  implicit none
  integer :: ir,l,i,npi
  real(dp) :: psu(ne,0:lmax,nieq),psd(ne,0:lmax,nieq),d
  real(dp),parameter :: pi=acos(-1.d0),pih=pi*0.5d0
  do ir=1,nieq
    do l=0,lmax
      if(l==0)then
        if(psu(1,0,ir)>pih) psu(1,0,ir)=psu(1,0,ir)-pi
      endif
      do i=2,ne
        npi=0
        d=psu(i,l,ir)-psu(i-1,l,ir)
        if(abs(d)>epsilon(1.d0)) npi=nint(d/pi)
        psu(i,l,ir)=psu(i,l,ir)-npi*pi
        !
        npi=0
        d=psd(i,l,ir)-psd(i-1,l,ir)
        if(abs(d)>epsilon(1.d0)) npi=nint(d/pi)
        psd(i,l,ir)=psd(i,l,ir)-npi*pi
      enddo
      if(l>0)then
        !both PS curves running close together.
        if(abs(psu(1,l,ir)-psd(1,l,ir))>pih)then
          if(psu(1,l,ir)<=0.d0) psd(1:ne,l,ir)=psd(1:ne,l,ir)-pi
          if(psu(1,l,ir)> 0.d0) psd(1:ne,l,ir)=psd(1:ne,l,ir)+pi
        endif
      endif
    enddo
  enddo
  return
  end subroutine PSnojump
  !--------------------------------------------------------------------    
  subroutine errVxc0(X,fitness)
  use param,only : dp
  use enrgy,only : ne,eev,v0ev
  implicit none
  integer  :: i
  real(dp) :: X(4),fitness,errv
  errv=0.d0
  do i=1,ne
    errv=errv+(v0ev(i)-max( x(1)+x(2)/sqrt(eev(i)+x(3)), X(4) ))**2
  enddo
  fitness=sqrt(errv/ne)
  return
  end subroutine errVxc0
  !----------------------------------------------------------------------
  subroutine Vxc0EincAprx(Vxc0,v0coef)
  !V0Einca is V0=a+b/sqrt(E+c) for E above minimum V0(E).
  use param,only : dp,outdir
  use enrgy,only : ne,eev
  use DifferentialEvolution,only : DEM,NP,method,itermax,strategy,refresh, &
                                   F_XC,F_CR,CR_XC,VTR,bestval,nfeval,itval
  implicit none
  integer :: i
  integer,parameter :: Dim_XC=4
  real(dp) :: XCmin(Dim_XC),XCmax(Dim_XC),bestmem_XC(Dim_XC)
  real(dp) :: Vxc0(ne),v0coef(4),v0tmp(ne),Vxc1
  real(dp),parameter :: rydb=13.60569172d0
  external errVxc0
  !
  NP=10*Dim_XC
  Vxc1 = minval(Vxc0)*rydb
  XCmax(1) = 0.5d0
  XCmin(1) = Vxc1
    XCmax(2) =  -2.d0
    XCmin(2) = -80.d0
      XCmax(3) = 20.d0
      XCmin(3) =  2.d0
        XCmax(4) = Vxc1+2.0d0
        XCmin(4) = Vxc1
  F_XC=0.5d0;  CR_XC=0.8d0
  method=(/0,1,0/)
  strategy=2; F_CR=0.8d0
  itermax=4000; refresh=-1
  VTR=0.0001d0
  call DEM(errVxc0,XCmin,XCmax,bestmem_XC,Dim_XC)
  if(bestval>=1.d+33)then
    write(61,'(a)')'errVxc0 failed, fitness too large, stop.'
    stop
  endif
  v0coef(1:4)=bestmem_XC(1:4)
  !logX of Vxc0EincAprx.
  write(61,910)'#calls,#improvements =',nfeval,itval
  write(61,921)'v0max  =',XCmax(1:4)
  write(61,921)'v0coef =',v0coef(1:4)
  write(61,921)'v0min  =',XCmin(1:4)
  write(61,'(a,f0.2)')'Vxc0(1)-Vxc0(ne) (eV) = ',(Vxc0(1)-Vxc0(ne))*rydb
  !Vxc0EincAprx plot.
  open(64,file=trim(outdir)//'/Vxc0EincAprx_v0coef',status='unknown')
  write(64,921)'v0coef =',v0coef(1:4)
  do i=1,ne
    v0tmp(i)=max( v0coef(1)+v0coef(2)/sqrt(eev(i)+v0coef(3)), v0coef(4) )
  enddo
  open(64,file=trim(outdir)//'/Vxc0EincAprx',status='unknown')
  write(64,930) (eev(i),v0tmp(i),i=1,ne)
  close(64)
  910 format(a,2(1x,i0),a,es9.2,1x,es9.2)
  921 format(a,4f9.4)
  930 format(2es14.6)
  return
  end subroutine Vxc0EincAprx
  !----------------------------------------------------------------------
  subroutine PStab(tag,delta,v0coef)
  use param,only : dp,outdir,SpinPS,lmax,idZA
  use space,only : nieq
  use enrgy,only : ne,eev, emv0
  implicit none
  character :: tag*2,pid*80,xid*80
  integer :: i,ir,l
  real(dp) :: delta(ne,0:lmax,nieq),v0coef(8)
  !
  !phase shift table PS.
  do ir=1,nieq
    if(tag=='sl')then
      pid='PS.'//trim(idZA(ir))//'.txt'
      xid='px.'//idZA(ir)
      open(40,file=trim(outdir)//'/'//trim(pid),status='replace')
!!!   open(41,file=trim(outdir)//'/'//trim(xid),status='replace')
      write(40,'(i0,1x,4(1x,f0.2),2x,a)') lmax+1,v0coef(1:4)
      do i=1,ne
!!!     Alex: below would have output incident E - muffin tin potential
!!!     changed to output incident E (changed emv0(i) to eev(i))
        write(40,910) eev(i),(delta(i,l,ir),l=0,lmax) !!! Line changed for ViPErLEED
!!!     Line above was before:   write(40,910) emv0(i),(delta(i,l,ir),l=0,lmax)
!!!     write(41,910) emv0(i),(delta(i,l,ir),l=0,min(lmax,9))
      enddo
      close(40)  !!!; close(41)
    elseif(tag=='su'.and.SpinPS=='y')then
      pid='PS.su.'//idZA(ir)
      open(40,file=trim(outdir)//'/'//trim(pid),status='replace')
      do i=1,ne
        write(40,910) emv0(i),(delta(i,l,ir),l=0,lmax)
      enddo
      close(40)
    elseif(tag=='sd'.and.SpinPS=='y')then
      pid='PS.sd.'//idZA(ir)
      open(40,file=trim(outdir)//'/'//trim(pid),status='replace')
      do i=1,ne
        write(40,910) emv0(i),(delta(i,l,ir),l=0,lmax)
      enddo
      close(40)
    endif
  enddo
  return
  910 format(f9.5,20(f9.5:))
  end subroutine PStab
  !----------------------------------------------------------------------
