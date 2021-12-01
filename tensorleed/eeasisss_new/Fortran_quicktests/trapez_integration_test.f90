! Program to quickly check integration of charge density file chgdenAu
! Author: Alexander M. Imre, October 2021

  program trapez_integration_test
        implicit none
        integer, parameter :: dp = kind(1.0d0)
        integer nx, i
        real(kind=dp) :: z, HFrmin, HFrmax,a
        real(dp),allocatable :: r(:),atrho4pir2(:), test(:), s(:)
        real :: result1, result2, result3
        ! names adapted from eeasisss_init.f90

        ! read file chgdenAu for testing into unit 50
        open(50,file='chgdenAu',status='old')
        read(50,*)
        read(50,*) z,nx,HFrmin,HFrmax ! nx is number of steps

        allocate(r(nx), atrho4pir2(nx), test(nx), s(nx))

        ! read in from file
        do i =1,nx
            read(50,*)r(i),atrho4pir2(i)
        end do
        close(50)


        write(*,*) "Expected charge: ", z
!        call trapez_unequal_dx(r,nx,atrho4pir2,1,nx,result1)
!        write(*,*) "Integrated charge, old: ", result1
        call trapezE(r,nx,1.3553317697208241d-2,atrho4pir2,1,nx,result2)
        write(*,*) "Integrated charge, tapezE: ", result2
        call trapezR(r,nx, atrho4pir2,1,nx,result3)
        write(*,*) "Integrated charge, trapez: ", result3

  end program trapez_integration_test


subroutine trapez_unequal_dx(rx,nx,f,j1,j2,result)
  ! in part based on trapez
  implicit none
  integer :: nx, j1, j2
  integer, parameter :: dp = kind(1.0d0)
  real(kind=dp) :: rx(nx),f(nx),s(nx-1), dx(nx-1), result
  ! initialize values to be safe
  s(:) = 0.0d0
  ! whole array operations
  dx(:) = rx(2:)-rx(:nx-1)
  s(:) = (f(:nx-1)+f(2:))
  result = sum(s*dx)*0.5d0
  return
  end subroutine trapez_unequal_dx


  subroutine trapezA(rx,nx,f,j1,j2,res)
  implicit none
  !trapez integrated on exponential radial grid.
  !Author: Alexander M. Imre, October 2021.
  integer :: nx,j1,j2
  integer, parameter :: dp = kind(1.0d0)
  real(dp) :: rx(nx),f(nx),s(nx-1),dx(nx-1)
  real :: res
  s(:)=0.d0
  dx(:)=(rx(2:)-rx(:nx-1))
  s(:)=(f(:nx-1)+f(2:))
  write(*,*) size(dx), size(s), size((f(:nx-1)+f(2:))), size(rx(2:)-rx(:nx-1))
  write(*,*) sum(s*dx)*0.5d0
  res=sum(s*dx)
  res = res * 0.5d0
  write(*,*) res
  return
  end subroutine trapezA
  !-------------------------------------------------------------------
subroutine trapezR(rx,nx,rho4pir2,h,k,q)
  implicit none
  !trapezoidal integral on exponential grid x = x(1)*exp[(i-1)*dx],
  !where rho4pir2(x)*dx = rho4pir2(x)*(dx/di)*di with dx/di = x(i)*dx.
  !trapez in
  !Abramowitz-Stegun, Handbook of Mathematical Functions, Sec. 25.4.4.
  integer :: nx,h,k
  integer, parameter :: dp = kind(1.0d0)
  real(dp) :: rx(nx),rho4pir2(nx),dx(nx-1)
  real :: q
  dx(:) = (rx(h+1:k)-rx(h:k-1))
  write(*,*) size(dx)
  q=(rho4pir2(h)*rx(h)*0.5d0 + sum(rho4pir2(h+1:k-1)*rx(h+1:k-1)) + &
     rho4pir2(k)*rx(k)*0.5d0)*log(rx(2)/rx(1))
  write(*,*) q/79d0
  write(*,*) "log(x2/x1) = ", log(rx(3)/rx(2))/2
  return
  end subroutine trapezR

  subroutine trapezE(rx,nx,dx,f, j1, j2,q)
    implicit none
    !trapez integrated on exponential radial grid.
    !Author: Alexander M. Imre, October 2021.
    integer :: nx, j1, j2, i
    integer, parameter :: dp = kind(1.0d0)
    real(dp) :: rx(nx),f(nx),s(nx-1),dx, step(nx-1)
    real q
    do i = j1,j2-1
      step(i)=rx(1)/exp(dx)*(exp((i+1)*dx)-exp((i)*dx))
    end do
    s(:)=0.d0
    s(:)=f(:nx-1)+f(2:)
    q=sum(s*step)*0.5d0
    return
    end subroutine trapezE