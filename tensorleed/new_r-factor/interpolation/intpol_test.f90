! A program to test out my new interpolation algorithm and compare it to other options

program intpol_test

    USE tenserleed_intpol
    !use bspline_module
    !use bspline_kinds_module, only: wp, ip
    use, intrinsic :: iso_fortran_env
    use lapack
    use utils
    use splines, only : spline3
    use interpolation
    implicit none
    integer, parameter :: fdp = REAL64
    ! Test data
    integer :: n_points, n_points_supersampled
    real(fdp), external :: test_data
    real(fdp), ALLOCATABLE :: x_data(:), y_data(:), new_y(:), new_y_comp(:), x_ss(:), y_ss(:), my_knots(:), y_new_deriv(:)
    !real(wp), ALLOCATABLE :: x(:), y(:), nx(:), ny(:)
    integer :: deg, i, j, new_n, ierr, n_my_knots
    real(fdp), ALLOCATABLE :: LHS(:,:), RHS(:), RHS_prep(:)
    integer :: nt, kl, ku, lhs_rows, lhs_cols, RHS_cols, nleft, nright
    integer, ALLOCATABLE :: ipiv(:)
    type(grid_pre_evaluation) :: pre_eval

    integer, ALLOCATABLE :: intervals(:)
    real(fdp), ALLOCATABLE :: deBoor_matrix(:,:), coeffs(:), deriv_deBoor_matrix(:,:)
    !integer :: info_values(info_size)

    ! For Bsplines module
    !integer(ip) :: iflag, np, k, one, zero
    !real(wp), ALLOCATABLE :: bcoef(:), knots(:), work(:)
    logical :: extrap

    ! for timing the calls
    real(fdp) :: start, finish, exec_time
    integer :: repeats

    ! for TensErLEED
    integer :: NGP, NB, ipr, ib, NE
    real(fdp), ALLOCATABLE :: A(:), E(:), X_TL(:), WORYT(:)
    real(fdp) :: EINCR

    ! Setup Data to interpolate
    repeats = 10000
    n_points = 600
    n_points_supersampled = n_points*5

    print*, "Interpolation Benchmark"
    print*, "N repeats :             ", repeats
    print*, "N points :              ", n_points
    print*, "N points supersampled : ", n_points_supersampled

    new_n = n_points_supersampled
    Allocate(x_data(n_points))
    ALLOCATE(y_data(n_points))
    ALLOCATE(x_ss(n_points_supersampled))
    ALLOCATE(new_y(n_points_supersampled), y_ss(n_points_supersampled), new_y_comp(n_points_supersampled))

    do i =1, n_points
        x_data(i) = 6d0/n_points*(i-1)! interval from -3 to 3
        y_data(i) = test_data(x_data(i))
    end do
    do i = 1, n_points_supersampled
        x_ss(i) = (6d0-0.1d0)/n_points_supersampled*(i-1) ! interval from -3 to 3
        y_ss(i) = test_data(x_ss(i)) 
    end do


    OPEN(10, file = "x_data.csv", STATUS='REPLACE')
    OPEN(20, file = "y_data.csv",STATUS='REPLACE')
    OPEN(30, file = "y_ss.csv",STATUS='REPLACE') !supersampled (real) y values
    OPEN(40, file = "x_ss.csv",STATUS='REPLACE') !supersampled x_values
    do i= 1,n_points
        write(10, *) x_data(i)
        write(20, *) y_data(i)
    end do
    do i = 1, n_points_supersampled
        write(30, *) y_ss(i)
        write(40, *) x_ss(i)
    end do
    CLOSE(10)
    CLOSE(20)
    Close(30)
    CLOSE(40)

    ! Test out interpolation as used in TensErLEED Refcalc
    ngp = n_points_supersampled+200 ! array size; bigger than n_points and n_points_supersampled
    NE = n_points
    EINCR = x_ss(2)-x_ss(1) ! Delta E value
    ib = 1 !beam number, dummy
    ipr = 0 ! write flag, dummy (!= 2)
    ALLOCATE(A(NGP), E(NGP), WORYT(NGP), X_TL(NGP))
    A = 0.0
    E = 0.0
    A(1:NE) = y_data
    E(1:NE) = x_data
    X_TL = 0.0
    WORYT = 0.0

    call cpu_time(start)
    do i=1, repeats
        CALL INTPOL(A,NE,E,NGP,EINCR,IPR,X_TL,WORYT,IB)
    end do
    call cpu_time(finish)
    exec_time = (finish - start)/repeats
    print*, "TensErLEED Interpolation:         ", exec_time, "s."!," For results see file: ", "y_out_TLEED.csv" 
    OPEN(99, file = "y_out_TLEED.csv", STATUS='REPLACE')
    do i=2, n_points_supersampled
        WRITE(99,*) A(i-1)
    end do
    CLOSE(99)

    ! Test out interpolation via my own module
    
    deg = 5
    new_y_comp = 1.0d0
    call cpu_time(start)
    ALLOCATE(intervals(n_points_supersampled))
    ALLOCATE(y_new_deriv(n_points_supersampled))

    call pre_evaluate_grid(x_data, n_points, x_ss, n_points_supersampled, deg, 1, pre_eval, ierr)
    !call pre_evaluate_deriv(2, n_points_supersampled, x_ss, info_values, my_knots, intervals, deriv_deBoor_matrix, ierr)
    do i=1, repeats
        call interpolate_fast(n_points, y_data, pre_eval, new_y_comp, coeffs, ierr)
        !call interpolate_deriv_fast(info_values, n_points_supersampled, intervals, coeffs, deriv_deBoor_matrix, y_new_deriv, ierr)
    end do
    call cpu_time(finish)
    exec_time = (finish - start)/repeats
    if (ierr.ne.0) then
        print*, "Error thrown", ierr
    end if
    print*, "My Interpolation 5th Order (new): ", exec_time, "s."!, " For results see file: ", "y_out_my_5.csv" 
    OPEN(99, file = "y_out_my_5.csv", STATUS='REPLACE')
    OPEN(98, file = "y_out_my_5_deriv.csv", STATUS='REPLACE')
    do i=1, n_points_supersampled
        WRITE(99,*) new_y_comp(i)
        WRITE(98,*) y_new_deriv(i)
    end do
    CLOSE(99)
    CLOSE(98)
    !DEALLOCATE(my_knots)
    !DEALLOCATE(intervals)
    !DEALLOCATE(LHS, RHS)
    !DEALLOCATE(coeffs)
    !DEALLOCATE(y_new_deriv)


    ! ! Test out interpolation with Fortran-Bsplines library
    ! deg = 3
    ! new_y = 1.0d0
    ! x = x_data
    ! y = y_data
    ! nx = x_ss
    
    ! k = deg+1
    ! np = n_points
    ! one = 1
    ! zero = 0
    ! extrap = .False.  
    ! call cpu_time(start)
    ! ALLOCATE(knots(np+k))
    ! ALLOCATE(bcoef(np))
    ! ALLOCATE(ny(n_points_supersampled))
    ! ALLOCATE(work(3*k))
    ! !CALL db1ink(x, k, np, y, 0, knots, bcoef, iflag
    ! do i = 1, repeats
    !     call db1ink(x,np,y,k,zero,knots,bcoef,iflag)
    !     do j=1,n_points_supersampled
    !         call db1val(nx(j), zero, knots, np, k, bcoef, ny(j), iflag, one, work, extrap)
    !     end do
    ! end do
    ! call cpu_time(finish)
    ! exec_time = (finish - start)/repeats
    ! print*, "F-Bspline Package 3rd Order:      ", exec_time, "s."!," For results see file: ", "y_out_bspline_pckg3.csv" 
    ! OPEN(99, file = "y_out_bsplines_pckg3.csv", STATUS='REPLACE')
    ! do i=1, n_points_supersampled
    !     WRITE(99,*) ny(i)
    ! end do
    ! CLOSE(99)
    ! DEALLOCATE(knots)
    ! DEALLOCATE(bcoef)
    ! DEALLOCATE(ny)
    ! DEALLOCATE(work)

    ! Test out interpolation with Fortran-Bsplines library
    ! deg = 5
    ! new_y = 1.0d0
    ! x = x_data
    ! y = y_data
    ! nx = x_ss
    
    ! k = deg+1
    ! np = n_points
    ! one = 1
    ! zero = 0
    ! extrap = .False.  
    ! call cpu_time(start)
    ! ALLOCATE(knots(np+k))
    ! ALLOCATE(bcoef(np))
    ! ALLOCATE(ny(n_points_supersampled))
    ! ALLOCATE(work(3*k))
    ! !CALL db1ink(x, k, np, y, 0, knots, bcoef, iflag
    ! do i = 1, repeats
    !     call db1ink(x,np,y,k,zero,knots,bcoef,iflag)
    !     do j=1,n_points_supersampled
    !         call db1val(nx(j), zero, knots, np, k, bcoef, ny(j), iflag, one, work, extrap)
    !     end do
    ! end do
    ! call cpu_time(finish)
    ! exec_time = (finish - start)/repeats
    ! print*, "F-Bspline Package 5th Order:      ", exec_time, "s." !," For results see file: ", "y_out_bspline_pckg3.csv" 
    ! OPEN(99, file = "y_out_bsplines_pckg5.csv", STATUS='REPLACE')
    ! do i=1, n_points_supersampled
    !     WRITE(99,*) ny(i)
    
    ! end do
    ! CLOSE(99)

!   !  Test out interpolation with Fortran-Utils Module
!    deg = 3
!    new_y = 1.0d0
!    call cpu_time(start)
!    do i=1, repeats
!        new_y = spline3(x_data,y_data,x_ss)
!    end do
!    call cpu_time(finish)
!    exec_time = (finish - start)/repeats
!    print*, "Cubic Spline from fortran-utils:  ", exec_time, "s."!," For results see file: ", "y_out_futils.csv" 
!    OPEN(99, file = "y_out_futils.csv", STATUS='REPLACE')
!    do i=1, n_points_supersampled
!        WRITE(99,*) new_y(i)
!    end do
!    CLOSE(99)

end program intpol_test


function test_data(x_in) result(y_out)
    use iso_fortran_env
    integer, PARAMETER :: dp = real64
    real(dp), INTENT(IN) :: x_in
    real(dp) :: y_out
    real(dp) :: a
    a = 5

    y_out = 2d0 + cos(a*(x_in-3.0d0))*exp(-(x_in-3.0d0)**2)
    !y_out = 2d0 + cos(x_in*8) + 3*exp(-abs(x_in-2.0d0)**2) + 0.2d0*sin(1.5*x_in) +min(tan(x_in)**2,2.0d0)
    return
end function test_data

