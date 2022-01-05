program intpol_test
    
    use, intrinsic :: iso_fortran_env

    use interpolation
    implicit none
    
    integer :: n_points, n_points_supersampled
    real(dp), ALLOCATABLE :: x_data(:), y_real(:), y_data(:)
    real(dp), ALLOCATABLE :: knots(:), work(:)
    integer :: n, n_knots, deg, m, ell, i

    integer :: kl, ku, ab_cols, ab_rows, offset
    real(dp), ALLOCATABLE :: ab(:,:)

    ! for interpolate test
    integer ierr
    real(dp), ALLOCATABLE :: c(:)

    n_points = 30
    n_points_supersampled = 300
    Allocate(x_data(n_points))
    ALLOCATE(y_data(n_points))
    ALLOCATE(y_real(n_points_supersampled))

    OPEN(10, file = "x_data.csv")
    OPEN(20, file = "y_data.csv")
    OPEN(30, file = "y_real.csv")
    do i= 1,n_points
        read(10, *) x_data(i)
        read(20, *) y_data(i)
    end do
    do i = 1, n_points_supersampled
        read(30, *) y_real(i)
    end do
    CLOSE(10)
    CLOSE(20)
    Close(30)

    n = n_points

    !print*, "x_data: ", x_data
    !print*, "y_data: ", y_data
    
    deg = 3
    ALLOCATE(knots(get_n_knots(deg,n)))
    call get_natural_knots(x_data, n, deg, knots, n_knots)


    kl = deg
    ku = deg
    ab_rows = 2*kl + ku +1
    ab_cols =n_knots - deg -1
    offset = 1

    print*, "Test interpolation"
    call interpolate_knots(x_data, y_data, n, knots, n_knots, deg, 0, c, ierr)
    print*, "Coeffs: ", c
end program intpol_test
