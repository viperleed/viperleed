program intpol_test
    
    use, intrinsic :: iso_fortran_env

    use interpolation!
    implicit none
    

    real(dp) :: x_data(6), data_out(30), y_data(6)
    real(dp), ALLOCATABLE :: knots(:), work(:)
    integer :: n, n_knots, deg, m, ell

    integer :: kl, ku, ab_cols, ab_rows, offset
    real(dp), ALLOCATABLE :: ab(:,:)

    ! for interpolate test
    integer ierr
    real(dp), ALLOCATABLE :: c(:)


    x_data = (/0,1,2,3,4,5/)
    n = size(x_data)
    y_data = (/0,1,2,3,-3,-1/)
    
    print*, "x_data: ", x_data
    print*, "y_data: ", y_data
    
    deg = 3
    ALLOCATE(knots(get_n_knots(deg,n)))
    call get_natural_knots(x_data, n, deg, knots, n_knots)


    kl = deg
    ku = deg
    ab_rows = 2*kl + ku +1
    ab_cols =n_knots - deg -1
    offset = 1
    allocate(AB(ab_rows, ab_cols))
    call build_colloc_matrix(x_data, n, knots, n_knots, deg, AB, ab_rows, ab_cols, offset)

    print*, "Test interpolation"
    call interpolate_knots(x_data, y_data, n, knots, n_knots, deg, 0, c, ierr)
    print*, "Coeffs: ", c
end program intpol_test
