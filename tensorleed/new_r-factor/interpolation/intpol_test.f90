program intpol_test
    
    use, intrinsic :: iso_fortran_env

    use interpolation!
    implicit none
    

    real(dp) :: x_data(11), data_out(30), y_data(11)
    real(dp), ALLOCATABLE :: knots(:)
    integer :: n, n_knots, deg
    x_data = (/0,1,2,3,4,5,6,7,8,9,10/)
    n = size(x_data)
    y_data = (/0,1,2,3,-3,-1,0,1,-2,3,0/)
    
    print*, "x_data: ", x_data
    print*, "y_data: ", y_data
    
    deg = 3
    print*, "get_n_knots: ", get_n_knots(deg, n)
    ALLOCATE(knots(get_n_knots(deg,n)))
    call get_natural_knots(x_data, n, deg, knots, n_knots)
    print*, "get_natural_knots:", knots
    print*, "find_interval: ", find_interval(knots, n, deg, 3.4d0, 1, 0)
    
    !print*, "de_Boor ", find_interval(knots, n, deg, 3.4, 1, 0)



end program intpol_test
