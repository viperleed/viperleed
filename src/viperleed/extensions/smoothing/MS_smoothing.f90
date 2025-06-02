module smoothing
implicit none

! f2py -m smoothing MS_smoothing.f90 -h smoothing.pyf --overwrite-signature --debug-capi
! f2py -c smoothing.pyf MS_smoothing.f90

real, parameter :: pi = 3.14159265358979323846
integer, parameter :: max_coeffs_rows = 4
integer, parameter :: max_degree = 10
! Coefficients for the MS1 filters, for obtaining a flat passband.
! The innermost arrays contain a, b, c for the fit
! kappa = a + b/(c - m) */

!data for 6th degree coefficient for flat passband
real(8), parameter :: correction_data_MS_6(3,1) = reshape((/0.001717576, 0.02437382, 1.64375/), (/3,1/))
! data for 8th degree coefficient for flat passband
real(8), parameter :: correction_data_MS_8(3,2) = reshape((/0.0043993373, 0.088211164, 2.359375,&
                                                             0.006146815, 0.024715371, 3.6359375/),&
                                                            (/3,2/))
! data for 10th degree coefficient for flat passband
real(8), parameter :: correction_data_MS_10(3,2) = reshape((/0.0011840032, 0.04219344, 2.746875,&
                                                             0.0036718843, 0.12780383, 2.7703125/),&
                                                            (/3,2/))

!data for 4th degree coefficient for flat passband
real(8), parameter :: correction_data_MS1_4(3,1) = reshape((/0.021944195, 0.050284006, 0.765625/), (/3,1/))
! data for 6th degree coefficient for flat passband
real(8), parameter :: correction_data_MS1_6(3,2) = reshape((/0.0018977303, 0.008476806, 1.2625,&
                                                             0.023064667, 0.13047926, 1.2265625/),&
                                                            (/3,2/))
! data for 8th degree coefficient for flat passband
real(8), parameter :: correction_data_MS1_8(3,3) = reshape((/0.0065903002, 0.057929456, 1.915625,&
                                                             0.0023234477, 0.010298849, 2.2726562,&
                                                             0.021046653, 0.16646601, 1.98125/),&
                                                            (/3,3/))

! data for 10th degree coefficient for flat passband
real(8), parameter :: correction_data_MS1_10(3,4) = reshape((/9.749618E-4, 0.0020742896, 3.74375,&
                                                             0.008975366, 0.09902466, 2.7078125,&
                                                             0.0024195414, 0.010064855, 3.296875,&
                                                             0.019185117, 0.18953617, 2.784961/),&
                                                            (/3,4/))

contains

subroutine MS_smoother(data, isMS1, degree, m, out_data)
    implicit none
    ! Input
    integer, intent(in) :: degree, m, isMS1
    real(8), intent(in) :: data(:)

    ! Output
    real(8), intent(out) :: out_data(size(data))

    ! Internal
    integer max_degree, m_min, radius, ncoeffs, fitlength, p
    real(8) :: coeffs(3*4) !TODO: use paramter max_coeffs_rows here
    real(8), allocatable :: extended_data(:), extended_smoothed(:), fit_weights(:), kernel(:)
    real(8) :: first_zero, beta

    ! Get coefficients from parameters above
    call get_coefficients(isMS1, degree, m, ncoeffs, coeffs)

    if (isMS1 == 1) then
        first_zero = (m+1)/(1+0.5d0*degree)
        beta = 0.65d0 + 0.35d0*exp(-0.55d0*(degree-4))
    else
        first_zero = (m+1)/(1.5d0+0.5d0*degree)
        beta = 0.7d0 + 0.14d0*exp(-0.6d0*(degree-4))
        fitlength = int(ceiling(first_zero*beta))
    end if

    fitlength = int(ceiling(first_zero*beta))
    allocate(fit_weights(fitlength))
    do p = 1, fitlength+1
        fit_weights(p) = (cos(0.5d0*pi/(first_zero*beta)*(p-1)))**2
    end do

    kernel = make_kernel(isMS1, degree, m, coeffs, ncoeffs) ! should be allocated automaticaly
    radius = size(kernel) - 1
    extended_data = extend_data(data, fit_weights, m)
    extended_smoothed = smooth_except_boundaries(extended_data, kernel)
    out_data = extended_smoothed(radius+1:radius+1+size(data))

    ! Below can be used for debugging; uncomment to print out arrays/values
    !write(*,*) "ncoeffs: ", ncoeffs
    !write(*,*) "coeffs: ", coeffs
    !write(*,*) "kernel: ", kernel
    !write(*,*) "fitlength: ", fitlength
    !write(*,*) "fit_weights: ", fit_weights
    !write(*,*) "extended_data: ", extended_data
    !write(*,*) "smooth_except_boundaies: ", extended_smoothed
    !write(*,*) "out: ", out_data
    return
end subroutine MS_smoother

pure function make_kernel(isMS1, degree, m, coeffs, ncoeffs) result(kernel)
    implicit none
    integer, intent(in) :: isMS1, degree, m, ncoeffs
    real(8), intent(in) :: coeffs(ncoeffs)

    real(8) kernel(m+1)
    real(8) sum, x, sinc_arg, k, decay
    integer nu, i, j

    sum = 0d0
    if (isMS1 == 1) then
        decay = 2 !MS: decay alpha =2: 13.5% at end without correction, 2sqrt2 sigma
        if (modulo((degree/2),2)==0) then
            nu = 2
        else
            nu = 1
        end if
    else
        decay = 4 !MS: decay alpha =2: 13.5% at end without correction, 2sqrt2 sigma
        ! nu not needed
    end if

    do i=0, m
            x = real(i)*(1d0/(m+1))
            if (isMS1 == 1) then
                sinc_arg = pi*0.5d0*(degree+2)*x ! depends on isMS1
            else
                sinc_arg = pi*0.5d0*(degree+4)*x ! depends on isMS1
            end if

            if (i==0) then
                k = 1 ! lim x->0 sin(x)/x = 1 (saved from 0/0 by l'Hopitals' rule)
            else
                k = sin(sinc_arg)/sinc_arg ! sinc = sin(x)/x
            end if
            if (isMS1 == 1) then
                do j=0, ncoeffs-1
                    k = k + coeffs(j+1)*x*sin((2*j+nu)*pi*x)
                end do
            else
                do j=0, ncoeffs-1
                    k = k + coeffs(j+1)*x*sin((j+1)*pi*x)
                end do
            end if
            k = k*(exp(-x**2*decay) + exp(-(x-2)**2*decay) + exp(-(x+2)**2*decay) - 2*exp(-decay) - exp(-9*decay))
            kernel(i+1) = k
            if (i>0) then
                sum = sum + 2*k !MS: off-center kernel elements appear twice
            else
                sum = sum + k
            end if
    end do


    !Finally normalize kernel elements
    kernel = kernel*(1/sum) ! Whole array operation
    return
end function make_kernel

subroutine get_coefficients(isMS1, degree, m, ncoeffs, coeffs)
    implicit none
    integer, intent(in) :: isMS1, degree, m

    integer, intent(out) :: ncoeffs
    real(8), intent(out) :: coeffs(3*4) !TODO: use paramter max_coeffs_rows here

    ! Internal
    real(8), allocatable :: corrForDeg(:,:)
    integer n_corr_lines, i
    real(8) :: abc(3), cm

    coeffs = 0
    if (isMS1 == 1) then ! MS1
        if (degree == 0) then
            n_corr_lines = 0
        elseif (degree == 2) then
            n_corr_lines = 0
        elseif (degree == 4) then
            n_corr_lines = 1
            allocate(corrForDeg(3,n_corr_lines))
            corrForDeg = correction_data_MS1_4
        elseif (degree == 6) then
            n_corr_lines = 2
            allocate(corrForDeg(3,n_corr_lines))
            corrForDeg = correction_data_MS1_6
        elseif (degree == 8) then
            n_corr_lines = 3
            allocate(corrForDeg(3,n_corr_lines))
            corrForDeg = correction_data_MS1_8
        elseif (degree == 10) then
            n_corr_lines = 4
            allocate(corrForDeg(3,n_corr_lines))
            corrForDeg = correction_data_MS1_10
        end if
    else ! MS
        if (degree == 0) then
            n_corr_lines = 0
        elseif (degree == 2) then
            n_corr_lines = 0
        elseif (degree == 4) then
            n_corr_lines = 0
        elseif (degree == 6) then
            n_corr_lines = 1
            allocate(corrForDeg(3,n_corr_lines))
            corrForDeg = correction_data_MS_6
        elseif (degree == 8) then
            n_corr_lines = 2
            allocate(corrForDeg(3,n_corr_lines))
            corrForDeg = correction_data_MS_8
        elseif (degree == 10) then
            n_corr_lines = 2
            allocate(corrForDeg(3,n_corr_lines))
            corrForDeg = correction_data_MS_10
        end if
    end if

    ncoeffs = n_corr_lines*3

    do i = 0, n_corr_lines-1
        abc = corrForDeg(:,i+1)
        cm = abc(3) - m
        coeffs(i+1) = abc(1) + abc(2)/(cm**3)
    end do


end subroutine get_coefficients

pure function extend_data(data, fit_weights, m) result(extended)
    implicit none
    integer, intent(in) :: m
    real(8), intent(in) :: data(:), fit_weights(:)

    !Output
    real(8), allocatable :: extended(:)

    integer fit_length, p
    real(8), allocatable :: lin_reg_x(:), lin_reg_y(:), lin_reg_weights(:)
    real(8) :: d_k(2)


    fit_length = int(min(size(fit_weights), size(data)))
    allocate(lin_reg_x(fit_length), lin_reg_y(fit_length), lin_reg_weights(fit_length))
    do concurrent (p = 1: fit_length)
        lin_reg_x(p) = p
        lin_reg_y(p) = data(p)
        lin_reg_weights(p) = fit_weights(p)
    end do
    d_k =linear_regression_weighted(lin_reg_x, lin_reg_y, lin_reg_weights)
    allocate(extended(size(data)+2*m))
    do p = 1, m
        extended(1+m-p) = d_k(1) - d_k(2)*(p-1)
    end do
    do concurrent (p = 1: fit_length)
        lin_reg_x(p) = p
        lin_reg_y(p) = data(size(data)-p+1) ! -p+1 rather than -p due to Fortran indices
        lin_reg_weights(p) = fit_weights(p) ! fit weights are of the form (cos(a*p))**2, which is symmertric thus p can be used as index
    end do
    d_k =linear_regression_weighted(lin_reg_x, lin_reg_y, lin_reg_weights)
    do p =1, m
        extended((size(data)+m+p)) = d_k(1) - d_k(2)*(p-1)
    end do
    extended(m+1: size(data)+m+1) = data
    return
end function extend_data

pure function smooth_except_boundaries(data, kernel) result(smoothed_data)
    implicit none
    real(8), intent(in) :: data(:), kernel(:)

    real(8) :: smoothed_data(size(data))

    integer radius, i, j

    smoothed_data = 0d0 ! initialization
    radius = size(kernel) -1
    do concurrent (i = radius: size(data)-radius-1)
        smoothed_data(i+1) = kernel(1)*data(i+1)
        do j = 1, size(kernel)-1
            smoothed_data(i+1) = smoothed_data(i+1) + kernel(j+1)*(data(i-j+1)+data(i+j+1))
        end do
    end do

end function  smooth_except_boundaries

pure function linear_regression_weighted(x, y, weights) result(d_k)
    implicit none
    real(8), intent(in) :: x(:), y(:), weights(:)

    !Returns
    real(8) d_k(2)

    !Internal
    real(8) :: sum_weights, sum_x, sum_y, sum_xy, sum_x2, sum_y2
    integer n

    n = size(x)
    sum_weights = sum(weights)
    sum_x = sum(x*weights)
    sum_y = sum(y*weights)
    sum_xy = sum(x*y*weights)
    sum_x2 = sum(x**2*weights)
    sum_y2 = sum(y**2*weights)

    ! Slope
    d_k(2) = (sum_xy-sum_x*sum_y*(1/sum_weights))/(sum_x2-sum_x**2*(1/sum_weights))
    ! Intercept
    d_k(1) = (sum_y-d_k(2)*sum_x)/sum_weights

end function linear_regression_weighted

end module smoothing
