module smoothing
implicit none

! f2py -m smoothing MS_smoothing.f90 -h smoothing.pyf --overwrite-signature --debug-capi
! f2py -c smoothing.pyf MS_smoothing.f90

real, parameter :: pi = 3.14159265358979323846
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
    integer, intent(in) :: degree, m, isMS1
    real(8), intent(in) :: data(:)

    real(8), intent(out) :: out_data(size(data))

    ! Internal
    integer max_degree, m_min, radius, ncoeffs
    real(8), allocatable :: extended_data(:), extended_smoothed(:), fit_weights(:), kernel(:), coeffs(:)


    max_degree = 10
    if ((degree<2).or.(degree>max_degree).or.(modulo(degree,2)==1)) then
        !invalid input
    end if
    m_min = int(degree/2)+2
    if (m < m_min) then
        ! problem, invalid input
    end if

    if (isMS1 == 1) then
        coeffs = get_coefficients_MS1(degree, m)
        ncoeffs = size(coeffs)
        fit_weights = make_fit_weights_MS1(degree, m) ! should be allocated automaticaly
        kernel = make_kernel_MS1(degree, m, coeffs, ncoeffs) ! should be allocated automaticaly
    else
        coeffs = get_coefficients_MS(degree, m)
        ncoeffs = size(coeffs)
        fit_weights = make_fit_weights_MS(degree, m) ! should be allocated automaticaly
        kernel = make_kernel_MS(degree, m, coeffs, ncoeffs) ! should be allocated automaticaly
    end if

    radius = size(kernel) - 1
    allocate(extended_data(size(data)+2*radius), extended_smoothed(size(data)+2*radius))
    extended_data = extend_data(data, fit_weights, m)
    extended_smoothed = smooth_except_boundaries(data, kernel)
    out_data = extended_smoothed(radius+1:radius+1+size(data))
    return
end subroutine MS_smoother


pure function make_kernel_MS(degree, m, coeffs, ncoeffs) result(kernel)
    implicit none
    integer, intent(in) :: degree, m, ncoeffs
    real(8), intent(in) :: coeffs(ncoeffs)

    real(8) kernel(m+1)
    real(8) sum, x, sinc_arg, k, decay
    integer nu, i, j

    sum = 0d0
    decay = 4 !MS: decay alpha =2: 13.5% at end without correction, 2sqrt2 sigma

    do i=1, m-1
        x = i*(1d0/(m+1))
        sinc_arg = pi*0.5d0*(degree+4)*x
        if (i==0) then
            k = 1 ! lim x->0 sin(x)/x = 1 (saved from 0/0 by l'Hopitals' rule)
        else
            k = sin(sinc_arg)/sinc_arg ! sinc = sin(x)/x
        end if
        do j=1, ncoeffs
            k = k + coeffs(j)*x*sin((j+1)*pi*x)
        end do
        k = k*(exp(-x**2*decay) + exp(-(x-2)**2*decay) + exp(-(x+2)**2*decay) - 2*exp(-decay) - exp(-9*decay))
        kernel(i) = k

        if (i>0) then
            sum = sum + 2*k !MS: off-center kernel elements appear twice
        else
            sum = sum + k
        end if
    end do
    !Finally normalize kernel elements
    kernel = kernel*(1/sum) ! Whole array operation
    return
end function make_kernel_MS

pure function make_kernel_MS1(degree, m, coeffs, ncoeffs) result(kernel)
    implicit none
    integer, intent(in) :: degree, m, ncoeffs
    real(8), intent(in) :: coeffs(ncoeffs)

    real(8) kernel(m+1)
    real(8) sum, x, sinc_arg, k, decay
    integer nu, i, j


    sum = 0d0
    decay = 2 !MS: decay alpha =2: 13.5% at end without correction, 2sqrt2 sigma

    if (modulo((degree/2),2)==0) then
        nu = 2
    else
        nu = 1
    end if

    do i=1, m-1
        x = i*(1d0/(m+1))
        sinc_arg = pi*0.5d0*(degree+2)*x
        if (i==0) then
            k = 1 ! lim x->0 sin(x)/x = 1 (saved from 0/0 by l'Hopitals' rule)
        else
            k = sin(sinc_arg)/sinc_arg ! sinc = sin(x)/x
        end if

        do j=1, ncoeffs
            k = k + coeffs(j)*x*sin((2*j+nu)*pi*x)
        end do

        k = k*(exp(-x**2*decay) + exp(-(x-2)**2*decay) + exp(-(x+2)**2*decay) - 2*exp(-decay) - exp(-9*decay))
        kernel(i) = k
        if (i>0) then
            sum = sum + 2*k !MS: off-center kernel elements appear twice
        else
            sum = sum + k
        end if
    end do
    !Finally normalize kernel elements
    kernel = kernel*(1/sum) ! Whole array operation
    return
end function make_kernel_MS1

pure function get_coefficients_MS(degree, m) result(coeffs)
    implicit none
    integer, intent(in) :: degree, m

    real(8), allocatable :: coeffs(:)

    ! Internal
    real(8), allocatable ::corrForDeg(:,:)
    integer n_coeffs, n_corr_lines, i
    real(8) :: abc(3), cm


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
    n_coeffs = n_corr_lines*3
    allocate(coeffs(n_coeffs))

    do i = 1, n_corr_lines
        abc = corrForDeg(:,i)
        cm = abc(3) - m
        coeffs(i) = abc(1) + abc(2)/(cm**3)
    end do

end function get_coefficients_MS

pure function get_coefficients_MS1(degree, m) result(coeffs)
    implicit none
    integer, intent(in) :: degree, m

    real(8), allocatable :: coeffs(:)

    ! Internal
    real(8), allocatable ::corrForDeg(:,:)
    integer n_coeffs, n_corr_lines, i
    real(8) :: abc(3), cm


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
    n_coeffs = n_corr_lines*3
    allocate(coeffs(n_coeffs))

    do i = 1, n_corr_lines
        abc = corrForDeg(:,i)
        cm = abc(3) - m
        coeffs(i) = abc(1) + abc(2)/(cm**3)
    end do

end function get_coefficients_MS1

pure function make_fit_weights_MS(degree, m) result(weights)
    implicit none
    integer, intent(in) :: degree, m

    real(8), allocatable :: weights(:)

    real first_zero, beta
    integer fitlength, p

    first_zero = (m+1)/(1+0.5d0*degree)
    beta = 0.65d0 + 0.35d0*exp(-0.55d0*(degree-4))
    fitlength = int(ceiling(first_zero*beta))
    allocate(weights(fitlength))
    do p = 0, fitlength
        weights(p) = sqrt(cos(0.5d0*pi/(first_zero*beta)*p))
    end do
end function make_fit_weights_MS

pure function make_fit_weights_MS1(degree, m) result(weights)
    implicit none
   integer, intent(in) :: degree, m

    real(8), allocatable :: weights(:)

    real first_zero, beta
    integer fitlength, p

    first_zero = (m+1)/(1.5d0+0.5d0*degree)
    beta = 0.7d0 + 0.14d0*exp(-0.6d0*(degree-4))
    fitlength = int(ceiling(first_zero*beta))
    allocate(weights(fitlength))
    do p = 0, fitlength
        weights(p) = sqrt(cos(0.5d0*pi/(first_zero*beta)*p))
    end do
end function make_fit_weights_MS1

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

    do concurrent (p = 0: fit_length-1)
        lin_reg_x(p) = p
        lin_reg_y(p) = data(p)
        lin_reg_weights(p) = fit_weights(p)
    end do
    d_k =linear_regression_weighted(lin_reg_x, lin_reg_y, lin_reg_weights)

    do p =1, m
        extended(m-p) = d_k(1) - d_k(2)*p
    end do

    do concurrent (p = 0: fit_length-1)
        lin_reg_x(p) = p
        lin_reg_y(p) = data(size(data)-p) ! -p rather than -(p+1) due to Fortran indices
        lin_reg_weights(p) = fit_weights(p) ! fit weights are of the form sqrt(cos(a*p)), which is symmertric thus p can be used as index
    end do
    d_k =linear_regression_weighted(lin_reg_x, lin_reg_y, lin_reg_weights)

    allocate(extended(size(data+2*m)))
    do p =1, m
        extended((size(data)+m-1+p)) = d_k(1) - d_k(2)*p
    end do
    extended(m: size(data)+m) = data
    return
end function extend_data

pure function smooth_except_boundaries(data, kernel) result(smoothed_data)
    implicit none
    real(8), intent(in) :: data(:), kernel(:)

    real(8) :: smoothed_data(size(data))

    integer radius, i, j


    radius = size(kernel) -1

    do concurrent (i = radius: size(data)-radius)
        smoothed_data(i) = kernel(1)*data(i)
        do j = 1, size(kernel)-1
            smoothed_data(i) = smoothed_data(i) + kernel(j)*(data(i-j)+data(i+j))
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
    d_k(2) = sum_xy-sum_x*sum_y*(1/sum_weights)/(sum_x2*sum_x**2*(1/sum_weights))
    ! Intercept
    d_k(1) = (sum_y-d_k(2)*sum_x)/sum_weights

end function linear_regression_weighted

end module smoothing
