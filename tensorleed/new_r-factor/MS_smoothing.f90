module smoothing

real, parameter :: pi = 3.14159265358979323846
! Coefficients for the MS1 filters, for obtaining a flat passband.
! The innermost arrays contain a, b, c for the fit
! kappa = a + b/(c - m) */

!data for 6th degree coefficient for flat passband
real(8), parameter :: correction_data_MS_6(3,1) =         (/0.001717576, 0.02437382, 1.64375/)
! data for 8th degree coefficient for flat passband
real(8), parameter :: correction_data_MS_8(3,2) = reshape((/0.0043993373, 0.088211164, 2.359375,&
                                                             0.006146815, 0.024715371, 3.6359375/),&
                                                            (/3,2/))
! data for 10th degree coefficient for flat passband
real(8), parameter :: correction_data_MS_10(3,2) = reshape((/0.0011840032, 0.04219344, 2.746875,&
                                                             0.0036718843, 0.12780383, 2.7703125/),&
                                                            (/3,2/))

!data for 4th degree coefficient for flat passband
real(8), parameter :: correction_data_MS1_4(3) =          (/0.021944195, 0.050284006, 0.765625/)
! data for 6th degree coefficient for flat passband
real(8), parameter :: correction_data_MS1_6(3,2) = reshape((/0.0018977303, 0.008476806, 1.2625,&
                                                             0.023064667, 0.13047926, 1.2265625/),&
                                                            (/3,2/))
! data for 8th degree coefficient for flat passband
real(8), parameter :: correction_data_MS1_10(3,2) = reshape((/0.0065903002, 0.057929456, 1.915625,&
                                                             0.0023234477, 0.010298849, 2.2726562,&
                                                             0.021046653, 0.16646601, 1.98125/),&
                                                            (/3,3/))

! data for 10th degree coefficient for flat passband
real(8), parameter :: correction_data_MS1_10(3,2) = reshape((/9.749618E-4, 0.0020742896, 3.74375,&
                                                             0.008975366, 0.09902466, 2.7078125,&
                                                             0.0024195414, 0.010064855, 3.296875,&
                                                             0.019185117, 0.18953617, 2.784961/),&
                                                            (/3,4/))

contains

subroutine MS_smoother(degree)
    integer degree

    ! Internal
    integer max_degree, m_min

    max_degree = 10
    if ((degree<2).or.(degree>max_degree).or.(modulo(degree,2)==1)) then
        !invalid input
    end if
    m_min = int(degree/2)+2
    if (m < m_min) then
        ! problem, invalid input
    end if
    kernel = make_kernel_MS(degree, m)
    fit_weigths = make_fit_weights_MS(degree, m)
    !public ModifiedSincSmoother(boolean isMS1, int degree, int m) {
    !    this.isMS1 = isMS1;
    !    this.degree = degree;
    !    if (degree < 2 || degree > MAX_DEGREE || (degree&0x1)!=0)
    !        throw new IllegalArgumentException("Invalid degree "+degree+"; only 2, 4, ... "+MAX_DEGREE+" supported");
    !    int mMin = isMS1 ? degree/2+1 : degree/2+2;
    !    if (m < mMin)  //kernel not wide enough for the wiggles of the sinc function
    !        throw new IllegalArgumentException("Invalid kernel half-width "+m+"; must be >= "+mMin);
    !    kernel = makeKernel(isMS1, degree, m);
    !    fitWeights = makeFitWeights(isMS1, degree, m);
    !}

end subroutine MS_smoother

pure function make_kernel_MS(degree, m, coeffs, n_coeffs) result(kernel)
    integer degree, m
    real coeffs(ncoeffs)

    real(8) kernel(m+1)
    real(8) sum, x, sinc_arg, k
    integer nu

    sum = 0d0
    decay = 4 !MS: decay alpha =2: 13.5% at end without correction, 2sqrt2 sigma

    do i=0, m
        x = i*(1d0/(m+1))
        sinc_arg = pi*0.5d0*(degree+4)*x
        if (i==0) then
            k = 1 ! lim x->0 sin(x)/x = 1 (saved from 0/0 by l'Hopitals' rule)
        else
            k = sin(sinc_arg)/sinc_arg ! sinc = sin(x)/x
        end if
        do j=0, ncoeffs
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

pure function make_kernel_MS1(degree, m, coeffs, n_coeffs) result(kernel)
    integer degree, m
    real coeffs(ncoeffs)

    real(8) kernel(m+1)
    real(8) sum, x, sinc_arg, k
    integer nu


    sum = 0d0
    decay = 2 !MS: decay alpha =2: 13.5% at end without correction, 2sqrt2 sigma

    if (modulo((degree/2),2)==0) then
        nu = 2
    else
        nu = 1
    end if

    do i=0, m
        x = i*(1d0/(m+1))
        sinc_arg = pi*0.5d0*(degree+2)*x
        if (i==0) then
            k = 1 ! lim x->0 sin(x)/x = 1 (saved from 0/0 by l'Hopitals' rule)
        else
            k = sin(sinc_arg)/sinc_arg ! sinc = sin(x)/x
        end if

        do j=0, ncoeffs
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
    integer, intent(in) :: degree, m

    real(8), intent(out), allocatable :: coeffs(:)

    ! Internal
    real(8), allocatable :: correctionData(:,:)
    real(8), allocatable ::corrForDeg(:,:)
    integer n_coeffs, n_corr_lines
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

    do i = 0, n_corr_lines
        abc = corrForDeg(i)
        cm = abc(2) - m
        coeffs(i) = abc(0) + abc(1)/(cm**3)
    end do

end function get_coefficients_MS

pure function get_coefficients_MS1(degree, m) result(coeffs)
    integer, intent(in) :: degree, m

    real(8), intent(out), allocatable :: coeffs(:)

    ! Internal
    real(8), allocatable :: correctionData(:,:)
    real(8), allocatable ::corrForDeg(:,:)
    integer n_coeffs, n_corr_lines
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

    do i = 0, n_corr_lines
        abc = corrForDeg(i)
        cm = abc(2) - m
        coeffs(i) = abc(0) + abc(1)/(cm**3)
    end do

end function get_coefficients_MS

pure function make_fit_weights_MS(degree, m) result(weights)
    integer, intent(in) :: degree, m

    real(8), allocatable :: weights(:)

    real first_zero, beta
    integer fitlength

    first_zero = (m+1)/(1+0.5d0*degree)
    beta = 0.65d0 + 0.35d0*exp(-0.55d0*(degree-4))
    fitlength = int(ceiling(first_zero*beta))
    allocate(weights(fitlength))
    do p = 0, fitlength
        weights(p) = sqrt(cos(0.5d0*pi/(first_zero*beta)*p))
    end do
end function make_fit_weights_MS

pure function make_fit_weights_MS1(degree, m) result(weights)
   integer, intent(in) :: degree, m

    real(8), allocatable :: weights(:)

    real first_zero, beta
    integer fitlength

    first_zero = (m+1)/(1.5d0+0.5d0*degree)
    beta = 0.7d0 + 0.14d0*exp(-0.6d0*(degree-4))
    fitlength = int(ceiling(first_zero*beta))
    allocate(weights(fitlength))
    do p = 0, fitlength
        weights(p) = sqrt(cos(0.5d0*pi/(first_zero*beta)*p))
    end do
end function make_fit_weights_MS1

pure function extend_data(data, fit_weights m, degree) result(extended)
    integer, intent(in) :: m, degree
    real(8), intent(in) :: data(:)

    real(8) extended_data(size(data+2*m))

    integer fit_length
    real(8), allocatable :: lin_reg_x(:), lin_reg_y(:), lin_reg_weights(:)


    fit_length = int(min(size(fit_weights), size(data)))
    allocate(lin_reg_x(fit_length), lin_reg_y(fit_length), lin_reg_weights(fit_length))

    ! Unfinished; continue here tomorrow
    do p = 0, fit_length
        lin_reg_x(p) = p
        lin_reg_y(p) = data(size(data)-1-p)
        lin_reg_weights(p) = fit_weights(p)
        ! Linear regression needs to be implemented
    end do

end function extend_data

pure function linear_regression_weighted(x, y, weights) result(d_k)
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
    d_k(1) = (sum_y-d_k(2)*sum_x)/sum_weighted

end function linear_regression

end module smoothing
