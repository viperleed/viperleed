! Module interpolation defined in file test_deriv_types.f90

subroutine f90wrap_info_values__get__deg(this, f90wrap_deg)
    use interpolation, only: info_values
    implicit none
    type info_values_ptr_type
        type(info_values), pointer :: p => NULL()
    end type info_values_ptr_type
    integer, intent(in)   :: this(2)
    type(info_values_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_deg
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_deg = this_ptr%p%deg
end subroutine f90wrap_info_values__get__deg

subroutine f90wrap_info_values__set__deg(this, f90wrap_deg)
    use interpolation, only: info_values
    implicit none
    type info_values_ptr_type
        type(info_values), pointer :: p => NULL()
    end type info_values_ptr_type
    integer, intent(in)   :: this(2)
    type(info_values_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_deg
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%deg = f90wrap_deg
end subroutine f90wrap_info_values__set__deg

subroutine f90wrap_info_values__get__n_knots(this, f90wrap_n_knots)
    use interpolation, only: info_values
    implicit none
    type info_values_ptr_type
        type(info_values), pointer :: p => NULL()
    end type info_values_ptr_type
    integer, intent(in)   :: this(2)
    type(info_values_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_n_knots
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_n_knots = this_ptr%p%n_knots
end subroutine f90wrap_info_values__get__n_knots

subroutine f90wrap_info_values__set__n_knots(this, f90wrap_n_knots)
    use interpolation, only: info_values
    implicit none
    type info_values_ptr_type
        type(info_values), pointer :: p => NULL()
    end type info_values_ptr_type
    integer, intent(in)   :: this(2)
    type(info_values_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_n_knots
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%n_knots = f90wrap_n_knots
end subroutine f90wrap_info_values__set__n_knots

subroutine f90wrap_info_values__get__nt(this, f90wrap_nt)
    use interpolation, only: info_values
    implicit none
    type info_values_ptr_type
        type(info_values), pointer :: p => NULL()
    end type info_values_ptr_type
    integer, intent(in)   :: this(2)
    type(info_values_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_nt
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_nt = this_ptr%p%nt
end subroutine f90wrap_info_values__get__nt

subroutine f90wrap_info_values__set__nt(this, f90wrap_nt)
    use interpolation, only: info_values
    implicit none
    type info_values_ptr_type
        type(info_values), pointer :: p => NULL()
    end type info_values_ptr_type
    integer, intent(in)   :: this(2)
    type(info_values_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_nt
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%nt = f90wrap_nt
end subroutine f90wrap_info_values__set__nt

subroutine f90wrap_info_values__get__kl(this, f90wrap_kl)
    use interpolation, only: info_values
    implicit none
    type info_values_ptr_type
        type(info_values), pointer :: p => NULL()
    end type info_values_ptr_type
    integer, intent(in)   :: this(2)
    type(info_values_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_kl
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_kl = this_ptr%p%kl
end subroutine f90wrap_info_values__get__kl

subroutine f90wrap_info_values__set__kl(this, f90wrap_kl)
    use interpolation, only: info_values
    implicit none
    type info_values_ptr_type
        type(info_values), pointer :: p => NULL()
    end type info_values_ptr_type
    integer, intent(in)   :: this(2)
    type(info_values_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_kl
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%kl = f90wrap_kl
end subroutine f90wrap_info_values__set__kl

subroutine f90wrap_info_values__get__ku(this, f90wrap_ku)
    use interpolation, only: info_values
    implicit none
    type info_values_ptr_type
        type(info_values), pointer :: p => NULL()
    end type info_values_ptr_type
    integer, intent(in)   :: this(2)
    type(info_values_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_ku
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_ku = this_ptr%p%ku
end subroutine f90wrap_info_values__get__ku

subroutine f90wrap_info_values__set__ku(this, f90wrap_ku)
    use interpolation, only: info_values
    implicit none
    type info_values_ptr_type
        type(info_values), pointer :: p => NULL()
    end type info_values_ptr_type
    integer, intent(in)   :: this(2)
    type(info_values_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_ku
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%ku = f90wrap_ku
end subroutine f90wrap_info_values__set__ku

subroutine f90wrap_info_values__get__nleft(this, f90wrap_nleft)
    use interpolation, only: info_values
    implicit none
    type info_values_ptr_type
        type(info_values), pointer :: p => NULL()
    end type info_values_ptr_type
    integer, intent(in)   :: this(2)
    type(info_values_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_nleft
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_nleft = this_ptr%p%nleft
end subroutine f90wrap_info_values__get__nleft

subroutine f90wrap_info_values__set__nleft(this, f90wrap_nleft)
    use interpolation, only: info_values
    implicit none
    type info_values_ptr_type
        type(info_values), pointer :: p => NULL()
    end type info_values_ptr_type
    integer, intent(in)   :: this(2)
    type(info_values_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_nleft
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%nleft = f90wrap_nleft
end subroutine f90wrap_info_values__set__nleft

subroutine f90wrap_info_values__get__nright(this, f90wrap_nright)
    use interpolation, only: info_values
    implicit none
    type info_values_ptr_type
        type(info_values), pointer :: p => NULL()
    end type info_values_ptr_type
    integer, intent(in)   :: this(2)
    type(info_values_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_nright
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_nright = this_ptr%p%nright
end subroutine f90wrap_info_values__get__nright

subroutine f90wrap_info_values__set__nright(this, f90wrap_nright)
    use interpolation, only: info_values
    implicit none
    type info_values_ptr_type
        type(info_values), pointer :: p => NULL()
    end type info_values_ptr_type
    integer, intent(in)   :: this(2)
    type(info_values_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_nright
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%nright = f90wrap_nright
end subroutine f90wrap_info_values__set__nright

subroutine f90wrap_info_values__get__LHS_rows(this, f90wrap_LHS_rows)
    use interpolation, only: info_values
    implicit none
    type info_values_ptr_type
        type(info_values), pointer :: p => NULL()
    end type info_values_ptr_type
    integer, intent(in)   :: this(2)
    type(info_values_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_LHS_rows
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_LHS_rows = this_ptr%p%LHS_rows
end subroutine f90wrap_info_values__get__LHS_rows

subroutine f90wrap_info_values__set__LHS_rows(this, f90wrap_LHS_rows)
    use interpolation, only: info_values
    implicit none
    type info_values_ptr_type
        type(info_values), pointer :: p => NULL()
    end type info_values_ptr_type
    integer, intent(in)   :: this(2)
    type(info_values_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_LHS_rows
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%LHS_rows = f90wrap_LHS_rows
end subroutine f90wrap_info_values__set__LHS_rows

subroutine f90wrap_info_values__get__LHS_cols(this, f90wrap_LHS_cols)
    use interpolation, only: info_values
    implicit none
    type info_values_ptr_type
        type(info_values), pointer :: p => NULL()
    end type info_values_ptr_type
    integer, intent(in)   :: this(2)
    type(info_values_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_LHS_cols
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_LHS_cols = this_ptr%p%LHS_cols
end subroutine f90wrap_info_values__get__LHS_cols

subroutine f90wrap_info_values__set__LHS_cols(this, f90wrap_LHS_cols)
    use interpolation, only: info_values
    implicit none
    type info_values_ptr_type
        type(info_values), pointer :: p => NULL()
    end type info_values_ptr_type
    integer, intent(in)   :: this(2)
    type(info_values_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_LHS_cols
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%LHS_cols = f90wrap_LHS_cols
end subroutine f90wrap_info_values__set__LHS_cols

subroutine f90wrap_info_values__get__RHS_cols(this, f90wrap_RHS_cols)
    use interpolation, only: info_values
    implicit none
    type info_values_ptr_type
        type(info_values), pointer :: p => NULL()
    end type info_values_ptr_type
    integer, intent(in)   :: this(2)
    type(info_values_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_RHS_cols
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_RHS_cols = this_ptr%p%RHS_cols
end subroutine f90wrap_info_values__get__RHS_cols

subroutine f90wrap_info_values__set__RHS_cols(this, f90wrap_RHS_cols)
    use interpolation, only: info_values
    implicit none
    type info_values_ptr_type
        type(info_values), pointer :: p => NULL()
    end type info_values_ptr_type
    integer, intent(in)   :: this(2)
    type(info_values_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_RHS_cols
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%RHS_cols = f90wrap_RHS_cols
end subroutine f90wrap_info_values__set__RHS_cols

subroutine f90wrap_info_values_initialise(this)
    use interpolation, only: info_values
    implicit none
    
    type info_values_ptr_type
        type(info_values), pointer :: p => NULL()
    end type info_values_ptr_type
    type(info_values_ptr_type) :: this_ptr
    integer, intent(out), dimension(2) :: this
    allocate(this_ptr%p)
    this = transfer(this_ptr, this)
end subroutine f90wrap_info_values_initialise

subroutine f90wrap_info_values_finalise(this)
    use interpolation, only: info_values
    implicit none
    
    type info_values_ptr_type
        type(info_values), pointer :: p => NULL()
    end type info_values_ptr_type
    type(info_values_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate(this_ptr%p)
end subroutine f90wrap_info_values_finalise

subroutine f90wrap_info_values_xinfo_size_array__array_getitem__items(f90wrap_this, f90wrap_i, itemsitem)
    
    use interpolation, only: info_size, info_values
    implicit none
    
    type info_values_xinfo_size_array
        type(info_values), dimension(info_size) :: items
    end type info_values_xinfo_size_array
    
    type info_values_xinfo_size_array_ptr_type
        type(info_values_xinfo_size_array), pointer :: p => NULL()
    end type info_values_xinfo_size_array_ptr_type
    type info_values_ptr_type
        type(info_values), pointer :: p => NULL()
    end type info_values_ptr_type
    integer, intent(in) :: f90wrap_this(2)
    type(info_values_xinfo_size_array_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_i
    integer, intent(out) :: itemsitem(2)
    type(info_values_ptr_type) :: items_ptr
    
    this_ptr = transfer(f90wrap_this, this_ptr)
    if (f90wrap_i < 1 .or. f90wrap_i > size(this_ptr%p%items)) then
        call f90wrap_abort("array index out of range")
    else
        items_ptr%p => this_ptr%p%items(f90wrap_i)
        itemsitem = transfer(items_ptr,itemsitem)
    endif
end subroutine f90wrap_info_values_xinfo_size_array__array_getitem__items

subroutine f90wrap_info_values_xinfo_size_array__array_setitem__items(f90wrap_this, f90wrap_i, itemsitem)
    
    use interpolation, only: info_size, info_values
    implicit none
    
    type info_values_xinfo_size_array
        type(info_values), dimension(info_size) :: items
    end type info_values_xinfo_size_array
    
    type info_values_xinfo_size_array_ptr_type
        type(info_values_xinfo_size_array), pointer :: p => NULL()
    end type info_values_xinfo_size_array_ptr_type
    type info_values_ptr_type
        type(info_values), pointer :: p => NULL()
    end type info_values_ptr_type
    integer, intent(in) :: f90wrap_this(2)
    type(info_values_xinfo_size_array_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_i
    integer, intent(in) :: itemsitem(2)
    type(info_values_ptr_type) :: items_ptr
    
    this_ptr = transfer(f90wrap_this, this_ptr)
    if (f90wrap_i < 1 .or. f90wrap_i > size(this_ptr%p%items)) then
        call f90wrap_abort("array index out of range")
    else
        items_ptr = transfer(itemsitem,items_ptr)
        this_ptr%p%items(f90wrap_i) = items_ptr%p
    endif
end subroutine f90wrap_info_values_xinfo_size_array__array_setitem__items

subroutine f90wrap_info_values_xinfo_size_array__array_len__items(f90wrap_this, f90wrap_n)
    
    use interpolation, only: info_size, info_values
    implicit none
    
    type info_values_xinfo_size_array
        type(info_values), dimension(info_size) :: items
    end type info_values_xinfo_size_array
    
    type info_values_xinfo_size_array_ptr_type
        type(info_values_xinfo_size_array), pointer :: p => NULL()
    end type info_values_xinfo_size_array_ptr_type
    type info_values_ptr_type
        type(info_values), pointer :: p => NULL()
    end type info_values_ptr_type
    integer, intent(out) :: f90wrap_n
    integer, intent(in) :: f90wrap_this(2)
    type(info_values_xinfo_size_array_ptr_type) :: this_ptr
    
    this_ptr = transfer(f90wrap_this, this_ptr)
    f90wrap_n = size(this_ptr%p%items)
end subroutine f90wrap_info_values_xinfo_size_array__array_len__items

subroutine f90wrap_info_values_xinfo_size_array_initialise(this)
    use interpolation, only: info_size, info_values
    implicit none
    
    type info_values_xinfo_size_array
        type(info_values), dimension(info_size) :: items
    end type info_values_xinfo_size_array
    
    type info_values_xinfo_size_array_ptr_type
        type(info_values_xinfo_size_array), pointer :: p => NULL()
    end type info_values_xinfo_size_array_ptr_type
    type(info_values_xinfo_size_array_ptr_type) :: this_ptr
    integer, intent(out), dimension(2) :: this
    allocate(this_ptr%p)
    this = transfer(this_ptr, this)
end subroutine f90wrap_info_values_xinfo_size_array_initialise

subroutine f90wrap_info_values_xinfo_size_array_finalise(this)
    use interpolation, only: info_size, info_values
    implicit none
    
    type info_values_xinfo_size_array
        type(info_values), dimension(info_size) :: items
    end type info_values_xinfo_size_array
    
    type info_values_xinfo_size_array_ptr_type
        type(info_values_xinfo_size_array), pointer :: p => NULL()
    end type info_values_xinfo_size_array_ptr_type
    type(info_values_xinfo_size_array_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate(this_ptr%p)
end subroutine f90wrap_info_values_xinfo_size_array_finalise

subroutine f90wrap_pack_values(deg, n_knots, nt, kl, ku, nleft, nright, lhs_rows, lhs_cols, rhs_cols, info)
    use interpolation, only: pack_values, info_size, info_values
    implicit none
    
    type info_values_xinfo_size_array
        type(info_values), dimension(info_size) :: items
    end type info_values_xinfo_size_array
    
    type info_values_xinfo_size_array_ptr_type
        type(info_values_xinfo_size_array), pointer :: p => NULL()
    end type info_values_xinfo_size_array_ptr_type
    integer, intent(in) :: deg
    integer, intent(in) :: n_knots
    integer, intent(in) :: nt
    integer, intent(in) :: kl
    integer, intent(in) :: ku
    integer, intent(in) :: nleft
    integer, intent(in) :: nright
    integer, intent(in) :: lhs_rows
    integer, intent(in) :: lhs_cols
    integer, intent(in) :: rhs_cols
    type(info_values_xinfo_size_array_ptr_type) :: info_ptr
    integer, intent(out), dimension(2) :: info
    allocate(info_ptr%p)
    call pack_values(deg=deg, n_knots=n_knots, nt=nt, kl=kl, ku=ku, nleft=nleft, nright=nright, LHS_rows=lhs_rows, &
        LHS_cols=lhs_cols, RHS_cols=rhs_cols, info=info_ptr%p%items)
    info = transfer(info_ptr, info)
end subroutine f90wrap_pack_values

subroutine f90wrap_unpack_values(info, deg, n_knots, nt, kl, ku, nleft, nright, lhs_rows, lhs_cols, rhs_cols)
    use interpolation, only: unpack_values, info_values
    implicit none
    
    type info_values_ptr_type
        type(info_values), pointer :: p => NULL()
    end type info_values_ptr_type
    type(info_values_ptr_type) :: info_ptr
    integer, intent(in), dimension(2) :: info
    integer, intent(out) :: deg
    integer, intent(out) :: n_knots
    integer, intent(out) :: nt
    integer, intent(out) :: kl
    integer, intent(out) :: ku
    integer, intent(out) :: nleft
    integer, intent(out) :: nright
    integer, intent(out) :: lhs_rows
    integer, intent(out) :: lhs_cols
    integer, intent(out) :: rhs_cols
    info_ptr = transfer(info, info_ptr)
    call unpack_values(info=info_ptr%p, deg=deg, n_knots=n_knots, nt=nt, kl=kl, ku=ku, nleft=nleft, nright=nright, &
        LHS_rows=lhs_rows, LHS_cols=lhs_cols, RHS_cols=rhs_cols)
end subroutine f90wrap_unpack_values

subroutine f90wrap_interpolation__get__sp(f90wrap_sp)
    use interpolation, only: interpolation_sp => sp
    implicit none
    integer, intent(out) :: f90wrap_sp
    
    f90wrap_sp = interpolation_sp
end subroutine f90wrap_interpolation__get__sp

subroutine f90wrap_interpolation__get__dp(f90wrap_dp)
    use interpolation, only: interpolation_dp => dp
    implicit none
    integer, intent(out) :: f90wrap_dp
    
    f90wrap_dp = interpolation_dp
end subroutine f90wrap_interpolation__get__dp

subroutine f90wrap_interpolation__get__qp(f90wrap_qp)
    use interpolation, only: interpolation_qp => qp
    implicit none
    integer, intent(out) :: f90wrap_qp
    
    f90wrap_qp = interpolation_qp
end subroutine f90wrap_interpolation__get__qp

subroutine f90wrap_interpolation__get__info_size(f90wrap_info_size)
    use interpolation, only: interpolation_info_size => info_size
    implicit none
    integer, intent(out) :: f90wrap_info_size
    
    f90wrap_info_size = interpolation_info_size
end subroutine f90wrap_interpolation__get__info_size

! End of module interpolation defined in file test_deriv_types.f90

