!> \file intarr_hashing.f90
!! \brief Module file for dictionary_t

!> Dictionary type that uses strings for the keys and values
!! Modified from fortran_hash_table by Pierre de Buyl
!! Downloaded 2020-12-04 from: https://github.com/pdebuyl/fortran_hash_table
!! See file intarr_hashing_LICENSE for licence.
!!
!! Design:
!!  - iahash hash function using array of integers
!!  - values are type real
!!  - There is no linked list nor pointers, only allocatable arrays for the dynamic data structure
!!  - set rewrites existing entries without complaining

module intarr_hashing
  implicit none

  private

  public :: dictionary_t

  !> Single entry in the dictionary
  type entry_t
     integer, allocatable :: key(:)
     real, allocatable :: value
  end type entry_t

  !> A bucket contains several entries
  type bucket_t
     type(entry_t), allocatable :: entries(:)
     integer :: current_size = 0
     integer :: current_idx = 0
   contains
     procedure :: find
  end type bucket_t

  !> The dictionary contains dict_size buckets (defined at run time)
  type dictionary_t
     type(bucket_t), allocatable :: buckets(:)
     integer :: dict_size = 0
   contains
!     procedure :: djb2
     procedure :: iahash
     procedure :: set
     procedure :: get
     procedure :: init
     procedure :: outdata
  end type dictionary_t

  integer, parameter :: BUCKET_EMPTY = -2
  integer, parameter :: BUCKET_ENTRY_NOT_FOUND = -4

contains

  !> djb2 hash function
  !!
  !! \param this the dictionary_t object
  !! \param s a string
  !!
  !! \return the hash value between 0 and dict_size-1
!  function djb2(this, s) result(r)
!    class(dictionary_t), intent(in) :: this
!    character(len=*), intent(in) :: s
!    integer :: r
!
!    integer :: i, l
!
!    l = len(s)
!
!    r = 5381
!
!    do i = 1, l
!       r = r*33 + ichar(s(i:i))
!    end do
!
!    r = modulo(r, this%dict_size)
!
!  end function djb2

  !> integer array hash function
  !!
  !! \param this the dictionary_t object
  !! \param s a string
  !!
  !! \return the hash value between 0 and dict_size-1
  function iahash(this, nk, k) result(r)
    class(dictionary_t), intent(in)   :: this
    integer, intent(in)                :: nk
    integer, dimension(nk), intent(in) :: k
    integer :: r

    integer :: i
    
    r = 33
    
    do i = 1, nk
       r = r*31 + k(i)
    end do
  
  r = modulo(r, this%dict_size)

  end function iahash

  !> Add or replace an entry in the dictionary
  !!
  !! \param this the dictionary_t object
  !! \param k the key
  !! \param v the value
  subroutine set(this, nk, k, v)
    class(dictionary_t), intent(inout) :: this
    integer, intent(in)                :: nk
    integer, dimension(nk), intent(in)  :: k
    real, intent(in)                   :: v

    type(bucket_t) :: tmp_bucket

    integer :: h, i, b_idx

    h = this%iahash(nk, k) + 1

    b_idx = this%buckets(h)%find(nk, k)

    if (b_idx == BUCKET_EMPTY) then
       ! allocate bucket for 1 entry
       ! also, means we can take the first entry
       allocate(this%buckets(h)%entries(1))
       this%buckets(h)%current_size = 1
       this%buckets(h)%current_idx = 1
       b_idx = 1
       this%buckets(h)%entries(1)%key = k
       this%buckets(h)%entries(1)%value = v
       ! the values are registered, exit
       return
    end if

    if (b_idx == BUCKET_ENTRY_NOT_FOUND) then
       ! copy and grow bucket entries
       
       allocate(tmp_bucket%entries(this%buckets(h)%current_size + 1))
       tmp_bucket%current_size = this%buckets(h)%current_size + 1
       tmp_bucket%current_idx = this%buckets(h)%current_idx + 1

       do i = 1, this%buckets(h)%current_size
          tmp_bucket%entries(i)%key = this%buckets(h)%entries(i)%key
          tmp_bucket%entries(i)%value = this%buckets(h)%entries(i)%value
       end do

       deallocate(this%buckets(h)%entries)
       allocate(this%buckets(h)%entries, source=tmp_bucket%entries)
       deallocate(tmp_bucket%entries)

       this%buckets(h)%current_size = tmp_bucket%current_size
       this%buckets(h)%current_idx = tmp_bucket%current_idx
       b_idx = this%buckets(h)%current_idx
    end if

    if (b_idx > 0) then
       this%buckets(h)%entries(b_idx)%key = k
       this%buckets(h)%entries(b_idx)%value = v
    end if

  end subroutine set

  !> Initialize a dictionary object
  !!
  !! \param this the dictionary_t object
  !! \param dict_size the size of the hash table
  subroutine init(this, dict_size)
    class(dictionary_t), intent(out) :: this
    integer, intent(in) :: dict_size

    allocate(this%buckets(dict_size))
    this%dict_size = dict_size

  end subroutine init

  !> Print the contents to file data.chem
  !!
  !! \param this the dictionary_t object
  subroutine outdata(this, fnum)
    class(dictionary_t), intent(in) :: this
	integer, intent(in) :: fnum

    integer :: i, j, s
    integer :: n
	character(len=10) :: file_id
	character(len=15) :: file_name
	
	write(file_id, '(i0)') fnum
	file_name = 'data' // trim(adjustl(file_id)) // '.chem'

	open (14, FILE=file_name, status='replace')
	write(14, *) "    R  | stored configurations"
	
    n = 0
    do i = 1, this%dict_size
       s = this%buckets(i)%current_idx
       if (s > 0) then
          do j = 1, s
             write(14, '(F7.4, " |", 500i4)') &
			   this%buckets(i)%entries(j)%value, &
			   this%buckets(i)%entries(j)%key
          end do
       end if
    end do
	
	close(14)

  end subroutine outdata

  !> Find the "in-bucket" index for a given key
  !!
  !! Negative return values correspond to module-defined return codes.
  !!
  !! \param this the bucket_t object
  !! \param k the key
  !!
  !! \return the index (1-based) of the key in the bucket or a return code
  function find(this, nk, k) result(r)
    class(bucket_t), intent(in)   :: this
    integer, intent(in)                :: nk
    integer, dimension(nk), intent(in)  :: k
    integer :: r

    integer :: i

    if (this%current_size == 0) then
       r = BUCKET_EMPTY
       return
    end if

    r = BUCKET_ENTRY_NOT_FOUND
    do i = 1, this%current_size
       if (all(this%entries(i)%key == k)) then
          r = i
          exit
       end if
    end do

  end function find

  !> Fetch an entry in the dictionary.
  !!
  !! \param this the dictionary_t object
  !! \param nk is the size of the k array
  !! \param k the key
  !!
  !! \return the value if found, 0.0 else
  function get(this, nk, k) result(r)
    class(dictionary_t), intent(in)    :: this
    integer, intent(in)                :: nk
    integer, dimension(nk), intent(in) :: k

    real    :: r

    integer :: h, b_idx

    h = this%iahash(nk, k) + 1

    b_idx = this%buckets(h)%find(nk, k)

    if ( (b_idx == BUCKET_EMPTY) .or. &
         (b_idx == BUCKET_ENTRY_NOT_FOUND) ) then
       r = 0.0
       return
    end if

    if (b_idx>0) then
       r = this%buckets(h)%entries(b_idx)%value
    end if
    
  end function get

end module intarr_hashing
