!Module whose procedures create, modify and study arrays.
module array_procedures
  implicit none

  private

  !List of doubles.
  type, public :: dbl_list
    double precision, dimension(:), allocatable :: array
  end type

  !List of integers.
  type, public :: int_list
    integer, dimension(:), allocatable :: array
  end type int_list

  !Pair of integers.
  type, public :: int_pair
    integer :: x, y
  end type int_pair

  interface add_item
    module procedure add_int, add_int_list, add_int_pair
  end interface add_item

  public add_item, del_int, pop

contains
  !Function that returns a list without the element placed at index.
  function pop(list, index) result(reduced_list)
    implicit none

    double precision, dimension(:), intent(in) :: list
    integer, intent(in) :: index

    double precision, dimension(:), allocatable :: reduced_list
    integer :: n

    n = size(list)
    if ((1 <= index).and.(index <= n)) then
      allocate(reduced_list(n-1))
      reduced_list = (/ list(:index-1), list(index+1:) /)
    else
      allocate(reduced_list(n))
      reduced_list = list
    end if
  end function pop


  !Subroutine that adds an integer to an array of integers.
  subroutine add_int(list, element)
    implicit none

    integer, dimension(:), allocatable, intent(inout) :: list
    integer, intent(in) :: element

    integer, dimension(:), allocatable :: copy_list
    integer :: isize

    if (allocated(list)) then
      isize = size(list)

      allocate(copy_list(1+isize))

      copy_list(1:isize) = list
      copy_list(1+isize) = element

      deallocate(list)
      allocate(list(1+isize))

      list = copy_list
      deallocate(copy_list)
    else
      allocate(list(1))
      list(1) = element
    end if
  end subroutine add_int


  !Subroutine that adds a list of integers to a list of lists.
  subroutine add_int_list(list, element)
    implicit none

    type(int_list), dimension(:), allocatable, intent(inout) :: list
    type(int_list), intent(in) :: element

    type(int_list), dimension(:), allocatable :: copy_list
    integer :: isize

    if (allocated(list)) then
      isize = size(list)

      allocate(copy_list(1+isize))

      copy_list(1:isize) = list
      copy_list(1+isize) = element

      deallocate(list)
      allocate(list(1+isize))

      list = copy_list
      deallocate(copy_list)
    else
      allocate(list(1))
      list(1) = element
    end if
  end subroutine add_int_list


  !Subroutine that adds an int_pair to a list of pairs.
  subroutine add_int_pair(list, element)
    implicit none

    type(int_pair), dimension(:), allocatable, intent(inout) :: list
    type(int_pair), intent(in) :: element

    type(int_pair), dimension(:), allocatable :: copy_list
    integer :: isize

    if (allocated(list)) then
      isize = size(list)

      allocate(copy_list(1+isize))

      copy_list(1:isize) = list
      copy_list(1+isize) = element

      deallocate(list)
      allocate(list(1+isize))

      list = copy_list
      deallocate(copy_list)
    else
      allocate(list(1))
      list(1) = element
    end if
  end subroutine add_int_pair


  !Subroutine that removes an integer from a list of integers.
  subroutine del_int(list, element)
    implicit none

    integer, dimension(:), allocatable, intent(inout) :: list
    integer, intent(in) :: element

    if (allocated(list)) list = pack(list, mask=list/=element)
  end subroutine del_int
end module array_procedures