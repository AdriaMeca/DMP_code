!Module whose procedures create, modify and study arrays.
module array_procedures
  implicit none

  private

  !List of doubles.
  type, public :: dbl_list
    double precision, dimension(:), allocatable :: array
  end type

  !Pair of integers.
  type, public :: int_int
    integer :: x, y
  end type int_int

  !List of integers.
  type, public :: int_list
    integer, dimension(:), allocatable :: array
  end type int_list
end module array_procedures