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
end module array_procedures