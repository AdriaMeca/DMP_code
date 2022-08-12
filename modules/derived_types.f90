!> Derived types used throughout the code.
!> Author: Adria Meca Montserrat.
!> Last modified date: 12/08/22.
module derived_types
  implicit none

  private

  type, public :: dbl_list
    double precision, allocatable :: array(:)
  end type dbl_list

  type, public :: int_list
    integer, allocatable :: array(:)
  end type int_list

  !> List of lists of integers.
  type, public :: int_llist
    type(int_list), allocatable :: time(:)
  end type int_llist

  type, public :: node
    integer, allocatable :: neighbors(:)
    integer, allocatable :: opposites(:)
  end type node

  !> Epidemiological parameters.
  type, public :: prm
    integer          :: t0      !> Observation time.
    double precision :: alpha   !> Probability that an S node becomes E by interacting with one of its E neighbors.
    double precision :: lambda  !> Probability that an S node becomes E by interacting with one of its I neighbors.
    double precision :: mu      !> Probability that an I node becomes R on its own.
    double precision :: nu      !> Probability that an E node becomes I on its own.
  end type prm
end module derived_types
