!> Procedures to create, modify and study arrays.
!> Author: Adria Meca Montserrat.
!> Last modified date: 12/08/22.
module array_procedures
  use derived_types, only: int_list

  implicit none

  private

  interface add
    module procedure add_dbl, add_int, add_int_list
  end interface add

  interface argsort
    module procedure dbl_argsort, int_argsort
  end interface argsort

  interface find
    module procedure find_int
  end interface find

  interface my_pack
    module procedure my_pack_int
  end interface my_pack

  public :: add, argsort, find, my_pack, pop, quicksort

contains
  !> Looks for the position of an integer in a list of integers.
  function find_int(list, element)
    integer, intent(in) :: element   !>
    integer, intent(in) :: list(:)   !>
    integer             :: find_int  !> Output.
    integer             :: i         !>

    find_int = 0
    do i = 1, size(list)
      if (list(i) == element) then
        find_int = i
        exit
      end if
    end do
  end function find_int


  !> Returns the indices that would sort a list of doubles in ascending order.
  function dbl_argsort(list)
    double precision, intent(in) :: list(:)                  !>
    double precision             :: copy_list(size(list))    !>
    integer                      :: dbl_argsort(size(list))  !> Output.
    integer                      :: i, j                     !>

    !> Insertion sort algorithm.
    copy_list = list
    dbl_argsort = [(i, i=1,size(list))]
    do i = 2, size(list)
      j = i
      do while ((j > 1) .and. (copy_list(j-1) > copy_list(j)))
        copy_list([j-1, j]) = copy_list([j, j-1])
        dbl_argsort([j-1, j]) = dbl_argsort([j, j-1])
        j = j - 1
      end do
    end do
  end function dbl_argsort


  !> Returns the indices that would sort a list of integers in ascending order.
  function int_argsort(list)
    integer, intent(in) :: list(:)                  !>
    integer             :: copy_list(size(list))    !>
    integer             :: i, j                     !>
    integer             :: int_argsort(size(list))  !> Output.

    !> Insertion sort algorithm.
    copy_list = list
    int_argsort = [(i, i=1,size(list))]
    do i = 2, size(list)
      j = i
      do while ((j > 1) .and. (copy_list(j-1) > copy_list(j)))
        copy_list([j-1, j]) = copy_list([j, j-1])
        int_argsort([j-1, j]) = int_argsort([j, j-1])
        j = j - 1
      end do
    end do
  end function int_argsort


  !> Returns a list without its i-th element.
  function pop(list, i)
    integer,                       intent(in) :: i        !>
    integer                                   :: n        !>
    double precision,              intent(in) :: list(:)  !>
    double precision, allocatable             :: pop(:)   !> Output.

    n = size(list)
    if ((1 <= i) .and. (i <= n)) then
      allocate(pop(n-1))
      pop = [list(:i-1), list(i+1:)]
    else
      allocate(pop(n))
      pop = list
    end if
  end function pop


  !> Adds a double to an array of doubles.
  subroutine add_dbl(list, element, lb_)
    integer,          optional,    intent(in)    :: lb_           !> Lower bound of 'list'.
    integer                                      :: isize, lb     !>
    double precision,              intent(in)    :: element       !>
    double precision, allocatable, intent(inout) :: list(:)       !>
    double precision, allocatable                :: copy_list(:)  !>

    if (allocated(list)) then
      lb = lbound(list, dim=1)

      !> We make a copy of the original list, adding the new element to it.
      isize = size(list)
      allocate(copy_list(lb:lb+isize))
      copy_list(lb:lb+isize-1) = list
      copy_list(lb+isize) = element

      !> We move the elements of the copied list to the original one.
      call move_alloc(copy_list, list)
    else
      lb = merge(lb_, 1, present(lb_))

      !> If the original list has no elements, we add the new element to it.
      allocate(list(lb:lb))
      list(lb) = element
    end if
  end subroutine add_dbl


  !> Adds an integer to an array of integers.
  subroutine add_int(list, element, lb_)
    integer,              intent(in)    :: element       !>
    integer, optional,    intent(in)    :: lb_           !> Lower bound of 'list'.
    integer, allocatable, intent(inout) :: list(:)       !>
    integer, allocatable                :: copy_list(:)  !>
    integer                             :: isize, lb     !>

    if (allocated(list)) then
      lb = lbound(list, dim=1)

      !> We make a copy of the original list, adding the new element to it.
      isize = size(list)
      allocate(copy_list(lb:lb+isize))
      copy_list(lb:lb+isize-1) = list
      copy_list(lb+isize) = element

      !> We move the elements of the copied list to the original one.
      call move_alloc(copy_list, list)
    else
      lb = merge(lb_, 1, present(lb_))

      !> If the original list has no elements, we add the new element to it.
      allocate(list(lb:lb))
      list(lb) = element
    end if
  end subroutine add_int


  !> Adds a list of integers to a list of lists (of integers).
  subroutine add_int_list(list, element, lb_)
    integer,        optional,    intent(in)    :: lb_           !> Lower bound of 'list'.
    integer                                    :: isize, lb     !>
    type(int_list),              intent(in)    :: element       !>
    type(int_list), allocatable, intent(inout) :: list(:)       !>
    type(int_list), allocatable                :: copy_list(:)  !>

    if (allocated(list)) then
      lb = lbound(list, dim=1)

      !> We make a copy of the original list, adding the new element to it.
      isize = size(list)
      allocate(copy_list(lb:lb+isize))
      copy_list(lb:lb+isize-1) = list
      copy_list(lb+isize) = element

      !> We move the elements of the copied list to the original one.
      call move_alloc(copy_list, list)
    else
      lb = merge(lb_, 1, present(lb_))

      !> If the original list has no elements, we add the new element to it.
      allocate(list(lb:lb))
      list(lb) = element
    end if
  end subroutine add_int_list


  !> Custom version of GFortran's PACK function that applies to lists of integers.
  subroutine my_pack_int(list, condition)
    integer, allocatable, intent(inout) :: list(:)                   !>
    integer, allocatable                :: copy_list(:)              !>
    integer                             :: i, j, new_size, old_size  !>
    logical,              intent(in)    :: condition(:)              !>

    if (.not. allocated(list)) then
      print *, 'An attempt was made to filter an unallocated list.'
      return
    end if

    !> We make a copy of the original list.
    old_size = size(list)
    allocate(copy_list(old_size))
    copy_list = list

    !> The new size of the original list is equal to the number of elements that
    !> meet the condition.
    new_size = count(condition)
    deallocate(list)
    allocate(list(new_size))

    !> We construct the new list.
    if (new_size > 0) then
      j = 1
      do i = 1, old_size
        if (condition(i)) then
          list(j) = copy_list(i)
          j = j + 1
        end if
      end do
    end if
  end subroutine my_pack_int


  !> Sorts two lists by increasing value of the former.
  recursive subroutine quicksort(energies, node_ids)
    integer,                       intent(inout) :: node_ids(:)                !>
    integer,          allocatable                :: n_eq(:), n_gt(:), n_lt(:)  !>
    integer                                      :: i, isize, n                !>
    double precision,              intent(inout) :: energies(:)                !>
    double precision, allocatable                :: e_eq(:), e_gt(:), e_lt(:)  !>
    double precision, parameter                  :: eps = 1.0d-12              !>
    double precision                             :: e, pivot                   !>

    !> If the list of energies has less than two elements, it is already sorted.
    isize = size(energies)
    if (isize > 1) then
      !> We arbitrarily choose a pivot.
      pivot = energies(1)

      !> We classify the energies according to whether they are equal to, less than
      !> or greater than the pivot. Nodes are classified in the same way.
      do i = 1, isize
        e = energies(i)
        n = node_ids(i)
        if (abs(e - pivot) < eps) then
          call add(e_eq, e)
          call add(n_eq, n)
        else if (e < pivot) then
          call add(e_lt, e)
          call add(n_lt, n)
        else if (e > pivot) then
          call add(e_gt, e)
          call add(n_gt, n)
        end if
      end do

      !> We make sure that the lists containing the elements that are smaller and
      !> greater than the pivot have been allocated to avoid memory errors.
      if (.not. allocated(e_lt)) allocate(e_lt(0), n_lt(0))
      if (.not. allocated(e_gt)) allocate(e_gt(0), n_gt(0))

      !> We recursively apply the subroutine to the lists containing the elements
      !> that are smaller and greater than the pivot.
      call quicksort(e_lt, n_lt)
      call quicksort(e_gt, n_gt)

      !> Finally, we construct the sorted lists.
      energies = [e_lt, e_eq, e_gt]
      node_ids = [n_lt, n_eq, n_gt]
    end if
  end subroutine quicksort
end module array_procedures
