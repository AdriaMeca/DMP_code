!> Procedures to create, modify and study arrays.
!> Author: Adria Meca Montserrat.
!> Last modified date: 30/09/22.
module array_procedures
  use derived_types, only: int_list

  implicit none

  private

  interface add
    module procedure add_dbl, add_int
  end interface add

  interface find
    module procedure find_int, find_str
  end interface find

  interface my_pack
    module procedure my_pack_int
  end interface my_pack

  public :: add, find, my_pack, pop, quicksort

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


  !> Looks for the position of a string in a list of strings.
  function find_str(list, element)
    integer                      :: find_str  !> Output.
    integer                      :: t         !>
    character(len=*), intent(in) :: element   !>
    character(len=*), intent(in) :: list(:)   !>

    find_str = size(list)
    do t = 1, find_str
      if (list(t) == element) then
        find_str = t
        exit
      end if
    end do
  end function find_str


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
  subroutine add_dbl(list, element)
    integer                                      :: isize         !>
    double precision,              intent(in)    :: element       !>
    double precision, allocatable, intent(inout) :: list(:)       !>
    double precision, allocatable                :: copy_list(:)  !>

    if (allocated(list)) then
      !> We make a copy of the original list, adding the new element to it.
      isize = size(list)
      allocate(copy_list(1+isize))
      copy_list(1:isize) = list
      copy_list(1+isize) = element

      !> We move the elements of the copied list to the original one.
      call move_alloc(copy_list, list)
    else
      !> If the original list has no elements, we add the new element to it.
      allocate(list(1))
      list(1) = element
    end if
  end subroutine add_dbl


  !> Adds an integer to an array of integers.
  subroutine add_int(list, element)
    integer,              intent(in)    :: element       !>
    integer, allocatable, intent(inout) :: list(:)       !>
    integer, allocatable                :: copy_list(:)  !>
    integer                             :: isize         !>

    if (allocated(list)) then
      !> We make a copy of the original list, adding the new element to it.
      isize = size(list)
      allocate(copy_list(1+isize))
      copy_list(1:isize) = list
      copy_list(1+isize) = element

      !> We move the elements of the copied list to the original one.
      call move_alloc(copy_list, list)
    else
      !> If the original list has no elements, we add the new element to it.
      allocate(list(1))
      list(1) = element
    end if
  end subroutine add_int


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
