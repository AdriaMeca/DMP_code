!Module whose procedures create, modify and study arrays.
!Author: Adrià Meca Montserrat.
!Last modified date: 16/05/22.
module array_procedures
  implicit none

  private

  !List of doubles.
  type, public :: dbl_list
    double precision, dimension(:), allocatable :: array
  end type dbl_list

  !List of integers.
  type, public :: int_list
    integer, dimension(:), allocatable :: array
  end type int_list

  !List of lists of integers.
  type, public :: int_list_list
    type(int_list), dimension(:), allocatable :: time
  end type int_list_list

  !Pair of integers.
  type, public :: int_pair
    integer :: x, y
  end type int_pair

  interface add_item
    module procedure add_dbl, add_int, add_int_list, add_int_pair
  end interface add_item

  interface idx_insertion_sort
    module procedure dbl_idx_insertion_sort, int_idx_insertion_sort
  end interface idx_insertion_sort

  public add_item, del_int, find_int, idx_insertion_sort, pop, quicksort

contains
  !Function that searches the index of an integer in a list of integers.
  function find_int(list, element) result(index)
    implicit none

    integer, dimension(:), intent(in) :: list
    integer, intent(in) :: element

    integer :: i, index

    index = 0
    do i = 1, size(list)
      if (list(i) == element) then
        index = i
        exit
      end if
    end do
  end function find_int


  !Function that returns the permutation of the indices that orders its
  !associated list of doubles.
  function dbl_idx_insertion_sort(list) result(indices)
    implicit none

    double precision, dimension(:), intent(in) :: list
    double precision, dimension(size(list)) :: copy_list
    integer, dimension(size(list)) :: indices

    integer :: i, j

    copy_list = list
    indices = [(i, i=1,size(list))]
    do i = 2, size(list)
      j = i
      do while (copy_list(j) < copy_list(j-1))
        copy_list([j-1, j]) = copy_list([j, j-1])
        indices([j-1, j]) = indices([j, j-1])
        j = j - 1
        if (j <= 1) exit
      end do
    end do
  end function dbl_idx_insertion_sort


  !Function that returns the permutation of the indices that orders its
  !associated list of integers.
  function int_idx_insertion_sort(list) result(indices)
    implicit none

    integer, dimension(:), intent(in) :: list
    integer, dimension(size(list)) :: copy_list, indices

    integer :: i, j

    copy_list = list
    indices = [(i, i=1,size(list))]
    do i = 2, size(list)
      j = i
      do while (copy_list(j) < copy_list(j-1))
        copy_list([j-1, j]) = copy_list([j, j-1])
        indices([j-1, j]) = indices([j, j-1])
        j = j - 1
        if (j <= 1) exit
      end do
    end do
  end function int_idx_insertion_sort


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
      reduced_list = [list(:index-1), list(index+1:)]
    else
      allocate(reduced_list(n))
      reduced_list = list
    end if
  end function pop


  !Subroutine that adds a double to an array of doubles.
  subroutine add_dbl(list, element)
    implicit none

    double precision, dimension(:), allocatable, intent(inout) :: list
    double precision, intent(in) :: element

    double precision, dimension(:), allocatable :: copy_list
    integer :: isize

    if (allocated(list)) then
      !We create a temporary list with one more element than the original.
      isize = size(list)
      allocate(copy_list(1+isize))

      !We copy the elements from the old list to the new one and add the
      !new element.
      copy_list(1:isize) = list
      copy_list(1+isize) = element

      !We deallocate the original list, make it bigger and transfer all
      !the elements to it, including the new one. Then, we deallocate
      !the temporary list because we do not need it anymore.
      deallocate(list)
      allocate(list(1+isize))
      list = copy_list
      deallocate(copy_list)
    else
      !If the original list has no elements, we pass the new item to it.
      allocate(list(1))
      list(1) = element
    end if
  end subroutine add_dbl


  !Subroutine that adds an integer to an array of integers.
  subroutine add_int(list, element)
    implicit none

    integer, dimension(:), allocatable, intent(inout) :: list
    integer, intent(in) :: element

    integer, dimension(:), allocatable :: copy_list
    integer :: isize

    if (allocated(list)) then
      !We create a temporary list with one more element than the original.
      isize = size(list)
      allocate(copy_list(1+isize))

      !We copy the elements from the old list to the new one and add the
      !new element.
      copy_list(1:isize) = list
      copy_list(1+isize) = element

      !We deallocate the original list, make it bigger and transfer all
      !the elements to it, including the new one. Then, we deallocate
      !the temporary list because we do not need it anymore.
      deallocate(list)
      allocate(list(1+isize))
      list = copy_list
      deallocate(copy_list)
    else
      !If the original list has no elements, we pass the new item to it.
      allocate(list(1))
      list(1) = element
    end if
  end subroutine add_int


  !Subroutine that adds a list of integers to a list of lists (of integers).
  subroutine add_int_list(list, element)
    implicit none

    type(int_list), dimension(:), allocatable, intent(inout) :: list
    type(int_list), intent(in) :: element

    type(int_list), dimension(:), allocatable :: copy_list
    integer :: isize

    if (allocated(list)) then
      !We create a temporary list with one more element than the original.
      isize = size(list)
      allocate(copy_list(1+isize))

      !We copy the elements from the old list to the new one and add the
      !new element.
      copy_list(1:isize) = list
      copy_list(1+isize) = element

      !We deallocate the original list, make it bigger and transfer all
      !the elements to it, including the new one. Then, we deallocate
      !the temporary list because we do not need it anymore.
      deallocate(list)
      allocate(list(1+isize))
      list = copy_list
      deallocate(copy_list)
    else
      !If the original list has no elements, we pass the new item to it.
      allocate(list(1))
      list(1) = element
    end if
  end subroutine add_int_list


  !Subroutine that adds an int_pair to a list of pairs (of integers).
  subroutine add_int_pair(list, element)
    implicit none

    type(int_pair), dimension(:), allocatable, intent(inout) :: list
    type(int_pair), intent(in) :: element

    type(int_pair), dimension(:), allocatable :: copy_list
    integer :: isize

    if (allocated(list)) then
      !We create a temporary list with one more element than the original.
      isize = size(list)
      allocate(copy_list(1+isize))

      !We copy the elements from the old list to the new one and add the
      !new element.
      copy_list(1:isize) = list
      copy_list(1+isize) = element

      !We deallocate the original list, make it bigger and transfer all
      !the elements to it, including the new one. Then, we deallocate
      !the temporary list because we do not need it anymore.
      deallocate(list)
      allocate(list(1+isize))
      list = copy_list
      deallocate(copy_list)
    else
      !If the original list has no elements, we pass the new item to it.
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


  !Recursive subroutine that sorts two lists in ascending order from the values
  !of the first.
  recursive subroutine quicksort(energies, node_ids)
    implicit none

    double precision, dimension(:), intent(inout) :: energies
    integer, dimension(:), intent(inout) :: node_ids

    double precision, dimension(:), allocatable :: e_eq, e_gt, e_lt
    double precision, parameter :: eps=1.0d-12
    double precision :: e, pivot
    integer, dimension(:), allocatable :: n_eq, n_gt, n_lt
    integer :: i, isize, n

    !If the size of the list of energies is less than or equal to 1 we do not
    !have to do anything to it because it is already sorted.
    isize = size(energies)
    if (isize > 1) then
      !We choose a pivot arbitrarily.
      pivot = energies(1)

      !We go through all the elements of the list of energies, classifying
      !them according to whether they are equal to, less than or greater
      !than the pivot; we classify node identifiers in the same way.
      do i = 1, isize
        e = energies(i)
        n = node_ids(i)
        if (abs(e-pivot) < eps) then
          call add_item(e_eq, e)
          call add_item(n_eq, n)
        else if (e < pivot) then
          call add_item(e_lt, e)
          call add_item(n_lt, n)
        else if (e > pivot) then
          call add_item(e_gt, e)
          call add_item(n_gt, n)
        end if
      end do

      !If the list containing the elements that are greater (smaller) than the
      !pivot has not been allocated, we allocate it to avoid memory errors.
      if (.not.allocated(e_lt)) allocate(e_lt(0), n_lt(0))
      if (.not.allocated(e_gt)) allocate(e_gt(0), n_gt(0))

      !We apply the subroutine recursively to the lists that contain the elements
      !that are smaller and greater than the pivot, respectively.
      call quicksort(e_lt, n_lt)
      call quicksort(e_gt, n_gt)

      !Finally, we construct the sorted lists using automatic reallocation.
      energies = [e_lt, e_eq, e_gt]
      node_ids = [n_lt, n_eq, n_gt]
    end if
  end subroutine quicksort
end module array_procedures
