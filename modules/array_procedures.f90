!Module whose procedures create, modify and study arrays.
!Author: Adrià Meca Montserrat.
!Last modified date: 25/05/22.
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
  type, public :: int_llist
    type(int_list), dimension(:), allocatable :: time
  end type int_llist

  !Pair of integers.
  type, public :: int_pair
    integer :: x, y
  end type int_pair

  !Interface that adds elements to lists.
  interface add
    module procedure add_dbl, add_int, add_int_list, add_int_pair
  end interface add

  !Interface that searches for elements in lists.
  interface find
    module procedure find_int
  end interface find

  !Interface that returns the indices that sort a given list in ascending order.
  interface idx_insertion_sort
    module procedure dbl_idx_insertion_sort, int_idx_insertion_sort
  end interface idx_insertion_sort

  !Interface that filters different lists given a condition.
  interface my_pack
    module procedure my_pack_int, my_pack_int_pair
  end interface my_pack

  public add, find, idx_insertion_sort, my_pack, pop, quicksort

contains
  !Function that searches the index of an integer in a list of integers.
  function find_int(list, element)
    implicit none

    !Input arguments.
    integer, dimension(:), intent(in) :: list
    integer, intent(in) :: element

    !Output arguments.
    integer :: find_int

    !Local variables.
    integer :: i

    find_int = 0
    do i = 1, size(list)
      if (list(i) == element) then
        find_int = i
        exit
      end if
    end do
  end function find_int


  !Function that returns the permutation of the indices that orders its
  !associated list of doubles.
  function dbl_idx_insertion_sort(list)
    implicit none

    !Input arguments.
    double precision, dimension(:), intent(in) :: list

    !Output arguments.
    integer, dimension(size(list)) :: dbl_idx_insertion_sort

    !Local variables.
    double precision, dimension(size(list)) :: copy_list

    integer :: i, j

    copy_list = list
    dbl_idx_insertion_sort = [(i, i=1,size(list))]
    do i = 2, size(list)
      j = i
      do while (copy_list(j) < copy_list(j-1))
        copy_list([j-1, j]) = copy_list([j, j-1])
        dbl_idx_insertion_sort([j-1, j]) = dbl_idx_insertion_sort([j, j-1])
        j = j - 1
        if (j <= 1) exit
      end do
    end do
  end function dbl_idx_insertion_sort


  !Function that returns the permutation of the indices that orders its
  !associated list of integers.
  function int_idx_insertion_sort(list)
    implicit none

    !Input arguments.
    integer, dimension(:), intent(in) :: list

    !Output arguments.
    integer, dimension(size(list)) :: int_idx_insertion_sort

    !Local variables.
    integer, dimension(size(list)) :: copy_list
    integer :: i, j

    copy_list = list
    int_idx_insertion_sort = [(i, i=1,size(list))]
    do i = 2, size(list)
      j = i
      do while (copy_list(j) < copy_list(j-1))
        copy_list([j-1, j]) = copy_list([j, j-1])
        int_idx_insertion_sort([j-1, j]) = int_idx_insertion_sort([j, j-1])
        j = j - 1
        if (j <= 1) exit
      end do
    end do
  end function int_idx_insertion_sort


  !Function that returns a list without the element placed at index.
  function pop(list, index)
    implicit none

    !Input arguments.
    double precision, dimension(:), intent(in) :: list

    integer, intent(in) :: index

    !Output arguments.
    double precision, dimension(:), allocatable :: pop

    !Local variables.
    integer :: n

    n = size(list)
    if ((1 <= index).and.(index <= n)) then
      allocate(pop(n-1))
      pop = [list(:index-1), list(index+1:)]
    else
      allocate(pop(n))
      pop = list
    end if
  end function pop


  !Subroutine that adds a double to an array of doubles.
  subroutine add_dbl(list, element)
    implicit none

    !Input/output arguments.
    double precision, dimension(:), allocatable, intent(inout) :: list
    double precision, intent(in) :: element

    !Local variables.
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

    !Input/output arguments.
    integer, dimension(:), allocatable, intent(inout) :: list
    integer, intent(in) :: element

    !Local variables.
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

    !Input/output arguments.
    type(int_list), dimension(:), allocatable, intent(inout) :: list
    type(int_list), intent(in) :: element

    !Local variables.
    integer :: isize

    type(int_list), dimension(:), allocatable :: copy_list

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

    !Input/output arguments.
    type(int_pair), dimension(:), allocatable, intent(inout) :: list
    type(int_pair), intent(in) :: element

    !Local variables.
    integer :: isize

    type(int_pair), dimension(:), allocatable :: copy_list

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


  !Custom version of the PACK function provided by GFortran that applies to
  !lists of integers.
  subroutine my_pack_int(list, condition)
    implicit none

    !Input/output arguments.
    integer, dimension(:), allocatable, intent(inout) :: list

    logical, dimension(:), intent(in) :: condition

    !Local variables.
    integer, dimension(:), allocatable :: copy_list
    integer :: i, j, new_size, old_size

    if (.not.allocated(list)) then
      print*, 'An attempt was made to filter an unallocated list.'
      return
    end if

    !We transfer all elements from list to copy_list.
    old_size = size(list)
    allocate(copy_list(old_size))
    copy_list = list

    !We compute the new size of the list based on the number of true elements
    !in the boolean array called condition.
    new_size = count(condition)
    deallocate(list)
    allocate(list(new_size))

    !If the new size is greater than zero, we pass all the elements that meet
    !the condition to our list.
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


  !Custom version of the PACK function provided by GFortran that applies to
  !lists of integer pairs.
  subroutine my_pack_int_pair(list, condition)
    implicit none

    !Input/output arguments.
    logical, dimension(:), intent(in) :: condition

    type(int_pair), dimension(:), allocatable, intent(inout) :: list

    !Local variables.
    integer :: i, j, new_size, old_size

    type(int_pair), dimension(:), allocatable :: copy_list

    if (.not.allocated(list)) then
      print*, 'An attempt was made to filter an unallocated list.'
      return
    end if

    !We transfer all elements from list to copy_list.
    old_size = size(list)
    allocate(copy_list(old_size))
    copy_list = list

    !We compute the new size of the list based on the number of true elements
    !in the boolean array called condition.
    new_size = count(condition)
    deallocate(list)
    allocate(list(new_size))

    !If the new size is greater than zero, we pass all the elements that meet
    !the condition to our list.
    if (new_size > 0) then
      j = 1
      do i = 1, old_size
        if (condition(i)) then
          list(j) = copy_list(i)
          j = j + 1
        end if
      end do
    end if
  end subroutine my_pack_int_pair


  !Recursive subroutine that sorts two lists in ascending order from the values
  !of the first.
  recursive subroutine quicksort(energies, node_ids)
    implicit none

    !Input/output arguments.
    double precision, dimension(:), intent(inout) :: energies

    integer, dimension(:), intent(inout) :: node_ids

    !Local variables.
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

      !If the list containing the elements that are greater (smaller) than the
      !pivot has not been allocated, we allocate it to avoid memory errors.
      if (.not.allocated(e_lt)) allocate(e_lt(0), n_lt(0))
      if (.not.allocated(e_gt)) allocate(e_gt(0), n_gt(0))

      !We apply the subroutine recursively to the lists that contain the elements
      !that are smaller and greater than the pivot, respectively.
      call quicksort(e_lt, n_lt)
      call quicksort(e_gt, n_gt)

      !Finally, we construct the sorted lists.
      energies = [e_lt, e_eq, e_gt]
      node_ids = [n_lt, n_eq, n_gt]
    end if
  end subroutine quicksort
end module array_procedures
