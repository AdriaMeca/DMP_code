!Module whose procedures create different networks.
!Author: Adria Meca Montserrat.
!Last modified date: 29/05/22.
module network_generation
  use array_procedures, only : add, int_list, int_pair, my_pack
  use random_number_generator, only : r1279

  implicit none

  private

  type, public :: node
    integer, dimension(:), allocatable :: neighbors, opposites
  end type node

  public PN, RRG

contains
  !Function that creates a Proximity Network (PN), i.e., a network in which
  !links are established based on the distance between nodes. A parameter l
  !controls the magnitude of the distances: larger values of l favor longer
  !connections. Here, l is implicitly present in the table of exponentials,
  !eij(idx) = exp(-dij(idx)/l), which is used to compute the probability of
  !connecting nodes i and j separated by a distance dij(idx).
  function PN(N, c, eij)
    implicit none

    !Input arguments.
    double precision, dimension(:), intent(in) :: eij

    integer, intent(in) :: c, N

    !Output arguments.
    type(node), dimension(N) :: PN

    !Local variables.
    double precision :: A, pij

    integer :: i, idx, j


    !We initialize the list of neighbors.
    do i = 1, N
      allocate(PN(i)%neighbors(0))
    end do

    !We calculate the normalization constant A.
    idx = 0
    A = 0.0d0
    do i = 1, N
      do j = i+1, N
        A = A + eij(idx+j-1)
      end do
      idx = idx + N - i - 1
    end do

    !Constant part of the probability of connecting any pair of nodes.
    pij = c * N / 2.0d0 / A

    !We connect the nodes.
    idx = 0
    do i = 1, N
      do j = i+1, N
        !We decide if i and j are connected or not.
        if (r1279() < pij*eij(idx+j-1)) then
          call add(PN(i)%neighbors, j)
          call add(PN(j)%neighbors, i)
        end if
      end do
      idx = idx + N - i - 1
    end do
  end function PN



  !Function that creates a Random Regular Graph (RRG), i.e., a network
  !whose nodes have the same number of neighbors c (degree).
  function RRG(N, c)
    implicit none

    !Input arguments.
    integer, intent(in) :: c, N

    !Output arguments.
    type(node), dimension(N) :: RRG

    !Local variables.
    double precision :: x

    integer, dimension(:), allocatable :: array_u, array_v
    integer :: i, idx_i, idx_j, idx_k, j, u, v

    type(int_list), dimension(N) :: G

    type(int_pair), dimension(:), allocatable :: array_w
    type(int_pair) :: link


    do while (.true.)
      allocate(array_u(N*c), array_v(N), array_w(0))

      !Step 1.
      do idx_i = 1, N
        allocate(G(idx_i)%array(c), RRG(idx_i)%neighbors(0))
        do idx_j = 1, c
          idx_k = idx_j + (idx_i-1)*c

          !We define N groups G {G(1), ..., G(N)} that contain c connectors.
          G(idx_i)%array(idx_j) = idx_k

          !We keep the N*c connectors (where N*c is an even number) in array_u.
          array_u(idx_k) = idx_k
        end do
        !We save the available groups in array_v.
        array_v(idx_i) = idx_i
      end do

      !Step 2.
      do while (size(array_v) > c)
        do while (.true.)
          !We choose two connectors i and j at random.
          i = array_u(1 + floor(size(array_u)*r1279()))
          j = array_u(1 + floor(size(array_u)*r1279()))

          !i and j belong to groups u and v, respectively.
          u = ceiling(i/real(c))
          v = ceiling(j/real(c))

          !If u and v are different and disconnected, we accept i and j.
          if ((u /= v).and.(all(RRG(u)%neighbors /= v))) exit
        end do

        !We remove i and j from array_u.
        call my_pack(array_u, (array_u/=i).and.(array_u/=j))

        !We remove i and j from their groups.
        call my_pack(G(u)%array, G(u)%array/=i)
        call my_pack(G(v)%array, G(v)%array/=j)

        !u and v are removed from array_v if they become empty.
        if (size(G(u)%array) == 0) call my_pack(array_v, array_v/=u)
        if (size(G(v)%array) == 0) call my_pack(array_v, array_v/=v)

        !We connect u and v.
        call add(RRG(u)%neighbors, v)
        call add(RRG(v)%neighbors, u)
      end do

      !We free array_u from memory because it is no longer needed.
      deallocate(array_u)

      !We save the remaining possible links {u, v} inside array_w.
      do idx_i = 1, size(array_v)
        do idx_j = idx_i+1, size(array_v)
          u = array_v(idx_i)
          v = array_v(idx_j)

          if (all(RRG(u)%neighbors /= v)) call add(array_w, int_pair(u, v))
        end do
      end do

      !Step 3.
      do while (size(array_w) > 0)
        do while (.true.)
          !We choose a link at random from array_w.
          link = array_w(1 + floor(size(array_w)*r1279()))

          u = link%x
          v = link%y

          !We accept u and v considering their size.
          x = dble(size(G(u)%array) * size(G(v)%array))
          if (r1279() < x/c**2) exit
        end do

        !We remove {u, v} from array_w.
        call my_pack(array_w, (array_w%x/=u).or.(array_w%y/=v))

        !We choose two connectors at random from u and v.
        i = G(u)%array(1 + floor(size(G(u)%array)*r1279()))
        j = G(v)%array(1 + floor(size(G(v)%array)*r1279()))

        !We remove i and j from their groups.
        call my_pack(G(u)%array, G(u)%array/=i)
        call my_pack(G(v)%array, G(v)%array/=j)

        !u and v are removed from array_v if they become empty. We also
        !remove links that contain u or v from array_w.
        if (size(G(u)%array) == 0) then
          call my_pack(array_v, array_v/=u)
          call my_pack(array_w, (array_w%x/=u).and.(array_w%y/=u))
        end if
        if (size(G(v)%array) == 0) then
          call my_pack(array_v, array_v/=v)
          call my_pack(array_w, (array_w%x/=v).and.(array_w%y/=v))
        end if

        !We connect u and v.
        call add(RRG(u)%neighbors, v)
        call add(RRG(v)%neighbors, u)
      end do

      !If there are no groups available, we are done.
      if (size(array_v) == 0) exit

      !If the attempt fails, we have to deallocate these arrays and start over.
      deallocate(array_v, array_w)
      do idx_i = 1, N
        deallocate(G(idx_i)%array, RRG(idx_i)%neighbors)
      end do
    end do
  end function RRG
end module network_generation
