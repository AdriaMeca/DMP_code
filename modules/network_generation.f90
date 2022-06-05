!Module whose procedures create different networks.
!Author: Adria Meca Montserrat.
!Last modified date: 05/06/22.
module network_generation
  use array_procedures, only : add, my_pack
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
  !connections.
  function PN(N, c, r, l)
    implicit none

    !Input arguments.
    double precision, dimension(:, :), intent(in) :: r
    double precision, intent(in) :: l

    integer, intent(in) :: c, N

    !Output arguments.
    type(node), dimension(N) :: PN

    !Local variables.
    double precision :: A, dij, pij, side, xij, yij

    integer :: i, j


    !Side of the 2-dimensional box where the nodes are confined.
    side = sqrt(dble(N))

    !We calculate the normalization constant A.
    A = 0.0d0
    do i = 1, N
      !We initialize the list of neighbors.
      allocate(PN(i)%neighbors(0))
      do j = i+1, N
        !Periodic Boundary Conditions (PBC).
        xij = min(abs(r(j, 1)-r(i, 1)), side-abs(r(j, 1)-r(i, 1)))
        yij = min(abs(r(j, 2)-r(i, 2)), side-abs(r(j, 2)-r(i, 2)))

        !Distance between nodes i and j.
        dij = sqrt(xij*xij + yij*yij)

        !We increment the value of A.
        A = A + exp(-dij/l)
      end do
    end do

    !Constant part of the probability of connecting any pair of nodes.
    pij = c * N / 2.0d0 / A

    !We connect the nodes.
    do i = 1, N
      do j = i+1, N
        !PBC.
        xij = min(abs(r(j, 1)-r(i, 1)), side-abs(r(j, 1)-r(i, 1)))
        yij = min(abs(r(j, 2)-r(i, 2)), side-abs(r(j, 2)-r(i, 2)))

        !Distance between nodes i and j.
        dij = sqrt(xij*xij + yij*yij)

        !We decide if i and j are connected or not.
        if (r1279() < pij*exp(-dij/l)) then
          call add(PN(i)%neighbors, j)
          call add(PN(j)%neighbors, i)
        end if
      end do
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
    integer, dimension(:), allocatable :: array_u, array_v
    integer :: i, idx_i, idx_j, idx_k, j, total, u, usize, v


    !Warning: N * c has to be an even number, otherwise the network cannot be
    !closed and we are stuck in an infinite loop.
    do while (.true.)
      allocate(array_u(N*c), array_v(N*c))

      !Step 1: we initialize the network and the arrays.
      do idx_i = 1, N
        allocate(RRG(idx_i)%neighbors(0))
        do idx_j = 1, c
          idx_k = idx_j + (idx_i-1)*c

          !We store all connectors in array_u.
          array_u(idx_k) = idx_k

          !Group to which connector 'idx_k' belongs.
          array_v(idx_k) = idx_i
        end do
      end do

      !Step 2: we establish the connections.
      total = 0
      usize = N * c
      do while (usize > 0)
        !We choose two connectors i and j at random.
        i = array_u(1 + floor(usize*r1279()))
        j = array_u(1 + floor(usize*r1279()))

        !i and j belong to groups u and v, respectively.
        u = array_v(i)
        v = array_v(j)

        !If u and v are equal or connected, we reject i and j.
        if ((u == v).or.(any(RRG(u)%neighbors == v))) then
          total = total + 1

          !If the counter reaches a certain value, we start over.
          if (total == usize * usize) exit
          !Otherwise, we keep looking for a suitable pair (i, j).
          cycle
        end if

        !If (i, j) is accepted, we reset the counter to zero.
        total = 0

        !We remove i and j from array_u.
        call my_pack(array_u, (array_u/=i).and.(array_u/=j))

        !We update the length of array_u.
        usize = usize - 2

        !We connect u and v.
        call add(RRG(u)%neighbors, v)
        call add(RRG(v)%neighbors, u)
      end do

      !If there are no connectors available, we are done.
      if (usize == 0) exit

      !If the attempt fails, we have to deallocate everything and start over.
      deallocate(array_u, array_v)
      do idx_i = 1, N
        deallocate(RRG(idx_i)%neighbors)
      end do
    end do
  end function RRG
end module network_generation
