!> Procedures that generate networks of different types.
!> Author: Adria Meca Montserrat.
!> Last modified date: 01/10/22.
!> Last reviewed date: 01/10/22.
module network_generation
  use array_procedures,        only: add, my_pack
  use derived_types,           only: node
  use network_properties,      only: internode_distance
  use random_number_generator, only: r1279

  implicit none

  private

  public :: PN, RRG

contains
  !> Generates a Proximity Network (PN) whose links are established based on the
  !> distance between nodes.
  function PN(N, c, r, l)
    integer,          intent(in) :: c            !> Parameter that controls connectivity.
    integer,          intent(in) :: N            !> Number of nodes.
    integer                      :: i, j         !>
    double precision, intent(in) :: l            !> Larger values of l allow longer link formation.
    double precision, intent(in) :: r(:, :)      !> Node positions.
    double precision             :: A, dij, pij  !>
    type(node)                   :: PN(N)        !> Output.

    !> We calculate the normalization constant A.
    A = 0.0d0
    do i = 1, N
      allocate(PN(i)%neighbors(0))

      do j = i+1, N
        !> Distance between nodes i and j.
        dij = internode_distance(r, i, j, pbc=.true.)

        !> We update A.
        A = A + exp(-dij/l)
      end do
    end do

    !> Constant part of the probability of connecting any pair of nodes.
    pij = N * c / 2.0d0 / A

    do i = 1, N
      do j = i+1, N
        !> Distance between nodes i and j.
        dij = internode_distance(r, i, j, pbc=.true.)

        !> We decide whether i and j are connected or not.
        if (r1279() < pij*exp(-dij/l)) then
          call add(PN(i)%neighbors, j)
          call add(PN(j)%neighbors, i)
        end if
      end do
    end do
  end function PN


  !> Generates a Random Regular Graph (RRG) whose nodes have the same degree.
  function RRG(N, c)
    integer,                 intent(in) :: c                                              !> Average degree.
    integer,                 intent(in) :: N                                              !> Number of nodes.
    integer,    allocatable             :: array_u(:), array_v(:)                         !>
    integer                             :: i, idx_i, idx_j, idx_k, j, total, u, usize, v  !>
    type(node)                          :: RRG(N)                                         !> Output.

    !> Warning: N * c must be an even number, otherwise the network cannot be closed
    !> and we are stuck in an infinite loop.
    do while (.true.)
      if (allocated(array_u)) deallocate(array_u); allocate(array_u(N*c))
      if (allocated(array_v)) deallocate(array_v); allocate(array_v(N*c))

      !> Step 1: we initialize the network and the arrays.
      do idx_i = 1, N
        if (allocated(RRG(idx_i)%neighbors)) deallocate(RRG(idx_i)%neighbors)
        allocate(RRG(idx_i)%neighbors(0))

        do idx_j = 1, c
          idx_k = idx_j + (idx_i-1)*c

          !> We store all connectors in 'array_u'.
          array_u(idx_k) = idx_k

          !> Group to which connector 'idx_k' belongs.
          array_v(idx_k) = idx_i
        end do
      end do

      !> Step 2: we establish the connections.
      total = 0
      usize = N * c
      do while (usize > 0)
        !> We choose two connectors i and j uniformly at random.
        i = array_u(1 + mod(int(usize*r1279()), usize))
        j = array_u(1 + mod(int(usize*r1279()), usize))

        !> i and j belong to groups u and v, respectively.
        u = array_v(i)
        v = array_v(j)

        !> If u and v are equal or connected, we reject i and j.
        if ((u == v) .or. any(RRG(u)%neighbors == v)) then
          total = total + 1

          !> If the counter reaches a certain value, we start over.
          if (total == usize * usize) exit

          !> Otherwise, we keep looking for a suitable pair (i, j).
          cycle
        end if

        !> If (i, j) is accepted, we reset the counter to zero.
        total = 0

        !> We remove i and j from 'array_u'.
        call my_pack(array_u, (array_u /= i) .and. (array_u /= j))

        !> We update the length of 'array_u'.
        usize = usize - 2

        !> We connect u and v.
        call add(RRG(u)%neighbors, v)
        call add(RRG(v)%neighbors, u)
      end do

      !> We exit the loop when the connectors run out. Otherwise, we start over.
      if (usize == 0) exit
    end do
  end function RRG
end module network_generation
