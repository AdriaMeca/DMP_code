!Module whose procedures create different networks.
!Author: Adrià Meca Montserrat.
!Last modified date: 14/05/22.
module network_generation
  use array_procedures, only : add_item, int_list, int_pair

  implicit none

  private

  type, public :: node
    character(len=1) :: state
    integer, dimension(:), allocatable :: neighbors, opposites
  end type node

  public proximity_network, random_regular_graph

contains
  !Function that creates a proximity network, i.e., a network in which
  !short connections are favored. The parameter l0 controls this behavior:
  !large values of l0 allow the presence of longer links, increasing the density.
  function proximity_network(N, c, l0) result(network)
    implicit none

    double precision, intent(in) :: l0
    integer, intent(in) :: c, N
    type(node), dimension(N) :: network

    double precision, dimension(N, 2) :: positions
    double precision, dimension(2) :: ri
    double precision :: dij, L, pij, random

    integer :: i, j

    L = sqrt(dble(N))
    pij = c / dble(N)

    !We initialize the positions of the nodes.
    do i = 1, N
      allocate(network(i)%neighbors(0))

      call random_number(ri)
      positions(i, :) = ri * L
    end do

    !We connect the nodes.
    do i = 1, N
      do j = i+1, N
        !Distance between nodes i and j.
        dij = sqrt(sum((positions(i, :)-positions(j, :))**2))

        !We decide if i and j are connected or not.
        call random_number(random)
        if (random < pij*exp(-dij/l0)) then
          call add_item(network(i)%neighbors, j)
          call add_item(network(j)%neighbors, i)
        end if
      end do
    end do
  end function proximity_network


  !Function that creates a random regular graph (RRG), i.e., a network
  !whose nodes have the same number of neighbors c (degree).
  function random_regular_graph(N, c) result(network)
    implicit none

    integer, intent(in) :: c, N
    type(node), dimension(N) :: network

    double precision :: r(2), x
    integer, dimension(:), allocatable :: array_u, array_v
    integer :: i, index_i, index_j, index_k, j, u, v

    type(int_list), dimension(N) :: G
    type(int_pair), dimension(:), allocatable :: array_w
    type(int_pair) :: link

    do while (.true.)
      allocate(array_u(N*c), array_v(N), array_w(0))

      !Step 1.
      do index_i = 1, N
        allocate(G(index_i)%array(c), network(index_i)%neighbors(0))
        do index_j = 1, c
          index_k = index_j + (index_i-1)*c

          !We define N groups G {G(1), ..., G(N)} that contain c points.
          G(index_i)%array(index_j) = index_k

          !We keep the N*c points (where N*c is an even number) in array_u.
          array_u(index_k) = index_k
        end do
        !We save the available groups in array_v.
        array_v(index_i) = index_i
      end do

      !Step 2.
      do while (size(array_v) > c)
        do while (.true.)
          call random_number(r)

          !We choose two points i and j at random.
          i = array_u(1 + floor(size(array_u)*r(1)))
          j = array_u(1 + floor(size(array_u)*r(2)))

          !i and j belong to groups u and v, respectively.
          u = ceiling(i/real(c))
          v = ceiling(j/real(c))

          !If u and v are different and disconnected, we accept i and j.
          if ((u /= v).and.(all(network(u)%neighbors /= v))) exit
        end do

        !We remove i and j from array_u.
        array_u = pack(array_u, mask=(array_u/=i).and.(array_u/=j))

        !We remove i and j from their groups.
        G(u)%array = pack(G(u)%array, mask=G(u)%array/=i)
        G(v)%array = pack(G(v)%array, mask=G(v)%array/=j)

        !u and v are removed from array_v if they become empty.
        if (size(G(u)%array) == 0) array_v = pack(array_v, mask=array_v/=u)
        if (size(G(v)%array) == 0) array_v = pack(array_v, mask=array_v/=v)

        !We connect u and v.
        call add_item(network(u)%neighbors, v)
        call add_item(network(v)%neighbors, u)
      end do

      !We save the remaining possible links {u, v} inside array_w.
      do index_i = 1, size(array_v)
        do index_j = index_i+1, size(array_v)
          u = array_v(index_i)
          v = array_v(index_j)

          if (all(network(u)%neighbors /= v)) call add_item(array_w, int_pair(u, v))
        end do
      end do

      !Step 3.
      do while (size(array_w) > 0)
        do while (.true.)
          call random_number(r)

          !We choose a link at random from array_w.
          link = array_w(1 + floor(size(array_w)*r(1)))

          u = link%x
          v = link%y

          !We accept u and v considering their size.
          x = dble(size(G(u)%array) * size(G(v)%array))
          if (r(2) < x/c**2) exit
        end do

        !We remove {u, v} from array_w.
        array_w = pack(array_w, mask=(array_w%x/=u).or.(array_w%y/=v))

        call random_number(r)

        !We choose two connections at random from u and v.
        i = G(u)%array(1 + floor(size(G(u)%array)*r(1)))
        j = G(v)%array(1 + floor(size(G(v)%array)*r(2)))

        !We remove i and j from their groups.
        G(u)%array = pack(G(u)%array, mask=G(u)%array/=i)
        G(v)%array = pack(G(v)%array, mask=G(v)%array/=j)

        !u and v are removed from array_v if they become empty. We also
        !remove links that contain u or v from array_w.
        if (size(G(u)%array) == 0) then
          array_v = pack(array_v, mask=array_v/=u)
          array_w = pack(array_w, mask=(array_w%x/=u).and.(array_w%y/=u))
        end if
        if (size(G(v)%array) == 0) then
          array_v = pack(array_v, mask=array_v/=v)
          array_w = pack(array_w, mask=(array_w%x/=v).and.(array_w%y/=v))
        end if

        !We connect u and v.
        call add_item(network(u)%neighbors, v)
        call add_item(network(v)%neighbors, u)
      end do

      !If there are no groups available, we are done.
      if (size(array_v) == 0) exit

      deallocate(array_u, array_v, array_w)
      do index_i = 1, N
        deallocate(G(index_i)%array, network(index_i)%neighbors)
      end do
    end do
  end function random_regular_graph
end module network_generation
