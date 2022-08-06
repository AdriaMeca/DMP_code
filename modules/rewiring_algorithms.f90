!> Procedures that rewire the links in a network over time.
!> Author: Adria Meca Montserrat.
!> Last modified date: 06/08/22.
module rewiring_algorithms
  use array_procedures, only: add, find, int_list, int_llist, my_pack
  use network_generation, only: node
  use random_number_generator, only: r1279

  implicit none

  private

  public :: rewiring

contains
  !> Records the link changes a network experiences over time and keeps track of
  !> its active connections at each time step.
  subroutine rewiring(model, network, history, indices, l, Q, r, c, t0)
    character(len=*), intent(in)    :: model                        !> Network type.
    double precision, intent(in)    :: l                            !> Parameter that controls locality.
    double precision, intent(in)    :: Q                            !> Rewiring probability.
    double precision, intent(in)    :: r(:, :)                      !> Node positions.
    integer,          intent(in)    :: c                            !> Parameter that controls connectivity.
    integer,          intent(in)    :: t0                           !> Observation time.
    integer                         :: altk, i, isize, k, ki, N, t  !>
    type(int_list)                  :: locations                    !> Active links of node i at time t.
    type(int_llist),  intent(inout) :: indices(:)                   !> Active links throughout the simulation.
    type(node),       intent(inout) :: history(:)                   !> Rewiring history.
    type(node),       intent(inout) :: network(:)                   !> Original network.

    !> Number of nodes.
    N = size(network)

    !> We initialize the rewiring history.
    do i = 1, N
      !> We remove the information from a hypothetical previous iteration.
      if (allocated(history(i)%neighbors)) deallocate(history(i)%neighbors)
      if (allocated(history(i)%opposites)) deallocate(history(i)%opposites)
      if (allocated(indices(i)%time)) deallocate(indices(i)%time)

      allocate(history(i)%neighbors(0), history(i)%opposites(0))
    end do

    !> We rewire links.
    do t = 0, t0
      if ((t > 0) .and. (Q > 0.0d0)) then
        select case (trim(model))
          case ('PN')
            call loc_rewiring(network, N, c, l, Q, r)
          case ('RRG')
            call std_rewiring(network, N, c, Q)
          case default
            call uni_rewiring(network, N, c, Q)
        end select
      end if

      !> We update 'history' and 'indices'.
      do i = 1, N
        allocate(locations%array(0))

        isize = size(network(i)%neighbors)
        do altk = 1, isize
          k = network(i)%neighbors(altk)

          !> We record every connection made throughout the simulation.
          ki = find(history(i)%neighbors, k)
          if (ki == 0) then
            call add(history(i)%neighbors, k)
            call add(history(k)%neighbors, i)
            call add(history(i)%opposites, size(history(k)%neighbors))
            call add(history(k)%opposites, size(history(i)%neighbors))

            call add(locations%array, size(history(i)%neighbors))
          else
            call add(locations%array, ki)
          end if
        end do

        !> We save the active links of node i at time t.
        call add(indices(i)%time, locations, 0)

        deallocate(locations%array)
      end do
    end do
  end subroutine rewiring


  !> Rewires links in a network based on their length to preserve locality.
  subroutine loc_rewiring(network, N, c, l, Q, r)
    double precision, intent(in)    :: l                                      !> Parameter that controls locality.
    double precision, intent(in)    :: Q                                      !> Rewiring probability.
    double precision, intent(in)    :: r(:, :)                                !> Node positions.
    double precision                :: A, cum, dij, pij, rng, side, xij, yij  !>
    integer,          intent(in)    :: c                                      !> Parameter that controls connectivity.
    integer,          intent(in)    :: N                                      !> Number of nodes.
    double precision                :: B(N)                                   !> Upper bounds of the towers.
    integer                         :: i, idx, j, k, ksize, m                 !>
    type(node),       intent(inout) :: network(:)                             !> Original network.

    !> Side of the two-dimensional box where the nodes are confined.
    side = sqrt(dble(N))

    !> We calculate the normalization constant A.
    A = 0.0d0
    do i = 1, N
      do j = i+1, N
        !> Periodic Boundary Conditions (PBC).
        xij = min(abs(r(j, 1)-r(i, 1)), side-abs(r(j, 1)-r(i, 1)))
        yij = min(abs(r(j, 2)-r(i, 2)), side-abs(r(j, 2)-r(i, 2)))

        !> Distance between nodes i and j.
        dij = sqrt(xij*xij + yij*yij)

        !> We update A.
        A = A + exp(-dij/l)
      end do
    end do

    !> Constant part of the probability of connecting any pair of nodes.
    pij = N * c / 2.0d0 / A

    !> We compute the upper bounds of the tower associated with each node.
    B = 0.0d0
    do i = 1, N
      do j = 1, N
        if (i == j) cycle

        !> We apply PBC and calculate the distance between nodes i and j.
        xij = min(abs(r(j, 1)-r(i, 1)), side-abs(r(j, 1)-r(i, 1)))
        yij = min(abs(r(j, 2)-r(i, 2)), side-abs(r(j, 2)-r(i, 2)))
        dij = sqrt(xij*xij + yij*yij)

        !> We update B(i).
        B(i) = B(i) + pij*exp(-dij/l)
      end do
    end do

    do idx = 1, N*c/4
      if (r1279() < Q) then
        do while (.true.)
          !> We choose a node i uniformly at random.
          i = 1 + mod(int(N*r1279()), N)

          !> We randomly choose a node j using the tower method.
          cum = 0.0d0
          rng = B(i) * r1279()
          do j = 1, N
            if (i == j) cycle

            !> We apply PBC and calculate the distance between nodes i and j.
            xij = min(abs(r(j, 1)-r(i, 1)), side-abs(r(j, 1)-r(i, 1)))
            yij = min(abs(r(j, 2)-r(i, 2)), side-abs(r(j, 2)-r(i, 2)))
            dij = sqrt(xij*xij + yij*yij)

            !> We update the cumulative probability.
            cum = cum + pij*exp(-dij/l)

            !> We have found j if 'rng' is less than the cumulative probability.
            if (rng < cum) exit
          end do
          !> We make sure that j is always less than N.
          j = merge(N, j, j>N)

          !> Temporary.
          if (i == j) print *, 'Ups!'

          !> If i and j are disconnected, we proceed.
          if (all(network(i)%neighbors /= j)) then
            !> We choose a node k uniformly at random.
            k = 1 + mod(int(N*r1279()), N)

            !> If k has at least one neighbor, we proceed.
            ksize = size(network(k)%neighbors)
            if (ksize > 0) then
              !> We choose a neighbor of k uniformly at random.
              m = network(k)%neighbors(1 + mod(int(ksize*r1279()), ksize))
              exit
            end if
          end if
        end do

        !> We remove the old connection.
        call my_pack(network(k)%neighbors, network(k)%neighbors /= m)
        call my_pack(network(m)%neighbors, network(m)%neighbors /= k)

        !> We add the new connection.
        call add(network(i)%neighbors, j)
        call add(network(j)%neighbors, i)
      end if
    end do
  end subroutine loc_rewiring


  !> Rewires the links in a network while preserving its degree distribution.
  subroutine std_rewiring(network, N, c, Q)
    double precision, intent(in)    :: Q                              !> Rewiring probability.
    integer,          intent(in)    :: c                              !> Parameter that controls connectivity.
    integer,          intent(in)    :: N                              !> Number of nodes.
    integer                         :: i, idx, isize, j, k, ksize, m  !>
    type(node),       intent(inout) :: network(:)                     !> Original network.

    do idx = 1, N*c/4
      if (r1279() < Q) then
        do while (.true.)
          !> We choose two nodes uniformly at random.
          i = 1 + mod(int(N*r1279()), N)
          k = 1 + mod(int(N*r1279()), N)

          !> We check if nodes i and k have at least one neighbor each.
          isize = size(network(i)%neighbors)
          ksize = size(network(k)%neighbors)
          if ((isize > 0) .and. (ksize > 0)) then
            !> If i and k are different and disconnected, we proceed.
            if ((i /= k) .and. all(network(i)%neighbors /= k)) then
              !> We choose a neighbor of each node uniformly at random.
              j = network(i)%neighbors(1 + mod(int(isize*r1279()), isize))
              m = network(k)%neighbors(1 + mod(int(ksize*r1279()), ksize))

              !> If j and m are different and disconnected, we proceed.
              if ((j /= m) .and. all(network(j)%neighbors /= m)) exit
            end if
          end if
        end do

        !> We remove the old connections.
        call my_pack(network(i)%neighbors, network(i)%neighbors /= j)
        call my_pack(network(j)%neighbors, network(j)%neighbors /= i)
        call my_pack(network(k)%neighbors, network(k)%neighbors /= m)
        call my_pack(network(m)%neighbors, network(m)%neighbors /= k)

        !> We add the new connections.
        call add(network(i)%neighbors, k)
        call add(network(k)%neighbors, i)
        call add(network(j)%neighbors, m)
        call add(network(m)%neighbors, j)
      end if
    end do
  end subroutine std_rewiring


  !> Uniformly rewires the links in a network.
  subroutine uni_rewiring(network, N, c, Q)
    double precision, intent(in)    :: Q                       !> Rewiring probability.
    integer,          intent(in)    :: c                       !> Parameter that controls connectivity.
    integer,          intent(in)    :: N                       !> Number of nodes.
    integer                         :: i, idx, j, jsize, k, m  !>
    type(node),       intent(inout) :: network(:)              !> Original network.

    do idx = 1, N*c/4
      if (r1279() < Q) then
        do while (.true.)
          !> We choose two nodes uniformly at random.
          i = 1 + mod(int(N*r1279()), N)
          k = 1 + mod(int(N*r1279()), N)

          !> If i and k are different and disconnected, we proceed.
          if ((i /= k) .and. all(network(i)%neighbors /= k)) then
            !> We choose a node uniformly at random.
            j = 1 + mod(int(N*r1279()), N)

            !> If j has at least one neighbor, we proceed.
            jsize = size(network(j)%neighbors)
            if (jsize > 0) then
              !> We choose a neighbor of j uniformly at random.
              m = network(j)%neighbors(1 + mod(int(jsize*r1279()), jsize))
              exit
            end if
          end if
        end do

        !> We remove the old connection.
        call my_pack(network(j)%neighbors, network(j)%neighbors /= m)
        call my_pack(network(m)%neighbors, network(m)%neighbors /= j)

        !> We add the new connection.
        call add(network(i)%neighbors, k)
        call add(network(k)%neighbors, i)
      end if
    end do
  end subroutine uni_rewiring
end module rewiring_algorithms
