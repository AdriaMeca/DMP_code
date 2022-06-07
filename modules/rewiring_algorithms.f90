!Module whose procedures modify the connections of a network, creating a history
!of the changes.
!Author: Adria Meca Montserrat.
!Last modified date: 07/06/22.
module rewiring_algorithms
  use array_procedures, only : add, find, int_list, int_llist, my_pack
  use network_generation, only : node
  use random_number_generator, only : r1279

  implicit none

  private

  public rewiring

contains
  !Subroutine that registers the changes that a network undergoes over time. It
  !also keeps track of its active connections at a certain time step.
  subroutine rewiring(model, network, history, indices, c, t0, r, l, Q)
    implicit none

    !Input arguments.
    character(len=*), intent(in) :: model

    double precision, dimension(:, :), intent(in) :: r
    double precision, intent(in) :: l, Q

    integer, intent(in) :: c, t0

    !Output arguments.
    type(node), dimension(:), intent(inout) :: history, network

    type(int_llist), dimension(:), intent(inout) :: indices

    !Local variables.
    integer :: altk, i, isize, k, ki, N, t

    type(int_list) :: locations


    !Number of nodes.
    N = size(network)

    !We initialize history.
    do i = 1, N
      !We remove the information from a hypothetical previous iteration.
      if (allocated(history(i)%neighbors)) deallocate(history(i)%neighbors)
      if (allocated(history(i)%opposites)) deallocate(history(i)%opposites)
      if (allocated(indices(i)%time)) deallocate(indices(i)%time)

      allocate(history(i)%neighbors(0), history(i)%opposites(0))
    end do

    do t = 1, t0
      !We modify the network connections.
      if ((t > 1).and.(Q > 0.0d0)) then
        select case (trim(model))
          case ('PN')
            call loc_rewiring(network, N, c, r, l, Q)
          case ('RRG')
            call std_rewiring(network, N, c, Q)
          case default
            call uni_rewiring(network, N, c, Q)
        end select
      end if

      !We update history and indices.
      do i = 1, N
        allocate(locations%array(0))
        isize = size(network(i)%neighbors)
        do altk = 1, isize
          k = network(i)%neighbors(altk)

          ki = find(history(i)%neighbors, k)
          !We record every connection made throughout the process.
          if (ki == 0) then
            call add(history(i)%neighbors, k)
            call add(history(k)%opposites, size(history(i)%neighbors))
            call add(history(k)%neighbors, i)
            call add(history(i)%opposites, size(history(k)%neighbors))

            call add(locations%array, size(history(i)%neighbors))
          else
            call add(locations%array, ki)
          end if
        end do

        !We save the active connections that node i has at time t.
        call add(indices(i)%time, locations)
        deallocate(locations%array)
      end do
    end do
  end subroutine rewiring



  !Subroutine that rewires the connections of a network based on their length,
  !trying to preserve locality.
  subroutine loc_rewiring(network, N, c, r, l, Q)
    implicit none

    !Input arguments.
    double precision, dimension(:, :), intent(in) :: r
    double precision, intent(in) :: l, Q

    integer, intent(in) :: c, N

    !Output arguments.
    type(node), dimension(:), intent(inout) :: network

    !Local variables.
    double precision, dimension(N) :: B
    double precision :: A, cum, dij, pij, rng, side, xij, yij
    integer :: i, idx, j, k, ksize, m


    !Side of the two-dimensional box where the nodes are confined.
    side = sqrt(dble(N))

    !Calculation of the normalization constant A.
    A = 0.0d0
    do i = 1, N
      do j = i+1, N
        !Periodic Boundary Conditions (PBC).
        xij = min(abs(r(j, 1)-r(i, 1)), side-abs(r(j, 1)-r(i, 1)))
        yij = min(abs(r(j, 2)-r(i, 2)), side-abs(r(j, 2)-r(i, 2)))
        !Distance between nodes i and j.
        dij = sqrt(xij*xij + yij*yij)

        !We update the value of A.
        A = A + exp(-dij/l)
      end do
    end do

    !Constant part of the probability of connecting any pair of nodes.
    pij = N * c / 2.0d0 / A

    !Calculation of the upper value B of the 'tower' associated with each node.
    B = 0.0d0
    do i = 1, N
      do j = 1, N
        !We apply PBC and calculate the distance between i and j.
        xij = min(abs(r(j, 1)-r(i, 1)), side-abs(r(j, 1)-r(i, 1)))
        yij = min(abs(r(j, 2)-r(i, 2)), side-abs(r(j, 2)-r(i, 2)))
        dij = sqrt(xij*xij + yij*yij)

        !We update B(i) for each i.
        B(i) = B(i) + pij*merge(exp(-dij/l), 0.0d0, dij>0.0d0)
      end do
    end do

    do idx = 1, N*c/4
      if (r1279() < Q) then
        do while (.true.)
          !We choose a node i uniformly at random.
          i = 1 + mod(int(N*r1279()), N)

          !We choose a node j at random using the tower method.
          cum = 0.0d0
          rng = B(i) * r1279()
          do j = 1, N
            !We apply PBC and calculate the distance between i and j.
            xij = min(abs(r(j, 1)-r(i, 1)), side-abs(r(j, 1)-r(i, 1)))
            yij = min(abs(r(j, 2)-r(i, 2)), side-abs(r(j, 2)-r(i, 2)))
            dij = sqrt(xij*xij + yij*yij)

            cum = cum + pij*merge(exp(-dij/l), 0.0d0, dij>0.0d0)

            !If rng is less than the cumulative probability, we have found j.
            if (rng < cum) exit
          end do
          !In principle j will always be less than N, but just in case.
          j = merge(N, j, j>N)

          !If i and j are different and disconnected, we proceed.
          if ((i /= j).and.all(network(i)%neighbors /= j)) then
            !We choose a node k uniformly at random.
            k = 1 + mod(int(N*r1279()), N)

            !If k has neighbors, we proceed.
            ksize = size(network(k)%neighbors)
            if (ksize > 0) then
              !We choose a neighbor of k uniformly at random.
              m = network(k)%neighbors(1 + mod(int(ksize*r1279()), ksize))
              exit
            end if
          end if
        end do

        !We remove the old connection.
        call my_pack(network(k)%neighbors, network(k)%neighbors/=m)
        call my_pack(network(m)%neighbors, network(m)%neighbors/=k)

        !We add the new connection.
        call add(network(i)%neighbors, j)
        call add(network(j)%neighbors, i)
      end if
    end do
  end subroutine loc_rewiring



  !Subroutine that rewires the connections of a network without modifying the
  !degree distribution of its nodes.
  subroutine std_rewiring(network, N, c, Q)
    implicit none

    !Input arguments.
    double precision, intent(in) :: Q

    integer, intent(in) :: c, N

    !Output arguments.
    type(node), dimension(:), intent(inout) :: network

    !Local variables.
    integer :: i, idx, isize, j, k, ksize, m


    do idx = 1, N*c/4
      if (r1279() < Q) then
        do while (.true.)
          !We choose two nodes uniformly at random.
          i = 1 + mod(int(N*r1279()), N)
          k = 1 + mod(int(N*r1279()), N)

          !We need to check if the nodes i and k have at least one neighbor.
          isize = size(network(i)%neighbors)
          ksize = size(network(k)%neighbors)
          if ((isize > 0).and.(ksize > 0)) then
            !If i and k are different and disconnected, we proceed.
            if ((i /= k).and.all(network(i)%neighbors /= k)) then
              !We choose a neighbor of each node uniformly at random.
              j = network(i)%neighbors(1 + mod(int(isize*r1279()), isize))
              m = network(k)%neighbors(1 + mod(int(ksize*r1279()), ksize))

              !If j and m are different and disconnected, we proceed.
              if ((j /= m).and.all(network(j)%neighbors /= m)) exit
            end if
          end if
        end do

        !We remove the previous connections.
        call my_pack(network(i)%neighbors, network(i)%neighbors/=j)
        call my_pack(network(j)%neighbors, network(j)%neighbors/=i)
        call my_pack(network(k)%neighbors, network(k)%neighbors/=m)
        call my_pack(network(m)%neighbors, network(m)%neighbors/=k)

        !We add the new connections.
        call add(network(i)%neighbors, k)
        call add(network(k)%neighbors, i)
        call add(network(j)%neighbors, m)
        call add(network(m)%neighbors, j)
      end if
    end do
  end subroutine std_rewiring



  !Subroutine that uniformly rewires the connections of a network.
  subroutine uni_rewiring(network, N, c, Q)
    implicit none

    !Input arguments.
    double precision, intent(in) :: Q

    integer, intent(in) :: c, N

    !Output arguments.
    type(node), dimension(:), intent(inout) :: network

    !Local variables.
    integer :: i, idx, j, jsize, k, m


    do idx = 1, N*c/4
      if (r1279() < Q) then
        do while (.true.)
          !We choose two nodes uniformly at random.
          i = 1 + mod(int(N*r1279()), N)
          k = 1 + mod(int(N*r1279()), N)

          !If i and k are different and disconnected, we proceed.
          if ((i /= k).and.all(network(i)%neighbors /= k)) then
            !We choose a node uniformly at random.
            j = 1 + mod(int(N*r1279()), N)

            !If j has neighbors, we proceed.
            jsize = size(network(j)%neighbors)
            if (jsize > 0) then
              !We choose a neighbor of j uniformly at random.
              m = network(j)%neighbors(1 + mod(int(jsize*r1279()), jsize))
              exit
            end if
          end if
        end do

        !We remove the old connection.
        call my_pack(network(j)%neighbors, network(j)%neighbors/=m)
        call my_pack(network(m)%neighbors, network(m)%neighbors/=j)

        !We add the new connection.
        call add(network(i)%neighbors, k)
        call add(network(k)%neighbors, i)
      end if
    end do
  end subroutine uni_rewiring
end module rewiring_algorithms
