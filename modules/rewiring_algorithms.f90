!Module whose procedures modify the connections of a network, creating a history
!of the changes.
!Author: Adria Meca Montserrat.
!Last modified date: 26/05/22.
module rewiring_algorithms
  use array_procedures, only : add, find, int_list, int_llist, my_pack
  use network_generation, only : node
  use random_number_generator, only : r1279

  implicit none

  private

  public rewiring

contains
  !Subroutine that records the changes that a network undergoes over time.
  !It also keeps track of the connections active at a certain time step.
  subroutine rewiring(type, network, history, indices, c, t0, Q)
    implicit none

    !Input arguments.
    character(len=*), intent(in) :: type

    double precision, intent(in) :: Q

    integer, intent(in) :: c, t0

    !Output arguments.
    type(node), dimension(:), intent(inout) :: history, network

    type(int_llist), dimension(:), intent(inout) :: indices

    !Local variables.
    integer :: altk, i, id_ki, isize, k, N, t

    type(int_list) :: locations


    !Number of nodes.
    N = size(network)

    !We initialize history.
    do i = 1, N
      allocate(history(i)%neighbors(0), history(i)%opposites(0))
    end do

    do t = 1, t0
      !If the rewiring probability Q is greater than 0.0, we modify the
      !network connections.
      if ((t > 1).and.(Q > 0.0d0)) then
        select case (type)
          case ('rrg')
            call std_rewiring(network, N, c, Q)
          case ('uni')
            call uni_rewiring(network, N, c, Q)
        end select
      end if

      !We update history and indices.
      do i = 1, N
        allocate(locations%array(0))
        isize = size(network(i)%neighbors)
        do altk = 1, isize
          k = network(i)%neighbors(altk)

          id_ki = find(history(i)%neighbors, k)
          !If the connection (i, k) has not occurred until this time step,
          !we save it in history.
          if (id_ki == 0) then
            call add(history(i)%neighbors, k)
            call add(history(k)%opposites, size(history(i)%neighbors))
            call add(history(k)%neighbors, i)
            call add(history(i)%opposites, size(history(k)%neighbors))

            call add(locations%array, size(history(i)%neighbors))
          else
            call add(locations%array, id_ki)
          end if
        end do

        !We save the active connections that node i has at this time step
        !inside indices.
        call add(indices(i)%time, locations)
        deallocate(locations%array)
      end do
    end do
  end subroutine rewiring



  !Subroutine that modifies the connections of a network using a standard rewiring
  !algorithm.
  subroutine std_rewiring(network, N, c, Q)
    implicit none

    !Input arguments.
    double precision, intent(in) :: Q

    integer, intent(in) :: c, N

    !Output arguments.
    type(node), dimension(:), intent(inout) :: network

    !Local variables.
    integer :: i, idx, j, k, m


    do idx = 1, N*c/4
      if (r1279() < Q) then
        do while (.true.)
          !We choose two random nodes from the network.
          i = 1 + floor(N*r1279())
          k = 1 + floor(N*r1279())

          !If i and k are different and disconnected, we proceed.
          if ((i /= k).and.(all(network(i)%neighbors /= k))) then
            !We choose two random neighbors of i and k, respectively.
            j = network(i)%neighbors(1 + floor(size(network(i)%neighbors)*r1279()))
            m = network(k)%neighbors(1 + floor(size(network(k)%neighbors)*r1279()))

            !If j and m are different and disconnected, we proceed.
            if ((j /= m).and.(all(network(j)%neighbors /= m))) exit
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



  !Subroutine that modifies the connections of a network using a uniform rewiring
  !algorithm.
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
          !We choose two random nodes from the network.
          i = 1 + floor(N*r1279())
          k = 1 + floor(N*r1279())

          !If i and k are different and disconnected, we proceed.
          if ((i /= k).and.(all(network(i)%neighbors /= k))) then
            !We choose a random node from the network.
            j = 1 + floor(N*r1279())

            !If j has neighbors we proceed.
            jsize = size(network(j)%neighbors)
            if (jsize > 0) then
              !We choose a random neighbor of j.
              m = network(j)%neighbors(1 + floor(jsize*r1279()))
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
