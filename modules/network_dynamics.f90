!Module whose procedures modify the connections of a network, creating a history of the changes.
module network_dynamics
  use array_procedures, only : add_item, del_int, find_int, int_list, int_list_list
  use network_generation, only : node
  use network_properties, only : mean_degree

  implicit none

  private

  public rewiring

contains
  !Subroutine that records the changes that a network undergoes over time. It also keeps
  !track of the connections active at a certain time step.
  subroutine rewiring(type, network, history, indices, q)
    implicit none

    character(len=*), intent(in) :: type
    double precision, intent(in) :: q
    type(node), dimension(:), intent(inout) :: history, network
    type(int_list_list), dimension(:), intent(inout) :: indices

    integer :: altk, c, counter=0, i, id_ki, isize, k, N
    type(int_list) :: locations

    N = size(network)
    c = nint(mean_degree(network))

    !If this is the first time we have called this subroutine, we initialize history.
    if (counter == 0) then
      do i = 1, N
        allocate(history(i)%neighbors(0), history(i)%opposites(0))
      end do
    !Otherwise, if the rewiring probability (q) is greater than 0, we modify the
    !network connections.
    else if ((counter > 0).and.(q > 0.0d0)) then
      select case (type)
        case default
          call std_rewiring(network, N, c, q)
      end select
    end if

    !We update history and indices.
    do i = 1, N
      allocate(locations%array(0))
      isize = size(network(i)%neighbors)
      do altk = 1, isize
        k = network(i)%neighbors(altk)

        id_ki = find_int(history(i)%neighbors, k)
        !If the connection (i, k) has not occurred until this time step,
        !we save it in history.
        if (id_ki == 0) then
          call add_item(history(i)%neighbors, k)
          call add_item(history(k)%opposites, size(history(i)%neighbors))
          call add_item(history(k)%neighbors, i)
          call add_item(history(i)%opposites, size(history(k)%neighbors))

          call add_item(locations%array, size(history(i)%neighbors))
        else
          call add_item(locations%array, id_ki)
        end if
      end do

      !We save the active connections that node i has at this time step inside indices.
      call add_item(indices(i)%time, locations)
      deallocate(locations%array)
    end do

    counter = counter + 1
  end subroutine rewiring
end module network_dynamics