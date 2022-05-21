!Module whose procedures study the properties of networks.
!Author: Adrià Meca Montserrat.
!Last modified date: 21/05/22.
module network_properties
  use network_generation, only : node

  implicit none

  private

  public edges, mean_degree

contains
  !Function that calculates the total number of connection points of a network.
  !In other words, it computes the sum over the degree of each node in a graph.
  function connection_points(network)
    implicit none

    !Input arguments.
    type(node), dimension(:), intent(in) :: network

    !Output arguments.
    integer :: connection_points

    !Local variables.
    integer :: i

    connection_points = 0
    do i = 1, size(network)
      connection_points = connection_points + size(network(i)%neighbors)
    end do
  end function connection_points


  !Function that calculates the number of edges of a network.
  function edges(network)
    implicit none

    !Input arguments.
    type(node), dimension(:), intent(in) :: network

    !Output arguments.
    integer :: edges

    edges = connection_points(network) / 2
  end function edges


  !Function that calculates the mean degree of a network.
  function mean_degree(network)
    implicit none

    !Input arguments.
    type(node), dimension(:), intent(in) :: network

    !Output arguments.
    double precision :: mean_degree

    mean_degree = connection_points(network) / dble(size(network))
  end function mean_degree
end module network_properties
