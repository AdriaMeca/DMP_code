!Module whose procedures study the properties of networks.
module network_properties
  use network_generation, only : node

  implicit none

  private

  public edges, mean_degree

contains
  !Function that calculates the mean degree of a network.
  double precision function mean_degree(network)
    implicit none

    type(node), dimension(:), intent(in) :: network

    integer :: N

    N = size(network)
    mean_degree = connection_points(network) / dble(N)
  end function mean_degree


  !Function that calculates the total number of connection points in a network.
  !In other words, it computes the sum over the degree of each node in a graph.
  integer function connection_points(network)
    implicit none

    type(node), dimension(:), intent(in) :: network

    integer :: i

    connection_points = 0
    do i = 1, size(network)
      connection_points = connection_points + size(network(i)%neighbors)
    end do
  end function connection_points


  !Function that calculates the number of edges of a network.
  integer function edges(network)
    implicit none

    type(node), dimension(:), intent(in) :: network

    edges = connection_points(network) / 2
  end function edges
end module network_properties
