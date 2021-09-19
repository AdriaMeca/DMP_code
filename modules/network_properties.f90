!Module whose procedures study the properties of networks.
module network_properties
  use network_generation, only : node

  implicit none

  private

  public mean_degree

contains
  !Function that calculates the mean degree of a network.
  double precision function mean_degree(network)
    implicit none

    type(node), dimension(:), intent(in) :: network

    integer :: i, N

    N = size(network)
    mean_degree = sum((/ (size(network(i)%neighbors), i=1,N) /)) / dble(N)
  end function mean_degree
end module network_properties