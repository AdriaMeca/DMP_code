!Module whose procedures study the properties of networks.
!Author: Adria Meca Montserrat.
!Last modified date: 26/05/22.
module network_properties
  use network_generation, only : node

  implicit none

  private

  public average_degree, draw_network, number_of_edges

contains
  !Function that calculates the average degree of a network.
  function average_degree(network)
    implicit none

    !Input arguments.
    type(node), dimension(:), intent(in) :: network

    !Output arguments.
    double precision :: average_degree


    average_degree = connection_points(network) / dble(size(network))
  end function average_degree



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
  function number_of_edges(network)
    implicit none

    !Input arguments.
    type(node), dimension(:), intent(in) :: network

    !Output arguments.
    integer :: number_of_edges


    number_of_edges = connection_points(network) / 2
  end function number_of_edges



  !Subroutine that writes the information needed for Gnuplot to draw a network.
  subroutine draw_network(network, positions)
    implicit none

    !Input arguments.
    double precision, dimension(:, :), intent(in) :: positions

    type(node), dimension(:), intent(in) :: network

    !Local variables.
    integer :: i, isize, k, ki, N


    !Number of nodes.
    N = size(network)

    !First, we write the positions of the nodes.
    do i = 1, N
      write(*, '(2es26.16)') positions(i, :)
    end do
    write(*, '(a)') achar(10)

    !In order for Gnuplot to draw a vector between nodes i and k, we have to
    !write (xi, yi) and (xk-xi, yk-yi).
    do i = 1, N
      isize = size(network(i)%neighbors)
      do ki = 1, isize
        k = network(i)%neighbors(ki)
        write(*, '(4es26.16)') positions(i, :), positions(k, :)-positions(i, :)
      end do
    end do
    write(*, '(a)') achar(10)
  end subroutine draw_network
end module network_properties
