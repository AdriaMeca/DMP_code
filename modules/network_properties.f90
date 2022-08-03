!> Procedures for studying the properties of networks.
!> Author: Adria Meca Montserrat.
!> Last modified date: 03/08/22.
module network_properties
  use network_generation, only: node

  implicit none

  private

  public :: average_degree, draw_network, number_of_edges

contains
  !> Calculates the average degree of a network.
  function average_degree(network)
    double precision             :: average_degree  !> Output.
    type(node),       intent(in) :: network(:)      !>

    average_degree = connection_points(network) / dble(size(network))
  end function average_degree


  !> Calculates the sum of the degrees of all nodes in a network.
  function connection_points(network)
    integer                :: connection_points  !> Output.
    integer                :: i                  !>
    type(node), intent(in) :: network(:)         !>

    connection_points = 0
    do i = 1, size(network)
      connection_points = connection_points + size(network(i)%neighbors)
    end do
  end function connection_points


  !> Calculates the number of edges in a network.
  function number_of_edges(network)
    integer                :: number_of_edges  !> Output.
    type(node), intent(in) :: network(:)       !>

    number_of_edges = connection_points(network) / 2
  end function number_of_edges


  !> Returns the information that Gnuplot needs to draw a network.
  subroutine draw_network(network, r)
    double precision, intent(in) :: r(:, :)             !> Node positions.
    integer                      :: i, isize, k, ki, N  !>
    type(node),       intent(in) :: network(:)          !>

    !> Number of nodes.
    N = size(network)

    !> We write the positions of the nodes.
    do i = 1, N
      write(*, '(2es26.16)') r(i, :)
    end do
    write(*, '(a)') achar(10)

    !> Writing (xi, yi) and (xk-xi, yk-yi) allows Gnuplot to draw a vector between
    !> nodes i and k.
    do i = 1, N
      isize = size(network(i)%neighbors)
      do ki = 1, isize
        k = network(i)%neighbors(ki)
        write(*, '(4es26.16)') r(i, :), r(k, :) - r(i, :)
      end do
    end do
    write(*, '(a)') achar(10)
  end subroutine draw_network
end module network_properties
