!> Procedures for studying the properties of networks.
!> Author: Adria Meca Montserrat.
!> Last modified date: 12/08/22.
module network_properties
  use derived_types, only: node

  implicit none

  private

  public :: average_degree, draw_network, internode_distance, number_of_edges

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


  !> Computes the distance between two nodes.
  function internode_distance(r, i, j, pbc)
    integer,          intent(in) :: i, j                !> Nodes of interest.
    integer                      :: N                   !>
    double precision, intent(in) :: r(:, :)             !> Node positions.
    double precision             :: internode_distance  !> Output.
    double precision             :: L, xij, yij         !>
    logical,          intent(in) :: pbc                 !> (F) T: (do not) use PBC to compute the distance.

    xij = r(j, 1) - r(i, 1)
    yij = r(j, 2) - r(i, 2)

    if (pbc) then
      !> Number of nodes.
      N = size(r, dim=1)

      !> Side of the two-dimensional box where the nodes are confined.
      L = sqrt(dble(N))

      !> Periodic Boundary Conditions (PBC).
      xij = min(abs(xij), L-abs(xij))
      yij = min(abs(yij), L-abs(yij))
    end if

    !> Distance between nodes i and j.
    internode_distance = sqrt(xij*xij + yij*yij)
  end function internode_distance


  !> Calculates the number of edges in a network.
  function number_of_edges(network)
    integer                :: number_of_edges  !> Output.
    type(node), intent(in) :: network(:)       !>

    number_of_edges = connection_points(network) / 2
  end function number_of_edges


  !> Returns the information that Gnuplot needs to draw a network.
  subroutine draw_network(network, r)
    integer                      :: i, k, ki, N  !>
    double precision, intent(in) :: r(:, :)      !> Node positions.
    type(node),       intent(in) :: network(:)   !>

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
      do ki = 1, size(network(i)%neighbors)
        k = network(i)%neighbors(ki)
        write(*, '(4es26.16)') r(i, :), r(k, :) - r(i, :)
      end do
    end do
    write(*, '(a)') achar(10)
  end subroutine draw_network
end module network_properties
