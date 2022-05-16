!Module that contains procedures that simulate the spread of epidemics on
!dynamic networks using Monte Carlo algorithms that impose different
!epidemiological models (SIR).
!Author: Adrià Meca Montserrat.
!Last modified date: 16/05/22.
module mc_simulations
  use array_procedures, only : int_list_list
  use network_generation, only : node
  use random_number_generator, only : r1279

  implicit none

  private

  public mc_sir

contains
  !Subroutine that spreads an infection that follows the SIR scheme through
  !a dynamic network using a Monte Carlo simulation.
  subroutine mc_sir(history, indices, origin, lambda, mu, t0, states)
    implicit none

    character(len=1), dimension(:), intent(out) :: states
    double precision, intent(in) :: lambda, mu
    integer, intent(in) :: origin, t0
    type(int_list_list), dimension(:), intent(in) :: indices
    type(node), dimension(:), intent(in) :: history

    character(len=1), dimension(size(history)) :: loc_states
    double precision :: p, r
    integer, dimension(:), allocatable :: nbrs
    integer :: i, isize, N, t

    N = size(history)

    !We initialize the state arrays.
    states = 'S'
    states(origin) = 'I'
    loc_states = states

    !Monte Carlo simulation for the SIR model.
    do t = 1, t0
      do i = 1, N
        !We compute the neighbors of node i at time t.
        isize = size(indices(i)%time(t)%array)
        allocate(nbrs(isize))
        nbrs = history(i)%neighbors(indices(i)%time(t)%array)

        !We apply the SIR rules at each t: a node I becomes R with probability mu
        !regardless of its neighbors; a node S becomes I with probability
        !lambda when it makes contact with one of its neighbors.
        r = r1279()
        if ((loc_states(i) == 'I').and.(r < mu)) then
          states(i) = 'R'
        else if (loc_states(i) == 'S') then
          p = 1.0d0 - product(1.0d0-lambda*merge(1.0d0, 0.0d0, loc_states(nbrs)=='I'))
          if (r < p) states(i) = 'I'
        end if

        deallocate(nbrs)
      end do
      !We update loc_states at the end of each time step.
      loc_states = states
    end do
  end subroutine mc_sir
end module mc_simulations
