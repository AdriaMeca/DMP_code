!Module that contains procedures that simulate the spread of epidemics on
!dynamic networks using Monte Carlo algorithms that impose different
!epidemiological models (SIR and SEIR).
!Author: Adrià Meca Montserrat.
!Last modified date: 19/05/22.
module mc_simulations
  use array_procedures, only : int_list_list
  use network_generation, only : node
  use random_number_generator, only : r1279

  implicit none

  private

  public mc_seir, mc_sir

contains
  !Subroutine that spreads an infection that follows the SEIR rules on a network
  !using a Monte Carlo simulation.
  subroutine mc_seir(history, indices, origin, alpha, lambda, mu, nu, t0, &
    states, realizations, tmp_mc_probs)
    implicit none

    !Input arguments.
    double precision, intent(in) :: alpha, lambda, mu, nu
    integer, intent(in) :: origin, realizations, t0
    type(int_list_list), dimension(:), intent(in) :: indices
    type(node), dimension(:), intent(in) :: history

    !Output arguments.
    character(len=1), dimension(:), intent(out) :: states
    double precision, dimension(:, :), intent(out) :: tmp_mc_probs

    !Local variables.
    double precision, dimension(t0, 4) :: mc_probs
    integer :: alti, i, N, t

    !Number of nodes.
    N = size(history)

    tmp_mc_probs = 0.0d0
    do alti = 1, realizations
      !We initialize the array of states.
      states = 'S'
      states(origin) = 'E'

      mc_probs = 0.0d0
      do t = 1, t0
        !We apply the MC kernel to modify the states using the SEIR rules.
        call mc_seir_kernel(history, indices, alpha, lambda, mu, nu, t, states)

        !We count the number of nodes that are in each state at time t.
        do i = 1, N
          select case (states(i))
            case ('S')
              mc_probs(t, 1) = mc_probs(t, 1) + 1.0d0
            case ('E')
              mc_probs(t, 2) = mc_probs(t, 2) + 1.0d0
            case ('I')
              mc_probs(t, 3) = mc_probs(t, 3) + 1.0d0
            case ('R')
              mc_probs(t, 4) = mc_probs(t, 4) + 1.0d0
          end select
        end do
      end do
      !We divide the above quantities by N, thus obtaining the time evolution
      !of the fraction of nodes that are in state X (i.e., S, E, I or R); we
      !compute this trajectory over different realizations and add them all.
      tmp_mc_probs = tmp_mc_probs + mc_probs/N
    end do
    !We divide the cumulative trajectory by the number of realizations, thus
    !obtaining the average trajectory.
    tmp_mc_probs = tmp_mc_probs / realizations
  end subroutine mc_seir


  !Kernel of the Monte Carlo simulation that applies the SEIR rules at each t.
  subroutine mc_seir_kernel(history, indices, alpha, lambda, mu, nu, t, states)
    implicit none

    !Input arguments.
    double precision, intent(in) :: alpha, lambda, mu, nu
    integer, intent(in) :: t
    type(int_list_list), dimension(:), intent(in) :: indices
    type(node), dimension(:), intent(in) :: history

    !Output arguments.
    character(len=1), dimension(:), intent(inout) :: states

    !Local variables.
    character(len=1), dimension(size(history)) :: loc_states
    double precision :: p1, p2, r
    integer, dimension(:), allocatable :: nbrs
    integer :: i, isize, N

    !Number of nodes.
    N = size(history)

    !We initialize an auxiliary array of states.
    loc_states = states
    do i = 1, N
      !We compute the neighbors of node i at time t.
      isize = size(indices(i)%time(t)%array)
      allocate(nbrs(isize))
      nbrs = history(i)%neighbors(indices(i)%time(t)%array)

      !We apply the SEIR rules at each t: an E (exposed) node becomes I with
      !probability nu, and an I node becomes R with probability mu, regardless
      !of their respective neighbors; a S node becomes E with probability alpha
      !when it makes contact with one of its E neighbors (modified SEIR); also,
      !it becomes E with probability lambda when it interacts with one of its I
      !neighbors.
      r = r1279()
      if ((loc_states(i) == 'E').and.(r < nu)) then
        states(i) = 'I'
      else if ((loc_states(i) == 'I').and.(r < mu)) then
        states(i) = 'R'
      else if (loc_states(i) == 'S') then
        p1 = 1.0d0 - product(1.0d0-alpha*merge(1.0d0, 0.0d0, loc_states(nbrs)=='E'))
        p2 = 1.0d0 - product(1.0d0-lambda*merge(1.0d0, 0.0d0, loc_states(nbrs)=='I'))
        if (r < p1+p2) states(i) = 'E'
      end if

      deallocate(nbrs)
    end do
  end subroutine mc_seir_kernel


  !Subroutine that spreads an infection that follows the SIR rules on a network
  !using a Monte Carlo simulation.
  subroutine mc_sir(history, indices, origin, lambda, mu, t0, states)
    implicit none

    !Input arguments.
    double precision, intent(in) :: lambda, mu
    integer, intent(in) :: origin, t0
    type(int_list_list), dimension(:), intent(in) :: indices
    type(node), dimension(:), intent(in) :: history

    !Output arguments.
    character(len=1), dimension(:), intent(out) :: states

    !Local variables.
    character(len=1), dimension(size(history)) :: loc_states
    double precision :: p, r
    integer, dimension(:), allocatable :: nbrs
    integer :: i, isize, N, t

    !Number of nodes.
    N = size(history)

    !We initialize the array of states.
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

        !We apply the SIR rules at each t: an I node becomes R with probability
        !mu regardless of its neighbors; a S node becomes I with probability
        !lambda when it makes contact with one of its I neighbors.
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
