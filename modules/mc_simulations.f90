!Module that contains procedures that simulate the spread of epidemics on
!dynamic networks using Monte Carlo algorithms that impose different
!epidemiological models (SIR and SEIR).
!Author: Adrià Meca Montserrat.
!Last modified date: 25/05/22.
module mc_simulations
  use array_procedures, only : int_llist
  use network_generation, only : node
  use random_number_generator, only : r1279

  implicit none

  private

  public mc_sim

contains
  !Subroutine that spreads an infection that follows the S(E)IR rules on a network
  !using a Monte Carlo simulation.
  subroutine mc_sim(type, history, indices, origins, alpha, lambda, mu, nu, t0, &
    states, realizations, tmp_mc_probs)
    implicit none

    !Input arguments.
    character(len=*), intent(in) :: type

    double precision, intent(in) :: alpha, lambda, mu, nu

    integer, dimension(:), intent(in) :: origins
    integer, intent(in) :: realizations, t0

    type(int_llist), dimension(:), intent(in) :: indices

    type(node), dimension(:), intent(in) :: history

    !Output arguments.
    character(len=1), dimension(:), intent(out) :: states

    double precision, dimension(:, :), intent(out) :: tmp_mc_probs

    !Local variables.
    character(len=1), dimension(size(history)) :: auxsts

    double precision, dimension(t0, 4) :: mc_probs
    double precision :: altalpha, altlambda, altmu, altnu, p1, p2, r

    integer, dimension(:), allocatable :: nbrs
    integer :: alti, i, isize, N, t

    !Number of nodes.
    N = size(history)

    !We choose the local epidemiological parameters based on the model in use.
    if (type == 'SIR') then
      altnu = mu
      altmu = 0.0d0
      altalpha = lambda
      altlambda = 0.0d0
    else
      altmu = mu
      altnu = nu
      altalpha = alpha
      altlambda = lambda
    end if

    tmp_mc_probs = 0.0d0
    do alti = 1, realizations
      !We initialize the array of states.
      states = 'S'
      states(origins) = 'E'
      !We initialize an auxiliary array of states.
      auxsts = states

      mc_probs = 0.0d0
      do t = 1, t0
        !We apply the MC core to modify the states using the S(E)IR rules.
        do i = 1, N
          !We compute the neighbors of node i at time t.
          isize = size(indices(i)%time(t)%array)
          allocate(nbrs(isize))
          nbrs = history(i)%neighbors(indices(i)%time(t)%array)

          !We apply the S(E)IR rules at each t: an E (exposed) node becomes I with
          !probability nu, and an I node becomes R with probability mu, regardless
          !of their respective neighbors; a S node becomes E with probability alpha
          !when it makes contact with one of its E neighbors (modified SEIR); also,
          !it becomes E with probability lambda when it interacts with one of its I
          !neighbors.
          r = r1279()
          if ((auxsts(i) == 'E').and.(r < altnu)) then
            states(i) = 'I'
          else if ((auxsts(i) == 'I').and.(r < altmu)) then
            states(i) = 'R'
          else if (auxsts(i) == 'S') then
            p1 = 1.0d0 - product(1.0d0-altalpha*merge(1, 0, auxsts(nbrs)=='E'))
            p2 = 1.0d0 - product(1.0d0-altlambda*merge(1, 0, auxsts(nbrs)=='I'))
            if (r < p1+p2) states(i) = 'E'
          end if

          deallocate(nbrs)
        end do
        !We update auxsts at the end of each time step.
        auxsts = states

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

    !If we are considering the SIR model, in the end we have to perform the
    !following exchanges: (PE, PI, PR) --> (0.0, PE, PI) and ('E', 'I') -->
    !('I', 'R').
    if (type == 'SIR') then
      tmp_mc_probs(:, 3:4) = tmp_mc_probs(:, 2:3)
      tmp_mc_probs(:, 2) = 0.0d0

      do i = 1, N
        if (states(i) == 'I') then
          states(i) = 'R'
        else if (states(i) == 'E') then
          states(i) = 'I'
        end if
      end do
    end if
  end subroutine mc_sim
end module mc_simulations
