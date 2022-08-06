!> Procedures that simulate the spread of epidemics on time-varying networks.
!> Author: Adria Meca Montserrat.
!> Last modified date: 06/08/22.
module mc_simulations
  use array_procedures, only: int_llist
  use network_generation, only: node
  use random_number_generator, only: r1279

  implicit none

  private

  public :: mc_sim

contains
  !> Spreads an S(E)IR epidemic on a time-varying network using Monte Carlo.
  subroutine mc_sim( &
    model,           &
    history,         &
    indices,         &
    origins,         &
    alpha_,          &
    lambda_,         &
    mu_,             &
    nu_,             &
    t0,              &
    tmp_states,      &
    repetitions,     &
    tmp_mc_probs     &
  )
    character(len=*),              intent(in)  :: model                      !> Epidemiological model.
    character(len=1),              intent(out) :: tmp_states(0:, :)          !> Node states.
    double precision,              intent(in)  :: alpha_, lambda_, mu_, nu_  !>
    double precision,              intent(out) :: tmp_mc_probs(:, :)         !> Temporal MC probabilities.
    double precision                           :: alpha, lambda, mu, nu      !> Epidemiological parameters.
    double precision                           :: p1, p2, r                  !>
    integer,                       intent(in)  :: origins(:)                 !> Patient zeros.
    integer,                       intent(in)  :: repetitions                !>
    integer,                       intent(in)  :: t0                         !> Observation time.
    integer,          allocatable              :: nbrs(:)                    !>
    integer                                    :: alti, i, isize, N, t       !>
    type(int_llist),               intent(in)  :: indices(:)                 !> Active links throughout the simulation.
    type(node),                    intent(in)  :: history(:)                 !> Rewiring history.
    character(len=1)                           :: auxsts(size(history))      !>
    character(len=1)                           :: states(size(history))      !>

    !> Number of nodes.
    N = size(history)

    !> We choose the local epidemiological parameters based on the model in use.
    if (trim(model) == 'SIR') then
      nu = mu_
      mu = 0.0d0
      alpha = lambda_
      lambda = 0.0d0
    else
      mu = mu_
      nu = nu_
      alpha = alpha_
      lambda = lambda_
    end if

    tmp_mc_probs = 0.0d0
    do alti = 1, repetitions
      !> We initialize the distribution of node states.
      states = 'S'
      states(origins) = 'E'

      !> We initialize an auxiliary array of states.
      auxsts = states

      tmp_states(0, :) = states

      !> We let the system evolve by simulating the S(E)IR rules with MC.
      do t = 1, t0
        do i = 1, N
          !> We retrieve the neighbors of node i at time t.
          isize = size(indices(i)%time(t)%array)
          allocate(nbrs(isize))
          nbrs = history(i)%neighbors(indices(i)%time(t)%array)

          r = r1279()
          if ((auxsts(i) == 'E') .and. (r < nu)) then
            states(i) = 'I'
          else if ((auxsts(i) == 'I') .and. (r < mu)) then
            states(i) = 'R'
          else if (auxsts(i) == 'S') then
            p1 = 1.0d0 - product(1.0d0-alpha*merge(1, 0, auxsts(nbrs) == 'E'))
            p2 = 1.0d0 - product(1.0d0-lambda*merge(1, 0, auxsts(nbrs) == 'I'))
            if (r < p1+p2) states(i) = 'E'
          end if

          deallocate(nbrs)
        end do

        !> We update auxsts at the end of each time step.
        auxsts = states

        tmp_states(t, :) = states

        !> We count the number of nodes that are in each state at time t.
        do i = 1, N
          select case (states(i))
            case ('S')
              tmp_mc_probs(t, 1) = tmp_mc_probs(t, 1) + 1.0d0
            case ('E')
              tmp_mc_probs(t, 2) = tmp_mc_probs(t, 2) + 1.0d0
            case ('I')
              tmp_mc_probs(t, 3) = tmp_mc_probs(t, 3) + 1.0d0
            case ('R')
              tmp_mc_probs(t, 4) = tmp_mc_probs(t, 4) + 1.0d0
          end select
        end do
      end do
    end do

    !> We divide the above quantities by N * repetitions to obtain the trajectory
    !> of the fraction of nodes in state X (i.e., S, E, I or R).
    tmp_mc_probs = tmp_mc_probs / N / repetitions

    !> For the SIR model, we must perform the following exchanges at the end of
    !> the simulation: (pe, pi, pr) --> (0.0d0, pe, pi) and (E, I) --> (I, R).
    if (trim(model) == 'SIR') then
      tmp_mc_probs(:, 3:4) = tmp_mc_probs(:, 2:3)
      tmp_mc_probs(:, 2) = 0.0d0

      do t = 0, t0
        do i = 1, N
          if (tmp_states(t, i) == 'I') then
            tmp_states(t, i) = 'R'
          else if (tmp_states(t, i) == 'E') then
            tmp_states(t, i) = 'I'
          end if
        end do
      end do
    end if
  end subroutine mc_sim
end module mc_simulations
