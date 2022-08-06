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
    patient_zeros,   &
    alpha_,          &
    lambda_,         &
    mu_,             &
    nu_,             &
    t0,              &
    states           &
  )
    character(len=*), intent(in)  :: model                      !> Epidemiological model.
    character(len=1), intent(out) :: states(0:, :)              !> Node states.
    double precision, intent(in)  :: alpha_, lambda_, mu_, nu_  !>
    double precision              :: alpha, lambda, mu, nu      !> Epidemiological parameters.
    double precision              :: p1, p2, r                  !>
    integer,          intent(in)  :: patient_zeros(:)           !>
    integer,          intent(in)  :: t0                         !> Observation time.
    integer                       :: altk, i, k, N, t           !>
    type(int_llist),  intent(in)  :: indices(:)                 !> Active links throughout the simulation.
    type(node),       intent(in)  :: history(:)                 !> Rewiring history.

    !> Number of nodes.
    N = size(history)

    !> We choose the local epidemiological parameters based on the model in use.
    select case (trim(model))
      case ('SIR')
        nu = mu_
        mu = 0.0d0
        alpha = lambda_
        lambda = 0.0d0
      case default
        mu = mu_
        nu = nu_
        alpha = alpha_
        lambda = lambda_
    end select

    !> We initialize the distribution of node states.
    states = 'S'
    states(:, patient_zeros) = 'E'

    !> We let the system evolve by simulating the S(E)IR rules with MC.
    do t = 0, t0
      do i = 1, N
        r = r1279()

        if ((states(t, i) == 'E') .and. (r < nu)) then
          states(t+1:, i) = 'I'
        else if ((states(t, i) == 'I') .and. (r < mu)) then
          states(t+1:, i) = 'R'
        else if (states(t, i) == 'S') then
          p1 = 1.0d0
          p2 = 1.0d0

          do altk = 1, size(indices(i)%time(t)%array)
            k = history(i)%neighbors(indices(i)%time(t)%array(altk))

            p1 = p1 * (1.0d0 - alpha*merge(1.0d0, 0.0d0, states(t, k) == 'E'))
            p2 = p2 * (1.0d0 - lambda*merge(1.0d0, 0.0d0, states(t, k) == 'I'))
          end do

          p1 = 1.0d0 - p1
          p2 = 1.0d0 - p2

          if (r < p1+p2) states(t+1:, i) = 'E'
        end if
      end do

      !> We update auxsts at the end of each time step.
      states(t, :) = states(t+1, :)
    end do

    !> For the SIR model, we must perform the following exchanges at the end of
    !> the simulation: (E, I) --> (I, R).
    if (trim(model) == 'SIR') then
      where (states == 'E')
        states = 'I'
      elsewhere (states == 'I')
        states = 'R'
      end where
    end if
  end subroutine mc_sim
end module mc_simulations
