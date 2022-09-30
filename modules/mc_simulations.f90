!> Procedures that simulate the spread of epidemics on time-varying networks.
!> Author: Adria Meca Montserrat.
!> Last modified date: 30/09/22.
module mc_simulations
  use derived_types,           only: int_list, node, prm
  use random_number_generator, only: r1279

  implicit none

  private

  public :: mc_sim

contains
  !> Monte Carlo simulation propagating an S(E)IR epidemic on a time-varying network.
  subroutine mc_sim(model, history, indices, patient_zeros, epi_params, states)
    integer,          intent(in)  :: patient_zeros(:)      !>
    integer                       :: altk, i, k, N, t, t0  !>
    character(len=*), intent(in)  :: model                 !> Epidemiological model.
    character(len=1), intent(out) :: states(:, 0:)         !> Distribution of node states.
    double precision              :: p1, p2, r             !>
    double precision              :: pa, pl, pm, pn        !>
    type(int_list),   intent(in)  :: indices(:, 0:)        !> Active links throughout the simulation.
    type(node),       intent(in)  :: history(:)            !> Rewiring history.
    type(prm),        intent(in)  :: epi_params            !> Epidemiological parameters.

    !> Number of nodes.
    N = size(history)

    !> Local observation time.
    t0 = epi_params%t0

    !> We choose the local epidemiological parameters based on the model in use.
    select case (trim(model))
      case ('SIR')
        pa = epi_params%lambda
        pl = 0.0d0
        pm = 0.0d0
        pn = epi_params%mu
      case default
        pa = epi_params%alpha
        pl = epi_params%lambda
        pm = epi_params%mu
        pn = epi_params%nu
    end select

    !> We initialize the distribution of node states.
    states = 'S'; states(patient_zeros, :) = 'E'

    !> We generate an epidemic by simulating the S(E)IR rules with MC up to t0.
    do t = 1, t0
      do i = 1, N
        r = r1279()

        select case (states(i, t))
          case ('S')
            p1 = 1.0d0
            p2 = 1.0d0

            do altk = 1, size(indices(i, t)%array)
              k = history(i)%neighbors(indices(i, t)%array(altk))

              p1 = p1 * (1.0d0 - pa*merge(1.0d0, 0.0d0, states(k, t) == 'E'))
              p2 = p2 * (1.0d0 - pl*merge(1.0d0, 0.0d0, states(k, t) == 'I'))
            end do

            p1 = 1.0d0 - p1
            p2 = 1.0d0 - p2

            if (r < p1+p2) states(i, t+1:) = 'E'
          case ('E')
            if (r < pn) states(i, t+1:) = 'I'
          case ('I')
            if (r < pm) states(i, t+1:) = 'R'
        end select
      end do

      !> We update the state of all nodes at the end of each time step.
      states(:, t) = states(:, t+1)
    end do

    !> For the SIR model, we must perform the following exchanges:
    if (trim(model) == 'SIR') then
      where (states == 'E')
        states = 'I'
      elsewhere (states == 'I')
        states = 'R'
      end where
    end if
  end subroutine mc_sim
end module mc_simulations
