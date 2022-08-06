!> Contains the dynamic message-passing (DMP) algorithm for the S(E)IR model on
!> time-varying networks.
!> Author: Adria Meca Montserrat.
!> Last modified date: 06/08/22.
module dmp_algorithms
  use array_procedures, only: dbl_list, int_llist, my_pack, pop
  use network_generation, only: node

  implicit none

  private

  public :: dmp

contains
  !> DMP(r) algorithm for the S(E)IR model.
  subroutine dmp(  &
    model,         &
    history,       &
    indices,       &
    nodes,         &
    patient_zeros, &
    alpha_,        &
    lambda_,       &
    mu_,           &
    nu_,           &
    t0,            &
    p              &
  )
    character(len=*), intent(in)  :: model                                         !> Epidemiological model.
    double precision, intent(in)  :: alpha_, lambda_, mu_, nu_                     !>
    double precision, intent(out) :: p(:, 0:, :)                                   !> DMP marginal probabilities.
    double precision              :: alpha, lambda, mu, nu                         !> Epidemiological parameters.
    double precision              :: theta_ki, xi                                  !>
    integer,          intent(in)  :: nodes(:)                                      !>
    integer,          intent(in)  :: patient_zeros(:)                              !> Patient zeros.
    integer,          intent(in)  :: t0                                            !> Observation time.
    integer                       :: alti, altk, hsize, i, ik, isize, k, ki, N, t  !>
    type(int_llist),  intent(in)  :: indices(:)                                    !> Active links throughout the simulation.
    type(node),       intent(in)  :: history(:)                                    !> Rewiring history.
    type(dbl_list)                :: nps(size(history))                            !>
    type(dbl_list)                :: ops(size(history))                            !>
    type(dbl_list)                :: phi(size(history))                            !>
    type(dbl_list)                :: psi(size(history))                            !>
    type(dbl_list)                :: theta(size(history))                          !>

    !> Number of nodes.
    N = size(history)

    !> We choose the local epidemiological parameters based on the model in use.
    select case (trim(model))
      case ('SIR')
        nu = mu_
        mu = 0.0d0
        xi = 0.0d0
        alpha = lambda_
        lambda = 0.0d0
      case default
        mu = mu_
        nu = nu_
        xi = nu_
        alpha = alpha_
        lambda = lambda_
    end select

    !> Initial conditions for the marginal probabilities.
    p(1, 0, :) = 1.0d0; p(1, 0, patient_zeros) = 0.0d0
    p(2, 0, :) = 1.0d0 - p(1, 0, :)
    p(3, 0, :) = 0.0d0
    p(4, 0, :) = 0.0d0

    !> Initial conditions for the DMP(r) messages.
    do i = 1, N
      hsize = size(history(i)%neighbors)

      allocate(ops(i)%array(hsize), nps(i)%array(hsize), phi(i)%array(hsize), &
        psi(i)%array(hsize), theta(i)%array(hsize))

      do ki = 1, hsize
        k = history(i)%neighbors(ki)

        ops(i)%array(ki) = p(1, 0, k)
        nps(i)%array(ki) = p(1, 0, k)
        phi(i)%array(ki) = p(3, 0, k)
        psi(i)%array(ki) = p(2, 0, k)
        theta(i)%array(ki) = 1.0d0
      end do
    end do

    do t = 1, t0
      !> We apply the DMP(r) iteration to compute the marginal probabilities that
      !> each node is in a given state at time t.
      do alti = 1, size(nodes)
        i = nodes(alti)

        isize = size(indices(i)%time(t)%array)
        if (isize > 0) then
          !> We update the thetas that enter node i at time t.
          do altk = 1, isize
            ki = indices(i)%time(t)%array(altk)
            theta(i)%array(ki) = theta(i)%array(ki) - alpha*psi(i)%array(ki)
            theta(i)%array(ki) = theta(i)%array(ki) - lambda*phi(i)%array(ki)
          end do

          !> We calculate the probability that node i is S at time t.
          p(1, t, i) = p(1, 0, i) * product(theta(i)%array)

          !> We update the npses that leave node i at time t.
          do ki = 1, size(history(i)%neighbors)
            k = history(i)%neighbors(ki)
            ik = history(i)%opposites(ki)
            theta_ki = theta(i)%array(ki)
            if (theta_ki > 0.0d0) then
              nps(k)%array(ik) = p(1, t, i) / theta_ki
            else
              nps(k)%array(ik) = p(1, 0, i) * product(pop(theta(i)%array, ki))
            end if
          end do
        end if
      end do

      do alti = 1, size(nodes)
        i = nodes(alti)

        !> We update the phis and psis that enter node i at time t.
        do altk = 1, size(indices(i)%time(t)%array)
          ki = indices(i)%time(t)%array(altk)
          psi(i)%array(ki) = (1.0d0-alpha) * psi(i)%array(ki)
          phi(i)%array(ki) = (1.0d0-lambda) * phi(i)%array(ki)
        end do
        phi(i)%array = (1.0d0-mu)*phi(i)%array + xi*psi(i)%array
        psi(i)%array = (1.0d0-nu)*psi(i)%array + ops(i)%array - nps(i)%array

        !> We update the opses that enter node i at time t.
        ops(i)%array = nps(i)%array

        !> We compute the other marginal probabilities at time t.
        p(4, t, i) = p(4, t-1, i) + mu*p(3, t-1, i)
        p(3, t, i) = p(3, t-1, i) - mu*p(3, t-1, i) + nu*p(2, t-1, i)

        p(2, t, i) = 1.0d0 - p(1, t, i) - p(3, t, i) - p(4, t, i)
      end do
    end do

    !> If we are considering the SIR model, in the end we have to perform the
    !> following exchanges: (pe, pi, pr) --> (0.0d0, pe, pi).
    if (trim(model) == 'SIR') p(3:4, :, :) = p(2:3, :, :); p(2, :, :) = 0.0d0
  end subroutine dmp
end module dmp_algorithms
