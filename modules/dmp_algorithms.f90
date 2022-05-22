!Module containing dynamic message-passing algorithms (DMP) for different
!epidemiological models (SIR and SEIR) in time-varying networks.
!Author: Adrià Meca Montserrat.
!Last modified date: 22/05/22.
module dmp_algorithms
  use array_procedures, only : dbl_list, int_list_list, pop
  use network_generation, only : node

  implicit none

  private

  public dmp

contains
  !DMP algorithm for the SIR and SEIR models.
  subroutine dmp(type, history, indices, origin, alpha, lambda, mu, nu, t0, &
    ps, pe, pi, pr, tmp_dmp_probs)
    implicit none

    !Input arguments.
    character(len=*), intent(in) :: type
    double precision, intent(in) :: alpha, lambda, mu, nu
    integer, intent(in) :: origin, t0
    type(int_list_list), dimension(:), intent(in) :: indices
    type(node), dimension(:), intent(in) :: history

    !Output arguments.
    double precision, dimension(:, :), intent(out) :: tmp_dmp_probs
    double precision, dimension(:), intent(out) :: ps, pe, pi, pr

    !Local variables.
    double precision, dimension(size(history)) :: ps0
    double precision :: altnu, theta_ki
    integer :: altk, hsize, i, ik, isize, k, ki, N, t
    type(dbl_list), dimension(size(history)) :: ops, nps, phi, psi, theta

    !Number of nodes.
    N = size(history)

    !Number that alternates the SEIR and SIR models.
    altnu = merge(nu, 0.0d0, type=='SEIR')

    !Initial conditions for the marginal probabilities.
    ps0 = 1.0d0
    ps0(origin) = 0.0d0

    ps = ps0
    pe = 1.0d0 - ps0
    pi = 0.0d0
    pr = 0.0d0

    !Initial conditions for the DMP messages.
    do i = 1, N
      hsize = size(history(i)%neighbors)

      allocate(ops(i)%array(hsize), nps(i)%array(hsize), phi(i)%array(hsize), &
        psi(i)%array(hsize), theta(i)%array(hsize))

      do ki = 1, hsize
        k = history(i)%neighbors(ki)

        ops(i)%array(ki) = ps(k)
        nps(i)%array(ki) = ps(k)
        phi(i)%array(ki) = pi(k)
        psi(i)%array(ki) = pe(k)
        theta(i)%array(ki) = 1.0d0
      end do
    end do

    tmp_dmp_probs = 0.0d0
    do t = 1, t0
      !We apply the DMP iteration to compute the marginal probabilities that
      !each node is in a given state at time t.
      do i = 1, N
        hsize = size(history(i)%neighbors)
        isize = size(indices(i)%time(t)%array)
        if (isize > 0) then
          !We update the thetas that enter node i at time t.
          do altk = 1, isize
            ki = indices(i)%time(t)%array(altk)
            theta(i)%array(ki) = theta(i)%array(ki) - alpha*psi(i)%array(ki) &
              - lambda*phi(i)%array(ki)
          end do

          !We calculate the probability that node i is S at time t.
          ps(i) = ps0(i) * product(theta(i)%array)

          !We update the npses that leave node i at time t.
          do ki = 1, hsize
            k = history(i)%neighbors(ki)
            ik = history(i)%opposites(ki)
            theta_ki = theta(i)%array(ki)
            if (theta_ki > 0.0d0) then
              nps(k)%array(ik) = ps(i) / theta_ki
            else
              nps(k)%array(ik) = ps0(i) * product(pop(theta(i)%array, ki))
            end if
          end do
        end if
      end do

      do i = 1, N
        !We update the phis and psis that enter node i at time t.
        isize = size(indices(i)%time(t)%array)
        do altk = 1, isize
          ki = indices(i)%time(t)%array(altk)
          psi(i)%array(ki) = (1.0d0-alpha) * psi(i)%array(ki)
          phi(i)%array(ki) = (1.0d0-lambda) * phi(i)%array(ki)
        end do
        phi(i)%array = (1.0d0-mu)*phi(i)%array + altnu*psi(i)%array
        psi(i)%array = (1.0d0-nu)*psi(i)%array + ops(i)%array - nps(i)%array

        !We update the opses that enter node i at time t.
        ops(i)%array = nps(i)%array

        !We compute the other marginal probabilities at time t.
        pr(i) = pr(i) + mu*pi(i)
        pi(i) = (1.0d0-mu)*pi(i) + nu*pe(i)
        pe(i) = 1.0d0 - ps(i) - pi(i) - pr(i)
      end do

      !We use the previous marginals to predict the time evolution of the
      !fraction of nodes that are in state X (i.e., S, E, I or R).
      do i = 1, N
        tmp_dmp_probs(t, 1) = tmp_dmp_probs(t, 1) + ps(i)
        tmp_dmp_probs(t, 2) = tmp_dmp_probs(t, 2) + pe(i)
        tmp_dmp_probs(t, 3) = tmp_dmp_probs(t, 3) + pi(i)
        tmp_dmp_probs(t, 4) = tmp_dmp_probs(t, 4) + pr(i)
      end do
    end do
    !We normalize the trajectories.
    tmp_dmp_probs = tmp_dmp_probs / N
  end subroutine dmp
end module dmp_algorithms
