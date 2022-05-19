!Module containing dynamic message-passing algorithms (DMP) for different
!epidemiological models (SIR and SEIR) in time-varying networks.
!Author: Adrià Meca Montserrat.
!Last modified date: 19/05/22.
module dmp_algorithms
  use array_procedures, only : dbl_list, int_list_list, pop
  use network_generation, only : node

  implicit none

  private

  public dmp_seir, dmp_sir

contains
  !DMP iteration for the SEIR model.
  subroutine dmp_seir(history, indices, origin, alpha, lambda, mu, nu, t0, &
    ps, pe, pi, pr, tmp_dmp_probs)
    implicit none

    !Input arguments.
    double precision, intent(in) :: alpha, lambda, mu, nu
    integer, intent(in) :: origin, t0
    type(int_list_list), dimension(:), intent(in) :: indices
    type(node), dimension(:), intent(in) :: history

    !Output arguments.
    double precision, dimension(:, :), intent(out) :: tmp_dmp_probs
    double precision, dimension(:), intent(out) :: ps, pe, pi, pr

    !Local variables.
    double precision, dimension(size(history)) :: ps0
    integer :: hsize, i, id_ki, k, N, t
    type(dbl_list), dimension(size(history)) :: ops, nps, phi, psi, theta

    !Number of nodes.
    N = size(history)

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

      do id_ki = 1, hsize
        k = history(i)%neighbors(id_ki)

        ops(i)%array(id_ki) = ps(k)
        nps(i)%array(id_ki) = ps(k)
        phi(i)%array(id_ki) = pi(k)
        psi(i)%array(id_ki) = pe(k)
        theta(i)%array(id_ki) = 1.0d0
      end do
    end do

    tmp_dmp_probs = 0.0d0
    do t = 1, t0
      !We apply the DMP algorithm to compute the marginal probabilities that
      !each node is in a given state at time t.
      call dmp_seir_kernel(history, indices, ops, nps, phi, psi, theta, &
        alpha, lambda, mu, nu, t, ps, ps0, pe, pi, pr)

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
  end subroutine dmp_seir


  !DMP algorithm for the SEIR model.
  subroutine dmp_seir_kernel(history, indices, ops, nps, phi, psi, theta, &
    alpha, lambda, mu, nu, t, ps, ps0, pe, pi, pr)
    implicit none

    !Input arguments.
    double precision, dimension(:), intent(in) :: ps0
    double precision, intent(in) :: alpha, lambda, mu, nu
    integer, intent(in) :: t
    type(int_list_list), dimension(:), intent(in) :: indices
    type(node), dimension(:), intent(in) :: history

    !Output arguments.
    double precision, dimension(:), intent(inout) :: ps, pe, pi, pr
    type(dbl_list), dimension(:), intent(inout) :: ops, nps, phi, psi, theta

    !Local variables.
    double precision :: theta_ki
    integer :: altk, hsize, i, id_ik, id_ki, isize, k, N

    !Number of nodes.
    N = size(history)

    do i = 1, N
      hsize = size(history(i)%neighbors)
      isize = size(indices(i)%time(t)%array)
      if (isize > 0) then
        !We update the thetas that enter node i at time t.
        do altk = 1, isize
          id_ki = indices(i)%time(t)%array(altk)
          theta(i)%array(id_ki) = theta(i)%array(id_ki) &
            - lambda*phi(i)%array(id_ki) - alpha*psi(i)%array(id_ki)
        end do

        !We calculate the probability that node i is S at time t.
        ps(i) = ps0(i) * product(theta(i)%array)

        !We update the npses that leave node i at time t.
        do id_ki = 1, hsize
          k = history(i)%neighbors(id_ki)
          id_ik = history(i)%opposites(id_ki)
          theta_ki = theta(i)%array(id_ki)
          if (theta_ki > 0.0d0) then
            nps(k)%array(id_ik) = ps(i) / theta_ki
          else
            nps(k)%array(id_ik) = ps0(i) * product(pop(theta(i)%array, id_ki))
          end if
        end do
      end if
    end do

    do i = 1, N
      !We update the phis and psis that enter node i at time t.
      isize = size(indices(i)%time(t)%array)
      do altk = 1, isize
        id_ki = indices(i)%time(t)%array(altk)
        phi(i)%array(id_ki) = (1.0d0-lambda) * phi(i)%array(id_ki)
      end do
      phi(i)%array = (1.0d0-mu)*phi(i)%array + nu*psi(i)%array
      do altk = 1, isize
        id_ki = indices(i)%time(t)%array(altk)
        phi(i)%array(id_ki) = phi(i)%array(id_ki) - nu*alpha*psi(i)%array(id_ki)
        psi(i)%array(id_ki) = (1.0d0-alpha) * psi(i)%array(id_ki)
      end do
      psi(i)%array = (1.0d0-nu)*psi(i)%array + ops(i)%array - nps(i)%array

      !We update the opses that enter node i at time t.
      ops(i)%array = nps(i)%array

      !We compute the other marginal probabilities at time t.
      pr(i) = pr(i) + mu*pi(i)
      pi(i) = (1.0d0-mu)*pi(i) + nu*pe(i)
      pe(i) = 1.0d0 - ps(i) - pi(i) - pr(i)
    end do
  end subroutine dmp_seir_kernel


  !DMP iteration for the SIR model.
  subroutine dmp_sir(history, indices, origin, lambda, mu, t0, ps, pi, pr)
    implicit none

    !Input arguments.
    double precision, intent(in) :: lambda, mu
    integer, intent(in) :: origin, t0
    type(int_list_list), dimension(:), intent(in) :: indices
    type(node), dimension(:), intent(in) :: history

    !Output arguments.
    double precision, dimension(:), intent(out) :: ps, pi, pr

    !Local variables.
    double precision, dimension(size(history)) :: ps0
    double precision :: theta_ki
    integer :: altk, hsize, i, id_ik, id_ki, isize, k, N, t
    type(dbl_list), dimension(size(history)) :: ops, nps, phi, theta

    !Number of nodes.
    N = size(history)

    !Initial conditions for the marginal probabilities.
    ps0 = 1.0d0
    ps0(origin) = 0.0d0

    ps = ps0
    pi = 1.0d0 - ps0
    pr = 0.0d0

    !Initial conditions for the DMP messages.
    do i = 1, N
      hsize = size(history(i)%neighbors)

      allocate(ops(i)%array(hsize), nps(i)%array(hsize), phi(i)%array(hsize), &
        theta(i)%array(hsize))

      do id_ki = 1, hsize
        k = history(i)%neighbors(id_ki)

        ops(i)%array(id_ki) = ps(k)
        nps(i)%array(id_ki) = ps(k)
        phi(i)%array(id_ki) = pi(k)
        theta(i)%array(id_ki) = 1.0d0
      end do
    end do

    !DMP algorithm for the SIR model.
    do t = 1, t0
      do i = 1, N
        hsize = size(history(i)%neighbors)
        isize = size(indices(i)%time(t)%array)
        if (isize > 0) then
          !We update the thetas that enter node i at time t.
          do altk = 1, isize
            id_ki = indices(i)%time(t)%array(altk)
            theta(i)%array(id_ki) = theta(i)%array(id_ki) - lambda*phi(i)%array(id_ki)
          end do

          !We calculate the probability that node i is S at time t.
          ps(i) = ps0(i) * product(theta(i)%array)

          !We update the npses that leave node i at time t.
          do id_ki = 1, hsize
            k = history(i)%neighbors(id_ki)
            id_ik = history(i)%opposites(id_ki)
            theta_ki = theta(i)%array(id_ki)
            if (theta_ki > 0.0d0) then
              nps(k)%array(id_ik) = ps(i) / theta_ki
            else
              nps(k)%array(id_ik) = ps0(i) * product(pop(theta(i)%array, id_ki))
            end if
          end do
        end if
      end do

      do i = 1, N
        !We update the phis that enter node i at time t.
        isize = size(indices(i)%time(t)%array)
        do altk = 1, isize
          id_ki = indices(i)%time(t)%array(altk)
          phi(i)%array(id_ki) = (1.0d0-lambda) * phi(i)%array(id_ki)
        end do
        phi(i)%array = (1.0d0-mu)*phi(i)%array + ops(i)%array - nps(i)%array

        !We update the opses that enter node i at time t.
        ops(i)%array = nps(i)%array

        !We compute the other marginal probabilities.
        pr(i) = pr(i) + mu*pi(i)
        pi(i) = 1.0d0 - ps(i) - pr(i)
      end do
    end do
  end subroutine dmp_sir
end module dmp_algorithms
