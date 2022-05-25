!Module containing dynamic message-passing algorithms (DMP) for different
!epidemiological models (SIR and SEIR) in time-varying networks.
!Author: Adrià Meca Montserrat.
!Last modified date: 25/05/22.
module dmp_algorithms
  use array_procedures, only : dbl_list, int_llist, my_pack, pop
  use network_generation, only : node

  implicit none

  private

  public dmp

contains
  !DMP algorithm for the SIR and SEIR models.
  subroutine dmp(type, restricted, history, indices, states, origins, &
    alpha, lambda, mu, nu, t0, ps, pe, pi, pr, tmp_dmp_probs)
    implicit none

    !Input arguments.
    character(len=1), dimension(:), intent(in) :: states
    character(len=*), intent(in) :: type

    double precision, intent(in) :: alpha, lambda, mu, nu

    integer, dimension(:), intent(in) :: origins
    integer, intent(in) :: t0

    logical, intent(in) :: restricted

    type(int_llist), dimension(:), intent(in) :: indices

    type(node), dimension(:), intent(in) :: history

    !Output arguments.
    double precision, dimension(:, :), intent(out) :: tmp_dmp_probs
    double precision, dimension(:), intent(out) :: ps, pe, pi, pr

    !Local variables.
    double precision, dimension(size(history)) :: ps0
    double precision :: altalpha, altlambda, altmu, altnu, altxi, theta_ki

    integer, dimension(:), allocatable :: non_susceptible
    integer :: alti, altk, gsize, hsize, i, ik, isize, k, ki, N, t

    type(dbl_list), dimension(size(history)) :: ops, nps, phi, psi, theta

    !Number of nodes.
    N = size(history)

    !We choose the local epidemiological parameters based on the model in use.
    if (type == 'SIR') then
      altnu = mu
      altmu = 0.0d0
      altxi = 0.0d0
      altalpha = lambda
      altlambda = 0.0d0
    else
      altmu = mu
      altnu = nu
      altxi = nu
      altalpha = alpha
      altlambda = lambda
    end if

    !Initial conditions for the marginal probabilities.
    ps0 = 1.0d0
    ps0(origins) = 0.0d0

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

    !In the restricted version of the DMP algorithm we only iterate over the
    !nodes whose state is different from S.
    allocate(non_susceptible(N))
    non_susceptible = [(i, i=1,N)]
    if (restricted) call my_pack(non_susceptible, states/='S')
    gsize = size(non_susceptible)

    tmp_dmp_probs = 0.0d0
    do t = 1, t0
      !We apply the DMP iteration to compute the marginal probabilities that
      !each node is in a given state at time t.
      do alti = 1, gsize
        i = non_susceptible(alti)

        hsize = size(history(i)%neighbors)
        isize = size(indices(i)%time(t)%array)
        if (isize > 0) then
          !We update the thetas that enter node i at time t.
          do altk = 1, isize
            ki = indices(i)%time(t)%array(altk)
            theta(i)%array(ki) = theta(i)%array(ki) - altalpha*psi(i)%array(ki) &
              - altlambda*phi(i)%array(ki)
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

      do alti = 1, gsize
        i = non_susceptible(alti)

        !We update the phis and psis that enter node i at time t.
        isize = size(indices(i)%time(t)%array)
        do altk = 1, isize
          ki = indices(i)%time(t)%array(altk)
          psi(i)%array(ki) = (1.0d0-altalpha) * psi(i)%array(ki)
          phi(i)%array(ki) = (1.0d0-altlambda) * phi(i)%array(ki)
        end do
        phi(i)%array = (1.0d0-altmu)*phi(i)%array + altxi*psi(i)%array
        psi(i)%array = (1.0d0-altnu)*psi(i)%array + ops(i)%array - nps(i)%array

        !We update the opses that enter node i at time t.
        ops(i)%array = nps(i)%array

        !We compute the other marginal probabilities at time t.
        pr(i) = pr(i) + altmu*pi(i)
        pi(i) = (1.0d0-altmu)*pi(i) + altnu*pe(i)
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

    !If we are considering the SIR model, in the end we have to perform the
    !following exchanges: (PE, PI, PR) --> (0.0, PE, PI).
    if (type == 'SIR') then
      tmp_dmp_probs(:, 3:4) = tmp_dmp_probs(:, 2:3)
      tmp_dmp_probs(:, 2) = 0.0d0

      pr = pi
      pi = pe
      pe = 0.0d0
    end if
  end subroutine dmp
end module dmp_algorithms
