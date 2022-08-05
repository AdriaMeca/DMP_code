!> Contains the dynamic message-passing (DMP) algorithm for the S(E)IR model on
!> time-varying networks.
!> Author: Adria Meca Montserrat.
!> Last modified date: 05/08/22.
module dmp_algorithms
  use array_procedures, only: dbl_list, int_llist, my_pack, pop
  use network_generation, only: node

  implicit none

  private

  public :: dmp

contains
  !> DMP(r) algorithm for the S(E)IR model.
  subroutine dmp( &
    model,        &
    restricted,   &
    history,      &
    indices,      &
    states,       &
    origins,      &
    alpha_,       &
    lambda_,      &
    mu_,          &
    nu_,          &
    t0,           &
    ps,           &
    pe,           &
    pi,           &
    pr,           &
    tmp_dmp_probs &
  )
    character(len=*),              intent(in)  :: model                                                !> Epidemiological model.
    double precision,              intent(in)  :: alpha_, lambda_, mu_, nu_                            !>
    double precision,              intent(out) :: ps(:), pe(:), pi(:), pr(:)                           !> DMP marginal probabilities.
    double precision,              intent(out) :: tmp_dmp_probs(:, :)                                  !> DMP trajectories.
    character(len=1),              intent(in)  :: states(:)                                            !> Node states.
    double precision                           :: alpha, lambda, mu, nu                                !> Epidemiological parameters.
    double precision                           :: ps0(size(history))                                   !>
    double precision                           :: theta_ki, xi                                         !>
    integer,                       intent(in)  :: origins(:)                                           !> Patient zeros.
    integer,                       intent(in)  :: t0                                                   !> Observation time.
    integer,          allocatable              :: non_susceptible(:)                                   !>
    integer                                    :: alti, altk, gsize, hsize, i, ik, isize, k, ki, N, t  !>
    logical,                       intent(in)  :: restricted                                           !>
    type(dbl_list)                             :: nps(size(history))                                   !>
    type(dbl_list)                             :: ops(size(history))                                   !>
    type(dbl_list)                             :: phi(size(history))                                   !>
    type(dbl_list)                             :: psi(size(history))                                   !>
    type(dbl_list)                             :: theta(size(history))                                 !>
    type(int_llist),               intent(in)  :: indices(:)                                           !> Active links throughout the simulation.
    type(node),                    intent(in)  :: history(:)                                           !> Rewiring history.

    !> Number of nodes.
    N = size(history)

    !> We choose the local epidemiological parameters based on the model in use.
    if (trim(model) == 'SIR') then
      nu = mu_
      mu = 0.0d0
      xi = 0.0d0
      alpha = lambda_
      lambda = 0.0d0
    else
      mu = mu_
      nu = nu_
      xi = nu_
      alpha = alpha_
      lambda = lambda_
    end if

    !> Initial conditions for the marginal probabilities.
    ps0 = 1.0d0
    ps0(origins) = 0.0d0

    ps = ps0
    pe = 1.0d0 - ps0
    pi = 0.0d0
    pr = 0.0d0

    !> Initial conditions for the DMP(r) messages.
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

    !> In the restricted version of the DMP algorithm we only iterate over the
    !> nodes whose state is different from S.
    allocate(non_susceptible(N))
    non_susceptible = [(i, i=1,N)]
    if (restricted) call my_pack(non_susceptible, states /= 'S')
    gsize = size(non_susceptible)

    tmp_dmp_probs = 0.0d0
    do t = 1, t0
      !> We apply the DMP(r) iteration to compute the marginal probabilities that
      !> each node is in a given state at time t.
      do alti = 1, gsize
        i = non_susceptible(alti)

        hsize = size(history(i)%neighbors)
        isize = size(indices(i)%time(t)%array)
        if (isize > 0) then
          !> We update the thetas that enter node i at time t.
          do altk = 1, isize
            ki = indices(i)%time(t)%array(altk)
            theta(i)%array(ki) = theta(i)%array(ki) - alpha*psi(i)%array(ki)
            theta(i)%array(ki) = theta(i)%array(ki) - lambda*phi(i)%array(ki)
          end do

          !> We calculate the probability that node i is S at time t.
          ps(i) = ps0(i) * product(theta(i)%array)

          !> We update the npses that leave node i at time t.
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

        !> We update the phis and psis that enter node i at time t.
        isize = size(indices(i)%time(t)%array)
        do altk = 1, isize
          ki = indices(i)%time(t)%array(altk)
          psi(i)%array(ki) = (1.0d0-alpha) * psi(i)%array(ki)
          phi(i)%array(ki) = (1.0d0-lambda) * phi(i)%array(ki)
        end do
        phi(i)%array = (1.0d0-mu)*phi(i)%array + xi*psi(i)%array
        psi(i)%array = (1.0d0-nu)*psi(i)%array + ops(i)%array - nps(i)%array

        !> We update the opses that enter node i at time t.
        ops(i)%array = nps(i)%array

        !> We compute the other marginal probabilities at time t.
        pr(i) = pr(i) + mu*pi(i)
        pi(i) = (1.0d0-mu)*pi(i) + nu*pe(i)
        pe(i) = 1.0d0 - ps(i) - pi(i) - pr(i)
      end do

      !> We use the previous marginals to predict the time evolution of the
      !> fraction of nodes that are in state X (i.e., S, E, I or R).
      do i = 1, N
        tmp_dmp_probs(t, 1) = tmp_dmp_probs(t, 1) + ps(i)
        tmp_dmp_probs(t, 2) = tmp_dmp_probs(t, 2) + pe(i)
        tmp_dmp_probs(t, 3) = tmp_dmp_probs(t, 3) + pi(i)
        tmp_dmp_probs(t, 4) = tmp_dmp_probs(t, 4) + pr(i)
      end do
    end do

    !> We normalize the trajectories.
    tmp_dmp_probs = tmp_dmp_probs / N

    !> If we are considering the SIR model, in the end we have to perform the
    !> following exchanges: (pe, pi, pr) --> (0.0d0, pe, pi).
    if (trim(model) == 'SIR') then
      tmp_dmp_probs(:, 3:4) = tmp_dmp_probs(:, 2:3)
      tmp_dmp_probs(:, 2) = 0.0d0

      pr = pi
      pi = pe
      pe = 0.0d0
    end if
  end subroutine dmp
end module dmp_algorithms
