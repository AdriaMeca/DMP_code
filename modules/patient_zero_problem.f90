!> Procedures that simulate the patient-zero problem.
!> Author: Adria Meca Montserrat.
!> Last modified date: 06/08/22.
module patient_zero_problem
  use array_procedures, only: add, find, int_llist, my_pack, quicksort
  use dmp_algorithms, only: dmp
  use mc_simulations, only: mc_sim
  use network_generation, only: node
  use random_number_generator, only: r1279

  implicit none

  private

  public :: pz_sim

contains
  !> Simulates the patient-zero problem: we use the DMP probabilities to compute
  !> the energies of the non-susceptible nodes of an epidemic. We rank the nodes
  !> by increasing value of their energy and return the rank the algorithm gives
  !> to the true patient zero.
  subroutine pz_sim( &
    model,           &
    dmpr,            &
    history,         &
    indices,         &
    seeds,           &
    states,          &
    origins,         &
    ranks,           &
    guess,           &
    gsize,           &
    alpha,           &
    lambda,          &
    mu,              &
    nu,              &
    t0               &
  )
    character(len=*),              intent(in)  :: model                                      !> Epidemiological model.
    character(len=1),              intent(out) :: states(0:, :)                              !> Node states.
    double precision,              intent(in)  :: alpha, lambda, mu, nu                      !> Epidemiological parameters.
    double precision, allocatable              :: energies(:)                                !> Node energies.
    double precision, parameter                :: small = 1.0d-300                           !>
    double precision                           :: joint                                      !>
    integer,                       intent(in)  :: seeds                                      !> Number of patient zeros.
    integer,                       intent(in)  :: t0                                         !> Observation time.
    double precision                           :: tmp_dmp_probs(t0, 4), tmp_mc_probs(t0, 4)  !>
    integer,                       intent(out) :: gsize(:), origins(:), ranks(:)             !>
    integer,          allocatable              :: non_susceptible(:)                         !>
    integer                                    :: counter, i, idx, N, node_i, total          !>
    integer                                    :: guess                                      !> Predicted patient zero.
    logical,                       intent(in)  :: dmpr                                       !>
    type(int_llist),               intent(in)  :: indices(:)                                 !> Active links throughout the simulation.
    type(node),                    intent(in)  :: history(:)                                 !> Rewiring history.
    double precision                           :: pe(size(history))                          !>
    double precision                           :: pi(size(history))                          !>
    double precision                           :: pr(size(history))                          !>
    double precision                           :: ps(size(history))                          !>

    !> Number of nodes.
    N = size(history)

    !> We neglect those instances where the epidemic does not spread.
    do while (.true.)
      !> We choose the origins at random.
      origins = [(1 + mod(int(N*r1279()), N), i=1,seeds)]

      !> We let the system evolve by simulating the S(E)IR rules with MC up to t0.
      call mc_sim(   &
        model,       &
        history,     &
        indices,     &
        origins,     &
        alpha,       &
        lambda,      &
        mu,          &
        nu,          &
        t0,          &
        states,      &
        1,           &
        tmp_mc_probs &
      )

      !> If the epidemic spread, we exit the loop.
      if (count(states(t0, :) /= 'S') > 1) exit
    end do

    !> Only non-susceptible nodes may have started the epidemic.
    gsize = 0
    counter = 0
    do i = 1, N
      select case (states(t0, i))
        case ('S')
          cycle
        case ('E')
          gsize(1) = gsize(1) + 1
        case ('I')
          gsize(2) = gsize(2) + 1
        case ('R')
          gsize(3) = gsize(3) + 1
      end select
      counter = counter + 1
      call add(non_susceptible, i)
    end do
    allocate(energies(counter))

    do idx = 1, counter
      node_i = non_susceptible(idx)

      !> We calculate the DMP(r) marginal probabilities that each node is in state
      !> X (i.e., S, E, I or R) at t0 given that node i is the patient zero.
      call dmp(        &
        model,         &
        dmpr,          &
        history,       &
        indices,       &
        states(t0, :), &
        [node_i],      &
        alpha,         &
        lambda,        &
        mu,            &
        nu,            &
        t0,            &
        ps,            &
        pe,            &
        pi,            &
        pr,            &
        tmp_dmp_probs  &
      )

      !> The probability that node i is the patient zero given the observation
      !> of the network at time t0 P(i|O) is proportional to P(O|i) (thanks to
      !> Bayes' theorem); we approximate the latter joint probability by a
      !> product of the DMP(r) marginals.
      total = 0
      joint = 1.0d0
      do i = 1, N
        select case (states(t0, i))
          case ('S')
            joint = joint * ps(i)
          case ('E')
            joint = joint * pe(i)
          case ('I')
            joint = joint * pi(i)
          case ('R')
            joint = joint * pr(i)
        end select

        !> It may happen that 'joint' becomes so small that GFortran converts its
        !> value to zero, losing all the information it carries. To avoid this,
        !> we multiply it by a huge constant every time it surpasses a certain
        !> threshold. Next, we use the properties of logarithms to compute the
        !> correct energy of a given node.
        if (abs(joint) < small) then
          joint = joint / small
          total = total + 1
        end if
      end do

      !> Sometimes joint <= 0, so to avoid having to deal with infinities we
      !> associate a large constant energy to the node in question.
      if (joint > 0.0d0) then
        energies(idx) = total*log(1/small) - log(joint)
      else
        energies(idx) = 1 / small
      end if
    end do

    !> We rank non-susceptible nodes by increasing value of their energy.
    call quicksort(energies, non_susceptible)

    !> We compute the ranks of the true patient zeros.
    ranks = [(find(non_susceptible, origins(i)) - 1, i=1,seeds)]

    !> We return the node that the algorithm believes to be the true patient zero.
    guess = non_susceptible(1)
  end subroutine pz_sim
end module patient_zero_problem
