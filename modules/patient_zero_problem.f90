!Module whose procedures simulate the patient zero problem in which we try to
!locate the node(s) that started an epidemic in a given network.
!Author: Adria Meca Montserrat.
!Last modified date: 16/06/22.
module patient_zero_problem
  use array_procedures, only : add, find, int_llist, my_pack, quicksort
  use dmp_algorithms, only : dmp
  use mc_simulations, only : mc_sim
  use network_generation, only : node
  use random_number_generator, only : r1279

  implicit none

  private

  public pz_sim

contains
  !Subroutine that simulates the problem of the patient zero: given a snapshot
  !of an epidemic, we rank the non-susceptible nodes by their energy, which is
  !calculated using the marginal DMP probabilities. Then, we extract the ranks
  !the algorithm gives to the original seeds and the size of the epidemic.
  subroutine pz_sim(model, dmpr, history, indices, seeds, states, origins, &
    ranks, guess, gsize, alpha, lambda, mu, nu, t0)
    implicit none

    !Input arguments.
    character(len=*), intent(in) :: model

    double precision, intent(in) :: alpha, lambda, mu, nu

    integer, intent(in) :: seeds, t0

    logical, intent(in) :: dmpr

    type(int_llist), dimension(:), intent(in) :: indices

    type(node), dimension(:), intent(in) :: history

    !Output arguments.
    character(len=1), dimension(:), intent(out) :: states

    integer, dimension(:), intent(out) :: gsize, origins, ranks
    integer :: guess

    !Local variables.
    double precision, dimension(t0, 4) :: tmp_dmp_probs, tmp_mc_probs
    double precision, dimension(:), allocatable :: energies
    double precision, dimension(size(history)) :: ps, pe, pi, pr
    double precision, parameter :: small=1.0d-300
    double precision :: joint

    integer, dimension(:), allocatable :: non_susceptible
    integer :: counter, i, idx, N, node_i, total


    !Number of nodes.
    N = size(history)

    !We neglect those cases where the epidemic does not spread.
    do while (.true.)
      !We choose the origins at random.
      origins = [(1 + mod(int(N*r1279()), N), i=1,seeds)]

      !We apply a Monte Carlo simulation to spread an infection on the network
      !using the rules of the S(E)IR model.
      call mc_sim(model, history, indices, origins, alpha, lambda, mu, nu, &
        t0, states, 1, tmp_mc_probs)

      !If the epidemic spread, we exit the loop.
      if (count(states /= 'S') > 1) exit
    end do

    !We look for nodes that are not S in the network at time t0 because they
    !are the ones that possibly started the epidemic.
    gsize = 0
    counter = 0
    do i = 1, N
      select case (states(i))
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

    !We loop over each of the non-susceptible nodes.
    do idx = 1, counter
      node_i = non_susceptible(idx)

      !We apply the DMP(r) algorithm to calculate the marginal probabilities
      !that each node is in a state X (i.e., S, (E), I or R) at time t0,
      !conditional on the fact that node i is the patient zero.
      call dmp(model, dmpr, history, indices, states, [node_i], alpha, &
        lambda, mu, nu, t0, ps, pe, pi, pr, tmp_dmp_probs)

      !The probability that node i is the patient zero given the observation
      !of the network at time t0 P(i|O) is proportional to P(O|i) (thanks to
      !Bayes' theorem); we approximate the later joint probability by a
      !product of the DMP(r) marginals.
      total = 0
      joint = 1.0d0
      do i = 1, N
        select case (states(i))
          case ('S')
            joint = joint * ps(i)
          case ('E')
            joint = joint * pe(i)
          case ('I')
            joint = joint * pi(i)
          case ('R')
            joint = joint * pr(i)
        end select

        !It may happen that 'joint' becomes so small that GFortran converts its
        !value to zero, losing all the information it carries. To avoid this,
        !we multiply it by a huge constant every time it surpasses a certain
        !threshold. Next, we use the properties of logarithms to compute the
        !correct energy of a given node.
        if (abs(joint) < small) then
          joint = joint / small
          total = total + 1
        end if
      end do

      !Sometimes joint <= 0, so to avoid having to deal with infinities we
      !associate a large constant energy to the node in question.
      if (joint > 0.0d0) then
        energies(idx) = total*log(1/small) - log(joint)
      else
        energies(idx) = 1 / small
      end if
    end do

    !We rank non-susceptible nodes by their energy in ascending order (0 <= r).
    call quicksort(energies, non_susceptible)
    !We compute the ranks of the true seeds.
    ranks = [(find(non_susceptible, origins(i))-1, i=1,seeds)]
    !We return the node that the algorithm believes to be the true patient zero.
    guess = non_susceptible(1)
  end subroutine pz_sim
end module patient_zero_problem
