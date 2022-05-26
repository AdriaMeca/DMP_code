!Simulation that aims to test the efficiency of the DMP(r) algorithm by studying
!the precision with which it locates the true patient zero(s) in a network at a
!certain time t as we vary the lambda parameter (with all the others fixed).
!Author: Adria Meca Montserrat.
!Last modified date: 26/05/22.
program dmp_lambda
  use array_procedures, only : find, int_llist, my_pack, quicksort
  use dmp_algorithms, only : dmp
  use mc_simulations, only : mc_sim
  use network_generation, only : node, RRG
  use random_number_generator, only : r1279, setr1279
  use rewiring_algorithms, only : rewiring

  implicit none

  !Variables.
  character(len=1), dimension(:), allocatable :: states
  character(len=:), allocatable :: model

  double precision, dimension(:, :), allocatable :: tmp_dmp_probs, tmp_mc_probs
  double precision, dimension(:), allocatable :: energies, ps, pe, pi, pr
  double precision, parameter :: a=5.0d-2, b=1.0d0-a, small=1.0d-300
  double precision :: alpha, avg_g, avg_prob, avg_prob2, avg_rank, avg_rank2, &
    err_prob, err_rank, joint, l0, lambda, mu, nu, Q, rank

  integer, dimension(:), allocatable :: non_susceptible, origins, ranks
  integer :: c, gsize, i, idx, ins, instances, length, N, node_i, point, points, &
    seeds, t0, total

  logical :: restr

  type(int_llist), dimension(:), allocatable :: indices

  type(node), dimension(:), allocatable :: history, network


  !Parameters.
  open(unit=10, file='params.txt')
    read(10, *) restr
    read(10, *) length; allocate(character(len=length) :: model)
    read(10, *) model
    read(10, *) seeds, points, instances
    read(10, *) c, l0, N
    read(10, *) alpha, mu, nu, t0, Q
  close(10)

  !We allocate the main arrays once we have loaded the parameters that tell us
  !their size.
  allocate(origins(seeds), ranks(seeds))
  allocate(tmp_dmp_probs(t0, 4), tmp_mc_probs(t0, 4))
  allocate(states(N), ps(N), pe(N), pi(N), pr(N), indices(N), history(N), &
    network(N))

  !We initialize the random number generator with a fixed seed so that we can
  !reproduce the same results.
  call setr1279(1)

  !We generate the network and apply the rewiring algorithm to it, which creates
  !the history and indices.
  network = RRG(N, c)
  call rewiring('rrg', network, history, indices, c, t0, Q)

  !We vary the lambda parameter.
  do point = 1, points
    lambda = a + (b-a)*(point-1)/(points-1)

    !Initialization of the average quantities.
    avg_g = 0.0d0
    avg_prob = 0.0d0
    avg_rank = 0.0d0
    avg_prob2 = 0.0d0
    avg_rank2 = 0.0d0

    !For each lambda, we run multiple instances.
    do ins = 1, instances
      !We neglect those cases where the epidemic does not spread.
      do while (.true.)
        !We choose the origins at random.
        origins = [(1 + floor(N*r1279()), i=1,seeds)]

        !We apply a Monte Carlo simulation to spread an infection on the network
        !using the rules of the S(E)IR model.
        call mc_sim(model, history, indices, origins, alpha, lambda, mu, nu, &
          t0, states, 1, tmp_mc_probs)

        !If the epidemic spread, we exit the loop.
        if (count(states /= 'S') > 1) exit
      end do

      !We look for nodes that are not S in the network at time t0 because they
      !are the ones that possibly started the epidemic.
      allocate(non_susceptible(N))
      non_susceptible = [(i, i=1,N)]
      call my_pack(non_susceptible, states/='S')
      gsize = size(non_susceptible)
      allocate(energies(gsize))

      !We loop over each of the non-susceptible nodes.
      do idx = 1, gsize
        node_i = non_susceptible(idx)

        !We apply the DMP(r) algorithm to calculate the marginal probabilities
        !that each node is in a state X (i.e., S, (E), I or R) at time t0,
        !conditional on the fact that node i is the patient zero.
        call dmp(model, restr, history, indices, states, [node_i], alpha, &
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
      !For multiple seeds, we use the average of their ranks as a measure of the
      !efficiency of the algorithm.
      rank = sum(ranks) / dble(gsize*seeds)

      !We count in how many instances the algorithm predicts that at least one
      !of the original seeds has a rank that falls within the bottom 1% of the
      !ranks.
      if (any(ranks <= N/100)) then
        avg_prob = avg_prob + 1.0d0
        avg_prob2 = avg_prob2 + 1.0d0
      end if

      !We accumulate the quantities of interest.
      avg_g = avg_g + dble(gsize)
      avg_rank = avg_rank + rank
      avg_rank2 = avg_rank2 + rank*rank

      !We deallocate the list of non-susceptible nodes and their energies for
      !the next iteration.
      deallocate(energies, non_susceptible)
    end do

    !Normalization of the main magitudes.
    avg_g = avg_g / instances
    avg_prob = avg_prob / instances
    avg_rank = avg_rank / instances
    avg_prob2 = avg_prob2 / instances
    avg_rank2 = avg_rank2 / instances

    !Uncertainties of the main quantities.
    err_prob = sqrt((avg_prob2-avg_prob**2)/dble(instances-1))
    err_rank = sqrt((avg_rank2-avg_rank**2)/dble(instances-1))

    write(*, '(6es26.16)') lambda, avg_g/N, avg_rank, err_rank, avg_prob, err_prob
  end do
end program dmp_lambda
