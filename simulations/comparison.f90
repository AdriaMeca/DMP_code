!Simulation that compares the marginal probabilities given by DMP and MC.
!Author: Adria Meca Montserrat.
!Last modified date: 16/06/22.
program comparison
  use array_procedures, only : idx_insertion_sort, int_llist
  use dmp_algorithms, only : dmp
  use mc_simulations, only : mc_sim
  use network_generation, only : node, PN, RRG
  use random_number_generator, only : ir1279, r1279, setr1279
  use rewiring_algorithms, only : rewiring

  implicit none

  !Variables.
  character(len=1), dimension(:), allocatable :: states
  character(len=9) :: graph, model

  double precision, dimension(:, :), allocatable :: probabilities, r, tmp_dmp, &
    tmp_mc
  double precision, dimension(:), allocatable :: pe_dmp, pi_dmp, pi_mc, pr_dmp, &
    pr_mc, ps_dmp, ps_mc
  double precision, dimension(6) :: averages, averages2, errors
  double precision :: alpha, l, lambda, mu, nu, Q

  integer, dimension(:), allocatable :: origins, p1, p2, p3
  integer, dimension(8) :: values
  integer :: alti, c, i, idx, instances, iseed, M, N, seeds, t0

  logical :: dmpr, rng

  type(int_llist), dimension(:), allocatable :: indices

  type(node), dimension(:), allocatable :: history, network


  !Parameters.
  open(unit=10, file='comparison.txt')
    read(10, *) rng
    read(10, *) dmpr
    read(10, *) c, l, N
    read(10, *) graph, model
    read(10, *) M, instances
    read(10, *) seeds, alpha, lambda, mu, nu, t0, Q
  close(10)

  allocate(probabilities(6, instances*N), r(N, 2), tmp_dmp(t0, 4), tmp_mc(t0, 4))

  allocate(history(N), indices(N), network(N), origins(seeds), pe_dmp(N), &
    pi_dmp(N), pi_mc(N), pr_dmp(N), pr_mc(N), ps_dmp(N), ps_mc(N), &
    p1(instances*N), p2(instances*N), p3(instances*N), states(N))

  !We initialize the random number generator.
  if (rng) then
    call date_and_time(values=values)
    iseed = sum(values)
  else
    iseed = 0
  end if
  call setr1279(iseed)

  do idx = 1, instances
    !We restart the random sequence.
    call setr1279(ir1279())

    !Initialization of the node positions.
    r = sqrt(dble(N)) * reshape([(r1279(), i=1,2*N)], [N, 2])

    !We choose the origins at random.
    origins = [(1 + mod(int(N*r1279()), N), i=1,seeds)]

    !We create a network of the chosen type.
    select case (trim(graph))
      case ('PN')
        network = PN(N, c, r, l)
      case ('RRG')
        network = RRG(N, c)
    end select

    !We rewire the connections of the network.
    call rewiring(graph, network, history, indices, c, t0, r, l, Q)

    !MC.
    ps_mc = 0.0d0
    pi_mc = 0.0d0
    pr_mc = 0.0d0
    do alti = 1, M
      call mc_sim(model, history, indices, origins, alpha, lambda, mu, nu, t0, &
        states, 1, tmp_mc)

      do i = 1, N
        select case (states(i))
          case ('S')
            ps_mc(i) = ps_mc(i) + 1.0d0
          case ('I')
            pi_mc(i) = pi_mc(i) + 1.0d0
          case ('R')
            pr_mc(i) = pr_mc(i) + 1.0d0
        end select
      end do
    end do
    ps_mc = ps_mc / M
    pi_mc = pi_mc / M
    pr_mc = pr_mc / M

    !DMP.
    call dmp(model, dmpr, history, indices, states, origins, alpha, lambda, &
      mu, nu, t0, ps_dmp, pe_dmp, pi_dmp, pr_dmp, tmp_dmp)

    !We save the probabilities of all the nodes.
    probabilities(1, 1+N*(idx-1):N*idx) = ps_mc
    probabilities(3, 1+N*(idx-1):N*idx) = pi_mc
    probabilities(5, 1+N*(idx-1):N*idx) = pr_mc
    probabilities(2, 1+N*(idx-1):N*idx) = ps_dmp
    probabilities(4, 1+N*(idx-1):N*idx) = pi_dmp
    probabilities(6, 1+N*(idx-1):N*idx) = pr_dmp
  end do

  !We compute the permutations of the indices that sort the MC probabilities in
  !ascending order (p1: S, p2: I, p3: R).
  p1 = idx_insertion_sort(probabilities(1, :))
  p2 = idx_insertion_sort(probabilities(3, :))
  p3 = idx_insertion_sort(probabilities(5, :))

  !We order the probabilities using the above indices.
  probabilities(1, :) = probabilities(1, p1)
  probabilities(2, :) = probabilities(2, p1)
  probabilities(3, :) = probabilities(3, p2)
  probabilities(4, :) = probabilities(4, p2)
  probabilities(5, :) = probabilities(5, p3)
  probabilities(6, :) = probabilities(6, p3)

  !Headers.
  write(*, '(a,a51,2a52)') '#', 'PS (MC & DMP)', 'PI (MC & DMP)', 'PR (MC & DMP)'

  !Probabilities of all nodes in a given number of instances.
  do i = 1, instances*N
    write(*, '(6es26.16)') probabilities(:, i)
  end do
  write(*, '(a)') achar(10)

  !Headers.
  write(*, '(a,a25,5a26)') '#', 'PS', 'dPS', 'PI', 'dPI', 'PR', 'dPR'

  !Average DMP probabilities.
  averages = 0.0d0
  averages2 = 0.0d0
  do i = 1, instances*N
    averages = averages + probabilities(:, i)
    averages2 = averages2 + probabilities(:, i)*probabilities(:, i)

    if (mod(i, instances) == 0) then
      averages = averages / instances
      averages2 = averages2 / instances

      errors = sqrt(abs(averages2-averages**2)/dble(instances-1))

      write(*, '(6es26.16)') averages(2), errors(2), averages(4), errors(4), &
        averages(6), errors(6)

      averages = 0.0d0
      averages2 = 0.0d0
    end if
  end do
end program comparison
