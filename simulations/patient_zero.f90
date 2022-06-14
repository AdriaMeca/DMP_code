!Simulation that aims to test the efficiency of the DMP(r) algorithm by studying
!the precision with which it locates the true patient zero(s) in a network at a
!certain time t as we vary one of the parameters, leaving all the others fixed.
!Author: Adria Meca Montserrat.
!Last modified date: 14/06/22.
program patient_zero
  use array_procedures, only : int_llist
  use network_generation, only : node, PN, RRG
  use patient_zero_problem, only : pz_sim
  use random_number_generator, only : ir1279, r1279, setr1279
  use rewiring_algorithms, only : rewiring

  implicit none

  !Variables.
  character(len=1), dimension(:), allocatable :: states
  character(len=9) :: graph, model, param, scale

  double precision, dimension(:, :), allocatable :: r
  double precision, dimension(3) :: avg_g
  double precision :: a, alpha, avg_norm, avg_norm2, avg_prob, avg_rank, &
    avg_rank2, avg_true, b, err_norm, err_rank, l, lambda, mu, norm, nu, &
    point, Q, rank

  integer, dimension(:), allocatable :: origins, ranks
  integer, dimension(8) :: values
  integer, dimension(3) :: gsize
  integer :: c, i, idx, instances, iseed, N, p, points, seeds, t0

  logical :: dmpr, rng

  type(int_llist), dimension(:), allocatable :: indices

  type(node), dimension(:), allocatable :: history, network


  !Parameters.
  open(unit=10, file='patient_zero.txt')
    read(10, *) rng
    read(10, *) dmpr
    read(10, *) c, l, N
    read(10, *) graph, model, param
    read(10, *) a, b, points, scale, instances
    read(10, *) seeds, alpha, lambda, mu, nu, t0, Q
  close(10)

  allocate(history(N), indices(N), network(N), r(N, 2), states(N))
  allocate(origins(seeds), ranks(seeds))

  !We initialize the random number generator.
  if (rng) then
    call date_and_time(values=values)
    iseed = sum(values)
  else
    iseed = 0
  end if
  call setr1279(iseed)

  !Initialization of the node positions.
  r = sqrt(dble(N)) * reshape([(r1279(), i=1,2*N)], [N, 2])

  !Headers.
  write(*, '(a,i0)') '# N=', N
  write(*, '(a,a25,8a26)') '#', 'x', 'I', 'R', 'r0', 'dr0', 'r0/G', 'd(r0/G)', &
    'P(r0=0)', 'P(r0<=r)'

  do p = 1, points
    !We restart the random sequence.
    call setr1279(ir1279())

    !The chosen parameter will vary differently depending on the selected scale.
    select case (trim(scale))
      case ('log')
        point = exp(log(a) + log(b/a)*(p-1)/(points-1))
      case default
        point = a + (b-a)*(p-1)/(points-1)
    end select

    !We modify the chosen parameter.
    select case (trim(param))
      case ('l')
        l = point
      case ('lambda')
        lambda = point
      case ('t0')
        t0 = int(point)
    end select

    !Initialization of the average quantities.
    avg_g = 0.0d0
    avg_norm = 0.0d0
    avg_prob = 0.0d0
    avg_rank = 0.0d0
    avg_true = 0.0d0
    avg_rank2 = 0.0d0
    avg_norm2 = 0.0d0

    do idx = 1, instances
      !We create a network of the chosen type.
      select case (trim(graph))
        case ('PN')
          network = PN(N, c, r, l)
        case ('RRG')
          network = RRG(N, c)
      end select

      !We rewire the connections of the network.
      call rewiring(graph, network, history, indices, c, t0, r, l, Q)

      !We compute the ranks of the seeds and the size of the epidemic.
      call pz_sim(model, dmpr, history, indices, seeds, states, origins, ranks, &
        gsize, alpha, lambda, mu, nu, t0)

      !We use the average of the rank of the seeds as a measure of performance.
      rank = sum(ranks) / dble(seeds)
      norm = rank / sum(gsize)

      !We count in how many instances the algorithm finds the true patient zero.
      if (any(ranks == 0)) then
        avg_true = avg_true + 1.0d0
      end if

      !We count in how many instances the algorithm predicts that at least one
      !of the original seeds has a rank that falls within the bottom 1% of the
      !ranks. If N is too large, we set the upper limit to 100.
      if (any(ranks <= min(100, N/100))) then
        avg_prob = avg_prob + 1.0d0
      end if

      !We accumulate the quantities of interest.
      avg_g = avg_g + dble(gsize)
      avg_norm = avg_norm + norm
      avg_rank = avg_rank + rank
      avg_rank2 = avg_rank2 + rank*rank
      avg_norm2 = avg_norm2 + norm*norm
    end do

    !Normalization of the main magitudes.
    avg_g = avg_g / instances
    avg_norm = avg_norm / instances
    avg_prob = avg_prob / instances
    avg_rank = avg_rank / instances
    avg_true = avg_true / instances
    avg_rank2 = avg_rank2 / instances
    avg_norm2 = avg_norm2 / instances

    !Uncertainty associated with the (normalized) rank.
    err_norm = sqrt(abs(avg_norm2-avg_norm**2)/dble(instances-1))
    err_rank = sqrt(abs(avg_rank2-avg_rank**2)/dble(instances-1))

    write(*, '(9es26.16)') point, avg_g(2:3), avg_rank, err_rank, avg_norm, &
      err_norm, avg_true, avg_prob
  end do
end program patient_zero
