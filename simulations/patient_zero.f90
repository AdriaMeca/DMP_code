!Simulation that aims to test the efficiency of the DMP(r) algorithm by studying
!the precision with which it locates the true patient zero(s) in a network at a
!certain time t as we vary one of the parameters, leaving all the others fixed.
!Author: Adria Meca Montserrat.
!Last modified date: 29/05/22.
program patient_zero
  use array_procedures, only : int_llist
  use network_generation, only : node, PN, RRG
  use patient_zero_problem, only : pz_sim
  use random_number_generator, only : r1279, setr1279
  use rewiring_algorithms, only : rewiring

  implicit none

  !Variables.
  character(len=9) :: graph, model, param, scale

  double precision, dimension(:, :), allocatable :: rij
  double precision :: a, alpha, avg_g, avg_prob, avg_prob2, avg_rank, avg_rank2, &
    b, err_prob, err_rank, l, lambda, mu, nu, point, Q, rank

  integer, dimension(:), allocatable :: ranks
  integer :: c, gsize, i, instances, N, p, points, seeds, t0

  logical :: dmpr

  type(int_llist), dimension(:), allocatable :: indices

  type(node), dimension(:), allocatable :: history, network


  !Parameters.
  open(unit=10, file='parameters.txt')
    read(10, *) dmpr
    read(10, *) c, l, N
    read(10, *) graph, model, param
    read(10, *) a, b, points, scale, instances
    read(10, *) seeds, alpha, lambda, mu, nu, t0, Q
  close(10)

  !We allocate the array containing the ranks of the seeds.
  allocate(ranks(seeds))

  !We initialize the random number generator.
  call setr1279(1)

  !Initialization of the node positions.
  if (trim(graph) == 'PN') then
    allocate(rij(N, 2))
    rij = sqrt(dble(N)) * reshape([(r1279(), i=1,2*N)], [N, 2])
  end if

  do p = 1, points
    !The chosen parameter will vary differently depending on the selected scale.
    select case (trim(scale))
      case ('log')
        point = exp(log(a) + log(b/a)*(p-1)/(points-1))
      case default
        point = a + (b-a)*(p-1)/(points-1)
    end select

    !We modify the chosen parameter.
    select case (trim(param))
      case ('N')
        N = int(point)
      case ('l')
        l = point
      case ('lambda')
        lambda = point
      case ('t0')
        t0 = int(point)
    end select

    !If we choose to vary N or l, we have to rebuild the network, history and
    !indices for each new value of the parameter.
    if ((trim(param) == 'N').and.allocated(network)) then
      deallocate(history, indices, network, rij)

      !Additionally, we have to reposition each node if N changes.
      allocate(rij(N, 2))
      rij = sqrt(dble(N)) * reshape([(r1279(), i=1,2*N)], [N, 2])
    else if ((trim(param) == 'l').and.allocated(network)) then
      deallocate(history, indices, network)
    end if

    if (.not.allocated(network)) then
      allocate(history(N), indices(N), network(N))

      !We create a network of the chosen type.
      select case (trim(graph))
        case ('PN')
          network = PN(N, c, rij, l)
        case ('RRG')
          network = RRG(N, c)
      end select

      !We rewire the network to obtain the history and indices.
      call rewiring(graph, network, history, indices, c, t0, Q)
    end if

    !Initialization of the average quantities.
    avg_g = 0.0d0
    avg_prob = 0.0d0
    avg_rank = 0.0d0
    avg_prob2 = 0.0d0
    avg_rank2 = 0.0d0

    do i = 1, instances
      !We compute the ranks of the seeds and the size of the epidemic.
      call pz_sim(model, dmpr, history, indices, seeds, ranks, gsize, alpha, &
        lambda, mu, nu, t0)

      !For multiple seeds, we use the average of their ranks as a measure of the
      !efficiency of the algorithm.
      rank = sum(ranks) / dble(gsize*seeds)

      !We count in how many instances the algorithm predicts that at least one
      !of the original seeds has a rank that falls within the bottom 1% of the
      !ranks. If N is too large, we set the upper limit to 100.
      if (any(ranks <= min(100, N/100))) then
        avg_prob = avg_prob + 1.0d0
        avg_prob2 = avg_prob2 + 1.0d0
      end if

      !We accumulate the quantities of interest.
      avg_g = avg_g + dble(gsize)
      avg_rank = avg_rank + rank
      avg_rank2 = avg_rank2 + rank*rank
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

    write(*, '(6es26.16)') point, avg_g/N, avg_rank, err_rank, avg_prob, err_prob
  end do
end program patient_zero
