!> Simulation that aims to test the efficiency of the DMP(r) algorithm by studying
!> the precision with which it locates the true patient zero(s) in a network at a
!> certain time t as we vary one of the parameters, leaving all the others fixed.
!> Author: Adria Meca Montserrat.
!> Last modified date: 06/08/22.
program instances
  use array_procedures, only: int_llist
  use network_generation, only: node, PN, RRG
  use patient_zero_problem, only: pz_sim
  use random_number_generator, only: ir1279, r1279, setr1279
  use rewiring_algorithms, only: rewiring

  implicit none

  character(len=1), allocatable :: states(:)                                              !>
  character(len=9)              :: graph, model, param, scale                             !>
  double precision, allocatable :: r(:, :)                                                !>
  double precision              :: a, alpha, b, l, lambda, mu, nu, point, Q, rank         !>
  integer,          allocatable :: origins(:), ranks(:)                                   !>
  integer                       :: c, guess, i, idx, ins, iseed, N, p, points, seeds, t0  !>
  integer                       :: gsize(3), values(8)                                    !>
  logical                       :: dmpr, rng                                              !>
  type(int_llist),  allocatable :: indices(:)                                             !>
  type(node),       allocatable :: history(:), network(:)                                 !>

  !> Parameters.
  open(unit=10, file='instances.txt')
    read(10, *) rng
    read(10, *) dmpr
    read(10, *) c, l, N
    read(10, *) graph, model, param
    read(10, *) a, b, points, scale, ins
    read(10, *) seeds, alpha, lambda, mu, nu, t0, Q
  close(10)

  allocate(r(N, 2))

  allocate(history(N), indices(N), network(N), origins(seeds), ranks(seeds), &
    states(N))

  !> We initialize the random number generator.
  if (rng) then
    call date_and_time(values=values)
    iseed = sum(values)
  else
    iseed = 0
  end if

  call setr1279(iseed)

  !> Initialization of the node positions.
  r = sqrt(dble(N)) * reshape([(r1279(), i=1,2*N)], [N, 2])

  do p = 1, points
    !> We restart the random sequence.
    call setr1279(ir1279())

    !> The chosen parameter will vary differently depending on the selected scale.
    if (points > 1) then
      select case (trim(scale))
        case ('log')
          point = exp(log(a) + log(b/a)*(p-1)/(points-1))
        case default
          point = a + (b-a)*(p-1)/(points-1)
      end select
    else
      point = lambda
    end if

    !> We modify the chosen parameter.
    select case (trim(param))
      case ('l')
        l = point
      case ('lambda')
        lambda = point
      case ('t0')
        t0 = int(point)
    end select

    do idx = 1, ins
      !> We create a network of the chosen type.
      select case (trim(graph))
        case ('PN')
          network = PN(N, c, r, l)
        case ('RRG')
          network = RRG(N, c)
      end select

      !> We rewire the network links.
      call rewiring(graph, network, history, indices, l, Q, r, c, t0)

      !> We compute the ranks of the seeds and the size of the epidemic.
      call pz_sim( &
        model,     &
        dmpr,      &
        history,   &
        indices,   &
        seeds,     &
        states,    &
        origins,   &
        ranks,     &
        guess,     &
        gsize,     &
        alpha,     &
        lambda,    &
        mu,        &
        nu,        &
        t0         &
      )

      !> We use the average of the rank of the seeds as a measure of performance.
      rank = sum(ranks) / dble(seeds)

      write(*, '(2es26.16,2i26)') point, rank, gsize(2:3)
    end do
  end do
end program instances
