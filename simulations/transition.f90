!Simulation that generates the phase diagram of the SIR model for a given set
!of epidemiological parameters.
!Author: Adria Meca Montserrat.
!Last modified date: 17/06/22.
program transition
  use array_procedures, only : int_llist
  use dmp_algorithms, only : dmp
  use network_generation, only : node, PN, RRG
  use random_number_generator, only : ir1279, r1279, setr1279
  use rewiring_algorithms, only : rewiring

  implicit none

  !Variables.
  character(len=1), dimension(:), allocatable :: states
  character(len=9) :: graph, model

  double precision, dimension(:, :, :), allocatable :: grid
  double precision, dimension(:, :), allocatable :: r, tmp_dmp_probs
  double precision, dimension(:), allocatable :: pe, pi, pr, ps
  double precision :: alpha, l, lambda, mu, nu, Q

  integer, dimension(:), allocatable :: origins
  integer, dimension(8) :: values
  integer :: c, i, idx, instances, iseed, N, points, seeds, t, t0, x, y

  logical :: dmpr, rng

  type(int_llist), dimension(:), allocatable :: indices

  type(node), dimension(:), allocatable :: history, network


  !Parameters.
  open(unit=10, file='transition.txt')
    read(10, *) rng
    read(10, *) dmpr
    read(10, *) c, l, N
    read(10, *) graph, model
    read(10, *) points, instances
    read(10, *) seeds, alpha, lambda, mu, nu, t0, Q
  close(10)

  allocate(grid(t0, points, points))

  allocate(r(N, 2), tmp_dmp_probs(t0, 4))

  allocate(history(N), indices(N), network(N), origins(seeds), pe(N), pi(N), &
    pr(N), ps(N), states(N))

  !We initialize the random number generator.
  if (rng) then
    call date_and_time(values=values)
    iseed = sum(values)
  else
    iseed = 0
  end if
  call setr1279(iseed)

  grid = 0.0d0
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

    do y = 1, points
      lambda = (y-1) / dble(points-1)

      do x = 1, points
        mu = (x-1) / dble(points-1)

        !DMP.
        call dmp(model, dmpr, history, indices, states, origins, alpha, &
          lambda, mu, nu, t0, ps, pe, pi, pr, tmp_dmp_probs)

        grid(:, x, y) = grid(:, x, y) + sum(tmp_dmp_probs(:, 2:4), dim=2)
      end do
    end do
  end do

  !Normalization of the grid.
  grid = grid / instances

  do t = 1, t0
    !Headers.
    write(*, '(a,i0)') '# t=', t
    write(*, '(a,a25,2a26)') '#', 'mu', 'lambda', 'Pq!=S'

    do y = 1, points
      lambda = (y-1) / dble(points-1)

      do x = 1, points
        mu = (x-1) / dble(points-1)

        write(*, '(3es26.16)') mu, lambda, grid(t, x, y)
      end do
      print*
    end do
    write(*, '(a)') achar(10)
  end do
end program transition
