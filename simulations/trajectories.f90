!Simulation that compares the trajectories given by DMP and MC.
!Author: Adria Meca Montserrat.
!Last modified date: 17/06/22.
program trajectories
  use array_procedures, only : int_llist
  use dmp_algorithms, only : dmp
  use mc_simulations, only : mc_sim
  use network_generation, only : node, PN, RRG
  use random_number_generator, only : ir1279, r1279, setr1279
  use rewiring_algorithms, only : rewiring

  implicit none

  !Variables.
  character(len=1), dimension(:), allocatable :: states
  character(len=9) :: graph, mode, model

  double precision, dimension(:, :), allocatable :: avg_dmp, avg_mc, grid, r, &
    tmp_dmp, tmp_mc
  double precision, dimension(:), allocatable :: pe, pi, pr, ps
  double precision :: alpha, l, lambda, mu, nu, Q

  integer, dimension(:), allocatable :: origins
  integer, dimension(8) :: values
  integer :: c, i, idx, instances, iseed, M, N, points, s, seeds, t, t0, x, y

  logical :: dmpr, rng

  type(int_llist), dimension(:), allocatable :: indices

  type(node), dimension(:), allocatable :: history, network


  !Parameters.
  open(unit=10, file='trajectories.txt')
    read(10, *) rng
    read(10, *) dmpr
    read(10, *) c, l, N
    read(10, *) graph, mode, model
    read(10, *) M, s, points, instances
    read(10, *) seeds, alpha, lambda, mu, nu, t0, Q
  close(10)

  allocate(avg_dmp(t0, 4), avg_mc(t0, 4), grid(points, points), r(N, 2), &
    tmp_dmp(t0, 4), tmp_mc(t0, 4))

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

  !Initialization of the average quantities.
  grid = 0.0d0
  avg_mc = 0.0d0
  avg_dmp = 0.0d0

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

    select case (trim(mode))
      !We compare the DMP and MC trajectories for a given pair (mu, lambda).
      case ('line')
        !DMP.
        call dmp(model, dmpr, history, indices, states, origins, alpha, lambda, &
          mu, nu, t0, ps, pe, pi, pr, tmp_dmp)

        !MC.
        call mc_sim(model, history, indices, origins, alpha, lambda, mu, nu, &
          t0, states, M, tmp_mc)

        avg_mc = avg_mc + tmp_mc
        avg_dmp = avg_dmp + tmp_dmp
      !We compute the difference between the DMP and MC trajectories for each
      !point (mu, lambda) in a given grid.
      case ('grid')
        do y = 1, points
          lambda = (y-1) / dble(points-1)

          do x = 1, points
            mu = (x-1) / dble(points-1)

            !DMP.
            call dmp(model, dmpr, history, indices, states, origins, alpha, &
              lambda, mu, nu, t0, ps, pe, pi, pr, tmp_dmp)

            !MC.
            call mc_sim(model, history, indices, origins, alpha, lambda, mu, &
              nu, t0, states, M, tmp_mc)

            grid(x, y) = grid(x, y) + sum(abs(tmp_dmp(:, s)-tmp_mc(:, s)))
          end do
        end do
    end select
  end do

  !Normalization of the main magnitudes.
  grid = grid / instances
  avg_mc = avg_mc / instances
  avg_dmp = avg_dmp / instances

  select case (trim(mode))
    case ('line')
      !Headers.
      write(*, '(a,a25,2a78)') '#', 't', 'S, I, R (MC)', 'S, I, R (DMP)'

      do t = 1, t0
        write(*, '(i26,6es26.16)') t, avg_mc(t, [1, 3, 4]), avg_dmp(t, [1, 3, 4])
      end do
    case ('grid')
      !Headers.
      write(*, '(a,a25,2a26)') '#', 'mu', 'lambda', 'difference'

      do y = 1, points
        lambda = (y-1) / dble(points-1)
        do x = 1, points
          mu = (x-1) / dble(points-1)
          write(*, '(3es26.16)') mu, lambda, grid(x, y)
        end do
        print*
      end do
  end select
end program trajectories
