!> Simulation that computes the DMP and MC marginal probabilities that each node
!> in a time-varying network is in states S, E, I or R at times t <= t0.
!> Author: Adria Meca Montserrat.
!> Last modified date: 01/10/22.
!> Last reviewed date: 01/10/22.
program probabilities
  use derived_types,           only: int_list, node, prm
  use dmp_algorithms,          only: dmp_alg
  use mc_simulations,          only: mc_sim
  use network_generation,      only: PN, RRG
  use random_number_generator, only: r1279, setr1279
  use rewiring_algorithms,     only: rewiring

  implicit none

  integer,          allocatable :: nodes(:)                                        !>
  integer                       :: alti, c, i, idx, instances, iseed, M, N, t, t0  !>
  integer                       :: patient_zeros(1)                                !>
  integer                       :: values(8)                                       !>
  character(len=1), allocatable :: states(:, :)                                    !> Distribution of node states.
  character(len=9)              :: graph, model                                    !>
  double precision, allocatable :: p_dmp(:, :, :), p_mc(:, :, :)                   !> Marginal probabilities.
  double precision, allocatable :: r(:, :)                                         !> Node positions.
  double precision              :: alpha, l, lambda, mu, nu, Q                     !>
  logical                       :: randomized                                      !>
  type(int_list),   allocatable :: indices(:, :)                                   !> Active links throughout the simulation.
  type(node),       allocatable :: history(:)                                      !> Rewiring history.
  type(node),       allocatable :: network(:)                                      !> Original network.
  type(prm)                     :: epi_params                                      !> Epidemiological parameters.

  !> Parameters.
  open(unit=10, file='parameters'//'.txt')
    read(10, *) M                          !> Number of MC repetitions.
    read(10, *) instances                  !> Number of instances.
    read(10, *) randomized                 !> T: Random seed; F: Fixed seed.
    read(10, *)                            !>

    read(10, *) graph                      !> Network type.
    read(10, *) c, l, N, Q                 !> Network properties.

    read(10, *) model                      !> Epidemiological model.
    read(10, *) alpha, lambda, mu, nu, t0  !> Epidemiological parameters.
  close(10)

  allocate(r(2, N))
  allocate(history(N), network(N), nodes(N))
  allocate(indices(N, 0:t0), states(N, 0:t0+1))
  allocate(p_dmp(N, 0:t0, 4), p_mc(N, 0:t0, 4))

  !> Nodes participating in the iteration of the DMP algorithm.
  nodes = [(i, i=1,N)]

  !> We initialize the epidemiological parameters.
  epi_params = prm(t0, alpha, lambda, mu, nu)

  !> We generate a random seed.
  call date_and_time(values=values); iseed = merge(sum(values), 0, randomized)

  !> We initialize the random number generator.
  call setr1279(iseed)

  do idx = 1, instances
    !> We create a network of the chosen type.
    select case (trim(graph))
      case ('PN')
        !> We initialize the node positions.
        r = sqrt(dble(N)) * reshape([(r1279(), i=1,2*N)], [2, N])

        network = PN(N, c, r, l)
      case ('RRG')
        network = RRG(N, c)
    end select

    !> We rewire the network links.
    call rewiring(graph, network, history, indices, l, Q, r, c, t0)

    !> We choose the patient zero(s) uniformly at random.
    patient_zeros = [(1 + mod(int(N*r1279()), N), i=1,size(patient_zeros))]

    !> We calculate the MC marginal probabilities that each node in a time-varying
    !> network is in states S, E, I or R at times t <= t0.
    p_mc = 0.0d0
    do alti = 1, M
      call mc_sim(model, history, indices, patient_zeros, epi_params, states)

      do t = 0, t0
        do i = 1, N
          select case (states(i, t))
            case ('S')
              p_mc(i, t, 1) = p_mc(i, t, 1) + 1.0d0
            case ('E')
              p_mc(i, t, 2) = p_mc(i, t, 2) + 1.0d0
            case ('I')
              p_mc(i, t, 3) = p_mc(i, t, 3) + 1.0d0
            case ('R')
              p_mc(i, t, 4) = p_mc(i, t, 4) + 1.0d0
          end select
        end do
      end do
    end do
    p_mc = p_mc / M

    !> We calculate the DMP marginal probabilities that each node in a time-varying
    !> network is in states S, E, I or R at times t <= t0.
    call dmp_alg(model, history, indices, N, nodes, patient_zeros, epi_params, p_dmp)

    !> We print the marginal probabilities given by both methods.
    do t = 0, t0
      do i = 1, N
        print *, t, p_mc(i, t, [1, 3, 4]), p_dmp(i, t, [1, 3, 4])
      end do
    end do
  end do
end program probabilities
