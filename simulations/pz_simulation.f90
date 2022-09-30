!> Simulation that tests the efficiency with which the DMP-based inference algorithm
!> locates the true patient zero of an epidemic.
!> Author: Adria Meca Montserrat.
!> Last modified date: 30/09/22.
program pz_simulation
  use array_procedures,        only: find
  use derived_types,           only: int_list, node, prm
  use network_generation,      only: PN, RRG
  use patient_zero_problem,    only: pz_problem
  use random_number_generator, only: r1279, setr1279
  use rewiring_algorithms,     only: rewiring

  implicit none

  integer,          allocatable :: non_S(:)                                  !> Non-susceptible nodes.
  integer                       :: alti, c, i, idx, instances, iseed, N, r0  !>
  integer                       :: patient_zeros(1)                          !>
  integer                       :: ppz, tpz                                  !>
  integer                       :: t0, t1, t2                                !>
  integer                       :: values(8)                                 !>
  character(len=1), allocatable :: states(:, :)                              !> Distribution of node states.
  character(len=9)              :: graph, model                              !>
  double precision, allocatable :: energies(:)                               !> Node energies.
  double precision, allocatable :: r(:, :)                                   !> Node positions.
  double precision, parameter   :: small = 1.0d-300                          !>
  double precision              :: alpha, e, l, lambda, mu, nu, Q            !>
  logical                       :: restricted, randomized                    !>
  type(int_list),   allocatable :: indices(:, :)                             !> Active links throughout the simulation.
  type(node),       allocatable :: history(:)                                !> Rewiring history.
  type(node),       allocatable :: network(:)                                !> Original network.
  type(prm)                     :: epi_params                                !> Epidemiological parameters.

  !> Parameters.
  open(unit=10, file='parameters'//'.txt')
    read(10, *)                            !>
    read(10, *) instances                  !> Number of instances.
    read(10, *) randomized                 !> T: Random seed; F: Fixed seed.
    read(10, *) restricted                 !> T: DMPr; F: DMP.

    read(10, *) graph                      !> Network type.
    read(10, *) c, l, N, Q                 !> Network properties.

    read(10, *) model                      !> Epidemiological model.
    read(10, *) alpha, lambda, mu, nu, t0  !> Epidemiological parameters.
  close(10)

  allocate(r(N, 2))
  allocate(states(0:t0+1, N))
  allocate(history(N), indices(N, 0:t0), network(N))

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
        r = sqrt(dble(N)) * reshape([(r1279(), i=1,2*N)], [N, 2])

        network = PN(N, c, r, l)
      case ('RRG')
        network = RRG(N, c)
    end select

    !> We rewire the network links.
    call rewiring(graph, network, history, indices, l, Q, r, c, t0)

    !> We create an epidemic and rank the non-susceptible nodes by increasing value
    !> of their energy.
    call pz_problem( &
      model,         &
      restricted,    &
      history,       &
      indices,       &
      patient_zeros, &
      epi_params,    &
      states,        &
      energies,      &
      non_S          &
    )

    !> Predicted and true patient zeros, respectively.
    ppz = non_S(1)
    tpz = patient_zeros(1)

    !> We compute the rank of the true patient zero.
    r0 = find(non_S, tpz) - 1

    !> We compute the time t1 at which the predicted patient zero got infected.
    t1 = find(states(:, ppz), 'I') - 1

    !> We compute the time t2 at which the true patient zero recovered.
    t2 = find(states(:, tpz), 'R') - 1

    !> Energies.
    open(unit=15, file='pz_simulation'//'_energies.dat', position='append')
      do alti = 1, min(size(non_S), 100)
        e = energies(alti); i = non_S(alti)

        if (e < 1/small) write(15, '(2i10,es26.16)') idx, i, e
      end do
    close(15)

    !> Output.
    write(*, '(12i10)')                                                        &
      tpz,                                                                     &
      ppz,                                                                     &
      r0,                                                                      &
      t1,                                                                      &
      t2,                                                                      &
      merge(1, 0, any(history(ppz)%neighbors(indices(ppz, t1)%array) == tpz)), &
      count(states(t0, :) == 'I'),                                             &
      count(states(t0, :) == 'R'),                                             &
      count(states(t1, :) == 'I'),                                             &
      count(states(t1, :) == 'R'),                                             &
      count(states(t2, :) == 'I'),                                             &
      count(states(t2, :) == 'R')
  end do
end program pz_simulation
