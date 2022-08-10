!> Simulation that tests the efficiency with which the DMP-based inference algorithm
!> locates the true patient zero of an epidemic.
!> Author: Adria Meca Montserrat.
!> Last modified date: 10/08/22.
program pz_simulation
  use array_procedures,        only: find
  use derived_types,           only: int_llist, node, prm
  use network_generation,      only: PN, RRG
  use patient_zero_problem,    only: pz_problem
  use random_number_generator, only: r1279, setr1279
  use rewiring_algorithms,     only: rewiring

  implicit none

  integer,          allocatable :: non_S(:)                                       !> Non-susceptible nodes.
  integer                       :: c, i, idx, instances, iseed, N, r0, t, t0, t1  !>
  integer                       :: patient_zeros(1)                               !>
  integer                       :: ppz, tpz                                       !> Predicted and true patient zeros, respectively.
  integer                       :: values(8)                                      !> Dummy array.
  character(len=1), allocatable :: states(:, :)                                   !> Distribution of node states.
  character(len=9)              :: graph, model                                   !>
  double precision, allocatable :: energies(:)                                    !> Node energies.
  double precision, allocatable :: r(:, :)                                        !> Node positions.
  double precision              :: alpha, l, lambda, mu, nu, Q                    !>
  logical                       :: restricted, randomized                         !>
  type(int_llist),  allocatable :: indices(:)                                     !> Active links throughout the simulation.
  type(node),       allocatable :: history(:)                                     !> Rewiring history.
  type(node),       allocatable :: network(:)                                     !> Original network.
  type(prm)                     :: epi_params                                     !> Epidemiological parameters.

  !> Parameters.
  open(unit=10, file='pz_simulation.txt')
    read(10, *) instances                  !> Number of instances.
    read(10, *) randomized                 !> T: Random seed; F: Fixed seed.
    read(10, *) restricted                 !> T: DMPr; F: DMP.

    read(10, *) graph                      !> Network type.
    read(10, *) c, l, N, Q                 !> Network Properties.

    read(10, *) model                      !> Epidemiological model.
    read(10, *) alpha, lambda, mu, nu, t0  !> Epidemiological parameters.
  close(10)

  allocate(r(N, 2))
  allocate(states(0:t0+1, N))
  allocate(history(N), indices(N), network(N))

  epi_params = prm(t0, alpha, lambda, mu, nu)

  !> We generate a random seed.
  if (randomized) then
    call date_and_time(values=values)
    iseed = sum(values)
  else
    iseed = 0
  end if

  !> We initialize the random number generator.
  call setr1279(iseed)

  !> We initialize the node positions.
  r = sqrt(dble(N)) * reshape([(r1279(), i=1,2*N)], [N, 2])

  !> Header.
  write(*, '(a,a5,9a7)') '# ', &
    'TPZ',                     &
    'PPZ',                     &
    'r0',                      &
    't1',                      &
    'c1',                      &
    'c2',                      &
    'I(t0)',                   &
    'R(t0)',                   &
    'I(t1)',                   &
    'R(t1)'

  do idx = 1, instances
    !> We create a network of the chosen type.
    select case (trim(graph))
      case ('PN')
        network = PN(N, c, r, l)
      case ('RRG')
        network = RRG(N, c)
    end select

    !> We rewire the network links.
    call rewiring(graph, network, history, indices, l, Q, r, c, t0)

    !> We create an epidemic and rank the non-susceptible nodes by their energy.
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

    !> We get the true patient zero.
    tpz = patient_zeros(1)

    !> We get the predicted patient zero.
    ppz = non_S(1)

    !> We compute the rank of the true patient zero.
    r0 = find(non_S, tpz) - 1

    !> We compute the time t1 at which the predicted patient zero got infected.
    do t = 0, t0
      if (states(t, ppz) == 'I') then
        t1 = t
        exit
      end if
    end do

    !> Output.
    write(*, '(10i7)')                                                              &
      tpz,                                                                          &
      ppz,                                                                          &
      r0,                                                                           &
      t1,                                                                           &
      merge(1, 0, any(history(ppz)%neighbors(indices(ppz)%time(t1)%array) == tpz)), &
      merge(1, 0, states(t1, tpz) == 'I'),                                          &
      count(states(t0, :) == 'I'),                                                  &
      count(states(t0, :) == 'R'),                                                  &
      count(states(t1, :) == 'I'),                                                  &
      count(states(t1, :) == 'R')
  end do
end program pz_simulation
