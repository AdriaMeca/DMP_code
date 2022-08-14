!> Simulation that tests the efficiency with which the DMP-based inference algorithm
!> locates the true patient zero of an epidemic.
!> Author: Adria Meca Montserrat.
!> Last modified date: 14/08/22.
program pz_simulation
  use array_procedures,        only: find
  use derived_types,           only: int_llist, node, prm
  use network_generation,      only: PN, RRG
  use patient_zero_problem,    only: pz_problem
  use random_number_generator, only: r1279, setr1279
  use rewiring_algorithms,     only: rewiring

  implicit none

  integer,          allocatable :: nbrs(:)                                              !> Neighbors of node i at time t.
  integer,          allocatable :: non_S(:)                                             !> Non-susceptible nodes.
  integer                       :: alti, c, i, idx, instances, iseed, N, r0, t, t0, t1  !>
  integer                       :: epi_size                                             !> Epidemic size.
  integer                       :: patient_zeros(1)                                     !>
  integer                       :: ppz, tpz                                             !>
  integer                       :: values(8)                                            !>
  character(len=1), allocatable :: states(:, :)                                         !> Distribution of node states.
  character(len=9)              :: graph, model                                         !>
  double precision, allocatable :: energies(:)                                          !> Node energies.
  double precision, allocatable :: r(:, :)                                              !> Node positions.
  double precision, parameter   :: small = 1.0d-300                                     !>
  double precision              :: alpha, e, l, lambda, mu, nu, Q                       !>
  logical                       :: restricted, randomized                               !>
  type(int_llist),  allocatable :: indices(:)                                           !> Active links throughout the simulation.
  type(node),       allocatable :: history(:)                                           !> Rewiring history.
  type(node),       allocatable :: network(:)                                           !> Original network.
  type(prm)                     :: epi_params                                           !> Epidemiological parameters.

  !> Parameters.
  open(unit=10, file='pz_simulation'//'.txt')
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
  allocate(history(N), indices(N), network(N))

  !> We initialize the epidemiological parameters.
  epi_params = prm(t0, alpha, lambda, mu, nu)

  !> We generate a random seed.
  call date_and_time(values=values); iseed = merge(sum(values), 0, randomized)

  !> We initialize the random number generator.
  call setr1279(iseed)

  !> We initialize the node positions.
  r = sqrt(dble(N)) * reshape([(r1279(), i=1,2*N)], [N, 2])

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
    do t = 0, t0
      if (states(t, ppz) == 'I') then
        t1 = t
        exit
      end if
    end do

    !> Monte Carlo history.
    open(unit=15, file='pz_simulation'//'_mc.dat', position='append')
      do t = 0, t0
        do i = 1, N
          allocate(nbrs(size(indices(i)%time(t)%array)))
          nbrs = history(i)%neighbors(indices(i)%time(t)%array)

          if (t == 0) then
            write(15, '(3i10,a10,999i10,:)') idx, t, i, states(t, i), nbrs
          else if (states(t-1, i) /= states(t, i)) then
            write(15, '(3i10,a10,999i10,:)') idx, t, i, states(t, i), nbrs
          else if (different_neighbors(history, indices, i, t)) then
            write(15, '(3i10,a10,999i10,:)') idx, t, i, states(t, i), nbrs
          end if

          deallocate(nbrs)
        end do
      end do
    close(15)

    !> Energies.
    open(unit=20, file='pz_simulation'//'_energies.dat', position='append')
      epi_size = size(non_S)
      do alti = 1, epi_size
        e = energies(alti); i = non_S(alti)

        if (e < 1/small) write(20, '(3i10,es26.16)') idx, epi_size, i, e
      end do
    close(20)

    !> Output.
    write(*, '(10i10)')                                                             &
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

contains
  !> Indicates whether the neighbors of node i have changed since the previous time step.
  function different_neighbors(history, indices, i, t)
    integer,         intent(in) :: i                    !> Node of interest.
    integer,         intent(in) :: t                    !> Time step.
    integer                     :: altk, isize, k       !>
    logical                     :: different_neighbors  !> Output.
    type(int_llist), intent(in) :: indices(:)           !> Active links throughout the simulation.
    type(node),      intent(in) :: history(:)           !> Rewiring history.

    !> Number of active links of node i at time t.
    isize = size(indices(i)%time(t)%array)

    different_neighbors = .false.
    if (isize /= size(indices(i)%time(t-1)%array)) then
      different_neighbors = .true.
    else
      do altk = 1, isize
        k = history(i)%neighbors(indices(i)%time(t)%array(altk))

        if (all(history(i)%neighbors(indices(i)%time(t-1)%array) /= k)) then
          different_neighbors = .true.
          exit
        end if
      end do
    end if
  end function different_neighbors
end program pz_simulation
