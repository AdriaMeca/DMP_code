!> Procedures that simulate the patient-zero problem.
!> Author: Adria Meca Montserrat.
!> Last modified date: 03/09/22.
module patient_zero_problem
  use array_procedures,        only: my_pack, quicksort
  use derived_types,           only: int_llist, node, prm
  use dmp_algorithms,          only: dmp_alg
  use mc_simulations,          only: mc_sim
  use random_number_generator, only: r1279

  implicit none

  private

  public :: pz_problem

contains
  !> Tries to locate the patient zero(s) of an epidemic.
  subroutine pz_problem( &
      model,             &
      restricted,        &
      history,           &
      indices,           &
      patient_zeros,     &
      epi_params,        &
      states,            &
      energies,          &
      non_S              &
  )
    integer,          allocatable, intent(out) :: non_S(:)                 !> Non-susceptible nodes.
    integer,                       intent(out) :: patient_zeros(:)         !>
    integer,          allocatable              :: nodes(:)                 !>
    integer                                    :: epi_size                 !> Epidemic size.
    integer                                    :: i, idx, j, N, t0, total  !>
    character(len=*),              intent(in)  :: model                    !> Epidemiological model.
    character(len=1),              intent(out) :: states(0:, :)            !> Distribution of node states.
    double precision, allocatable, intent(out) :: energies(:)              !> Node energies.
    double precision, allocatable              :: p(:, :, :)               !> DMP marginal probabilities.
    double precision, parameter                :: small = 1.0d-300         !>
    double precision                           :: joint                    !> P(O|i).
    logical,                       intent(in)  :: restricted               !> T: DMPr; F: DMP.
    type(int_llist),               intent(in)  :: indices(:)               !> Active links throughout the simulation.
    type(node),                    intent(in)  :: history(:)               !> Rewiring history.
    type(prm),                     intent(in)  :: epi_params               !> Epidemiological parameters.

    !> Number of nodes.
    N = size(history)

    !> Local observation time.
    t0 = epi_params%t0

    allocate(nodes(N))
    allocate(non_S(N))
    allocate(p(4, 0:t0, N))

    do while (.true.)
      !> We choose the patient zero(s) uniformly at random.
      patient_zeros = [(1 + mod(int(N*r1279()), N), i=1,size(patient_zeros))]

      !> We create an epidemic by simulating the S(E)IR rules with MC up to t0.
      call mc_sim(model, history, indices, patient_zeros, epi_params, states)

      !> If the epidemic spread, we exit the loop. Otherwise, we start over.
      if (count(states(t0, :) /= 'S') > 1) exit
    end do

    !> Only non-susceptible nodes may have started the epidemic.
    nodes = [(i, i=1,N)]
    non_S = nodes; call my_pack(non_S, states(t0, :) /= 'S')

    !> Each non-susceptible node has an energy that quantifies its likelihood of
    !> being the true patient zero.
    epi_size = size(non_S)
    allocate(energies(epi_size))

    do idx = 1, epi_size
      i = non_S(idx)

      !> We calculate the DMP marginal probabilities that each node j in a time-
      !> varying network is in states S, E, I or R at time t0, given that node i
      !> is the patient zero.
      if (.not. restricted) then
        call dmp_alg(model, history, indices, N, nodes, [i], epi_params, p)
      else
        call dmp_alg(model, history, indices, N, non_S, [i], epi_params, p)
      end if

      !> The probability that node i is the patient zero given the epidemic O is
      !> proportional to P(O|i), which we obtain by multiplying the DMP marginal
      !> probabilities (mean-field type approximation).
      total = 0
      joint = 1.0d0
      do j = 1, N
        select case (states(t0, j))
          case ('S')
            joint = joint * p(1, t0, j)
          case ('E')
            joint = joint * p(2, t0, j)
          case ('I')
            joint = joint * p(3, t0, j)
          case ('R')
            joint = joint * p(4, t0, j)
        end select

        !> 'joint' may become so small that GFortran treats it as zero. To avoid
        !> this, we multiply 'joint' by a huge constant each time it surpasses a
        !> certain threshold.
        if (abs(joint) < small) then
          joint = joint / small
          total = total + 1
        end if
      end do

      !> We compute the energy of each non-susceptible node. When joint <= 0, we
      !> give that particular node a large constant energy to avoid dealing with
      !> infinities.
      if (joint > 0.0d0) then
        energies(idx) = total*log(1/small) - log(joint)
      else
        energies(idx) = 1 / small
      end if
    end do

    deallocate(nodes)
    deallocate(p)

    !> We rank non-susceptible nodes by increasing value of their energy.
    call quicksort(energies, non_S)
  end subroutine pz_problem
end module patient_zero_problem
