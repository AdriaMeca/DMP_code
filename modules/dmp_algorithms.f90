!> Procedures for computing the marginal probabilities that each node in a time-
!> varying network is in states S, E, I or R at times t <= t0.
!> Author: Adria Meca Montserrat.
!> Last modified date: 24/08/22.
module dmp_algorithms
  use array_procedures, only: my_pack, pop
  use derived_types,    only: dbl_list, int_llist, node, prm

  implicit none

  private

  public :: dmp_alg

contains
  !> Uses the DMP algorithm to compute the marginal probabilities that each node
  !> in a time-varying network is in states S, E, I or R at times t <= t0.
  subroutine dmp_alg( &
      model,          &
      history,        &
      indices,        &
      N,              &
      nodes,          &
      patient_zeros,  &
      epi_params,     &
      p               &
  )
    integer,          intent(in)  :: N                                  !> Number of nodes.
    integer,          intent(in)  :: nodes(:)                           !> Nodes participating in the iteration.
    integer,          intent(in)  :: patient_zeros(:)                   !>
    integer                       :: alti, altk, i, ik, k, ki, t, t0    !>
    integer                       :: hsize, isize, nsize                !>
    character(len=*), intent(in)  :: model                              !> Epidemiological model.
    double precision, intent(out) :: p(:, 0:, :)                        !> DMP marginal probabilities.
    double precision              :: pa, pl, pm, pn, xi                 !>
    type(dbl_list)                :: mf(N), mn(N), mo(N), mp(N), mt(N)  !> DMP messages.
    type(int_llist),  intent(in)  :: indices(:)                         !> Active links throughout the simulation.
    type(node),       intent(in)  :: history(:)                         !> Rewiring history.
    type(prm),        intent(in)  :: epi_params                         !> Epidemiological parameters.

    !> Local observation time.
    t0 = epi_params%t0

    !> We choose the local epidemiological parameters based on the model in use.
    select case (trim(model))
      case ('SIR')
        pa = epi_params%lambda
        pl = 0.0d0
        pm = 0.0d0
        pn = epi_params%mu
        xi = 0.0d0
      case default
        pa = epi_params%alpha
        pl = epi_params%lambda
        pm = epi_params%mu
        pn = epi_params%nu
        xi = epi_params%nu
    end select

    !> Initial condition for the marginal probability of being S.
    p(1, 0, :) = 1.0d0; p(1, 0, patient_zeros) = 0.0d0

    !> Initial conditions for the marginal probabilities of being E, I or R, respectively.
    p(2, 0, :) = 1.0d0 - p(1, 0, :)
    p(3, 0, :) = 0.0d0
    p(4, 0, :) = 0.0d0

    !> Initial conditions for the DMP messages.
    do i = 1, N
      hsize = size(history(i)%neighbors)

      allocate(mf(i)%array(hsize))
      allocate(mn(i)%array(hsize))
      allocate(mo(i)%array(hsize))
      allocate(mp(i)%array(hsize))
      allocate(mt(i)%array(hsize))

      do ki = 1, hsize
        k = history(i)%neighbors(ki)

        mf(i)%array(ki) = p(3, 0, k)
        mn(i)%array(ki) = p(1, 0, k)
        mo(i)%array(ki) = p(1, 0, k)
        mp(i)%array(ki) = p(2, 0, k)
        mt(i)%array(ki) = 1.0d0
      end do
    end do

    !> We compute the DMP marginal probabilities at times t <= t0.
    nsize = size(nodes)
    do t = 1, t0
      do alti = 1, nsize
        i = nodes(alti)

        !> Number of active links of node i at time t.
        isize = size(indices(i)%time(t)%array)

        if (isize > 0) then
          !> We update the thetas that enter node i at time t.
          do altk = 1, isize
            ki = indices(i)%time(t)%array(altk)

            mt(i)%array(ki) = mt(i)%array(ki) - pa*mp(i)%array(ki) - pl*mf(i)%array(ki)
          end do

          !> We calculate the probability that node i is S at time t.
          p(1, t, i) = p(1, 0, i) * product(mt(i)%array)

          !> We update the npses that leave node i at time t.
          do ki = 1, size(history(i)%neighbors)
            k = history(i)%neighbors(ki); ik = history(i)%opposites(ki)

            if (mt(i)%array(ki) > 0.0d0) then
              mn(k)%array(ik) = p(1, t, i) / mt(i)%array(ki)
            else
              mn(k)%array(ik) = p(1, 0, i) * product(pop(mt(i)%array, ki))
            end if
          end do
        else
          p(1, t, i) = p(1, t-1, i)
        end if
      end do

      do alti = 1, nsize
        i = nodes(alti)

        !> We update the phis and psis that enter node i at time t.
        do altk = 1, size(indices(i)%time(t)%array)
          ki = indices(i)%time(t)%array(altk)

          mf(i)%array(ki) = (1.0d0-pl) * mf(i)%array(ki)
          mp(i)%array(ki) = (1.0d0-pa) * mp(i)%array(ki)
        end do
        mf(i)%array = (1.0d0-pm) * mf(i)%array + mp(i)%array * xi
        mp(i)%array = (1.0d0-pn) * mp(i)%array + mo(i)%array - mn(i)%array

        !> We update the opses that enter node i at time t.
        mo(i)%array = mn(i)%array

        !> We compute the probabilities that node i is R or I at time t, respectively.
        p(4, t, i) = p(4, t-1, i) + pm*p(3, t-1, i)
        p(3, t, i) = p(3, t-1, i) - pm*p(3, t-1, i) + pn*p(2, t-1, i)

        !> We compute the probability that node i is E at time t.
        p(2, t, i) = 1.0d0 - p(1, t, i) - p(3, t, i) - p(4, t, i)
      end do
    end do

    !> For the SIR model, we must perform the following exchanges:
    if (trim(model) == 'SIR') p(3:4, :, :) = p(2:3, :, :); p(2, :, :) = 0.0d0
  end subroutine dmp_alg
end module dmp_algorithms
