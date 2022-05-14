!Module containing dynamic message-passing algorithms (DMP) for different
!epidemiological models (SIR and SEIR) in time-varying networks.
!Author: Adrià Meca Montserrat.
!Last modified date: 14/05/22.
module dmp_algorithms
  use array_procedures, only : dbl_list, int_list_list, pop
  use network_generation, only : node

  implicit none

  private

  public dmp_sir

contains
  !DMP algorithm for the SIR model.
  subroutine dmp_sir(history, indices, origin, lambda, mu, t0, ps, pi, pr)
    implicit none

    double precision, dimension(:), intent(out) :: ps, pi, pr
    double precision, intent(in) :: lambda, mu
    integer, intent(in) :: origin, t0
    type(int_list_list), dimension(:), intent(in) :: indices
    type(node), dimension(:), intent(in) :: history

    double precision, dimension(size(history)) :: ps0
    double precision :: theta_ki
    integer :: altk, hsize, i, id_ik, id_ki, isize, k, N, t
    type(dbl_list), dimension(size(history)) :: ops, nps, phi, theta

    N = size(history)

    !Initial conditions for the marginal probabilities.
    ps0 = 1.0d0
    ps0(origin) = 0.0d0

    ps = ps0
    pi = 1.0d0 - ps0
    pr = 0.0d0

    !Initial conditions for the DMP messages.
    do i = 1, N
      hsize = size(history(i)%neighbors)

      allocate(ops(i)%array(hsize), nps(i)%array(hsize), phi(i)%array(hsize), &
        theta(i)%array(hsize))

      do id_ki = 1, hsize
        k = history(i)%neighbors(id_ki)

        ops(i)%array(id_ki) = ps(k)
        nps(i)%array(id_ki) = ps(k)
        phi(i)%array(id_ki) = pi(k)
        theta(i)%array(id_ki) = 1.0d0
      end do
    end do

    !DMP iteration for the SIR model.
    do t = 1, t0
      do i = 1, N
        hsize = size(history(i)%neighbors)
        isize = size(indices(i)%time(t)%array)
        if (isize > 0) then
          do altk = 1, isize
            id_ki = indices(i)%time(t)%array(altk)
            theta(i)%array(id_ki) = theta(i)%array(id_ki) - lambda*phi(i)%array(id_ki)
          end do
          ps(i) = ps0(i) * product(theta(i)%array)
          do id_ki = 1, hsize
            k = history(i)%neighbors(id_ki)
            id_ik = history(i)%opposites(id_ki)
            theta_ki = theta(i)%array(id_ki)
            if (theta_ki > 0.0d0) then
              nps(k)%array(id_ik) = ps(i) / theta_ki
            else
              nps(k)%array(id_ik) = ps0(i) * product(pop(theta(i)%array, id_ki))
            end if
          end do
        end if
      end do

      do i = 1, N
        isize = size(indices(i)%time(t)%array)
        do altk = 1, isize
          id_ki = indices(i)%time(t)%array(altk)
          phi(i)%array(id_ki) = (1.0d0-lambda) * phi(i)%array(id_ki)
        end do
        phi(i)%array = (1.0d0-mu)*phi(i)%array + ops(i)%array - nps(i)%array
        ops(i)%array = nps(i)%array
        pr(i) = pr(i) + mu*pi(i)
        pi(i) = 1.0d0 - ps(i) - pr(i)
      end do
    end do
  end subroutine dmp_sir
end module dmp_algorithms
