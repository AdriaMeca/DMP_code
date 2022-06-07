!Module whose procedures generate random numbers.
!Authors: Matteo Palassini and Adria Meca Montserrat.
!Last modified date: 07/06/22.
module random_number_generator
  implicit none

  private

  !Global variables.
  double precision, parameter :: inv_maxint=1.0d0/2147483647.0d0

  integer, dimension(0:2047) :: index1, index2, irand
  integer :: ioffset

  public ir1279, r1279, setr1279

contains
  !Function that returns random integers from the uniform distribution over the
  !interval [0, maxint].
  function ir1279()
    implicit none

    !Output arguments.
    integer :: ir1279


    ioffset = iand(ioffset+1, 2047)
    irand(ioffset) = irand(index1(ioffset)) * irand(index2(ioffset))
    ir1279 = ishft(irand(ioffset), -1)
  end function ir1279



  !Function that returns random doubles from the uniform distribution over the
  !interval [0.0, 1.0].
  function r1279()
    implicit none

    !Output arguments.
    double precision :: r1279


    r1279 = ir1279() * inv_maxint
  end function r1279



  !L'Ecuyer's long-period random number generator (> 2x10^{18}) with Bays-Durham
  !shuffle and additional safeguards. Returns a uniform random deviation between
  !0.0 and 1.0 (excluding endpoint values). Call with idum a negative integer to
  !initialize; thereafter, do not alter idum between successive deviations in a
  !sequence. The variable rnmx should approximate the largest float value that
  !is less than 1. This function has been taken from the book: Numerical Recipes
  !in Fortran 77, and modified so that gfortran stops complaining so much.
  function ran2(idum)
    implicit none

    !Input arguments.
    integer, intent(inout) :: idum

    !Output arguments.
    double precision :: ran2

    !Local variables.
    double precision, parameter :: am=1.0d0/2147483563.0d0, eps=1.2d-7, &
      rnmx=1-eps

    integer, parameter :: im1=2147483563, im2=2147483399, imm1=im1-1, &
      ia1=40014, ia2=40692, iq1=53668, iq2=52774, ir1=12211, ir2=3791, &
      ntab=32, ndiv=67108862
    integer, save :: idum2=123456789, iv(ntab)=ntab*0, iy=0
    integer :: idx1, idx2, j, k


    !Initialize.
    if (idum <= 0) then
      !Be sure to prevent idum = 0.
      idum = max(-idum, 1)
      idum2 = idum

      !Load the shuffle table (alternative version).
      idx1 = 1; idx2 = ntab
      do while (idx2 >= 1)
        k = idum / iq1
        idum = ia1*(idum-k*iq1) - k*ir1
        if (idum < 0) idum = idum + im1

        !Eight warm-ups.
        if (idx1 <= 8) then
          idx1 = idx1 + 1
          cycle
        end if

        iv(idx2) = idum
        idx2 = idx2 - 1
      end do
      iy = iv(1)
    end if

    !Start here when not initializing.
    k = idum / iq1
    !Compute idum = mod(ia1*idum, im1) without overflows by Schrage's method.
    idum = ia1*(idum-k*iq1) - k*ir1
    if (idum < 0) idum = idum + im1

    k = idum2 / iq2
    !Compute idum2 = mod(ia2*idum2, im2) likewise.
    idum2 = ia2*(idum2-k*iq2) - k*ir2
    if (idum2 < 0) idum2 = idum2 + im2

    !Will be in the range 1:ntab.
    j = 1 + iy/ndiv
    !Here idum is shuffled, idum and idum2 are combined to generate output.
    iy = iv(j) - idum2
    iv(j) = idum
    if (iy < 1) iy = iy + imm1

    !Because users do not expect endpoint values.
    ran2 = min(am*iy, rnmx)
  end function ran2



  !Subroutine that initializes the Lagged Fibonacci random number generator.
  subroutine setr1279(iseed)
    implicit none

    !Input arguments.
    integer, intent(in) :: iseed

    !Local variables.
    integer, parameter :: nbitm1=31
    integer :: ibit, ispoke, one_bit, localseed


    !Initialize ioffset. This will be increased by (1 mod 2048) for
    !each random number which is called.
    ioffset = 0

    !Set up the two arrays which give the locations of the two random
    !numbers which will be multiplied to get the new random number.
    do ispoke = 0, 2047
      index1(ispoke) = iand(ispoke-1279, 2047)
      index2(ispoke) = iand(ispoke-418, 2047)
    end do

    !Set up the initial array of 2048 integer random numbers.
    !Each bit is separately initialized using the function ran2.
    localseed = -abs(iseed)

    !Matteo: I have substituted lshift with ishft.
    do ispoke = 0, 2047
      irand(ispoke) = 0
      do ibit = 0, nbitm1
        one_bit = 0
        if (ran2(localseed) > 0.5d0) one_bit = 1
        irand(ispoke) = ior(irand(ispoke), ishft(one_bit, ibit))
      end do
      irand(ispoke) = 2*irand(ispoke) + 1
    end do
  end subroutine setr1279
end module random_number_generator
