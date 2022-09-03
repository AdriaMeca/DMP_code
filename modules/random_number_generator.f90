!> Procedures for generating random numbers.
!> Authors: Matteo Palassini and Adria Meca Montserrat.
!> Last modified date: 03/09/22.
module random_number_generator
  implicit none

  private

  !> Global variables.
  integer                     :: index1(0:2047), index2(0:2047), irand(0:2047)  !>
  integer                     :: ioffset                                        !>
  double precision, parameter :: inv_maxint = 1.0d0 / 2147483647.0d0            !>

  public :: ir1279, r1279, setr1279

contains
  !> Returns random doubles uniformly distributed over [0.0d0, 1.0d0].
  function r1279()
    integer          :: rngint  !>
    double precision :: r1279   !> Output.

    call ir1279(rngint)

    r1279 = rngint * inv_maxint
  end function r1279


  !> Returns a uniform random deviation between 0.0d0 and 1.0d0 (excluding endpoints).
  !> Adapted from the book 'Numerical recipes in Fortran 77'.
  function ran2(idum)
    integer,                     intent(inout) :: idum                                     !>
    integer,          parameter                :: ia1 = 40014, ia2 = 40692                 !>
    integer,          parameter                :: im1 = 2147483563, imm1 = im1 - 1         !>
    integer,          parameter                :: im2 = 2147483399                         !>
    integer,          parameter                :: iq1 = 53668, iq2 = 52774                 !>
    integer,          parameter                :: ir1 = 12211, ir2 = 3791                  !>
    integer,          parameter                :: ntab = 32, ndiv = 67108862               !>
    integer,          save                     :: idum2 = 123456789, iv(ntab) = 0, iy = 0  !>
    integer                                    :: idx1, idx2, j, k                         !>
    double precision, parameter                :: am = 1.0d0 / 2147483563.0d0              !>
    double precision, parameter                :: eps = 1.2d-7, rnmx = 1.0d0 - eps         !>
    double precision                           :: ran2                                     !> Output.

    !> Initialization.
    if (idum <= 0) then
      !> We prevent idum = 0.
      idum = max(-idum, 1)
      idum2 = idum

      !> We load the shuffle table (alternative version).
      idx1 = 1
      idx2 = ntab
      do while (idx2 >= 1)
        k = idum / iq1
        idum = ia1*(idum-k*iq1) - k*ir1
        if (idum < 0) idum = idum + im1

        !> Eight warm-ups.
        if (idx1 <= 8) then
          idx1 = idx1 + 1
          cycle
        end if

        iv(idx2) = idum
        idx2 = idx2 - 1
      end do
      iy = iv(1)
    end if

    !> We compute idum = mod(ia1*idum, im1) without overflows by Schrage's method.
    k = idum / iq1
    idum = ia1*(idum-k*iq1) - k*ir1
    if (idum < 0) idum = idum + im1

    !> We compute idum2 = mod(ia2*idum2, im2) likewise.
    k = idum2 / iq2
    idum2 = ia2*(idum2-k*iq2) - k*ir2
    if (idum2 < 0) idum2 = idum2 + im2

    !> j lies within the interval [1, ntab].
    j = 1 + iy/ndiv

    !> We shuffle 'idum', which is combined with 'idum2' to generate output.
    iy = iv(j) - idum2; iv(j) = idum
    if (iy < 1) iy = iy + imm1

    !> This is done because users do not expect endpoint values.
    ran2 = min(am*iy, rnmx)
  end function ran2


  !> Returns random integers uniformly distributed over [0, maxint].
  subroutine ir1279(rngint)
    integer, intent(out) :: rngint  !> Random integer that lies within the interval [0, maxint].

    ioffset = iand(ioffset+1, 2047)
    irand(ioffset) = irand(index1(ioffset)) * irand(index2(ioffset))

    rngint = ishft(irand(ioffset), -1)
  end subroutine ir1279


  !> Initializes the Lagged Fibonacci random number generator.
  subroutine setr1279(iseed)
    integer,            intent(in) :: iseed                             !>
    integer, parameter             :: nbitm1 = 31                       !>
    integer                        :: ibit, ispoke, one_bit, localseed  !>

    !> We initialize 'ioffset', which will increase by (1 mod 2048) for each call.
    ioffset = 0

    !> We set up the arrays that give the locations of the two random numbers which
    !> will be multiplied to produce a new random number.
    do ispoke = 0, 2047
      index1(ispoke) = iand(ispoke-1279, 2047)
      index2(ispoke) = iand(ispoke-418, 2047)
    end do

    !> We set up the initial array of 2048 random integers; each bit is separately
    !> initialized using the function 'ran2'.
    localseed = -abs(iseed)
    do ispoke = 0, 2047
      irand(ispoke) = 0
      do ibit = 0, nbitm1
        one_bit = 0
        if (ran2(localseed) > 0.5d0) one_bit = 1

        !> Matteo: I have substituted 'lshift' with 'ishft'.
        irand(ispoke) = ior(irand(ispoke), ishft(one_bit, ibit))
      end do
      irand(ispoke) = 2*irand(ispoke) + 1
    end do
  end subroutine setr1279
end module random_number_generator
