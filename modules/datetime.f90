!Module whose procedures measure the time it takes to execute a piece of code.
module datetime
  implicit none

  private

  public total_time

contains
  !Function that calculates the time lapse between t1 and t2, writing it in the
  !most suitable units of time.
  function total_time(t1, t2) result(message)
    implicit none

    real, intent(in) :: t1, t2

    integer, parameter :: case1=60, case2=case1*60, case3=case2*24
    character(len=19) :: message
    character(len=11) :: header='# CPU_TIME='
    real :: time

    time = t2 - t1
    select case (floor(time))
      case (:case1-1)
        write(message, 100) header, time, 's.'
      case (case1:case2-1)
        write(message, 100) header, time/case1, 'm.'
      case (case2:case3-1)
        write(message, 100) header, time/case2, 'h.'
      case (case3:)
        write(message, 100) header, time/case3, 'd.'
    end select
    100 format (a11,f5.2,a3)
  end function total_time
end module datetime
