    module MiscUtils
    implicit none

    INTERFACE PresentDefault
    module procedure PresentDefault_S, PresentDefault_L, PresentDefault_I, PresentDefault_R, PresentDefault_D
    END INTERFACE PresentDefault

    INTERFACE IfThenElse
    module procedure IfThenElse_S, IfThenElse_L, IfThenElse_I, IfThenElse_R, IfThenElse_D
    END INTERFACE IfThenElse

    integer, parameter :: TTimer_dp = Kind(1.d0)
    Type TTimer
        real(TTimer_dp) start_time
    contains
    procedure Start => TTimer_Start
    procedure Time => TTimer_Time
    procedure WriteTime => TTimer_WriteTime    
    end type TTimer

    contains

    function TimerTime()
    real(TTimer_dp) time
    real(TTimer_dp) :: TimerTime
#ifdef MPI
    TimerTime = MPI_WTime()
#else
    call cpu_time(time)
    TimerTime=  time
#endif
    end function TimerTime

    subroutine TTimer_Start(this)
    class(TTimer) :: this
    this%start_time = TimerTime()
    end subroutine TTimer_Start

    real(TTimer_dp) function TTimer_Time(this)
    class(TTimer) :: this
    TTimer_Time =  TimerTime() - this%start_time
    end function TTimer_Time

    subroutine TTimer_WriteTime(this,Msg, start)
    class(TTimer) :: this
    character(LEN=*), intent(in), optional :: Msg
    real(TTimer_dp), optional :: start
    real(TTimer_dp) T

    if (present(start)) then
        T=start
    else
        T=this%start_time
    end if

    if (present(Msg)) then
        write (*,*) trim(Msg)//': ', TimerTime() - T
    end if
    if (.not. present(start)) this%start_time = TimerTime()

    end subroutine TTimer_WriteTime


    function PresentDefault_S(default, S) result(Sout)
    character(LEN=*), intent(in), target :: default
    character(LEN=*), intent(in), target, optional :: S
    character(LEN=:), pointer :: Sout

    if (present(S)) then
        SOut => S
    else
        SOut => default
    end if
    end function PresentDefault_S

    function PresentDefault_L(default, S) result(Sout)
    logical, intent(in) :: default
    logical, intent(in), optional :: S
    logical :: Sout

    if (present(S)) then
        SOut = S
    else
        SOut = default
    end if
    end function PresentDefault_L

    function PresentDefault_I(default, S) result(Sout)
    integer, intent(in) :: default
    integer, intent(in), optional :: S
    integer :: Sout

    if (present(S)) then
        SOut = S
    else
        SOut = default
    end if
    end function PresentDefault_I

    function PresentDefault_R(default, S) result(Sout)
    real, intent(in) :: default
    real, intent(in), optional :: S
    real:: Sout

    if (present(S)) then
        SOut = S
    else
        SOut = default
    end if
    end function PresentDefault_R

    function PresentDefault_D(default, S) result(Sout)
    double precision, intent(in) :: default
    double precision, intent(in), optional :: S
    double precision:: Sout

    if (present(S)) then
        SOut = S
    else
        SOut = default
    end if
    end function PresentDefault_D


    function IfThenElse_S(flag, either, or) result(IfThenElse)
    logical, intent(in) :: flag
    character(LEN=:), pointer :: IfThenElse
    character(LEN=*), target :: either, or

    if (flag) then
        IfThenElse => either
    else
        IfThenElse => or
    end if

    end function

    function IfThenElse_L(flag, either, or) result(IfThenElse)
    logical, intent(in) :: flag
    logical :: IfThenElse, either, or

    if (flag) then
        IfThenElse = either
    else
        IfThenElse = or
    end if

    end function

    function IfThenElse_I(flag, either, or) result(IfThenElse)
    logical, intent(in) :: flag
    integer :: IfThenElse, either, or

    if (flag) then
        IfThenElse = either
    else
        IfThenElse = or
    end if

    end function

    function IfThenElse_R(flag, either, or) result(IfThenElse)
    logical, intent(in) :: flag
    real :: IfThenElse, either, or

    if (flag) then
        IfThenElse = either
    else
        IfThenElse = or
    end if

    end function
    function IfThenElse_D(flag, either, or) result(IfThenElse)
    logical, intent(in) :: flag
    double precision :: IfThenElse, either, or

    if (flag) then
        IfThenElse = either
    else
        IfThenElse = or
    end if

    end function


    end module MiscUtils
