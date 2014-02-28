    module MiscUtils
    implicit none

    INTERFACE PresentDefault
    module procedure PresentDefault_S, PresentDefault_L, PresentDefault_I
    END INTERFACE PresentDefault

    contains

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
    logical, intent(in), target :: default
    logical, intent(in), target, optional :: S
    logical :: Sout

    if (present(S)) then
        SOut = S
    else
        SOut = default
    end if
    end function PresentDefault_L
    
    
    function PresentDefault_I(default, S) result(Sout)
    integer, intent(in), target :: default
    integer, intent(in), target, optional :: S
    integer :: Sout

    if (present(S)) then
        SOut = S
    else
        SOut = default
    end if
    end function PresentDefault_I

    end module MiscUtils