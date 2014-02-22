    module StringUtils
    implicit none

    INTERFACE CONCAT
    module procedure concat_s, concat_s_n
    END INTERFACE


    INTERFACE RealToStr
    module procedure SingleToStr, DoubleToStr
    END INTERFACE RealToStr

    contains

    function GetParamCount()
    integer GetParamCount

    GetParamCount = command_argument_count() 

    end function GetParamCount

    function GetParam(i)
    character(LEN=:), allocatable :: GetParam
    integer, intent(in) :: i
    character(LEN=:), allocatable :: tmp
    integer l

    if (GetParamCount() < i) then
        GetParam = ''
    else
        call get_command_argument(i,length=l)
        allocate(character(l)::tmp)
        call get_command_argument(i,value=tmp)
        GetParam = trim(tmp)
    end if

    end function GetParam

    function StringStarts(S, substring) result(OK)
    character(LEN=*), intent(in) :: S, substring
    logical OK

    OK = S(1:min(len(S),len_trim(substring)))==substring

    end function

    subroutine StringReplace(FindS, RepS, S)
    character(LEN=*), intent(in) :: FindS, RepS
    character(LEN=:), allocatable, intent(inout) :: S
    integer i

    i = index(S,FindS)
    if (i>0) then
        S = S(1:i-1)//trim(RepS)//S(i+len_trim(FindS):len_trim(S))
    end if

    end subroutine StringReplace

    function numcat(S, num)
    character(LEN=*) S
    character(LEN=1024) numstr
    character(LEN=:), allocatable :: numcat
    integer num

    write (numstr, *) num
    numcat = trim(S) // trim(adjustl(numstr))
    !OK, so can probably do with with a format statement too... 
    end function numcat

    function IntToStr(I, minlen)
    integer , intent(in) :: I
    character(LEN=:), allocatable :: IntToStr
    integer, intent(in), optional :: minlen
    integer n
    character (LEN=128) :: form, tmp

    if (present(minlen)) then
        n = minlen
        if (I<0) n=n+1
        form = concat('(I',n,'.',minlen,')')
        write (tmp,form) i
        IntToStr = trim(tmp)
    else
        write (tmp,*) i
        IntToStr = trim(adjustl(tmp))
    end if

    end function IntToStr

    function StrToInt(S)
    integer :: StrToInt
    character(LEN=*), intent(in) :: S

    read(S,*) StrToInt

    end function StrToInt

    function concat_s(S1,S2,S3,S4,S5,S6,S7,S8) result(outstr)
    character(LEN=*), intent(in) :: S1, S2
    character(LEN=*), intent(in) , optional :: S3, S4, S5, S6,S7,S8
    character(LEN = 1000) concat
    character(LEN=:), allocatable :: outstr

    concat = trim(S1) // S2
    if (present(S3)) then
        concat = trim(concat) // S3
        if (present(S4)) then
            concat = trim(concat) // S4
            if (present(S5)) then
                concat = trim(concat) // S5
                if (present(S6)) then
                    concat = trim(concat) // S6
                    if (present(S7)) then
                        concat = trim(concat) // S7
                        if (present(S8)) then
                            concat = trim(concat) // S8
                        end if
                    end if    
                end if
            end if
        end if
    end if
    outstr = trim(concat)

    end function concat_s

    function concat_s_n(SS1,N2,SS3,N4,SS5,N6,SS7,N8,SS9,N10,SS11) result(outstr)
    character(LEN=*), intent(in) :: SS1
    integer, intent(in) :: N2
    character(LEN=*), intent(in) , optional :: SS3, SS5, SS7, SS9,SS11
    integer, intent(in), optional ::N4,N6,N8, N10
    character(LEN = 1000) concat
    character(LEN=:), allocatable :: outstr

    concat = trim(SS1) //trim(IntToStr(N2))
    if (present(SS3)) then
        concat = trim(concat) // SS3
        if (present(N4)) then
            concat = trim(concat) // trim(IntToStr(N4))
            if (present(SS5)) then
                concat = trim(concat) // SS5
                if (present(N6)) then
                    concat = trim(concat) // trim(intToStr(N6))
                    if (present(SS7)) then
                        concat = trim(concat) // SS7
                        if (present(N8)) then
                            concat = trim(concat) // trim(intToStr(N8))
                            if (present(SS9)) then
                                concat = trim(concat) // SS9
                                if (present(N10)) then
                                    concat = trim(concat) // trim(intToStr(N10))
                                    if (present(SS11)) then
                                        concat = trim(concat) // SS11
                                    end if
                                end if
                            end if       
                        end if
                    end if
                end if
            end if
        end if
    end if
    outstr = trim(concat)

    end function concat_s_n


    function DoubleToStr(R, figs)
    double precision, intent(in) :: R
    integer, intent(in), optional :: figs
    character(LEN=:), allocatable :: DoubleToStr

    DoubleToStr = SingleToStr(real(R),figs)

    end function DoubleToStr

    function SingleToStr(R, figs)
    real, intent(in) :: R
    integer, intent(in), optional :: figs
    character(LEN=:), allocatable :: SingleToStr
    character(LEN=30) tmp

    if (abs(R)>=0.001 .or. R==0.) then
        write (tmp,'(f12.6)') R

        tmp = adjustl(tmp)
        if (present(figs)) then
            SingleToStr = tmp(1:figs)
        else
            SingleToStr = tmp(1:6)
        end if
    else
        if (present(figs)) then
            write (tmp,trim(numcat('(E',figs))//'.2)') R
        else
            write (tmp,'(G9.2)') R
        end if
        SingleToStr = trim(adjustl(tmp))
    end if

    end function SingleToStr

    subroutine WriteFormatInts(unit, formatst, i1,i2,i3,i4)
    integer, intent(in) :: unit
    character(LEN=*), intent(in) :: formatst
    integer, intent(in) :: i1
    integer, intent(in),optional :: i2,i3,i4
    character(LEN=:), allocatable :: S

    S = formatst
    call StringReplace('%u', IntToStr(i1), S)
    if (present(i2)) call StringReplace('%u', IntToStr(i2), S)
    if (present(i3)) call StringReplace('%u', IntToStr(i3), S)
    if (present(i4)) call StringReplace('%u', IntToStr(i4), S)

    write(unit,'(a)') trim(S)

    end subroutine WriteFormatInts


    end module StringUtils
