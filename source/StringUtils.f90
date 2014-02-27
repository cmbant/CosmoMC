    module StringUtils
    implicit none


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
    character(LEN=:), allocatable :: numcat
    integer num

    numcat = concat(S,num)
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

    subroutine StringAppend(S,X)
    character(LEN=:), allocatable :: S
    class(*) X

    if (.not. allocated(S)) S=''
    select type (X)
    type is (character(LEN=*))
        S = S // trim(X)
    type is (integer)
        S = S // IntToStr(X)
    type is (real)
        S = S // RealToStr(X)
    type is (double precision)
        S=S //RealToStr(X)
        class default
        stop 'StringAppend: Unknown type'
    end select
    end subroutine

    function concat(S1,S2,S3,S4,S5,S6,S7,S8) result(outstr)
    character(LEN=*), intent(in) :: S1
    class(*), intent(in) :: S2
    class(*), intent(in) , optional :: S3, S4, S5, S6,S7,S8
    character(LEN=:), allocatable :: outstr

    outstr=S1
    call StringAppend(outstr,S2)
    if (present(S3)) then
        call StringAppend(outstr,S3)
        if (present(S4)) then
            call StringAppend(outstr,S4)
            if (present(S5)) then
                call StringAppend(outstr,S5)
                if (present(S6)) then
                    call StringAppend(outstr,S6)
                    if (present(S7)) then
                        call StringAppend(outstr,S7)
                        if (present(S8)) then
                            call StringAppend(outstr,S8)
                        end if
                    end if    
                end if
            end if
        end if
    end if

    end function concat


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

    function SubNextFormat(S, X) result(OK)
    character(LEN=:), allocatable :: S
    class(*) X
    logical OK
    integer ix, P
    character c

    P=1 
    do
        ix=scan(S(P:),'%')
        OK = ix/=0 .and. ix < len(S)
        if (.not. OK) return
        c = S(ix+P:ix+P)
        if (c=='%') then
            P=P+Ix+1
        else
            exit
        end if
    end do

    select type (X)
    type is (integer)
        if (c/='u') stop 'Wrong format for type'
        call StringReplace('%u', IntToStr(X), S)
    type is (Character(LEN=*))
        if (c/='s') stop 'Wrong format for type'
        call StringReplace('%s', X, S)
        class default
        stop 'Unsupported format type'
    end select

    end function SubNextFormat

    function FormatString(formatst, i1,i2,i3,i4,i5,i6) result(S)
    character(LEN=*), intent(in) :: formatst
    class(*), intent(in) :: i1
    class(*), intent(in),optional :: i2,i3,i4,i5,i6
    character(LEN=:), allocatable :: S
    logical OK
    !Note that this routine is incomplete and very simple (so buggy in complex cases)
    !(should not substitute on the previously substituted string, etc, etc..)

    S = formatst
    OK = SubNextFormat(S, i1)
    if (OK .and. present(i2)) OK = SubNextFormat(S, i2)
    if (OK .and. present(i3)) OK = SubNextFormat(S, i3)
    if (OK .and. present(i4)) OK = SubNextFormat(S, i4)
    if (OK .and. present(i5)) OK = SubNextFormat(S, i5)
    if (OK .and. present(i6)) OK = SubNextFormat(S, i6)
    if (.not. OK) stop 'FormatString: Wrong number or kind of formats in string'
    call StringReplace('%%', '%', S)

    end function FormatString


    end module StringUtils
