    module FileUtils
    use AmlUtils
    implicit none
    !Utils using F2008 features

    contains

    subroutine FileUtils_Error(aname, msg, errormsg)
    character(LEN=*), intent(IN) :: aname, msg
    character(LEN=*), intent(IN), optional :: errormsg

    if (present(errormsg)) then
        call MpiStop(trim(errormsg)//' : '//trim(aname))
    else
        call MpiStop(trim(msg)//' : '//trim(aname))
    end if

    end subroutine 

    function OpenNewTxtFile(aname, errormsg) result(aunit)
    character(LEN=*), intent(IN) :: aname
    character(LEN=*), intent(IN), optional :: errormsg

    integer :: aunit

    aunit = OpenNewFile(aname,'formatted',errormsg)

    end function OpenNewTxtFile

    function OpenNewFile(aname, mode,errormsg) result(aunit)
    character(LEN=*), intent(IN) :: aname,mode
    integer status,aunit
    character(LEN=*), intent(IN), optional :: errormsg

    open(file=aname,form=mode,status='old', action='read', newunit=aunit, iostat=status)
    if (status/=0) call FileUtils_Error(aname, 'File not found', errormsg)

    end function OpenNewFile

    function CreateNewTxtFile(aname,errormsg) result(aunit)
    character(LEN=*), intent(IN) :: aname
    integer :: aunit
    character(LEN=*), intent(IN), optional :: errormsg

    aunit = CreateNewFile(aname,'formatted', errormsg)

    end function CreateNewTxtFile


    function CreateNewFile(aname, mode,errormsg) result(aunit)
    character(LEN=*), intent(IN) :: aname,mode
    integer :: aunit, status
    character(LEN=*), intent(IN), optional :: errormsg

    open(file=aname,form=mode,status='replace', newunit=aunit, iostat=status)
    if (status/=0) call FileUtils_Error(aname, 'Error creating file', errormsg)

    end function CreateNewFile


    function ReadLine(aunit, InLine) result(OK)
    integer, intent(IN) :: aunit
    character(LEN=*) :: InLine
    logical :: OK
    integer status

    read (aunit,'(a)',iostat=status) InLine
    OK = status==0

    end function ReadLine

    end module FileUtils