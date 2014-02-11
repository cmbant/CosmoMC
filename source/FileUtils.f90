    module FileUtils
    implicit none
    !Utils using F2008 features

    integer, parameter :: max_line_length = 1024*64

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


    function FileExists(aname)
    character(LEN=*), intent(IN) :: aname
    logical FileExists

    inquire(file=aname, exist = FileExists)

    end function FileExists


    function OpenNewTxtFile(aname, errormsg) result(aunit)
    character(LEN=*), intent(IN) :: aname
    character(LEN=*), intent(IN), optional :: errormsg

    integer :: aunit

    aunit = OpenNewFile(aname,'formatted',errormsg)

    end function OpenNewTxtFile

    function OpenNewFile(aname, mode,errormsg) result(aunit)
    character(LEN=*), intent(IN) :: aname
    character(LEN=*), intent(IN), optional :: mode
    integer status,aunit
    character(LEN=*), intent(IN), optional :: errormsg
    character(LEN=:), allocatable :: amode

    if (present(mode)) then
        amode=mode
    else
        amode = 'unformatted'
    end if
    open(file=aname,form=amode,status='old', action='read', newunit=aunit, iostat=status)
    if (status/=0) call FileUtils_Error(aname, 'File not found', errormsg)

    end function OpenNewFile

    function CreateNewTxtFile(aname,errormsg) result(aunit)
    character(LEN=*), intent(IN) :: aname
    integer :: aunit
    character(LEN=*), intent(IN), optional :: errormsg

    aunit = CreateNewFile(aname,'formatted', errormsg)

    end function CreateNewTxtFile


    function CreateNewFile(aname, mode,errormsg) result(aunit)
    character(LEN=*), intent(IN) :: aname
    character(LEN=*), intent(in), optional :: mode
    character(LEN=*), intent(IN), optional :: errormsg
    integer :: aunit, status
    character(LEN=:), allocatable :: amode

    if (present(mode)) then
        amode=mode
    else
        amode = 'unformatted'
    end if
    open(file=aname,form=amode,status='replace', newunit=aunit, iostat=status)
    if (status/=0) call FileUtils_Error(aname, 'Error creating file', errormsg)

    end function CreateNewFile


    function  CreateOpenNewTxtFile(aname, append) result(aunit)
    character(LEN=*), intent(IN) :: aname
    integer :: aunit
    logical, optional, intent(in) :: append
    logical A

    if (present(append)) then
        A=append
    else
        A = .false.
    endif

    aunit = CreateOpenNewFile(aname,'formatted',A)

    end function CreateOpenNewTxtFile

    function CreateOpenNewFile(aname, mode, append) result(aunit)
    character(LEN=*), intent(IN) :: aname
    character(LEN=*), intent(IN), optional :: mode
    integer :: aunit
    logical, optional, intent(in) :: append
    logical A
    integer status
    character(LEN=:), allocatable :: amode

    if (present(append)) then
        A=append
    else
        A = .false.
    endif
    if (present(mode)) then
        amode=mode
    else
        amode='unformatted'
    end if

    if (A) then
        open(newunit=aunit,file=aname,form=amode,status='unknown', iostat=status, position='append')
    else
        open(newunit=aunit,file=aname,form=amode,status='replace', iostat=status)
    end if
    if (status/=0) call MpiStop('Error creatinging or opening '//trim(aname))

    end function CreateOpenNewFile


    function ReadLine(aunit, InLine) result(OK)
    integer, intent(IN) :: aunit
    character(LEN=:), allocatable, optional :: InLine
    character(LEN=max_line_length) :: InS
    logical :: OK
    integer status

    read (aunit,'(a)',iostat=status) InS
    OK = status==0
    if (OK .and. present(InLine)) InLine = trim(InS)

    end function ReadLine


    function TxtNumberColumns(InLine) result(n)
    character(LEN=*) :: InLine
    integer n,i
    logical isNum    

    n=0
    isNum=.false.
    do i=1, len_trim(InLIne)
        if (verify(InLine(i:i),'-+eE.0123456789') == 0) then
            if (.not. IsNum) n=n+1
            IsNum=.true.
        else
            IsNum=.false.     
        end if
    end do

    end function TxtNumberColumns

    function TxtColumns(InLine) result(n)
    character(LEN=*) :: InLine
    integer n,i
    logical isNum

    n=0
    isNum=.false.
    do i=1, len_trim(InLine)
        if (InLine(i:i) > char(32)) then
            if (.not. IsNum) n=n+1
            IsNum=.true.
        else
            IsNum=.false.
        end if
    end do

    end function TxtColumns

    function FileColumns(aunit) result(n)
    integer, intent(in) :: aunit
    integer n
    character(LEN=:), allocatable :: InLine

    if (ReadLine(aunit, InLine)) then
        n = TxtNumberColumns(InLine)
    else
        n=0
    end if
    rewind aunit

    end function FileColumns

    function FileLines(aunit) result(n)
    integer, intent(in) :: aunit
    integer n

    n=0
    do while (ReadLine(aunit))
        n = n+1
    end do
    rewind aunit

    end function FileLines


    function TopCommentLine(aname) result(res)
    character(LEN=*), intent(IN) :: aname
    integer file_id 
    character(LEN=:), allocatable :: res

    file_id = OpenNewTxtFile(aname)
    res=''
    do while (res == '' .and. ReadLine(file_id,res))
    end do
    If (res(1:1)/='#') then
        res = ''
    end if
    close(file_id)

    end function TopCommentLine


    function TxtFileColumns(aname) result(n)
    character(LEN=*), intent(IN) :: aname
    integer n, file_id 

    file_id = OpenNewTxtFile(aname)
    n = FileColumns(file_id)
    close(file_id)

    end function TxtFileColumns


    function LastFileLine(aname)
    character(LEN=*), intent(IN) :: aname
    character(LEN=:), allocatable :: LastFileLine
    integer file_id

    LastFileLine = ''
    file_id= OpenNewTxtFile(aname)
    do while (ReadLine(file_id, LastFileLine))
    end do
    close(file_id)

    end function LastFileLine


    function ExtractFilePath(aname)
    character(LEN=*), intent(IN) :: aname
    character(LEN=:), allocatable :: ExtractFilePath
    integer len, i

    len = len_trim(aname)
    do i = len, 1, -1
        if (aname(i:i)=='/') then
            ExtractFilePath = aname(1:i)
            return
        end if
    end do
    ExtractFilePath = ''

    end function ExtractFilePath

    function ExtractFileExt(aname)
    character(LEN=*), intent(IN) :: aname
    character(LEN=:), allocatable :: ExtractFileExt
    integer len, i

    len = len_trim(aname)
    do i = len, 1, -1
        if (aname(i:i)=='/') then
            ExtractFileExt = ''
            return
        else if (aname(i:i)=='.') then
            ExtractFileExt= aname(i:len)   
            return
        end if
    end do
    ExtractFileExt = ''

    end function ExtractFileExt


    function ExtractFileName(aname)
    character(LEN=*), intent(IN) :: aname
    character(LEN=:), allocatable :: ExtractFileName
    integer len, i

    len = len_trim(aname)
    do i = len, 1, -1
        if (aname(i:i)=='/') then
            ExtractFileName = aname(i+1:len)
            return
        end if
    end do
    ExtractFileName = trim(aname)

    end function ExtractFileName

    function ChangeFileExt(aname,ext)
    character(LEN=*), intent(IN) :: aname,ext
    character(LEN=:), allocatable :: ChangeFileExt
    integer len, i

    len = len_trim(aname)
    do i = len, 1, -1
        if (aname(i:i)=='.') then
            ChangeFileExt = aname(1:i) // trim(ext)
            return
        end if
    end do
    ChangeFileExt = trim(aname) // '.' // trim(ext)

    end function ChangeFileExt

    end module FileUtils
