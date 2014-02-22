    module FileUtils
    use MpiUtils
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

    function FileExists(aname)
    character(LEN=*), intent(IN) :: aname
    logical FileExists

    inquire(file=aname, exist = FileExists)

    end function FileExists

    function DirectoryExists(aname)
    character(LEN=*), intent(IN) :: aname
    logical DirectoryExists

    !Does not work in some cases
    inquire(file=CheckTrailingSlash(aname)//'/.', exist=DirectoryExists)

    end function DirectoryExists


    function FileSize(name) result(fsize)
    integer fsize
    character(LEN=*), intent(in)::name

    inquire(file=name, size=fsize)

    end function FileSize

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


    function ReadLine(aunit, InLine, trimmed) result(OK)
    integer, intent(IN) :: aunit
    character(LEN=:), allocatable, optional :: InLine
    logical, intent(in), optional :: trimmed
    integer, parameter :: line_buf_len= 1024*4
    character(LEN=line_buf_len) :: InS
    logical :: OK, set
    integer status, size

    OK = .false.
    set = .true.
    do
        read (aunit,'(a)',advance='NO',iostat=status, size=size) InS
        OK = .not. IS_IOSTAT_END(status)
        if (.not. OK) return
        if (present(InLine)) then
            if (set) then
                InLine = InS(1:size)
                set=.false.
            else
                InLine = InLine // InS(1:size)
            end if
        end if
        if (IS_IOSTAT_EOR(status)) exit
    end do
    if (present(trimmed) .and. present(InLine)) then
        if (trimmed) InLine = trim(adjustl(InLine))
    end if

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


    function CheckTrailingSlash(aname)
    character(LEN=*), intent(in) :: aname
    character(LEN=:), allocatable :: CheckTrailingSlash
    integer alen
    character, parameter :: win_slash = char(92)

    alen = len_trim(aname)
    if (aname(alen:alen) /= win_slash .and. aname(alen:alen) /= '/') then
        CheckTrailingSlash = trim(aname)//'/'
    else
        CheckTrailingSlash = aname(1:alen)
    end if

    end function CheckTrailingSlash


    subroutine DeleteFile(aname)
    character(LEN=*), intent(IN) :: aname
    integer file_id, status

    if (FileExists(aname)) then
        open(newunit = file_id, file = aname, iostat=status)
        if (status/=0) return
        close(unit = file_id, status = 'DELETE')
    end if

    end subroutine DeleteFile


    subroutine FlushFile(aunit)
    integer, intent(IN) :: aunit

    flush(aunit)

    end subroutine FlushFile

    end module FileUtils
