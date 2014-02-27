    module FileUtils
    use MpiUtils
    use MiscUtils
    use StringUtils
    use, intrinsic :: ISO_FORTRAN_ENV, only : INT64
    implicit none
    !Utils using F2008 features
    private

    integer, parameter :: File_size_int = KIND(INT64)
    character(LEN=*), parameter :: mode_formatted = 'formatted'
    character(LEN=*), parameter :: mode_binary = 'unformatted'
    character(LEN=*), parameter :: access_formatted = 'sequential'
    character(LEN=*), parameter :: access_binary = 'stream'

    Type :: TFileStream
        integer :: unit = 0
        character(LEN=:), allocatable :: FileName
    contains
    procedure :: Open => TFileStream_Open
    procedure :: OpenFile => TFileStream_OpenFile
    procedure :: Flush => TFileStream_Flush
    procedure :: CreateFile
    procedure :: CreateOpenFile
    procedure :: Position => TFileStream_Position
    procedure :: Size => TFileStream_Size
    procedure :: Close => TFileStream_Close
    procedure :: Rewind => TFileStream_Rewind
    procedure :: Error
    procedure :: Opened
    procedure :: WriteTrim
    procedure, private :: ReadStringFunc
    procedure, private :: ReadStringSub
    procedure, private :: ReadItemSub
    procedure, private :: ReadItemFunc
    procedure, private :: ReadItems
    procedure, private :: ReadArray
    procedure, private :: ReadArray2
    procedure, private :: ReadArray2Func
    procedure, private :: ReadArrayFunc
    procedure, private :: WriteItemSub
    procedure, private :: WriteArray
    procedure, private :: WriteArray2
    procedure, private :: WriteItems
    procedure, private  :: WriteSizedArray1
    procedure, private  :: WriteSizedArray2
    procedure, private  :: ReadSizedArray_R
    procedure, private  :: ReadSizedArray_D
    procedure, private  :: ReadSizedArray_I
    generic :: Write => WriteItems, WriteArray, WriteArray2
    generic :: Read => ReadItems, ReadArray, ReadArray2
    generic :: ReadItem => ReadItemFunc, ReadArrayFunc, ReadArray2Func
    generic :: ReadString => ReadStringSub
    generic :: ReadStringItem => ReadStringFunc
    generic :: WriteSizedArray => WriteSizedArray1,WriteSizedArray2
    generic :: ReadSizedArray => ReadSizedArray_R,ReadSizedArray_D,ReadSizedArray_I 
    final :: TFileStream_Free
    end type

    Type, extends(TFileStream) :: TBinaryFile
    end type


    Type, extends(TFileStream) :: TTextFile
        character(LEN=20) :: RealFormat = '(*(E17.7))'
        character(LEN=20) :: IntegerFormat = '(*(I10))'
        logical :: AdvanceDefault = .true.
    contains
    procedure :: Open => OpenTxtFile
    procedure :: CreateFile => CreateTxtFile
    procedure :: CreateOpenFile => CreateOpenTxtFile
    procedure :: ReadLine
    procedure :: ReadLineSkipEmptyAndComments
    procedure :: NewLine
    procedure :: SkipLines
    procedure :: Lines
    procedure :: Columns
    procedure :: WriteLeftAligned
    procedure :: WriteItemTxt
    procedure :: WriteArrayTxt
    procedure :: WriteInLineItems
    procedure :: WriteInLineTrim
    procedure :: WriteFormat
    procedure :: WriteTrim => WriteTrimTxt
    procedure, private :: ReadStringSub => ReadStringTxt
    procedure, private :: ReadItemTxt
    procedure, private :: ReadArrayTxt
    procedure, private :: WriteInLineItem
    procedure, private :: WriteInLineArray
    procedure, private :: DefaultAdvance
    procedure, private :: WriteItem => TTextFile_WriteItem
    procedure, private :: WriteArray => TTextFile_WriteArray
    procedure, private :: WriteArray2 => WriteArray2Txt
    procedure, private :: ReadItemSub => ReadItemTxt
    procedure, private :: ReadArray => ReadArrayTxt
    procedure, private :: ReadArray2 => ReadArray2Txt
    procedure, private :: WriteItems => WriteItemsTxt
    generic :: WriteInLine => WriteInLineItem, WriteInLineArray
    end type

    !Functions on filenames and text
    !File instance below acts as a namespace
    Type TFile
    contains
    procedure, nopass :: TxtNumberColumns
    procedure, nopass :: TxtColumns
    procedure, nopass :: TxtFileColumns
    procedure, nopass :: LastLine => LastFileLine
    procedure, nopass :: Size => FileSize
    procedure, nopass :: Exists => FileExists
    procedure, nopass :: ExtractName => ExtractFileName
    procedure, nopass :: ExtractPath => ExtractFilePath
    procedure, nopass :: ChangeExt => ChangeFileExt
    procedure, nopass :: CheckTrailingSlash
    procedure, nopass :: ExtractExt => ExtractFileExt
    procedure, nopass :: Delete => DeleteFile
    procedure, nopass :: ReadTextMatrix
    procedure, nopass :: ReadTextVector
    procedure, nopass :: WriteTextVector
    end type

    type(TFile), save :: File

    public TFileStream, TBinaryFile, TTextFile, File, File_size_int
    contains

    function FileExists(aname)
    character(LEN=*), intent(IN) :: aname
    logical FileExists

    inquire(file=aname, exist = FileExists)

    end function FileExists

    function FileSize(name) result(fsize)
    integer(file_size_int) fsize
    character(LEN=*), intent(in)::name

    inquire(file=name, size=fsize)

    end function FileSize

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


    function Opened(this)
    class(TFileStream) :: this
    logical Opened
    Opened = this%unit /=0
    end function Opened

    subroutine Error(this, msg, errormsg)
    class(TFileStream) :: this
    character(LEN=*), intent(IN) :: msg
    character(LEN=*), intent(IN), optional :: errormsg

    if (present(errormsg)) then
        call MpiStop(trim(errormsg)//' : '//this%FileName )
    else
        call MpiStop(trim(msg)//' : '// this%FileName )
    end if

    end subroutine

    subroutine TFileStream_Close(this)
    class(TFileStream) :: this
    if (this%unit/=0) close(this%unit)
    this%unit=0
    end subroutine TFileStream_Close

    subroutine TFileStream_Free(this)
    Type(TFileStream) :: this

    call this%Close()

    end subroutine TFileStream_Free

    subroutine TFileStream_Flush(this)
    class(TFileStream) :: this
    flush(this%unit)

    end subroutine TFileStream_Flush

    subroutine TFileStream_Rewind(this)
    class(TFileStream) :: this
    rewind(this%unit)
    end subroutine TFileStream_Rewind

    function TFileStream_Position(this)
    class(TFileStream) :: this
    integer(file_size_int) TFileStream_Position

    inquire(this%unit, pos=TFileStream_Position)

    end function TFileStream_Position


    function TFileStream_Size(this)
    class(TFileStream) :: this
    integer(file_size_int) TFileStream_Size

    if (this%unit /=0) then
        inquire(this%unit, size=TFileStream_Size)
    else if (allocated(this%FileName)) then
        TFileStream_Size = File%Size(this%FileName)
    else
        call this%Error('File not defined for size')
    end if

    end function TFileStream_Size

    subroutine TFileStream_Open(this, aname, errormsg, status)
    class(TFileStream) :: this
    character(LEN=*), intent(IN) :: aname
    character(LEN=*), intent(IN), optional :: errormsg
    integer, intent(out), optional :: status

    call this%OpenFile(aname,errormsg=errormsg, status=status)

    end subroutine TFileStream_Open

    subroutine TFileStream_OpenFile(this, aname, mode, errormsg, forwrite, append, status)
    class(TFileStream) :: this
    character(LEN=*), intent(IN) :: aname
    character(LEN=*), intent(IN), optional :: mode
    character(LEN=*), intent(IN), optional :: errormsg
    integer, intent(out), optional :: status
    logical, intent(in), optional :: forwrite, append
    character(LEN=:), allocatable :: amode, access, state, action, pos
    integer out_status

    call this%Close()

    this%FileName = trim(aname)

    amode = PresentDefault(mode_binary, mode)
    if (amode == mode_formatted) then
        access = access_formatted
    else
        access = access_binary
    end if
    if (PresentDefault(.false., forwrite)) then
        state = 'replace'
        action = 'readwrite'
    else
        state = 'old'
        action = 'read'
    end if
    if (PresentDefault(.false., append) .and. FileExists(aname)) then
        pos = 'append'
        state = 'old'
        action = 'readwrite'
    else
        pos='asis'
    end if

    open(file=aname,form=amode,status=state, action=action, newunit=this%unit, &
    & iostat=out_status, position =pos,  access=access)
    if (present(status)) then
        status=out_status
        if (out_status/=0) this%unit = 0
    else
        if (out_status/=0) then
            if (state == 'replace') then
                call this%Error('Error creating file', errormsg)
            else
                call this%Error('File not found', errormsg)
            end if
            this%unit = 0
        end if
    end if
    end subroutine TFileStream_OpenFile

    subroutine CreateFile(this,aname, errormsg)
    class(TFileStream) :: this
    character(LEN=*), intent(IN) :: aname
    character(LEN=*), intent(IN), optional :: errormsg

    call this%OpenFile(aname, errormsg=errormsg, forwrite =.true.)

    end subroutine CreateFile

    subroutine CreateOpenFile(this, aname, append, errormsg)
    class(TFileStream) :: this
    character(LEN=*), intent(IN) :: aname
    logical, optional, intent(in) :: append
    character(LEN=*), intent(IN), optional :: errormsg

    call this%OpenFile(aname, forwrite =.true., append=append, errormsg=errormsg)

    end subroutine CreateOpenFile


    subroutine ReadStringSub(this, S, OK) 
    class(TFileStream) :: this
    character(LEN=:), allocatable :: S
    logical, optional :: OK
    logical isOK
    integer i, status

    read(this%unit, iostat=status) i
    isOK = status==0
    if (isOK) then
        if (allocated(S)) deallocate(S)
        allocate(character(LEN=i)::S)
        read(this%unit, iostat=status) S
        isOK = status==0
    end if
    if (present(OK)) OK = isOK
    end subroutine ReadStringSub

    function ReadStringFunc(this) result(S)
    class(TFileStream) :: this
    character(LEN=:), allocatable :: S
    call this%ReadString(S)
    end function ReadStringFunc

    subroutine ReadItemSub(this, R, OK)
    class(TFileStream) :: this
    class(*), intent(out) :: R
    logical, optional :: OK
    logical res
    integer status

    select type(R)
    type is (real)
        read(this%unit, iostat=status) R
    type is (double precision)
        read(this%unit, iostat=status) R
    type is (integer)
        read(this%unit, iostat=status) R
    type is (logical)
        read(this%unit, iostat=status) R
    type is (character(LEN=*))
        read(this%unit, iostat=status) R
        class default
        call this%Error('Unknown type to read')
    end select

    res = status==0
    if (status/=0 .and. (.not. IS_IOSTAT_END(status) .or. .not. present(OK))) &
    & call this%Error('Error reading item')
    if (present(OK)) OK = res
    end subroutine ReadItemSub


    function ReadItemFunc(this, R) result(res)
    class(TFileStream) :: this
    class(*), intent(out) :: R
    logical :: res

    call this%ReadItemSub(R, res)
    end function ReadItemFunc

    subroutine ReadArray(this, R, n, OK)
    class(TFileStream) :: this
    class(*) :: R(1:)
    integer, intent(in), optional :: n
    logical, optional :: OK
    integer status

    select type(R)
    type is (real)
        read(this%unit, iostat=status) R(1:PresentDefault(size(R),n))
    type is (double precision)
        read(this%unit, iostat=status) R(1:PresentDefault(size(R),n))
    type is (integer)
        read(this%unit, iostat=status) R(1:PresentDefault(size(R),n))
    type is (logical)
        read(this%unit, iostat=status) R(1:PresentDefault(size(R),n))
        class default
        call this%Error('Unknown type to read')
    end select
    if (status/=0 .and. (.not. IS_IOSTAT_END(status) .or. .not. present(OK))) &
    & call this%Error('Error reading item')
    if (present(OK)) OK = status==0
    end subroutine ReadArray


    function ReadArrayFunc(this, R, n) result(res)
    class(TFileStream) :: this
    class(*) :: R(1:)
    integer, intent(in), optional :: n
    logical :: res

    call this%ReadArray(R,n,res)

    end function ReadArrayFunc


    function ReadArray2Func(this, R) result(res)
    class(TFileStream) :: this
    class(*) :: R(:,:)
    logical :: res
    call this%ReadArray2(R,res)
    end function ReadArray2Func


    subroutine ReadArray2(this, R, OK)
    class(TFileStream) :: this
    class(*) :: R(:,:)
    logical, optional :: OK
    integer status

    select type(R)
    type is (real)
        read(this%unit, iostat=status) R
    type is (double precision)
        read(this%unit, iostat=status) R
    type is (integer)
        read(this%unit, iostat=status) R
    type is (logical)
        read(this%unit, iostat=status) R
        class default
        call this%Error('Unknown type to read')
    end select
    if (status/=0 .and. (.not. IS_IOSTAT_END(status) .or. .not. present(OK))) &
    & call this%Error('Error reading item')
    if (present(OK)) OK = status==0
    end subroutine ReadArray2


    subroutine ReadSizedArray_R(this, R) 
    class(TFileStream) :: this
    Real, allocatable :: R(:)
    integer sz

    call this%Read(sz)
    if (allocated(R)) deallocate(R)
    allocate(R(sz))
    call this%ReadArray(R)

    end subroutine ReadSizedArray_R


    subroutine ReadSizedArray_D(this, R) 
    class(TFileStream) :: this
    double precision, allocatable :: R(:)
    integer sz

    call this%Read(sz)
    if (allocated(R)) deallocate(R)
    allocate(R(sz))
    call this%ReadArray(R)

    end subroutine ReadSizedArray_D

    subroutine ReadSizedArray_I(this, R) 
    class(TFileStream) :: this
    integer, allocatable :: R(:)
    integer sz

    call this%Read(sz)
    if (allocated(R)) deallocate(R)
    allocate(R(sz))
    call this%ReadArray(R)

    end subroutine ReadSizedArray_I

    subroutine ReadSizedArray2_R(this, R) 
    class(TFileStream) :: this
    real, allocatable :: R(:,:)
    integer sz1, sz2

    call this%Read(sz1)
    call this%Read(sz2)
    if (allocated(R)) deallocate(R)
    allocate(R(sz1,sz2))
    call this%ReadArray2(R)

    end subroutine ReadSizedArray2_R

    subroutine ReadSizedArray2_D(this, R) 
    class(TFileStream) :: this
    double precision, allocatable :: R(:,:)
    integer sz1, sz2

    call this%Read(sz1)
    call this%Read(sz2)
    if (allocated(R)) deallocate(R)
    allocate(R(sz1,sz2))
    call this%ReadArray2(R)

    end subroutine ReadSizedArray2_D

    subroutine ReadSizedArray2_I(this, R) 
    class(TFileStream) :: this
    integer, allocatable :: R(:,:)
    integer sz1, sz2

    call this%Read(sz1)
    call this%Read(sz2)
    if (allocated(R)) deallocate(R)
    allocate(R(sz1,sz2))
    call this%ReadArray2(R)

    end subroutine ReadSizedArray2_I

    subroutine WriteItemSub(this, R)
    class(TFileStream) :: this
    class(*), intent(in) :: R

    select type(R)
    type is (real)
        Write(this%unit) R
    type is (double precision)
        Write(this%unit) R
    type is (integer)
        Write(this%unit) R
    type is (logical)
        Write(this%unit) R
    type is (character(LEN=*))
        Write(this%unit) len(R)
        Write(this%unit) R
        class default
        call this%Error('Unknown type to Write')
    end select

    end subroutine WriteItemSub

    subroutine WriteTrim(this, S)
    class(TFileStream) :: this
    character(LEN=*), intent(in) :: S
    call this%WriteItemSub(trim(S))
    end subroutine WriteTrim

    subroutine WriteSizedArray1(this, R, n) 
    class(TFileStream) :: this
    class(*), intent(in) :: R(1:)
    integer, intent(in), optional :: n
    integer sz

    sz=PresentDefault(size(R),n)
    call this%Write(sz)
    call this%WriteArray(R,n)

    end subroutine WriteSizedArray1

    subroutine WriteSizedArray2(this, R) 
    class(TFileStream) :: this
    class(*), intent(in) :: R(:,:)

    call this%Write(size(R,dim=1))
    call this%Write(size(R,dim=2))
    call this%WriteArray2(R)

    end subroutine WriteSizedArray2

    subroutine WriteArray(this, R, n) 
    class(TFileStream) :: this
    class(*), intent(in) :: R(1:)
    integer, intent(in), optional :: n
    integer sz

    sz=PresentDefault(size(R),n)
    select type(R)
    type is (real)
        Write(this%unit) R(1:sz)
    type is (double precision)
        Write(this%unit) R(1:sz)
    type is (integer)
        Write(this%unit) R(1:sz)
    type is (logical)
        Write(this%unit) R(1:sz)
        class default
        call this%Error('Unknown type to Write')
    end select

    end subroutine WriteArray

    subroutine WriteArray2(this, R) 
    class(TFileStream) :: this
    class(*), intent(in) :: R(:,:)

    select type(R)
    type is (real)
        Write(this%unit) R
    type is (double precision)
        Write(this%unit) R
    type is (integer)
        Write(this%unit) R
    type is (logical)
        Write(this%unit) R
        class default
        call this%Error('Unknown type to Write')
    end select

    end subroutine WriteArray2

    subroutine WriteItems(this, S1, S2,S3,S4,S5,S6)
    class(TFileStream) :: this
    class(*), intent(in) :: S1
    class(*), intent(in), optional :: S2, S3,S4,S5,S6

    call this%WriteItemSub(S1)
    if (present(S2)) call this%WriteItemSub(S2)
    if (present(S3)) call this%WriteItemSub(S3)
    if (present(S4)) call this%WriteItemSub(S4)
    if (present(S5)) call this%WriteItemSub(S5)
    if (present(S6)) call this%WriteItemSub(S6)

    end subroutine WriteItems

    subroutine ReadItems(this, S1, S2,S3,S4,S5,S6,OK)
    class(TFileStream) :: this
    class(*) S1
    class(*), optional :: S2,S3,S4,S5,S6
    logical, optional :: OK

    call this%ReadItemSub(S1)
    if (present(S2)) call this%ReadItemSub(S2)
    if (present(S3)) call this%ReadItemSub(S3)
    if (present(S4)) call this%ReadItemSub(S4)
    if (present(S5)) call this%ReadItemSub(S5)
    if (present(S6)) call this%ReadItemSub(S6)

    end subroutine ReadItems

    !Text unformatted files

    subroutine OpenTxtFile(this, aname, errormsg, status)
    class(TTextFile) :: this
    character(LEN=*), intent(IN) :: aname
    character(LEN=*), intent(IN), optional :: errormsg
    integer, intent(out), optional :: status

    call this%TFileStream%OpenFile(aname,'formatted',errormsg, status=status)

    end subroutine OpenTxtFile

    subroutine CreateTxtFile(this, aname,errormsg)
    class(TTextFile) :: this
    character(LEN=*), intent(IN) :: aname
    character(LEN=*), intent(IN), optional :: errormsg

    call this%OpenFile(aname, mode='formatted', forwrite =.true., errormsg=errormsg)

    end subroutine CreateTxtFile

    subroutine CreateOpenTxtFile(this,aname, append, errormsg)
    class(TTextFile) :: this
    character(LEN=*), intent(IN) :: aname
    logical, optional, intent(in) :: append
    character(LEN=*), intent(IN), optional :: errormsg

    call this%OpenFile(aname, mode='formatted', forwrite =.true., append=append, errormsg=errormsg)

    end subroutine CreateOpenTxtFile

    function ReadLine(this, InLine, trimmed) result(OK)
    class(TTextFile) :: this
    character(LEN=:), allocatable, optional :: InLine
    logical, intent(in), optional :: trimmed
    integer, parameter :: line_buf_len= 1024*4
    character(LEN=line_buf_len) :: InS
    logical :: OK, set
    integer status, size

    OK = .false.
    set = .true.
    do
        read (this%unit,'(a)',advance='NO',iostat=status, size=size) InS
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

    function ReadLineSkipEmptyAndComments(this, InLine, comment) result(OK)
    class(TTextFile) :: this
    logical :: OK
    character(LEN=:), allocatable :: InLine
    character(LEN=:), allocatable, optional, intent(out) :: comment

    if (present(comment)) comment=''
    do
        OK = this%ReadLine(InLine, trimmed=.true.)
        if (.not. OK) return
        if (InLine=='') cycle
        if (InLine(1:1)/='#') then
            return
        else
            if (present(comment)) comment = trim(InLine(2:))
        end if
    end do

    end function ReadLineSkipEmptyAndComments

    function SkipLines(this, n) result(OK)
    class(TTextFile) :: this
    integer, intent(in) :: n
    logical OK
    integer ix

    do ix = 1, n
        if (.not. this%ReadLine()) then
            OK = .false.
            return
        end if
    end do
    OK = .true.

    end function SkipLines

    function Columns(this) result(n)
    class(TTextFile) :: this
    integer n
    character(LEN=:), allocatable :: InLine

    if (this%ReadLine(InLine)) then
        n = File%TxtNumberColumns(InLine)
    else
        n=0
    end if
    call this%Rewind()

    end function Columns

    function Lines(this) result(n)
    class(TTextFile) :: this
    integer n

    n=0
    do while (this%ReadLine())
        n = n+1
    end do
    call this%Rewind()

    end function Lines

    function DefaultAdvance(this,advance)
    class(TTextFile) :: this
    logical, intent(in), optional :: advance
    character(3) :: DefaultAdvance

    if (PresentDefault(this%AdvanceDefault,advance)) then
        DefaultAdvance='YES'
    else
        DefaultAdvance='NO'
    end if

    end function

    subroutine NewLine(this)
    class(TTextFile) :: this
    call this%Write('')
    end subroutine

    subroutine WriteLeftAligned(this, Form, str)
    class(TTextFile) :: this
    character(LEN=*) str, Form
    Character(LEN=max(len(str),128)) tmp

    tmp = str
    write(this%unit,form, advance='NO') tmp

    end subroutine WriteLeftAligned


    subroutine WriteInLineItems(this, S1, S2,S3,S4,S5,S6)
    class(TTextFile) :: this
    class(*), intent(in) :: S1
    class(*), intent(in), optional :: S2, S3,S4,S5,S6

    call this%WriteInLine(S1)
    if (present(S2)) call this%WriteInLine(S2)
    if (present(S3)) call this%WriteInLine(S3)
    if (present(S4)) call this%WriteInLine(S4)
    if (present(S5)) call this%WriteInLine(S5)
    if (present(S6)) call this%WriteInLine(S6)

    end subroutine WriteInLineItems

    subroutine WriteItemsTxt(this, S1, S2,S3,S4,S5,S6)
    class(TTextFile) :: this
    class(*), intent(in) :: S1
    class(*), intent(in), optional :: S2, S3,S4,S5,S6

    call this%WriteInLineItems(S1, S2,S3,S4,S5,S6)
    call this%NewLine()

    end subroutine WriteItemsTxt

    subroutine WriteInLineItem(this,str,form)
    class(TTextFile) :: this
    class(*), intent(in) :: str
    character(LEN=*), intent(in), optional :: form
    call this%WriteItemTxt(str,form,.false.)
    end subroutine WriteInLineItem

    subroutine WriteInLineArray(this,str,form,n)
    class(TTextFile) :: this
    class(*), intent(in) :: str(:)
    character(LEN=*), intent(in), optional :: form
    integer, intent(in), optional :: n

    call this%WriteArrayTxt(str,form,.false.,number=n)
    end subroutine WriteInLineArray

    subroutine TTextFile_WriteItem(this,R)
    class(TTextFile) :: this
    class(*), intent(in) :: R
    call this%WriteItemTxt(R,advance=.true.)
    end subroutine TTextFile_WriteItem

    subroutine TTextFile_WriteArray(this,R,n)
    class(TTextFile) :: this
    class(*), intent(in) :: R(:)
    integer, intent(in), optional :: n
    call this%WriteArrayTxt(R,number=n,advance=.true.)
    end subroutine TTextFile_WriteArray

    subroutine WriteItemTxt(this, str, form, advance)
    class(TTextFile) :: this
    class(*), intent(in) :: str
    character(LEN=*), intent(in), optional :: form
    logical, intent(in), optional :: advance
    character(LEN=3) :: Ad

    Ad = this%DefaultAdvance(advance)
    select type(str)
    type is (character(LEN=*))
        if (Ad=='YES') then
            write(this%unit,PresentDefault('(a)',form)) trim(str)
        else
            write(this%unit,PresentDefault('(a)',form), advance=Ad) str
        end if
    type is (real)
        write(this%unit,PresentDefault(this%RealFormat,form), advance=Ad) str
    type is (double precision)
        write(this%unit,PresentDefault(this%RealFormat,form), advance=Ad) str
    type is (integer)
        write(this%unit,PresentDefault(this%IntegerFormat,form), advance=Ad) str
    type is (logical)
        write(this%unit,PresentDefault('(L2)',form), advance=Ad) str
        class default
        call this%Error('unknown type to write')
    end select

    end subroutine WriteItemTxt

    subroutine WriteArrayTxt(this, str, form, advance, number)
    class(TTextFile) :: this
    class(*), intent(in) :: str(:)
    character(LEN=*), intent(in), optional :: form
    logical, intent(in), optional :: advance
    integer, intent(in), optional :: number
    integer n
    character(LEN=3) :: Ad

    Ad = this%DefaultAdvance(advance)
    n = PresentDefault(size(str),number)
    select type(str)
    type is (character(LEN=*))
        write(this%unit,PresentDefault('(a)',form), advance=Ad) str(1:n)
    type is (real)
        write(this%unit,PresentDefault(this%RealFormat,form), advance=Ad) str(1:n)
    type is (double precision)
        write(this%unit,PresentDefault(this%RealFormat,form), advance=Ad) str(1:n)
    type is (integer)
        write(this%unit,PresentDefault(this%IntegerFormat,form), advance=Ad) str(1:n)
    type is (logical)
        write(this%unit,PresentDefault('(*(L2))',form), advance=Ad) str(1:n)
        class default
        call this%Error('unknown type to write')
    end select

    end subroutine WriteArrayTxt

    subroutine WriteArray2Txt(this, R)
    class(TTextFile) :: this
    class(*), intent(in) :: R(:,:)
    integer i

    do i=1, size(R,1)
        call this%WriteArrayTxt(R(i,:))
    end do

    end subroutine WriteArray2Txt


    subroutine WriteTrimTxt(this, S)
    class(TTextFile) :: this
    character(LEN=*), intent(in) :: S

    call this%WriteItemTxt(S, advance=.true.)

    end subroutine WriteTrimTxt

    subroutine WriteInLineTrim(this, string)
    class(TTextFile) :: this
    character(LEN=*), intent(in) :: string

    call this%WriteItemTxt(trim(string), advance=.false.)
    end subroutine WriteInLineTrim


    subroutine ReadArray2Txt(this, R, OK)
    class(TTextFile) :: this
    class(*) ::  R(:,:)
    logical, optional :: OK
    integer i

    do i=1, size(R,1)
        call this%ReadArrayTxt(R(i,:), OK=OK)
        if (present(OK)) then
            if (.not. OK) return
        end if
    end do

    end subroutine ReadArray2Txt


    subroutine ReadItemTxt(this, R, OK)
    class(TTextFile) :: this
    class(*), intent(out) :: R
    logical, optional :: OK

    integer status
    select type(R)
    type is (character(LEN=*))
        Read(this%unit,'(a)', iostat=status) R
    type is (real)
        Read(this%unit,*, iostat=status) R
    type is (double precision)
        Read(this%unit,*, iostat=status) R
    type is (integer)
        Read(this%unit,*, iostat=status) R
        class default
        call this%Error('unknown type to Read')
    end select
    if (status/=0 .and. (.not. IS_IOSTAT_END(status) .or. .not. present(OK))) &
    & call this%Error('Error reading item')
    if (present(OK)) OK = status==0

    end subroutine ReadItemTxt

    subroutine ReadArrayTxt(this, R, n, OK)
    class(TTextFile) :: this
    class(*) :: R(1:)
    integer, intent(in), optional :: n
    logical, optional :: OK
    integer status

    select type(R)
    type is (real)
        Read(this%unit,*, iostat=status) R(1:PresentDefault(size(R),n))
    type is (double precision)
        Read(this%unit,*, iostat=status) R(1:PresentDefault(size(R),n))
    type is (integer)
        Read(this%unit,*, iostat=status) R(1:PresentDefault(size(R),n))
        class default
        call this%Error('unknown type to Read')
    end select

    if (status/=0 .and. (.not. IS_IOSTAT_END(status) .or. .not. present(OK))) &
    & call this%Error('Error reading item')
    if (present(OK)) OK = status==0

    end subroutine ReadArrayTxt

    subroutine ReadStringTxt(this, S, OK)
    class(TTextFile) :: this
    character(LEN=:), allocatable :: S
    logical, optional :: OK
    logical isOK

    isOK = this%ReadLine(S)
    if (present(OK)) OK = isOK

    end subroutine  ReadStringTxt

    subroutine WriteFormat(this, formatst, i1,i2,i3,i4,i5,i6)
    class(TTextFile) :: this
    character(LEN=*), intent(in) :: formatst
    class(*), intent(in) :: i1
    class(*), intent(in),optional :: i2,i3,i4,i5,i6

    write(this%unit,'(a)') FormatString(formatst,i1,i2,i3,i4,i5,i6)

    end subroutine WriteFormat

    !Misc functions

    function TopCommentLine(aname) result(res)
    character(LEN=*), intent(IN) :: aname
    character(LEN=:), allocatable :: res
    Type(TTextFile) :: F

    call F%Open(aname)
    res=''
    do while (res == '' .and. F%ReadLine(res))
    end do
    If (res(1:1)/='#') then
        res = ''
    end if
    call F%Close()

    end function TopCommentLine


    function TxtFileColumns(aname) result(n)
    character(LEN=*), intent(IN) :: aname
    integer n
    Type(TTextFile) :: F

    call F%Open(aname)
    n = F%Columns()
    call F%Close()

    end function TxtFileColumns


    function LastFileLine(aname)
    character(LEN=*), intent(IN) :: aname
    character(LEN=:), allocatable :: LastFileLine
    Type(TTextFile) :: F

    LastFileLine = ''
    call F%Open(aname)
    do while (F%ReadLine(LastFileLine))
    end do
    call F%Close()
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


    subroutine ReadTextVector(aname, vec, n)
    character(LEN=*), intent(IN) :: aname
    integer, intent(in) :: n
    class(*), intent(out) :: vec(n)
    integer j
    Type(TTextFile) :: F

    call F%Open(aname)
    do j=1,n
        select type(vec)
        type is (real)
            if (.not. F%ReadItem(vec(j))) call F%Error('vector file is the wrong size')
        type is (double precision)
            if (.not. F%ReadItem(vec(j))) call F%Error('vector file is the wrong size')
            class default
            stop 'wrong type for vector'
        end select
    end do
    call F%Close()

    end subroutine ReadTextVector

    subroutine WriteTextVector(aname, vec, n)
    character(LEN=*), intent(IN) :: aname
    integer, intent(in) :: n
    class(*), intent(in) :: vec(n)
    integer j
    Type(TTextFile) :: F

    call F%CreateFile(aname)
    do j=1,n
        select type(vec)
        type is (real)
            call F%Write(vec(j))
        type is (double precision)
            call F%Write(vec(j))
            class default
            stop 'wrong type for vector'
        end select
    end do
    call F%Close()

    end subroutine WriteTextVector

    subroutine ReadTextMatrix(aname, mat, inm,inn)
    character(LEN=*), intent(IN) :: aname
    integer, intent(in), optional :: inm,inn
    real(kind(1.d0)), intent(out) :: mat(:,:)
    integer j,k, status, n,m
    real tmp
    Type(TTextFile) :: F

    m = PresentDefault(size(mat,dim=1),inm)
    n = PresentDefault(size(mat,dim=2),inn)
    call F%Open(aname)

    do j=1,m
        read (F%unit,*, iostat=status) mat(j,1:n)
        if (status/=0) exit
    end do
    if (status/=0) then
        call F%Rewind()  !Try other possible format
        do j=1,m
            do k=1,n
                read (F%unit,*, iostat=status) mat(j,k)
                if (status/=0) call F%Error( 'matrix file is the wrong size')
            end do
        end do
    end if

    read (F%unit,*, iostat=status) tmp
    if (status==0) call F%Error( 'matrix file is the wrong size (too big)')

    call F%Close()

    end subroutine ReadTextMatrix

    end module FileUtils
