    !Module to read in name/value pairs from a file, with each line of the form line 'name = value'
    !Should correctly interpret FITS headers
    !Antony Lewis (http://cosmologist.info/). Released to the public domain.
    !2014 using Fortran 2003/2008

    module IniObjects
    use FileUtils
    use StringUtils
    use MiscUtils
    implicit none
    private

    character(LEN=0), target :: Empty_String = ''
    integer, parameter :: Ini_Enumeration_Len = 64

    type TNameValue
        character(LEN=:), allocatable :: Name
        character(LEN=:), allocatable :: Value
    end type TNameValue

    type TNameValue_pointer
        Type(TNameValue), pointer :: P
    end type TNameValue_pointer

    Type TNameValueList
        integer :: Count = 0
        integer :: Delta = 128
        integer :: Capacity = 0
        logical :: ignoreDuplicates = .false.
        logical :: AllowDuplicateKeys = .false.
        !Use pointer so arrays can be moved around quickly without deep copies
        type(TNameValue_pointer), dimension(:), allocatable :: Items
    contains
    procedure :: Init => TNameValueList_Init
    procedure :: Clear => TNameValueList_Clear
    procedure :: ValueOf => TNameValueList_ValueOf
    procedure :: IndexOf => TNameValueList_IndexOf
    procedure :: HasKey => TNameValueList_HasKey
    procedure :: SetCapacity => TNameValueList_SetCapacity
    procedure :: Delete => TNameValueList_Delete
    procedure :: Error => TNameValueList_Error
    procedure :: FailStop => TNameValueList_FailStop
    procedure :: AddString => TNameValueList_Add
    procedure :: Name => TNameValueList_Name
    procedure :: Value => TNameValueList_Value
    procedure :: Override => TNameValueList_Override
    procedure :: TNameValueList_AddDouble
    procedure :: TNameValueList_AddReal
    procedure :: TNameValueList_AddInt
    procedure :: TNameValueList_AddLogical
    generic :: Add => AddString, TNameValueList_AddDouble, TNameValueList_AddReal, TNameValueList_AddInt, TNameValueList_AddLogical
    final :: TNameValueList_Free
    end Type TNameValueList

    Type, extends(TNameValueList) :: TIniFile
        logical :: SlashComments = .false.
        logical :: Echo_Read = .false.
        logical :: Fail_on_not_found = .false.
        character(LEN=:), allocatable :: Original_filename
        Type(TNameValueList) :: ReadValues
    contains
    procedure :: Open => Ini_Open
    procedure :: Open_Fromlines => Ini_Open_Fromlines
    procedure :: Close => Ini_Close
    procedure :: Read_String => Ini_Read_String
    procedure :: Read_String_Default => Ini_Read_String_Default
    procedure :: Read_String_Array => Ini_Read_String_Array
    procedure :: Read_Int_Array => Ini_Read_Int_Array
    procedure :: Read_Int => Ini_Read_Int
    procedure :: Read_Double => Ini_Read_Double
    procedure :: Read_Double_Array => Ini_Read_Double_Array
    procedure :: Read_Real => Ini_Read_Real
    procedure :: Read_Real_Array => Ini_Read_Real_Array
    procedure :: Read_Logical => Ini_Read_Logical
    procedure :: Read_Enumeration => Ini_Read_Enumeration
    procedure :: Read_Enumeration_List => Ini_Read_Enumeration_List
    procedure :: TestEqual => Ini_TestEqual
    procedure :: SaveReadValues => Ini_SaveReadValues
    procedure :: Key_To_Arraykey => Ini_Key_To_Arraykey
    procedure :: EnumerationValue => Ini_EnumerationValue
    procedure, private :: NameValue_AddLine => Ini_NameValue_AddLine
    procedure, private :: EmptyCheckDefault => Ini_EmptyCheckDefault
    procedure, private :: Ini_Read_Real_Change
    procedure, private :: Ini_Read_Double_Change
    procedure, private :: Ini_Read_Int_Change
    procedure, private :: Ini_Read_Logical_Change
    procedure, private :: Ini_Read_String_Change
    generic :: Read => Ini_Read_Real_Change, Ini_Read_Double_Change, Ini_Read_Int_Change, &
        & Ini_Read_Logical_Change,Ini_Read_String_Change
    end Type TIniFile


    public TNameValueList, TIniFile, Ini_Enumeration_Len
    contains

    subroutine TNameValueList_Init(L, ignoreDuplicates)
    class(TNameValueList) :: L
    logical, intent(in), optional :: ignoreDuplicates

    L%Count = 0
    L%Capacity = 0
    L%Delta = 128
    L%ignoreDuplicates = DefaultFalse(ignoreDuplicates)

    end subroutine TNameValueList_Init

    subroutine TNameValueList_Clear(L)
    class(TNameValueList) :: L
    integer i, status

    do i=L%count,1,-1
        deallocate (L%Items(i)%P, stat = status)
    end do
    deallocate (L%Items, stat = status)
    call L%Init()

    end subroutine TNameValueList_Clear

    subroutine TNameValueList_Free(L)
    Type(TNameValueList) :: L
    call L%Clear()
    end subroutine TNameValueList_Free

    function TNameValueList_ValueOf(L, AName) result(AValue)
    class(TNameValueList), intent(in) :: L
    character(LEN=*), intent(in) :: AName
    character(LEN=:), pointer :: AValue
    integer i

    i = L%IndexOf(AName)
    if (i/=-1) then
        AValue => L%Items(i)%P%Value
    else
        AValue => Empty_String
    end if

    end function TNameValueList_ValueOf

    function TNameValueList_Name(L, i) result(Name)
    class(TNameValueList), intent(in) :: L
    integer i
    character(LEN=:), pointer :: Name

    Name => L%Items(i)%P%Name

    end function TNameValueList_Name

    function TNameValueList_Value(L, i) result(Value)
    class(TNameValueList), intent(in) :: L
    integer i
    character(LEN=:), pointer :: Value

    Value => L%Items(i)%P%Value

    end function TNameValueList_Value


    function TNameValueList_IndexOf(L, AName) result (AValue)
    class(TNameValueList), intent(in) :: L
    character(LEN=*), intent(in) :: AName
    integer :: AValue
    integer i

    do i=1, L%Count
        if (L%Items(i)%P%Name == AName) then
            AValue = i
            return
        end if
    end do
    AValue = -1

    end function TNameValueList_IndexOf

    function TNameValueList_HasKey(L, AName) result (AValue)
    class(TNameValueList), intent(in) :: L
    character(LEN=*), intent(in) :: AName
    logical :: AValue

    AValue = L%IndexOf(AName) /= -1

    end function TNameValueList_HasKey
    
    subroutine TNameValueList_Override(L, Settings, only_if_exists)
    class(TNameValueList) :: L
    class(TNameValueList), intent(in) :: Settings
    logical, intent(in), optional :: only_if_exists
    integer i, ix
    character(LEN=:), pointer :: name
    
    do i=1, Settings%Count
        name => Settings%Name(i)
        ix = L%IndexOf(name)
        if (ix/=-1) then
            L%Items(i)%P%Value = Settings%Value(i)
        elseif (.not. DefaultFalse(only_if_exists)) then
            call L%Add(name, Settings%Value(i))
        end if
    end do

    end subroutine TNameValueList_Override

    subroutine TNameValueList_Add(L, AName, AValue, only_if_undefined)
    class(TNameValueList) :: L
    character(LEN=*), intent(in) :: AName, AValue
    logical, optional, intent(in) :: only_if_undefined
    logical isDefault

    isDefault = DefaultTrue(only_if_undefined)

    if ((.not. L%AllowDuplicateKeys .or. isDefault) .and. L%HasKey(AName)) then
        if (L%ignoreDuplicates .or. isDefault) return
        call L%Error('duplicate key name',AName)
    end if
    if (L%Count == L%Capacity) call L%SetCapacity(L%Capacity + L%Delta)
    L%Count = L%Count + 1
    allocate(L%Items(L%Count)%P)
    L%Items(L%Count)%P%Name = trim(adjustl(AName))
    L%Items(L%Count)%P%Value = trim(adjustl(AValue))

    end subroutine TNameValueList_Add

    subroutine TNameValueList_AddReal(L, AName, AValue)
    class(TNameValueList) :: L
    character(LEN=*), intent(in) :: AName
    real, intent(in) :: AValue
    character(LEN=32) tmp

    write(tmp,*) AValue
    call L%AddString(AName, Tmp)

    end subroutine TNameValueList_AddReal

    subroutine TNameValueList_AddDouble(L, AName, AValue)
    class(TNameValueList) :: L
    character(LEN=*), intent(in) :: AName
    double precision, intent(in) :: AValue
    character(LEN=32) tmp

    write(tmp,*) AValue
    call L%AddString(AName, Tmp)

    end subroutine TNameValueList_AddDouble


    subroutine TNameValueList_AddInt(L, AName, AValue)
    class(TNameValueList) :: L
    character(LEN=*), intent(in) :: AName
    integer, intent(in) :: AValue
    character(LEN=32) tmp

    write(tmp,*) AValue
    call L%AddString(AName, Tmp)

    end subroutine TNameValueList_AddInt


    subroutine TNameValueList_AddLogical(L, AName, AValue)
    class(TNameValueList) :: L
    character(LEN=*), intent(in) :: AName
    logical, intent(in) :: AValue
    character(LEN=32) tmp

    write(tmp,*) AValue
    call L%AddString(AName, Tmp)

    end subroutine TNameValueList_AddLogical


    subroutine TNameValueList_SetCapacity(L, C)
    class(TNameValueList) :: L
    integer C
    type(TNameValue_pointer), dimension(:), allocatable :: TmpItems

    if (L%Count > 0) then
        if (C < L%Count) call L%Error('TNameValueList_SetCapacity, smaller than Count')
        allocate(TmpItems(C))
        TmpItems(:L%Count) = L%Items(:L%Count)
        call move_alloc(TmpItems, L%Items)
    else
        allocate(L%Items(C))
    end if
    L%Capacity = C

    end subroutine TNameValueList_SetCapacity

    subroutine TNameValueList_Delete(L, i)
    class(TNameValueList) :: L
    integer, intent(in) :: i

    deallocate(L%Items(i)%P)
    if (L%Count > 1) L%Items(i:L%Count-1) = L%Items(i+1:L%Count)
    L%Count = L%Count -1

    end subroutine TNameValueList_Delete

    subroutine TNameValueList_FailStop(L)
    class(TNameValueList) :: L

    stop
    end subroutine TNameValueList_FailStop

    subroutine TNameValueList_Error(L, Msg, Key)
    class(TNameValueList) :: L
    character(LEN=*) :: Msg
    character(LEN=*), optional :: Key

    if (present(Key)) then
        write(*,*) 'Error for key "'//trim(Key)//'" : '//Msg
    else
        write(*,*) 'Error :'//Msg
    end if
    call L%FailStop()

    end subroutine TNameValueList_Error

    subroutine Ini_EmptyCheckDefault(Ini, Key, Default)
    class(TIniFile) :: Ini
    character(LEN=*), intent(in) :: Key
    class(*), optional, intent(in) :: Default

    if (.not. present(Default)) call Ini%Error('missing key',Key)

    end subroutine Ini_EmptyCheckDefault



    subroutine Ini_NameValue_AddLine(Ini,AInLine,only_if_undefined)
    class(TIniFile) :: Ini
    character (LEN=*), intent(IN) :: AInLine
    integer EqPos, slashpos, lastpos
    logical, optional, intent(in) :: only_if_undefined
    logical isDefault
    character (LEN=len(AInLine)) :: AName, S, InLine

    isDefault = DefaultFalse(only_if_undefined)

    InLine=trim(adjustl(AInLine))
    EqPos = scan(InLine,'=')
    if (EqPos/=0 .and. InLine(1:1)/='#' .and. .not. StringStarts(InLine,'COMMENT')) then
        AName = trim(InLine(1:EqPos-1))

        S = adjustl(InLine(EqPos+1:))
        if (Ini%SlashComments) then
            slashpos=scan(S,'/')
            if (slashpos /= 0) then
                S  = S(1:slashpos-1)
            end if
        end if
        lastpos=len_trim(S)
        if (lastpos>1) then
            if (S(1:1)=='''' .and. S(lastpos:lastpos)=='''') then
                S = S(2:lastpos-1)
            end if
        end if
        call Ini%Add(AName, S,only_if_undefined = isDefault)
    end if

    end subroutine Ini_NameValue_AddLine

    recursive subroutine Ini_Open(Ini, filename, error, slash_comments, append,only_if_undefined)
    class(TIniFile) :: Ini
    character (LEN=*), intent(IN) :: filename
    logical, intent(OUT), optional :: error
    logical, optional, intent(IN) :: slash_comments
    logical, optional, intent(in) :: append, only_if_undefined
    character (LEN=:), allocatable :: IncludeFile
    character(LEN=:), allocatable :: InLine
    integer lastpos, i, status
    Type(TNameValueList) IncudeFiles, DefaultValueFiles
    logical FileExists, isDefault
    Type(TTextFile) :: F


    isDefault = DefaultFalse(only_if_undefined)

    if (.not. DefaultFalse(append)) then
        call Ini%TNameValueList%Init()
        call Ini%ReadValues%Init(.true.)
        Ini%Original_filename = filename
    end if

    Ini%SlashComments = DefaultFalse(slash_comments)

    call F%Open(filename, status=status)
    if (status/=0) then
        if (present(error)) then
            error=.true.
            return
        else
            call Ini%Error('Ini_Open, file not found: '//trim(filename))
        end if
    end if

    call IncudeFiles%Init()
    call DefaultValueFiles%Init()

    do
        if (.not. F%ReadLineSkipEmptyAndComments(InLine)) exit
        if (InLine == 'END') exit
        if (StringStarts(InLine,'INCLUDE(')) then
            lastpos = scan(InLine,')')
            if (lastpos/=0) then
                call IncudeFiles%Add(InLine(9:lastpos-1),'')
            else
                call Ini%Error('Ini_Open, error in INCLUDE line: '//trim(filename))
            end if
        elseif (StringStarts(InLine,'DEFAULT(')) then
            !Settings to read in as defaults, overridden by subsequent re-definitions
            lastpos = scan(InLine,')')
            if (lastpos/=0) then
                call DefaultValueFiles%Add(InLine(9:lastpos-1),'')
            else
                call Ini%Error('Ini_Open, error in DEFAULT line: '//trim(filename))
            end if
        elseif (InLine /= '') then
            call Ini%NameValue_AddLine(InLine, only_if_undefined=isDefault)
        end if
    end do

    call F%Close()
    if (present(error)) error=.false.

    do i=1, IncudeFiles%Count
        if (DefaultFalse(error)) exit
        IncludeFile= IncudeFiles%Items(i)%P%Name
        inquire(file=IncludeFile, exist = FileExists)
        if (.not. FileExists) then
            IncludeFile= File%ExtractPath(filename)//trim(IncludeFile)
            inquire(file=IncludeFile, exist = FileExists)
            if (.not. FileExists) then
                call Ini%Error('Ini_Open, error in INCLUDE file not found: '//trim(IncudeFiles%Items(i)%P%Name))
            end if
        end if
        call Ini%Open(IncludeFile, error, slash_comments, append=.true.,only_if_undefined=isDefault)
    end do
    do i=1, DefaultValueFiles%Count
        if (DefaultFalse(error)) exit
        IncludeFile=DefaultValueFiles%Items(i)%P%Name
        inquire(file=IncludeFile, exist = FileExists)
        if (.not. FileExists) then
            IncludeFile=trim(File%ExtractPath(filename))//trim(IncludeFile)
            inquire(file=IncludeFile, exist = FileExists)
            if (.not. FileExists) then
                call Ini%Error('Ini_Open, error in DEFAULT file not found: '//trim(DefaultValueFiles%Items(i)%P%Name))
            end if
        end if
        call Ini%Open(IncludeFile, error, slash_comments, append=.true., only_if_undefined=.true.)
    end do
    call IncudeFiles%Clear()
    call DefaultValueFiles%Clear()

    end subroutine Ini_Open

    subroutine Ini_Open_Fromlines(Ini, Lines, NumLines, slash_comments)
    class(TIniFile) :: Ini
    integer, intent(IN) :: NumLines
    character (LEN=*), dimension(NumLines), intent(IN) :: Lines
    logical, intent(IN) :: slash_comments
    integer i

    call Ini%TNameValueList%Init()
    call Ini%ReadValues%Init(.true.)

    Ini%SlashComments = slash_comments

    do i=1,NumLines
        call Ini%NameValue_AddLine(Lines(i))
    end do

    end  subroutine Ini_Open_Fromlines


    subroutine Ini_Close(Ini)
    class(TIniFile) :: Ini

    call Ini%TNameValueList%Clear()
    call Ini%ReadValues%Clear()

    end  subroutine Ini_Close

    function Ini_Read_String(Ini, Key, NotFoundFail) result(AValue)
    class(TIniFile) :: Ini
    character (LEN=*), intent(IN) :: Key
    logical, optional, intent(IN) :: NotFoundFail
    character(LEN=:), pointer :: AValue

    AValue => Ini%ValueOf(Key)

    if (AValue/='') then
        call  Ini%ReadValues%Add(Key, AValue)
        if (Ini%Echo_Read) write (*,*) trim(Key)//' = ',trim(AValue)
        return
    end if
    if (PresentDefault(Ini%fail_on_not_found, NotFoundFail)) then
        call Ini%Error('key not found',Key)
    end if

    end function Ini_Read_String

    function Ini_Read_String_Default(Ini, Key, Default, AllowBlank, EnvDefault) result(AValue)
    !Returns from Ini first; if not found tries environment variable if EnvDefault is true,
    !If not found returns Default (unless not present, in which case gives error)
    !AllBlank specifies whether an empty string counts as a valid result to give instead of Default
    class(TIniFile) :: Ini
    character (LEN=*), intent(IN) :: Key
    character (LEN=*), intent(IN), optional ::Default
    logical, intent(in), optional :: AllowBlank
    logical, optional, intent(IN) :: EnvDefault
    character(LEN=:), allocatable :: AValue
    logical is_present

    if (Ini%HasKey(Key)) then
        AValue = Ini%Read_String(Key, .false.)
        if (AValue/='' .or. DefaultFalse(AllowBlank)) return
    else
        AValue=''
    end if
    if (DefaultFalse(EnvDefault)) then
        AValue = GetEnvironmentVariable(Key,is_present)
        if (DefaultFalse(AllowBlank) .and. is_present) return
    end if
    if (AValue=='') then
        if (.not. present(Default)) call Ini%Error('key not found',Key)
        AValue = Default
    end if
    call  Ini%ReadValues%Add(Key, AValue)
    if (Ini%Echo_Read) write (*,*) trim(Key)//' = ',trim(AValue)

    end function Ini_Read_String_Default

    function Ini_TestEqual(Ini, Key, value, EmptyOK) result(OK)
    class(TIniFile) :: Ini
    character (LEN=*), intent(IN) :: Key, value
    logical, intent(in), optional :: EmptyOK
    character(LEN=:), pointer :: val
    logical :: OK

    val => Ini%ValueOf(Key)
    OK = val == value
    if (.not. OK .and. present(EmptyOK)) then
        OK = EmptyOK .and. val==''
    end if

    end function Ini_TestEqual

    function Ini_Key_To_Arraykey(Ini,Key, index)  result(AValue)
    class(TIniFile), intent(in) :: Ini
    character (LEN=*), intent(IN) :: Key
    integer, intent(in) :: index
    character(LEN=:), allocatable :: AValue
    character(LEN=32) :: numstr

    write (numstr,*) index
    numstr=adjustl(numstr)
    AValue = trim(Key) // '(' // trim(numStr) // ')'

    end function Ini_Key_To_Arraykey

    function Ini_Read_String_Array(Ini, Key, index, NotFoundFail) result(AValue)
    class(TIniFile) :: Ini
    integer, intent(in) :: index
    character (LEN=*), intent(IN) :: Key
    logical, optional, intent(IN) :: NotFoundFail
    character(LEN=:), pointer :: AValue
    character(LEN=:), allocatable :: ArrayKey

    ArrayKey = Ini%Key_To_Arraykey(Key,index)
    AValue => Ini%Read_String(ArrayKey, NotFoundFail)

    end function Ini_Read_String_Array

    function Ini_Read_Int_Array(Ini,Key, index, Default)
    !Reads Key(1), Key(2), etc.
    class(TIniFile) :: Ini
    integer Ini_Read_Int_Array
    integer, optional, intent(IN) :: Default
    integer, intent(in) :: index
    character(LEN=*), intent(IN) :: Key
    character(LEN=:), allocatable :: ArrayKey

    ArrayKey = Ini%Key_To_Arraykey(Key,index)
    Ini_Read_Int_Array = Ini%Read_Int(ArrayKey, Default)

    end function Ini_Read_Int_Array

    function Ini_Read_Int(Ini, Key, Default, min, max, OK)
    class(TIniFile) :: Ini
    integer Ini_Read_Int
    integer, optional, intent(IN) :: Default, min, max
    logical, intent(out), optional :: OK
    character(LEN=*), intent(IN) :: Key
    character(LEN=:), pointer :: S
    integer status

    S => Ini%Read_String(Key,.not. present(Default))
    if (S == '') then
        call Ini%EmptyCheckDefault(Key,Default)
        Ini_Read_Int = Default
        call  Ini%ReadValues%Add(Key, Default)
    else
        if (verify(trim(S),'-+0123456789') /= 0) then
            status=1
            if (present(OK)) then
                Ini_Read_Int=-1
                OK = .false.
                return
            end if
        else
            read (S,*, iostat=status) Ini_Read_Int
        end if
        if (status/=0) call Ini%Error('error reading integer',Key)
        if (present(max)) then
            if (Ini_Read_Int > max) call Ini%Error('value > max',Key)
        end if
        if (present(min)) then
            if (Ini_Read_Int < min) call Ini%Error('value < min',Key)
        end if
    end if
    if (present(OK)) OK = .true.

    end function Ini_Read_Int

    integer function Ini_EnumerationValue(Ini, S, Names)
    class(TIniFile) :: Ini
    character(LEN=*) :: S
    character(LEN=Ini_Enumeration_Len), intent(in) :: Names(:)
    integer i

    do i=1, size(Names)
        if (S==Names(i)) then
            Ini_EnumerationValue = i
            return
        end if
    end do
    Ini_EnumerationValue = -1

    end function Ini_EnumerationValue

    function Ini_Read_Enumeration(Ini, Key, Names, Default)
    class(TIniFile) :: Ini
    character(LEN=*), intent(in) :: Key
    character(LEN=Ini_Enumeration_Len), intent(in) :: Names(:)
    integer, optional, intent(in) :: Default
    integer Ini_Read_Enumeration
    character(LEN=:), pointer :: S
    logical OK
    integer i

    Ini_Read_Enumeration = Ini%Read_Int(Key, Default, OK = OK)
    if (OK) then
        if (Ini_Read_Enumeration<1 .or. Ini_Read_Enumeration> size(Names)) &
            & call Ini%Error('enumeration value not valid',Key)
    else
        S => Ini%ValueOf(Key)
        Ini_Read_Enumeration = Ini%EnumerationValue(S, Names)
        if (Ini_Read_Enumeration<0) call Ini%Error('"'//S//'" enumeration name not recognised',Key)
    end if

    end function Ini_Read_Enumeration

    subroutine Ini_Read_Enumeration_List(Ini, Key, Names, Enums, nvalues, max_enums, Default)
    class(TIniFile) :: Ini
    character(LEN=*), intent(in) :: Key
    character(LEN=Ini_Enumeration_Len), intent(in) :: Names(:)
    integer, allocatable, intent(out) :: Enums(:)
    character(LEN=*), optional, intent(in) :: Default
    integer, intent(in), optional :: nvalues, max_enums
    integer, parameter :: defmax = 128
    integer :: maxvalues
    integer, allocatable :: Values(:)
    character(LEN=Ini_Enumeration_Len+8) part
    character(LEN=:), allocatable :: InLine
    integer i, slen, pos, n, status

    InLine = Ini%Read_String_Default(Key, Default)
    pos = 1
    slen = len_trim(InLine)
    n=0
    maxvalues = PresentDefault(PresentDefault(defmax, max_enums), nvalues)
    allocate(Values(maxvalues))
    do
        do while (pos <= slen)
            if (IsWhiteSpace(InLine(pos:pos))) then
                pos = pos+1
            else
                exit
            endif
        end do
        read(InLine(pos:), *, iostat=status) part
        if (status/=0) exit
        pos = pos + len_trim(part)
        i= Ini%EnumerationValue(part, Names)
        if (i<0) call Ini%Error('"'//part//'" enumeration name not recognised',Key)
        n=n+1
        if (n > defmax) call Ini%Error('More than maximum enumeration values', Key)
        values(n) = i
        if (n == maxvalues) exit
    end do
    if (present(nvalues)) then
        if (n==1) then
            !Fill whole array with same value
            allocate(Enums(nvalues), source= values(1))
            return
        elseif (n/=nvalues) then
            call Ini%Error('Wrong number of enumeration values', Key)
        end if
    end if
    allocate(Enums(n), source= values(:n))

    end subroutine Ini_Read_Enumeration_List

    function Ini_Read_Double(Ini,Key, Default, min, max)
    class(TIniFile) :: Ini
    double precision Ini_Read_Double
    double precision, optional, intent(IN) :: Default, min,max
    character (LEN=*), intent(IN) :: Key
    character(LEN=:), pointer :: S
    integer status

    S => Ini%Read_String(Key,.not. present(Default))
    if (S == '') then
        call Ini%EmptyCheckDefault(Key,Default)
        Ini_Read_Double = Default
        call  Ini%ReadValues%Add(Key, Default)
    else
        read (S,*, iostat=status) Ini_Read_Double
        if (status/=0) call Ini%Error('error reading double',Key)
    end if
    if (present(max)) then
        if (Ini_Read_Double > max) call Ini%Error('value > max',Key)
    end if
    if (present(min)) then
        if (Ini_Read_Double < min) call Ini%Error('value < min',Key)
    end if

    end function Ini_Read_Double

    function Ini_Read_Double_Array(Ini,Key, index, Default, min, max)
    !Reads Key(1), Key(2), etc.
    class(TIniFile) :: Ini
    double precision Ini_Read_Double_Array
    double precision, optional, intent(IN) :: Default, min, max
    integer, intent(in) :: index
    character(LEN=*), intent(IN) :: Key
    character(LEN=:), allocatable :: ArrayKey

    ArrayKey = Ini%Key_To_Arraykey(Key,index)
    Ini_Read_Double_Array = Ini%Read_Double(ArrayKey, Default, min, max)

    end function Ini_Read_Double_Array


    function Ini_Read_Real(Ini,Key, Default, min, max)
    class(TIniFile) :: Ini
    real Ini_Read_Real
    real, optional, intent(IN) :: Default, min, max
    character(LEN=*), intent(IN) :: Key
    character(LEN=:), pointer :: S
    integer status

    S => Ini%Read_String(Key,.not. present(Default))
    if (S == '') then
        call Ini%EmptyCheckDefault(Key,Default)
        Ini_Read_Real = Default
        call  Ini%ReadValues%Add(Key, Default)
    else
        read (S,*, iostat=status) Ini_Read_Real
        if (status/=0) call Ini%Error('error reading real',Key)
    end if
    if (present(max)) then
        if (Ini_Read_Real > max) call Ini%Error('value > max',Key)
    end if
    if (present(min)) then
        if (Ini_Read_Real < min) call Ini%Error('value < min',Key)
    end if

    end function Ini_Read_Real


    function Ini_Read_Real_Array(Ini,Key, index, Default, min, max)
    !Reads Key(1), Key(2), etc.
    class(TIniFile) :: Ini
    real Ini_Read_Real_Array
    real, optional, intent(IN) :: Default, min, max
    integer, intent(in) :: index
    character(LEN=*), intent(IN) :: Key
    character(LEN=:), allocatable :: ArrayKey

    ArrayKey = Ini%Key_To_Arraykey(Key,index)
    Ini_Read_Real_Array = Ini%Read_Real(ArrayKey, Default,min,max)

    end function Ini_Read_Real_Array

    function Ini_Read_Logical(Ini, Key, Default)
    class(TIniFile) :: Ini
    logical Ini_Read_Logical
    logical, optional, intent(IN) :: Default
    character(LEN=*), intent(IN) :: Key
    character(LEN=:), pointer :: S
    integer status

    S => Ini%Read_String(Key,.not. present(Default))
    if (S == '') then
        call Ini%EmptyCheckDefault(Key,Default)
        Ini_Read_Logical = Default
        call  Ini%ReadValues%Add(Key, Default)
    else
        if (verify(trim(S),'10TF') /= 0) then
            status=1
        else
            read (S,*, iostat=status) Ini_Read_Logical
        end if
        if (status/=0) call Ini%Error('error reading logical',Key)
    end if

    end function Ini_Read_Logical


    subroutine Ini_SaveReadValues(Ini, afile)
    class(TIniFile) :: Ini
    character(LEN=*), intent(in) :: afile
    integer i
    Type(TTextFile) :: F

    call F%CreateFile(afile)

    do i=1, Ini%ReadValues%Count
        call F%Write(Ini%ReadValues%Items(i)%P%Name, ' = ', Ini%ReadValues%Items(i)%P%Value)
    end do

    call F%Close()

    end subroutine Ini_SaveReadValues


    subroutine Ini_Read_Int_Change(Ini, Key, Current, min, max)
    class(TIniFile) :: Ini
    character(LEN=*), intent(IN) :: Key
    integer, intent(inout) :: Current
    integer, optional, intent(IN) :: min, max
    Current = Ini%Read_Int(Key,Current,min,max)
    end subroutine Ini_Read_Int_Change

    subroutine Ini_Read_Double_Change(Ini,Key, Current, min, max)
    class(TIniFile) :: Ini
    character(LEN=*), intent(IN) :: Key
    double precision, intent(inout) :: Current
    double precision, optional, intent(IN) :: min,max
    Current = Ini%Read_Double(Key, Current, min, max)
    end subroutine Ini_Read_Double_Change

    subroutine Ini_Read_Real_Change(Ini,Key, Current, min, max)
    class(TIniFile) :: Ini
    character(LEN=*), intent(IN) :: Key
    real, intent(inout) :: Current
    real, optional, intent(IN) :: min,max
    Current = Ini%Read_Real(Key, Current, min, max)
    end subroutine Ini_Read_Real_Change

    subroutine Ini_Read_Logical_Change(Ini,Key, Current)
    class(TIniFile) :: Ini
    character(LEN=*), intent(IN) :: Key
    logical, intent(inout) :: Current
    Current = Ini%Read_Logical(Key, Current)
    end subroutine Ini_Read_Logical_Change

    subroutine Ini_Read_String_Change(Ini,Key, Current)
    class(TIniFile) :: Ini
    character(LEN=*), intent(IN) :: Key
    character(LEN=:), allocatable, intent(inout) :: Current

    Current = Ini%Read_String_Default(Key, Current)

    end subroutine Ini_Read_String_Change

    end module IniObjects
