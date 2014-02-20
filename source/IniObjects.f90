    !Module to read in name/value pairs from a file, with each line of the form line 'name = value'
    !Should correctly interpret FITS headers
    !Antony Lewis (http://cosmologist.info/). Released to the public domain.
    !2014 using Fortran 2003

    module IniObjects
    use FileUtils
    implicit none
    public

    integer, parameter :: Ini_max_string_len = 1024*4

    character(LEN=0), target :: Empty_String = ''

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
    procedure :: TNameValueList_AddDouble
    procedure :: TNameValueList_AddReal
    procedure :: TNameValueList_AddInt
    procedure :: TNameValueList_AddLogical
    generic :: Add => AddString, TNameValueList_AddDouble, TNameValueList_AddReal, TNameValueList_AddInt, TNameValueList_AddLogical
    end Type TNameValueList

    Type, extends(TNameValueList) :: TIniFile
        logical :: SlashComments = .false.
        logical :: Echo_Read = .false.
        logical :: Fail_on_not_found = .false.
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
    procedure :: SaveReadValues => Ini_SaveReadValues
    procedure :: Key_To_Arraykey => Ini_Key_To_Arraykey
    procedure :: NameValue_Add => Ini_NameValue_Add
    procedure :: EmptyCheckDefault => Ini_EmptyCheckDefault
    procedure, nopass :: ExtractFilePath => Ini_ExtractFilePath
    procedure :: Ini_Read_Real_Change
    procedure :: Ini_Read_Double_Change
    procedure :: Ini_Read_Int_Change
    procedure :: Ini_Read_Logical_Change
    generic :: Read => Ini_Read_Real_Change, Ini_Read_Double_Change, Ini_Read_Int_Change, Ini_Read_Logical_Change
    end Type TIniFile

    contains

    subroutine TNameValueList_Init(L, ignoreDuplicates)
    class(TNameValueList) :: L
    logical, intent(in), optional :: ignoreDuplicates

    L%Count = 0
    L%Capacity = 0
    L%Delta = 128
    L%ignoreDuplicates = .false.
    if (present(ignoreDuplicates)) L%ignoreDuplicates=ignoreDuplicates

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

    function TNameValueList_ValueOf(L, AName) result(AValue)
    class(TNameValueList), intent(in) :: L
    character(LEN=*), intent(in) :: AName
    character(LEN=:), pointer :: AValue
    integer i

    do i=1, L%Count
        if (L%Items(i)%P%Name == AName) then
            AValue => L%Items(i)%P%Value
            return
        end if
    end do
    AValue => Empty_String

    end function TNameValueList_ValueOf


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

    subroutine TNameValueList_Add(L, AName, AValue, only_if_undefined)
    class(TNameValueList) :: L
    character(LEN=*), intent(in) :: AName, AValue
    logical, optional, intent(in) :: only_if_undefined
    logical isDefault

    if (present(only_if_undefined)) then
        isDefault = only_if_undefined
    else
        isDefault = .true.
    end if
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
        write(*,*) 'Error for key '//trim(Key)//' : '//Msg
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



    subroutine Ini_NameValue_Add(Ini,AInLine,only_if_undefined)
    class(TIniFile) :: Ini
    character (LEN=*), intent(IN) :: AInLine
    integer EqPos, slashpos, lastpos
    logical, optional, intent(in) :: only_if_undefined
    logical isDefault
    character (LEN=len(AInLine)) :: AName, S, InLine

    if (present(only_if_undefined)) then
        isDefault = only_if_undefined
    else
        isDefault = .false.
    end if
    InLine=trim(adjustl(AInLine))
    EqPos = scan(InLine,'=')
    if (EqPos/=0 .and. InLine(1:1)/='#' .and. InLine(1:min(len(InLine),7)) /= 'COMMENT' ) then
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
        call Ini%TNameValueList%Add(AName, S,only_if_undefined = isDefault)
    end if

    end subroutine Ini_NameValue_Add

    function Ini_ExtractFilePath(aname)
    character(LEN=*), intent(IN) :: aname
    character(LEN=:), allocatable :: Ini_ExtractFilePath
    integer len, i

    len = len_trim(aname)
    do i = len, 1, -1
        if (aname(i:i)=='/') then
            Ini_ExtractFilePath = aname(1:i)
            return
        end if
    end do
    Ini_ExtractFilePath = ''

    end function Ini_ExtractFilePath

    recursive subroutine Ini_Open(Ini, filename, error, slash_comments, append,only_if_undefined)
    class(TIniFile) :: Ini
    character (LEN=*), intent(IN) :: filename
    integer :: unit_id
    logical, intent(OUT), optional :: error
    logical, optional, intent(IN) :: slash_comments
    logical, optional, intent(in) :: append, only_if_undefined
    character (LEN=Ini_max_string_len) :: IncludeFile
    character(LEN=:), allocatable :: InLine
    integer lastpos, i, status
    Type(TNameValueList) IncudeFiles, DefaultValueFiles
    logical doappend, FileExists, isDefault

    if (present(append)) then
        doappend=append
    else
        doappend=.false.
    end if

    if (present(only_if_undefined)) then
        isDefault = only_if_undefined
    else
        isDefault = .false.
    end if

    if (.not. doappend) then
        call Ini%TNameValueList%Init()
        call Ini%ReadValues%Init(.true.)
    end if

    if (present(slash_comments)) then
        Ini%SlashComments = slash_comments
    else
        Ini%SlashComments = .false.
    end if

    open(newunit=unit_id,file=filename,form='formatted',status='old', iostat=status)
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
        if (.not. ReadLine(unit_id,InLine)) exit
        if (InLine == 'END') exit
        if (InLine(1:min(8,len(InLine))) == 'INCLUDE(') then
            lastpos = scan(InLine,')')
            if (lastpos/=0) then
                call IncudeFiles%Add(InLine(9:lastpos-1),'')
            else
                call Ini%Error('Ini_Open, error in INCLUDE line: '//trim(filename))
            end if
        elseif (InLine(1:min(8,len(InLine))) == 'DEFAULT(') then
            !Settings to read in as defaults, overridden by subsequent re-definitions
            lastpos = scan(InLine,')')
            if (lastpos/=0) then
                call DefaultValueFiles%Add(InLine(9:lastpos-1),'')
            else
                call Ini%Error('Ini_Open, error in DEFAULT line: '//trim(filename))
            end if
        elseif (InLine /= '') then
            call Ini%NameValue_Add(InLine, only_if_undefined=isDefault)
        end if
    end do

    close(unit_id)
    if (present(error)) error=.false.

    do i=1, IncudeFiles%Count
        if (present(error)) then
            if (error) exit
        end if
        IncludeFile= IncudeFiles%Items(i)%P%Name
        inquire(file=IncludeFile, exist = FileExists)
        if (.not. FileExists) then
            IncludeFile= Ini%ExtractFilePath(filename)//trim(IncludeFile)
            inquire(file=IncludeFile, exist = FileExists)
            if (.not. FileExists) then
                call Ini%Error('Ini_Open, error in INCLUDE file not found: '//trim(IncudeFiles%Items(i)%P%Name))
            end if
        end if
        call Ini%Open(IncludeFile, error, slash_comments, append=.true.,only_if_undefined=isDefault)
    end do
    do i=1, DefaultValueFiles%Count
        if (present(error)) then
            if (error) exit
        end if
        IncludeFile=DefaultValueFiles%Items(i)%P%Name
        inquire(file=IncludeFile, exist = FileExists)
        if (.not. FileExists) then
            IncludeFile=trim(Ini%ExtractFilePath(filename))//trim(IncludeFile)
            inquire(file=IncludeFile, exist = FileExists)
            if (.not. FileExists) then
                call Ini%Error('Ini_Open, error in DEFAULT file not found: '//trim(IncudeFiles%Items(i)%P%Name))
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
        call Ini%NameValue_Add(Lines(i))
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
    if (present(NotFoundFail)) then
        if (NotFoundFail) then
            call Ini%Error('key not found',Key)
        end if
    else if (Ini%fail_on_not_found) then
        call Ini%Error('key not found',Key)
    end if

    end function Ini_Read_String

    function Ini_Read_String_Default(Ini, Key, Default, AllowBlank) result(AValue)
    class(TIniFile) :: Ini
    character (LEN=*), intent(IN) :: Key, Default
    character(LEN=:), allocatable :: AValue
    logical, intent(in), optional :: AllowBlank

    if (Ini%HasKey(Key)) then
        AValue = Ini%Read_String(Key, .false.)
        if (present(AllowBlank)) then
            if (AllowBlank) return
        end if
        if (AValue=='') AValue = Default
    else
        AValue = Default
    end if

    end function Ini_Read_String_Default


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

    function Ini_Read_Int(Ini, Key, Default, min, max)
    class(TIniFile) :: Ini
    integer Ini_Read_Int
    integer, optional, intent(IN) :: Default, min, max
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

    end function Ini_Read_Int

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
    integer :: unit_id
    integer i

    unit_id = CreateNewTxtFile(afile)

    do i=1, Ini%ReadValues%Count
        write (unit_id,'(a)') Ini%ReadValues%Items(i)%P%Name // ' = ' &
        & //Ini%ReadValues%Items(i)%P%Value
    end do

    close(unit_id)

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

    end module IniObjects