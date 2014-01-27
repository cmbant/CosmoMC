    !Module to read in name/value pairs from a file, with each line of the form line 'name = value'
    !Should correctly interpret FITS headers
    !Antony Lewis (http://cosmologist.info/). Released to the public domain.
    !2014 using Fortran 2003

    module IniObjects
    use FileUtils
    implicit none
    public

    integer, parameter :: Ini_max_name_len = 128

    integer, parameter :: Ini_max_string_len = 1024
    logical :: Ini_fail_on_not_found = .false.

    logical :: Ini_Echo_Read = .false.

    type TNameValue
        !no known way to make character string pointers..
        character(Ini_max_name_len)  :: Name
        character(Ini_max_string_len):: Value
    end type TNameValue

    type TNameValue_pointer
        Type(TNameValue), pointer :: P
    end type TNameValue_pointer

    Type TNameValueList
        integer Count
        integer Delta
        integer Capacity
        logical ignoreDuplicates
        logical :: Ini_AllowDuplicateKeys = .false.
        type(TNameValue_pointer), dimension(:), allocatable :: Items
    contains
    procedure :: Init => TNameValueList_Init
    procedure :: Clear => TNameValueList_Clear
    procedure :: ValueOf => TNameValueList_ValueOf
    procedure :: IndexOf => TNameValueList_IndexOf
    procedure :: HasKey => TNameValueList_HasKey
    procedure :: Add => TNameValueList_Add
    procedure :: SetCapacity => TNameValueList_SetCapacity
    procedure :: Delete => TNameValueList_Delete
    end Type TNameValueList

    Type TIniFile
        logical :: SlashComments = .false.
        Type(TNameValueList) :: L, ReadValues
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
    procedure :: HasKey => Ini_HasKey
    procedure :: Key_To_Arraykey => Ini_Key_To_Arraykey
    procedure :: NameValue_Add => Ini_NameValue_Add
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
    call TNameValueList_Init(L)

    end subroutine TNameValueList_Clear

    subroutine TNameValueList_ValueOf(L, AName, AValue)
    class(TNameValueList), intent(in) :: L
    character(LEN=*), intent(in) :: AName
    character(LEN=*), intent(out) :: AValue
    integer i

    do i=1, L%Count
        if (L%Items(i)%P%Name == AName) then
            AValue = L%Items(i)%P%Value
            return
        end if
    end do
    AValue = ''

    end subroutine TNameValueList_ValueOf


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
    if ((.not. L%Ini_AllowDuplicateKeys .or. isDefault) .and. L%HasKey(AName)) then
        if (L%ignoreDuplicates .or. isDefault) return
        write (*,*) 'IniFile,TNameValueList_Add: duplicate key name in file: '//trim(AName)
        stop
    end if
    if (L%Count == L%Capacity) call L%SetCapacity(L%Capacity + L%Delta)
    L%Count = L%Count + 1
    allocate(L%Items(L%Count)%P)
    L%Items(L%Count)%P%Name = AName
    L%Items(L%Count)%P%Value = AValue

    end subroutine TNameValueList_Add

    subroutine TNameValueList_SetCapacity(L, C)
    class(TNameValueList) :: L
    integer C
    type(TNameValue_pointer), dimension(:), pointer :: TmpItems

    if (L%Count > 0) then
        if (C < L%Count) stop 'TNameValueList_SetCapacity: smaller than Count'
        allocate(TmpItems(L%Count))
        TmpItems = L%Items(1:L%Count)
        deallocate(L%Items)
        allocate(L%Items(C))
        L%Items(1:L%Count) = TmpItems
        deallocate(TmpItems)
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
    if (EqPos/=0 .and. InLine(1:1)/='#' .and. InLine(1:7) /= 'COMMENT' ) then
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
        call TNameValueList_Add(Ini%L, AName, S,only_if_undefined = isDefault )
    end if

    end subroutine Ini_NameValue_Add

    function Ini_ExtractFilePath(aname)
    character(LEN=*), intent(IN) :: aname
    character(LEN=Ini_max_string_len) Ini_ExtractFilePath
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

    recursive subroutine Ini_Open(Ini, filename, &
    error, slash_comments, append,only_if_undefined)
    class(TIniFile) :: Ini
    character (LEN=*), intent(IN) :: filename
    integer :: unit_id
    logical, intent(OUT) :: error
    logical, optional, intent(IN) :: slash_comments
    logical, optional, intent(in) :: append, only_if_undefined
    character (LEN=Ini_max_string_len) :: InLine, IncludeFile
    integer lastpos, i
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
        call Ini%L%Init()
        call Ini%ReadValues%Init(.true.)
    end if

    call IncudeFiles%Init()
    call DefaultValueFiles%Init()

    if (present(slash_comments)) then
        Ini%SlashComments = slash_comments
    else
        Ini%SlashComments = .false.
    end if

    open(newunit=unit_id,file=filename,form='formatted',status='old', err=500)

    do
        if (.not. ReadLine(unit_id,InLine)) exit
        if (InLine == 'END') exit
        if (InLine(1:8) == 'INCLUDE(') then
            lastpos = scan(InLine,')')
            if (lastpos/=0) then
                call TNameValueList_Add(IncudeFiles, trim(adjustl(InLine(9:lastpos-1))),'')
            else
                stop 'Ini_Open: error in INCLUDE line'
            end if
        elseif (InLine(1:8) == 'DEFAULT(') then
            !Settings to read in as defaults, overridden by subsequent re-definitions
            lastpos = scan(InLine,')')
            if (lastpos/=0) then
                call TNameValueList_Add(DefaultValueFiles, trim(adjustl(InLine(9:lastpos-1))),'')
            else
                stop 'Ini_Open: error in DEFAULT line'
            end if
        elseif (InLine /= '') then
            call Ini%NameValue_Add(InLine, only_if_undefined=isDefault)
        end if
    end do

    close(unit_id)
    error=.false.

    do i=1, IncudeFiles%Count
        if (error) exit
        IncludeFile=IncudeFiles%Items(i)%P%Name
        inquire(file=IncludeFile, exist = FileExists)
        if (.not. FileExists) then
            IncludeFile=trim(Ini_ExtractFilePath(filename))//trim(IncludeFile)
            inquire(file=IncludeFile, exist = FileExists)
            if (.not. FileExists) then
                write(*,*) 'Ini_Open: INCLUDE file not found: '//trim(IncudeFiles%Items(i)%P%Name)
                stop
            end if
        end if
        call Ini%Open(IncludeFile, error, slash_comments, append=.true.,only_if_undefined=isDefault)
    end do
    do i=1, DefaultValueFiles%Count
        if (error) exit
        IncludeFile=DefaultValueFiles%Items(i)%P%Name
        inquire(file=IncludeFile, exist = FileExists)
        if (.not. FileExists) then
            IncludeFile=trim(Ini_ExtractFilePath(filename))//trim(IncludeFile)
            inquire(file=IncludeFile, exist = FileExists)
            if (.not. FileExists) then
                write(*,*) 'Ini_Open: DEFAULT file not found:' //trim(DefaultValueFiles%Items(i)%P%Name)
                stop
            end if
        end if
        call Ini%Open(IncludeFile, error, slash_comments, append=.true., only_if_undefined=.true.)
    end do

    call TNameValueList_Clear(IncudeFiles)
    call TNameValueList_Clear(DefaultValueFiles)

    return

500 error=.true.
    call TNameValueList_Clear(IncudeFiles)
    call TNameValueList_Clear(DefaultValueFiles)

    end subroutine Ini_Open

    subroutine Ini_Open_Fromlines(Ini, Lines, NumLines, slash_comments)
    class(TIniFile) :: Ini

    integer, intent(IN) :: NumLines
    character (LEN=*), dimension(NumLines), intent(IN) :: Lines
    logical, intent(IN) :: slash_comments
    integer i

    call TNameValueList_Init(Ini%L)
    call TNameValueList_Init(Ini%ReadValues, .true.)

    Ini%SlashComments = slash_comments

    do i=1,NumLines
        call Ini%NameValue_Add(Lines(i))
    end do

    end  subroutine Ini_Open_Fromlines


    subroutine Ini_Close(Ini)
    class(TIniFile) :: Ini

    call TNameValueList_Clear(Ini%L)
    call TNameValueList_Clear(Ini%ReadValues)

    end  subroutine Ini_Close

    function Ini_Read_String(Ini, Key, NotFoundFail) result(AValue)
    class(TIniFile) :: Ini
    character (LEN=*), intent(IN) :: Key
    logical, optional, intent(IN) :: NotFoundFail
    character(LEN=Ini_max_string_len) :: AValue

    call TNameValueList_ValueOf(Ini%L, Key, AValue)

    if (AValue/='') then
        call  TNameValueList_Add(Ini%ReadValues, Key, AValue)
        if (Ini_Echo_Read) write (*,*) trim(Key)//' = ',trim(AValue)
        return
    end if
    if (present(NotFoundFail)) then
        if (NotFoundFail) then
            write(*,*) 'key not found : '//trim(Key)
            stop
        end if
    else if (Ini_fail_on_not_found) then
        write(*,*) 'key not found : '//trim(Key)
        stop
    end if

    end function Ini_Read_String

    function Ini_Read_String_Default(Ini, Key, Default, AllowBlank) result(AValue)
    class(TIniFile) :: Ini
    character (LEN=*), intent(IN) :: Key, Default
    character(LEN=Ini_max_string_len) :: AValue
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

    function Ini_HasKey(Ini, Key) result(AValue)
    class(TIniFile), intent(in) :: Ini
    character (LEN=*), intent(IN) :: Key
    logical AValue

    Avalue = TNameValueList_HasKey(Ini%L, Key)

    end function Ini_HasKey

    function Ini_Key_To_Arraykey(Ini,Key, index)  result(AValue)
    class(TIniFile), intent(in) :: Ini
    character (LEN=*), intent(IN) :: Key
    integer, intent(in) :: index
    character(LEN=Ini_max_string_len) :: AValue
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
    character(LEN=Ini_max_string_len) :: AValue
    character(LEN=Ini_max_string_len) :: ArrayKey

    ArrayKey = Ini%Key_To_Arraykey(Key,index)
    AValue = Ini%Read_String(ArrayKey, NotFoundFail)

    end function Ini_Read_String_Array

    function Ini_Read_Int_Array(Ini,Key, index, Default)
    !Reads Key(1), Key(2), etc.
    class(TIniFile) :: Ini
    integer Ini_Read_Int_Array
    integer, optional, intent(IN) :: Default
    integer, intent(in) :: index
    character (LEN=*), intent(IN) :: Key
    character(LEN=Ini_max_string_len) :: ArrrayKey

    ArrrayKey = Ini%Key_To_Arraykey(Key,index)
    Ini_Read_Int_Array = Ini%Read_Int(ArrrayKey, Default)

    end function Ini_Read_Int_Array

    function Ini_Read_Int(Ini, Key, Default)
    class(TIniFile) :: Ini
    integer Ini_Read_Int
    integer, optional, intent(IN) :: Default
    character  (LEN=*), intent(IN) :: Key
    character(LEN=Ini_max_string_len) :: S

    S = Ini%Read_String(Key,.not. present(Default))
    if (S == '') then
        if (.not. present(Default)) then
            write(*,*) 'no value for key: '//Key
            stop
        end if
        Ini_Read_Int = Default
        write (S,*) Default
        call  TNameValueList_Add(Ini%ReadValues, Key, S)
    else
        if (verify(trim(S),'-+0123456789') /= 0) goto 10
        read (S,*, err = 10) Ini_Read_Int
    end if
    return
10  write (*,*) 'error reading integer for key: '//Key
    stop

    end function Ini_Read_Int

    function Ini_Read_Double(Ini,Key, Default)
    class(TIniFile) :: Ini
    double precision Ini_Read_Double
    double precision, optional, intent(IN) :: Default
    character (LEN=*), intent(IN) :: Key
    character(LEN=Ini_max_string_len) :: S

    S = Ini%Read_String(Key,.not. present(Default))
    if (S == '') then
        if (.not. present(Default)) then
            write(*,*) 'no value for key: '//Key
            stop
        end if
        Ini_Read_Double = Default
        write (S,*) Default

        call  TNameValueList_Add(Ini%ReadValues, Key, S)
    else
        read (S,*, err=10) Ini_Read_Double
    end if

    return

10  write (*,*) 'error reading double for key: '//Key
    stop

    end function Ini_Read_Double

    function Ini_Read_Double_Array(Ini,Key, index, Default)
    !Reads Key(1), Key(2), etc.
    class(TIniFile) :: Ini
    double precision Ini_Read_Double_Array
    double precision, optional, intent(IN) :: Default
    integer, intent(in) :: index
    character (LEN=*), intent(IN) :: Key
    character(LEN=Ini_max_string_len) ::  ArrrayKey

    ArrrayKey = Ini%Key_To_Arraykey(Key,index)
    Ini_Read_Double_Array = Ini%Read_Double(ArrrayKey, Default)

    end function Ini_Read_Double_Array


    function Ini_Read_Real(Ini,Key, Default)
    class(TIniFile) :: Ini
    real Ini_Read_Real
    real, optional, intent(IN) :: Default
    character (LEN=*), intent(IN) :: Key
    character(LEN=Ini_max_string_len) :: S

    S = Ini%Read_String(Key,.not. present(Default))
    if (S == '') then
        if (.not. present(Default)) then
            write(*,*) 'no value for key: '//Key
            stop
        end if
        Ini_Read_Real = Default
        write (S,*) Default
        call  TNameValueList_Add(Ini%ReadValues, Key, S)
    else
        read (S,*, err=10) Ini_Read_Real
    end if

    return

10  write (*,*) 'error reading double for key: '//Key
    stop

    end function Ini_Read_Real


    function Ini_Read_Real_Array(Ini,Key, index, Default)
    !Reads Key(1), Key(2), etc.
    class(TIniFile) :: Ini
    real Ini_Read_Real_Array
    real, optional, intent(IN) :: Default
    integer, intent(in) :: index
    character (LEN=*), intent(IN) :: Key
    character(LEN=Ini_max_string_len) :: ArrrayKey

    ArrrayKey = Ini%Key_To_Arraykey(Key,index)
    Ini_Read_Real_Array = Ini%Read_Real(ArrrayKey, Default)

    end function Ini_Read_Real_Array

    function Ini_Read_Logical(Ini, Key, Default)
    class(TIniFile) :: Ini
    logical Ini_Read_Logical
    logical, optional, intent(IN) :: Default
    character  (LEN=*), intent(IN) :: Key

    character(LEN=Ini_max_string_len) :: S

    S = Ini%Read_String(Key,.not. present(Default))
    if (S == '') then
        if (.not. present(Default)) then
            write(*,*) 'no value for key: '//Key
            stop
        end if
        Ini_Read_Logical = Default
        write (S,*) Default

        call  TNameValueList_Add(Ini%ReadValues, Key, S)
    else
        if (verify(trim(S),'10TF') /= 0) goto 10
        read (S,*, err = 10) Ini_Read_Logical
    end if

    return

10  write (*,*) 'error reading logical for key: '//Key
    stop
    end function Ini_Read_Logical


    subroutine Ini_SaveReadValues(Ini, afile)
    class(TIniFile) :: Ini
    character(LEN=*), intent(in) :: afile
    integer :: unit_id
    integer i

    unit_id = CreateNewTxtFile(afile)

    do i=1, Ini%ReadValues%Count
        write (unit_id,'(a)') trim(Ini%ReadValues%Items(i)%P%Name) // ' = ' &
        //trim(Ini%ReadValues%Items(i)%P%Value)
    end do

    close(unit_id)

    end subroutine Ini_SaveReadValues

    end module IniObjects
