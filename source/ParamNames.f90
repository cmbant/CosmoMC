    module ParamNames
    use AMlUtils
    implicit none

    integer, parameter :: ParamNames_maxlen = 128

    Type TParamNames
        integer :: nnames =0
        integer :: num_MCMC = 0
        integer :: num_derived = 0
        character(LEN=ParamNames_maxlen), dimension(:), pointer ::  name=> null()
        character(LEN=ParamNames_maxlen), dimension(:), pointer ::  label => null()
        character(LEN=ParamNames_maxlen), dimension(:), pointer ::  comment => null()
        logical, dimension(:), pointer ::  is_derived
    end Type TParamNames

    contains

    function IsWhiteSpace(C)
    character, intent(in) :: C
    logical IsWhiteSpace

    IsWhiteSpace = (C==' ') .or. (C==char(9)) 

    end function IsWhiteSpace


    function ParamNames_ParseLine(Names,InLine,n) result(res)
    type(TParamNames) :: Names
    character(LEN=*) :: InLine
    integer n
    logical res
    integer pos, len

    len = len_trim(InLIne)
    pos =1
    do while (pos < len .and. IsWhiteSpace(InLIne(pos:pos))) 
        pos = pos+1
    end do 
    read(InLine(pos:), *, end=400, err=400) Names%name(n)
    pos = pos + len_trim(Names%name(n))
    do while (pos < len .and. IsWhiteSpace(InLIne(pos:pos))) 
        pos = pos+1
    end do 
    Names%label(n) = trim(adjustl(InLine(pos:len))) 
    pos = scan(Names%label(n),'#')
    if (pos/=0) then
        Names%comment(n) = Names%label(n)(pos+1: len_trim(Names%label(n)))    
        Names%label(n) = Names%label(n)(1:pos-1)
    else
        Names%comment(n) = ''         
    endif 
    pos = scan(Names%label(n),char(9))
    if (pos/=0) Names%label(n) = Names%label(n)(1:pos-1)      
    Names%name(n) = trim(adjustl(Names%name(n)))
    len = len_trim( Names%name(n) )
    if (Names%name(n)(len:len)=='*') then 
        Names%name(n)(len:len)=' ' 
        Names%is_derived(n) = .true.
    else
        Names%is_derived(n) = .false. 
    end if
    res = .true.
    return
400 res=.false.
    return

    end function ParamNames_ParseLine

    subroutine ParamNames_Alloc(Names,n)
    Type(TParamNames) :: Names
    integer,intent(in) :: n

    allocate(Names%name(n))
    allocate(Names%label(n)) 
    allocate(Names%comment(n))
    allocate(Names%is_derived(n))
    Names%nnames = n
    Names%is_derived = .false. 
    Names%num_MCMC = 0
    Names%num_derived = 0
    Names%name = '' 
    Names%comment=''
    Names%label=''

    end subroutine ParamNames_Alloc
    
    subroutine ParamNames_dealloc(Names)
        Type(TParamNames) :: Names
        integer i
        if (associated(Names%name)) &
         deallocate(Names%name,Names%label,Names%comment,Names%is_derived)
      
    end subroutine ParamNames_dealloc
    
    subroutine ParamNames_Init(Names, filename)
    Type(TParamNames) :: Names
    character(Len=*), intent(in) :: filename
    integer handle,n
    character (LEN=ParamNames_maxlen*3) :: InLine

    handle = new_file_unit()
    call OpenTxtFile(filename,handle)
    n = FileLines(handle)
    call ParamNames_Alloc(Names,n)

    n=0
    do 
        read (handle,'(a)',end=500) InLine
        if (trim(InLine)=='') cycle
        n=n+1
        if (.not. ParamNames_ParseLine(Names,InLine,n)) then
            call MpiStop(concat('ParamNames_Init: error parsing line: ',n))
        end if
    end do

500 call CloseFile(handle) 

    Names%nnames = n
    Names%num_derived = count(Names%is_derived)
    Names%num_MCMC = Names%nnames - Names%num_derived 


    end subroutine ParamNames_Init

    subroutine ParamNames_AssignItem(Names, Names2,n,i)
    Type(TParamNames), target :: Names, Names2
    integer n, i

    Names%name(n) = Names2%name(i)
    Names%label(n) = Names2%label(i)
    Names%comment(n) = Names2%comment(i)
    Names%is_derived(n) = Names2%is_derived(i)

    end subroutine ParamNames_AssignItem


    subroutine ParamNames_Add(Names, Names2)
    Type(TParamNames), target :: Names, Names2, Tmp
    integer n,i, newold, derived
    Type(TParamNames),pointer :: P

    n=0
    do i=1, names2%nnames
        if (ParamNames_index(Names,Names2%name(i))==-1) then
            n=n+1
        end if
    end do
    if (n==0) return

    call ParamNames_Alloc(Tmp,Names%nnames + n)
    Tmp%nnames = 0
    do derived=0,1
        P=> Names
        do newold=0,1
            do i=1, P%nnames
                if (ParamNames_index(Tmp,P%name(i))==-1) then
                    if (derived==0 .and. .not. P%is_derived(i) .or.derived==1 .and. P%is_derived(i) ) then
                        Tmp%nnames = Tmp%nnames + 1
                        call ParamNames_AssignItem(Tmp, P, Tmp%nnames , i)
                    end if
                end if
            end do
            P=> Names2
        enddo
    end do
    if (Tmp%nnames/= Names%nnames + n) stop 'ParamNames_Add: duplicate parameters?'

    Tmp%num_derived = count(Tmp%is_derived)
    Tmp%num_MCMC= Tmp%nnames-Tmp%num_derived
    call ParamNames_dealloc(Names)
    Names=Tmp

    end subroutine ParamNames_Add

    subroutine ParamNames_SetLabels(Names,filename)
    Type(TParamNames) :: Names, LabNames
    character(Len=*), intent(in) :: filename
    integer i,ix

    call ParamNames_init(LabNames,filename)
    do i=1, LabNames%nnames
        ix = ParamNames_index(Names, LabNames%name(i))
        if (ix/=-1) then
            Names%label(ix) = LabNames%label(i)
        end if
    end do

    end subroutine ParamNames_SetLabels

    function ParamNames_index(Names,name) result(ix)
    Type(TParamNames) :: Names
    character(len=*), intent(in) :: name
    integer ix,i

    do i=1,Names%nnames
        if (Names%name(i) == name) then
            ix = i
            return
        end if
    end do
    ix = -1

    end function ParamNames_index


    function ParamNames_label(Names,name) result(lab)
    Type(TParamNames) :: Names
    character(len=*), intent(in) :: name
    character(len = ParamNames_maxlen) lab
    integer ix

    ix = ParamNames_index(Names,name)
    if (ix>0) then
        lab = Names%label(ix)
    else
        lab = ''
    end if 

    end function ParamNames_label

    function ParamNames_name(Names,ix) result(name)
    Type(TParamNames) :: Names
    character(len=ParamNames_maxlen)  :: name
    integer, intent(in) :: ix

    if (ix <= Names%nnames) then
        name = Names%name(ix)
    else
        name = ''
    end if  

    end function ParamNames_name


    subroutine ParamNames_ReadIndices(Names,InLine, params, num)
    Type(TParamNames) :: Names
    character(LEN=*), intent(in) :: InLine
    integer, intent(out) :: params(*)
    integer  :: num
    character(LEN=ParamNames_maxlen) part
    integer param,len,ix, pos, max_num
    integer, parameter :: unknown_num = 1024

    if (num==0) return
    len = len_trim(InLine)
    pos = 1
    if (num==-1) then
        max_num = unknown_num 
    else
        max_num = num
    end if
    do param = 1, max_num
        do while (pos < len .and. IsWhiteSpace(InLine(pos:pos)) ) 
            pos = pos+1
        end do 
        read(InLine(pos:), *, end=400, err=400) part
        if (max_num == unknown_num) num = param
        pos = pos + len_trim(part)
        ix = ParamNames_index(Names,part)
        if (ix>0) then
            params(param) = ix
        else
            if (verify(trim(part),'0123456789') /= 0) then
                call MpiStop( 'ParamNames: Unknown parameter name '//trim(part))
            end if
            read(part,*) params(param)
        end if
    end do 
    return
400 if (max_num==unknown_num) return
    call MpiStop('ParamNames: Not enough names or numbers - '//trim(InLine))

    end subroutine ParamNames_ReadIndices

    function ParamNames_AsString(Names, i, want_comment) result(line)
    Type(TParamNames) :: Names
    integer, intent(in) :: i
    logical ,intent(in), optional :: want_comment
    character(LEN=ParamNames_maxlen*3) Line
    logical wantCom

    if (present(want_comment)) then
        wantCom = want_comment
    else
        wantCom = .false.
    end if

    if (i> Names%nnames) call MpiStop('ParamNames_AsString: index out of range')
    Line = trim(Names%name(i))
    if (Names%is_derived(i))Line = concat(Line,'*')
    Line =  trim(Line)//char(9)//trim(Names%label(i))
    if (wantCom .and. Names%comment(i)/='') then
        Line = trim(Line)//char(9)//'#'//trim(Names%comment(i))
    end if

    end function ParamNames_AsString

    subroutine ParamNames_WriteFile(Names, fname)
    Type(TParamNames) :: Names
    character(LEN=*), intent(in) :: fname
    integer :: unit
    integer i

    unit = new_file_unit()
    call CreateTxtFile(fname,unit)

    do i=1, Names%nnames
        write(unit,*) trim(ParamNames_AsString(Names,i))
    end do   

    call CloseFile(unit)

    end subroutine ParamNames_WriteFile


    function ParamNames_NameOrNumber(Names,ix) result(name)
    Type(TParamNames) :: Names
    character(len=ParamNames_maxlen)  :: name
    integer, intent(in) :: ix

    name = ParamNames_name(Names,ix) 
    if (name == '') name = IntToStr(ix)

    end function ParamNames_NameOrNumber

    subroutine ParamNames_WriteMatlab(Names,  unit, headObj)
    Type(TParamNames) :: Names
    character(len=ParamNames_maxlen) name
    character(len=*), intent(in) :: headObj
    integer :: unit
    integer i

    do i=1, Names%nnames
        name = ParamNames_name(Names,i) 
        if (name /= '') then
            write(unit,'(a)', advance='NO') trim(headObj)//trim(name)//'= struct(''n'','''//trim(name) &
            //''',''i'','//trim(intToStr(i))//',''label'','''//trim(Names%label(i))//''',''isDerived'','
            if (Names%is_derived(i)) then
                write(unit,'(a)') 'true);'
            else
                write(unit,'(a)') 'false);'
            endif
        end if
    end do

    end subroutine ParamNames_WriteMatlab

    function ParamNames_ReadIniForParam(Names,Ini,Key, param) result(input)
    ! read Key[name] or Keyn where n is the parameter number
    use IniFile
    Type(TParamNames) :: Names
    Type(TIniFile) :: Ini
    character(LEN=*), intent(in) :: Key
    integer, intent(in) :: param
    character(LEN=128) input

    input = ''
    if (Names%nnames>0) then
        input = Ini_Read_String_File(Ini,trim(key)//'['//trim(Names%name(param))//']')
    end if
    if (input=='') then
        input = Ini_Read_String_File(Ini,trim(Key)//trim(IntToStr(param)))
    end if

    end function ParamNames_ReadIniForParam 

    function ParamNames_HasReadIniForParam(Names,Ini,Key, param) result(B)
    ! read Key[name] or Keyn where n is the parameter number
    use IniFile
    Type(TParamNames) :: Names
    Type(TIniFile) :: Ini
    character(LEN=*), intent(in) :: Key
    integer, intent(in) :: param
    logical B

    B = .false.
    if (Names%nnames>0) then
        B = Ini_HasKey_File(Ini,trim(key)//'['//trim(Names%name(param))//']')
    end if
    if (.not. B) then
        B = Ini_HasKey_File(Ini,trim(Key)//trim(IntToStr(param)))
    end if

    end function ParamNames_HasReadIniForParam 



    end module ParamNames
