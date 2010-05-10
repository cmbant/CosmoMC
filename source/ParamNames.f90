module ParamNames
 use AMlUtils
 implicit none
 
 integer, parameter :: ParamNames_maxlen = 128
 
  Type TParamNames
    integer nnames, num_MCMC, num_derived
    character(LEN=ParamNames_maxlen), dimension(:), pointer ::  name
    character(LEN=ParamNames_maxlen), dimension(:), pointer ::  label    
    logical, dimension(:), pointer ::  is_derived
  end Type TParamNames

contains

function IsWhiteSpace(C)
 character, intent(in) :: C
 logical IsWhiteSpace
 
 IsWhiteSpace = (C==' ') .or. (C==char(9)) 

end function IsWhiteSpace

subroutine ParamNames_Init(Names, filename)
 Type(TParamNames) :: Names
 character(Len=*), intent(in) :: filename
 integer handle,n, len, pos
 character (LEN=ParamNames_maxlen*3) :: InLine
  
  handle = new_file_unit()
  call OpenTxtFile(filename,handle)
  n = FileLines(handle)
  allocate(Names%name(n))
  allocate(Names%label(n)) 
  allocate(Names%is_derived(n))
  Names%is_derived = .false. 
  Names%num_MCMC = 0
  Names%num_derived = 0

  Names%name = '' 
  n=0
  do 
      read (handle,'(a)',end=400) InLine
      len = len_trim(InLIne)
      n=n+1
      pos =1
      do while (pos < len .and. IsWhiteSpace(InLIne(pos:pos))) 
       pos = pos+1
      end do 
      read(InLine(pos:), *, end=400, err=400) Names%name(n)
      pos = pos + len_trim(Names%name(n))
      do while (pos < len .and. IsWhiteSpace(InLIne(pos:pos))) 
       pos = pos+1
      end do 
      Names%label(n) = trim(adjustl(InLine(pos:ParamNames_maxlen))) 
      pos = scan(Names%label(n),'#')
      if (pos/=0) Names%label(n) = Names%label(n)(1:pos-1)
      pos = scan(Names%label(n),char(9))
      if (pos/=0) Names%label(n) = Names%label(n)(1:pos-1)      
      Names%name(n) = trim(adjustl(Names%name(n)))
      len = len_trim( Names%name(n) )
      if (Names%name(n)(len:len)=='*') then 
       Names%name(n)(len:len)=' ' 
       Names%is_derived(n) = .true.
       Names%num_derived = Names%num_derived + 1
      else
       Names%num_MCMC =  Names%num_MCMC + 1
      end if 
  end do

400  call CloseFile(handle) 
   Names%nnames = n

end subroutine ParamNames_Init

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
 integer i
 
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
     do while (pos < len .and. InLine(pos:pos)==' ') 
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

subroutine ParamNames_WriteFile(Names, fname)
 Type(TParamNames) :: Names
 character(LEN=*), intent(in) :: fname
 integer :: unit
 integer i
 character(LEN=ParamNames_maxlen) nm

 unit = new_file_unit()
 call CreateTxtFile(fname,unit)
   
   do i=1, Names%nnames
     nm = trim(Names%name(i))
     if (Names%is_derived(i)) nm = concat(nm,'*')
     write(unit,*) trim(nm)//char(9)//trim(Names%label(i))
   end do   

 call CloseFile(unit)
   
end subroutine ParamNames_WriteFile

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
