
module IO
!Module to wrap output functions
!Replace file e.g. to store output in a pipline database (indexed by filename)
!Al May09
!Valid handles must be /=0
use AmlUtils
implicit none

 

contains


   subroutine IO_Ini_Load(Ini,InputFile, bad)
     use IniFile
     Type(TIniFile) :: Ini
     character(LEN=*), intent(in) :: InputFile
     logical bad
     integer file_unit
     
     file_unit = new_file_unit()
     call Ini_Open(InputFile, file_unit, bad, .false.)
     call ClearFileUnit(file_unit)
  end subroutine IO_Ini_Load

  function IO_OpenChainForRead(name) result(handle)
   character(len=*), intent(in) :: name
   integer handle

    handle = new_file_unit()
    call OpenTxtFile(name,handle)
  
  end function IO_OpenChainForRead

  function IO_OpenDataForRead(name) result(handle)
   character(len=*), intent(in) :: name
   integer handle

   !assumes handle is a fortran unit; If want to get from database then probably to file here
   
    handle = new_file_unit()
    call OpenFile(name,handle,'unformatted')
  
  end function IO_OpenDataForRead

  function IO_OutputOpenForWrite(name, append, isLogFile) result(handle)
   character(len=*), intent(in) :: name
   logical, intent(in), optional :: append, isLogFile
   logical app
   integer handle
   
   if (present(append)) then
     app = append
   else
     app = .false.
   end if
    
   handle = new_file_unit()
   call CreateOpenTxtFile(name,handle,app)
  
  end function IO_OutputOpenForWrite
        
   function IO_DataOpenForWrite(name, append) result(handle)
   !e.g. cached C_l data
   character(len=*), intent(in) :: name
   logical, intent(in), optional :: append
   logical app
   integer handle
   !assumes handle is a fortran unit; If want to put in database then probably convert file when handle closed
   
   if (present(append)) then
     app = append
   else
     app = .false.
   end if
    
   handle = new_file_unit()
   call CreateOpenFile(name,handle,'unformatted',app)
  
  end function IO_DataOpenForWrite      
        
        
  subroutine IO_Close(handle, IsLogFile)
   integer, intent(in) :: handle
   logical, intent(in), optional :: isLogFile
   
    call CloseFile(handle)
   
  end subroutine IO_Close      
  
  subroutine IO_DataCloseWrite(handle)
   integer, intent(in) :: handle
    !modify e.g. to grab data from temporary file 
   
    call CloseFile(handle)
   
  end subroutine IO_DataCloseWrite      
  
  function IO_Exists(name) result(res)
   logical res
   character(LEN=*), intent(in)::name
   
   res = FileExists(name)
   
  end function IO_Exists

 subroutine IO_OutputChainRow(handle, mult, like, values, nvalues)
  integer, intent(in) :: handle
  real mult, like, values(:)
  integer, intent(in), optional :: nvalues
  character(LEN =128) fmt
  integer n
  
  if (present(nvalues)) then
   n = nvalues
  else
   n = size(values)
  end if 

   fmt = trim(numcat('(2E16.7,',n))//'E16.7)'
   write (handle,fmt) mult, like, values(1:n)

   if (flush_write) call FlushFile(handle)

 end subroutine IO_OutputChainRow
 
 subroutine IO_WriteLog(handle, S)
   integer, intent(in) :: handle
   character(LEN=*), intent(in) :: S
   
   write(handle,*) trim(S)
   if (flush_write) call FlushFile(handle) 
 
 end subroutine IO_WriteLog

 function IO_ReadChainRow(handle, mult, like, values, nvalues) result(OK)
  logical OK
  integer, intent(in) :: handle
  real, intent(out) :: mult, like, values(:)
  integer, intent(in), optional :: nvalues
  integer n
  
    if (present(nvalues)) then
    n = nvalues
    else
    n = size(values)
    end if 

    read(handle, *, end=100,err=100) mult, like, values(1:n)    
    OK = .true.
    return   
100 OK = .false. 

 end function IO_ReadChainRow

 subroutine IO_ReadLastChainParams(name, mult, like, values, nvalues)
  character(LEN=*), intent(in) :: name
  real, intent(out) :: mult, like, values(:)
  integer, intent(in), optional :: nvalues
  character(LEN=5000) :: InLine
  integer n
  
  if (present(nvalues)) then
   n = nvalues
  else
   n = size(values)
  end if 
       
  InLine = LastFileLine(name)
  read(InLine, *) mult, like, values(1:n)

 end subroutine IO_ReadLastChainParams

 subroutine IO_OutputParamNames(Names, fname)
   use ParamNames
   Type(TParamNames) :: Names
   character(len=*), intent(in) :: fname
      
   call ParamNames_WriteFile(Names,fname)

 end subroutine IO_OutputParamNames
 

end module IO