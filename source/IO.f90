
module IO
!Module to wrap output functions
!Replace file e.g. to store output in a pipline database (indexed by filename)
!Valid handles must be /=0
use AmlUtils
use settings
use MatrixUtils
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

  function IO_OpenChainForRead(name, OK) result(handle)
   character(len=*), intent(in) :: name
   integer handle
   logical, optional, intent(out) :: OK

    handle = new_file_unit()
    if (.not. present(OK)) then  
     call OpenTxtFile(name,handle)
    else
    
     open(unit=handle,file=name,form='formatted',status='old', err=28)
     OK=.true.
     return
28   OK=.false.
   
    endif     
    
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

 subroutine IO_WriteProposeMatrix(pmat, prop_mat, comment)
   real pmat(:,:)
   character(LEN=*), intent(in) :: prop_mat
   character(LEN=*), optional, intent(in) :: comment
   
   if (present(comment)) then
      call Matrix_write(prop_mat,pmat,.true.,commentline=comment)
   else
      call Matrix_write(prop_mat,pmat,.true.)
   endif 
   
 end subroutine IO_WriteProposeMatrix

 subroutine IO_ReadProposeMatrix(pmat, prop_mat)
   use ParamNames
   real pmat(:,:)
   character(LEN=1024), intent(in) :: prop_mat
   real,allocatable :: tmpMat(:,:)
   integer i,y
   integer file_id 
   character(LEN=4096) :: InLine
   integer num, cov_params(256) 
    
    file_id = new_file_unit()
    call OpenTxtFile(prop_mat, file_id)
    InLine=''
    do while (InLine == '') 
     read(file_id,'(a)', end = 10) InLine
    end do
    InLine = adjustl(InLine)
    If (InLine(1:1)=='#') then
     !Have paramnames to identify
      InLine = InLine(2:len_trim(InLine))
      num=-1
      call ParamNames_ReadIndices(NameMapping,InLine, cov_params, num)
      pmat=0
      y=0
      do
       read(file_id,'(a)', end = 20) InLine
       if (InLine/='') then
         y=y+1
         read(InLine,*,end=20) pmat(cov_params(1:num),cov_params(y))
         if (y==num) exit
       end if
      end do 
      call CloseFile(file_id) 
      return
20   call mpiStop('ReadProposeMatrix: wrong number of rows/columns in .covmat')
        
    end if
10  call CloseFile(file_id)  
 
   i=TxtNumberColumns(InLine)
   if (i==num_params) then
    call ReadMatrix(prop_mat,pmat, num_params, num_params)
   else if (i==num_real_params) then
    allocate(tmpMat(num_real_params,num_real_params))
    call ReadMatrix(prop_mat,tmpmat, num_real_params, num_real_params)
    pmat=0
    pmat(1:num_real_params,1:num_real_params) = tmpMat
    deallocate(tmpMat)
   else
    call MpiStop('Propose matrix the wrong size: '//trim(prop_mat))
   end if
   
end subroutine IO_ReadProposeMatrix


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

 function IO_SkipChainRows(handle,nrows) result(OK)
   integer, intent(in):: handle,nrows
   logical OK
   integer ix
   character(LEN=10000) InLine  
   
      do ix = 1, nrows
        read (handle,'(a)',end=1) InLine
      end do
      OK = .true.
      return
 1    OK = .false.      
      
 end function IO_SkipChainRows

 function IO_ReadChainRow(handle, mult, like, values, nvalues, chainOK, samples_chains) result(OK)
  !Returns OK=false if end of file or if not enough values on each line, otherwise OK = true
  !Returns chainOK = false if bad line or NaN, chainOK=false and OK=true for NaN (continue reading)
  logical OK
  integer, intent(in) :: handle
  real, intent(out) :: mult, like, values(:)
  integer, intent(in), optional :: nvalues
  logical, optional, intent(out) :: ChainOK
  logical, optional, intent(in) :: samples_chains
  logical samples_are_chains
  integer n
   character(LEN=10000) InLine  
  
    if (present(nvalues)) then
    n = nvalues
    else
    n = size(values)
    end if 
   
    if (present(samples_chains)) then
     samples_are_chains=samples_chains
     else
     samples_are_chains = .true.
    endif

    if (present(ChainOK)) chainOK = .true.

    read (handle,'(a)',end=100) InLine
    if (SCAN (InLine, 'N') /=0) then
       OK = .true.
       if (present(ChainOK)) chainOK = .false.
       return
    end if
    
    if (samples_are_chains) then
     read(InLine, *, end=100,err=110) mult, like, values(1:n)    
    else
     mult=1
     like=1
     read(InLine, *, end=100,err=110) values(1:n)    
    end if
    OK = .true.
    return   
100 OK = .false. 
    return
110 OK = .false.
    if (present(ChainOK)) chainOK = .false.

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
      
   call ParamNames_WriteFile(Names,trim(fname)//'.paramnames')

 end subroutine IO_OutputParamNames

 subroutine IO_ReadParamNames(Names, in_root)
   use ParamNames
   Type(TParamNames) :: Names
   character(LEN=*), intent(in) :: in_root
   character(LEN=Ini_max_string_len) infile
      
        infile = trim(in_root) // '.paramnames'
        if (FileExists(infile)) then
             call ParamNames_Init(Names,infile)
        end if         

 end subroutine IO_ReadParamNames


 
 function IO_ReadChainRows(in_root, chain_ix,chain_num, ignorerows, nrows, &
        ncols,max_rows,coldata,samples_are_chains) result(OK)
  !OK = false if chain not found or not enough samples
  character(LEN=*), intent(in) :: in_root
  integer,intent(in) :: chain_ix, chain_num
  integer, intent(in) :: max_rows, ignorerows
  integer, intent(in) :: ncols
  real(KIND(1.d0)), intent(inout) :: coldata(ncols,0:max_rows) !(col_index, row_index)
  logical, intent(in) :: samples_are_chains
  integer, intent(inout) :: nrows
  logical OK
  real invars(1:ncols)
  logical chainOK
  integer chain_handle, row_start
  character(LEN=Ini_max_string_len) infile, numstr

       row_start=nrows
       if (chain_num == 0) then
            infile = trim(in_root) // '.txt'
         else
            write (numstr,*) chain_ix
            infile = trim(in_root) //'_'//trim(adjustl(numstr))// '.txt'
         end if

         write (*,*) 'reading ' // trim(infile)

         chain_handle = IO_OpenChainForRead(infile, chainOK)
         if (.not. chainOK) then
          write (*,'(" chain ",1I4," missing")') chain_ix
          OK = .false.
          return
         end if       
 
        if (ignorerows >=1) then
           if (.not. IO_SkipChainRows(chain_handle,ignorerows)) then
            call IO_Close(chain_handle) 
            OK = .false.
            return
           end if
          end if
          
          OK = .true.
          do
            if (.not. IO_ReadChainRow(chain_handle, invars(1), invars(2), &
                    invars(3:),ncols-2,chainOK,samples_are_chains)) then
             if (.not. chainOK) then
              write (*,*) 'error reading line ', nrows -row_start + ignorerows ,' - skipping rest of file'
             endif 
             call IO_Close(chain_handle)  
             return 
            else
             if (.not. chainOK) then
               write (*,*) 'WARNING: skipping line with probable NaN'
               cycle 
             end if 
            end if 
   
           coldata(1:ncols, nrows) = invars(1:ncols)
           nrows = nrows + 1
           if (nrows > max_rows) stop 'need to increase max_rows'
         
          end do 

 end function IO_ReadChainRows

 subroutine IO_OutputMargeStats(Names, froot,num_vars,num_contours, contours,contours_str, &
           cont_lines, colix, mean, sddev, has_limits_bot, has_limits_top, labels, force_twotail)
        use ParamNames
       Type(TParamNames) :: Names
       character(LEN=*), intent(in) :: froot
       integer, intent(in) :: num_vars, num_contours
       logical,intent(in) :: force_twotail, has_limits_bot(*),has_limits_top(*)
       real, intent(in) :: mean(*), sddev(*), contours(*), cont_lines(:,:,:)
       character(LEN=*), intent(in) :: contours_str
       integer,intent(in) :: colix(*)
       character(LEN=128) labels(*), tag
        
       integer i,j,file_id
         
         file_id = new_file_unit()
         open(unit=file_id,file=trim(froot)//'.margestats',form='formatted',status='replace')
          write(file_id,'(a)',advance='NO') 'param  mean           sddev          '       
          do j=1, num_contours
           write(file_id,'(a)',advance='NO') trim(concat('lower',j))//'         '//trim(concat('upper',j))//'         '
          end do
          write(file_id,'(a)') ''
         
          do j=1, num_vars
             write(file_id,'(1I5,2E15.7)', advance='NO') colix(j)-2, mean(j), sddev(j)
             do i=1, num_contours
               write(file_id,'(2E15.7)',advance='NO') cont_lines(j,1:2,i)
             end do
             write(file_id,'(a)') '   '//trim(labels(colix(j)))
          end do
          write (file_id,*) ''
 
          write (file_id,'(a)') 'Limits are: ' // trim(contours_str)
          if (.not. force_twotail) then
              do j=1, num_vars
              if (has_limits_bot(colix(j)).and. has_limits_top(colix(j))) then
                  tag='no tail'
                elseif (has_limits_bot(colix(j))) then
                 tag= '> one tail' 
                elseif (has_limits_top(colix(j))) then 
                 tag = '< one tail' 
                else
                 tag = 'two tail'
                end if
               write(file_id,'(1I5," ",1A20," ",1A12)', advance='NO')  &
                    colix(j)-2, ParamNames_name(Names,colix(j)-2),tag
               write(file_id,'(a)') trim(labels(colix(j)))

              end do
          else
              write (file_id,*) 'All limits are two tail'
          end if

         call CloseFile(file_id)
          
 end subroutine IO_OutputMargeStats

end module IO