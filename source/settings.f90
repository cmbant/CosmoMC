module settings
  use AMLutils
  use Random
  use IniFile
  implicit none

  real :: AccuracyLevel = 1.
  !Set to >1 to use CAMB etc on higher accuracy settings. 
  !Does not affect MCMC (except making it all slower)

!num_norm is number of calibration and normalization parameters (fast)
!num_initpower is the number of initial power parameters
!num_hard is the number of 'hard' parameters

  integer, parameter :: num_hard =7 
  integer, parameter :: num_initpower = 3 
  integer, parameter :: num_norm = 3 
     !Should be 2 or larger (scale for the scalar Cl, and the T/S ratio

  logical :: use_fast_slow = .true.
    !Set to false if using a slow likelihood function so no there's point is treating
    !'fast' parameters differently (in fact, doing so will make performance worse)

  integer :: oversample_fast = 0
  integer, parameter :: sampling_metropolis = 1, sampling_slice = 2, sampling_fastslice =3, &
         sampling_slowgrid = 4,  sampling_multicanonical = 5,  sampling_wang_landau = 6

  integer :: sampling_method = sampling_metropolis
 
  !scale of the proposal distribution in units of the posterior standard deviation
  real    :: propose_scale  = 2.4 

!The rest are set up automatically

  real, parameter :: cl_norm = 1e-10 !units for As

  integer, parameter :: num_fast_params = num_initpower + num_norm

  integer, parameter :: num_params = num_norm + num_initpower + num_hard
  integer, parameter :: index_initpower = num_hard+1
  integer, parameter :: index_norm = index_initpower + num_initpower
  integer, dimension(:), allocatable :: params_used,fast_params_used
  integer num_params_used, num_fast, num_slow

  integer :: num_threads = 0
  integer :: instance = 0
  integer :: MPIchains = 1, MPIrank = 0

  logical :: Use_LSS = .true.

  integer :: logfile_unit  = 0
  integer :: outfile_unit = 0 
  integer :: indepfile_unit = 0
  integer :: slow_proposals = 0
  integer :: output_lines = 0

  real, parameter :: logZero = 1e30
  character (LEN =120) :: FileChangeIni = '', FileChangeIniAll = ''


contains


 subroutine CheckParamChangeF(F)
   character(LEN=*), intent(in) ::  F
   logical bad, doexit

   if (F /= '') then

       call Ini_Open(F, tmp_file_unit, bad, .false.)
       if (bad) return
       Ini_fail_on_not_found = .false.
       doexit = (Ini_Read_Int('exit',0) == 1) 
       FeedBack = Ini_Read_Int('feedback',Feedback)
       num_threads = Ini_Read_Int('num_threads',num_threads)
       call Ini_Close
       if (F== FileChangeIni) call DeleteFile(FileChangeini)
       if (doexit) stop

   end if

 end subroutine CheckParamChangeF

 subroutine CheckParamChange

   call CheckParamChangeF(FileChangeIni)
   if (FileChangeIni/=FileChangeIniAll) call CheckParamChangeF(FileChangeIniAll)

 end subroutine CheckParamChange


 subroutine ReadMatrix(aname, mat, m,n)
   character(LEN=*), intent(IN) :: aname
   integer, intent(in) :: m,n
   real, intent(out) :: mat(m,n)
   integer j,k
   real tmp

   if (Feedback > 0) write(*,*) 'reading: '//trim(aname)
   call OpenTxtFile(aname, tmp_file_unit)

   do j=1,m
      read (tmp_file_unit,*, end = 200, err=100) mat(j,1:n)
   end do
   goto 120

100 rewind(tmp_file_unit)  !Try other possible format
   do j=1,m 
    do k=1,n
      read (tmp_file_unit,*, end = 200) mat(j,k)
    end do
   end do

120 read (tmp_file_unit,*, err = 150, end =150) tmp
   goto 200

150 close(tmp_file_unit)
    return

 200 write (*,*) 'matrix file '//trim(aname)//' is the wrong size'
     stop

 end subroutine ReadMatrix

  subroutine ReadVector(aname, vec, n)
   character(LEN=*), intent(IN) :: aname
   integer, intent(in) :: n
   real, intent(out) :: vec(n)
   integer j

   if (Feedback > 0) write(*,*) 'reading: '//trim(aname)

   call OpenTxtFile(aname, tmp_file_unit)

   do j=1,n
      read (tmp_file_unit,*, end = 200) vec(j)
   end do


    close(tmp_file_unit)
    return

 200 write (*,*) 'vector file '//trim(aname)//' is the wrong size'
     stop

 end subroutine ReadVector

 subroutine WriteVector(aname, vec, n)
   character(LEN=*), intent(IN) :: aname
   integer, intent(in) :: n
   real, intent(in) :: vec(n)
   integer j
  
   call CreateTxtFile(aname, tmp_file_unit)

   do j=1,n
      write (tmp_file_unit,'(1E15.6)') vec(j)
   end do

   close(tmp_file_unit)
  
 end subroutine WriteVector

  subroutine Matrix_Write(aname, mat, forcetable)
   character(LEN=*), intent(in) :: aname
   real, intent(in) :: mat(:,:)
   logical, intent(in), optional :: forcetable
   integer i,k
   character(LEN=50) fmt
   integer shp(2)
   logical WriteTab

   shp = shape(mat)
   WriteTab = shp(2)<=50
   if (present(forcetable)) then
     if (forcetable) WriteTab = .true.
    end if

   call CreateTxtFile(aname, tmp_file_unit)
   fmt = trim(numcat('(',shp(2)))//'E15.6)'
   do i=1, shp(1)
     if (.not. WriteTab) then
      do k=1, shp(2)
       write (tmp_file_unit, '(1E15.6)') mat(i,k)
      end do
     else
      write (tmp_file_unit, fmt) mat(i,1:shp(2))
     end if
   end do

   close(tmp_file_unit)

 end subroutine Matrix_Write


 subroutine WriteSqMatrix(aname, mat,n)
   character(LEN=*), intent(in) :: aname
   integer, intent(in) :: n
   real, intent(in) :: mat(n,n)
   integer i
   character(LEN=50) fmt

   call CreateTxtFile(aname, tmp_file_unit)
   fmt = trim(numcat('(',n))//'E15.6)'
   do i=1, n
      write (tmp_file_unit, fmt) mat(i,1:n)
   end do

   close(tmp_file_unit)

 end subroutine WriteSqMatrix

       subroutine Diagonalize(m, diag, n)
          !Does m = U diag U^T, returning U in m
          integer, intent(in) :: n
          real m(n,n), diag(n)
          integer ierr, tmpsize
          real tmp(3*n**2)

          tmpsize = 3*n**2
          call SSYEV('V','U',n,m,n,diag,tmp,tmpsize,ierr) !evalues and vectors of symmetric matrix
        
       end subroutine Diagonalize


  subroutine Matrix_Inverse(M)
   !This should not be used in real situations, but useful for quick testing
     real, intent(inout):: M(:,:)
     real w(Size(M,DIM=1)),tmp(Size(M,DIM=1),Size(M,DIM=1))
     integer i, n
     real sigmas(Size(M,DIM=1))

     n=Size(M,DIM=1)
     if (n.eq.1) then
       M=1/M
       return
     end if
     if (n<0) then
       write(*,*) 'ERROR: n<0 in Matrix_Inverse!!'
       return
     end if

     if (Size(M,DIM=2)/=n) stop 'Matrix_Inverse: non-square matrix'

     do i=1, n
       if (M(i,i) < 1e-15) stop 'Matrix_Inverse: small or negative diagonal element'
       sigmas(i) = sqrt(M(i,i))
       M(i,:) = M(i,:) /sigmas(i)
       M(:,i) = M(:,i) /sigmas(i)
     end do
     call Diagonalize(M,w,n)      
     do i=1, n
        tmp(i,:) = M(:,i)/w(i)
     end do
     M = matmul(M,tmp)     
     do i=1, n
       M(i,:) = M(i,:) /sigmas(i)
       M(:,i) = M(:,i) /sigmas(i)
     end do

  end subroutine Matrix_Inverse



end module settings

