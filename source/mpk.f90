!Module storing observed matter power spectrum datasets, their points and window functions
!and routines for computing the likelihood

!This code is based on that in cmbdata.f90 
!and on Sam Leach's incorporation of Max Tegmark's SDSS code

!Originally SLB Sept 2004
!AL April 2006: added covariance matrix support (following 2df 2005)

module mpk
use settings
use cmbtypes
implicit none

 Type mpkdataset
    logical :: use_set
    integer :: num_mpk_points_use ! total number of points used (ie. max-min+1)
    integer :: num_mpk_kbands_use ! total number of kbands used (ie. max-min+1)
    character(LEN=20) :: name
    real, pointer, dimension(:,:) :: N_inv
    real, pointer, dimension(:,:) :: mpk_W, mpk_invcov
    real, pointer, dimension(:) :: mpk_P, mpk_sdev, mpk_k
   !for Q and A see e.g. astro-ph/0501174, astro-ph/0604335
    logical :: Q_marge
    real :: Q_mid, Q_sigma, Ag
   end Type mpkdataset

 integer :: num_mpk_datasets = 0
 Type(mpkdataset) mpkdatasets(10)

 !Note all units are in k/h here
 
  logical :: use_mpk = .false.
  
contains 

  subroutine ReadmpkDataset(gname)   
    character(LEN=*), intent(IN) :: gname
    character(LEN=120) :: kbands_file, measurements_file, windows_file, cov_file
    Type (mpkdataset) :: mset

    integer i,iopb
    real keff,klo,khi,beff
    integer :: num_mpk_points_full ! actual number of bandpowers in the infile
    integer :: num_mpk_kbands_full ! actual number of k positions " in the infile
    integer :: max_mpk_points_use ! in case you don't want the smallest scale modes (eg. sdss)
    integer :: min_mpk_points_use ! in case you don't want the largest scale modes
    integer :: max_mpk_kbands_use ! in case you don't want to calc P(k) on the smallest scales (will truncate P(k) to zero here!)
    integer :: min_mpk_kbands_use ! in case you don't want to calc P(k) on the largest scales (will truncate P(k) to zero here!)
    real, dimension(:,:), allocatable :: mpk_Wfull, mpk_covfull
    real, dimension(:), allocatable :: mpk_kfull, mpk_fiducial

    character(80) :: dummychar
    logical bad
 
    num_mpk_datasets = num_mpk_datasets + 1
    if (num_mpk_datasets > 10) stop 'too many datasets'

    call Ini_Open(gname, 1, bad, .false.)
    if (bad) then
      write (*,*)  'Error opening dataset file '//trim(gname)
      stop
    end if

    mset%name = Ini_Read_String('name') 
    mset%use_set =.true.
    if (Feedback > 0) write (*,*) 'reading: '//trim(mset%name)
    num_mpk_points_full = Ini_Read_Int('num_mpk_points_full',0)
    if (num_mpk_points_full.eq.0) write(*,*) ' ERROR: parameter num_mpk_points_full not set'
    num_mpk_kbands_full = Ini_Read_Int('num_mpk_kbands_full',0)
    if (num_mpk_kbands_full.eq.0) write(*,*) ' ERROR: parameter num_mpk_kbands_full not set'
    min_mpk_points_use = Ini_Read_Int('min_mpk_points_use',1)
    min_mpk_kbands_use = Ini_Read_Int('min_mpk_kbands_use',1)
    max_mpk_points_use = Ini_Read_Int('max_mpk_points_use',num_mpk_points_full)
    max_mpk_kbands_use = Ini_Read_Int('max_mpk_kbands_use',num_mpk_kbands_full)
    mset%num_mpk_points_use = max_mpk_points_use - min_mpk_points_use +1
    mset%num_mpk_kbands_use = max_mpk_kbands_use - min_mpk_kbands_use +1

    allocate(mpk_Wfull(num_mpk_points_full,num_mpk_kbands_full))
    allocate(mpk_kfull(num_mpk_kbands_full))
    allocate(mset%mpk_P(mset%num_mpk_points_use))
    allocate(mset%mpk_sdev(mset%num_mpk_points_use))  ! will need to replace with the covmat
    allocate(mset%mpk_k(mset%num_mpk_kbands_use))
    allocate(mset%mpk_W(mset%num_mpk_points_use,mset%num_mpk_kbands_use))
    allocate(mpk_fiducial(mset%num_mpk_points_use))

    kbands_file  = Ini_Read_String('kbands_file')
    call ReadVector(kbands_file,mpk_kfull,num_mpk_kbands_full)
    mset%mpk_k(1:mset%num_mpk_kbands_use)=mpk_kfull(min_mpk_kbands_use:max_mpk_kbands_use) 
    if (Feedback > 1) then 
       write(*,*) 'reading: ',mset%name,' data'
       write(*,*) 'Using kbands windows between',mset%mpk_k(1),' < k/h < ',mset%mpk_k(mset%num_mpk_kbands_use)      
    endif

    measurements_file  = Ini_Read_String('measurements_file')
    call OpenTxtFile(measurements_file, tmp_file_unit)
    mset%mpk_P=0.
    read (tmp_file_unit,*)dummychar
    read (tmp_file_unit,*)dummychar
    do i= 1, (min_mpk_points_use-1)
       read (tmp_file_unit,*, iostat=iopb) keff,klo,khi,beff,beff,beff
    end do
    if (Feedback > 1 .and. min_mpk_points_use>1) write(*,*) 'Not using bands with keff=  ',keff,' or below'
    do i =1, mset%num_mpk_points_use
       read (tmp_file_unit,*, iostat=iopb) keff,klo,khi,mset%mpk_P(i),mset%mpk_sdev(i),mpk_fiducial(i)
    end do
    close(tmp_file_unit) 
    if (Feedback > 1) write(*,*) 'bands truncated at keff=  ',keff
    
    windows_file  = Ini_Read_String('windows_file')
    if (windows_file.eq.'') write(*,*) 'ERROR: mpk windows_file not specified'
    call ReadMatrix(windows_file,mpk_Wfull,num_mpk_points_full,num_mpk_kbands_full)
    mset%mpk_W(1:mset%num_mpk_points_use,1:mset%num_mpk_kbands_use)= &
       mpk_Wfull(min_mpk_points_use:max_mpk_points_use,min_mpk_kbands_use:max_mpk_kbands_use)
    

    cov_file  = Ini_Read_String('cov_file')
    if (cov_file /= '') then
     allocate(mpk_covfull(num_mpk_points_full,num_mpk_points_full))
     call ReadMatrix(cov_file,mpk_covfull,num_mpk_points_full,num_mpk_points_full)
     allocate(mset%mpk_invcov(mset%num_mpk_points_use,mset%num_mpk_points_use))
     mset%mpk_invcov=  mpk_covfull(min_mpk_points_use:max_mpk_points_use,min_mpk_points_use:max_mpk_points_use)
     call Matrix_Inverse(mset%mpk_invcov)
     deallocate(mpk_covfull)
    else
     nullify(mset%mpk_invcov)
    end if

    mset%Q_marge = Ini_Read_Logical('Q_marge',.false.)
    if (mset%Q_marge) then
     mset%Q_mid = Ini_Read_Real('Q_mid')
     mset%Q_sigma = Ini_Read_Real('Q_sigma')
     mset%Ag = Ini_Read_Real('Ag', 1.4)
    end if 
    if (iopb.ne.0) then
       stop 'Error reading mpk file'
    endif
 
   call Ini_Close

   deallocate(mpk_Wfull, mpk_kfull,mpk_fiducial)

   mpkdatasets(num_mpk_datasets) = mset
 
  end subroutine ReadmpkDataset

 
  function LSS_mpklike(Theory,mset) result(LnLike)
   Type (mpkdataset) :: mset
   Type (CosmoTheory) Theory
   real LnLike
   real, dimension(:), allocatable :: mpk_Pth, mpk_lin
   real, dimension(:), allocatable :: w
   real, dimension(:), allocatable :: mpk_WPth
   real :: covdat(mset%num_mpk_points_use), covth(mset%num_mpk_points_use)
   real :: normV, Ag, Q, minchisq
   real tmp
   integer :: i, iQ
   logical :: do_marge
   integer, parameter :: nQ=6
   real :: dQ = 0.4
   real chisq(-nQ:nQ)
   real calweights(-nQ:nQ)
 
   allocate(mpk_lin(mset%num_mpk_kbands_use) ,mpk_Pth(mset%num_mpk_kbands_use))
   allocate(w(mset%num_mpk_points_use))
   allocate(mpk_WPth(mset%num_mpk_points_use))
   chisq = 0

   if (.not. mset%use_set) then
      LnLike = 0
      return
   end if

   ! won't actually want to do this multiple times for multiple galaxy pk data sets?..
   do i=1, mset%num_mpk_kbands_use
      mpk_lin(i)=MatterPowerAt(Theory,mset%mpk_k(i))
   end do


  do_marge = mset%Q_Marge
  if (mset%Q_sigma==0) do_marge = .false.

  do iQ=-nQ,nQ
    Q = mset%Q_mid +iQ*mset%Q_sigma*dQ 
 
   if (mset%Q_marge) then
      mpk_Pth=mpk_lin*(1+Q*mset%mpk_k**2)/(1+mset%Ag*mset%mpk_k)
   else 
      mpk_Pth = mpk_lin
   end if

   ! SLB stuff
!   if (mset%mpk_sdev(1).eq.mset%mpk_P(1)) then ! the data is dummy data, so output Pk for simulating data
!      write(*,*) ' P(k) points and errors are the same. '
!      write(*,*) ' So just writing out theoretical P(k) and exiting'
!      do i=1, mset%num_mpk_kbands_use
!	     write(10,*) mset%mpk_k(i), mpk_Pth(i)      ! for simulating data from
!      end do
!      write(*,*) ' Wrote ',mset%num_mpk_kbands_use,' theory P(k) values to fort.10'
!      stop
!   end if
  
    mpk_WPth = matmul(mset%mpk_W,mpk_Pth)

   !with analytic marginalization over normalization nuisance (flat prior on b^2)
   !See appendix F of cosmomc paper
   
    if (associated(mset%mpk_invcov)) then
     covdat = matmul(mset%mpk_invcov,mset%mpk_P)
     covth = matmul(mset%mpk_invcov,mpk_WPth)
     normV = sum(mpk_WPth*covth)
     chisq(iQ) = sum(mset%mpk_P*covdat)  - sum(mpk_WPth*covdat)**2/normV  + log(normV)  

    else
     !with analytic marginalization over normalization nuisance (flat prior on b^2)
     w=1/(mset%mpk_sdev**2)
     normV = sum(mpk_WPth*mpk_WPth*w)
     tmp=sum(mpk_WPth*mset%mpk_P*w)/normV ! avoid subtracting one large number from another
     chisq(iQ) = sum(mset%mpk_P*(mset%mpk_P - mpk_WPth*tmp)*w)  + log(normV)
    end if

   if (do_marge) then
     calweights(iQ) = exp(-(iQ*dQ)**2/2)
   else 
     LnLike = chisq(iQ)/2
    exit
   end if
   
   end do


   !without analytic marginalization
 !! chisq = sum((mset%mpk_P(:) - mpk_WPth(:))**2*w) ! uncommented for debugging purposes

   if (do_marge) then
    minchisq=minval(chisq)
    LnLike = sum(exp(-(chisq-minchisq)/2)*calweights)/sum(calweights)
     if (LnLike == 0) then
     LnLike = LogZero
    else
     LnLike =  -log(LnLike) + minchisq/2
    end if
   end if
   if (Feedback>1) write(*,*) 'mpk chi-sq:', LnLike*2

   if (LnLike > 1e8) then
       write(*,*) 'Chisq is huge, maybe there is a problem? chisq=',chisq
   end if

   deallocate(mpk_Pth,mpk_lin)
   deallocate(mpk_WPth)
   
 end function LSS_mpklike


 function LSSLnLike(CMB, Theory)
   Type (CMBParams) CMB
   Type (CosmoTheory) Theory
   real LSSLnLike
   integer i
   real tot(num_mpk_datasets)

  do i=1, num_mpk_datasets
     if (mpkdatasets(i)%name == 'twodf') then
        stop 'twodf no longer supported - use data/2df_2005.dataset'
     else
      tot(i) = LSS_mpklike(Theory,mpkdatasets(i))
     end if
  end do
  LSSLnLike = SUM(tot) 
  
 end function LSSLnLike


end module 


