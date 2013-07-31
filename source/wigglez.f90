!Module storing observed matter power spectrum datasets, their points and window functions
!and routines for computing the likelihood

!This code is based on that in cmbdata.f90 
!and on Sam Leach's incorporation of Max Tegmark's SDSS code
!
!Originally SLB Sept 2004
!AL April 2006: added covariance matrix support (following 2df 2005)
!LV_06 : incorporation of LRG DR4 from Tegmark et al . astroph/0608632 
!AL: modified LV SDSS to do Q and b^2 or b^2*Q marge internally as for 2df
!BR09: added model LRG power spectrum.
!AL Oct 20: switch to Ini_Read_xxx_File; fortran compatibility changes

!WiggleZ Matter power spectrum likelihood module.  Format is based upon mpk.f90
!DP & JD 2013 For compatibility with the latest version of CosmoMC (March2013)


module wigglezinfo
!David Parkinson 12th March 2012
  use settings
  use cmbtypes
  use Precision
!use CMB_Cls

  implicit none

  integer, dimension(4) :: izwigglez
  real(dl), parameter :: z0 = 0.d0, za = 0.22d0, zb = 0.41d0, zc = 0.6d0, zd = 0.78d0

!in CAMB: 5=now (z=0), 4=a, 3=b, 2=c, 1=d; opposite order in matter_power
!! now generalized indices iz0, iza, izb, izc, izd
  real(dl), dimension(4) :: zeval, zweight, sigma2BAOfid, sigma2BAO

  real(dl) om_val, ol_val, ok_check, wval  ! passed in from CMBparams CMB

 ! power spectra evaluated at GiggleZ fiducial cosmological theory 
 real, dimension(num_matter_power,4) :: power_hf_fid
 !make allocatable to avoid compile-time range errors when matter_power_lnzsteps<4

contains

  subroutine GiggleZinfo_init(redshift)
    integer :: iopb, i, ios, iz
    real(dl) :: kval, power_nl
    real(mcp) redshift
    logical save

    zeval(1) = za
    zeval(2) = zb
    zeval(3) = zc
    zeval(4) = zd
    
    iz = 0
    do i=1,4
       if(abs(redshift-zeval(i)).le.0.001) iz = i
    enddo

    !! first read in everything needed from the CAMB output files.
    iopb = 0 !! check later if there was an error

    if(iz.eq.1) then
       open(unit=tmp_file_unit,file=trim(DataDir)//'gigglezfiducialmodel_matterpower_a.dat',form='formatted',err=500, iostat=ios)
    else if(iz.eq.2) then
       open(unit=tmp_file_unit,file=trim(DataDir)//'gigglezfiducialmodel_matterpower_b.dat',form='formatted',err=500, iostat=ios)
    else if(iz.eq.3) then
       open(unit=tmp_file_unit,file=trim(DataDir)//'gigglezfiducialmodel_matterpower_c.dat',form='formatted',err=500, iostat=ios)
    else if(iz.eq.4) then
       open(unit=tmp_file_unit,file=trim(DataDir)//'gigglezfiducialmodel_matterpower_d.dat',form='formatted',err=500, iostat=ios)
    else 
       call MpiStop('could not indentify redshift')       
    endif

    do i = 1, num_matter_power
       read (tmp_file_unit,*,iostat=iopb) kval, power_nl
       power_hf_fid(i,iz) = power_nl
    end do
    close(tmp_file_unit)


500 if(ios .ne. 0) stop 'Unable to open file'
  if(iopb .ne. 0) stop 'Error reading model or fiducial theory files.'    
    
  end subroutine GiggleZinfo_init

! HARD CODING OF POLYNOMIAL FITS TO FOUR REDSHIFT BINS.
  subroutine GiggleZtoICsmooth(k,fidpolys)
    real(mcp), intent(in) :: k
    real(mcp) :: fidz_0, fidz_1, fidz_2, fidz_3, fidz_4
    real(mcp), dimension(4), intent(out) :: fidpolys


    fidz_1 = (4.619d0 - 13.7787d0*k + 58.941d0*k**2 - 175.24d0*k**3 + 284.321d0*k**4 - 187.284d0*k**5)
    fidz_2 = (4.63079d0 - 12.6293d0*k + 42.9265d0*k**2 - 91.8068d0*k**3 + 97.808d0*k**4 - 37.633d0*k**5)
    fidz_3 = (4.69659d0 - 12.7287d0*k + 42.5681d0*k**2 - 89.5578d0*k**3 + 96.664d0*k**4 - 41.2564*k**5)
    fidz_4 = (4.6849d0 - 13.4747d0*k + 53.7172d0*k**2 - 145.832d0*k**3 + 216.638d0*k**4 - 132.782*k**5)
    
    
    fidpolys(1) = 10**fidz_1
    fidpolys(2) = 10**fidz_2
    fidpolys(3) = 10**fidz_3
    fidpolys(4) = 10**fidz_4
    return

  end subroutine GiggleZtoICsmooth

  subroutine fill_GiggleZTheory(Theory, minkh, dlnkh,z)
    Type(TheoryPredictions) Theory
    real(mcp), intent(in) :: minkh, dlnkh
    real(mcp), intent(in) :: z
    real(mcp) :: logmink, xi, kval, expval,  nlrat
    real(mcp), dimension(4) :: fidpolys
    real(mcp) y, dz, matter_power_dlnz
    real(mcp) pk_hf, holdval
    integer :: i,iz, iz2, ik
    character(len=32) fname
    logmink = log(minkh)

    iz = 0
    do i=1,4
       if(abs(z-zeval(i)).le.0.001) iz = i
    enddo

    do ik=1,num_matter_power
       xi = logmink + dlnkh*(ik-1)
       kval = exp(xi)
       Theory%WiggleZPk(ik) = 0.
       pk_hf = Theory%matter_power(ik,izwigglez(iz))
   
       call GiggleZtoICsmooth(kval,fidpolys)

       holdval = pk_hf*fidpolys(iz)/power_hf_fid(ik,iz)
       
       Theory%WiggleZPk(ik) = Theory%WiggleZPk(ik) + holdval     

    end do

  end subroutine fill_GiggleZTheory


end module wigglezinfo


module wigglez
use precision
use settings
use cmbtypes
use likelihood
use wigglezinfo
implicit none

type, extends(CosmologyLikelihood) :: WiggleZLikelihood
  logical :: use_set
  integer :: num_mpk_points_use ! total number of points used (ie. max-min+1)
  integer :: num_mpk_kbands_use ! total number of kbands used (ie. max-min+1)
  integer :: num_regions_used   ! total number of wigglez regions being used
! 1st index always refers to the region
! so mpk_P(1,:) is the Power spectrum in the first active region
  real(mcp), pointer, dimension(:,:,:) :: mpk_W, mpk_invcov
  real(mcp), pointer, dimension(:,:) :: mpk_P
  real(mcp), pointer, dimension(:) :: mpk_k
  logical :: use_scaling !as SDSS_lrgDR3
  logical, pointer, dimension(:) :: regions_active
  logical use_gigglez
  !for Q and A see e.g. astro-ph/0501174, astro-ph/0604335
  logical :: Q_marge, Q_flat
  real(mcp):: Q_mid, Q_sigma, Ag
  real(mcp) :: redshift ! important to know
  contains
   procedure :: LogLike => WiggleZ_Lnlike
 end type WiggleZLikelihood

 integer, parameter :: max_num_wigglez_regions = 7
 !Note all units are in k/h here
 integer, parameter :: mpk_d = kind(1.d0)

contains 

  subroutine WiggleZLikelihood_Add(LikeList, Ini)
    use IniFile
    use settings
    class(LikelihoodList) :: LikeList
    Type(TIniFile) :: ini
    Type(WiggleZLikelihood), pointer :: like
    integer nummpksets, i
    
    
    use_wigglez_mpk = (Ini_Read_Logical_File(Ini, 'use_wigglez_mpk',.false.))
    
    if(.not. use_wigglez_mpk) return
    
    nummpksets = Ini_Read_Int('mpk_wigglez_numdatasets',0)
    do i= 1, nummpksets
      allocate(like)
      call ReadWiggleZDataset(like, ReadIniFileName(Ini,numcat('wigglez_dataset',i)) )
      like%LikelihoodType = 'MPK'
      like%needs_powerspectra = .true.
      call LikeList%Add(like)
    end do
    if (Feedback>1) write(*,*) 'read WiggleZ MPK datasets'
      
  end subroutine WiggleZLikelihood_Add
  
  subroutine ReadWiggleZDataset(wmset,gname)   
! this will be called once for each redshift bin
    use wigglezinfo
    use MatrixUtils
    implicit none
    type(WiggleZLikelihood) wmset
    character(LEN=*), intent(IN) :: gname
    character(LEN=Ini_max_string_len) :: kbands_file, measurements_file, windows_file, cov_file


    integer i,iopb,i_regions
    real(mcp) keff,klo,khi,beff
    integer :: num_mpk_points_full ! actual number of bandpowers in the infile
    integer :: num_mpk_kbands_full ! actual number of k positions " in the infile
    integer :: max_mpk_points_use ! in case you don't want the smallest scale modes (eg. sdss)
    integer :: min_mpk_points_use ! in case you don't want the largest scale modes
    integer :: max_mpk_kbands_use ! in case you don't want to calc P(k) on the smallest scales (will truncate P(k) to zero here!)
    integer :: min_mpk_kbands_use ! in case you don't want to calc P(k) on the largest scales (will truncate P(k) to zero here!)
    real(mcp), dimension(:,:,:), allocatable :: mpk_Wfull, mpk_covfull
    real(mcp), dimension(:), allocatable :: mpk_kfull  
    real(mcp), dimension(:,:), allocatable :: invcov_tmp
    character(len=64) region_string
    character(80) :: dummychar
    character z_char
    integer iz,count
    integer :: file_unit
    logical bad
    Type(TIniFile) :: Ini
    
    iopb = 0
    
    file_unit = new_file_unit()
    call Ini_Open_File(Ini, gname, file_unit, bad, .false.)
    if (bad) then
      write (*,*)  'Error opening dataset file '//trim(gname)
      stop
    end if

#ifndef WIGZ
       call MpiStop('mpk: edit makefile to have "EXTDATA = WIGZ" to inlude WiggleZ data')
#endif

 

! we've opened the file, now we have to break it up ourselves


    wmset%name = Ini_Read_String_File(Ini,'name') 
    wmset%redshift = Ini_Read_Real_File(Ini,'redshift',0.0)
    if(wmset%redshift.eq.0.0) then
       call MpiStop('mpk: failed  to read in WiggleZ redshift')
    end if
    
    Ini_fail_on_not_found = .false.
    wmset%use_set =.true.
    if (Feedback > 0) write (*,*) 'reading: '//trim(wmset%name)
    num_mpk_points_full = Ini_Read_Int_File(Ini,'num_mpk_points_full',0)
    if (num_mpk_points_full.eq.0) write(*,*) ' ERROR: parameter num_mpk_points_full not set'
    num_mpk_kbands_full = Ini_Read_Int_File(Ini,'num_mpk_kbands_full',0)
    if (num_mpk_kbands_full.eq.0) write(*,*) ' ERROR: parameter num_mpk_kbands_full not set'
    min_mpk_points_use = Ini_Read_Int_File(Ini,'min_mpk_points_use',1)
    min_mpk_kbands_use = Ini_Read_Int_File(Ini,'min_mpk_kbands_use',1)
    max_mpk_points_use = Ini_Read_Int_File(Ini,'max_mpk_points_use',num_mpk_points_full)
    max_mpk_kbands_use = Ini_Read_Int_File(Ini,'max_mpk_kbands_use',num_mpk_kbands_full)

! region 1 = 9h
! region 2 = 11h
! region 3 = 15h
! region 4 = 22h
! region 5 = 0h
! region 6 = 1h
! region 7 = 3h

    allocate(wmset%regions_active(max_num_wigglez_regions))
    do i_regions=1,7
       if(i_regions.eq.1) then
          region_string = 'Use_9-hr_region'
       else if(i_regions.eq.2) then
          region_string = 'Use_11-hr_region'
       else if(i_regions.eq.3) then
          region_string = 'Use_15-hr_region'
       else if(i_regions.eq.4) then
          region_string = 'Use_22-hr_region'
       else if(i_regions.eq.5) then
          region_string = 'Use_1-hr_region'
       else if(i_regions.eq.6) then
          region_string = 'Use_3-hr_region'
       else if(i_regions.eq.7) then
          region_string = 'Use_0-hr_region'
       endif
       wmset%regions_active(i_regions) =  Ini_Read_Logical_File(Ini,region_string,.false.)
    enddo

!  ... work out how many regions are being used

    wmset%num_regions_used = 0
    do i_regions = 1,max_num_wigglez_regions
       if(wmset%regions_active(i_regions)) wmset%num_regions_used = wmset%num_regions_used + 1
    enddo

    if(wmset%num_regions_used.eq.0) then
       print*, wmset%name
       call MpiStop('mpk: no regions begin used in this data set')
    endif

    wmset%num_mpk_points_use = max_mpk_points_use - min_mpk_points_use +1
    wmset%num_mpk_kbands_use = max_mpk_kbands_use - min_mpk_kbands_use +1

    if(allocated(mpk_kfull)) deallocate(mpk_kfull)
    allocate(mpk_kfull(num_mpk_kbands_full))
    allocate(wmset%mpk_P(wmset%num_regions_used,wmset%num_mpk_points_use))
    allocate(wmset%mpk_k(wmset%num_mpk_kbands_use))
    allocate(wmset%mpk_W(wmset%num_regions_used,wmset%num_mpk_points_use,wmset%num_mpk_kbands_use))



    kbands_file  = ReadIniFileName(Ini,'kbands_file')
    call ReadVector(kbands_file,mpk_kfull,num_mpk_kbands_full)
    wmset%mpk_k(:)=mpk_kfull(min_mpk_kbands_use:max_mpk_kbands_use) 
    if (Feedback > 1) then 
       write(*,*) 'reading: '//trim(wmset%name)//' data'
       write(*,*) 'Using kbands windows between',real(wmset%mpk_k(1)),' < k/h < ',real(wmset%mpk_k(wmset%num_mpk_kbands_use))      
    endif
    if  (wmset%mpk_k(1) < matter_power_minkh) then
       write (*,*) 'WARNING: k_min in '//trim(wmset%name)//'less than setting in cmbtypes.f90'
       write (*,*) 'all k<matter_power_minkh will be set to matter_power_minkh' 
    end if

    measurements_file  = ReadIniFileName(Ini,'measurements_file')
    call OpenTxtFile(measurements_file, tmp_file_unit)
    wmset%mpk_P=0.
    count = 0
    do i_regions =1,7
       if(wmset%regions_active(i_regions)) then
          count = count+1
          read (tmp_file_unit,*) dummychar
          read (tmp_file_unit,*) dummychar
          do i= 1, (min_mpk_points_use-1)
             read (tmp_file_unit,*, iostat=iopb) keff,klo,khi,beff,beff,beff
          end do

          if (Feedback > 1 .and. min_mpk_points_use>1) write(*,*) 'Not using bands with keff=  ',real(keff),&
               ' or below in region', i_regions
          do i =1, wmset%num_mpk_points_use
             read (tmp_file_unit,*, iostat=iopb) keff,klo,khi,wmset%mpk_P(count,i),beff,beff
          end do
          ! NB do something to get to the end of the list
          do i=1, num_mpk_points_full-wmset%num_mpk_points_use-min_mpk_points_use+1
             read (tmp_file_unit,*, iostat=iopb) klo,klo,khi,beff,beff,beff
             if(iopb.ne.0) stop
          end do
       else
          read (tmp_file_unit,*) dummychar
          read (tmp_file_unit,*) dummychar
          do i=1,50
             read (tmp_file_unit,*, iostat=iopb) klo,klo,khi,beff,beff,beff
             if(iopb.ne.0) stop
          enddo
       endif
    enddo
    close(tmp_file_unit)
    if (Feedback > 1) write(*,*) 'bands truncated at keff=  ',real(keff)

    allocate(mpk_Wfull(max_num_wigglez_regions,num_mpk_points_full,num_mpk_kbands_full))
    windows_file  = ReadIniFileName(Ini,'windows_file')
    if (windows_file.eq.'') write(*,*) 'ERROR: WiggleZ mpk windows_file not specified'
    call ReadWiggleZMatrices(windows_file,mpk_Wfull,max_num_wigglez_regions,num_mpk_points_full,num_mpk_kbands_full)
    count = 0
    do i_regions=1,max_num_wigglez_regions
       if(wmset%regions_active(i_regions)) then
          count = count + 1
          wmset%mpk_W(count,1:wmset%num_mpk_points_use,1:wmset%num_mpk_kbands_use)= &
               mpk_Wfull(i_regions,min_mpk_points_use:max_mpk_points_use,min_mpk_kbands_use:max_mpk_kbands_use)
       endif
    enddo
    
    deallocate(mpk_Wfull)
!    deallocate(mpk_kfull)
 
    cov_file  = ReadIniFileName(Ini,'cov_file')
    if (cov_file /= '') then
       allocate(mpk_covfull(max_num_wigglez_regions,num_mpk_points_full,num_mpk_points_full))
       allocate(invcov_tmp(wmset%num_mpk_points_use,wmset%num_mpk_points_use))
! ... read the entire covraiance matrix in, then decide which regions we want...
       call ReadWiggleZMatrices(cov_file,mpk_covfull,max_num_wigglez_regions,num_mpk_points_full,num_mpk_points_full)
       allocate(wmset%mpk_invcov(wmset%num_regions_used,wmset%num_mpk_points_use,wmset%num_mpk_points_use))
       count = 0
       do i_regions=1,max_num_wigglez_regions
          if(wmset%regions_active(i_regions)) then
             count = count + 1
! ... the covariance matrix has two indices for the different k-values, and another one for the region...       
!             wmset%mpk_invcov(count,1:wmset%num_mpk_points_use,1:wmset%num_mpk_points_use)=  &
             invcov_tmp(:,:) = &
                  mpk_covfull(i_regions,min_mpk_points_use:max_mpk_points_use,min_mpk_points_use:max_mpk_points_use)
!             call Matrix_Inverse(wmset%mpk_invcov(count,:,:))
             call Matrix_Inverse(invcov_tmp)
             wmset%mpk_invcov(count,1:wmset%num_mpk_points_use,1:wmset%num_mpk_points_use) = invcov_tmp(:,:)
          endif
       enddo
       deallocate(mpk_covfull)
       deallocate(invcov_tmp)
    else
       nullify(wmset%mpk_invcov)
    end if

    if (iopb.ne.0) then
       stop 'Error reading WiggleZ mpk file'
    endif

    wmset%use_scaling = Ini_Read_Logical_File(Ini,'use_scaling',.false.)
    wmset%use_gigglez = Ini_Read_Logical_File(Ini,'Use_gigglez',.false.)

    if(wmset%use_gigglez) then
      call GiggleZinfo_init(wmset%redshift)
    endif

!... IMPORTANT: need to talk to CB about Q marginalisation ...
    wmset%Q_marge = Ini_Read_Logical_File(Ini,'Q_marge',.false.)
    if (wmset%Q_marge) then
       wmset%Q_flat = Ini_Read_Logical_File(Ini,'Q_flat',.false.)
       if (.not. wmset%Q_flat) then
          !gaussian prior on Q
          wmset%Q_mid = Ini_Read_Real_File(Ini,'Q_mid')
          wmset%Q_sigma = Ini_Read_Real_File(Ini,'Q_sigma')
       end if
       wmset%Ag = Ini_Read_Real_File(Ini,'Ag', 1.4)
    end if 

    call Ini_Close_File(Ini)
    call ClearFileUnit(file_unit)

    if (Feedback > 1) write(*,*) 'read: '//trim(wmset%name)//' data'
 
 end subroutine ReadWiggleZDataset
 
 subroutine ReadWiggleZMatrices(aname,mat,num_regions,m,n)
! suborutine to read all the matrices from each of the different regions, enclosed in one file

   implicit none
   character(LEN=*), intent(IN) :: aname
   integer, intent(in) :: m,n,num_regions
   real(mcp), intent(out) :: mat(num_regions,m,n)
   integer j,k,i_region
   real(mcp) tmp
   character(LEN=64) dummychar



   if (Feedback > 1) write(*,*) 'reading: '//trim(aname)
   call OpenTxtFile(aname, tmp_file_unit)
   do i_region=1,num_regions
      read (tmp_file_unit,*, end = 200, err=100) dummychar
      do j=1,m
         read (tmp_file_unit,*, end = 200, err=100) mat(i_region,j,1:n)
      enddo
   enddo
   goto 120

100 write(*,*) 'matrix file '//trim(aname)//' is the wrong size',i_region,j,n,mat(num_regions,m,n)
   stop

120 read (tmp_file_unit,*, err = 200, end =200) tmp
   goto 200

   
200 close(tmp_file_unit)
   return
     

 end subroutine ReadWiggleZMatrices
 
 function WiggleZ_LnLike(like,CMB,Theory,DataParams) ! LV_06 added CMB here
   Class(CMBParams) CMB
   Class(WiggleZLikelihood) :: like
   Class(TheoryPredictions) Theory
   real(mcp) :: DataParams(:)   
   real(mcp) :: WiggleZ_LnLike, LnLike
   real(mcp), dimension(:), allocatable :: mpk_Pth, mpk_k2,mpk_lin,k_scaled !LV_06 added for LRGDR4
   real(mcp), dimension(:), allocatable :: mpk_WPth, mpk_WPth_k2
   real(mcp), dimension(:), allocatable :: diffs, step1
   real(mcp), dimension(:), allocatable :: Pk_delta_delta,Pk_delta_theta,Pk_theta_theta
   real(mcp), dimension(:), allocatable :: damp1, damp2, damp3
   real(mcp) :: covdat(like%num_mpk_points_use)
   real(mcp) :: covth(like%num_mpk_points_use)
   real(mcp) :: covth_k2(like%num_mpk_points_use)
   real(mcp), dimension(:), allocatable :: mpk_WPth_large, covdat_large, covth_large, mpk_Pdata_large
   integer imin,imax
   real :: normV, Q, minchisq
   real(mcp) :: a_scl  !LV_06 added for LRGDR4
   integer :: i, iQ,ibias,ik,j,iz
   logical :: do_marge
   integer, parameter :: nQ=6
   integer, parameter :: nbias = 100
   real(mcp) b0, bias_max, bias_step,bias,old_chisq,beta_val,kval,xi
   real(mcp) :: tmp, dQ = 0.4
   real(mcp), dimension(:), allocatable :: chisq(:)
   real(mcp) calweights(-nQ:nQ)
   real(mcp) vec2(2),Mat(2,2)
   real(mcp) final_term, b_out
   real(mcp) z,omk_fid, omv_fid,w_fid
   integer i_region
   character(len=32) fname

   If(Feedback > 1) print*, 'Calling WiggleZ likelihood routines'
   allocate(mpk_lin(like%num_mpk_kbands_use),mpk_Pth(like%num_mpk_kbands_use))
   allocate(mpk_WPth(like%num_mpk_points_use))
   allocate(k_scaled(like%num_mpk_kbands_use))!LV_06 added for LRGDR4 !! IMPORTANT: need to check k-scaling
!   allocate(like%num_mpk_points_use))
   allocate(diffs(like%num_mpk_points_use),step1(like%num_mpk_points_use))
   allocate(Pk_delta_delta(like%num_mpk_kbands_use),Pk_delta_theta(like%num_mpk_kbands_use))
   allocate(Pk_theta_theta(like%num_mpk_kbands_use))
   allocate(damp1(like%num_mpk_kbands_use),damp2(like%num_mpk_kbands_use),damp3(like%num_mpk_kbands_use))

   
   allocate(chisq(-nQ:nQ))
   If(CMB%wa/=0) call MpiStop('WiggleZ MPK module not compatible with wa /= 0')
   chisq = 0

   if (.not. like%use_set) then
      LnLike = 0
      return
   end if

   ! won't actually want to do this multiple times for multiple galaxy pk data sets?..
   omk_fid = 0.d0
   omv_fid = 0.705d0
   w_fid = -1.d0
   z = 1.d0*dble(like%redshift) ! accuracy issues
   if(like%use_scaling) then
     call compute_scaling_factor_z(z,dble(CMB%omk),dble(CMB%omv),dble(CMB%w),omk_fid,omv_fid,w_fid,a_scl)
   else
     a_scl = 1
   end if
   
   iz = 0
   do i=1,4
      if(abs(z-zeval(i)).le.0.001) iz = i
   enddo
   if(iz.eq.0) call MpiStop('could not indentify redshift')

   if(like%use_gigglez) then
      call fill_GiggleZTheory(Theory,matter_power_minkh,matter_power_dlnkh,z)
   endif

   do i=1, like%num_mpk_kbands_use 
      ! It could be that when we scale the k-values, the lowest bin drops off the bottom edge
      !Errors from using matter_power_minkh at lower end should be negligible
         k_scaled(i)=max(matter_power_minkh,like%mpk_k(i)*a_scl)
         if(like%use_gigglez) then
            mpk_lin(i) = WiggleZPowerAt(Theory,k_scaled(i))/a_scl**3
         else
            mpk_lin(i)=MatterPowerAt_zbin(Theory,k_scaled(i),izwigglez(iz))/a_scl**3
         endif
   end do
   
    do_marge = like%Q_Marge
    if (do_marge .and. like%Q_flat) then
       !Marginalize analytically with flat prior on b^2 and b^2*Q
       !as recommended by Max Tegmark for SDSS
       allocate(mpk_k2(like%num_mpk_kbands_use))
       allocate(mpk_WPth_k2(like%num_mpk_points_use))
       
       Mat(:,:) = 0.d0
       vec2(:) = 0.d0
       final_term = 0.d0
       do i_region=1,like%num_regions_used
          mpk_Pth(:)=mpk_lin(:)/(1+like%Ag*k_scaled)
          mpk_k2(:)=mpk_Pth(:)*k_scaled(:)**2
          
          
          mpk_WPth(:) = matmul(like%mpk_W(i_region,:,:),mpk_Pth(:))
          mpk_WPth_k2(:) = matmul(like%mpk_W(i_region,:,:),mpk_k2(:))
          
          covdat(:) = matmul(like%mpk_invcov(i_region,:,:),like%mpk_P(i_region,:))
          covth(:) = matmul(like%mpk_invcov(i_region,:,:),mpk_WPth(:))
          covth_k2(:) = matmul(like%mpk_invcov(i_region,:,:),mpk_WPth_k2(:))
          
          Mat(1,1) = Mat(1,1) + sum(covth(:)*mpk_WPth(:))
          Mat(2,2) = Mat(2,2) + sum(covth_k2(:)*mpk_WPth_k2(:))
          Mat(1,2) = Mat(1,2) + sum(covth(:)*mpk_WPth_k2(:))
          Mat(2,1) = Mat(1,2)

          vec2(1) = vec2(1) + sum(covdat(:)*mpk_WPth(:))
          vec2(2) = vec2(2) + sum(covdat(:)*mpk_WPth_k2(:))
          final_term = final_term + sum(like%mpk_P(i_region,:)*covdat(:))
       enddo
       LnLike = log( Mat(1,1)*Mat(2,2)-Mat(1,2)**2)
       call inv_mat22(Mat)
       !          LnLike = (sum(mset%mpk_P*covdat) - sum(vec2*matmul(Mat,vec2)) + LnLike ) /2
       LnLike = (final_term - sum(vec2*matmul(Mat,vec2)) + LnLike ) /2

       deallocate(mpk_k2,mpk_WPth_k2)      
    else
       if (like%Q_sigma==0) do_marge = .false.
       ! ... sum the chi-squared contributions for all regions first
       chisq(:) = 0.d0
       old_chisq = 1.d30
       if(feedback > 1) print*, "starting analytic marginalisation over bias"
       allocate(mpk_Pdata_large(like%num_mpk_points_use*like%num_regions_used))
       allocate(mpk_WPth_large(like%num_mpk_points_use*like%num_regions_used))
       allocate(covdat_large(like%num_mpk_points_use*like%num_regions_used))       
       allocate(covth_large(like%num_mpk_points_use*like%num_regions_used))
       normV = 0.d0
       do iQ=-nQ,nQ
          Q = like%Q_mid +iQ*like%Q_sigma*dQ 
          if (like%Q_marge) then
             mpk_Pth(:)=mpk_lin(:)*(1+Q*k_scaled(:)**2)/(1+like%Ag*k_scaled(:))
          else 
             mpk_Pth(:) = mpk_lin(:)
          end if
          do i_region=1,like%num_regions_used
             imin = (i_region-1)*like%num_mpk_points_use+1
             imax = i_region*like%num_mpk_points_use
             mpk_WPth(:) = matmul(like%mpk_W(i_region,:,:),mpk_Pth(:))
             mpk_Pdata_large(imin:imax) = like%mpk_P(i_region,:)
             mpk_WPth_large(imin:imax) = mpk_WPth(:) 
             
             !with analytic marginalization over normalization nuisance (flat prior on b^2)
             !See appendix F of cosmomc paper
             
             covdat_large(imin:imax) = matmul(like%mpk_invcov(i_region,:,:),like%mpk_P(i_region,:))
             covth_large(imin:imax) = matmul(like%mpk_invcov(i_region,:,:),mpk_WPth(:))
          enddo
          normV = normV + sum(mpk_WPth_large*covth_large)
          b_out =  sum(mpk_WPth_large*covdat_large)/sum(mpk_WPth_large*covth_large)
          if(Feedback.ge.2) print*, "Bias value:", b_out
          chisq(iQ) = sum(mpk_Pdata_large*covdat_large)  - sum(mpk_WPth_large*covdat_large)**2/normV!  + log(normV)
                   
          if (do_marge) then
             calweights(iQ) = exp(-(iQ*dQ)**2/2)
          else 
             LnLike = chisq(iQ)/2
             exit
          end if
          
       end do
       deallocate(covdat_large,covth_large,mpk_Pdata_large,mpk_WPth_large)
    
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

    end if !not analytic over Q
   WiggleZ_LnLike=LnLike   
   if (Feedback>1) write(*,*) 'WiggleZ mpk chi-sq:', LnLike*2
   
   if (LnLike > 1e8) then
      write(*,*) 'Chisq is huge, maybe there is a problem? chisq=',chisq
   end if
   
   deallocate(mpk_Pth,mpk_lin)
   deallocate(mpk_WPth,k_scaled)!,w)
   deallocate(chisq)
   
 end function WiggleZ_LnLike 
 
  

 subroutine inv_mat22(M)
    real(mcp) M(2,2), Minv(2,2), det

    det = M(1,1)*M(2,2)-M(1,2)*M(2,1)
    Minv(1,1)=M(2,2)
    Minv(2,2) = M(1,1)
    Minv(1,2) = - M(2,1)
    Minv(2,1) = - M(1,2)
    M = Minv/det

 end subroutine inv_mat22

!-----------------------------------------------------------------------------
!LV added to include lrg DR4

!DP added an input zeff value to compute_scaling_factor
subroutine compute_scaling_factor_z(z,Ok,Ol,w,Ok0,Ol0,w0,a)
  ! a = dV for z=0.35 relative to its value for flat Om=0.25 model.
  ! This is the factor by which the P(k) measurement would shift 
  ! sideways relative to what we got for this fiducial flat model.
  ! * a = (a_angular**2 * a_radial)**(1/3)
  ! * a_angular = comoving distance to z=0.35 in Mpc/h relative to its value for flat Om=0.25 model
  !     dA = (c/H)*eta = (2997.92458 Mpc/h)*eta, so we care only about 
  !     eta scaling, not h scaling.
  !     For flat w=-1 models, a ~ (Om/0.25)**(-0.065)
  !     For the LRG mean redshift z=0.35, the power law fit 
  !    dA(z,Om= 0.3253 (Om/0.25)^{-0.065}c H_0^{-1} is quite good within 
  !    our range of interest,
  !     accurate to within about 0.1% for 0.2<Om<0.3.
  ! * a_radial = 1/H(z) relative to its value for flat Om=0.25 model
  implicit none
  real(mpk_d) Or, Om, Ok, Ol, w, Ok0, Om0, Ol0, w0, z, eta, eta0, Hrelinv, Hrelinv0, tmp
  real(mpk_d) a_radial, a_angular
  real(mcp) a
  !Or= 0.0000415996*(T_0/2.726)**4 / h**2
  Or= 0! Radiation density totally negligible at  z < 0.35
  Om= 1-Ok-Ol-Or
  !!!z  = 0.35  !!edited by Beth 21-11-08: change to zeff of Will's LRG sample.
 
  Hrelinv= 1/sqrt(Ol*(1+z)**(3*(1+w)) + Ok*(1+z)**2 + Om*(1+z)**3 + Or*(1+z)**4)
!  write(*,*) Ok,Ol,w  
  call compute_z_eta(Or,Ok,Ol,w,z,eta)
  tmp = sqrt(abs(Ok))
  if (Ok.lt.-1.d-6) eta = sin(tmp*eta)/tmp
  if (Ok.gt.1d-6)   eta = (exp(tmp*eta)-exp(-tmp*eta))/(2*tmp) ! sinh(tmp*eta)/tmp
!  Ok0= 0
!  Ol0= 0.75
!  w0= -1
  Om0= 1-Ok0-Ol0-Or
  call compute_z_eta(Or,Ok0,Ol0,w0,z,eta0)
  Hrelinv0= 1/sqrt(Ol0*(1+z)**(3*(1+w0)) + Ok0*(1+z)**2 + Om0*(1+z)**3 + Or*(1+z)**4)
  !a_angular = (Om/0.25)**(-0.065) * (-w*Otot)**0.14 ! Approximation based on Taylor expansion
  a_angular = eta/eta0
  a_radial= Hrelinv/Hrelinv0
  a=  (a_angular**2 * a_radial)**(1/3.d0)
  !write(*,'(9f10.5)') Ok,Ol,w,a,a_radial,a_angular,(Om/0.25)**(-0.065) * (-w*(1-Ok))**0.14
  !write(*,'(9f10.5)') Ok,Ol,w,a,a_radial**(2/3.d0),a_angular**(4/3.d0),((Om/0.25)**(-0.065) * (-w*(1-Ok))**0.14)**(4/3.d0)
  !!! BR09 -- in previous version, scale factor is applied in the wrong direction!  So take a = 1/a to fix it.
  a = 1.0/a
end subroutine compute_scaling_factor_z

subroutine eta_demo
  implicit none
  real(mpk_d) Or, Ok, Ol, w, h, z, eta
  h  = 0.7
  Ok = 0
  Ol = 0.7
  Or = 0.0000416/h**2 
  w  = -1
  z  = 1090
  call compute_z_eta(Or,Ok,Ol,w,z,eta)
!  print *,'eta.............',eta
!  print *,'dlss in Gpc.....',(2.99792458/h)*eta
end subroutine eta_demo

!INTERFACE
logical function nobigbang2(Ok,Ol,w)
  ! Test if we're in the forbidden zone where the integrand blows up
  ! (where's H^2 < 0 for some a between 0 and 1).
  ! The function f(a) = Omega_m + Omega_k*a + Omega_l*a**(-3*w)
  ! can have at most one local minimum, at (Ok/(3*w*Ol))**(-1/(1+3*w)), 
  ! so simply check if f(a)<0 there or at the endpoints a=0, a=1.
  ! g(0) = Omega_m - Omega_l*a**(-3*w) < 0 if w > 0 & Omega_k > 1
  !                                     or if w = 0 & Omega_l < 1       
  ! g(1) = Omega_m + Omega_k + Omega_l = 1 > 0
  implicit none
  real(mpk_d) Ok, Ol, w, Om, tmp, a, epsilon
  integer failure
  failure = 0
  epsilon = 0
  !epsilon = 0.04  ! Numerical integration fails even before H^2 goes negative.
  Om = 1.d0 - Ok - Ol
  if (w*Ol.ne.0) then
     tmp = Ok/(3*w*Ol)
     if ((tmp.gt.0).and.(1+3*w.ne.0)) then ! f'(0)=0 for some a>0
        a = tmp**(-1/(1+3*w))
        if (a.lt.1) then
           if (Om + Ok*a + Ol*a**(-3*w).lt.epsilon) failure = 1
        end if
     end if
  end if
  if ((w.eq.0).and.(Ok.gt.1)) failure = 2
  if ((w.gt.0).and.(Ol.lt.0)) failure = 3
  nobigbang2 = (failure.gt.0)
  if (failure.gt.0) print *,'Big Bang failure mode ',failure
  return
end function nobigbang2
!END INTERFACE

real(mpk_d) function eta_integrand(a)
  implicit none
  real(mpk_d) Or, Ok, Ox, w
  common/eta/Or, Ok, Ox, w
  real(mpk_d) a, Om
  ! eta = int (H0/H)dz = int (H0/H)(1+z)dln(1+z) = int (H0/H)/a dlna = int (H0/H)/a^2 da = 
  ! Integrand = (H0/H)/a^2
  ! (H/H0)**2 = Ox*a**(-3*(1+w)) + Ok/a**2 + Om/a**3 + Or/a**4 
  if (a.eq.0.d0) then 
     eta_integrand = 0.d0
  else
     Om = 1.d0 - Or - Ok - Ox
     eta_integrand = 1.d0/sqrt(Ox*a**(1-3*w) + Ok*a**2 + Om*a + Or)
  end if
  return
end function eta_integrand

subroutine eta_z_integral(Omega_r,Omega_k,Omega_x,w_eos,z,eta)
  ! Computes eta as a function
  ! of the curvature Omega_k, the dark energy density Omega_x
  ! and its equation of state w.
  implicit none
  real(mpk_d) Or, Ok, Ox, w
  common/eta/Or, Ok, Ox, w
  real(mpk_d) Omega_r, Omega_k,Omega_x,w_eos, z, eta, epsabs, epsrel, amin, amax!, eta_integrand
  Or = Omega_r
  Ok = Omega_k
  Ox = Omega_x
  w  = w_eos
  epsabs  = 0
  epsrel  = 1.d-10
  amin= 1/(1+z)
  amax= 1
  call qromb2(eta_integrand,amin,amax,epsabs,epsrel,eta)
  return
end subroutine eta_z_integral

subroutine compute_z_eta(Or,Ok,Ox,w,z,eta)
  ! Computes the conformal distance eta(z)
  implicit none
  real(mpk_d) Or, Ok, Ox, w, z, eta
!  logical nobigbang2
  if (nobigbang2(Ok,Ox,w)) then
     print *,'No big bang, so eta undefined if z>zmax.'
     eta = 99 
  else
     call eta_z_integral(Or,Ok,Ox,w,z,eta) 
     ! print *,'Or, Ok, Ox, w, z, H_0 t_0...',Or, Ok, Ox, w, eta
  end if
  return
end subroutine compute_z_eta


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! num rec routines
!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE qromb2(func,a,b,epsabs,epsrel,ss)
! The numerical recipes routine, but modified so that is decides
! it's done when either the relative OR the absolute accuracy has been attained.
! The old version used relative errors only, so it always failed when
! when the integrand was near zero.
! epsabs = epsrel = 1e-6 are canonical choices.
  INTEGER JMAX,JMAXP,K,KM
  real(mpk_d) a,b,func,ss,epsabs,epsrel
  EXTERNAL func
  PARAMETER (JMAX=20, JMAXP=JMAX+1, K=5, KM=K-1)
                                !    USES polint,trapzd
  INTEGER j
  real(mpk_d) dss,h(JMAXP),s(JMAXP)
  h(1)=1.d0
  do j=1,JMAX
     call trapzd(func,a,b,s(j),j)
     if (j.ge.K) then
        call polint(h(j-KM),s(j-KM),K,0.d0,ss,dss)
        if (abs(dss).le.epsrel*abs(ss)) return
        if (abs(dss).le.epsabs) return
     endif
     s(j+1)=s(j)
     h(j+1)=0.25d0*h(j)
  ENDDO
  print *,'Too many steps in qromb'
      
  RETURN 
END SUBROUTINE qromb2
  
SUBROUTINE polint(xa,ya,n,x,y,dy) ! From Numerical Recipes
  INTEGER n,NMAX
  real(mpk_d) dy,x,y,xa(n),ya(n)
  PARAMETER (NMAX=10)
  INTEGER i,m,ns
  real(mpk_d) den,dif,dift,ho,hp,w,c(NMAX),d(NMAX)
  ns=1
  dif=abs(x-xa(1))
  do  i=1,n
     dift=abs(x-xa(i))
     if (dift.lt.dif) then
        ns=i
        dif=dift
     endif
     c(i)=ya(i)
     d(i)=ya(i)
  enddo
  y=ya(ns)
  ns=ns-1
  do  m=1,n-1
     do  i=1,n-m
        ho=xa(i)-x
        hp=xa(i+m)-x
        w=c(i+1)-d(i)
        den=ho-hp
        if(den.eq.0.) then
           print*, 'failure in polint'
           stop
        endif
        den=w/den
        d(i)=hp*den
        c(i)=ho*den
     enddo
     if (2*ns.lt.n-m)then
        dy=c(ns+1)
     else
        dy=d(ns)
        ns=ns-1
     endif
     y=y+dy
  enddo
  return
END SUBROUTINE polint
      
SUBROUTINE trapzd(func,a,b,s,n) ! From Numerical Recipes
  INTEGER n
  real(mpk_d) a,b,s,func
  EXTERNAL func
  INTEGER it,j
  real(mpk_d) del,sum,tnm,x
  if (n.eq.1) then
     s=0.5*(b-a)*(func(a)+func(b))
  else
     it=2**(n-2)
     tnm=it
     del=(b-a)/tnm
     x=a+0.5*del
     sum=0.
     do  j=1,it
        sum=sum+func(x)
        x=x+del
     enddo
     s=0.5*(s+(b-a)*sum/tnm)
  endif
  return
END SUBROUTINE trapzd

end module
