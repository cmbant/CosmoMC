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


module LRGinfo
use settings
use cmbtypes
use Precision
use lrggettheory
!use CMB_Cls

implicit none

!! these are the LRG redshift subsample weights.
real(dl), parameter :: w0 = 0.0d0, wNEAR = 0.395d0, wMID = 0.355d0, wFAR = 0.250d0

!in CAMB: 4=now (z=0), 3=NEAR, 2=MID, 1=FAR; opposite order in matter_power
!! now generalized indices iz0lrg, izNEARlrg, izMIDlrg, izFARlrg
real(dl), dimension(4) :: zeval, zweight, sigma2BAOfid, sigma2BAO

real(dl) om_val, ol_val, ok_check, wval  ! passed in from CMBparams CMB

! power spectra evaluated at fiducial cosmological theory (WMAP5 recommended values)
real, allocatable :: ratio_power_nw_nl_fid(:,:)
!real,dimension(num_matter_power,matter_power_lnzsteps) :: ratio_power_nw_nl_fid
!make allocatable to avoid compile-time range errors when matter_power_lnzsteps<4

contains

subroutine LRGinfo_init()
  integer :: iopb, i, ios
  real(dl) :: omegakdummy,omegavdummy,wdummy,getabstransferscalefiddummy
  real(dl) :: kval, plin, psmooth, rationwhalofit

  !!BR09 only needed for LRGs, so only 4 redshifts no matter what matter_power_lnzsteps is
  allocate(ratio_power_nw_nl_fid(num_matter_power,4))

  sigma2BAOfid(1) = 1.0e-5  !! don't do any smearing at z=0; this won't be used anyway.
  sigma2BAOfid(2) = sigma2BAONEAR
  sigma2BAOfid(3) = sigma2BAOMID
  sigma2BAOfid(4) = sigma2BAOFAR

  zeval(1) = z0
  zeval(2) = zNEAR
  zeval(3) = zMID
  zeval(4) = zFAR

  zweight(1) = w0
  zweight(2) = wNEAR
  zweight(3) = wMID
  zweight(4) = wFAR

  !! first read in everything needed from the CAMB output files.
  iopb = 0 !! check later if there was an error

  open(unit=tmp_file_unit,file=trim(DataDir)//'new_lrgdr7fiducialmodel_matterpowerzNEAR.dat',form='formatted',err=500, iostat=ios)
  read (tmp_file_unit,*,iostat=iopb) getabstransferscalefiddummy, omegakdummy,omegavdummy,wdummy
  do i = 1, num_matter_power
    read (tmp_file_unit,*,iostat=iopb) kval, plin, psmooth, rationwhalofit
    ratio_power_nw_nl_fid(i,2) = rationwhalofit
  end do
  close(tmp_file_unit)

  open(unit=tmp_file_unit,file=trim(DataDir)//'new_lrgdr7fiducialmodel_matterpowerzMID.dat',form='formatted',err=500, iostat=ios)
  read (tmp_file_unit,*,iostat=iopb) getabstransferscalefiddummy,omegakdummy,omegavdummy,wdummy
  do i = 1, num_matter_power
    read (tmp_file_unit,*,iostat=iopb) kval, plin, psmooth, rationwhalofit
    ratio_power_nw_nl_fid(i,3) = rationwhalofit
  end do
  close(tmp_file_unit)

  open(unit=tmp_file_unit,file=trim(DataDir)//'new_lrgdr7fiducialmodel_matterpowerzFAR.dat',form='formatted',err=500,iostat=ios)
  read (tmp_file_unit,*,iostat=iopb) getabstransferscalefiddummy,omegakdummy,omegavdummy,wdummy
  do i = 1, num_matter_power
    read (tmp_file_unit,*,iostat=iopb) kval, plin, psmooth, rationwhalofit
    ratio_power_nw_nl_fid(i,4) = rationwhalofit
  end do
  close(tmp_file_unit)

500 if(ios .ne. 0) stop 'Unable to open file'
  if(iopb .ne. 0) stop 'Error reading model or fiducial theory files.'
end subroutine LRGinfo_init

! HARD CODING OF POLYNOMIAL FITS TO NEAR, MID, FAR SUBSAMPLES.
subroutine LRGtoICsmooth(k,fidpolys)
  real(dl), intent(in) :: k
  real(dl) :: fidNEAR, fidMID, fidFAR
  real(dl), dimension(2:4), intent(out) :: fidpolys

  if(k < 0.194055d0) then !!this is where the two polynomials are equal
    fidNEAR = (1.0d0 - 0.680886d0*k + 6.48151d0*k**2)
  else
    fidNEAR = (1.0d0 - 2.13627d0*k + 21.0537d0*k**2 - 50.1167d0*k**3 + 36.8155d0*k**4)*1.04482d0
  end if

  if(k < 0.19431) then
    fidMID = (1.0d0 - 0.530799d0*k + 6.31822d0*k**2)
  else
    fidMID = (1.0d0 - 1.97873d0*k + 20.8551d0*k**2 - 50.0376d0*k**3 + 36.4056d0*k**4)*1.04384
  end if

  if(k < 0.19148) then
    fidFAR = (1.0d0 - 0.475028d0*k + 6.69004d0*k**2)
  else
    fidFAR = (1.0d0 - 1.84891d0*k + 21.3479d0*k**2 - 52.4846d0*k**3 + 38.9541d0*k**4)*1.03753
  end if
  fidpolys(2) = fidNEAR
  fidpolys(3) = fidMID
  fidpolys(4) = fidFAR
end subroutine LRGtoICsmooth

subroutine fill_LRGTheory(Theory, minkh, dlnkh)
  Type(CosmoTheory) Theory
  real, intent(in) :: minkh, dlnkh
  real(dl) :: logmink, xi, kval, expval, psmear, nlrat
  real(dl), dimension(2:4) :: fidpolys, holdval

  integer :: iz, ik, matterpowerindx

  do iz = 1, 4
    sigma2BAO(iz) = sigma2BAOfid(iz)
  end do

  logmink = log(minkh)
  do ik=1,num_matter_power
   xi = logmink + dlnkh*(ik-1)
   kval = exp(xi)
   Theory%finalLRGtheoryPk(ik) = 0.
   do iz = 2,4
    if(iz == 2) matterpowerindx = izNEARlrg
    if(iz == 3) matterpowerindx = izMIDlrg
    if(iz == 4) matterpowerindx = izFARlrg
    expval = exp(-kval**2*sigma2BAO(iz)*0.5)
    psmear = (Theory%matter_power(ik,matterpowerindx))*expval + (Theory%mpk_nw(ik,matterpowerindx))*(1.0-expval)
    psmear = psmear*powerscaletoz0(iz)
    nlrat = (Theory%mpkrat_nw_nl(ik,matterpowerindx))/(ratio_power_nw_nl_fid(ik,iz))
    call LRGtoICsmooth(kval,fidpolys)
    holdval(iz) = zweight(iz)*psmear*nlrat*fidpolys(iz)
    Theory%finalLRGtheoryPk(ik) = Theory%finalLRGtheoryPk(ik) + holdval(iz)
   end do

  end do

end subroutine fill_LRGTheory

end module LRGinfo

module wigglezinfo
!David Parkinson 12th March 2012
  use settings
  use cmbtypes
  use Precision
!use lrggettheory

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
!real,dimension(num_matter_power,matter_power_lnzsteps) :: ratio_power_nw_nl_fid
!make allocatable to avoid compile-time range errors when matter_power_lnzsteps<4
  logical :: use_wigz10 = .false.

contains

  subroutine GiggleZinfo_init(redshift)
    integer :: iopb, i, ios, iz
    real(dl) :: kval, power_nl
    real redshift
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
    real(dl), intent(in) :: k
    real(dl) :: fidz_0, fidz_1, fidz_2, fidz_3, fidz_4
    real(dl), dimension(4), intent(out) :: fidpolys

    !This is z=0, put here for the sake of completeness, not to be used
!    fidz_0 = (4.54679d0 - 11.2413d0*k + 30.0297d0*k**2 - 34.9125d0*k**3)

    fidz_1 = (4.619d0 - 13.7787d0*k + 58.941d0*k**2 - 175.24d0*k**3 + 284.321d0*k**4 - 187.284d0*k**5)
    fidz_2 = (4.63079d0 - 12.6293d0*k + 42.9265d0*k**2 - 91.8068d0*k**3 + 97.808d0*k**4 - 37.633d0*k**5)
    fidz_3 = (4.69659d0 - 12.7287d0*k + 42.5681d0*k**2 - 89.5578d0*k**3 + 96.664d0*k**4 - 41.2564*k**5)
    fidz_4 = (4.6849d0 - 13.4747d0*k + 53.7172d0*k**2 - 145.832d0*k**3 + 216.638d0*k**4 - 132.782*k**5)

!    fidz_1 = (4.46156d0 - 9.41743d0*k + 18.6255d0*k**2 - 15.555d0*k**3)
    
!    fidz_2 = (4.50568d0 - 9.41886d0*k + 18.2398d0*k**2 - 14.9585d0*k**3)

!    fidz_3 = (4.5837d0 - 9.80682d0*k + 19.777d0*k**2 - 16.8924d0*k**3)
    
!    fidz_4 = (4.5447d0 - 9.51283d0*k + 18.9695d0*k**2 - 16.0172d0*k**3)
    
    
    fidpolys(1) = 10**fidz_1
    fidpolys(2) = 10**fidz_2
    fidpolys(3) = 10**fidz_3
    fidpolys(4) = 10**fidz_4
    return

  end subroutine GiggleZtoICsmooth

  subroutine fill_GiggleZTheory(Theory, minkh, dlnkh,z)
    Type(CosmoTheory) Theory
    real, intent(in) :: minkh, dlnkh
    real(dl), intent(in) :: z
    real(dl) :: logmink, xi, kval, expval,  nlrat
    real(dl), dimension(4) :: fidpolys
    real(dl) y, dz, matter_power_dlnz
    real(dl) pk_hf, holdval
    integer :: i,iz, iz2, ik
    character(len=32) fname
    logmink = log(minkh)

    iz = 0
    do i=1,4
       if(abs(z-zeval(i)).le.0.001) iz = i
    enddo

!    write(fname,'(a,i2.2,a)') 'pk_hf_',int(z*100.),'.dat'
!    open(unit=60,file=fname,status='unknown')
!    write(fname,'(a,i2.2,a)') 'pk_final_',int(z*100.d0),'.dat'
!    open(unit=61,file=fname,status='unknown')
!    write(fname,'(a,i2.2,a)') 'pk_poly_',int(z*100.),'.dat'
!    open(unit=62,file=fname,status='unknown')


    do ik=1,num_matter_power
       xi = logmink + dlnkh*(ik-1)
       kval = exp(xi)
       Theory%finalLRGtheoryPk(ik) = 0.
       pk_hf = Theory%matter_power(ik,izwigglez(iz))
       
       ! matter_power from Theory is non-linear ! DRP 22-June-2012
!       write(60,*) kval, pk_hf

       call GiggleZtoICsmooth(kval,fidpolys)
!       write(61,*) kval, real(fidpolys(iz))

       holdval = pk_hf*fidpolys(iz)/power_hf_fid(ik,iz)
       ! we're using the finalLRGtheory data structure, even though these aren't actually LRGs
       Theory%finalLRGtheoryPk(ik) = Theory%finalLRGtheoryPk(ik) + holdval
!       write(62,*) kval, Theory%finalLRGtheoryPk(ik), holdval       

    end do
!    close(60)
!    close(61)
!    close(62)

  end subroutine fill_GiggleZTheory


end module wigglezinfo


module mpk
use precision
use settings
use cmbtypes
use LRGinfo
implicit none

 Type mpkdataset
    logical :: use_set
    integer :: num_mpk_points_use ! total number of points used (ie. max-min+1)
    integer :: num_mpk_kbands_use ! total number of kbands used (ie. max-min+1)
    character(LEN=20) :: name
    real, pointer, dimension(:,:) :: N_inv
    real, pointer, dimension(:,:) :: mpk_W, mpk_invcov
    real, pointer, dimension(:) :: mpk_P, mpk_sdev, mpk_k
    real, pointer, dimension(:) :: mpk_zerowindowfxn
    real, pointer, dimension(:) :: mpk_zerowindowfxnsubtractdat
    real :: mpk_zerowindowfxnsubdatnorm !!the 0th entry in windowfxnsubtract file
    logical :: use_scaling !as SDSS_lrgDR3
   !for Q and A see e.g. astro-ph/0501174, astro-ph/0604335
    logical :: Q_marge, Q_flat
    real :: Q_mid, Q_sigma, Ag
   end Type mpkdataset

 integer :: num_mpk_datasets = 0
 Type(mpkdataset) mpkdatasets(10)

 Type wigglez_mpkdataset
    logical :: use_set
    integer :: num_mpk_points_use ! total number of points used (ie. max-min+1)
    integer :: num_mpk_kbands_use ! total number of kbands used (ie. max-min+1)
    integer :: num_regions_used   ! total number of wigglez regions being used
    character(LEN=20) :: name
! 1st index always refers to the region
! so mpk_P(1,:) is the Power spectrum in the first active region
    real, pointer, dimension(:,:,:) :: mpk_W, mpk_invcov
    real, pointer, dimension(:,:) :: mpk_P
    real, pointer, dimension(:) :: mpk_k
!    real, pointer, dimension(:,:) :: mpk_zerowindowfxn
!    real, pointer, dimension(:,:) :: mpk_zerowindowfxnsubtractdat
!    real :: mpk_zerowindowfxnsubdatnorm !!the 0th entry in windowfxnsubtract file
    logical :: use_scaling !as SDSS_lrgDR3
    logical, pointer, dimension(:) :: regions_active
    logical use_jennings,use_simpledamping,use_gigglez
   !for Q and A see e.g. astro-ph/0501174, astro-ph/0604335
    logical :: Q_marge, Q_flat
    real :: Q_mid, Q_sigma, Ag, damp_sigv
    real :: redshift ! important to know
 end Type wigglez_mpkdataset

 integer :: num_wigglez_mpk_datasets = 0
 Type(wigglez_mpkdataset) wmpkdatasets(10)

 integer, parameter :: max_num_wigglez_regions = 7

 !Note all units are in k/h here
 
  integer, parameter :: mpk_d = kind(1.d0)
 
  logical :: use_mpk = .false.
  real, parameter :: Pi_num = 3.1415926535 
  ! constants describing the allowed a1,a2 regions.
  ! must check the functions below before changing these, because the shape of the space may change!

  integer, parameter :: wp = selected_real_kind(11,99)

  !!these are the 'nonconservative' nuisance parameter bounds
  !!real(dl), parameter :: k1 = 0.1d0, k2 = 0.2d0, s1 = 0.02d0, s2 = 0.05d0, a1maxval = 0.5741d0
  real(dl), parameter :: k1 = 0.1d0, k2 = 0.2d0, s1 = 0.04d0, s2 = 0.10d0, a1maxval = 1.1482d0
  integer, parameter :: nptsa1 = 41, nptsa2 = 41, nptstot = 325  
  !! but total number of points to evaluate is much smaller than 41**2 because lots of the space 
  !is not allowed by the s1,s2 constraints.
 
   ! only want to compute these once. 
   real(dl), dimension(nptstot) :: a1list, a2list

  real sigma_v, k_eval, growth_eval, growth_rate_eval, Omega_mat_growth, Omega_lam_growth, Omega_k_growth

contains 

  subroutine mpk_SetTransferRedshifts(redshifts)
   real, intent(inout) :: redshifts(*)
   !input is default log z spacing; can change here; check for consistency with other (e.g. lya)
         
   !Note internal ordering in CAMB is the opposite to that used in cosmomc transfer arrays (as here)
   !first index here must be redshift zero
   
      if(use_dr7lrg .and. matter_power_lnzsteps < 4) &
       call MpiStop('For LRGs matter_power_lnzsteps should be set to at least 4 (hardcoded in cmbtypes)')

       if (matter_power_lnzsteps==1 .or. .not. use_dr7lrg) return

       !! assigning indices to LRG NEAR, MID, FAR.  If you want to reorder redshifts, just change here.
       iz0lrg = 1  !! we use the z=0 output to normalize things; this is already assumed index 1 elsewhere 
                    !(like in calculation of sigma8).
       izNEARlrg = 2
       izMIDlrg = 3
       izFARlrg = 4
       redshifts(izNEARlrg) = zNEAR
       redshifts(izMIDlrg) = zMID
       redshifts(izFARlrg) = zFAR  
       if(iz0lrg /= 1) then
          redshifts(iz0lrg) = 0.0d0
       else
          if(redshifts(1) > 0.001) call MpiStop('redshifts(1) should be at z=0!')
       endif
         
  end subroutine mpk_SetTransferRedshifts

 subroutine wigglez_sdss_SetTransferRedshifts(redshifts)
   use wigglezinfo
   implicit none
   real, intent(inout) :: redshifts(matter_power_lnzsteps)
   integer :: zix_dum,extra_iz
   real red_dum
   !input is default log z spacing; can change here; check for consistency with other (e.g. lya)
         
   if(use_dr7lrg .and. matter_power_lnzsteps < 8) &
        call MpiStop('For combined LRGs and WiggleZ matter_power_lnzsteps should be set to at least 8 (hardcoded in cmbtypes)')

! wigglez redshifts: z0 = 0.d0, za = 0.22d0, zb = 0.41d0, zc = 0.6d0, zd = 0.78d0
! sdss lrg redshifts: z0 = 0.0d0, zNEAR = 0.235d0, zMID = 0.342d0, zFAR = 0.421d0   
!   if (use_dr7lrg) then
      iz0lrg = 1
      izwigglez(1) = 2
      izNEARlrg = 3
      izMIDlrg = 4
      izwigglez(2) = 5
      izFARlrg = 6
      izwigglez(3) = 7
      izwigglez(4) = 8
      redshifts(izNEARlrg) = zNEAR
      redshifts(izMIDlrg) = zMID
      redshifts(izFARlrg) = zFAR  
      redshifts(izwigglez(1)) = za
      redshifts(izwigglez(2)) = zb
      redshifts(izwigglez(3)) = zc
      redshifts(izwigglez(4)) = zd
!   endif
   redshifts(iz0lrg) = 0.0d0

   return
 end subroutine wigglez_sdss_SetTransferRedshifts
  
  subroutine ReadmpkDataset(gname)   
    use MatrixUtils
    character(LEN=*), intent(IN) :: gname
    character(LEN=Ini_max_string_len) :: kbands_file, measurements_file, windows_file, cov_file
    !! added for the LRG window function subtraction
    character(LEN=Ini_max_string_len) :: zerowindowfxn_file, zerowindowfxnsubtractdat_file

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

    real, dimension(:), allocatable :: mpk_zerowindowfxnfull
    real, dimension(:), allocatable :: mpk_zerowindowfxnsubfull

    character(80) :: dummychar
    logical bad
    Type(TIniFile) :: Ini
    integer file_unit
    
 
    num_mpk_datasets = num_mpk_datasets + 1
    if (num_mpk_datasets > 10) stop 'too many datasets'
    file_unit = new_file_unit()
    call Ini_Open_File(Ini, gname, file_unit, bad, .false.)
    if (bad) then
      write (*,*)  'Error opening dataset file '//trim(gname)
      stop
    end if

    mset%name = Ini_Read_String_File(Ini,'name') 
    Ini_fail_on_not_found = .false.
    mset%use_set =.true.
    if (Feedback > 0) write (*,*) 'reading: '//trim(mset%name)
    num_mpk_points_full = Ini_Read_Int_File(Ini,'num_mpk_points_full',0)
    if (num_mpk_points_full.eq.0) write(*,*) ' ERROR: parameter num_mpk_points_full not set'
    num_mpk_kbands_full = Ini_Read_Int_File(Ini,'num_mpk_kbands_full',0)
    if (num_mpk_kbands_full.eq.0) write(*,*) ' ERROR: parameter num_mpk_kbands_full not set'
    min_mpk_points_use = Ini_Read_Int_File(Ini,'min_mpk_points_use',1)
    min_mpk_kbands_use = Ini_Read_Int_File(Ini,'min_mpk_kbands_use',1)
    max_mpk_points_use = Ini_Read_Int_File(Ini,'max_mpk_points_use',num_mpk_points_full)
    max_mpk_kbands_use = Ini_Read_Int_File(Ini,'max_mpk_kbands_use',num_mpk_kbands_full)
    mset%num_mpk_points_use = max_mpk_points_use - min_mpk_points_use +1
    mset%num_mpk_kbands_use = max_mpk_kbands_use - min_mpk_kbands_use +1

    allocate(mpk_Wfull(num_mpk_points_full,num_mpk_kbands_full))
    allocate(mpk_kfull(num_mpk_kbands_full))
    allocate(mset%mpk_P(mset%num_mpk_points_use))
    allocate(mset%mpk_sdev(mset%num_mpk_points_use))  ! will need to replace with the covmat
    allocate(mset%mpk_k(mset%num_mpk_kbands_use))
    allocate(mset%mpk_W(mset%num_mpk_points_use,mset%num_mpk_kbands_use))
    allocate(mset%mpk_zerowindowfxn(mset%num_mpk_kbands_use))
    allocate(mset%mpk_zerowindowfxnsubtractdat(mset%num_mpk_points_use))
    allocate(mpk_fiducial(mset%num_mpk_points_use))
    allocate(mpk_zerowindowfxnsubfull(num_mpk_points_full+1)) 
      !!need to add 1 to get the normalization held in the first (really zeroth) entry
    allocate(mpk_zerowindowfxnfull(num_mpk_kbands_full))

    kbands_file  = ReadIniFileName(Ini,'kbands_file')
    call ReadVector(kbands_file,mpk_kfull,num_mpk_kbands_full)
    mset%mpk_k(1:mset%num_mpk_kbands_use)=mpk_kfull(min_mpk_kbands_use:max_mpk_kbands_use) 
    if (Feedback > 1) then 
       write(*,*) 'reading: ',mset%name,' data'
       write(*,*) 'Using kbands windows between',mset%mpk_k(1),' < k/h < ',mset%mpk_k(mset%num_mpk_kbands_use)      
    endif
    if  (mset%mpk_k(1) < matter_power_minkh) then
       write (*,*) 'WARNING: k_min in '//trim(mset%name)//'less than setting in cmbtypes.f90'
       write (*,*) 'all k<matter_power_minkh will be set to matter_power_minkh' 
    end if

    measurements_file  = ReadIniFileName(Ini,'measurements_file')
    call OpenTxtFile(measurements_file, tmp_file_unit)
    mset%mpk_P=0.
    read (tmp_file_unit,*) dummychar
    read (tmp_file_unit,*) dummychar
    do i= 1, (min_mpk_points_use-1)
       read (tmp_file_unit,*, iostat=iopb) keff,klo,khi,beff,beff,beff
    end do
    if (Feedback > 1 .and. min_mpk_points_use>1) write(*,*) 'Not using bands with keff=  ',keff,' or below'
    do i =1, mset%num_mpk_points_use
       read (tmp_file_unit,*, iostat=iopb) keff,klo,khi,mset%mpk_P(i),mset%mpk_sdev(i),mpk_fiducial(i)
    end do
    close(tmp_file_unit) 
    if (Feedback > 1) write(*,*) 'bands truncated at keff=  ',keff
    
    windows_file  = ReadIniFileName(Ini,'windows_file')
    if (windows_file.eq.'') write(*,*) 'ERROR: mpk windows_file not specified'
    call ReadMatrix(windows_file,mpk_Wfull,num_mpk_points_full,num_mpk_kbands_full)
    mset%mpk_W(1:mset%num_mpk_points_use,1:mset%num_mpk_kbands_use)= &
       mpk_Wfull(min_mpk_points_use:max_mpk_points_use,min_mpk_kbands_use:max_mpk_kbands_use)
    
 
    if (mset%name == 'lrg_2009') then
#ifndef DR71RG
        call MpiStop('mpk: edit makefile to have "EXTDATA = LRG" to inlude LRGs')
#else
        use_dr7lrg = .true.
        zerowindowfxn_file  = ReadIniFileName(Ini,'zerowindowfxn_file')
    
        print *, 'trying to read this many points', num_mpk_kbands_full
        if (zerowindowfxn_file.eq.'') write(*,*) 'ERROR: mpk zerowindowfxn_file not specified'
        call ReadVector(zerowindowfxn_file,mpk_zerowindowfxnfull,num_mpk_kbands_full)
        mset%mpk_zerowindowfxn(1:mset%num_mpk_kbands_use) = mpk_zerowindowfxnfull(min_mpk_kbands_use:max_mpk_kbands_use)
        zerowindowfxnsubtractdat_file  = ReadIniFileName(Ini,'zerowindowfxnsubtractdat_file')
        if (zerowindowfxnsubtractdat_file.eq.'') write(*,*) 'ERROR: mpk zerowindowfxnsubtractdat_file not specified'
        call ReadVector(zerowindowfxnsubtractdat_file,mpk_zerowindowfxnsubfull,num_mpk_points_full+1)
        mset%mpk_zerowindowfxnsubtractdat(1:mset%num_mpk_points_use) = &
         mpk_zerowindowfxnsubfull(min_mpk_points_use+1:max_mpk_points_use+1)
        mset%mpk_zerowindowfxnsubdatnorm = mpk_zerowindowfxnsubfull(1)
#endif
    end if

    cov_file  = ReadIniFileName(Ini,'cov_file')
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

    mset%use_scaling = Ini_Read_Logical_File(Ini,'use_scaling',.false.)

    mset%Q_marge = Ini_Read_Logical_File(Ini,'Q_marge',.false.)
    if (mset%Q_marge) then
     mset%Q_flat = Ini_Read_Logical_File(Ini,'Q_flat',.false.)
     if (.not. mset%Q_flat) then
      !gaussian prior on Q
      mset%Q_mid = Ini_Read_Real_File(Ini,'Q_mid')
      mset%Q_sigma = Ini_Read_Real_File(Ini,'Q_sigma')
     end if
     mset%Ag = Ini_Read_Real_File(Ini,'Ag', 1.4)
    end if 
    if (iopb.ne.0) then
       stop 'Error reading mpk file'
    endif
 
   call Ini_Close_File(Ini)
   call ClearFileUnit(file_unit)
   
   deallocate(mpk_Wfull, mpk_kfull,mpk_fiducial)

   mpkdatasets(num_mpk_datasets) = mset
 
  if (mset%name == 'lrg_2009') call LSS_LRG_mpklike_init()

  end subroutine ReadmpkDataset


  subroutine wigglez_ReadmpkDataset(gname)   
! this will be called once for each redshift bin
    use wigglezinfo
    use MatrixUtils
    implicit none
    character(LEN=*), intent(IN) :: gname
    character(LEN=Ini_max_string_len) :: kbands_file, measurements_file, windows_file, cov_file
    !! added for the LRG window function subtraction
    character(LEN=Ini_max_string_len) :: zerowindowfxn_file, zerowindowfxnsubtractdat_file

    Type (wigglez_mpkdataset) :: wmset

    integer i,iopb,i_regions
    real keff,klo,khi,beff
    integer :: num_mpk_points_full ! actual number of bandpowers in the infile
    integer :: num_mpk_kbands_full ! actual number of k positions " in the infile
    integer :: max_mpk_points_use ! in case you don't want the smallest scale modes (eg. sdss)
    integer :: min_mpk_points_use ! in case you don't want the largest scale modes
    integer :: max_mpk_kbands_use ! in case you don't want to calc P(k) on the smallest scales (will truncate P(k) to zero here!)
    integer :: min_mpk_kbands_use ! in case you don't want to calc P(k) on the largest scales (will truncate P(k) to zero here!)
    real, dimension(:,:,:), allocatable :: mpk_Wfull, mpk_covfull
    real, dimension(:), allocatable :: mpk_kfull  !, mpk_fiducial
    real, dimension(:,:), allocatable :: invcov_tmp
!    real, dimension(:), allocatable :: mpk_zerowindowfxnfull
!    real, dimension(:), allocatable :: mpk_zerowindowfxnsubfull
    character(len=64) region_string
    character(80) :: dummychar
    character z_char
    integer iz,count
    integer :: file_unit
    logical bad
    Type(TIniFile) :: Ini
    integer, parameter :: tmp_file_unit2=51
    
    iopb = 0
    num_wigglez_mpk_datasets = num_wigglez_mpk_datasets + 1
    if (num_wigglez_mpk_datasets > 10) stop 'too many datasets'
    file_unit = new_file_unit()
    call Ini_Open_File(Ini, gname, file_unit, bad, .false.)
    if (bad) then
      write (*,*)  'Error opening dataset file '//trim(gname)
      stop
    end if

#ifndef WIGZ
       call MpiStop('mpk: edit makefile to have "EXTDATA = WIGZ" to inlude WiggleZ data')
#else
       use_wigz10 = .true.
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
!    if(associated(wmset%mpk_P)) nullify(wmset%mpk_P)
    allocate(wmset%mpk_P(wmset%num_regions_used,wmset%num_mpk_points_use))
!    if(associated(wmset%mpk_k)) deallocate(wmset%mpk_k)
    allocate(wmset%mpk_k(wmset%num_mpk_kbands_use))
!    if(associated(wmset%mpk_W)) deallocate(wmset%mpk_W)
    allocate(wmset%mpk_W(wmset%num_regions_used,wmset%num_mpk_points_use,wmset%num_mpk_kbands_use))
!    allocate(mset%mpk_zerowindowfxn(mset%num_mpk_kbands_use))
!    allocate(mset%mpk_zerowindowfxnsubtractdat(mset%num_mpk_points_use))
!    allocate(mpk_fiducial(mset%num_mpk_points_use))
!    allocate(mpk_zerowindowfxnsubfull(num_mpk_points_full+1)) 
      !!need to add 1 to get the normalization held in the first (really zeroth) entry
!    allocate(mpk_zerowindowfxnfull(num_mpk_kbands_full))


    kbands_file  = ReadIniFileName(Ini,'kbands_file')
    call ReadVector(kbands_file,mpk_kfull,num_mpk_kbands_full)
    wmset%mpk_k(:)=mpk_kfull(min_mpk_kbands_use:max_mpk_kbands_use) 
    if (Feedback > 1) then 
       write(*,*) 'reading: ',wmset%name,' data'
       write(*,*) 'Using kbands windows between',wmset%mpk_k(1),' < k/h < ',wmset%mpk_k(wmset%num_mpk_kbands_use)      
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

          if (Feedback > 1 .and. min_mpk_points_use>1) write(*,*) 'Not using bands with keff=  ',keff,&
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
    if (Feedback > 1) write(*,*) 'bands truncated at keff=  ',keff

    allocate(mpk_Wfull(max_num_wigglez_regions,num_mpk_points_full,num_mpk_kbands_full))
    windows_file  = ReadIniFileName(Ini,'windows_file')
    if (windows_file.eq.'') write(*,*) 'ERROR: mpk windows_file not specified'
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
       stop 'Error reading mpk file'
    endif

    wmset%use_scaling = Ini_Read_Logical_File(Ini,'use_scaling',.false.)
    wmset%use_jennings = Ini_Read_Logical_File(Ini,'Use_jennings',.false.)
    wmset%use_simpledamping = Ini_Read_Logical_File(Ini,'Use_simpledamp',.false.)
    wmset%use_gigglez = Ini_Read_Logical_File(Ini,'Use_gigglez',.false.)

    if(wmset%use_jennings.and.wmset%use_simpledamping) then
         call MpiStop('mpk: cannot use both jennings formula and simple damping')
      endif

    if(wmset%use_gigglez.and.wmset%use_simpledamping) then
         call MpiStop('mpk: cannot use both jennings formula and gigglez scaling')
      endif

    if(wmset%use_jennings.and.wmset%use_gigglez) then
         call MpiStop('mpk: cannot use both gigglez scaling and simple damping')
      endif

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

    if(wmset%use_simpledamping) then ! need to read in sigma v
       wmset%damp_sigv = Ini_Read_Real_File(Ini,'sigmav')
    endif

    call Ini_Close_File(Ini)
    call ClearFileUnit(file_unit)


    wmpkdatasets(num_wigglez_mpk_datasets) = wmset

    if (Feedback > 1) write(*,*) 'read: ',wmpkdatasets(num_wigglez_mpk_datasets)%name,' data'
 
 end subroutine Wigglez_ReadmpkDataset
 
  function LSS_mpklike(Theory,mset,CMB) result(LnLike) ! LV_06 added CMB here
   Type (mpkdataset) :: mset
   Type (CosmoTheory) Theory
   Type(CMBparams) CMB     !LV_06 added for LRGDR4
   real LnLike
   real, dimension(:), allocatable :: mpk_Pth, mpk_k2,mpk_lin,k_scaled !LV_06 added for LRGDR4
   real, dimension(:), allocatable :: w
   real, dimension(:), allocatable :: mpk_WPth, mpk_WPth_k2
   real :: covdat(mset%num_mpk_points_use), covth(mset%num_mpk_points_use),  covth_k2(mset%num_mpk_points_use)
   real :: normV, Q, minchisq
   real :: a_scl  !LV_06 added for LRGDR4
   integer :: i, iQ
   logical :: do_marge
   integer, parameter :: nQ=6
   real :: tmp, dQ = 0.4
   real chisq(-nQ:nQ)
   real calweights(-nQ:nQ)
   real vec2(2),Mat(2,2)
   real(mpk_d) z, omk_fid,omv_fid, w_fid
   allocate(mpk_lin(mset%num_mpk_kbands_use) ,mpk_Pth(mset%num_mpk_kbands_use))
   allocate(mpk_WPth(mset%num_mpk_points_use))
   allocate(k_scaled(mset%num_mpk_kbands_use))!LV_06 added for LRGDR4
   allocate(w(mset%num_mpk_points_use))

   chisq = 0

   if (.not. mset%use_set) then
      LnLike = 0
      return
   end if

   ! won't actually want to do this multiple times for multiple galaxy pk data sets?..
   z = 0.
   omk_fid = 0.d0
   omv_fid = 0.75d0
   w_fid = -1.d0
   IF(mset%use_scaling) then   
      call compute_scaling_factor(z,dble(CMB%omk),dble(CMB%omv),dble(CMB%w),omk_fid,omv_fid,w_fid,a_scl)      
   else
     a_scl = 1
   end if
      

   do i=1, mset%num_mpk_kbands_use 
     !Errors from using matter_power_minkh at lower end should be negligible
         k_scaled(i)=max(matter_power_minkh,a_scl*mset%mpk_k(i))
         mpk_lin(i)=MatterPowerAt(Theory,k_scaled(i))/a_scl**3
   end do
      
      
    do_marge = mset%Q_Marge
    if (do_marge .and. mset%Q_flat) then
        !Marginalize analytically with flat prior on b^2 and b^2*Q
        !as recommended by Max Tegmark for SDSS
          allocate(mpk_k2(mset%num_mpk_kbands_use))
          allocate(mpk_WPth_k2(mset%num_mpk_points_use))
 
          mpk_Pth=mpk_lin/(1+mset%Ag*k_scaled)
          mpk_k2=mpk_Pth*k_scaled**2
          mpk_WPth = matmul(mset%mpk_W,mpk_Pth)
          mpk_WPth_k2 = matmul(mset%mpk_W,mpk_k2)

          if (associated(mset%mpk_invcov)) then
            covdat = matmul(mset%mpk_invcov,mset%mpk_P)
            covth = matmul(mset%mpk_invcov,mpk_WPth)
            covth_k2 = matmul(mset%mpk_invcov,mpk_WPth_k2)
          else
            w=1/(mset%mpk_sdev**2)
            covdat = mset%mpk_P*w
            covth = mpk_WPth*w
            covth_k2 = mpk_WPth_k2*w     
          end if
           
          Mat(1,1) = sum(covth*mpk_WPth)
          Mat(2,2) = sum(covth_k2*mpk_WPth_k2)
          Mat(1,2) = sum(covth*mpk_WPth_k2)
          Mat(2,1) = Mat(1,2)
          LnLike = log( Mat(1,1)*Mat(2,2)-Mat(1,2)**2)
          call inv_mat22(Mat)
          vec2(1) = sum(covdat*mpk_WPth)
          vec2(2) = sum(covdat*mpk_WPth_k2)
          LnLike = (sum(mset%mpk_P*covdat) - sum(vec2*matmul(Mat,vec2)) + LnLike ) /2

          deallocate(mpk_k2,mpk_WPth_k2)      
    else

      if (mset%Q_sigma==0) do_marge = .false.
      
      do iQ=-nQ,nQ
         Q = mset%Q_mid +iQ*mset%Q_sigma*dQ 
         
         if (mset%Q_marge) then
            mpk_Pth=mpk_lin*(1+Q*k_scaled**2)/(1+mset%Ag*k_scaled)
         else 
            mpk_Pth = mpk_lin
         end if
         
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

   end if !not analytic over Q
      
   if (Feedback>1) write(*,*) 'mpk chi-sq:', LnLike*2
   
   if (LnLike > 1e8) then
      write(*,*) 'Chisq is huge, maybe there is a problem? chisq=',chisq
   end if
   
   deallocate(mpk_Pth,mpk_lin)
   deallocate(mpk_WPth,k_scaled,w)
   
 end function LSS_mpklike


 function LSSLnLike(CMB, Theory)
   use wigglezinfo
   Type (CMBParams) CMB
   Type (CosmoTheory) Theory
   real LSSLnLike
   integer i
   real tot(num_mpk_datasets)
   real tot_wigglez(num_wigglez_mpk_datasets)

  do i=1, num_mpk_datasets
     if (mpkdatasets(i)%name == 'twodf') then
        stop 'twodf no longer supported - use data/2df_2005.dataset'
     else if (mpkdatasets(i)%name == 'lrg_2009') then
        tot(i) = LSS_LRG_mpklike(Theory,mpkdatasets(i),CMB)
     else
      tot(i) = LSS_mpklike(Theory,mpkdatasets(i),CMB) !LV_06 added CMB here
     end if
  end do

  if(use_wigz10) then
#ifndef WIGZ
     call MpiStop('mpk: edit makefile to have "EXTDATA = WIGZ" to inlude WiggleZ data')
#else
     do i=1, num_wigglez_mpk_datasets
        tot_wigglez(i) = WiggleZ_mpklike(Theory,wmpkdatasets(i),CMB) 
     end do
#endif
  endif

  LSSLnLike = SUM(tot)
  
#ifdef WIGZ
  LSSLnLike = LSSLnLike + SUM(tot_wigglez)
#endif  
 end function LSSLnLike

 subroutine inv_mat22(M)
    real M(2,2), Minv(2,2), det

    det = M(1,1)*M(2,2)-M(1,2)*M(2,1)
    Minv(1,1)=M(2,2)
    Minv(2,2) = M(1,1)
    Minv(1,2) = - M(2,1)
    Minv(2,1) = - M(1,2)
    M = Minv/det

 end subroutine inv_mat22

!-----------------------------------------------------------------------------
!LV added to include lrg DR4

subroutine compute_scaling_factor(z,Ok,Ol,w,Ok0,Ol0,w0,a)
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
  real a
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
end subroutine compute_scaling_factor

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
            


!! added by Beth Reid for LRG P(k) analysis

  function a2maxpos(a1val) result(a2max)
    real(dl), intent(in) :: a1val
    real(dl) a2max
    a2max = -1.0d0
    if (a1val <= min(s1/k1,s2/k2)) then
      a2max = min(s1/k1**2 - a1val/k1, s2/k2**2 - a1val/k2)
    end if
  end function a2maxpos

  function a2min1pos(a1val) result(a2min1)
    real(dl), intent(in) :: a1val
    real(dl) a2min1
    a2min1 = 0.0d0
    if(a1val <= 0.0d0) then
      a2min1 = max(-s1/k1**2 - a1val/k1, -s2/k2**2 - a1val/k2, 0.0d0)
    end if
  end function a2min1pos

  function a2min2pos(a1val) result(a2min2)
    real(dl), intent(in) :: a1val
    real(dl) a2min2
    a2min2 = 0.0d0
    if(abs(a1val) >= 2.0d0*s1/k1 .and. a1val <= 0.0d0)  then
      a2min2 = a1val**2/s1*0.25d0
    end if
  end function a2min2pos

  function a2min3pos(a1val) result(a2min3)
    real(dl), intent(in) :: a1val
    real(dl) a2min3
    a2min3 = 0.0d0
    if(abs(a1val) >= 2.0d0*s2/k2 .and. a1val <= 0.0d0)  then
      a2min3 = a1val**2/s2*0.25d0
    end if
  end function a2min3pos

  function a2minfinalpos(a1val) result(a2minpos)
    real(dl), intent(in) :: a1val
    real(dl) a2minpos
    a2minpos = max(a2min1pos(a1val),a2min2pos(a1val),a2min3pos(a1val))
  end function a2minfinalpos

  function a2minneg(a1val) result(a2min)
    real(dl), intent(in) :: a1val
    real(dl) a2min
    if (a1val >= max(-s1/k1,-s2/k2)) then
      a2min = max(-s1/k1**2 - a1val/k1, -s2/k2**2 - a1val/k2)
    else
      a2min = 1.0d0
    end if
  end function a2minneg

  function a2max1neg(a1val) result(a2max1)
    real(dl), intent(in) :: a1val
    real(dl) a2max1
    if(a1val >= 0.0d0) then
      a2max1 = min(s1/k1**2 - a1val/k1, s2/k2**2 - a1val/k2, 0.0d0)
    else
      a2max1 = 0.0d0
    end if
  end function a2max1neg

  function a2max2neg(a1val) result(a2max2)
    real(dl), intent(in) :: a1val
    real(dl) a2max2
    a2max2 = 0.0d0
    if(abs(a1val) >= 2.0d0*s1/k1 .and. a1val >= 0.0d0)  then
      a2max2 = -a1val**2/s1*0.25d0
    end if
  end function a2max2neg

  function a2max3neg(a1val) result(a2max3)
    real(dl), intent(in) :: a1val
    real(dl) a2max3
    a2max3 = 0.0d0
    if(abs(a1val) >= 2.0d0*s2/k2 .and. a1val >= 0.0d0)  then
      a2max3 = -a1val**2/s2*0.25d0
    end if
  end function a2max3neg

  function a2maxfinalneg(a1val) result(a2maxneg)
    real(dl), intent(in) :: a1val
    real(dl) a2maxneg
    a2maxneg = min(a2max1neg(a1val),a2max2neg(a1val),a2max3neg(a1val))
  end function a2maxfinalneg


function testa1a2(a1val, a2val) result(testresult)
    real(dl), intent(in) :: a1val,a2val
    logical :: testresult

    real(dl) :: kext, diffval
    testresult = .true.

    ! check if there's an extremum; either a1val or a2val has to be negative, not both
    kext = -a1val/2.0d0/a2val
    diffval = abs(a1val*kext + a2val*kext**2)
    if(kext > 0.0d0 .and. kext <= k1 .and. diffval > s1) testresult = .false.
    if(kext > 0.0d0 .and. kext <= k2 .and. diffval > s2) testresult = .false.

    if (abs(a1val*k1 + a2val*k1**2) > s1) testresult = .false.
    if (abs(a1val*k2 + a2val*k2**2) > s2) testresult = .false.

end function testa1a2

!! copying LSS_mpklike above.
!! points_use is how many points to use in the likelihood calculation; 
!!kbands_use is how many points you need to have a theory for in order to convolve the theory with the window function.


! this subroutine fills in the a1 and a2 values only once.
 subroutine LSS_LRG_mpklike_init()
   real(dl) :: a1val, a2val
   real(dl) :: da1, da2  ! spacing of numerical integral over nuisance params.
   integer :: countcheck = 0
   integer :: i, j
   !! this is just for checking the 'theory' curves for fiducial model
   real :: fidLnLike
   type(CosmoTheory) :: temptheory
   type(CMBparams) :: tempCMB

   da1 = a1maxval/(nptsa1/2)
   da2 = a2maxpos(-a1maxval)/(nptsa2/2)
   do i = -nptsa1/2, nptsa1/2
      do j = -nptsa2/2, nptsa2/2
         a1val = da1*i
         a2val = da2*j

         if ((a2val >= 0.0d0 .and. a2val <= a2maxpos(a1val) .and. a2val >= a2minfinalpos(a1val)) .or. &
     & (a2val <= 0.0d0 .and. a2val <= a2maxfinalneg(a1val) .and. a2val >= a2minneg(a1val))) then
            if(testa1a2(a1val,a2val) .eqv. .false.)  then
               print *,'Failed a1, a2: ',a1val,a2val
               if (a2val >= 0.0d0) print *,'pos', a2maxpos(a1val), a2minfinalpos(a1val)
               if (a2val <= 0.0d0) print *,'neg', a2maxfinalneg(a1val), a2minneg(a1val)
               stop
            end if
            countcheck = countcheck + 1
            if(countcheck > nptstot) then
               print *, 'countcheck > nptstot failure.'
               stop
            end if
            a1list(countcheck) = a1val
            a2list(countcheck) = a2val
            !print *, countcheck, a1list(countcheck), a2list(countcheck)
         end if
      end do
   end do
   if(countcheck .ne. nptstot) then
     print *, 'countcheck issue', countcheck, nptstot
     stop
   end if
  call LRGinfo_init()
 end subroutine LSS_LRG_mpklike_init


 function LSS_LRG_mpklike(Theory,mset,CMB) result(LnLike)  ! LV_06 added CMB here
   Type (mpkdataset) :: mset
   Type (CosmoTheory) Theory
   Type(CMBparams) CMB                  !LV_06 added for LRGDR4
   real LnLike
   integer :: i
   real, dimension(:), allocatable :: mpk_raw, mpk_Pth, mpk_Pth_k, mpk_Pth_k2, k_scaled
   real, dimension(:), allocatable :: mpk_WPth, mpk_WPth_k, mpk_WPth_k2
   real :: covdat(mset%num_mpk_points_use), covth(mset%num_mpk_points_use), &
          & covth_k(mset%num_mpk_points_use), covth_k2(mset%num_mpk_points_use), & 
          & covth_zerowin(mset%num_mpk_points_use)

   real, dimension(nptstot) :: chisq, chisqmarg  !! minus log likelihood list
   real :: minchisq,maxchisq,deltaL

   real(dl) :: a1val, a2val, zerowinsub
   real :: sumDD, sumDT, sumDT_k, sumDT_k2, sumTT,& 
     &  sumTT_k, sumTT_k2, sumTT_k_k, sumTT_k_k2, sumTT_k2_k2, &
     &  sumDT_tot, sumTT_tot, &
     &  sumDT_zerowin, sumTT_zerowin, sumTT_k_zerowin, sumTT_k2_zerowin, sumTT_zerowin_zerowin

   real :: sumzerow_Pth, sumzerow_Pth_k, sumzerow_Pth_k2
   real(mpk_d) :: z, omk_fid, omv_fid, w_fid
   real :: a_scl      !LV_06 added for LRGDR4

   real(wp) :: temp1,temp2,temp3
   real :: temp4

   !! added for no marg
   integer :: myminchisqindx
   real :: currminchisq, currminchisqmarg, minchisqtheoryamp, chisqnonuis
   real :: minchisqtheoryampnonuis, minchisqtheoryampminnuis
   real(dl), dimension(2) :: myerfval

   call fill_LRGTheory(Theory,matter_power_minkh,matter_power_dlnkh)
   allocate(mpk_raw(mset%num_mpk_kbands_use) ,mpk_Pth(mset%num_mpk_kbands_use))
   allocate(mpk_Pth_k(mset%num_mpk_kbands_use) ,mpk_Pth_k2(mset%num_mpk_kbands_use))
   allocate(mpk_WPth(mset%num_mpk_points_use),mpk_WPth_k(mset%num_mpk_points_use),mpk_WPth_k2(mset%num_mpk_points_use))
   allocate(k_scaled(mset%num_mpk_kbands_use))!LV_06 added for LRGDR4

   chisq = 0

   if (.not. mset%use_set) then
      LnLike = 0
      return
   end if
   z = zeffDR7
   omk_fid = 0.d0
   omv_fid = 0.75d0
   w_fid = -1.d0
   IF(mset%use_scaling) then
!      call compute_scaling_factor(z,dble(CMB%omk),dble(CMB%omv),dble(CMB%w),a_scl)
      call compute_scaling_factor(z,dble(CMB%omk),dble(CMB%omv),dble(CMB%w),omk_fid,omv_fid,w_fid,a_scl) 
      !! this step now applied in compute_scaling_factor
      !! this fixes the bug most easily !!
      !!a_scl = 1.0d0/a_scl
   else
     a_scl = 1
     stop 'use_scaling should be set to true for the LRGs!'
   end if

   do i=1, mset%num_mpk_kbands_use
         k_scaled(i)=max(matter_power_minkh,a_scl*mset%mpk_k(i))
         mpk_raw(i)=LRGPowerAt(Theory,k_scaled(i))/a_scl**3
   end do

   mpk_Pth = mpk_raw

   mpk_Pth_k = mpk_Pth*k_scaled
   mpk_Pth_k2 = mpk_Pth*k_scaled**2
   mpk_WPth = matmul(mset%mpk_W,mpk_Pth)
   mpk_WPth_k = matmul(mset%mpk_W,mpk_Pth_k)
   mpk_WPth_k2 = matmul(mset%mpk_W,mpk_Pth_k2)

   sumzerow_Pth = sum(mset%mpk_zerowindowfxn*mpk_Pth)/mset%mpk_zerowindowfxnsubdatnorm
   sumzerow_Pth_k = sum(mset%mpk_zerowindowfxn*mpk_Pth_k)/mset%mpk_zerowindowfxnsubdatnorm
   sumzerow_Pth_k2 = sum(mset%mpk_zerowindowfxn*mpk_Pth_k2)/mset%mpk_zerowindowfxnsubdatnorm


   covdat = matmul(mset%mpk_invcov,mset%mpk_P)
   covth = matmul(mset%mpk_invcov,mpk_WPth)
   covth_k = matmul(mset%mpk_invcov,mpk_WPth_k)
   covth_k2 = matmul(mset%mpk_invcov,mpk_WPth_k2)
   covth_zerowin = matmul(mset%mpk_invcov,mset%mpk_zerowindowfxnsubtractdat)

   sumDD = sum(mset%mpk_P*covdat)
   sumDT = sum(mset%mpk_P*covth)
   sumDT_k = sum(mset%mpk_P*covth_k)
   sumDT_k2 = sum(mset%mpk_P*covth_k2)
   sumDT_zerowin = sum(mset%mpk_P*covth_zerowin)

   sumTT = sum(mpk_WPth*covth)
   sumTT_k = sum(mpk_WPth*covth_k)
   sumTT_k2 = sum(mpk_WPth*covth_k2)
   sumTT_k_k = sum(mpk_WPth_k*covth_k)
   sumTT_k_k2 = sum(mpk_WPth_k*covth_k2)
   sumTT_k2_k2 = sum(mpk_WPth_k2*covth_k2)
   sumTT_zerowin = sum(mpk_WPth*covth_zerowin)
   sumTT_k_zerowin = sum(mpk_WPth_k*covth_zerowin)
   sumTT_k2_zerowin = sum(mpk_WPth_k2*covth_zerowin)
   sumTT_zerowin_zerowin = sum(mset%mpk_zerowindowfxnsubtractdat*covth_zerowin)

   currminchisq = 1000.0d0
   do i=1,nptstot
     a1val = a1list(i)
     a2val = a2list(i)
     zerowinsub = -(sumzerow_Pth + a1val*sumzerow_Pth_k + a2val*sumzerow_Pth_k2)

     sumDT_tot = sumDT + a1val*sumDT_k + a2val*sumDT_k2 + zerowinsub*sumDT_zerowin
     sumTT_tot = sumTT + a1val**2.0d0*sumTT_k_k + a2val**2.0d0*sumTT_k2_k2 + & 
                 & zerowinsub**2.0d0*sumTT_zerowin_zerowin &
       & + 2.0d0*a1val*sumTT_k + 2.0d0*a2val*sumTT_k2 + 2.0d0*a1val*a2val*sumTT_k_k2 &
       & + 2.0d0*zerowinsub*sumTT_zerowin + 2.0d0*zerowinsub*a1val*sumTT_k_zerowin &
       & + 2.0d0*zerowinsub*a2val*sumTT_k2_zerowin
     minchisqtheoryamp = sumDT_tot/sumTT_tot
     chisq(i) = sumDD - 2.0d0*minchisqtheoryamp*sumDT_tot + minchisqtheoryamp**2.0d0*sumTT_tot
#ifdef DR71RG
     myerfval(1) = sumDT_tot/2.0d0/sqrt(sumTT_tot)
     call geterf(myerfval)
     chisqmarg(i) = sumDD - sumDT_tot**2.0d0/sumTT_tot & 
         & + log(sumTT_tot) &
         & - 2.0*log(1.0d0 + myerfval(2))
#else
     !!leave out the erf term, just to get it to compile.  This should never run.
     chisqmarg(i) = sumDD - sumDT_tot**2.0d0/sumTT_tot & 
         & + log(sumTT_tot)
     if(0 .eq. 0) stop 'Logic problem.  Shouldnt be here.'
#endif
!this should always be here, but we're using gsl to call erf, so this function is only available if gsl is installed.
     if(i == 1 .or. chisq(i) < currminchisq) then
        myminchisqindx = i
        currminchisq = chisq(i)
        currminchisqmarg = chisqmarg(i)
        minchisqtheoryampminnuis = minchisqtheoryamp
     end if
     if(i == int(nptstot/2)+1) then
        chisqnonuis = chisq(i)
        minchisqtheoryampnonuis = minchisqtheoryamp
        if(abs(a1val) > 0.001 .or. abs(a2val) > 0.001) then
           print *, 'ahhhh! violation!!', a1val, a2val
        end if
     end if

   end do

! numerically marginalize over a1,a2 now using values stored in chisq
   minchisq = minval(chisqmarg)
   maxchisq = maxval(chisqmarg)

   LnLike = sum(exp(-(chisqmarg-minchisq)/2.0d0)/(nptstot*1.0d0))
   if(LnLike == 0) then
     LnLike = LogZero
   else
     LnLike = -log(LnLike) + minchisq/2.0d0
   end if
 deltaL = (maxchisq - minchisq)*0.5
 if(Feedback > 1) print *,'LRG P(k) LnLike = ',LnLike

   deallocate(mpk_raw, mpk_Pth)
   deallocate(mpk_Pth_k, mpk_Pth_k2)
   deallocate(mpk_WPth, mpk_WPth_k, mpk_WPth_k2)
   deallocate(k_scaled) 

 end function LSS_LRG_mpklike

  function WiggleZ_mpklike(Theory,wmset,CMB) result(LnLike) 
    use wigglezinfo
   Type (wigglez_mpkdataset) :: wmset
   Type (CosmoTheory) Theory
   Type(CMBparams) CMB     
   real LnLike
   real, dimension(:), allocatable :: mpk_Pth, mpk_k2,mpk_lin,k_scaled !LV_06 added for LRGDR4
!   real, dimension(:), allocatable :: w
   real, dimension(:), allocatable :: mpk_WPth, mpk_WPth_k2
   real, dimension(:), allocatable :: diffs, step1
   real, dimension(:), allocatable :: Pk_delta_delta,Pk_delta_theta,Pk_theta_theta
   real, dimension(:), allocatable :: damp1, damp2, damp3
   real :: covdat(wmset%num_mpk_points_use)
   real :: covth(wmset%num_mpk_points_use)
   real :: covth_k2(wmset%num_mpk_points_use)
   real, dimension(:), allocatable :: mpk_WPth_large, covdat_large, covth_large, mpk_Pdata_large
   integer imin,imax
   real :: normV, Q, minchisq
   real :: a_scl  !LV_06 added for LRGDR4
   integer :: i, iQ,ibias,ik,j,iz
   logical :: do_marge
   integer, parameter :: nQ=6
   integer, parameter :: nbias = 100
   real b0, bias_max, bias_step,bias,old_chisq,beta_val,kval,xi
   real :: tmp, dQ = 0.4
   real, dimension(:), allocatable :: chisq(:)
   real calweights(-nQ:nQ)
   real vec2(2),Mat(2,2)
   real final_term, b_out
   real(dl) z,omk_fid, omv_fid,w_fid
   integer i_region
   character(len=32) fname

   If(Feedback > 1) print*, 'Calling WiggleZ likelihood routines'
   allocate(mpk_lin(wmset%num_mpk_kbands_use),mpk_Pth(wmset%num_mpk_kbands_use))
   allocate(mpk_WPth(wmset%num_mpk_points_use))
   allocate(k_scaled(wmset%num_mpk_kbands_use))!LV_06 added for LRGDR4 !! IMPORTANT: need to check k-scaling
!   allocate(wmset%num_mpk_points_use))
   allocate(diffs(wmset%num_mpk_points_use),step1(wmset%num_mpk_points_use))
   allocate(Pk_delta_delta(wmset%num_mpk_kbands_use),Pk_delta_theta(wmset%num_mpk_kbands_use))
   allocate(Pk_theta_theta(wmset%num_mpk_kbands_use))
   allocate(damp1(wmset%num_mpk_kbands_use),damp2(wmset%num_mpk_kbands_use),damp3(wmset%num_mpk_kbands_use))

   if(wmset%use_jennings.or.wmset%use_simpledamping) then
      allocate(chisq(nbias))
      b0 = 0.0
      bias_max = 2.0
      bias_step = (bias_max-b0)/real(nbias)
   else 
      allocate(chisq(-nQ:nQ))
   endif
   chisq = 0

   if (.not. wmset%use_set) then
      LnLike = 0
      return
   end if

   ! won't actually want to do this multiple times for multiple galaxy pk data sets?..
   omk_fid = 0.d0
   omv_fid = 0.705d0
   w_fid = -1.d0
   z = 1.d0*dble(wmset%redshift) ! accuracy issues
   IF(wmset%use_scaling) then
!      call compute_scaling_factor(dble(z),dble(CMB%omk),dble(CMB%omv),dble(CMB%w),a_scl)
      call compute_scaling_factor(z,dble(CMB%omk),dble(CMB%omv),dble(CMB%w),omk_fid,omv_fid,w_fid,a_scl)
   else
     a_scl = 1
   end if

   iz = 0
   do i=1,4
      if(abs(z-zeval(i)).le.0.001) iz = i
   enddo
   if(iz.eq.0) call MpiStop('could not indentify redshift')

   if(wmset%use_gigglez) then
!      sigma_v = sqrt(func_sigma_v(num_matter_power,matter_power_minkh,matter_power_dlnkh,real(z*1.),&
!           CMB%omk,CMB%omv,Theory%matter_power(:,izwigglez(iz))))
!      if(feedback.ge.2) print*, 'redshift', z, 'sigma v', sigma_v
      call fill_GiggleZTheory(Theory,matter_power_minkh,matter_power_dlnkh,z)
!     call fill_GiggleZTheory(Theory,matter_power_minkh,matter_power_dlnkh,z,CMB%h,sigma_v)
   endif

   do i=1, wmset%num_mpk_kbands_use 
      ! It could be that when we scale the k-values, the lowest bin drops off the bottom edge
     !Errors from using matter_power_minkh at lower end should be negligible
         k_scaled(i)=max(matter_power_minkh,wmset%mpk_k(i)*a_scl)
         if(wmset%use_gigglez) then
            mpk_lin(i) = LRGPowerAt(Theory,k_scaled(i))/a_scl**3
         else
            mpk_lin(i)=MatterPowerAt_zbin(Theory,k_scaled(i),izwigglez(iz))/a_scl**3
         endif
   end do
   
    do_marge = wmset%Q_Marge
    if (do_marge .and. wmset%Q_flat) then
       !Marginalize analytically with flat prior on b^2 and b^2*Q
       !as recommended by Max Tegmark for SDSS
       allocate(mpk_k2(wmset%num_mpk_kbands_use))
       allocate(mpk_WPth_k2(wmset%num_mpk_points_use))
       
       Mat(:,:) = 0.d0
       vec2(:) = 0.d0
       final_term = 0.d0
       do i_region=1,wmset%num_regions_used
          mpk_Pth(:)=mpk_lin(:)/(1+wmset%Ag*k_scaled)
          mpk_k2(:)=mpk_Pth(:)*k_scaled(:)**2
          
          
          mpk_WPth(:) = matmul(wmset%mpk_W(i_region,:,:),mpk_Pth(:))
          mpk_WPth_k2(:) = matmul(wmset%mpk_W(i_region,:,:),mpk_k2(:))
          
          covdat(:) = matmul(wmset%mpk_invcov(i_region,:,:),wmset%mpk_P(i_region,:))
          covth(:) = matmul(wmset%mpk_invcov(i_region,:,:),mpk_WPth(:))
          covth_k2(:) = matmul(wmset%mpk_invcov(i_region,:,:),mpk_WPth_k2(:))
          
          Mat(1,1) = Mat(1,1) + sum(covth(:)*mpk_WPth(:))
          Mat(2,2) = Mat(2,2) + sum(covth_k2(:)*mpk_WPth_k2(:))
          Mat(1,2) = Mat(1,2) + sum(covth(:)*mpk_WPth_k2(:))
          Mat(2,1) = Mat(1,2)

          vec2(1) = vec2(1) + sum(covdat(:)*mpk_WPth(:))
          vec2(2) = vec2(2) + sum(covdat(:)*mpk_WPth_k2(:))
          final_term = final_term + sum(wmset%mpk_P(i_region,:)*covdat(:))
       enddo
       LnLike = log( Mat(1,1)*Mat(2,2)-Mat(1,2)**2)
       call inv_mat22(Mat)
       !          LnLike = (sum(mset%mpk_P*covdat) - sum(vec2*matmul(Mat,vec2)) + LnLike ) /2
       LnLike = (final_term - sum(vec2*matmul(Mat,vec2)) + LnLike ) /2

       deallocate(mpk_k2,mpk_WPth_k2)      
    else
       if (wmset%Q_sigma==0) do_marge = .false.
       ! ... sum the chi-squared contributions for all regions first
       chisq(:) = 0.d0
       old_chisq = 1.d30
       if(wmset%use_jennings.or.wmset%use_simpledamping) then
          if(feedback > 1) print*, "starting direct marginalisation over bias"
          do_marge = .false.
          ! ... need to numerically marginalise over bias ...
          mpk_Pth(:) = mpk_lin(:)
          if(wmset%use_jennings) then
             call angle_averaged_pk(wmset%num_mpk_kbands_use,k_scaled,mpk_Pth,wmset%redshift,CMB%omk,CMB%omv&
                  ,Pk_delta_delta,Pk_delta_theta,Pk_theta_theta,Theory)
          else if(wmset%use_simpledamping) then
             call compute_damping_terms(wmset%num_mpk_kbands_use,k_scaled,mpk_Pth,wmset%redshift,wmset%damp_sigv,&
                  CMB%omk,CMB%omv,damp1,damp2,damp3)
          endif
          do ibias = 1,nbias
             bias = b0+real(ibias)*bias_step
             if(wmset%use_jennings) then
                mpk_Pth(:) = bias**2*Pk_delta_delta(:) + bias*Pk_delta_theta(:) + Pk_theta_theta(:) 
             else if(wmset%use_simpledamping) then
                beta_val = growth_rate(real(z),real(CMB%omk),real(CMB%omv))/bias
                mpk_Pth(:) = mpk_lin(:)*(damp1(:) + 2.d0*beta_val*damp2(:) + beta_val**2+damp3(:))
             else
                mpk_Pth(:) = bias**2*mpk_lin(:)
             endif
!             if(ibias.eq.12) then
!                write(fname,'(a,i2.2,a)') 'mpk_damped_',int(wmset%redshift*100),'.dat'
!                open(unit=21,file=TRIM(fname),status='unknown')
             !               write(fname,'(a,i2.2,a)') 'mpk_jennings_',int(wmset%redshift*100),'.dat'
             !               open(unit=22,file=TRIM(fname),status='unknown')
             
!                do i=1,wmset%num_mpk_kbands_use
!                   write(21,*) k_scaled(i), mpk_Pth(i), damp1(i), 2.d0*beta_val*damp2(i), beta_val**2*damp3(i)
             !                  write(21,*) k_scaled(i), mpk_lin(i) 
             !                  write(22,*) k_scaled(i),  Pk_delta_delta(i) + Pk_delta_theta(i) + Pk_theta_theta(i) 
!                enddo
!                close(21)
             !               close(22)
!             endif
             
             do i_region=1,wmset%num_regions_used
                mpk_WPth(:) = matmul(wmset%mpk_W(i_region,:,:),mpk_Pth(:))
                diffs(:) = wmset%mpk_P(i_region,:)-mpk_WPth(:)
                step1(:) = matmul(wmset%mpk_invcov(i_region,:,:),diffs(:))
                chisq(ibias) = chisq(ibias) +  dot_product(diffs,step1)
                if(chisq(ibias).lt.old_chisq) then
                   old_chisq = chisq(ibias)
                endif
             enddo
          enddo
          minchisq = minval(chisq)
          LnLike = sum(exp(-(chisq-minchisq)/2))/real(nbias)
          if (LnLike == 0) then
             LnLike = LogZero
          else
             LnLike =  -log(LnLike) + minchisq/2
          end if
          deallocate(diffs,step1)
          deallocate(Pk_delta_delta,Pk_delta_theta,Pk_theta_theta)
       else
         if(feedback > 1) print*, "starting analytic marginalisation over bias"
         allocate(mpk_Pdata_large(wmset%num_mpk_points_use*wmset%num_regions_used))
         allocate(mpk_WPth_large(wmset%num_mpk_points_use*wmset%num_regions_used))
         allocate(covdat_large(wmset%num_mpk_points_use*wmset%num_regions_used))       
         allocate(covth_large(wmset%num_mpk_points_use*wmset%num_regions_used))
         normV = 0.d0
         do iQ=-nQ,nQ
            Q = wmset%Q_mid +iQ*wmset%Q_sigma*dQ 
            if (wmset%Q_marge) then
               mpk_Pth(:)=mpk_lin(:)*(1+Q*k_scaled(:)**2)/(1+wmset%Ag*k_scaled(:))
            else 
               mpk_Pth(:) = mpk_lin(:)
            end if
            do i_region=1,wmset%num_regions_used
               imin = (i_region-1)*wmset%num_mpk_points_use+1
               imax = i_region*wmset%num_mpk_points_use
               mpk_WPth(:) = matmul(wmset%mpk_W(i_region,:,:),mpk_Pth(:))
               mpk_Pdata_large(imin:imax) = wmset%mpk_P(i_region,:)
               mpk_WPth_large(imin:imax) = mpk_WPth(:) 
               
               !with analytic marginalization over normalization nuisance (flat prior on b^2)
               !See appendix F of cosmomc paper
               
               !         if (associated(mset%mpk_invcov)) then
               covdat_large(imin:imax) = matmul(wmset%mpk_invcov(i_region,:,:),wmset%mpk_P(i_region,:))
               covth_large(imin:imax) = matmul(wmset%mpk_invcov(i_region,:,:),mpk_WPth(:))
            enddo
            normV = normV + sum(mpk_WPth_large*covth_large)
            b_out =  sum(mpk_WPth_large*covdat_large)/sum(mpk_WPth_large*covth_large)
            if(Feedback.ge.2) print*, "Bias value:", b_out
            chisq(iQ) = sum(mpk_Pdata_large*covdat_large)  - sum(mpk_WPth_large*covdat_large)**2/normV!  + log(normV)
            !         else
            
            !            w=1/(mset%mpk_sdev**2)
            !            normV = sum(mpk_WPth*mpk_WPth*w)
            !            tmp=sum(mpk_WPth*mset%mpk_P*w)/normV ! avoid subtracting one large number from another
            !            chisq(iQ) = sum(mset%mpk_P*(mset%mpk_P - mpk_WPth*tmp)*w)  + log(normV)
            !         end if
            
            if (do_marge) then
               calweights(iQ) = exp(-(iQ*dQ)**2/2)
            else 
               LnLike = chisq(iQ)/2
               exit
            end if
            
         end do
         deallocate(covdat_large,covth_large,mpk_Pdata_large,mpk_WPth_large)
      endif
      !without analytic marginalization
      !! chisq = sum((mset%mpk_P(:) - mpk_WPth(:))**2*w) ! uncommented for debugging purposes
       if (do_marge) then
          if(.not.wmset%use_jennings) then
             minchisq=minval(chisq)
             LnLike = sum(exp(-(chisq-minchisq)/2)*calweights)/sum(calweights)
             if (LnLike == 0) then
                LnLike = LogZero
             else
                LnLike =  -log(LnLike) + minchisq/2
             end if
          endif
       end if

    end if !not analytic over Q
      
   if (Feedback>1) write(*,*) 'mpk chi-sq:', LnLike*2
   
   if (LnLike > 1e8) then
      write(*,*) 'Chisq is huge, maybe there is a problem? chisq=',chisq
   end if
   
   deallocate(mpk_Pth,mpk_lin)
   deallocate(mpk_WPth,k_scaled)!,w)
   deallocate(chisq)
   
 end function WiggleZ_mpklike


 subroutine ReadWiggleZMatrices(aname,mat,num_regions,m,n)
! suborutine to read all the matrices from each of the different regions, enclosed in one file

   implicit none
   character(LEN=*), intent(IN) :: aname
   integer, intent(in) :: m,n,num_regions
   real, intent(out) :: mat(num_regions,m,n)
   integer j,k,i_region
   real tmp
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

! Include effects of angle averaged power spectrum
! and velocity damping
! assuming 

 subroutine angle_averaged_pk(nkbands,k,Pk_input,z,omega_k,omega_lambda,Pk_dd_z,Pk_dt_z,Pk_tt_z,Theory)
! only works with linear power spectrum and linear bias
! does integral over mu
   use wigglezinfo
   implicit none
   integer, intent(in) :: nkbands
   real, intent(in) :: z
   real, intent(in), dimension(nkbands) :: k
   real, intent(in), dimension(nkbands) :: Pk_input ! density Pk at redshift z
   real, intent(out), dimension(nkbands) :: Pk_dd_z, Pk_dt_z, Pk_tt_z
   Type(CosmoTheory), intent(in) :: Theory
   real, dimension(nkbands) :: Pk_linear_0
   real omega_k, omega_lambda, growth_0, growth_eval2, growth2_0
   integer ik, iz,i 
   real, dimension(:), allocatable :: k_dummy, P_dd_dummy, P_lin_dummy
! coefficients for angle averaged power spectra
   real, dimension(nkbands) :: A0, A2, A4
! density-velocity and velocity velcoity power spectra
   real, dimension(nkbands) :: P_delta_theta, P_theta_theta
! coffeicients for delta_theta and theta_theta power spectra
! from Jennings et al 2010
   real, dimension(0:3) :: alpha_delta_theta, alpha_theta_theta
! mimimum and maximum mu
   real(dl) :: mu_min, mu_max
   real  scaling_coeff
   real(dl), parameter :: eps=1.d-4
   real(dl) rombint
   external rombint
   character(len=27) fname
   if(feedback.ge.2) print*, 'computing angle averaged pk'

   do ik=1,nkbands
      Pk_linear_0(ik) = MatterPowerAt(Theory,k(ik)) ! linear/halofit pk at z=0
   enddo 

! use fitting formula to estimate density-velocity diveregence and velocity dv-velocity divergence power spec
   alpha_delta_theta(0) = -12288.7
   alpha_delta_theta(1) = 1.43
   alpha_delta_theta(2) = 1367.7
   alpha_delta_theta(3) = 1.54
   P_delta_theta(:) = ((alpha_delta_theta(0)*sqrt(Pk_linear_0(:)))+(alpha_delta_theta(1)*(Pk_linear_0(:)**2)))&
        /(alpha_delta_theta(2)+(alpha_delta_theta(3)*Pk_linear_0(:)))
   alpha_theta_theta(0) = -12462.1
   alpha_theta_theta(1) = 0.839
   alpha_theta_theta(2) = 1446.6
   alpha_theta_theta(3) = 0.806
   P_theta_theta(:) = ((alpha_theta_theta(0)*sqrt(Pk_linear_0(:)))+(alpha_theta_theta(1)*(Pk_linear_0(:)**2)))&
        /(alpha_theta_theta(2)+(alpha_theta_theta(3)*Pk_linear_0(:)))

!   do ik=1,nkbands
!      write(90,*) k(ik), Pk_linear(ik), P_delta_theta(ik), P_theta_theta(ik)
!   enddo


! need to scale P_delta_theta and P_theta_theta to new redshift

   growth_eval = growth(z,omega_k,omega_lambda)
   growth_eval2 = growth2(z,omega_k,omega_lambda)
   growth_0 = 1.d0 !growth(0.0)
!   growth2_0 = growth2(0.0,omega_k,omega_lambda)
!   scaling_coeff = (growth_eval+growth_eval**2+growth_eval**3)/(growth_0+growth_0**2+growth_0**3)
   scaling_coeff = (growth_0+growth_0**2+growth_0**3)/(growth_eval+growth_eval**2+growth_eval**3)
   if(feedback.ge.2) print*, 'scaling coefficient', scaling_coeff
   P_delta_theta(:) = (P_delta_theta(:) -Pk_linear_0(:))/(scaling_coeff**2) + Pk_input(:)
   P_theta_theta(:) = (P_theta_theta(:) -Pk_linear_0(:))/(scaling_coeff**2) + Pk_input(:)


!   write(fname,'(a,i2.2,a)') 'P_delta_theta_',int(z*100),'.dat'
!   open(unit=25,file=TRIM(fname),status='unknown')
!   write(fname,'(a,i2.2,a)') 'P_theta_theta_',int(z*100),'.dat'
!   open(unit=26,file=TRIM(fname),status='unknown')
!   do ik=1,nkbands
!      write(25,*) k(ik), P_delta_theta(ik)
!      write(26,*) k(ik), P_theta_theta(ik)
!   enddo
!   close(25)
!   close(26)
! now compute sigma_v
   
   iz = 0
   do i=1,5
      if(abs(z-zeval(i)).le.0.001) iz = i
   enddo

   sigma_v = sqrt(func_sigma_v(num_matter_power,matter_power_minkh,matter_power_dlnkh,real(z*1.),&
           omega_k,omega_lambda,Theory%matter_power(:,iz)))
   if(feedback.ge.2) print*, 'redshift', z, 'sigma v', sigma_v, 'growth eval', growth_eval, 'growth eval2', growth_eval2 

! now compute A coefficients



   mu_min = -1.
   mu_max = 1.
   do ik=1,nkbands
      k_eval = k(ik)
      growth_rate_eval = growth_rate(z,omega_k,omega_lambda)
      A0(ik) = real(rombint(func_A0,mu_min,mu_max,eps))
      A2(ik) = real(rombint(func_A2,mu_min,mu_max,eps))
      A4(ik) = real(rombint(func_A4,mu_min,mu_max,eps))
   enddo
   A0(:) = A0(:)/2.0
   A2(:) = A2(:)/2.0
   A4(:) = A4(:)/2.0

!   do ik=1,nkbands
!      write(91,*) k(ik), A0(ik), A2(ik), A4(ik)
!   enddo
! fold the A coffeficents and growth rate into the power spec before returning
!   do ik=1,nkbands
!      write(92,*) k(ik), Pk_input(ik), P_delta_theta(ik), P_theta_theta(ik)
!   enddo


   Pk_dd_z(:) = Pk_input(:)*A0(:)
   Pk_dt_z(:) = growth_rate(z,omega_k,omega_lambda)*A2(:)*P_delta_theta(:)
   Pk_tt_z(:) = growth_rate(z,omega_k,omega_lambda)**2*A4(:)*P_theta_theta(:)


! returns angle-averaged Pks at redshift z
   return
 end subroutine angle_averaged_pk

 subroutine compute_damping_terms(nkbands,k,Pk_input,z,sigv,omega_k,omega_lambda,damp1,damp2,damp3)
   ! does integral over mu
   implicit none
   integer, intent(in) :: nkbands
   real, intent(in) :: z,sigv
   real, intent(in), dimension(nkbands) :: k
   real, intent(in), dimension(nkbands) :: Pk_input ! density Pk at redshift z
   real, intent(out), dimension(nkbands) :: damp1,damp2,damp3
   integer ik
   real omega_k, omega_lambda
   real(dl) mu_min, mu_max
   real(dl), parameter :: eps=1.d-4
   real(dl) rombint
   external rombint

   sigma_v = sigv
   mu_min = 0.
   mu_max = 1.

   do ik=1,nkbands
      k_eval = k(ik)
      damp1(ik) = real(rombint(func_damp1_int,mu_min,mu_max,eps))
      damp2(ik) = real(rombint(func_damp2_int,mu_min,mu_max,eps))
      damp3(ik) = real(rombint(func_damp3_int,mu_min,mu_max,eps))
    enddo

   return

 end subroutine compute_damping_terms

 


 function func_sigma_v(nk,kmin,dlnk,z,ok,ol,Pk_dd)
   ! function to compute sigma_v

   integer nk, ik
   real, dimension(nk) :: k, Ptt, Pk_dd, ln_k
   real kmin, dlnk
   real func_sigma_v, dk,z,ok,ol, a
   character(len=27) fname
   Omega_k_growth = ok
   Omega_lam_growth = ol
   Omega_mat_growth = 1.0-ok-ol
! we'll do this by simpsons rule
   func_sigma_v = 0.0

   do ik=1,nk
     k(ik) = kmin*exp(dlnk*real(ik-1))
   enddo
   a = 1.0/(1.0+z)
   Ptt(:) = Pk_dd(:)*growth_rate(a,ok,ol)**2

! Simpson's rule
!   do ik=2,nk-1
!      func_sigma_v = func_sigma_v + (1./6.)*(ln_k(ik)-ln_k(ik-1))*(Ptt(ik-1)+4.*Ptt(ik)+Ptt(ik+1))
!   enddo
! 5degree newton-cotes
   do ik=3,nk-1
      func_sigma_v = func_sigma_v + (1./24.)*(k(ik)-k(ik-1))*(11.0*Ptt(ik-2)+Ptt(ik-1)+Ptt(ik)+11.*Ptt(ik+1))
   enddo

   func_sigma_v = func_sigma_v*2.0/(3.0*(2.0*Pi_num)**2)
   return
 end function func_sigma_v


 function func_A0(mu)
   ! function to be integrated to give A0

   implicit none
   real(dl) func_A0, mu

   func_A0 = exp(-(growth_rate_eval*mu*k_eval*sigma_v)**2)

   return
 end function func_A0

 function func_A2(mu)
   ! function to be integrated to give A2

   implicit none
   real(dl) func_A2, mu

   func_A2 = mu**2*exp(-(growth_rate_eval*mu*k_eval*sigma_v)**2)

   return
 end function func_A2


 function func_A4(mu)
   ! function to be integrated to give A2

   implicit none
   real(dl) func_A4, mu

   func_A4 = mu**4*exp(-(growth_rate_eval*mu*k_eval*sigma_v)**2)

   return
 end function func_A4


 function func_damp1_int(mu)
   ! function to be integrated to give damp1
 
   implicit none
   real(dl) func_damp1_int, mu

   func_damp1_int = 1.d0/(1.d0+(mu*k_eval*(sigma_v/100.d0))**2)

   return
 end function func_damp1_int

function func_damp2_int(mu)
   ! function to be integrated to give damp1

   implicit none
   real(dl) func_damp2_int, mu

   func_damp2_int = mu**2/(1.d0+(mu*k_eval*(sigma_v/100.d0))**2)

   return
 end function func_damp2_int

function func_damp3_int(mu)
   ! function to be integrated to give damp1

   implicit none
   real(dl) func_damp3_int, mu

   func_damp3_int = mu**4/(1.d0+(mu*k_eval*(sigma_v/100.d0))**2)

   return
 end function func_damp3_int



 function growth_rate(a,omega_k,omega_lambda)
   ! f=dlnD/dlna

   real growth_rate, z,a,omega_k,omega_lambda
   real, parameter :: gamma = 0.545

   z = 1.0/a - 1.0
   Omega_mat_growth = 1.0 - omega_k - omega_lambda
   Omega_k_growth = omega_k
   Omega_lam_growth = omega_lambda
   growth_rate = omegamz(z)**gamma
   return
 end function growth_rate

 function func_growth_int(lna)
   ! function to be integrated to give growth factor

   real(dl) func_growth_int
   real(dl) z,a, oma, lna
   real(dl), parameter :: gamma = 0.545
   a = exp(lna)
   z = 1.d0/a - 1.d0
   oma = (Omega_mat_growth*(1.d0+z)**3)/(Omega_lam_growth+Omega_mat_growth*(1.d0+z)**3+Omega_k_growth*(1.d0+z)**2)
   func_growth_int = oma**gamma
    return
 end function func_growth_int

 function omegamz(z)
   ! Omega_matter as a function of redshift
   real omegamz,z
   
   omegamz = Omega_mat_growth*(1.+z)**3/(Omega_lam_growth+Omega_mat_growth*(1.+z)**3+Omega_k_growth*(1.+z)**2)
   return
 end function omegamz

 function growth(z,ok,ol)
   ! growth factor for a given redshift
   real z,growth,ok,ol, gv0
   real(dl) a, lna, lna_end
   real(dl), parameter :: eps=1.d-4
   real(dl) rombint
   external rombint

   a = 1.d0/(1.d0+z)
   lna = log(a)
   lna_end = log(1.d0)
   Omega_k_growth = ok
   Omega_lam_growth = ol
   Omega_mat_growth = 1.0-ok-ol
   growth = rombint(func_growth_int,lna_end,lna,eps)
   growth = exp(growth)
!   gv0 = rombint(func_growth_int,lna_end,1.d-2,eps)
!   gv0 = 1.0*exp(gv0)
!   growth = growth/gv0
   return
 end function growth


 function growth_rate2(a)
   real(dl) growth_rate2, E, a
   E = sqrt(Omega_lam_growth+Omega_mat_growth/a**3+Omega_k_growth/a**2)
   
   growth_rate2 = 1.0/((a*E)**3)
   return
 end function growth_rate2

 function growth2(z,ok,ol)
   
   real z, growth2, gv0
   real(dl) amax
   real E, ok, ol
   real, parameter :: eps=1.d-4
   real(dl) rombint
   external rombint

   Omega_k_growth = ok
   Omega_lam_growth = ol
   Omega_mat_growth = 1.0-ok-ol
   amax = 1.d0/(1.d0+z)
   E = sqrt(Omega_lam_growth+Omega_mat_growth*(1.+z)**3+Omega_k_growth*(1.+z)**2)
   growth2 = rombint(growth_rate2,1.d-10,amax,eps)
   gv0 = rombint(growth_rate2,1.d-10,1.d0,eps)
   growth2 = E*growth2/gv0
   return
 end function growth2
 
   function MatterPowerAt_zbin(T,kh,iz)
     !get matter power spectrum today at kh = k/h by interpolation from stored values
     real, intent(in) :: kh
     Type(CosmoTheory) T
     real MatterPowerAt_zbin
     real x, d
     integer i,iz
   
     x = log(kh/matter_power_minkh) / matter_power_dlnkh
     if (x < 0 .or. x >= num_matter_power-1) then
        write (*,*) ' k/h out of bounds in MatterPowerAt (',kh,')'
        stop 
     end if
     i = int(x)
     d = x - i
     MatterPowerAt_zbin = exp(log(T%matter_power(i+1,iz))*(1-d) &
       + log(T%matter_power(i+2,iz))*d)
     !Just do linear interpolation in logs for now..
     !(since we already cublic-spline interpolated to get the stored values)
     !Assume matter_power_lnzsteps is at redshift zero
   end function MatterPowerAt_zbin



end module
