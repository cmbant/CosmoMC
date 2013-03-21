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
real(mcp), allocatable :: ratio_power_nw_nl_fid(:,:)
!real(mcp),dimension(num_matter_power,matter_power_lnzsteps) :: ratio_power_nw_nl_fid
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

  open(unit=tmp_file_unit,file=trim(DataDir)//'lrgdr7fiducialmodel_matterpowerzNEAR.dat',form='formatted',err=500, iostat=ios)
  read (tmp_file_unit,*,iostat=iopb) getabstransferscalefiddummy, omegakdummy,omegavdummy,wdummy
  do i = 1, num_matter_power
    read (tmp_file_unit,*,iostat=iopb) kval, plin, psmooth, rationwhalofit
    ratio_power_nw_nl_fid(i,2) = rationwhalofit
  end do
  close(tmp_file_unit)

  open(unit=tmp_file_unit,file=trim(DataDir)//'lrgdr7fiducialmodel_matterpowerzMID.dat',form='formatted',err=500, iostat=ios)
  read (tmp_file_unit,*,iostat=iopb) getabstransferscalefiddummy,omegakdummy,omegavdummy,wdummy
  do i = 1, num_matter_power
    read (tmp_file_unit,*,iostat=iopb) kval, plin, psmooth, rationwhalofit
    ratio_power_nw_nl_fid(i,3) = rationwhalofit
  end do
  close(tmp_file_unit)

  open(unit=tmp_file_unit,file=trim(DataDir)//'lrgdr7fiducialmodel_matterpowerzFAR.dat',form='formatted',err=500,iostat=ios)
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
  Type(TheoryPredictions) Theory
  real(mcp), intent(in) :: minkh, dlnkh
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
    nlrat = (Theory%mpkrat_nw_nl(ik,matterpowerindx))/(ratio_power_nw_nl_fid(ik,matterpowerindx))
    call LRGtoICsmooth(kval,fidpolys)
    holdval(iz) = zweight(iz)*psmear*nlrat*fidpolys(iz)
    Theory%finalLRGtheoryPk(ik) = Theory%finalLRGtheoryPk(ik) + holdval(iz)
   end do

  end do

end subroutine fill_LRGTheory

end module

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
    real(mcp), pointer, dimension(:,:) :: N_inv
    real(mcp), pointer, dimension(:,:) :: mpk_W, mpk_invcov
    real(mcp), pointer, dimension(:) :: mpk_P, mpk_sdev, mpk_k
    real(mcp), pointer, dimension(:) :: mpk_zerowindowfxn
    real(mcp), pointer, dimension(:) :: mpk_zerowindowfxnsubtractdat
    real(mcp) :: mpk_zerowindowfxnsubdatnorm !!the 0th entry in windowfxnsubtract file
    logical :: use_scaling !as SDSS_lrgDR3
   !for Q and A see e.g. astro-ph/0501174, astro-ph/0604335
    logical :: Q_marge, Q_flat
    real(mcp) :: Q_mid, Q_sigma, Ag
   end Type mpkdataset

 integer :: num_mpk_datasets = 0
 Type(mpkdataset) mpkdatasets(10)

 !Note all units are in k/h here

  integer, parameter :: mpk_d = kind(1.d0)

  logical :: use_mpk = .false.

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

contains

  subroutine mpk_SetTransferRedshifts(redshifts)
   real(mcp), intent(inout) :: redshifts(*)
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

  subroutine ReadmpkDataset(gname)
    use MatrixUtils
    character(LEN=*), intent(IN) :: gname
    character(LEN=Ini_max_string_len) :: kbands_file, measurements_file, windows_file, cov_file
    !! added for the LRG window function subtraction
    character(LEN=Ini_max_string_len) :: zerowindowfxn_file, zerowindowfxnsubtractdat_file

    Type (mpkdataset) :: mset

    integer i,iopb
    real(mcp) keff,klo,khi,beff
    integer :: num_mpk_points_full ! actual number of bandpowers in the infile
    integer :: num_mpk_kbands_full ! actual number of k positions " in the infile
    integer :: max_mpk_points_use ! in case you don't want the smallest scale modes (eg. sdss)
    integer :: min_mpk_points_use ! in case you don't want the largest scale modes
    integer :: max_mpk_kbands_use ! in case you don't want to calc P(k) on the smallest scales (will truncate P(k) to zero here!)
    integer :: min_mpk_kbands_use ! in case you don't want to calc P(k) on the largest scales (will truncate P(k) to zero here!)
    real(mcp), dimension(:,:), allocatable :: mpk_Wfull, mpk_covfull
    real(mcp), dimension(:), allocatable :: mpk_kfull, mpk_fiducial

    real(mcp), dimension(:), allocatable :: mpk_zerowindowfxnfull
    real(mcp), dimension(:), allocatable :: mpk_zerowindowfxnsubfull

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


  function LSS_mpklike(Theory,mset,CMB) result(LnLike) ! LV_06 added CMB here
   Type (mpkdataset) :: mset
   Type (TheoryPredictions) Theory
   Type(CMBparams) CMB     !LV_06 added for LRGDR4
   real(mcp) LnLike
   real(mcp), dimension(:), allocatable :: mpk_Pth, mpk_k2,mpk_lin,k_scaled !LV_06 added for LRGDR4
   real(mcp), dimension(:), allocatable :: w
   real(mcp), dimension(:), allocatable :: mpk_WPth, mpk_WPth_k2
   real(mcp) :: covdat(mset%num_mpk_points_use), covth(mset%num_mpk_points_use),  covth_k2(mset%num_mpk_points_use)
   real(mcp) :: normV, Q, minchisq
   real(mcp) :: a_scl  !LV_06 added for LRGDR4
   integer :: i, iQ
   logical :: do_marge
   integer, parameter :: nQ=6
   real(mcp) :: tmp, dQ = 0.4
   real(mcp) chisq(-nQ:nQ)
   real(mcp) calweights(-nQ:nQ)
   real(mcp) vec2(2),Mat(2,2)

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

   IF(mset%use_scaling) then
      call compute_scaling_factor(dble(CMB%omk),dble(CMB%omv),dble(CMB%w),a_scl)
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
   Type (CMBParams) CMB
   Type (TheoryPredictions) Theory
   real(mcp) LSSLnLike
   integer i
   real(mcp) tot(num_mpk_datasets)

  do i=1, num_mpk_datasets
     if (mpkdatasets(i)%name == 'twodf') then
        stop 'twodf no longer supported - use data/2df_2005.dataset'
     else if (mpkdatasets(i)%name == 'lrg_2009') then
        tot(i) = LSS_LRG_mpklike(Theory,mpkdatasets(i),CMB)
     else
      tot(i) = LSS_mpklike(Theory,mpkdatasets(i),CMB) !LV_06 added CMB here
     end if
  end do
  LSSLnLike = SUM(tot)

 end function LSSLnLike

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

subroutine compute_scaling_factor(Ok,Ol,w,a)
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
  z = zeffDR7
  Hrelinv= 1/sqrt(Ol*(1+z)**(3*(1+w)) + Ok*(1+z)**2 + Om*(1+z)**3 + Or*(1+z)**4)
!  write(*,*) Ok,Ol,w
call compute_z_eta(Or,Ok,Ol,w,z,eta)
  tmp = sqrt(abs(Ok))
  if (Ok.lt.-1.d-6) eta = sin(tmp*eta)/tmp
  if (Ok.gt.1d-6)   eta = (exp(tmp*eta)-exp(-tmp*eta))/(2*tmp) ! sinh(tmp*eta)/tmp
  Ok0= 0
  Ol0= 0.75
  w0= -1
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
   real(mcp) :: fidLnLike
   type(TheoryPredictions) :: temptheory
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
   Type (TheoryPredictions) Theory
   Type(CMBparams) CMB                  !LV_06 added for LRGDR4
   real(mcp) LnLike
   integer :: i
   real(mcp), dimension(:), allocatable :: mpk_raw, mpk_Pth, mpk_Pth_k, mpk_Pth_k2, k_scaled
   real(mcp), dimension(:), allocatable :: mpk_WPth, mpk_WPth_k, mpk_WPth_k2
   real(mcp) :: covdat(mset%num_mpk_points_use), covth(mset%num_mpk_points_use), &
          & covth_k(mset%num_mpk_points_use), covth_k2(mset%num_mpk_points_use), &
          & covth_zerowin(mset%num_mpk_points_use)

   real(mcp), dimension(nptstot) :: chisq, chisqmarg  !! minus log likelihood list
   real(mcp) :: minchisq,maxchisq,deltaL

   real(dl) :: a1val, a2val, zerowinsub
   real(mcp) :: sumDD, sumDT, sumDT_k, sumDT_k2, sumTT,&
     &  sumTT_k, sumTT_k2, sumTT_k_k, sumTT_k_k2, sumTT_k2_k2, &
     &  sumDT_tot, sumTT_tot, &
     &  sumDT_zerowin, sumTT_zerowin, sumTT_k_zerowin, sumTT_k2_zerowin, sumTT_zerowin_zerowin

   real(mcp) :: sumzerow_Pth, sumzerow_Pth_k, sumzerow_Pth_k2

   real(mcp) :: a_scl      !LV_06 added for LRGDR4

   real(wp) :: temp1,temp2,temp3
   real(mcp) :: temp4

   !! added for no marg
   integer :: myminchisqindx
   real(mcp) :: currminchisq, currminchisqmarg, minchisqtheoryamp, chisqnonuis
   real(mcp) :: minchisqtheoryampnonuis, minchisqtheoryampminnuis
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

   IF(mset%use_scaling) then
      call compute_scaling_factor(dble(CMB%omk),dble(CMB%omv),dble(CMB%w),a_scl)
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

end module
