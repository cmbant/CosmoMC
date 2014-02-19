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

    !JD 03/08/2013 fixed compute_scaling_factor and associated functions
    !to work with w_a/=0

    !JD 09/13: Replaced compute_scaling_factor routines with routines that use CAMB's
    !          built in D_V function.
    
    !JD 02/14  Moved common MPK functions to power_spec.f90 implemented AL's 
    !          Calculator_Cosmology functions.

    module wigglezinfo
    !David Parkinson 12th March 2012
    use settings
    use cmbtypes
    use CosmoTheory

    implicit none

    logical :: use_wigglez_mpk = .false.  !DP for WiggleZ MPK
    real(mcp), parameter :: z0 = 0.d0, za = 0.22d0, zb = 0.41d0, zc = 0.6d0, zd = 0.78d0
    real(mcp), dimension(4) :: zeval, zweight, sigma2BAOfid, sigma2BAO

    !settings power spectra evaluated at GiggleZ fiducial cosmological theory
    integer, parameter :: gigglez_num_matter_power = 500
    real(mcp), parameter :: gigglez_minkh =  0.999e-4
    real(mcp), parameter :: gigglez_dlnkh =  0.024
    !power spectra evaluated at GiggleZ fiducial cosmological theory
    real(mcp), dimension(gigglez_num_matter_power,4) :: power_hf_fid
    !for GiggleZ adjusted WiggleZ Power spectrum
    real(mcp), dimension(gigglez_num_matter_power) :: WiggleZPk
    contains

    subroutine GiggleZinfo_init(redshift)
    integer :: iopb, i, iz
    real(mcp) :: kval, power_nl
    real(mcp) redshift
    integer tmp_file_unit
    character(LEN=:), allocatable :: fname

    iz = 0
    do i=1,4
        if(abs(redshift-zeval(i)).le.0.001) iz = i
    enddo

    !! first read in everything needed from the CAMB output files.
    iopb = 0 !! check later if there was an error

    if(iz.eq.1) then
        fname = 'gigglezfiducialmodel_matterpower_a.dat'
    else if(iz.eq.2) then
        fname = 'gigglezfiducialmodel_matterpower_b.dat'
    else if(iz.eq.3) then
        fname = 'gigglezfiducialmodel_matterpower_c.dat'
    else if(iz.eq.4) then
        fname = 'gigglezfiducialmodel_matterpower_d.dat'
    else
        call MpiStop('could not indentify redshift')
    endif
    tmp_file_unit= OpenNewTxtFile(DataDir//fname)
    do i = 1, gigglez_num_matter_power
        read (tmp_file_unit,*,iostat=iopb) kval, power_nl
        if(iopb .ne. 0) stop 'Error reading model or fiducial theory files.'
        power_hf_fid(i,iz) = power_nl
    end do
    close(tmp_file_unit)

    end subroutine GiggleZinfo_init

    ! HARD CODING OF POLYNOMIAL FITS TO FOUR REDSHIFT BINS.
    subroutine GiggleZtoICsmooth(k,fidpolys)
    real(mcp), intent(in) :: k
    real(mcp) :: fidz_1, fidz_2, fidz_3, fidz_4
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

    subroutine fill_GiggleZTheory(Theory,zbin)
    class(TCosmoTheoryPredictions) Theory
    integer, intent(in) :: zbin
    real(mcp) :: xi, kval
    real(mcp), dimension(4) :: fidpolys
    real(mcp) pk_hf, holdval
    integer :: i,iz,ik

    iz = 0
    do i=1,4
        if(abs(Theory%redshifts(zbin)-zeval(i)).le.0.001) iz = i
    end do

    WiggleZPk = 0.

    do ik=1,gigglez_num_matter_power
        xi = log(gigglez_minkh) + gigglez_dlnkh*(ik-1)
        kval = exp(xi)
        pk_hf = Theory%MatterPowerAt_zbin(kval,zbin,.true.)

        call GiggleZtoICsmooth(kval,fidpolys)

        holdval = pk_hf*fidpolys(iz)/power_hf_fid(ik,iz)

        WiggleZPk(ik) = WiggleZPk(ik) + holdval
    end do

    end subroutine fill_GiggleZTheory

    function WiggleZPowerAt(kh)
    !get WiggleZ matter power spectrum today at kh = k/h by interpolation from stored values
    real(mcp), intent(in) :: kh
    real(mcp) WiggleZPowerAt
    real(mcp) x, d
    integer i

    x = log(kh/gigglez_minkh)/gigglez_dlnkh
    if (x < 0 .or. x >= gigglez_num_matter_power-1) then
        write (*,*) ' k/h out of bounds in WiggleZPowerAt (',kh,')'
        call MpiStop('')
    end if

    i = int(x)
    d = x - i
    WiggleZPowerAt = exp(log(WiggleZPk(i+1))*(1-d) + log(WiggleZPk(i+2))*d)
    !Just do linear interpolation in logs for now..
    !(since we already cublic-spline interpolated to get the stored values)
    end function WiggleZPowerAt

    end module wigglezinfo


    module wigglez
    use settings
    use cmbtypes
    use CosmoTheory
    use likelihood
    use wigglezinfo
    implicit none

    type, extends(DatasetFileLikelihood) :: TWiggleZCommon
    contains
    procedure :: ReadIni => TWiggleZCommon_ReadIni
    end type TWiggleZCommon

    type, extends(CosmologyLikelihood) :: WiggleZLikelihood
        logical :: use_set
        ! 1st index always refers to the region
        ! so mpk_P(1,:) is the Power spectrum in the first active region
        real(mcp), pointer, dimension(:,:,:) :: mpk_W, mpk_invcov
        real(mcp), pointer, dimension(:,:) :: mpk_P
        real(mcp), pointer, dimension(:) :: mpk_k
        !JD 09/13 added DV_fid so we can use camb routines to calculate a_scl
        real(mcp):: DV_fid
        real(mcp):: redshift ! important to know
    contains
    procedure :: LogLike => WiggleZ_Lnlike
    procedure :: ReadIni => WiggleZ_ReadIni
    end type WiggleZLikelihood

    type(TWiggleZCommon), save, target :: WiggleZCommon

    integer, parameter :: max_num_wigglez_regions = 7
    !Note all units are in k/h here
    integer, parameter :: mpk_d = kind(1.d0)

    !JD 09/13  Moved a bunch of stuff here so we only set it once and settings are
    !common across all used WiggleZ datasets
    integer :: num_mpk_points_use ! total number of points used (ie. max-min+1)
    integer :: num_mpk_kbands_use ! total number of kbands used (ie. max-min+1)
    integer :: num_regions_used   ! total number of wigglez regions being used

    integer :: num_mpk_points_full ! actual number of bandpowers in the infile
    integer :: num_mpk_kbands_full ! actual number of k positions " in the infile
    integer :: max_mpk_points_use ! in case you don't want the smallest scale modes (eg. sdss)
    integer :: min_mpk_points_use ! in case you don't want the largest scale modes
    integer :: max_mpk_kbands_use ! in case you don't want to calc P(k) on the smallest scales (will truncate P(k) to zero here!)
    integer :: min_mpk_kbands_use ! in case you don't want to calc P(k) on the largest scales (will truncate P(k) to zero here!)

    logical, pointer, dimension(:) :: regions_active

    logical :: use_scaling !as SDSS_lrgDR3 !JD 09/13 now using CAMB functions for a_scl

    logical :: use_gigglez = .false.
    logical :: nonlinear_wigglez = .false.

    !for Q and A see e.g. astro-ph/0501174, astro-ph/0604335
    logical :: Q_marge, Q_flat
    real(mcp):: Q_mid, Q_sigma, Ag



    contains

    subroutine WiggleZLikelihood_Add(LikeList, Ini)
    use IniObjects
    use settings
    use wigglezinfo
    class(LikelihoodList) :: LikeList
    class(TSettingIni) :: ini
    Type(WiggleZLikelihood), pointer :: like
    integer nummpksets, i

    use_wigglez_mpk = (Ini%Read_Logical('use_wigglez_mpk',.false.))

    if(.not. use_wigglez_mpk) return

    use_gigglez = Ini%Read_Logical('Use_gigglez',.false.)
    nonlinear_wigglez = Ini%Read_Logical('nonlinear_wigglez',.false.)

    call WiggleZCommon%ReadDatasetFile(Ini%ReadFileName('wigglez_common_dataset'))
    WiggleZCommon%LikelihoodType = 'MPK'

    nummpksets = Ini%Read_Int('mpk_wigglez_numdatasets',0)
    do i= 1, nummpksets
        allocate(like)
        like%LikelihoodType = 'MPK'
        like%needs_powerspectra = .true.
        like%needs_exact_z = .true.
        like%num_z = 1
        call like%ReadDatasetFile(Ini%ReadFileName(numcat('wigglez_dataset',i)))
        like%CommonData=> WiggleZCommon
        call LikeList%Add(like)
    end do
    if (Feedback>1) write(*,*) 'read WiggleZ MPK datasets'

    end subroutine WiggleZLikelihood_Add


    subroutine TWiggleZCommon_ReadIni(like,Ini)
    use IniObjects
    use wigglezinfo
    class(TWiggleZCommon) :: like
    class(TSettingIni) :: ini
    character(len=64) region_string
    integer i_regions

    zeval(1) = za
    zeval(2) = zb
    zeval(3) = zc
    zeval(4) = zd


    num_mpk_points_full = Ini%Read_Int('num_mpk_points_full',0)
    if (num_mpk_points_full.eq.0) write(*,*) ' ERROR: parameter num_mpk_points_full not set'
    num_mpk_kbands_full = Ini%Read_Int('num_mpk_kbands_full',0)
    if (num_mpk_kbands_full.eq.0) write(*,*) ' ERROR: parameter num_mpk_kbands_full not set'

    min_mpk_points_use = Ini%Read_Int('min_mpk_points_use',1)
    min_mpk_kbands_use = Ini%Read_Int('min_mpk_kbands_use',1)
    max_mpk_points_use = Ini%Read_Int('max_mpk_points_use',num_mpk_points_full)
    max_mpk_kbands_use = Ini%Read_Int('max_mpk_kbands_use',num_mpk_kbands_full)

    ! region 1 = 9h
    ! region 2 = 11h
    ! region 3 = 15h
    ! region 4 = 22h
    ! region 5 = 0h
    ! region 6 = 1h
    ! region 7 = 3h

    allocate(regions_active(max_num_wigglez_regions))
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
        regions_active(i_regions) =  Ini%Read_Logical(region_string,.false.)
    enddo

    !  ... work out how many regions are being used
    num_regions_used = 0
    do i_regions = 1,max_num_wigglez_regions
        if(regions_active(i_regions)) num_regions_used = num_regions_used + 1
    enddo

    if(num_regions_used.eq.0) then
        call MpiStop('WiggleZ_mpk: no regions being used in this data set')
    endif

    num_mpk_points_use = max_mpk_points_use - min_mpk_points_use +1
    num_mpk_kbands_use = max_mpk_kbands_use - min_mpk_kbands_use +1

    use_scaling = Ini%Read_Logical('use_scaling',.false.)

    if(use_gigglez .and. .not. nonlinear_wigglez) then
        write(*,*) 'ERROR!:  GiggleZ non-linear prescription only available'
        write(*,*) '         when setting nonlinear_wigglez = T in WiggleZ_MPK.ini'
        call MPIstop()
    end if

    if(.not. use_gigglez .and. nonlinear_wigglez)then
        write(*,*)'WARNING! Using non-linear model in WiggleZ module without'
        write(*,*)'GiggleZ prescription.  This method may not be as accurate.'
        write(*,*)'See arXiv:1210.2130 for details.'
    end if

    Q_marge = Ini%Read_Logical('Q_marge',.false.)
    if (Q_marge) then
        Q_flat = Ini%Read_Logical('Q_flat',.false.)
        if (.not. Q_flat) then
            !gaussian prior on Q
            Q_mid = Ini%Read_Real('Q_mid')
            Q_sigma = Ini%Read_Real('Q_sigma')
        end if
        Ag = Ini%Read_Real('Ag', 1.4)
    end if

    end subroutine TWiggleZCommon_ReadIni

    subroutine WiggleZ_ReadIni(like,Ini)
    ! this will be called once for each redshift bin
    use wigglezinfo
    use MatrixUtils
    implicit none
    class(WiggleZLikelihood) like
    class(TSettingIni) :: Ini
    character(LEN=:), allocatable :: kbands_file, measurements_file, windows_file, cov_file
    integer i,iopb,i_regions
    real(mcp) keff,klo,khi,beff
    real(mcp), dimension(:,:,:), allocatable :: mpk_Wfull, mpk_covfull
    real(mcp), dimension(:), allocatable :: mpk_kfull
    real(mcp), dimension(:,:), allocatable :: invcov_tmp
    character(80) :: dummychar
    integer count,tmp_file_unit

    iopb = 0

    allocate(like%exact_z(like%num_z))
    allocate(like%exact_z_index(like%num_z))
    like%exact_z(1) = Ini%Read_Double('redshift',0.d0)

    if(like%exact_z(1).eq.0.0) then
        call MpiStop('mpk: failed  to read in WiggleZ redshift')
    end if

    like%use_set =.true.
    if (Feedback > 0) write (*,*) 'reading: '//trim(like%name)

    if(allocated(mpk_kfull)) deallocate(mpk_kfull)
    allocate(mpk_kfull(num_mpk_kbands_full))
    allocate(like%mpk_P(num_regions_used,num_mpk_points_use))
    allocate(like%mpk_k(num_mpk_kbands_use))
    allocate(like%mpk_W(num_regions_used,num_mpk_points_use,num_mpk_kbands_use))

    kbands_file  = Ini%ReadFileName('kbands_file')
    call ReadVector(kbands_file,mpk_kfull,num_mpk_kbands_full)
    like%mpk_k(:)=mpk_kfull(min_mpk_kbands_use:max_mpk_kbands_use)
    if (Feedback > 1) then
        write(*,*) 'reading: '//trim(like%name)//' data'
        write(*,*) 'Using kbands windows between',real(like%mpk_k(1)),' < k/h < ',real(like%mpk_k(num_mpk_kbands_use))
    endif
    !if  (like%mpk_k(1) < matter_power_minkh) then
    !    write (*,*) 'WARNING: k_min in '//trim(like%name)//'less than setting in cmbtypes.f90'
    !    write (*,*) 'all k<matter_power_minkh will be set to matter_power_minkh'
    !end if

    measurements_file  = Ini%ReadFileName('measurements_file')
    tmp_file_unit = OpenNewTxtFile(measurements_file)
    like%mpk_P=0.
    count = 0
    do i_regions =1,7
        if(regions_active(i_regions)) then
            count = count+1
            read (tmp_file_unit,*) dummychar
            read (tmp_file_unit,*) dummychar
            do i= 1, (min_mpk_points_use-1)
                read (tmp_file_unit,*, iostat=iopb) keff,klo,khi,beff,beff,beff
            end do

            if (Feedback > 1 .and. min_mpk_points_use>1) write(*,*) 'Not using bands with keff=  ',real(keff),&
            ' or below in region', i_regions
            do i =1, num_mpk_points_use
                read (tmp_file_unit,*, iostat=iopb) keff,klo,khi,like%mpk_P(count,i),beff,beff
            end do
            ! NB do something to get to the end of the list
            do i=1, num_mpk_points_full-num_mpk_points_use-min_mpk_points_use+1
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
    windows_file  = Ini%ReadFileName('windows_file')
    if (windows_file.eq.'') write(*,*) 'ERROR: WiggleZ mpk windows_file not specified'
    call ReadWiggleZMatrices(windows_file,mpk_Wfull,max_num_wigglez_regions,num_mpk_points_full,num_mpk_kbands_full)
    count = 0
    do i_regions=1,max_num_wigglez_regions
        if(regions_active(i_regions)) then
            count = count + 1
            like%mpk_W(count,1:num_mpk_points_use,1:num_mpk_kbands_use)= &
            mpk_Wfull(i_regions,min_mpk_points_use:max_mpk_points_use,min_mpk_kbands_use:max_mpk_kbands_use)
        endif
    enddo

    deallocate(mpk_Wfull)
    !    deallocate(mpk_kfull)

    cov_file  = Ini%ReadFileName('cov_file')
    if (cov_file /= '') then
        allocate(mpk_covfull(max_num_wigglez_regions,num_mpk_points_full,num_mpk_points_full))
        allocate(invcov_tmp(num_mpk_points_use,num_mpk_points_use))
        ! ... read the entire covraiance matrix in, then decide which regions we want...
        call ReadWiggleZMatrices(cov_file,mpk_covfull,max_num_wigglez_regions,num_mpk_points_full,num_mpk_points_full)
        allocate(like%mpk_invcov(num_regions_used,num_mpk_points_use,num_mpk_points_use))
        count = 0
        do i_regions=1,max_num_wigglez_regions
            if(regions_active(i_regions)) then
                count = count + 1
                ! ... the covariance matrix has two indices for the different k-values, and another one for the region...
                !             like%mpk_invcov(count,1:num_mpk_points_use,1:num_mpk_points_use)=  &
                invcov_tmp(:,:) = &
                mpk_covfull(i_regions,min_mpk_points_use:max_mpk_points_use,min_mpk_points_use:max_mpk_points_use)
                !             call Matrix_Inverse(like%mpk_invcov(count,:,:))
                call Matrix_Inverse(invcov_tmp)
                like%mpk_invcov(count,1:num_mpk_points_use,1:num_mpk_points_use) = invcov_tmp(:,:)
            endif
        enddo
        deallocate(mpk_covfull)
        deallocate(invcov_tmp)
    else
        nullify(like%mpk_invcov)
    end if

    !JD 09/13 Read in fiducial D_V for use when calculating a_scl
    if(use_scaling) then
        like%DV_fid = Ini%Read_Double('DV_fid',-1.d0)
        if(like%DV_fid == -1.d0) then
            write(*,*)'ERROR: use_scaling = T and no DV_fid given '
            write(*,*)'       for dataset '//trim(like%name)//'.'
            write(*,*)'       Please check your .dataset files.'
            call MPIstop()
        end if
    end if

    if(use_gigglez) then
        call GiggleZinfo_init(like%exact_z(1))
    endif

    like%kmax = 0.8
    if(nonlinear_wigglez) then
        like%needs_nonlinear_pk = .true.
        like%kmax=1.2
    end if

    if (Feedback > 1) write(*,*) 'read: '//trim(like%name)//' data'

    if (iopb.ne.0) then
        stop 'Error reading WiggleZ mpk file'
    endif

    end subroutine WiggleZ_ReadIni

    subroutine ReadWiggleZMatrices(aname,mat,num_regions,m,n)
    ! suborutine to read all the matrices from each of the different regions, enclosed in one file

    implicit none
    character(LEN=*), intent(IN) :: aname
    integer, intent(in) :: m,n,num_regions
    real(mcp), intent(out) :: mat(num_regions,m,n)
    integer j,i_region,tmp_file_unit
    real(mcp) tmp
    character(LEN=64) dummychar

    if (Feedback > 1) write(*,*) 'reading: '//trim(aname)
    tmp_file_unit= OpenNewTxtFile(aname)
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
    Class(TCosmoTheoryPredictions) Theory
    real(mcp) :: DataParams(:)
    real(mcp) :: WiggleZ_LnLike, LnLike
    real(mcp), dimension(:), allocatable :: mpk_Pth, mpk_k2,mpk_lin,k_scaled !LV_06 added for LRGDR4
    real(mcp), dimension(:), allocatable :: mpk_WPth, mpk_WPth_k2
    real(mcp), dimension(:), allocatable :: diffs, step1
    real(mcp), dimension(:), allocatable :: Pk_delta_delta,Pk_delta_theta,Pk_theta_theta
    real(mcp), dimension(:), allocatable :: damp1, damp2, damp3
    real(mcp) :: covdat(num_mpk_points_use)
    real(mcp) :: covth(num_mpk_points_use)
    real(mcp) :: covth_k2(num_mpk_points_use)
    real(mcp), dimension(:), allocatable :: mpk_WPth_large, covdat_large, covth_large, mpk_Pdata_large
    integer imin,imax
    real :: normV, Q, minchisq
    real(mcp) :: a_scl  !LV_06 added for LRGDR4
    integer :: i, iQ,iz
    logical :: do_marge
    integer, parameter :: nQ=6
    integer, parameter :: nbias = 100
    real(mcp) old_chisq
    real(mcp) :: dQ = 0.4
    real(mcp), dimension(:), allocatable :: chisq(:)
    real(mcp) calweights(-nQ:nQ)
    real(mcp) vec2(2),Mat(2,2)
    real(mcp) final_term, b_out
    real(mcp) z
    integer i_region

    If(Feedback > 1) print*, 'Calling WiggleZ likelihood routines'
    allocate(mpk_lin(num_mpk_kbands_use),mpk_Pth(num_mpk_kbands_use))
    allocate(mpk_WPth(num_mpk_points_use))
    allocate(k_scaled(num_mpk_kbands_use))!LV_06 added for LRGDR4 !! IMPORTANT: need to check k-scaling
    !   allocate(num_mpk_points_use))
    allocate(diffs(num_mpk_points_use),step1(num_mpk_points_use))
    allocate(Pk_delta_delta(num_mpk_kbands_use),Pk_delta_theta(num_mpk_kbands_use))
    allocate(Pk_theta_theta(num_mpk_kbands_use))
    allocate(damp1(num_mpk_kbands_use),damp2(num_mpk_kbands_use),damp3(num_mpk_kbands_use))


    allocate(chisq(-nQ:nQ))

    chisq = 0

    if (.not. like%use_set) then
        LnLike = 0
        return
    end if

    z = 1.d0*dble(like%exact_z(1)) ! accuracy issues

    !JD 09/13 new compute_scaling_factor functions
    if(use_scaling) then
        call compute_scaling_factor(z,CMB,Theory,like%DV_fid,a_scl)
    else
        a_scl = 1
    end if

    iz = 0
    do i=1,4
        if(abs(z-zeval(i)).le.0.001) iz = i
    enddo
    if(iz.eq.0) call MpiStop('could not indentify redshift')

    if(abs(z-Theory%redshifts(like%exact_z_index(1)))>1.d-3)then
        write(*,*)'ERROR: WiggleZ redshift does not match the value stored'
        write(*,*)'       in the Theory%redshifts array.'
        call MpiStop()
    end if

    if(use_gigglez) then
        call fill_GiggleZTheory(Theory, like%exact_z_index(1))
    endif

    do i=1, num_mpk_kbands_use
        ! It could be that when we scale the k-values, the lowest bin drops off the bottom edge
        !Errors from using matter_power_minkh at lower end should be negligible
        k_scaled(i)=max(exp(Theory%log_kh(1)),like%mpk_k(i)*a_scl)
        if(use_gigglez) then
            mpk_lin(i) = WiggleZPowerAt(k_scaled(i))/a_scl**3
        else if(nonlinear_wigglez)then
            mpk_lin(i)=Theory%MatterPowerAt_zbin(k_scaled(i),like%exact_z_index(1),.true.)/a_scl**3
        else
            mpk_lin(i)=Theory%MatterPowerAt_zbin(k_scaled(i),like%exact_z_index(1))/a_scl**3
        endif
    end do

    do_marge = Q_marge
    if (do_marge .and. Q_flat) then
        !Marginalize analytically with flat prior on b^2 and b^2*Q
        !as recommended by Max Tegmark for SDSS
        allocate(mpk_k2(num_mpk_kbands_use))
        allocate(mpk_WPth_k2(num_mpk_points_use))

        Mat(:,:) = 0.d0
        vec2(:) = 0.d0
        final_term = 0.d0
        do i_region=1,num_regions_used
            mpk_Pth(:)=mpk_lin(:)/(1+Ag*k_scaled)
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
        if (Q_sigma==0) do_marge = .false.
        ! ... sum the chi-squared contributions for all regions first
        chisq(:) = 0.d0
        old_chisq = 1.d30
        if(feedback > 1) print*, "starting analytic marginalisation over bias"
        allocate(mpk_Pdata_large(num_mpk_points_use*num_regions_used))
        allocate(mpk_WPth_large(num_mpk_points_use*num_regions_used))
        allocate(covdat_large(num_mpk_points_use*num_regions_used))
        allocate(covth_large(num_mpk_points_use*num_regions_used))
        normV = 0.d0
        do iQ=-nQ,nQ
            Q = Q_mid +iQ*Q_sigma*dQ
            if (Q_marge) then
                mpk_Pth(:)=mpk_lin(:)*(1+Q*k_scaled(:)**2)/(1+Ag*k_scaled(:))
            else
                mpk_Pth(:) = mpk_lin(:)
            end if
            do i_region=1,num_regions_used
                imin = (i_region-1)*num_mpk_points_use+1
                imax = i_region*num_mpk_points_use
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
    if (Feedback>1) write(*,'("WiggleZ bin ",I0," MPK Likelihood = ",F10.5)')iz,LnLike

    if (LnLike > 1e8) then
        write(*,'("WARNING: WiggleZ bin",I0," Likelihood is huge!")')iz
        write(*,'("         Maybe there is a problem? Likelihood = ",F10.5)')LnLike
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
    ! JD 09/13: Replaced compute_scaling_factor routines so we use
    !           D_V calculations from CAMB.  New routines below
    
    subroutine compute_scaling_factor(z,CMB,Theory,DV_fid,a_scl)
    implicit none
    Class(TCosmoTheoryPredictions) Theory
    Class(CMBParams) CMB
    real(mcp), intent(in) :: z, DV_fid
    real(mcp), intent(out) :: a_scl

    !We use H_0*D_V because we dont care about scaling of h since
    !k is in units of h/Mpc
    !a_scl = CMB%H0*Theory%Config%Calculator%BAO_D_v(z)/DV_fid
    call MPIstop('Need to get BAO_D_V passed properly')
    !Like in original code, we need to apply a_scl in the correct direction
    a_scl = 1.0_mcp/a_scl
    end subroutine compute_scaling_factor
    
    end module wigglez
