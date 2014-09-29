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
    use CosmologyTypes
    use CosmoTheory
    use Interpolation
    implicit none

    real(mcp), parameter :: za = 0.22d0, zb = 0.41d0, zc = 0.6d0, zd = 0.78d0
    real(mcp), dimension(4) :: zeval

    !settings power spectra evaluated at GiggleZ fiducial cosmological theory
    !power spectra evaluated at GiggleZ fiducial cosmological theory
    Type(TCubicSpline), allocatable ::  GiggleZPK(:)
    integer, parameter :: GiggleZ_numk = 500
    integer, parameter :: GiggleZ_numz = 4
    contains

    subroutine GiggleZinfo_init()
    integer :: iopb, ik, iz
    real(mcp),allocatable :: kval(:), power_nl(:)
    character(LEN=:), allocatable :: fname
    Type(TTextFile) :: F

    allocate(GiggleZPK(GiggleZ_numz))
    allocate(kval(GiggleZ_numk))
    allocate(power_nl(GiggleZ_numk))

    do iz=1,GiggleZ_numz
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
        end if
        call F%Open(DataDir//fname)
        do ik=1, GiggleZ_numk
            read (F%unit,*,iostat=iopb)kval(ik),power_nl(ik)
            if(iopb .ne. 0) stop 'Error reading model or fiducial theory files.'
            !JD PK arrays store log(k); we choose to store log(PK) for interpolation
        end do
        call F%Close()
        call GiggleZPK(iz)%Init(kval,log(power_nl),GiggleZ_numk)
    end do

    end subroutine GiggleZinfo_init

    ! HARD CODING OF POLYNOMIAL FITS TO FOUR REDSHIFT BINS.
    function GiggleZtoICsmooth(k,zbin)
    real(mcp), intent(in) :: k
    integer, intent(in) :: zbin
    real(mcp) :: GiggleZtoICsmooth
    real(mcp) :: fidz

    if(zbin==1)then
        fidz = (4.619d0 - 13.7787d0*k + 58.941d0*k**2 - 175.24d0*k**3 + 284.321d0*k**4 - 187.284d0*k**5)
    else if(zbin==2)then
        fidz = (4.63079d0 - 12.6293d0*k + 42.9265d0*k**2 - 91.8068d0*k**3 + 97.808d0*k**4 - 37.633d0*k**5)
    else if(zbin==3)then
        fidz = (4.69659d0 - 12.7287d0*k + 42.5681d0*k**2 - 89.5578d0*k**3 + 96.664d0*k**4 - 41.2564*k**5)
    else if(zbin==4)then
        fidz = (4.6849d0 - 13.4747d0*k + 53.7172d0*k**2 - 145.832d0*k**3 + 216.638d0*k**4 - 132.782*k**5)
    end if
    GiggleZtoICsmooth = 10._mcp**fidz

    end function GiggleZtoICsmooth

    !Calculate GiggleZ adjusted WiggleZ Power spectrum
    subroutine WiggleZPowerAt(kh,zbin,PK)
    real(mcp), intent(in) :: kh
    integer, intent(in) :: zbin
    real(mcp), intent(inout) :: PK

    PK = PK*GiggleZtoICsmooth(kh,zbin)/exp(GiggleZPK(zbin)%value(kh))

    end subroutine WiggleZPowerAt

    end module wigglezinfo


    module wigglez
    use settings
    use CosmologyTypes
    use CosmoTheory
    use likelihood
    use wigglezinfo
    use MPK_Common
    use MatrixUtils
    implicit none
    private

    type, extends(TDatasetFileLikelihood) :: TWiggleZCommon
    contains
    procedure :: ReadIni => TWiggleZCommon_ReadIni
    end type TWiggleZCommon

    type, extends(TCosmologyPKLikelihood) :: WiggleZLikelihood
        ! 1st index always refers to the region
        ! so mpk_P(1,:) is the Power spectrum in the first active region
        Type(TPKLikelihoodCommon), dimension(:), allocatable :: PKData
        !Redshift bin that is being probed by a particular call
        integer :: zbin
    contains
    procedure :: LogLike => WiggleZ_Lnlike
    procedure :: ReadIni => WiggleZ_ReadIni
    end type WiggleZLikelihood

    type(TWiggleZCommon), save, target :: WiggleZCommon

    integer, parameter :: max_num_wigglez_regions = 7

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

    logical, allocatable, dimension(:) :: regions_active

    logical :: use_scaling !as SDSS_lrgDR3 !JD 09/13 now using CAMB functions for a_scl

    logical :: use_gigglez = .false.
    logical :: nonlinear_wigglez = .false.
    logical :: use_wigglez_mpk = .false.  !DP for WiggleZ MPK

    !for Q and A see e.g. astro-ph/0501174, astro-ph/0604335
    logical :: Q_marge, Q_flat
    real(mcp):: Q_mid, Q_sigma, Ag

    public TWiggleZCommon, WiggleZLikelihood_Add
    contains

    subroutine WiggleZLikelihood_Add(LikeList, Ini)
    class(TLikelihoodList) :: LikeList
    class(TSettingIni) :: ini
    Type(WiggleZLikelihood), pointer :: this
    integer nummpksets, i

    use_wigglez_mpk = (Ini%Read_Logical('use_wigglez_mpk',.false.))

    if(.not. use_wigglez_mpk) return

    use_gigglez = Ini%Read_Logical('Use_gigglez',.false.)
    nonlinear_wigglez = Ini%Read_Logical('nonlinear_wigglez',.false.)

    call WiggleZCommon%ReadDatasetFile(Ini%ReadFileName('wigglez_common_dataset'))
    WiggleZCommon%LikelihoodType = 'MPK'

    nummpksets = Ini%Read_Int('mpk_wigglez_numdatasets',0)
    do i= 1, nummpksets
        allocate(this)
        this%LikelihoodType = 'MPK'
        this%needs_powerspectra = .true.
        this%needs_exact_z = .true.
        this%num_z = 1
        this%needs_nonlinear_pk = nonlinear_wigglez
        call this%ReadDatasetFile(Ini%ReadFileName(numcat('wigglez_dataset',i)))
        this%CommonData=> WiggleZCommon
        call LikeList%Add(this)
    end do
    if (Feedback>1) write(*,*) 'read WiggleZ MPK data sets'

    end subroutine WiggleZLikelihood_Add


    subroutine TWiggleZCommon_ReadIni(this,Ini)
    class(TWiggleZCommon) :: this
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

    if(use_gigglez) then
        call GiggleZinfo_init()
    endif

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

    subroutine WiggleZ_ReadIni(this,Ini)
    ! this will be called once for each redshift bin
    class(WiggleZLikelihood) this
    class(TSettingIni) :: Ini
    character(LEN=:), allocatable :: kbands_file, measurements_file, windows_file, cov_file
    integer i,iopb,i_regions
    real(mcp) keff,klo,khi,beff
    real(mcp), dimension(:,:,:), allocatable :: mpk_Wfull, mpk_covfull
    real(mcp), dimension(:), allocatable :: mpk_kfull
    real(mcp), dimension(:,:), allocatable :: invcov_tmp
    character(80) :: dummychar
    integer count
    Type(TTextFile) :: F

    iopb = 0

    allocate(this%exact_z(this%num_z))
    this%exact_z(1) = Ini%Read_Double('redshift',0.d0)

    if(this%exact_z(1).eq.0.0) then
        call MpiStop('WiggleZMPK: failed to read in WiggleZ redshift')
    end if

    this%zbin = 0
    do i=1,4
        if(abs(this%exact_z(1)-zeval(i)).le.0.001) this%zbin = i
    enddo

    if(this%zbin.eq.0) call MpiStop('WiggleZMPK: could not indentify redshift')

    if (Feedback > 0) write (*,*) 'reading: '//trim(this%name)

    allocate(this%PKdata(num_regions_used))

    if(allocated(mpk_kfull)) deallocate(mpk_kfull)
    allocate(mpk_kfull(num_mpk_kbands_full))

    kbands_file  = Ini%ReadFileName('kbands_file')
    call File%ReadTextVector(kbands_file,mpk_kfull,num_mpk_kbands_full)
    if (Feedback > 1) then
        write(*,*) 'reading: '//trim(this%name)//' data'
        write(*,*) 'Using kbands windows between',real(mpk_kfull(min_mpk_kbands_use)),&
        ' < k/h < ',real(mpk_kfull(max_mpk_kbands_use))
    endif

    measurements_file  = Ini%ReadFileName('measurements_file')
    call F%Open(measurements_file)
    count = 0
    do i_regions =1,7
        if(regions_active(i_regions)) then
            count = count+1
            allocate(this%PKdata(count)%mpk_P(num_mpk_points_use))
            allocate(this%PKdata(count)%mpk_k(num_mpk_kbands_use))
            allocate(this%PKdata(count)%mpk_W(num_mpk_points_use,num_mpk_kbands_use))
            this%PKdata(count)%mpk_k(:)=mpk_kfull(min_mpk_kbands_use:max_mpk_kbands_use)
            this%PKdata(count)%mpk_P=0.

            read (F%unit,*) dummychar
            read (F%unit,*) dummychar
            do i= 1, (min_mpk_points_use-1)
                read (F%unit,*, iostat=iopb) keff,klo,khi,beff,beff,beff
            end do

            if (Feedback > 1 .and. min_mpk_points_use>1) write(*,*) 'Not using bands with keff=  ',real(keff),&
            ' or below in region', i_regions
            do i =1, num_mpk_points_use
                read (F%unit,*, iostat=iopb) keff,klo,khi,this%PKdata(count)%mpk_P(i),beff,beff
            end do
            ! NB do something to get to the end of the list
            do i=1, num_mpk_points_full-num_mpk_points_use-min_mpk_points_use+1
                read (F%unit,*, iostat=iopb) klo,klo,khi,beff,beff,beff
                if(iopb.ne.0) stop
            end do
        else
            read (F%unit,*) dummychar
            read (F%unit,*) dummychar
            do i=1,50
                read (F%unit,*, iostat=iopb) klo,klo,khi,beff,beff,beff
                if(iopb.ne.0) stop
            enddo
        endif
    enddo
    call F%Close()
    if (Feedback > 1) write(*,*) 'bands truncated at keff=  ',real(keff)

    allocate(mpk_Wfull(max_num_wigglez_regions,num_mpk_points_full,num_mpk_kbands_full))
    windows_file  = Ini%ReadFileName('windows_file')
    if (windows_file.eq.'') write(*,*) 'ERROR: WiggleZ mpk windows_file not specified'
    call ReadWiggleZMatrices(windows_file,mpk_Wfull,max_num_wigglez_regions,num_mpk_points_full,num_mpk_kbands_full)
    count = 0
    do i_regions=1,max_num_wigglez_regions
        if(regions_active(i_regions)) then
            count = count + 1
            this%PKdata(count)%mpk_W(1:num_mpk_points_use,1:num_mpk_kbands_use)= &
            mpk_Wfull(i_regions,min_mpk_points_use:max_mpk_points_use,min_mpk_kbands_use:max_mpk_kbands_use)
        endif
    enddo

    cov_file  = Ini%ReadFileName('cov_file')
    if (cov_file /= '') then
        allocate(mpk_covfull(max_num_wigglez_regions,num_mpk_points_full,num_mpk_points_full))
        allocate(invcov_tmp(num_mpk_points_use,num_mpk_points_use))
        ! ... read the entire covraiance matrix in, then decide which regions we want...
        call ReadWiggleZMatrices(cov_file,mpk_covfull,max_num_wigglez_regions,num_mpk_points_full,num_mpk_points_full)
        count = 0
        do i_regions=1,max_num_wigglez_regions
            if(regions_active(i_regions)) then
                count = count + 1
                allocate(this%PKdata(count)%mpk_invcov(num_mpk_points_use,num_mpk_points_use))
                invcov_tmp(:,:) = &
                mpk_covfull(i_regions,min_mpk_points_use:max_mpk_points_use,min_mpk_points_use:max_mpk_points_use)
                call Matrix_Inverse(invcov_tmp)
                this%PKdata(count)%mpk_invcov(1:num_mpk_points_use,1:num_mpk_points_use) = invcov_tmp(:,:)
            endif
        enddo
    end if

    !JD 09/13 Read in fiducial D_V for use when calculating a_scl
    if(use_scaling) then
        this%DV_fid = Ini%Read_Double('DV_fid',-1.d0)
        if(this%DV_fid == -1.d0) then
            write(*,*)'ERROR: use_scaling = T and no DV_fid given '
            write(*,*)'       for dataset '//trim(this%name)//'.'
            write(*,*)'       Please check your .dataset files.'
            call MPIstop()
        end if
    end if

    if(this%needs_nonlinear_pk) then
        this%kmax=1.2_mcp
    else
        this%kmax=0.8_mcp
    end if

    if (Feedback > 1) write(*,*) 'read: '//trim(this%name)//' data'

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
    integer j,i_region
    real(mcp) tmp
    character(LEN=64) dummychar
    type(TTextFile) :: F

    if (Feedback > 1) write(*,*) 'reading: '//trim(aname)
    call F%Open(aname)
    do i_region=1,num_regions
        read (F%unit,*, end = 200, err=100) dummychar
        do j=1,m
            read (F%unit,*, end = 200, err=100) mat(i_region,j,1:n)
        enddo
    enddo
    goto 120

100 write(*,*) 'matrix file '//trim(aname)//' is the wrong size',i_region,j,n,mat(num_regions,m,n)
    stop

120 read (F%unit,*, err = 200, end =200) tmp
    goto 200


200 call F%Close()
    return

    end subroutine ReadWiggleZMatrices

    function WiggleZ_LnLike(this,CMB,Theory,DataParams) ! LV_06 added CMB here
    Class(CMBParams) CMB
    Class(WiggleZLikelihood) :: this
    Class(TCosmoTheoryPredictions), target :: Theory
    Type (TCosmoTheoryPK), pointer :: PK
    real(mcp) :: DataParams(:)
    real(mcp) :: WiggleZ_LnLike, LnLike
    real(mcp), dimension(:), allocatable :: mpk_Pth, mpk_k2,mpk_lin,k_scaled !LV_06 added for LRGDR4
    real(mcp), dimension(:), allocatable :: mpk_WPth, mpk_WPth_k2
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

    allocate(chisq(-nQ:nQ))

    if (this%needs_nonlinear_pk) then
        if(.not. allocated(Theory%NL_MPK))then
            write(*,*) 'ERROR: Your Theory%NL_MPK derived type is not initialized. Make sure you are'
            write(*,*) '       calling a SetPk routine and filling your power spectra.'
            call MPIstop()
        end if
        PK=>Theory%NL_MPK
    else
        if(.not. allocated(Theory%MPK))then
            write(*,*) 'ERROR: Your Theory%MPK derived type is not initialized. Make sure you are'
            write(*,*) '       calling a SetPk routine and filling your power spectra.'
            call MPIstop()
        end if
        PK=>Theory%MPK
    end if

    chisq = 0

    z = this%exact_z(1)

    if(abs(z-PK%y(this%exact_z_index(1)))>1.d-3)then
        write(*,*)'ERROR: WiggleZ redshift does not match the value stored'
        write(*,*)'       in the PK%y array.'
        call MpiStop()
    end if

    !JD 09/13 new compute_scaling_factor functions
    if(use_scaling) then
        a_scl = this%compute_scaling_factor(z,CMB)
    else
        a_scl = 1
    end if

    do i=1, num_mpk_kbands_use
        ! It could be that when we scale the k-values, the lowest bin drops off the bottom edge
        !Errors from using matter_power_minkh at lower end should be negligible
        k_scaled(i)=max(exp(PK%x(1)),this%PKdata(1)%mpk_k(i)*a_scl)
        mpk_lin(i) = PK%PowerAt(k_scaled(i),this%exact_z(1))/a_scl**3
        if(use_gigglez) call WiggleZPowerAt(k_scaled(i),this%zbin,mpk_lin(i))
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


            mpk_WPth(:) = matmul(this%PKdata(i_region)%mpk_w(:,:),mpk_Pth(:))
            mpk_WPth_k2(:) = matmul(this%PKdata(i_region)%mpk_w(:,:),mpk_k2(:))

            covdat(:) = matmul(this%PKdata(i_region)%mpk_invcov(:,:),this%PKdata(i_region)%mpk_p(:))
            covth(:) = matmul(this%PKdata(i_region)%mpk_invcov(:,:),mpk_WPth(:))
            covth_k2(:) = matmul(this%PKdata(i_region)%mpk_invcov(:,:),mpk_WPth_k2(:))

            Mat(1,1) = Mat(1,1) + sum(covth(:)*mpk_WPth(:))
            Mat(2,2) = Mat(2,2) + sum(covth_k2(:)*mpk_WPth_k2(:))
            Mat(1,2) = Mat(1,2) + sum(covth(:)*mpk_WPth_k2(:))
            Mat(2,1) = Mat(1,2)

            vec2(1) = vec2(1) + sum(covdat(:)*mpk_WPth(:))
            vec2(2) = vec2(2) + sum(covdat(:)*mpk_WPth_k2(:))
            final_term = final_term + sum(this%PKdata(i_region)%mpk_p(:)*covdat(:))
        enddo
        LnLike = log( Mat(1,1)*Mat(2,2)-Mat(1,2)**2)
        call Matrix_Inverse(Mat)
        !          LnLike = (sum(mset%mpk_P*covdat) - sum(vec2*matmul(Mat,vec2)) + LnLike ) /2
        LnLike = (final_term - sum(vec2*matmul(Mat,vec2)) + LnLike ) /2
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
                mpk_WPth(:) = matmul(this%PKdata(i_region)%mpk_w(:,:),mpk_Pth(:))
                mpk_Pdata_large(imin:imax) = this%PKdata(i_region)%mpk_p(:)
                mpk_WPth_large(imin:imax) = mpk_WPth(:)

                !with analytic marginalization over normalization nuisance (flat prior on b^2)
                !See appendix F of cosmomc paper

                covdat_large(imin:imax) = matmul(this%PKdata(i_region)%mpk_invcov(:,:),this%PKdata(i_region)%mpk_p(:))
                covth_large(imin:imax) = matmul(this%PKdata(i_region)%mpk_invcov(:,:),mpk_WPth(:))
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
    if (Feedback>1) write(*,'("WiggleZ bin ",I0," MPK Likelihood = ",F10.5)')this%zbin,LnLike

    if (LnLike > 1e8) then
        write(*,'("WARNING: WiggleZ bin",I0," Likelihood is huge!")')this%zbin
        write(*,'("         Maybe there is a problem? Likelihood = ",F10.5)')LnLike
    end if

    end function WiggleZ_LnLike

    end module wigglez
