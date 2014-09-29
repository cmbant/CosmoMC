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

    !JD 29/07/2013 removed LRG stuff
    !JD 03/08/2013 fixed compute_scaling_factor and associated functions
    !to work with w_a/=0

    !JD 09/13: Replaced compute_scaling_factor routines with routines that use CAMB's
    !          built in D_V function.

    !JD 02/14  CosmoTheory changes;  Added MPK_Common

    module MPK_Common
    use settings
    use CosmologyTypes
    use CosmoTheory
    use Calculator_Cosmology
    use Likelihood_Cosmology
    implicit none
    private

    Type TPKLikelihoodCommon
        real(mcp), allocatable, dimension(:,:) :: mpk_W, mpk_invcov
        real(mcp), allocatable, dimension(:) :: mpk_P, mpk_sdev, mpk_k
    end type TPKLikelihoodCommon

    type, extends(TCosmoCalcLikelihood) :: TCosmologyPKLikelihood
        real(mcp) DV_fid   !Fiducial D_V
    contains
    procedure :: compute_scaling_factor
    end type TCosmologyPKLikelihood

    public TCosmologyPKLikelihood, TPKLikelihoodCommon
    contains

    function compute_scaling_factor(this,z,CMB)
    Class(TCosmologyPKLikelihood) this
    Class(CMBParams) CMB
    real(mcp), intent(in) :: z
    real(mcp) :: compute_scaling_factor

    !We use H_0*D_V because we dont care about scaling of h since
    !k is in units of h/Mpc
    compute_scaling_factor = this%DV_fid/(CMB%H0*this%Calculator%BAO_D_v(z))

    end function compute_scaling_factor

    end module MPK_Common

    module mpk
    use settings
    use CosmologyTypes
    use CosmoTheory
    use likelihood
    use MatrixUtils
    use MPK_Common
    implicit none
    private

    type, extends(TCosmologyPKLikelihood) :: MPKLikelihood
        logical :: use_set
        integer :: num_mpk_points_use ! total number of points used (ie. max-min+1)
        integer :: num_mpk_kbands_use ! total number of kbands used (ie. max-min+1)
        logical :: use_scaling !as SDSS_lrgDR3
        !for Q and A see e.g. astro-ph/0501174, astro-ph/0604335
        logical :: Q_marge, Q_flat
        real(mcp) :: Q_mid, Q_sigma, Ag
        Type(TPKLikelihoodCommon) :: PKData
    contains
    procedure :: LogLike => MPK_Lnlike
    procedure :: ReadIni => MPK_ReadIni
    end type MPKLikelihood

    logical :: use_mpk = .false.

    public use_mpk, MPKLikelihood, MPKLikelihood_Add

    contains

    subroutine MPKLikelihood_Add(LikeList, Ini)
    use IniObjects
    use settings
    class(TLikelihoodList) :: LikeList
    class(TSettingIni) :: ini
    type(MPKLikelihood), pointer :: this
    integer nummpksets, i

    use_mpk = (Ini%Read_Logical('use_mpk',.false.))

    if(.not. use_mpk) return

    nummpksets = Ini%Read_Int('mpk_numdatasets',0)
    do i= 1, nummpksets
        allocate(this)
        this%LikelihoodType = 'MPK'
        this%needs_powerspectra = .true.
        this%needs_exact_z = .true.
        this%num_z = 1
        this%needs_nonlinear_pk = Ini%Read_Logical(numcat('mpk_dataset_nonlinear',i),.false.)
        call this%ReadDatasetFile(Ini%ReadFileName(numcat('mpk_dataset',i)) )
        call LikeList%Add(this)
    end do
    if (Feedback>1 .and. nummpksets>0) write(*,*) 'read MPK data sets'

    end subroutine MPKLikelihood_Add

    subroutine MPK_ReadIni(this,Ini)
    use MatrixUtils
    class(MPKLikelihood) this
    class(TSettingIni) :: Ini
    character(LEN=:), allocatable :: kbands_file, measurements_file, windows_file, cov_file

    integer i,iopb,i_conflict
    real(mcp) keff,klo,khi,beff
    integer :: num_mpk_points_full ! actual number of bandpowers in the infile
    integer :: num_mpk_kbands_full ! actual number of k positions " in the infile
    integer :: max_mpk_points_use ! in case you don't want the smallest scale modes (eg. sdss)
    integer :: min_mpk_points_use ! in case you don't want the largest scale modes
    integer :: max_mpk_kbands_use ! in case you don't want to calc P(k) on the smallest scales (will truncate P(k) to zero here!)
    integer :: min_mpk_kbands_use ! in case you don't want to calc P(k) on the largest scales (will truncate P(k) to zero here!)
    real(mcp), dimension(:,:), allocatable :: mpk_Wfull, mpk_covfull
    real(mcp), dimension(:), allocatable :: mpk_kfull, mpk_fiducial
    Type(TTextFile) :: F

    character(80) :: dummychar

    if (this%name == 'twodf') stop 'twodf no longer supported - use data/2df_2005.dataset'
    if (this%name == 'lrg_2009') stop 'lrg_2009 no longer supported'
    this%use_set =.true.
    if (Feedback > 0) write (*,*) 'reading: '//trim(this%name)
    num_mpk_points_full = Ini%Read_Int('num_mpk_points_full',0)
    if (num_mpk_points_full.eq.0) write(*,*) ' ERROR: parameter num_mpk_points_full not set'
    num_mpk_kbands_full = Ini%Read_Int('num_mpk_kbands_full',0)
    if (num_mpk_kbands_full.eq.0) write(*,*) ' ERROR: parameter num_mpk_kbands_full not set'
    min_mpk_points_use = Ini%Read_Int('min_mpk_points_use',1)
    min_mpk_kbands_use = Ini%Read_Int('min_mpk_kbands_use',1)
    max_mpk_points_use = Ini%Read_Int('max_mpk_points_use',num_mpk_points_full)
    max_mpk_kbands_use = Ini%Read_Int('max_mpk_kbands_use',num_mpk_kbands_full)
    this%num_mpk_points_use = max_mpk_points_use - min_mpk_points_use +1
    this%num_mpk_kbands_use = max_mpk_kbands_use - min_mpk_kbands_use +1

    allocate(mpk_Wfull(num_mpk_points_full,num_mpk_kbands_full))
    allocate(mpk_kfull(num_mpk_kbands_full))
    allocate(this%PKData%mpk_P(this%num_mpk_points_use))
    allocate(this%PKData%mpk_sdev(this%num_mpk_points_use))  ! will need to replace with the covmat
    allocate(this%PKData%mpk_k(this%num_mpk_kbands_use))
    allocate(this%PKData%mpk_W(this%num_mpk_points_use,this%num_mpk_kbands_use))
    allocate(mpk_fiducial(this%num_mpk_points_use))

    kbands_file  = Ini%ReadFileName('kbands_file')
    call File%ReadTextVector(kbands_file,mpk_kfull,num_mpk_kbands_full)
    this%PKData%mpk_k(1:this%num_mpk_kbands_use)=mpk_kfull(min_mpk_kbands_use:max_mpk_kbands_use)
    if (Feedback > 1) then
        write(*,*) 'reading: '//trim(this%name)//' data'
        write(*,*) 'Using kbands windows between',real(this%PKData%mpk_k(1)),' < k/h < ', &
        & real(this%PKData%mpk_k(this%num_mpk_kbands_use))
    endif
    !if  (this%PKData%mpk_k(1) < matter_power_minkh) then
    !    write (*,*) 'WARNING: k_min in '//trim(this%name)//'less than setting in CosmologyTypes.f90'
    !    write (*,*) 'all k<matter_power_minkh will be set to matter_power_minkh'
    !end if

    measurements_file  = Ini%ReadFileName('measurements_file')
    call F%Open(measurements_file)
    this%PKData%mpk_P=0.
    read (F%unit,*) dummychar
    read (F%unit,*) dummychar
    do i= 1, (min_mpk_points_use-1)
        read (F%unit,*, iostat=iopb) keff,klo,khi,beff,beff,beff
    end do
    if (Feedback > 1 .and. min_mpk_points_use>1) write(*,*) 'Not using bands with keff=  ',real(keff),' or below'
    do i =1, this%num_mpk_points_use
        read (F%unit,*, iostat=iopb) keff,klo,khi,this%PKData%mpk_P(i),this%PKData%mpk_sdev(i),mpk_fiducial(i)
    end do
    call F%Close()
    if (Feedback > 1) write(*,*) 'bands truncated at keff=  ',real(keff)

    windows_file  = Ini%ReadFileName('windows_file')
    if (windows_file.eq.'') write(*,*) 'ERROR: mpk windows_file not specified'
    call File%ReadTextMatrix(windows_file,mpk_Wfull,num_mpk_points_full,num_mpk_kbands_full)
    this%PKData%mpk_W(1:this%num_mpk_points_use,1:this%num_mpk_kbands_use)= &
    mpk_Wfull(min_mpk_points_use:max_mpk_points_use,min_mpk_kbands_use:max_mpk_kbands_use)


    cov_file  = Ini%ReadFileName('cov_file')
    if (cov_file /= '') then
        allocate(mpk_covfull(num_mpk_points_full,num_mpk_points_full))
        call File%ReadTextMatrix(cov_file,mpk_covfull,num_mpk_points_full,num_mpk_points_full)
        allocate(this%PKData%mpk_invcov(this%num_mpk_points_use,this%num_mpk_points_use))
        this%PKData%mpk_invcov=  mpk_covfull(min_mpk_points_use:max_mpk_points_use,min_mpk_points_use:max_mpk_points_use)
        call Matrix_Inverse(this%PKData%mpk_invcov)
    end if

    this%use_scaling = Ini%Read_Logical('use_scaling',.false.)


    !JD 09/13 Read in fiducial D_V and redshift for use when calculating a_scl
    if(this%use_scaling) then
        !DV_fid should be in units CMB%H0*BAO_D_v(z)
        this%DV_fid = Ini%Read_Double('DV_fid',-1.d0)
        if(this%DV_fid == -1.d0) then
            write(*,*)'ERROR: use_scaling = T and no DV_fid given '
            write(*,*)'       for dataset '//trim(this%name)//'.'
            write(*,*)'       Please check your .dataset files.'
            call MPIstop()
        end if
    end if

    !JD 10/13  New settings for new mpk array handling
    allocate(this%exact_z(this%num_z))
    this%exact_z(1) = Ini%Read_Double('redshift',0.35d0)

    if(this%needs_nonlinear_pk) then
        this%kmax=1.2_mcp
    else
        this%kmax=0.8_mcp
    end if

    this%Q_marge = Ini%Read_Logical('Q_marge',.false.)
    if (this%Q_marge) then
        this%Q_flat = Ini%Read_Logical('Q_flat',.false.)
        if (.not. this%Q_flat) then
            !gaussian prior on Q
            this%Q_mid = Ini%Read_Real('Q_mid')
            this%Q_sigma = Ini%Read_Real('Q_sigma')
        end if
        this%Ag = Ini%Read_Real('Ag', 1.4)
    end if

    if (iopb.ne.0) then
        stop 'Error reading mpk file'
    endif

    end subroutine MPK_ReadIni


    function MPK_LnLike(this,CMB,Theory,DataParams) ! LV_06 added CMB here
    Class(CMBParams) CMB
    Class(MPKLikelihood) :: this
    Class(TCosmoTheoryPredictions), target :: Theory
    Type(TCosmoTheoryPK), pointer :: PK
    real(mcp) :: DataParams(:)
    real(mcp) MPK_LnLike,LnLike
    real(mcp), dimension(:), allocatable :: mpk_Pth, mpk_k2,mpk_lin,k_scaled !LV_06 added for LRGDR4
    real(mcp), dimension(:), allocatable :: w
    real(mcp), dimension(:), allocatable :: mpk_WPth, mpk_WPth_k2
    real(mcp) :: covdat(this%num_mpk_points_use), covth(this%num_mpk_points_use),  covth_k2(this%num_mpk_points_use)
    real(mcp) :: normV, Q, minchisq
    real(mcp) :: a_scl  !LV_06 added for LRGDR4
    integer :: i, iQ
    logical :: do_marge
    integer, parameter :: nQ=6
    real(mcp) :: tmp, dQ = 0.4
    real(mcp) chisq(-nQ:nQ)
    real(mcp) calweights(-nQ:nQ)
    real(mcp) vec2(2),Mat(2,2)

    allocate(mpk_lin(this%num_mpk_kbands_use) ,mpk_Pth(this%num_mpk_kbands_use))
    allocate(mpk_WPth(this%num_mpk_points_use))
    allocate(k_scaled(this%num_mpk_kbands_use))!LV_06 added for LRGDR4
    allocate(w(this%num_mpk_points_use))

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

    if (.not. this%use_set) then
        LnLike = 0
        return
    end if

    !JD 09/13 new compute_scaling_factor functions
    if(this%use_scaling) then
        a_scl = this%compute_scaling_factor(this%exact_z(1),CMB)
    else
        a_scl = 1
    end if

    if(abs(this%exact_z(1)-PK%y(this%exact_z_index(1)))>1.d-3)then
        write(*,*)'ERROR: MPK redshift does not match the value stored'
        write(*,*)'       in the PK%y array.'
        call MpiStop()
    end if

    do i=1, this%num_mpk_kbands_use
        !Errors from using matter_power_minkh at lower end should be negligible
        k_scaled(i)=max(exp(PK%x(1)),a_scl*this%PKData%mpk_k(i))
        mpk_lin(i)=PK%PowerAt(k_scaled(i),this%exact_z(1))/a_scl**3
    end do


    do_marge = this%Q_Marge
    if (do_marge .and. this%Q_flat) then
        !Marginalize analytically with flat prior on b^2 and b^2*Q
        !as recommended by Max Tegmark for SDSS
        allocate(mpk_k2(this%num_mpk_kbands_use))
        allocate(mpk_WPth_k2(this%num_mpk_points_use))

        mpk_Pth=mpk_lin/(1+this%Ag*k_scaled)
        mpk_k2=mpk_Pth*k_scaled**2
        mpk_WPth = matmul(this%PKData%mpk_W,mpk_Pth)
        mpk_WPth_k2 = matmul(this%PKData%mpk_W,mpk_k2)

        if (allocated(this%PKData%mpk_invcov)) then
            covdat = matmul(this%PKData%mpk_invcov,this%PKData%mpk_P)
            covth = matmul(this%PKData%mpk_invcov,mpk_WPth)
            covth_k2 = matmul(this%PKData%mpk_invcov,mpk_WPth_k2)
        else
            w=1/(this%PKData%mpk_sdev**2)
            covdat = this%PKData%mpk_P*w
            covth = mpk_WPth*w
            covth_k2 = mpk_WPth_k2*w
        end if

        Mat(1,1) = sum(covth*mpk_WPth)
        Mat(2,2) = sum(covth_k2*mpk_WPth_k2)
        Mat(1,2) = sum(covth*mpk_WPth_k2)
        Mat(2,1) = Mat(1,2)
        LnLike = log( Mat(1,1)*Mat(2,2)-Mat(1,2)**2)
        call Matrix_Inverse(Mat)
        vec2(1) = sum(covdat*mpk_WPth)
        vec2(2) = sum(covdat*mpk_WPth_k2)
        LnLike = (sum(this%PKData%mpk_P*covdat) - sum(vec2*matmul(Mat,vec2)) + LnLike ) /2
    else
        if (this%Q_sigma==0) do_marge = .false.

        do iQ=-nQ,nQ
            Q = this%Q_mid +iQ*this%Q_sigma*dQ

            if (this%Q_marge) then
                mpk_Pth=mpk_lin*(1+Q*k_scaled**2)/(1+this%Ag*k_scaled)
            else
                mpk_Pth = mpk_lin
            end if

            mpk_WPth = matmul(this%PKData%mpk_W,mpk_Pth)

            !with analytic marginalization over normalization nuisance (flat prior on b^2)
            !See appendix F of cosmomc paper

            if (allocated(this%PKData%mpk_invcov)) then
                covdat = matmul(this%PKData%mpk_invcov,this%PKData%mpk_P)
                covth = matmul(this%PKData%mpk_invcov,mpk_WPth)
                normV = sum(mpk_WPth*covth)
                chisq(iQ) = sum(this%PKData%mpk_P*covdat)  - sum(mpk_WPth*covdat)**2/normV  + log(normV)
            else
                !with analytic marginalization over normalization nuisance (flat prior on b^2)
                w=1/(this%PKData%mpk_sdev**2)
                normV = sum(mpk_WPth*mpk_WPth*w)
                tmp=sum(mpk_WPth*this%PKData%mpk_P*w)/normV ! avoid subtracting one large number from another
                chisq(iQ) = sum(this%PKData%mpk_P*(this%PKData%mpk_P - mpk_WPth*tmp)*w)  + log(normV)
            end if

            if (do_marge) then
                calweights(iQ) = exp(-(iQ*dQ)**2/2)
            else
                LnLike = chisq(iQ)/2
                exit
            end if
        end do

        !without analytic marginalization
        !! chisq = sum((this%PKData%mpk_P(:) - mpk_WPth(:))**2*w) ! uncommented for debugging purposes

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

    if (Feedback>1) write(*,*)trim(this%name)//' MPK Likelihood = ', LnLike

    if (LnLike > 1e8) then
        write(*,*)'WARNING: '//trim(this%name)//' MPK Likelihood is huge!'
        write(*,*)'          Maybe there is a problem? Likelihood = ',LnLike
    end if

    MPK_LnLike=LnLike

    end function MPK_LnLike

    end module mpk
