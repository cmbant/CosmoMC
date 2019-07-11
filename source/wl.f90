    ! Module for galaxy weak lensing, galaxy-galaxy and galaxy auto
    ! e.g. for DES 1 YR
    !AL 2018, following exactly the same approximations as in the DES papers
    !(can only use Weyl potential for lensing)
    ! MR 2019 update to use Weyl potential for galaxy-lensing cross

    module wl

    use settings
    use CosmologyTypes
    use CosmoTheory
    use Calculator_Cosmology
    use Likelihood_Cosmology
    use Interpolation

    implicit none

    private

    integer, parameter :: measurement_xip = 1, measurement_xim = 2, &
        measurement_gammat = 3, measurement_wtheta = 4
    character(LEN=Ini_Enumeration_Len), parameter :: measurement_names(4) = &
        [character(Ini_Enumeration_Len):: 'xip', 'xim', 'gammat', 'wtheta']

    logical,  parameter :: WL_timing = .false.

    type, extends(TCosmoCalcLikelihood) :: WLLikelihood
        real(mcp), allocatable :: invcov(:,:)
        integer :: num_z_bins !lensing sources
        integer :: num_gal_bins !lensing galaxies
        integer :: num_theta_bins
        real(mcp), allocatable  :: theta_bins(:), theta_bin_radians(:)
        real(mcp), allocatable :: z_bins(:)
        integer :: num_z_p ! Source distribution p(z,bin)
        Type(TCubicSpline), allocatable :: p_sp(:), pgal_sp(:)
        real(mcp), allocatable, dimension(:) :: z_p
        integer :: nmeasurement_types
        integer, allocatable :: measurement_types(:), num_type(:)
        integer, allocatable :: used_measurement_types(:)
        logical :: want_type(size(measurement_names))
        real, allocatable :: data_selection(:,:,:,:)
        integer :: num_used
        integer, allocatable :: used_indices(:), used_items(:,:)
        integer, allocatable :: bin_pairs(:,:,:)
        real(mcp), allocatable :: corr_data(:,:,:,:)

        real(mcp) :: ah_factor ! factor to rescale covariance
        integer :: intrinsic_alignment_model
        logical :: use_non_linear ! Whether to use non-linear corrections
        logical :: use_weyl !Wether to get lensing directly from the Weyl potential

        real(mcp), private, allocatable :: data_vector(:) !derived based on cuts
        real(mcp), private, allocatable :: corr_theory(:,:,:,:)
        real(mcp), private, allocatable :: ls_bessel(:)
        integer, private, allocatable :: ls_cl(:)
        real(mcp), private, allocatable :: j0s(:,:), j2s(:,:), j4s(:,:)
        integer, private :: first_theta_bin_used
        integer :: lmax = 50000
        real(mcp) :: acc = 1._mcp !accuracy parameter

    contains

    procedure :: LogLike => WL_LnLike
    procedure :: ReadIni => WL_ReadIni
    procedure :: WriteLikelihoodData => WL_WriteLikelihoodData
    procedure, private :: make_vector
    procedure, private :: calc_theory
    procedure, private :: cl2corr
    procedure, private :: init_bessel_integration

    end type WLLikelihood

    integer, parameter :: intrinsic_alignment_none=1, intrinsic_alignment_DES1YR=2
    character(LEN=Ini_Enumeration_Len), parameter :: intrinsic_alignments(2) = &
        [character(Ini_Enumeration_Len)::'none','DES1YR']

    public WLLikelihood, WLLikelihood_Add

    contains

    subroutine WLLikelihood_Add(LikeList, Ini)
    class(TLikelihoodList) :: LikeList
    class(TSettingIni) :: ini
    Type(WLLikelihood), pointer :: this
    Type(TSettingIni) :: DataSets, OverrideSettings
    integer i
    logical :: nonlinear, useweyl

    !Written generally, but currently only supports DES parameters
    if (Ini%Read_Logical('use_WL',.false.)) then
        nonlinear = Ini%Read_Logical('wl_use_non_linear',.true.)
        useweyl = Ini%Read_Logical('wl_use_weyl',.false.)
        call Ini%TagValuesForName('wl_dataset', DataSets)
        do i= 1, DataSets%Count
            call Ini%SettingValuesForTagName('wl_dataset',DataSets%Name(i),OverrideSettings)
            allocate(this)
            this%needs_nonlinear_pk = nonlinear
            this%use_non_linear = nonlinear
            this%use_weyl = useweyl
            call this%ReadDatasetFile(DataSets%Value(i),OverrideSettings)
            call Ini%Read(Ini%NamedKey('wl_dataset_speed',DataSets%Name(i)),this%speed)
            this%LikelihoodType = 'WL'
            this%tag = DataSets%Name(i)
            this%needs_powerspectra = .true.
            this%needs_weylpower = useweyl
            call LikeList%Add(this)
        end do
        if (Feedback>1) write(*,*) 'read WL data sets'
    end if

    end subroutine WLLikelihood_Add

    subroutine WL_ReadIni(this, Ini)
    use MatrixUtils
    class(WLLikelihood) this
    class(TSettingIni) :: Ini
    Type(TTextFile) :: F
    real(mcp), allocatable :: nz_source(:,:)
    character(LEN=:), allocatable :: InLine
    integer bin1, bin2, maxbin, theta_bin
    real(mcp) theta_range(2)
    character(Ini_Enumeration_Len) tp
    integer status
    integer cov_ix
    character(LEN=:), allocatable :: measurements_format
    integer, allocatable :: used_indices(:), used_items(:,:)
    real(mcp), allocatable :: wl_cov(:,:)
    real(mcp) :: theta, dat, x
    integer lastbin1, lastbin2
    integer i, j, b, maxused, this_type
    integer, allocatable :: ls_tmp(:)
    real(mcp), allocatable :: p(:)

    if (Feedback > 0) write (*,*) 'reading WL data set: '//trim(this%name)

    measurements_format = Ini%Read_String('measurements_format',NotFoundFail=.true.)
    IF (measurements_format /= 'DES') call MpiStop('WL: unknown or old measurements_format')

    this%num_z = Ini%Read_Int('nz_wl',100)
    this%max_z = Ini%Read_Double('max_z',0.d0)

    this%num_z_bins = Ini%Read_Int('num_z_bins')
    this%num_gal_bins = Ini%Read_Int('num_gal_bins', 0)
    maxbin = max(this%num_z_bins, this%num_gal_bins)

    this%acc = Ini%Read_Double('acc',this%acc)
    this%lmax = Ini%Read_int('lmax',this%lmax)

    call File%LoadTxt(Ini%ReadRelativeFilename('nz_file'), nz_source)
    this%num_z_p = size(nz_source(:,2)) + 2
    allocate(this%z_p(this%num_z_p))
    this%z_p(1:this%num_z_p-2) = nz_source(:,2)
    !end with zero, and allow range to extend a bit because of marginalizing over mean
    this%z_p(this%num_z_p-1) = 2*this%z_p(this%num_z_p-2) - this%z_p(this%num_z_p-3)
    this%z_p(this%num_z_p) = 3*this%z_p(this%num_z_p-2) - 2*this%z_p(this%num_z_p-3)
    this%max_z = max(this%max_z, maxval(this%z_p))
    allocate(p(this%num_z_p))
    allocate(this%P_sp(this%num_z_bins))
    do i=1,this%num_z_bins
        p(1:this%num_z_p-2) = nz_source(:,4+i-1)
        p(this%num_z_p-1:this%num_z_p) = 0
        call this%P_sp(i)%Init(this%z_p, p)
    end do
    deallocate(nz_source)

    if (this%num_gal_bins > 0) then
        call File%LoadTxt(Ini%ReadRelativeFilename('nz_gal_file'), nz_source)
        if (size(nz_source(:,2)) /= this%num_z_p-2) call MpiStop('wl assumes windows used same bins')
        if (any(nz_source(:,2) /= this%z_p(1:this%num_z_p-2))) &
            call MpiStop('wl assumes windows used same bins')
        allocate(this%Pgal_sp(this%num_gal_bins))
        do i=1,this%num_gal_bins
            p(1:this%num_z_p-2) = nz_source(:,4+i-1)
            p(this%num_z_p-1:this%num_z_p) = 0
            call this%Pgal_sp(i)%Init(this%z_p, p)
        end do
        deallocate(nz_source)
    end if
    deallocate(p)

    call File%LoadTxt(Ini%ReadRelativeFilename('theta_bins_file'), this%theta_bins)
    this%num_theta_bins = Ini%Read_Int('num_theta_bins',size(this%theta_bins))
    if (size(this%theta_bins) /= this%num_theta_bins ) error stop 'size mismatch in theta_bins_file'

    allocate(this%theta_bin_radians(this%num_theta_bins))
    this%theta_bin_radians = this%theta_bins / 60 * pi/ 180
    !Above is workaround for gfortran bug
    !allocate(this%theta_bin_radians, source=this%theta_bins / 60 * pi/ 180)

    this%kmax = Ini%Read_Double('kmax')
    this%ah_factor = Ini%Read_Double('ah_factor',1.0d0)
    call File%LoadTxt(Ini%ReadRelativeFilename('cov_file'),wl_cov)

    this%intrinsic_alignment_model = &
        Ini%Read_Enumeration('intrinsic_alignment_model', &
        intrinsic_alignments,intrinsic_alignment_DES1YR)

    call Ini%Read_Enumeration_List('data_types',measurement_names, &
        this%measurement_types)
    this%nmeasurement_types = size(this%measurement_types)
    if (Ini%HasKey('used_data_types')) then
        call Ini%Read_Enumeration_List('used_data_types',measurement_names, &
            this%used_measurement_types)
    else
        allocate(this%used_measurement_types, source = this%measurement_types)
    end if
    this%want_type = .false.
    this%want_type(this%used_measurement_types) = .true.

    call this%loadParamNames(Ini%ReadRelativeFileName('nuisance_params',NotFoundFail=.true.))

    allocate(this%data_selection(size(measurement_names),maxbin,maxbin,2))
    this%data_selection = -1
    call F%Open(Ini%ReadRelativeFilename('data_selection'))
    do while (F%ReadLineSkipEmptyAndComments(InLine))
        read(InLine, *, iostat=status) tp, bin1, bin2, theta_range(:)
        if (status/= 0) call MpiStop('WL: Error reading data_selection: ' //InLine)
        i = Ini%EnumerationValue(tp, measurement_names)
        if (i==-1) call MpiStop('data_selection has unknown measurement type')
        if (bin1 < 1 .or. bin1 > maxbin .or. bin2<1 .or. bin2 > maxbin) &
            call MpiStop('data_selection: invalid bin')
        if (this%want_type(i)) then
            this%data_selection(i, bin1, bin2, :) = theta_range
        end if
    end do
    call F%Close

    cov_ix = 0
    this%num_used = 0
    maxused = this%nmeasurement_types*this%num_theta_bins*maxbin**2
    allocate(used_indices(maxused))
    allocate(used_items(maxused,4))
    allocate(this%num_type(this%nmeasurement_types))
    allocate(this%bin_pairs(2, maxbin**2, this%nmeasurement_types))
    allocate(this%corr_data(this%num_theta_bins,maxbin, maxbin, this%nmeasurement_types))

    this%first_theta_bin_used = this%num_theta_bins
    do i = 1, this%nmeasurement_types
        this%num_type(i)=0
        this_type = this%measurement_types(i)
        lastbin1= 0
        lastbin2=0
        call F%Open(Ini%ReadRelativeFilename('measurements['// &
            trim(measurement_names(this_type))//']'))
        do while (F%ReadLineSkipEmptyAndComments(InLine))
            read(InLine, *, iostat=status) bin1, bin2, theta_bin, dat
            if (status/= 0) call MpiStop('WL: Error reading measurements ' &
                //measurement_names(this_type))
            if (theta_bin <1 .or. theta_bin > this%num_theta_bins) &
                call MpiStop('WL: invalid theta bin: '//InLine)
            cov_ix = cov_ix + 1
            if (lastbin1/=bin1 .or. lastbin2/=bin2) then
                this%num_type(i) =this%num_type(i)+1
                this%bin_pairs(1, this%num_type(i),i) = bin1
                this%bin_pairs(2, this%num_type(i),i) = bin2
                lastbin1 = bin1
                lastbin2 = bin2
            end if
            this%corr_data(theta_bin, bin1, bin2, i) = dat
            if (this%want_type(this_type)) then
                theta_range = this%data_selection(this_type, bin1, bin2, 1:2)
                theta = this%theta_bins(theta_bin)
                if (theta>=theta_range(1) .and. theta <= theta_range(2)) then
                    this%num_used = this%num_used + 1
                    used_indices(this%num_used) = cov_ix
                    used_items(this%num_used,1) = i
                    used_items(this%num_used,2) = bin1
                    used_items(this%num_used,3) = bin2
                    used_items(this%num_used,4) = theta_bin
                    this%first_theta_bin_used = min(theta_bin, this%first_theta_bin_used)
                end if
            end if
        end do
        call F%Close
    end do
    if (cov_ix /= size(wl_cov, dim=1) .or. &
        cov_ix /= size(wl_cov, dim=2)) call MpiStop('WL: cov size does not match data size')

    allocate(this%used_items, source = used_items(1:this%num_used,:))
    allocate(this%used_indices, source = used_indices(1:this%num_used))

    allocate(this%invcov, source = wl_cov(this%used_indices,this%used_indices))
    call Matrix_Inverse(this%invcov)

    allocate(this%data_vector(this%num_used))
    call this%make_vector(this%corr_data, this%data_vector)

    call this%init_bessel_integration()
    !Get ell for calculating C_L. Linear then log.
    b=0
    allocate(ls_tmp(this%lmax))
    do i=2, 100 -int(4/this%acc), max(1,int(4/this%acc))
        b=b+1
        ls_tmp(b) = i
    end do
    i=0
    do while (ls_tmp(b) < this%lmax)
        b=b+1
        ls_tmp(b) = nint(100*exp(0.1266*i/this%acc))
        i=i+1
    end do
    allocate(this%ls_cl, source = ls_tmp(1:b))

    end subroutine WL_ReadIni

    subroutine WL_WriteLikelihoodData(this,Theory,DataParams, root)

    implicit none

    class(WLLikelihood) :: this
    class(TTheoryPredictions) :: Theory
    real(mcp), intent(in) :: DataParams(:)
    character(LEN=*), intent(in) :: root
    real(mcp), allocatable :: corr_theory(:,:,:,:)
    type(TTextFile) F
    integer :: i, j, k, tp, type_ix

    F%IntegerFormat = '(*(I6))'
    ! create the output file:
    call F%CreateFile( trim(root)//'_'//trim(this%getTag())//'.theory' )
    ! write the header with the comment:
    call F%WriteInLine('#   theta')
    do type_ix = 1, this%nmeasurement_types
        tp = this%measurement_types(type_ix)
        do j = 1, this%num_type(type_ix)
            call F%WriteInLine( '   '//trim(measurement_names(tp))&
                &//trim(integer_to_string( this%bin_pairs(1,j,type_ix) ))&
                &//trim(integer_to_string( this%bin_pairs(2,j,type_ix) )) )
        end do
    end do
    call F%NewLine()
    ! write the theory prediciton:
    do i=1, this%num_theta_bins
        call F%WriteInLine(this%theta_bins(i))
        do type_ix = 1, this%nmeasurement_types
            tp = this%measurement_types(type_ix)
            do j = 1, this%num_type(type_ix)
                call F%WriteInLine( this%corr_theory(i, this%bin_pairs(1,j,type_ix), this%bin_pairs(2,j,type_ix), tp ) )
            end do
        end do
        call F%NewLine()
    end do
    ! close file:
    call F%Close()

    contains

    ! helper to convert number to string:
    function integer_to_string( number )
    implicit none
    integer, intent(in) :: number               !< Input integer number
    character(10)       :: integer_to_string    !< Output string with the number
    write( integer_to_string, '(i10)' ) number
    integer_to_string = TRIM(ADJUSTL( integer_to_string ))
    end function integer_to_string

    end subroutine WL_WriteLikelihoodData

    function WL_LnLike(this, CMB, Theory, DataParams)
    use MatrixUtils
    Class(WLLikelihood) :: this
    Class(CMBParams) CMB
    Class(TCosmoTheoryPredictions), target :: Theory
    real(mcp) :: DataParams(:)
    real(mcp) WL_LnLike
    real(mcp) vec(this%num_used)

    if (allocated(this%corr_theory) ) deallocate(this%corr_theory)
    allocate(this%corr_theory, source = this%corr_data*0)
    call this%calc_theory(CMB,Theory, this%corr_theory, DataParams)
    call this%make_vector(this%corr_theory, vec)

    vec = vec - this%data_vector
    WL_LnLike = Matrix_QuadForm(this%invcov,vec) / 2

    end function WL_LnLike

    subroutine init_bessel_integration(this)
    Class(WLLikelihood) :: this
    real(mcp) dlog, tmp0, tmp2, tmp4, x
    integer n, ix, ell, i, j, ell_last
    integer, allocatable :: ell_sum_min(:), ell_sum_max(:), bigell(:), dell(:)

    !Get array of roughly log-spaced ls_bessel to sample in the C_L
    n = int(500*this%acc)
    dlog = (log(real(this%lmax)) - log(1.))/n
    allocate(dell(n))
    ix = 0
    ell_last = 1
    do i=1, n
        ell = int(exp(i*dlog))
        if (ell /= ell_last) then
            ix = ix+1
            dell(ix) = ell-ell_last
            ell_last = ell
        end if
    end do
    allocate(this%ls_bessel(ix))
    allocate(ell_sum_min(ix), ell_sum_max(ix))
    ell = 2
    do i=1, ix
        this%ls_bessel(i) = (2*ell+dell(i) -1. )/2
        ell_sum_min(i) = ell
        ell_sum_max(i) = ell + dell(i) -1
        ell = ell+ dell(i)
    end do
    !Calculate average of Bessels over each ell bin range
    allocate(this%j0s(size(this%ls_bessel), this%first_theta_bin_used:this%num_theta_bins))
    allocate(this%j2s(size(this%ls_bessel), this%first_theta_bin_used:this%num_theta_bins))
    allocate(this%j4s(size(this%ls_bessel), this%first_theta_bin_used:this%num_theta_bins))

    do i = 1, size(this%ls_bessel)
        do j = this%first_theta_bin_used, this%num_theta_bins
            tmp0 = 0
            tmp2 = 0
            tmp4 = 0
            do ell = ell_sum_min(i),ell_sum_max(i)
                x = ell * this%theta_bin_radians(j)
                tmp0 = tmp0 + ell*Bessel_J0(x)
                tmp2 = tmp2 + ell*Bessel_JN(2,x)
                tmp4 = tmp4 + ell*Bessel_JN(4,x)
            end do
            this%j0s(i,j) = tmp0/(2*pi)
            this%j2s(i,j) = tmp2/(2*pi)
            this%j4s(i,j) = tmp4/(2*pi)
        end do
    end do

    !allocate(this%ls_bessel, &
    !    source = real((/ (i, i = 2, this%lmax, this%dl_bessel) /), mcp))
    !!Precompute bessels
    !allocate(this%j0s(size(this%ls_bessel), this%first_theta_bin_used:this%num_theta_bins))
    !allocate(this%j2s(size(this%ls_bessel), this%first_theta_bin_used:this%num_theta_bins))
    !allocate(this%j4s(size(this%ls_bessel), this%first_theta_bin_used:this%num_theta_bins))
    !do i=this%first_theta_bin_used,this%num_theta_bins
    !    do j=1, size(this%ls_bessel)
    !        x = this%ls_bessel(j)*this%theta_bin_radians(i)
    !        this%j0s(j,i) = this%ls_bessel(j)*Bessel_J0(x)
    !        this%j2s(j,i) = this%ls_bessel(j)*Bessel_JN(2,x)
    !        this%j4s(j,i) = this%ls_bessel(j)*Bessel_JN(4,x)
    !    end do
    !end do
    !this%j0s = this%j0s * this%dl_bessel / (2 * pi)
    !this%j2s = this%j2s * this%dl_bessel / (2 * pi)
    !this%j4s = this%j4s * this%dl_bessel / (2 * pi)

    end subroutine init_bessel_integration

    subroutine make_vector(this, corr, vec)
    Class(WLLikelihood) :: this
    real(mcp), intent(in) :: corr(:,:,:,:)
    real(mcp), intent(out) :: vec(this%num_used)
    integer i
    integer type_ix, f1, f2, theta_bin

    do i=1, this%num_used
        type_ix = this%used_items(i,1)
        f1 = this%used_items(i,2)
        f2 = this%used_items(i,3)
        theta_bin = this%used_items(i,4)
        vec(i) = corr(theta_bin, f1, f2, type_ix)
    end do

    end subroutine make_vector

    subroutine calc_theory(this,CMB,Theory,corrs, DataParams)
    use Interpolation
    use ArrayUtils
    Class(WLLikelihood) :: this
    Class(CMBParams) CMB
    Class(TCosmoTheoryPredictions), target :: Theory
    real(mcp), intent(out) :: corrs(:,:,:,:)
    real(mcp), intent(in) :: DataParams(:)
    type(TCosmoTheoryPK), pointer :: PK, WPK, MWPK
    real(mcp) h, omm
    real(mcp), allocatable :: chis(:), dchis(:), Hs(:), D_growth(:)
    real(mcp) zshift
    real(mcp) Alignment_z(this%num_z_p), fac(this%num_z_p)
    real(mcp), allocatable :: qs(:,:), n_chi(:,:), qgal(:,:)
    real(mcp), allocatable :: cl_kappa(:,:,:)
    real(mcp), allocatable :: cl_w(:,:,:), cl_cross(:,:,:)
    integer xim_index
    integer i,b, ii, j, f1, f2, tp, type_ix, ix
    real(mcp) kh, khmax, khmin, fac2, cltmp

    real(mcp) bin_bias(this%num_gal_bins)
    real(mcp) shear_calibration_parameters(this%num_z_bins)
    real(mcp) intrinsic_alignment_A,  intrinsic_alignment_alpha,&
        intrinsic_alignment_z0
    real(mcp) source_photoz_errors(this%num_z_bins)
    real(mcp) lens_photoz_errors(this%num_gal_bins)

    real(mcp) :: tmparr(size(this%ls_cl))
    real(mcp) :: kharr(this%num_z_p),zarr(this%num_z_p), powers(this%num_z_p), wpowers(this%num_z_p), mwpowers(this%num_z_p), tmp(this%num_z_p), wtmp(this%num_z_p), mwtmp(this%num_z_p)
    real(mcp) :: time

    time= TimerTime()

    bin_bias = DataParams(1:this%num_gal_bins)
    i = this%num_gal_bins
    shear_calibration_parameters = DataParams(i+1:i+this%num_z_bins)
    i = i + this%num_z_bins
    intrinsic_alignment_A=DataParams(i+1)
    intrinsic_alignment_alpha=DataParams(i+2)
    intrinsic_alignment_z0=DataParams(i+3)
    i=i+3
    lens_photoz_errors = DataParams(i+1:i+this%num_gal_bins)
    i = i + this%num_gal_bins
    source_photoz_errors = DataParams(i+1:i+this%num_z_bins)
    if (this%use_non_linear) then
        PK  => Theory%NL_MPK
        if ( this%use_weyl ) then
            WPK  => Theory%NL_MPK_WEYL
            MWPK => Theory%NL_MPK_WEYL_CROSS
        end if
    else
        PK  => Theory%MPK
        if ( this%use_weyl ) then
            WPK  => Theory%MPK_WEYL
            MWPK => Theory%MPK_WEYL_CROSS
        end if
    end if

    h = CMB%H0/100
    omm = CMB%omdm+CMB%omb

    allocate(chis(this%num_z_p), dchis(this%num_z_p))
    call this%Calculator%ComovingRadialDistanceArr(this%z_p, chis, this%num_z_p)
    dchis(1) = (chis(2) + chis(1))/2
    dchis(this%num_z_p) = chis(this%num_z_p) - chis(this%num_z_p-1)
    dchis(2:this%num_z_p-1) = (chis(3:this%num_z_p) - chis(1:this%num_z_p-2))/2
    allocate(Hs(this%num_z_p))
    allocate(D_growth(this%num_z_p))

    !$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i)
    do i =1 , this%num_z_p
        Hs(i) = this%Calculator%Hofz(this%z_p(i))
        D_growth(i)= Theory%MPK%PowerAt(0.01d0,this%z_p(i))
    end do
    !$OMP END PARALLEL DO

    D_growth = sqrt(D_growth/Theory%MPK%PowerAt(0.01d0,0.d0))

    Alignment_z = intrinsic_alignment_A * ((1 + this%z_p) / &
        (1 + intrinsic_alignment_z0)) ** intrinsic_alignment_alpha &
        * 0.0134 / D_growth

    allocate(qs(this%num_z_p, this%num_z_bins))
    allocate(qgal(this%num_z_p, this%num_gal_bins))
    allocate(n_chi(this%num_z_p, this%num_z_bins))

    !Get lensing source qs and galaxy (lens) qgal, including intrinsic alignment model for qs
    !$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(b, zshift)
    do i=1, this%num_z_p
        !Neglecting any change in normalization of n due to shifting (?)
        do b = 1, this%num_z_bins
            zshift = this%z_p(i)- source_photoz_errors(b)
            if (zshift <this%z_p(1) .or. zshift > this%z_p(this%num_z_p)) then
                n_chi(i,b) =0
            else
                n_chi(i,b) = Hs(i)*this%P_sp(b)%Value(zshift)
            end if
        end do
        do b = 1, this%num_gal_bins
            zshift = this%z_p(i)- lens_photoz_errors(b)
            if (zshift <this%z_p(1) .or. zshift > this%z_p(this%num_z_p)) then
                qgal(i,b) =0
            else
                qgal(i,b) = Hs(i)*this%Pgal_sp(b)%Value(zshift)*bin_bias(b)
            end if
        end do
    end do
    !$OMP END PARALLEL DO
    !FIRSTPRIVATE is a workaround for ifort issues on some machines
    !$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i), FIRSTPRIVATE(fac)
    do b = 1, this%num_z_bins
        fac = dchis*n_chi(:,b)
        do i=1, this%num_z_p
            qs(i,b) = dot_product(fac(i:this%num_z_p),(1 - chis(i) / chis(i:this%num_z_p)))
        end do
        if (this%intrinsic_alignment_model == intrinsic_alignment_DES1YR) then
            qs(:,b) = qs(:,b) - Alignment_z * n_chi(:,b) / (chis * (1 + this%z_p) * 3 * h**2 * (1e5 / const_c) ** 2 / 2)
        end if
        if (this%use_weyl) then
            qs(:,b) = qs(:,b) * chis
        else
            qs(:,b) = qs(:,b) * (3/2._mcp * omm * h**2 * (1e5 / const_c) ** 2) * chis * (1 + this%z_p)
        end if
    end do
    !$OMP END PARALLEL DO

    if (WL_timing)  print *, 'time 1', TimerTime() - time
    time= TimerTime()

    !Get C_kappa
    khmin = exp(PK%x(1))
    khmax = exp(PK%x(PK%nx))
    allocate(cl_kappa(size(this%ls_cl),this%num_z_bins,this%num_z_bins))
    allocate(cl_w(size(this%ls_cl),this%num_gal_bins,this%num_gal_bins))
    allocate(cl_cross(size(this%ls_cl),this%num_gal_bins,this%num_z_bins))
    cl_kappa=0
    cl_w=0
    cl_cross=0

    fac = dchis/chis**2
    !FIRSTPRIVATE is a workaround for ifort issues on some machines
    !$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(j,kh, type_ix, tp, f1, f2, cltmp, ii, ix), &
    !$OMP FIRSTPRIVATE(kharr, zarr, powers, wpowers, mwpowers, tmp, wtmp, mwtmp )
    do i=1, size(this%ls_cl)
        ix =0
        do j = 1, this%num_z_p
            kh= (this%ls_cl(i) + 0.5) / chis(j)/h
            if (kh >= khmin .and. kh <= khmax) then
                ix = ix +1
                zarr(ix) = this%z_p(j)
                kharr(ix) = kh
            end if
        end do
        call PK%PowerAtArr (kharr, zarr, ix, powers )
        if ( this%use_weyl ) then
            call WPK%PowerAtArr(kharr, zarr, ix, wpowers)
            call MWPK%PowerAtArr(kharr, zarr, ix, mwpowers)
        end if
        ix=0
        do j = 1, this%num_z_p
            kh = (this%ls_cl(i) + 0.5) / chis(j)/h
            if (kh >= khmin .and. kh <= khmax) then
                ix = ix+1
                tmp(j)  = fac(j)*powers(ix)/h**3
                if ( this%use_weyl ) then
                    wtmp(j)  = fac(j)*wpowers(ix)
                    mwtmp(j) = -fac(j)*mwpowers(ix)
                else
                    wtmp(j)  = tmp(j)
                    mwtmp(j) = tmp(j)
                end if
            else
                tmp(j)   = 0
                wtmp(j)  = 0
                mwtmp(j) = 0
            end if
        end do
        do type_ix = 1, this%nmeasurement_types
            tp = this%measurement_types(type_ix)
            if (tp==measurement_xim .or. .not. this%want_type(tp)) cycle !assume get from xip (cl_kappa)
            do j=1, this%num_type(type_ix)
                f1 = this%bin_pairs(1,j,type_ix)
                f2 = this%bin_pairs(2,j,type_ix)
                if (tp==measurement_xip) then
                    cltmp = 0
                    do ii = 1, this%num_z_p
                        cltmp = cltmp + wtmp(ii)*qs(ii,f1)*qs(ii,f2)
                    end do
                    cl_kappa(i,f1,f2) = cltmp
                else if (tp==measurement_wtheta) then
                    cltmp = 0
                    do ii = 1, this%num_z_p
                        cltmp = cltmp + tmp(ii)*(qgal(ii,f1)*qgal(ii,f2))
                    end do
                    cl_w(i,f1,f2)=cltmp
                else if (tp==measurement_gammat) then
                    cltmp = 0
                    do ii = 1, this%num_z_p
                        cltmp = cltmp + mwtmp(ii)*(qgal(ii,f1)*qs(ii,f2))
                    end do
                    cl_cross(i,f1,f2)=cltmp
                end if
            end do
        end do
    end do
    !$OMP END PARALLEL DO

    if (WL_timing) print *, 'time 2', TimerTime() - time
    time= TimerTime()

    xim_index = IndexOf(measurement_xim,this%measurement_types, &
        this%nmeasurement_types)
    do type_ix = 1, this%nmeasurement_types
        tp = this%measurement_types(type_ix)
        if (tp==measurement_xim .or. .not. this%want_type(tp)) cycle
        !$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(j)
        do j=1, this%num_type(type_ix)
            call this%cl2corr(tp, cl_kappa, cl_w, cl_cross, corrs, type_ix, &
                xim_index,j, shear_calibration_parameters)
        end do
        !$OMP END PARALLEL DO
    end do
    if (WL_timing) print *, 'time 3', TimerTime() - time

    end subroutine calc_theory

    subroutine cl2corr(this, tp, cl_kappa, cl_w, cl_cross, corrs, &
        type_ix, xim_index, j, shear_calibration_parameters)
    !This only a subroutine to work around issues with OPENMP
    class(WLLikelihood) :: this
    integer, intent(in) :: tp, type_ix, xim_index, j
    real(mcp), intent(in) :: shear_calibration_parameters(*)
    real(mcp), intent(in) :: cl_kappa(:,:,:), cl_w(:,:,:), cl_cross(:,:,:)
    real(mcp), intent(inout) :: corrs(:,:,:,:)
    integer f1, f2
    real(mcp) cl_bessels(size(this%ls_bessel)), fac2
    Type(TCubicSpline) :: CL_sp

    !Note j0s, j2s and j4s already contain L Delta_L/(2*pi) factor
    f1 = this%bin_pairs(1,j,type_ix)
    f2 = this%bin_pairs(2,j,type_ix)
    if (tp==measurement_xip) then
        call CL_sp%Init(this%ls_cl, cl_kappa(:,f1,f2), size(this%ls_cl))
        call CL_sp%Array(this%ls_bessel, cl_bessels)
        fac2 =  (1 + shear_calibration_parameters(f1)) &
            * ( 1 + shear_calibration_parameters(f2))
        corrs(this%first_theta_bin_used:,f1,f2,type_ix) = matmul(cl_bessels, this%j0s) *fac2
        if (xim_index/=0) &
            corrs(this%first_theta_bin_used:,f1,f2,xim_index) = matmul(cl_bessels, this%j4s) *fac2
    else if (tp == measurement_gammat) then
        call CL_sp%Init(this%ls_cl, cl_cross(:,f1,f2), size(this%ls_cl))
        call CL_sp%Array(this%ls_bessel, cl_bessels)
        fac2 = ( 1 + shear_calibration_parameters(f2))
        corrs(this%first_theta_bin_used:,f1,f2,type_ix) = matmul(cl_bessels, this%j2s) *fac2
    else if (tp ==  measurement_wtheta) then
        call CL_sp%Init(this%ls_cl, cl_w(:,f1,f2), size(this%ls_cl))
        call CL_sp%Array(this%ls_bessel, cl_bessels)
        corrs(this%first_theta_bin_used:,f1,f2,type_ix) = matmul(cl_bessels, this%j0s)
    end if

    end subroutine cl2corr

    end module wl


