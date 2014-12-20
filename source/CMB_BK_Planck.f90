    module BK_planck
    !BICEP, Keck, Planck B mode likelihood
    use CMBlikes
    use CosmologyTypes
    use FileUtils
    private

    real(mcp), parameter :: T_CMB = 2.7255_mcp     ! CMB temperature
    real(mcp), parameter :: h = 6.62606957e-34_mcp ! Planck's constant
    real(mcp), parameter :: kB = 1.3806488e-23_mcp ! Boltzmann constant
    real(mcp), parameter :: Ghz_Kelvin = h/kB*1e9_mcp

    Type TBandpass
        real(mcp), allocatable :: R(:,:)
        real(mcp), allocatable :: dnu(:)
    end Type TBandpass

    Type, extends(TCMBLikes) :: TBK_planck
        real(mcp) :: T_dust = 19.6_mcp
        Type(TBandpass), allocatable :: Bandpasses(:)
    contains
    procedure :: ReadIni => TBK_planck_ReadIni
    procedure :: AddForegrounds => TBK_planck_AddForegrounds
    procedure :: ReadBandpass => TBK_planck_Read_Bandpass
    end Type TBK_planck

    public TBK_planck
    contains

    subroutine TBK_planck_ReadIni(this, Ini)
    class(TBK_planck) :: this
    class(TSettingIni) :: Ini
    character(LEN=:), allocatable :: fname
    integer i

    !Read all standard parameters
    call this%TCMBLikes%ReadIni(Ini)
    this%has_foregrounds = .true.
    !Set up nuisance parameters
    call this%loadParamNames(Ini%ReadFileName('nusiance_params',relative=.true.,NotFoundFail=.true.))

    !Override dust temperature if wanted
    call Ini%Read('T_dust',this%T_dust)

    !Load in the bandpass files for each map
    allocate(this%Bandpasses(this%nmaps_required))
    do i = 1, this%nmaps_required
        fname = Ini%ReadFileName('bandpass['//this%map_order%Item(i)//']',relative = .true., NotFoundFail=.true.)
        call this%ReadBandpass(fname, this%Bandpasses(i))
    end do


    end subroutine TBK_planck_ReadIni

    subroutine TBK_planck_Read_Bandpass(this, fname, Bandpass)
    class(TBK_planck) :: this
    character(LEN=*), intent(in) :: fname
    real(mcp), pointer :: nu(:)
    Type(TBandpass), target :: Bandpass
    integer i, n

    call File%LoadTxt(fname, Bandpass%R, n)
    nu => Bandpass%R(:,1)
    allocate(Bandpass%dnu(n))
    Bandpass%dnu(1) = nu(2) -nu(1)
    do i=2, n-1
        Bandpass%dnu(i) = (nu(i+1)-nu(i-1))/2
    end do
    Bandpass%dnu(n) = nu(n) - nu(n-1)

    end subroutine TBK_planck_Read_Bandpass

    ! Calculates greybody scaling of dust signal defined at 353 GHz
    ! to specified bandpass.
    subroutine DustScaling(beta,Tdust,bandpass,fdust)
    real(mcp), intent(in) :: beta
    real(mcp), intent(in) :: Tdust
    Type(TBandpass), intent(in) :: bandpass
    real(mcp), intent(out) :: fdust
    real(mcp) :: gb_int = 0 ! Integrate greybody scaling.
    real(mcp) :: th_int = 0 ! Integrate thermodynamic temperature conversion.
    real(mcp) :: nu0 = 353  ! Pivot frequency for dust (353 GHz).
    real(mcp) :: gb0          ! Greybody scaling at pivot.
    real(mcp) :: th0          ! Thermodynamic temperature conversion at pivot.

    ! Integrate greybody scaling and thermodynamic temperature conversion
    ! across experimental bandpass.
    gb_int = sum( bandpass%dnu * bandpass%R(:,2) * bandpass%R(:,1)**(3+beta) &
        / (exp(Ghz_Kelvin*bandpass%R(:,1)/Tdust) - 1))
    th_int = sum( bandpass%dnu * bandpass%R(:,2) * bandpass%R(:,1)**4 * exp(Ghz_Kelvin*bandpass%R(:,1)/T_CMB) &
        / (exp(Ghz_Kelvin*bandpass%R(:,1)/T_CMB) - 1)**2)

    ! Calculate values at pivot frequency.
    gb0 = nu0**(3+beta) / (exp(Ghz_Kelvin*nu0/Tdust) - 1)
    th0 = nu0**4 * exp(Ghz_Kelvin*nu0/T_CMB) / (exp(Ghz_Kelvin*nu0/T_CMB) - 1)**2
    ! Calculate dust scaling.
    fdust = (gb_int / gb0) / (th_int / th0)

    end subroutine DustScaling


    subroutine TBK_planck_AddForegrounds(this,Cls,DataParams)
    class(TBK_planck) :: this
    class(TMapCrossPowerSpectrum), intent(inout) :: Cls(:,:)
    real(mcp), intent(in) :: DataParams(:)
    real(mcp) :: Adust, Async, alphadust, beta
    real(mcp), parameter :: fsync_B2=1.027449,fsync_P217=0.572670
    real(mcp), parameter :: fsync_P353=0.585792
    ! real(mcp), parameter :: fdust_B2=0.046010,fdust_P217=0.144359
    ! real(mcp), parameter :: fdust_P353=1.130129
    ! real(mcp), parameter :: fdust(3) = [fdust_B2,fdust_P217,fdust_P353]
    real(mcp) :: fdust(this%nmaps_required)
    real(mcp), parameter :: fsync(3) = [fsync_B2,fsync_P217,fsync_P353]
    real(mcp) :: dust, sync
    integer i,j,l

    Adust = DataParams(1)
    Async = DataParams(2)
    alphadust = DataParams(3)
    beta = DataParams(4)

    do i=1, this%nmaps_required
        call DustScaling(beta,this%T_dust,this%Bandpasses(i),fdust(i))
        write(*,*) "dust scaling ", this%map_order%Item(i), fdust(i)
    end do

    do i=1, this%nmaps_required
        do j=1, i
            associate(CL=> Cls(i,j))
                !dust = fdust(CL%map_i)*fdust(CL%map_j)
                dust = fdust(i)*fdust(j)
                sync = fsync(CL%map_i)*fsync(CL%map_j)
                do l=this%pcl_lmin,this%pcl_lmax
                    CL%CL(l) = CL%CL(l)+ (Adust)*(l/80.)**(alphadust)*dust+(Async)*(l/80.)**(-0.6)*sync
                end do
            end associate
        end do
    end do

    end subroutine TBK_planck_AddForegrounds

    end module BK_planck
