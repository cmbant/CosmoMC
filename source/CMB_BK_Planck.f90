    module BK_planck
    !BICEP, Keck, Planck B mode likelihood
    use CMBlikes
    use CosmologyTypes
    use FileUtils
    implicit none
    private

    real(mcp), parameter :: T_CMB = 2.7255_mcp     ! CMB temperature
    real(mcp), parameter :: h = 6.62606957e-34_mcp ! Planck's constant
    real(mcp), parameter :: kB = 1.3806488e-23_mcp ! Boltzmann constant
    real(mcp), parameter :: Ghz_Kelvin = h/kB*1e9_mcp

    Type TBandpass
        real(mcp), allocatable :: R(:,:)
        real(mcp), allocatable :: dnu(:)
        real(mcp) :: th353, th150
    end Type TBandpass

    Type, extends(TCMBLikes) :: TBK_planck
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
    call this%loadParamNames(Ini%ReadFileName('nuisance_params',relative=.true.,NotFoundFail=.true.))

    !Load in the bandpass files for each map
    allocate(this%Bandpasses(this%nmaps_required))
    do i = 1, this%nmaps_required
        fname = Ini%ReadFileName('bandpass['//this%used_map_order%Item(i)//']',relative = .true., NotFoundFail=.true.)
        call this%ReadBandpass(fname, this%Bandpasses(i))
    end do

    end subroutine TBK_planck_ReadIni

    subroutine TBK_planck_Read_Bandpass(this, fname, Bandpass)
    class(TBK_planck) :: this
    character(LEN=*), intent(in) :: fname
    real(mcp), pointer :: nu(:)
    Type(TBandpass), target :: Bandpass
    integer i, n
    real(mcp) :: th_int, nu0, th0

    call File%LoadTxt(fname, Bandpass%R, n)
    nu => Bandpass%R(:,1)
    allocate(Bandpass%dnu(n))
    Bandpass%dnu(1) = nu(2) -nu(1)
    do i=2, n-1
        Bandpass%dnu(i) = (nu(i+1)-nu(i-1))/2
    end do
    Bandpass%dnu(n) = nu(n) - nu(n-1)

    ! Calculate thermodynamic temperature conversion between this bandpass
    ! and pivot frequencies 353 GHz (used for dust) and 150 GHz (used for
    ! sync).
    th_int = sum( Bandpass%dnu * Bandpass%R(:,2) * Bandpass%R(:,1)**4 * exp(Ghz_Kelvin*bandpass%R(:,1)/T_CMB) &
        / (exp(Ghz_Kelvin*bandpass%R(:,1)/T_CMB) - 1)**2)
    nu0=353
    th0 = nu0**4 * exp(Ghz_Kelvin*nu0/T_CMB) / (exp(Ghz_Kelvin*nu0/T_CMB) - 1)**2
    Bandpass%th353 = th_int / th0
    nu0=150
    th0 = nu0**4 * exp(Ghz_Kelvin*nu0/T_CMB) / (exp(Ghz_Kelvin*nu0/T_CMB) - 1)**2
    Bandpass%th150 = th_int / th0

    end subroutine TBK_planck_Read_Bandpass

    ! Calculates greybody scaling of dust signal defined at 353 GHz
    ! to specified bandpass.
    subroutine DustScaling(beta,Tdust,bandpass,fdust)
    real(mcp), intent(in) :: beta
    real(mcp), intent(in) :: Tdust
    Type(TBandpass), intent(in) :: bandpass
    real(mcp), intent(out) :: fdust
    real(mcp), parameter :: nu0 = 353._mcp  ! Pivot frequency for dust (353 GHz).
    real(mcp) :: gb_int  ! Integrate greybody scaling.
    real(mcp) :: gb0     ! Greybody scaling at pivot.

    ! Integrate greybody scaling and thermodynamic temperature conversion
    ! across experimental bandpass.
    gb_int = sum( bandpass%dnu * bandpass%R(:,2) * bandpass%R(:,1)**(3+beta) &
        / (exp(Ghz_Kelvin*bandpass%R(:,1)/Tdust) - 1))

    ! Calculate values at pivot frequency.
    gb0 = nu0**(3+beta) / (exp(Ghz_Kelvin*nu0/Tdust) - 1)

    ! Calculate dust scaling.
    fdust = (gb_int / gb0) / bandpass%th353

    end subroutine DustScaling

    ! Calculates power-law scaling of synchrotron signal defined at 150 GHz
    ! to specified bandpass.
    subroutine SyncScaling(beta,bandpass,fsync)
    real(mcp), intent(in) :: beta
    Type(TBandpass), intent(in) :: bandpass
    real(mcp), intent(out) :: fsync
    real(mcp), parameter :: nu0 = 150.0_mcp ! Pivot frequency for sync (150 GHz).
    real(mcp) :: pl_int  ! Integrate power-law scaling.
    real(mcp) :: pl0     ! Power-law scaling at pivot.

    ! Integrate power-law scaling and thermodynamic temperature conversion
    ! across experimental bandpass.
    pl_int = sum( bandpass%dnu * bandpass%R(:,2) * bandpass%R(:,1)**(2+beta))

    ! Calculate values at pivot frequency.
    pl0 = nu0**(2+beta)

    ! Calculate dust scaling.
    fsync = (pl_int / pl0) / bandpass%th150

    end subroutine SyncScaling

    subroutine TBK_planck_AddForegrounds(this,Cls,DataParams)
    class(TBK_planck) :: this
    class(TMapCrossPowerSpectrum), target, intent(inout) :: Cls(:,:)
    class(TMapCrossPowerSpectrum), pointer :: CL
    real(mcp), intent(in) :: DataParams(:)
    real(mcp) :: Adust, Async, alphadust, betadust, Tdust
    real(mcp) :: alphasync, betasync, dustsync_corr
    real(mcp) :: fdust(this%nmaps_required)
    real(mcp) :: fsync(this%nmaps_required)
    real(mcp) :: dust, sync, dustsync
    real(mcp) :: EEtoBB_dust,EEtoBB_sync
    integer i,j,l
    real(mcp) :: lpivot = 80.0_mcp

    Adust = DataParams(1)
    Async = DataParams(2)
    alphadust = DataParams(3)
    betadust = DataParams(4)
    Tdust = DataParams(5)
    alphasync = DataParams(6)
    betasync = DataParams(7)
    dustsync_corr = DataParams(8)
    EEtoBB_dust = DataParams(9)
    EEtoBB_sync = DataParams(10)

    do i=1, this%nmaps_required
        call DustScaling(betadust,Tdust,this%Bandpasses(i),fdust(i))
        call SyncScaling(betasync, this%Bandpasses(i), fsync(i))
    end do

    do i=1, this%nmaps_required
        do j=1, i
            CL=> Cls(i,j)
            dust = fdust(i)*fdust(j)
            sync = fsync(i)*fsync(j)
            dustsync = fdust(i)*fsync(j) + fsync(i)*fdust(j)
            If (CL%theory_i==2 .and. CL%theory_j==2) then
                ! EE spectrum: multiply foregrounds by EE/BB ratio
                dust = dust * EEtoBB_dust
                sync = sync * EEtoBB_sync
                dustsync = dustsync * sqrt(EEtoBB_dust*EEtoBB_sync)
            end if

            if ((CL%theory_i==2 .and. CL%theory_j==2) .or. (CL%theory_i==3 .and. CL%theory_j==3)) then
                ! Only add foregrounds to EE or BB.
                do l=this%pcl_lmin,this%pcl_lmax
                    CL%CL(l) = CL%CL(l) + &
                        dust*Adust*(l/lpivot)**(alphadust) + &
                        sync*Async*(l/lpivot)**(alphasync) + &
                        dustsync_corr*dustsync*sqrt(Adust*Async)*(l/lpivot)**((alphadust+alphasync)/2)
                end do
            end if
        end do
    end do

    end subroutine TBK_planck_AddForegrounds

    end module BK_planck
