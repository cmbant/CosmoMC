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
        real(mcp) :: th_dust, th_sync
    end Type TBandpass

    Type, extends(TCMBLikes) :: TBK_planck
        Type(TBandpass), allocatable :: Bandpasses(:)
        real(mcp) :: fpivot_dust, fpivot_sync
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

    !Assign foreground pivot frequencies.
    !Defaults are 353 GHz for dust, 23 GHz for sync.
    this%fpivot_dust = Ini%Read_Double('fpivot_dust', 353.0_mcp)
    this%fpivot_sync = Ini%Read_Double('fpivot_sync', 23.0_mcp)

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
    nu0 = this%fpivot_dust
    th0 = nu0**4 * exp(Ghz_Kelvin*nu0/T_CMB) / (exp(Ghz_Kelvin*nu0/T_CMB) - 1)**2
    Bandpass%th_dust = th_int / th0
    nu0 = this%fpivot_sync
    th0 = nu0**4 * exp(Ghz_Kelvin*nu0/T_CMB) / (exp(Ghz_Kelvin*nu0/T_CMB) - 1)**2
    Bandpass%th_sync = th_int / th0

    end subroutine TBK_planck_Read_Bandpass

    ! Calculates greybody scaling of dust signal defined at 353 GHz
    ! to specified bandpass.
    subroutine DustScaling(beta,Tdust,bandpass,nu0,fdust)
    real(mcp), intent(in) :: beta
    real(mcp), intent(in) :: Tdust
    Type(TBandpass), intent(in) :: bandpass
    real(mcp), intent(in) :: nu0 ! Pivot frequency
    real(mcp), intent(out) :: fdust
    real(mcp) :: gb_int  ! Integrate greybody scaling.
    real(mcp) :: gb0     ! Greybody scaling at pivot.

    ! Integrate greybody scaling and thermodynamic temperature conversion
    ! across experimental bandpass.
    gb_int = sum( bandpass%dnu * bandpass%R(:,2) * bandpass%R(:,1)**(3+beta) &
        / (exp(Ghz_Kelvin*bandpass%R(:,1)/Tdust) - 1))

    ! Calculate values at pivot frequency.
    gb0 = nu0**(3+beta) / (exp(Ghz_Kelvin*nu0/Tdust) - 1)

    ! Calculate dust scaling.
    fdust = (gb_int / gb0) / bandpass%th_dust

    end subroutine DustScaling

    ! Calculates power-law scaling of synchrotron signal defined at 150 GHz
    ! to specified bandpass.
    subroutine SyncScaling(beta,bandpass,nu0,fsync)
    real(mcp), intent(in) :: beta
    Type(TBandpass), intent(in) :: bandpass
    real(mcp), intent(in) :: nu0 ! Pivot frequency
    real(mcp), intent(out) :: fsync
    real(mcp) :: pl_int  ! Integrate power-law scaling.
    real(mcp) :: pl0     ! Power-law scaling at pivot.

    ! Integrate power-law scaling and thermodynamic temperature conversion
    ! across experimental bandpass.
    pl_int = sum( bandpass%dnu * bandpass%R(:,2) * bandpass%R(:,1)**(2+beta))

    ! Calculate values at pivot frequency.
    pl0 = nu0**(2+beta)

    ! Calculate sync scaling.
    fsync = (pl_int / pl0) / bandpass%th_sync

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
    real(mcp) dustpow(this%pcl_lmin:this%pcl_lmax)
    real(mcp) syncpow(this%pcl_lmin:this%pcl_lmax)
    real(mcp) dustsyncpow(this%pcl_lmin:this%pcl_lmax)

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
        call DustScaling(betadust,Tdust,this%Bandpasses(i),this%fpivot_dust,fdust(i))
        call SyncScaling(betasync, this%Bandpasses(i),this%fpivot_sync,fsync(i))
    end do
    do l=this%pcl_lmin,this%pcl_lmax
        dustpow(l) = Adust * (l / lpivot) ** alphadust
        syncpow(l) = Async * (l / lpivot) ** alphasync
        dustsyncpow(l) = dustsync_corr * sqrt(Adust * Async) * (l / lpivot) ** ((alphadust + alphasync) / 2)
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
                    CL%CL(l) = CL%CL(l) + dust * dustpow(l) + sync * syncpow(l) + dustsync * dustsyncpow(l)
                end do
            end if
        end do
    end do

    end subroutine TBK_planck_AddForegrounds

    end module BK_planck
