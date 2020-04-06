    module BK_planck
    !BICEP, Keck, Planck B mode likelihood
    use CMBlikes
    use CosmologyTypes
    use FileUtils
    implicit none
    private

    real(mcp), parameter :: T_CMB = 2.72548_mcp    ! CMB temperature
    real(mcp), parameter :: h = 6.62606957e-34_mcp ! Planck's constant
    real(mcp), parameter :: kB = 1.3806488e-23_mcp ! Boltzmann constant
    real(mcp), parameter :: Ghz_Kelvin = h/kB*1e9_mcp

    Type TBandpass
        real(mcp), allocatable :: R(:,:)
        real(mcp), allocatable :: dnu(:)
        real(mcp) :: th_dust, th_sync
        real(mcp) :: nu_bar
    end Type TBandpass

    Type, extends(TCMBLikes) :: TBK_planck
        Type(TBandpass), allocatable :: Bandpasses(:)
        real(mcp) :: fpivot_dust, fpivot_sync
        real(mcp) :: fpivot_dust_decorr(2), fpivot_sync_decorr(2)
        character(LEN=:), allocatable :: lform_dust_decorr, lform_sync_decorr
    contains
    procedure :: ReadIni => TBK_planck_ReadIni
    procedure :: AddForegrounds => TBK_planck_AddForegrounds
    procedure :: ReadBandpass => TBK_planck_Read_Bandpass
    procedure :: Decorrelation
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

    ! Assign pivot frequencies for foreground decorrelation model.
    this%fpivot_dust_decorr(1) = Ini%Read_Double_Array('fpivot_dust_decorr', 1, 217.0_mcp)
    this%fpivot_dust_decorr(2) = Ini%Read_Double_Array('fpivot_dust_decorr', 2, 353.0_mcp)
    this%fpivot_sync_decorr(1) = Ini%Read_Double_Array('fpivot_sync_decorr', 1, 23.0_mcp)
    this%fpivot_sync_decorr(2) = Ini%Read_Double_Array('fpivot_sync_decorr', 2, 33.0_mcp)

    ! Functional form of ell scaling for foreground decorrelation model.
    this%lform_dust_decorr = Ini%Read_String_Default('lform_dust_decorr', 'flat')
    this%lform_sync_decorr = Ini%Read_String_Default('lform_sync_decorr', 'flat')

    ! Load in the bandpass files for each map
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
    Bandpass%dnu(1) = nu(2) - nu(1)
    do i=2, n-1
        Bandpass%dnu(i) = (nu(i+1) - nu(i-1))/2
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

    ! Calculate bandpass center-of-mass (i.e. mean frequency).
    Bandpass%nu_bar = sum(Bandpass%dnu * Bandpass%R(:,1) * Bandpass%R(:,2)) / &
        sum(Bandpass%dnu * Bandpass%R(:,2))

    end subroutine TBK_planck_Read_Bandpass

    ! Calculates greybody scaling of dust signal defined at 353 GHz
    ! to specified bandpass.
    subroutine DustScaling(beta,Tdust,bandpass,nu0,bandcenter_err,fdust)
    real(mcp), intent(in) :: beta
    real(mcp), intent(in) :: Tdust
    Type(TBandpass), intent(in) :: bandpass
    real(mcp), intent(in) :: nu0 ! Pivot frequency
    real(mcp), intent(in) :: bandcenter_err
    real(mcp), intent(out) :: fdust
    real(mcp) :: gb_int  ! Integrate greybody scaling.
    real(mcp) :: gb0     ! Greybody scaling at pivot.
    real(mcp) :: th_err  ! Conversion factor error due to bandcenter error.
    real(mcp) :: gb_err  ! Greybody scaling error due to bandcenter error.

    ! Integrate greybody scaling and thermodynamic temperature conversion
    ! across experimental bandpass.
    gb_int = sum( bandpass%dnu * bandpass%R(:,2) * bandpass%R(:,1)**(3+beta) &
        / (exp(Ghz_Kelvin*bandpass%R(:,1)/Tdust) - 1))

    ! Calculate values at pivot frequency.
    gb0 = nu0**(3+beta) / (exp(Ghz_Kelvin*nu0/Tdust) - 1)

    ! Add correction for band center error
    if (bandcenter_err /= 1.) then
        th_err = (bandcenter_err)**4 * &
            exp(Ghz_Kelvin * bandpass%nu_bar * (bandcenter_err - 1) / T_CMB) * &
            (exp(Ghz_Kelvin * bandpass%nu_bar / T_CMB) - 1)**2 / &
            (exp(Ghz_Kelvin * bandpass%nu_bar * bandcenter_err / T_CMB) - 1)**2
        gb_err = (bandcenter_err)**(3+beta) * &
            (exp(Ghz_Kelvin * bandpass%nu_bar / Tdust) - 1) / &
            (exp(Ghz_Kelvin * bandpass%nu_bar * bandcenter_err / Tdust) - 1)
    else
        th_err = 1.0_mcp
        gb_err = 1.0_mcp
    end if

    ! Calculate dust scaling.
    fdust = (gb_int / gb0) / bandpass%th_dust * (gb_err / th_err)

    end subroutine DustScaling

    ! Calculates power-law scaling of synchrotron signal defined at 150 GHz
    ! to specified bandpass.
    subroutine SyncScaling(beta,bandpass,nu0,bandcenter_err,fsync)
    real(mcp), intent(in) :: beta
    Type(TBandpass), intent(in) :: bandpass
    real(mcp), intent(in) :: nu0 ! Pivot frequency
    real(mcp), intent(in) :: bandcenter_err
    real(mcp), intent(out) :: fsync
    real(mcp) :: pl_int  ! Integrate power-law scaling.
    real(mcp) :: pl0     ! Power-law scaling at pivot.
    real(mcp) :: th_err  ! Conversion factor error due to bandcenter error.
    real(mcp) :: pl_err  ! Power-law scaling error due to bandcenter error.

    ! Integrate power-law scaling and thermodynamic temperature conversion
    ! across experimental bandpass.
    pl_int = sum( bandpass%dnu * bandpass%R(:,2) * bandpass%R(:,1)**(2+beta))

    ! Calculate values at pivot frequency.
    pl0 = nu0**(2+beta)

    ! Add correction for band center error
    if (bandcenter_err /= 1.) then
        th_err = (bandcenter_err)**4 * &
            exp(Ghz_Kelvin * bandpass%nu_bar * (bandcenter_err - 1) / T_CMB) * &
            (exp(Ghz_Kelvin * bandpass%nu_bar / T_CMB) - 1)**2 / &
            (exp(Ghz_Kelvin * bandpass%nu_bar * bandcenter_err / T_CMB) - 1)**2
        pl_err = (bandcenter_err)**(2+beta)
    else
        th_err = 1.0_mcp
        pl_err = 1.0_mcp
    end if

    ! Calculate sync scaling.
    fsync = (pl_int / pl0) / bandpass%th_sync * (pl_err / th_err)

    end subroutine SyncScaling

    ! Calculate factor by which foreground (dust or sync) power is decreased
    ! for a cross-spectrum between two different frequencies.
    subroutine Decorrelation(this, Delta, nu0, nu1, nupivot, l, lform, Deltap)
    class(TBK_planck) :: this
    real(mcp), intent(in) :: Delta, nu0, nu1
    real(mcp), intent(in) :: nupivot(:)
    integer l
    character(len=*), intent(in) :: lform
    real(mcp), intent(out) :: Deltap
    real(mcp) :: lpivot = 80.0_mcp
    real(mcp) :: scl_nu, scl_ell

    ! Decorrelation scales as log^2(nu0/nu1)
    scl_nu = (log(nu0 / nu1)**2) / (log(nupivot(1) / nupivot(2))**2)
    ! Functional form for ell scaling is specified in .dataset file.
    select case (lform)
    case ("flat")
        scl_ell = 1.0_mcp
    case ("lin")
        scl_ell = l / lpivot
    case ("quad")
        scl_ell = (l / lpivot)**2
        case default
        scl_ell = 1.0_mcp
    end select

    ! Even for small cval, correlation can become negative for sufficiently large frequency
    ! difference or ell value (with linear or quadratic scaling).
    ! Following Vansyngel et al, A&A, 603, A62 (2017), we use an exponential function to
    ! remap the correlation coefficient on to the range [0,1].
    ! We symmetrically extend this function to (non-physical) correlation coefficients
    ! greater than 1 -- this is only used for validation tests of the likelihood model.
    ! Value returned corresponds to the "re-mapped" decorrelation parameter, denoted as
    ! $\Delta'_d$ in Appendix F of the BK15 paper (equations F4 and F5)
    if (Delta > 1) then
        ! If using a physical prior for Delta, then this scenario should never happen.
        Deltap = 2.0_mcp - exp(log(2.0_mcp - Delta) * scl_nu * scl_ell)
    else
        ! This is for physically-relevant values of Delta.
        Deltap = exp(log(Delta) * scl_nu * scl_ell)
    end if

    end subroutine Decorrelation

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
    real(mcp) :: EEtoBB_dust, EEtoBB_sync
    real(mcp) :: Delta_dust, Delta_sync    ! Dust/sync decorrelation model parameters
    real(mcp) :: Deltap_dust, Deltap_sync  ! Remapped dust/sync decorrelation
    integer i,j,l
    real(mcp) :: lpivot = 80.0_mcp
    real(mcp) :: bandcenter_err(this%nmaps_required)
    real(mcp) dustpow(this%pcl_lmin:this%pcl_lmax)
    real(mcp) syncpow(this%pcl_lmin:this%pcl_lmax)
    real(mcp) dustsyncpow(this%pcl_lmin:this%pcl_lmax)
    logical :: need_sync_decorr, need_dust_decorr

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
    Delta_dust = DataParams(11)
    Delta_sync = DataParams(12)

    ! Calculate dust and sync scaling for each map.
    do i=1, this%nmaps_required
        !Read and assign values to band center error params
        if (index(this%used_map_order%Item(i) , '95') > 0) then
            bandcenter_err(i)= DataParams(13)+ DataParams(14) + 1.
        else if (index(this%used_map_order%Item(i),'150') > 0) then
            bandcenter_err(i)= DataParams(13)+ DataParams(15) + 1.
        else if (index(this%used_map_order%Item(i), '220') > 0) then
            bandcenter_err(i)= DataParams(13)+ DataParams(16) + 1.
        else
            bandcenter_err(i)= 1.
        end if

        call DustScaling(betadust,Tdust,this%Bandpasses(i),this%fpivot_dust,bandcenter_err(i),fdust(i))
        call SyncScaling(betasync,this%Bandpasses(i),this%fpivot_sync,bandcenter_err(i),fsync(i))
    end do

    ! Calculate dust/sync/corr angular power spectra at pivot frequencies.
    do l=this%pcl_lmin,this%pcl_lmax
        dustpow(l) = Adust * (l / lpivot) ** alphadust
        syncpow(l) = Async * (l / lpivot) ** alphasync
        dustsyncpow(l) = dustsync_corr * sqrt(Adust * Async) * (l / lpivot) ** ((alphadust + alphasync) / 2)
    end do

    ! Only calculate foreground decorrelation if necessary.
    need_dust_decorr = abs(Delta_dust-1) > 1d-5
    need_sync_decorr = abs(Delta_sync-1) > 1d-5

    ! Loop over all auto and cross spectra
    do i=1, this%nmaps_required
        do j=1, i
            CL=> Cls(i,j)

            ! Only add foregrounds to EE or BB.
            if ((CL%theory_i==2 .and. CL%theory_j==2) .or. (CL%theory_i==3 .and. CL%theory_j==3)) then

                ! Calculate dust/sync/corr scaling for this spectrum.
                dust = fdust(i)*fdust(j)
                sync = fsync(i)*fsync(j)
                dustsync = fdust(i)*fsync(j) + fsync(i)*fdust(j)
                if (CL%theory_i==2 .and. CL%theory_j==2) then
                    ! EE spectrum: multiply foregrounds by EE/BB ratio
                    dust = dust * EEtoBB_dust
                    sync = sync * EEtoBB_sync
                    dustsync = dustsync * sqrt(EEtoBB_dust*EEtoBB_sync)
                end if

                do l=this%pcl_lmin,this%pcl_lmax
                    ! Calculate correlation factors for dust and sync.
                    if ((need_dust_decorr) .and. (i /= j)) then
                        call this%Decorrelation(Delta_dust,this%Bandpasses(i)%nu_bar*bandcenter_err(i),&
                            this%Bandpasses(j)%nu_bar*bandcenter_err(j),this%fpivot_dust_decorr, l, &
                            this%lform_dust_decorr, Deltap_dust)
                    else
                        ! No dust decorrelation for auto-spectra.
                        Deltap_dust = 1.0_mcp
                    end if
                    if ((need_sync_decorr) .and. (i /= j)) then
                        call this%Decorrelation(Delta_sync,this%Bandpasses(i)%nu_bar*bandcenter_err(i),&
                            this%Bandpasses(j)%nu_bar*bandcenter_err(j),this%fpivot_sync_decorr, l, &
                            this%lform_sync_decorr, Deltap_sync)
                    else
                        ! No sync decorrelation for auto-spectra.
                        Deltap_sync = 1.0_mcp
                    end if

                    ! Add foreground model to theory spectrum.
                    ! NOTE: Decorrelation is not implemented for the dust/sync correlated component.
                    !       In BK15, we never turned on correlation and decorrelation parameters
                    !       simultaneously.
                    CL%CL(l) = CL%CL(l) + dust * dustpow(l) * Deltap_dust + &
                        sync * syncpow(l) * Deltap_sync + dustsync * dustsyncpow(l)
                end do
            end if
        end do
    end do

    end subroutine TBK_planck_AddForegrounds

    end module BK_planck
