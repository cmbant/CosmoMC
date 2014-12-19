    module BK_planck
    !BICEP, Keck, Planck B mode likelihood
    use CMBlikes
    use CosmologyTypes
    use FileUtils
    private

    real(mcp), parameter :: T_CMB = 2.7255_mcp

    Type TBandpass
        real(mcp), allocatable :: R(:,:) 
        real(mcp), allocatable :: dnu(:)
    end Type TBandpass

    Type, extends(TCMBLikes) :: TBK_planck
        real(mcp) :: T_dust = 19.6
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


    subroutine TBK_planck_AddForegrounds(this,Cls,DataParams)
    class(TBK_planck) :: this
    class(TMapCrossPowerSpectrum), intent(inout) :: Cls(:,:)
    real(mcp), intent(in) :: DataParams(:)
    real(mcp) :: Adust, Async, alphadust, beta
    real(mcp), parameter :: fsync_B2=1.027449,fsync_P217=0.572670
    real(mcp), parameter :: fsync_P353=0.585792
    real(mcp), parameter :: fdust_B2=0.046010,fdust_P217=0.144359
    real(mcp), parameter :: fdust_P353=1.130129
    real(mcp), parameter :: fdust(3) = [fdust_B2,fdust_P217,fdust_P353]
    real(mcp), parameter :: fsync(3) = [fsync_B2,fsync_P217,fsync_P353]
    real(mcp) :: dust, sync
    integer i,j,l

    Adust = DataParams(1)
    Async = DataParams(2)
    alphadust = DataParams(3)
    beta = DataParams(4)
    
    !TODO map beta into acutal dust, sync deleting above hard coded

    do i=1, this%nmaps_required
        do j=1, i
            associate(CL=> Cls(i,j))
                dust = fdust(CL%map_i)*fdust(CL%map_j)
                sync = fsync(CL%map_i)*fsync(CL%map_j)
                do l=this%pcl_lmin,this%pcl_lmax
                    CL%CL(l) = CL%CL(l)+ (Adust)*(l/80.)**(alphadust)*dust+(Async)*(l/80.)**(-0.6)*sync
                end do
            end associate
        end do
    end do

    end subroutine TBK_planck_AddForegrounds

    end module BK_planck