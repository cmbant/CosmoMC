    module BK_planck
    !BICEP, Keck, Planck B mode likelihood
    use CMBlikes
    use CosmologyTypes

    Type, extends(TCMBLikes) :: TBK_planck
    contains
    procedure :: ReadIni => TBK_planck_ReadIni
    procedure :: AddForegrounds => TBK_planck_AddForegrounds
    end Type TBK_planck


    contains

    subroutine TBK_planck_ReadIni(this, Ini)
    class(TBK_planck) :: this
    class(TSettingIni) :: Ini

    call this%TCMBLikes%ReadIni(Ini)
    this%has_foregrounds = .true.
    call this%loadParamNames(Ini%ReadFileName('nusiance_params',relative=.true.,NotFoundFail=.true.))

    end subroutine TBK_planck_ReadIni

    subroutine TBK_planck_AddForegrounds(this,Cls,DataParams)
    class(TBK_planck) :: this
    class(TMapCrossPowerSpectrum), intent(inout) :: Cls(:,:)
    real(mcp), intent(in) :: DataParams(:)
    real(mcp) :: Adust, Async, alphadust
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