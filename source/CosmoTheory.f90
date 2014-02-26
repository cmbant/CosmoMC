    module CosmoTheory
    use settings
    use cmbtypes
    use GeneralTypes
    use likelihood
    use Interpolation, only : spline, SPLINE_DANGLE
    implicit none
    private

    Type TCosmoTheoryPK
        !everything is a function of k/h
        integer   ::  num_k, num_z
        real(mcp), dimension(:), allocatable :: log_kh
        real(mcp), dimension(:), allocatable :: redshifts
        real(mcp), dimension(:,:), allocatable :: matter_power, ddmatter_power
        !whether power is stored as log(P_k)
        logical :: islog
    contains
    !MPK stuff
    procedure :: IOPK_GetSplines
    procedure :: InitPK
    procedure :: ClearPK
    procedure :: PowerAt_zbin
    procedure :: PowerAt
    procedure :: PowerAt_Z
    end Type TCosmoTheoryPK

    Type, extends(TTheoryPredictions) :: TCosmoTheoryPredictions
        real(mcp) cl(lmax,num_cls_tot)
        !TT, TE, EE (BB) + other C_l (e.g. lensing)  in that order
        real(mcp) sigma_8
        real(mcp) tensor_ratio_r10, tensor_ratio_02
        integer numderived
        real(mcp) derived_parameters(max_derived_parameters)

        type(TCosmoTheoryPK), pointer :: MPK
        type(TCosmoTheoryPK), pointer :: NL_MPK
    contains
    procedure :: ClsFromTheoryData
    procedure :: WriteTextCls
    !Inherited overrides
    procedure :: WriteTheory => TCosmoTheoryPredictions_WriteTheory
    procedure :: ReadTheory => TCosmoTheoryPredictions_ReadTheory
    procedure :: WriteBestFitData => TCosmoTheoryPredictions_WriteBestFitData
    end Type TCosmoTheoryPredictions

    public TCosmoTheoryPredictions, TCosmoTheoryPK
    contains

    subroutine InitPK(PK,num_k,num_z,islog)
    class(TCosmoTheoryPK) :: PK
    integer, intent(in) :: num_k, num_z
    logical, intent(in) :: islog

    PK%num_k = num_k
    PK%num_z = num_z
    PK%islog = islog
    call PK%ClearPK()
    allocate(PK%log_kh(num_k))
    allocate(PK%matter_power(num_k,num_z))
    allocate(PK%ddmatter_power(num_k,num_z))
    allocate(PK%redshifts(num_z))

    end subroutine InitPK

    subroutine ClearPK(PK)
    class(TCosmoTheoryPK) :: PK

    if(allocated(PK%log_kh))deallocate(PK%log_kh)
    if(allocated(PK%matter_power))deallocate(PK%matter_power)
    if(allocated(PK%ddmatter_power))deallocate(PK%ddmatter_power)
    if(allocated(PK%redshifts))deallocate(PK%redshifts)

    end subroutine ClearPK

    subroutine IOPK_GetSplines(PK)
    class(TCosmoTheoryPK) PK
    integer :: zix

    do zix=1, PK%num_z
        call spline(PK%log_kh,PK%matter_power(:,zix),PK%num_k,SPLINE_DANGLE,&
        SPLINE_DANGLE,PK%ddmatter_power(:,zix))
    end do

    end subroutine IOPK_GetSplines

    function PowerAt_zbin(PK, kh, itf, isAt_z) result(outpower)
    !Get matter power spectrum at particular k/h by interpolation
    class(TCosmoTheoryPK) :: PK
    integer, intent(in) :: itf
    real (mcp), intent(in) :: kh
    logical, optional, intent(in) :: isAt_z
    logical :: wantforAt_z
    real(mcp) :: logk
    integer klo,khi
    real(mcp) outpower, dp
    real(mcp) ho,a0,b0
    real(mcp), dimension(2) :: matpower, ddmat
    integer, save :: i_last = 1

    if(.not. allocated(PK%log_kh)) then
        write(*,*) 'ERROR:  PowerAt_zbin least one of your MPK arrays is not initialized:'
        write(*,*) '        Make sure you are calling a SetPk and filling your power spectra.'
        write(*,*) '        This error could also mean you are doing importance sampling'
        write(*,*) '        and need to turn on redo_pk.'
        call MPIstop()
    end if

    if(present(isAt_z)) then
        wantforAt_z = isAt_z
    else
        wantforAt_z = .false.
    end if

    logk = log(kh)
    if (logk < PK%log_kh(1)) then
        matpower=PK%matter_power(1:2,itf)
        dp = (matpower(2)-matpower(1))/(PK%log_kh(2)-PK%log_kh(1))
        outpower = matpower(1) + dp*(logk-PK%log_kh(1))
    else if (logk > PK%log_kh(PK%num_k)) then
        !Do dodgy linear extrapolation on assumption accuracy of result won't matter
        matpower=PK%matter_power(PK%num_k-1:PK%num_k,itf)
        dp = (matpower(2)-matpower(1))/(PK%log_kh(PK%num_k)-PK%log_kh(PK%num_k-1))
        outpower = matpower(2) + dp*(logk-PK%log_kh(PK%num_k))
    else
        klo=min(i_last,PK%num_k)
        do while (PK%log_kh(klo) > logk)
            klo=klo-1
        end do
        do while (PK%log_kh(klo+1)< logk)
            klo=klo+1
        end do
        i_last =klo
        khi=klo+1

        matpower=PK%matter_power(klo:khi,itf)
        ddmat = PK%ddmatter_power(klo:khi,itf)

        ho=PK%log_kh(khi)-PK%log_kh(klo)
        a0=(PK%log_kh(khi)-logk)/ho
        b0=(logk-PK%log_kh(klo))/ho

        outpower = a0*matpower(1)+b0*matpower(2)+((a0**3-a0)*ddmat(1) &
        + (b0**3-b0)*ddmat(2))*ho**2/6
    end if

    if (PK%islog .and. .not. wantforAt_z) outpower = exp(outpower)

    end function PowerAt_zbin

    function PowerAt(PK, kh) result(outpower)
    !Get matter power spectrum at particular k/h by interpolation
    class(TCosmoTheoryPK) :: PK
    real (mcp), intent(in) :: kh
    real(mcp) outpower

    outpower = PK%PowerAt_zbin(kh,1)

    end function PowerAt

    function PowerAt_Z(PK, kh, z) result(outpower)
    !Get matter power spectrum at particular k/h by interpolation
    class(TCosmoTheoryPK) :: PK
    real (mcp), intent(in) :: kh, z
    integer zlo, zhi, iz, itf
    real(mcp) outpower
    real(mcp) ho,a0,b0
    real(mcp), dimension(4) :: matpower, ddmat, zvec
    integer, save :: zi_last = 1

    if(.not. allocated(PK%log_kh)) then
        write(*,*) 'ERROR:  PowerAt_Z least one of your MPK arrays is not initialized:'
        write(*,*) '        Make sure you are calling a SetPk and filling your power spectra.'
        write(*,*) '        This error could also mean you are doing importance sampling'
        write(*,*) '        and need to turn on redo_pk.'
        call MPIstop()
    end if

    if(z>PK%redshifts(PK%num_z)) then
        write (*,*) ' z out of bounds in PowerAt_Z (',z,')'
        call MPIstop()
    end if

    zlo=min(zi_last,PK%num_z)
    do while (PK%redshifts(zlo) > z)
        zlo=zlo-1
    end do
    do while (PK%redshifts(zlo+1)< z)
        zlo=zlo+1
    end do
    zi_last=zlo
    zhi=zlo+1

    if(zlo==1)then
        zvec = 0
        matpower = 0
        ddmat = 0
        iz = 2
        zvec(2:4)=PK%redshifts(zlo:zhi+1)
        do itf=zlo, zhi+1
            matpower(iz) = PK%PowerAt_zbin(kh,itf,.true.)
            iz=iz+1
        end do
        call spline(zvec(2:4),matpower(2:4),3,SPLINE_DANGLE,SPLINE_DANGLE,ddmat(2:4))
    else if(zhi==PK%num_z)then
        zvec = 0
        matpower = 0
        ddmat = 0
        iz = 1
        zvec(1:3)=PK%redshifts(zlo-1:zhi)
        do itf=zlo-1, zhi
            matpower(iz) = PK%PowerAt_zbin(kh,itf,.true.)
            iz=iz+1
        end do
        call spline(zvec(1:3),matpower(1:3),3,SPLINE_DANGLE,SPLINE_DANGLE,ddmat(1:3))
    else
        iz = 1
        zvec(:)=PK%redshifts(zlo-1:zhi+1)
        do itf=zlo-1, zhi+1
            matpower(iz) = PK%PowerAt_zbin(kh,itf,.true.)
            iz=iz+1
        end do
        call spline(zvec,matpower,4,SPLINE_DANGLE,SPLINE_DANGLE,ddmat)
    end if

    ho=zvec(3)-zvec(2)
    a0=(zvec(3)-z)/ho
    b0=(z-zvec(2))/ho

    outpower = a0*matpower(2)+b0*matpower(3)+((a0**3-a0)*ddmat(2) &
    +(b0**3-b0)*ddmat(3))*ho**2/6

    if(PK%islog) outpower = exp(outpower)

    end function PowerAt_Z

    subroutine ClsFromTheoryData(T, Cls)
    class(TCosmoTheoryPredictions) T
    real(mcp) Cls(lmax,num_cls_tot)

    Cls(2:lmax,1:num_clsS) =T%cl(2:lmax,1:num_clsS)
    if (num_cls>3 .and. num_ClsS==3) Cls(2:lmax,num_cls)=0
    if (num_cls_ext>0) then
        Cls(2:lmax,num_cls+1:num_cls_tot) =T%cl(2:lmax,num_clsS+1:num_clsS+num_cls_ext)
    end if

    end subroutine ClsFromTheoryData

    subroutine WriteTextCls(T,aname)
    class(TCosmoTheoryPredictions) T
    character (LEN=*), intent(in) :: aname
    integer l
    real(mcp) Cls(lmax,num_cls_tot), nm
    character(LEN=*), parameter :: fmt = '(1I6,*(E15.5))'
    Type(TTextFile) :: F

    call F%CreateFile(aname)
    call T%ClsFromTheoryData(Cls)
    do l = 2, lmax
        nm = 2*pi/(l*(l+1))
        if (num_cls_ext > 0) then
            write (F%unit,fmt) l, cls(l,1:num_cls)/nm, cls(l,num_cls+1:num_cls_tot)
        else
            write (F%unit,fmt) l, cls(l,:)/nm
        end if
    end do
    call F%Close()

    end subroutine WriteTextCls

    subroutine TCosmoTheoryPredictions_WriteTheory(T, F)
    Class(TCosmoTheoryPredictions) T
    class(TFileStream) :: F
    integer, parameter :: varcount = 0
    integer tmp(varcount)
    logical, save :: first = .true.

    if (first .and. new_chains) then
        first = .false.
        write(F%unit) use_LSS, compute_tensors, get_sigma8
        write(F%unit) lmax, lmax_tensor, num_cls, num_cls_ext
        write(F%unit) varcount
        write(F%unit) tmp(1:varcount)
    end if

    write(F%unit) T%numderived
    write(F%unit) T%derived_parameters(1:T%numderived)
    write(F%unit) T%cl(2:lmax,1:num_cls)
    if (num_cls_ext>0) write(F%unit) T%cl(2:lmax,num_cls+1:num_cls_tot)

    if (compute_tensors) then
        write(F%unit) T%tensor_ratio_02, T%tensor_ratio_r10
    end if

    if (get_sigma8 .or. use_LSS) write(F%unit) T%sigma_8

    if (use_LSS) then
        write(F%unit) T%MPK%num_k, T%MPK%num_z
        write(F%unit) T%MPK%log_kh
        write(F%unit) T%MPK%redshifts
        write(F%unit) T%MPK%matter_power
        write(F%unit) use_nonlinear
        if(use_nonlinear) write(F%unit) T%NL_MPK%matter_power
    end if

    end subroutine TCosmoTheoryPredictions_WriteTheory

    subroutine TCosmoTheoryPredictions_ReadTheory(T, F)
    Class(TCosmoTheoryPredictions) T
    class(TFileStream) :: F
    integer unused
    logical, save :: first = .true.
    logical, save :: has_sigma8, has_LSS, has_tensors
    integer, save :: almax, almaxtensor, anumcls, anumclsext, tmp(1)
    !JD 02/14 new variables for handling new pk arrays
    integer :: num_k, num_z, stat
    logical  has_nonlinear

    if (first) then
        first = .false.
        read(F%unit) has_LSS, has_tensors, has_sigma8
        read(F%unit) almax, almaxtensor, anumcls, anumclsext
        if (almax > lmax) call MpiStop('ReadTheory: reading file with larger lmax')
        if (anumcls /= num_cls) call MpiStop('ReadTheory: reading file with different Cls')
        if (anumclsext /= num_cls_ext) call MpiStop('ReadTheory: reading file with different ext Cls')
        read(F%unit) unused
        read(F%unit) tmp(1:unused)
    end if

    T%cl = 0
    T%derived_parameters=0
    read(F%unit) T%numderived
    read(F%unit) T%derived_parameters(1:T%numderived)
    read(F%unit) T%cl(2:almax,1:anumcls)
    if (anumclsext >0) read(F%unit) T%cl(2:almax,num_cls+1:num_cls+anumclsext)

    if (has_tensors) then
        read(F%unit) T%tensor_ratio_02, T%tensor_ratio_r10
    end if

    if (has_sigma8 .or. has_LSS) read(F%unit) T%sigma_8
    if (has_LSS) then
        read(F%unit) num_k, num_z
        if(.not. associated(T%MPK)) allocate(T%MPK)
        call T%MPK%InitPK(num_k,num_z,.true.)
        read(F%unit) T%MPK%log_kh
        read(F%unit) T%MPK%redshifts
        read(F%unit, iostat=stat) T%MPK%matter_power
        call T%MPK%IOPK_GetSplines()
        read(F%unit)has_nonlinear
        if(has_nonlinear) then
            if(.not. associated(T%NL_MPK)) allocate(T%NL_MPK)
            T%NL_MPK=T%MPK
            read(F%unit)T%NL_MPK%matter_power
            call T%NL_MPK%IOPK_GetSplines()
            if(.not. use_nonlinear) then
                write(*,*)"WARNING:  ReadTheory - Your data files have nonlinear power spectra,"
                write(*,*)"          but you are not using them. Be careful that this"
                write(*,*)"          is what you intended."
            end if
        end if
    end if

    end subroutine TCosmoTheoryPredictions_ReadTheory

    subroutine TCosmoTheoryPredictions_WriteBestFitData(Theory,fnameroot)
    class(TCosmoTheoryPredictions) Theory
    character(LEN=*), intent(in) :: fnameroot

    if (use_CMB) call Theory%WriteTextCls(fnameroot //'.bestfit_cl')

    end subroutine TCosmoTheoryPredictions_WriteBestFitData

    end module CosmoTheory
