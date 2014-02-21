    module CosmoTheory
    use settings
    use cmbtypes
    use GeneralTypes
    use likelihood
    use Interpolation, only : spline, SPLINE_DANGLE
    implicit none
    
    !JD work in progress
    !Type TCosmoTheoryPK
    !    integer   ::  num_k
    !    real(mcp), dimension(:), allocatable :: log_kh
    !    real(mcp), dimension(:), allocatable :: redshifts
    !    real(mcp), dimension(:,:), allocatable :: PK, ddPK
    !    logical :: islog
    !contains
    !MPK stuff
    !procedure :: IOTheory_GetSplines
    !procedure :: InitPK
    !procedure :: MatterPowerAt_zbin
    !procedure :: MatterPowerAt
    !procedure :: MatterPowerAt_Z
    !end Type TCosmoTheoryPK
    
    Type, extends(TTheoryPredictions) :: TCosmoTheoryPredictions
        real(mcp) cl(lmax,num_cls_tot)
        !TT, TE, EE (BB) + other C_l (e.g. lensing)  in that order
        real(mcp) sigma_8
        real(mcp) tensor_ratio_r10, tensor_ratio_02
        integer numderived
        real(mcp) derived_parameters(max_derived_parameters)

        !everything is a function of k/h
        integer   ::  num_k
        real(mcp), dimension(:), allocatable :: log_kh
        real(mcp), dimension(:), allocatable :: redshifts
        !matpower is log(P_k)
        real(mcp), dimension(:,:), allocatable :: matter_power, ddmatter_power
        real(mcp), dimension(:,:), allocatable :: nlmatter_power, ddnlmatter_power
        !JD work in progress
        !type(TCosmoTheoryPK) MPK
        !type(TCosmoTheoryPK) NL_MPK
        
        
    contains
    procedure :: ClsFromTheoryData
    procedure :: WriteTextCls
    !MPK stuff
    procedure :: IOTheory_GetSplines
    procedure :: InitPK
    procedure :: MatterPowerAt_zbin
    procedure :: MatterPowerAt
    procedure :: MatterPowerAt_Z
    !Inherited overrides
    procedure :: WriteTheory => TCosmoTheoryPredictions_WriteTheory
    procedure :: ReadTheory => TCosmoTheoryPredictions_ReadTheory
    procedure :: WriteBestFitData => TCosmoTheoryPredictions_WriteBestFitData
    end Type TCosmoTheoryPredictions

    contains
    
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
    character(LEN=80) fmt
    integer unit

    unit=CreateNewTxtFile(aname)
    call T%ClsFromTheoryData(Cls)
    fmt = concat('(1I6,',num_cls_tot,'E15.5)')
    do l = 2, lmax
        nm = 2*pi/(l*(l+1))
        if (num_cls_ext > 0) then
            write (unit,fmt) l, cls(l,1:num_cls)/nm, cls(l,num_cls+1:num_cls_tot)
        else
            write (unit,fmt) l, cls(l,:)/nm
        end if
    end do
    close(unit)

    end subroutine WriteTextCls

    subroutine InitPK(Theory, num_k, num_z)
    class(TCosmoTheoryPredictions) :: Theory
    integer, intent(in) :: num_k, num_z

    if(allocated(Theory%log_kh))deallocate(Theory%log_kh)
    if(allocated(Theory%matter_power))deallocate(Theory%matter_power)
    if(allocated(Theory%ddmatter_power))deallocate(Theory%ddmatter_power)
    if(allocated(Theory%redshifts))deallocate(Theory%redshifts)
    allocate(Theory%log_kh(num_k))
    allocate(Theory%matter_power(num_k,num_z))
    allocate(Theory%ddmatter_power(num_k,num_z))
    allocate(Theory%redshifts(num_z))
    if(use_nonlinear) then
        if(allocated(Theory%nlmatter_power))deallocate(Theory%nlmatter_power)
        if(allocated(Theory%ddnlmatter_power))deallocate(Theory%ddnlmatter_power)
        allocate(Theory%nlmatter_power(num_k,num_z))
        allocate(Theory%ddnlmatter_power(num_k,num_z))
    end if

    end subroutine InitPK

    subroutine IOTheory_GetSplines(Theory)
    class(TCosmoTheoryPredictions) Theory
    integer :: zix,num_k,nz

    num_k = Theory%num_k
    nz = size(Theory%redshifts)

    do zix=1, nz
        call spline(Theory%log_kh,Theory%matter_power(:,zix),num_k,SPLINE_DANGLE,&
        SPLINE_DANGLE,Theory%ddmatter_power(:,zix))

        if(use_nonlinear .and. allocated(Theory%nlmatter_power))&
        call spline(Theory%log_kh,Theory%nlmatter_power(:,zix),num_k,SPLINE_DANGLE,&
        SPLINE_DANGLE,Theory%ddnlmatter_power(:,zix))
    end do

    end subroutine IOTheory_GetSplines
    
    function MatterPowerAt_zbin(Theory, kh, itf, NNL) result(outpower)
    !Get matter power spectrum at particular k/h by interpolation
    class(TCosmoTheoryPredictions) :: Theory
    integer, intent(in) :: itf
    real (mcp), intent(in) :: kh
    logical, optional, intent(in) :: NNL
    logical :: NL
    real(mcp) :: logk
    integer klo,khi
    real(mcp) outpower, dp
    real(mcp) ho,a0,b0
    real(mcp), dimension(2) :: matpower, ddmat
    integer, save :: i_last = 1

    if(.not. allocated(Theory%log_kh)) then
        write(*,*) 'MPK arrays are not initialized:' 
        write(*,*) 'Make sure you are calling SetPk and filling your power spectra'
        call MPIstop()
    end if
    
    if(present(NNL))then
        NL = NNL
    else
        NL = .false.
    end if

    !AL commenting as not currrently defined
    !if(NL .and. .not. allocated(Theory%nlmatter_power)) then
    !    write(*,*)"You are asking for a nonlinear MPK without having initialized nlmatter_power"
    !    write(*,*)"Most likely you are doing importance sampling and need to turn on redo_pk"
    !    call MPIstop()
    !end if

    logk = log(kh)
    if (logk < Theory%log_kh(1)) then
        if( NL ) then
            matpower=Theory%nlmatter_power(1:2,itf)
        else
            matpower=Theory%matter_power(1:2,itf)
        end if
        dp = (matpower(2)-matpower(1))/(Theory%log_kh(2)-Theory%log_kh(1))
        outpower = matpower(1) + dp*(logk-Theory%log_kh(1))
    else if (logk > Theory%log_kh(Theory%num_k)) then
        !Do dodgy linear extrapolation on assumption accuracy of result won't matter
        if( NL ) then
            matpower=Theory%nlmatter_power(Theory%num_k-1:Theory%num_k,itf)
        else
            matpower=Theory%matter_power(Theory%num_k-1:Theory%num_k,itf)
        end if
        dp = (matpower(2)-matpower(1))/(Theory%log_kh(Theory%num_k)-Theory%log_kh(Theory%num_k-1))
        outpower = matpower(2) + dp*(logk-Theory%log_kh(Theory%num_k))
    else
        klo=min(i_last,Theory%num_k)
        do while (Theory%log_kh(klo) > logk)
            klo=klo-1
        end do
        do while (Theory%log_kh(klo+1)< logk)
            klo=klo+1
        end do
        i_last =klo
        khi=klo+1

        if( NL ) then
            matpower=Theory%nlmatter_power(klo:khi,itf)
            ddmat = Theory%ddnlmatter_power(klo:khi,itf)
        else
            matpower=Theory%matter_power(klo:khi,itf)
            ddmat = Theory%ddmatter_power(klo:khi,itf)
        end if

        ho=Theory%log_kh(khi)-Theory%log_kh(klo)
        a0=(Theory%log_kh(khi)-logk)/ho
        b0=1-a0

        outpower = a0*matpower(1)+b0*matpower(2)+((a0**3-a0)*ddmat(1) &
        + (b0**3-b0)*ddmat(2))*ho**2/6
    end if

    outpower = exp(max(-30._mcp,outpower))

    end function MatterPowerAt_zbin

    function MatterPowerAt(Theory, kh, NNL) result(outpower)
    !Get matter power spectrum at particular k/h by interpolation
    class(TCosmoTheoryPredictions) :: Theory
    real (mcp), intent(in) :: kh
    logical, optional, intent(in) :: NNL
    logical :: NL
    real(mcp) outpower

    if(present(NNL))then
        NL = NNL
    else
        NL = .false.
    end if

    outpower = Theory%MatterPowerAt_zbin(kh,1,NL)

    end function MatterPowerAt

    function MatterPowerAt_Z(Theory, kh, z, NNL) result(outpower)
    !Get matter power spectrum at particular k/h by interpolation
    class(TCosmoTheoryPredictions) :: Theory
    real (mcp), intent(in) :: kh, z
    logical, optional, intent(in) :: NNL
    logical :: NL
    integer zlo, zhi, iz, itf, nz
    real(mcp) outpower
    real(mcp) ho,a0,b0
    real(mcp), dimension(4) :: matpower, ddmat, zvec
    integer, save :: zi_last = 1
    
    if(.not. allocated(Theory%log_kh)) then
        write(*,*) 'MPK arrays are not initialized:' 
        write(*,*) 'Make sure you are calling SetPk and filling your power spectra'
        call MPIstop()
    end if

    nz = size(Theory%redshifts)
    
    if(present(NNL))then
        NL = NNL
    else
        NL = .false.
    end if

    if(z>Theory%redshifts(nz)) then
        write (*,*) ' z out of bounds in MatterPowerAt_Z (',z,')'
        call MPIstop()
    end if

    zlo=min(zi_last,size(Theory%redshifts))
    do while (Theory%redshifts(zlo) > z)
        zlo=zlo-1
    end do
    do while (Theory%redshifts(zlo+1)< z)
        zlo=zlo+1
    end do
    zi_last=zlo
    zhi=zlo+1

    if(zlo==1)then
        zvec = 0
        matpower = 0
        ddmat = 0
        iz = 2
        zvec(2:4)=Theory%redshifts(zlo:zhi+1)
        do itf=zlo, zhi+1
            matpower(iz) = log(Theory%MatterPowerAt_zbin(kh,itf,NL))
            iz=iz+1
        end do
        call spline(zvec(2:4),matpower(2:4),3,SPLINE_DANGLE,SPLINE_DANGLE,ddmat(2:4))
    else if(zhi==nz)then
        zvec = 0
        matpower = 0
        ddmat = 0
        iz = 1
        zvec(1:3)=Theory%redshifts(zlo-1:zhi)
        do itf=zlo-1, zhi
            matpower(iz) = log(Theory%MatterPowerAt_zbin(kh,itf,NL))
            iz=iz+1
        end do
        call spline(zvec(1:3),matpower(1:3),3,SPLINE_DANGLE,SPLINE_DANGLE,ddmat(1:3))
    else
        iz = 1
        zvec(:)=Theory%redshifts(zlo-1:zhi+1)
        do itf=zlo-1, zhi+1
            matpower(iz) = log(Theory%MatterPowerAt_zbin(kh,itf,NL))
            iz=iz+1
        end do
        call spline(zvec,matpower,4,SPLINE_DANGLE,SPLINE_DANGLE,ddmat)
    end if

    ho=zvec(3)-zvec(2)
    a0=(zvec(3)-z)/ho
    b0=(z-zvec(2))/ho

    outpower = a0*matpower(2)+b0*matpower(3)+((a0**3-a0)*ddmat(2) &
    +(b0**3-b0)*ddmat(3))*ho**2/6

    outpower = exp(max(-30._mcp,outpower))

    end function MatterPowerAt_Z
        
    subroutine TCosmoTheoryPredictions_WriteTheory(T, unit)
    Class(TCosmoTheoryPredictions) T
    integer, intent(in) :: unit
    integer, parameter :: varcount = 1
    integer tmp(varcount)
    logical, save :: first = .true.

    if (first .and. new_chains) then
        first = .false.
        write(unit) use_LSS, compute_tensors
        write(unit) lmax, lmax_tensor, num_cls, num_cls_ext
        write(unit) varcount
        if (get_sigma8) then
            tmp(1)=1
        else
            tmp(1)=0
        end if
        write(unit) tmp(1:varcount)
    end if

    write(unit) T%numderived
    write(unit) T%derived_parameters(1:T%numderived)
    write(unit) T%cl(2:lmax,1:num_cls)
    if (num_cls_ext>0) write(unit) T%cl(2:lmax,num_cls+1:num_cls_tot)

    if (compute_tensors) then
        write(unit) T%tensor_ratio_02, T%tensor_ratio_r10
    end if

    if (get_sigma8 .or. use_LSS) write(unit) T%sigma_8

    if (use_LSS) then
        write(unit) T%num_k, size(T%redshifts)
        write(unit) T%log_kh
        write(unit) T%redshifts
        write(unit) T%matter_power
        write(unit) use_nonlinear
        if(use_nonlinear) write(unit) T%nlmatter_power
    end if

    end subroutine TCosmoTheoryPredictions_WriteTheory

    subroutine TCosmoTheoryPredictions_ReadTheory(T, unit)
    Class(TCosmoTheoryPredictions) T
    integer, intent(in) :: unit
    integer unused
    logical, save :: first = .true.
    logical, save :: has_sigma8, has_LSS, has_tensors
    integer, save :: almax, almaxtensor, anumcls, anumclsext, tmp(1)
    logical, save :: planck1_format
    real(mcp) old_matterpower(74,1)
    !JD 10/13 new variables for handling new pk arrays
    integer :: num_z, stat
    logical  has_nonlinear

    if (first) then
        first = .false.
        read(unit) has_LSS, has_tensors
        read(unit) almax, almaxtensor, anumcls, anumclsext
        if (almax > lmax) call MpiStop('ReadTheory: reading file with larger lmax')
        if (anumcls /= num_cls) call MpiStop('ReadTheory: reading file with different Cls')
        if (anumclsext /= num_cls_ext) call MpiStop('ReadTheory: reading file with different ext Cls')
        read(unit) unused
        planck1_format = unused==0
        if (unused>0) read(unit) tmp(1:unused)
        if (.not. planck1_format) has_sigma8 = tmp(1)==1
    end if

    T%cl = 0
    T%derived_parameters=0
    read(unit) T%numderived
    read(unit) T%derived_parameters(1:T%numderived)
    read(unit) T%cl(2:almax,1:anumcls)
    if (anumclsext >0) read(unit) T%cl(2:almax,num_cls+1:num_cls+anumclsext)

    if (has_tensors) then
        read(unit) T%tensor_ratio_02, T%tensor_ratio_r10
    end if

    if (planck1_format) then
        if (has_LSS) then
            read(unit) T%sigma_8, old_matterpower
            !discard unused mpk array
        end if
    else
        if (has_sigma8 .or. has_LSS) read(unit) T%sigma_8
        if (has_LSS) then
            read(unit) T%num_k, num_z
            call T%InitPK(T%num_k,num_z)
            read(unit) T%log_kh
            read(unit) T%redshifts
            read(unit, iostat=stat) T%matter_power
            if (IS_IOSTAT_END(stat) .and. use_nonlinear) then
                deallocate(T%nlmatter_power,T%ddnlmatter_power)
                write(*,*)"ReadTheory:  You want a nonlinear MPK but but your datafile is "
                write(*,*)"in an old format does not include one."
                write(*,*)"Make sure you set redo_pk = T or the program will fail"
            else
                read(unit)has_nonlinear
                if(has_nonlinear .and. use_nonlinear) then
                    read(unit)T%nlmatter_power
                else if(has_nonlinear .and. .not. use_nonlinear) then
                    if(allocated(T%nlmatter_power))deallocate(T%nlmatter_power)
                    if(allocated(T%ddnlmatter_power))deallocate(T%ddnlmatter_power)
                    allocate(T%nlmatter_power(T%num_k,num_z))
                    allocate(T%ddnlmatter_power(T%num_k,num_z))
                    write(*,*)"Your data files have nonlinear power spectra, but you are not using"
                    write(*,*)"nonlinear power spectra.  Be careful that this is what you intended"
                    read(unit)T%nlmatter_power
                else 
                    deallocate(T%nlmatter_power,T%ddnlmatter_power)
                    write(*,*)"ReadTheory:  You want a nonlinear MPK but but your datafile does not include one."
                    write(*,*)"Make sure you set redo_pk = T or the program will fail"
                end if
            end if
            call T%IOTheory_GetSplines()
        end if
    end if

    end subroutine TCosmoTheoryPredictions_ReadTheory

    subroutine TCosmoTheoryPredictions_WriteBestFitData(Theory,fnameroot)
    class(TCosmoTheoryPredictions) Theory
    character(LEN=*), intent(in) :: fnameroot

    if (use_CMB) call Theory%WriteTextCls(fnameroot //'.bestfit_cl')

    end subroutine TCosmoTheoryPredictions_WriteBestFitData

    end module CosmoTheory