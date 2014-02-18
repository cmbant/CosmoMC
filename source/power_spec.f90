    module powerspec
    use settings
    use cmbtypes
    use precision
    use Interpolation, only : spline, SPLINE_DANGLE
    implicit none

    contains

    function MatterPowerAt_zbin(PK, kh, itf, NNL) result(outpower)
    !Get matter power spectrum at particular k/h by interpolation
    class(TCosmoTheoryPredictions) :: PK
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

    if(.not. allocated(Pk%log_kh) then
        write(*,*) 'MPK arrays are not initialized:' 
        write(*,*) 'Make sure you are calling SetPk and filling your power spectra'
        call MPIstop()
    end if
    
    if(present(NNL))then
        NL = NNL
    else
        NL = .false.
    end if
    
    if(NL .and. .not allocated(Theory%nlmatter_power)) then
        write(*,*)"You are asking for a nonlinear MPK without having initialized nlmatter_power"
        write(*,*)"Most likely you are doing importance sampling and need to turn on redo_pk"
        call MPIstop()
    end if

    logk = log(kh)
    if (logk < PK%log_kh(1)) then
        if( NL ) then
            matpower=PK%nlmatter_power(1:2,itf)
        else
            matpower=PK%matter_power(1:2,itf)
        end if
        dp = (matpower(2)-matpower(1))/(PK%log_kh(2)-PK%log_kh(1))
        outpower = matpower(1) + dp*(logk-PK%log_kh(1))
    else if (logk > PK%log_kh(PK%num_k)) then
        !Do dodgy linear extrapolation on assumption accuracy of result won't matter
        if( NL ) then
            matpower=PK%nlmatter_power(PK%num_k-1:PK%num_k,itf)
        else
            matpower=PK%matter_power(PK%num_k-1:PK%num_k,itf)
        end if
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

        if( NL ) then
            matpower=PK%nlmatter_power(klo:khi,itf)
            ddmat = PK%ddnlmatter_power(klo:khi,itf)
        else
            matpower=PK%matter_power(klo:khi,itf)
            ddmat = PK%ddmatter_power(klo:khi,itf)
        end if

        ho=PK%log_kh(khi)-PK%log_kh(klo)
        a0=(PK%log_kh(khi)-logk)/ho
        b0=1-a0

        outpower = a0*matpower(1)+b0*matpower(2)+((a0**3-a0)*ddmat(1) &
        + (b0**3-b0)*ddmat(2))*ho**2/6
    end if

    outpower = exp(max(-30._dl,outpower))

    end function MatterPowerAt_zbin

    function MatterPowerAt(PK, kh, NNL) result(outpower)
    !Get matter power spectrum at particular k/h by interpolation
    class(TCosmoTheoryPredictions) :: PK
    real (mcp), intent(in) :: kh
    logical, optional, intent(in) :: NNL
    logical :: NL
    real(mcp) outpower

    if(present(NNL))then
        NL = NNL
    else
        NL = .false.
    end if

    outpower = MatterPowerAt_zbin(PK,kh,1,NL)

    end function MatterPowerAt

    function MatterPowerAt_Z(PK, kh, z, NNL) result(outpower)
    !Get matter power spectrum at particular k/h by interpolation
    class(TCosmoTheoryPredictions) :: PK
    real (mcp), intent(in) :: kh, z
    logical, optional, intent(in) :: NNL
    logical :: NL
    integer zlo, zhi, iz, itf, nz
    real(mcp) outpower
    real(mcp) ho,a0,b0
    real(mcp), dimension(4) :: matpower, ddmat, zvec
    integer, save :: zi_last = 1
    
    if(.not. allocated(Pk%log_kh) then
        write(*,*) 'MPK arrays are not initialized:' 
        write(*,*) 'Make sure you are calling SetPk and filling your power spectra'
        call MPIstop()
    end if

    nz = size(PK%redshifts)
    
    if(present(NNL))then
        NL = NNL
    else
        NL = .false.
    end if

    if(z>PK%redshifts(nz)) then
        write (*,*) ' z out of bounds in MatterPowerAt_Z (',z,')'
        call MPIstop()
    end if

    zlo=min(zi_last,size(PK%redshifts))
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
            matpower(iz) = log(MatterPowerAt_zbin(PK,kh,itf,NL))
            iz=iz+1
        end do
        call spline(zvec(2:4),matpower(2:4),3,SPLINE_DANGLE,SPLINE_DANGLE,ddmat(2:4))
    else if(zhi==nz)then
        zvec = 0
        matpower = 0
        ddmat = 0
        iz = 1
        zvec(1:3)=PK%redshifts(zlo-1:zhi)
        do itf=zlo-1, zhi
            matpower(iz) = log(MatterPowerAt_zbin(PK,kh,itf,NL))
            iz=iz+1
        end do
        call spline(zvec(1:3),matpower(1:3),3,SPLINE_DANGLE,SPLINE_DANGLE,ddmat(1:3))
    else
        iz = 1
        zvec(:)=PK%redshifts(zlo-1:zhi+1)
        do itf=zlo-1, zhi+1
            matpower(iz) = log(MatterPowerAt_zbin(PK,kh,itf,NL))
            iz=iz+1
        end do
        call spline(zvec,matpower,4,SPLINE_DANGLE,SPLINE_DANGLE,ddmat)
    end if

    ho=zvec(3)-zvec(2)
    a0=(zvec(3)-z)/ho
    b0=(z-zvec(2))/ho

    outpower = a0*matpower(2)+b0*matpower(3)+((a0**3-a0)*ddmat(2) &
    +(b0**3-b0)*ddmat(3))*ho**2/6

    outpower = exp(max(-30._dl,outpower))

    end function MatterPowerAt_Z

    end module powerspec


