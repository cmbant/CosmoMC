    module powerspec
    use settings
    use cmbtypes
    use precision
    use ModelData
    use Transfer
    implicit none

    contains

    subroutine Theory_GetMatterPowerData(MTrans, Theory, in)
    !Get total matter power spectrum in units of (h Mpc^{-1})^3 ready for interpolation.
    !Here there definition is < Delta^2(x) > = 1/(2 pi)^3 int d^3k P_k(k)
    !We are assuming that Cls are generated so any baryonic wiggles are well sampled and that matter power
    !sepctrum is generated to beyond the CMB k_max
    Type(MatterTransferData), intent(in) :: MTrans
    Type(TheoryPredictions) Theory
    Type(MatterPowerData) :: Cosmo_PK
    integer, intent(in) :: in
    integer nz,zix

    Theory%num_k = MTrans%num_q_trans
    nz = num_power_redshifts

    call InitPK(Theory, Theory%num_k,nz)

    do zix=1,nz
        call Transfer_GetMatterPowerData(MTrans,Cosmo_PK,1,CP%Transfer%PK_redshifts_index(nz-zix+1))
        if(zix==1) Theory%log_kh=Cosmo_Pk%log_kh
        Theory%redshifts(zix) = CP%Transfer%PK_redshifts(nz-zix+1)
        Theory%matter_power(:,zix) = Cosmo_PK%matpower(:,1)
        Theory%ddmatter_power(:,zix) = Cosmo_PK%ddmat(:,1)
        if(use_nonlinear) then
            call MatterPowerdata_MakeNonlinear(Cosmo_PK)
            Theory%nlmatter_power(:,zix) = Cosmo_PK%matpower(:,1)
            Theory%ddnlmatter_power(:,zix) = Cosmo_PK%ddmat(:,1)
        end if
        call MatterPowerdata_Free(Cosmo_PK)
    end do

    end subroutine Theory_GetMatterPowerData

    function MatterPowerAt_zbin(PK, kh, itf, NNL) result(outpower)
    !Get matter power spectrum at particular k/h by interpolation
    Type(TheoryPredictions) :: PK
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

    if(present(NNL))then
        NL = NNL
    else
        NL = .false.
    end if

    logk = log(kh)
    if (logk < PK%log_kh(1)) then
        if( NL ) then
            matpower=PK%nlmatter_power(1:2,itf)
        else
            matpower=PK%nlmatter_power(1:2,itf)
        end if
        dp = (matpower(2)-matpower(1))/(PK%log_kh(2)-PK%log_kh(1))
        outpower = matpower(1) + dp*(logk-PK%log_kh(1))
    else if (logk > PK%log_kh(PK%num_k)) then
        !Do dodgy linear extrapolation on assumption accuracy of result won't matter
        if( NL ) then
            matpower=PK%nlmatter_power(PK%num_k-1:PK%num_k,itf)
        else
            matpower=PK%nlmatter_power(PK%num_k-1:PK%num_k,itf)
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
    Type(TheoryPredictions) :: PK
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
    Type(TheoryPredictions) :: PK
    real (mcp), intent(in) :: kh, z
    logical, optional, intent(in) :: NNL
    logical :: NL
    integer zlo, zhi, iz, itf
    real(mcp) outpower
    real(mcp) ho,a0,b0
    real(mcp), dimension(4) :: matpower, ddmat, zvec
    integer, save :: zi_last = 1

    if(present(NNL))then
        NL = NNL
    else
        NL = .false.
    end if

    if(z>PK%redshifts(size(PK%redshifts))) then
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
        iz = 2
        zvec(2:4)=PK%redshifts(zlo:zhi+1)
        do itf=zlo, zhi+1
            matpower(iz) = log(MatterPowerAt_zbin(PK,kh,itf,NL))
            iz=iz+1
        end do
        call spline_double(zvec(2:4),matpower(2:4),4,ddmat(2:4))
    else
        iz = 1
        zvec(:)=PK%redshifts(zlo-1:zhi+1)
        do itf=zlo-1, zhi+1
            matpower(iz) = log(MatterPowerAt_zbin(PK,kh,itf,NL))
            iz=iz+1
        end do
        call spline_double(zvec,matpower,4,ddmat)
    end if

    ho=zvec(3)-zvec(2)
    a0=(zvec(3)-z)/ho
    b0=(z-zvec(2))/ho

    outpower = a0*matpower(2)+b0*matpower(3)+((a0**3-a0)*ddmat(2) &
    +(b0**3-b0)*ddmat(3))*ho**2/6

    outpower = exp(max(-30._dl,outpower))

    end function MatterPowerAt_Z

    end module powerspec


