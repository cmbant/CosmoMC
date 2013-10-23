
    module minimize
    use ParamDef
    use CalcLike
    use Powell_ConstrainedOptimize
    use MatrixUtils
    use MonteCarlo
    implicit none
    private

    Type(ParamSet), save :: MinParams

    real(mcp) :: max_like_radius = 0.01_mcp
    integer :: max_like_iterations  = 10000
    integer :: minimization_points_factor = 2
    real(mcp) :: minimize_loglike_tolerance = 0._mcp
    logical :: minimize_separate_fast = .true.
    integer :: minimize_mcmc_refine_num = 20
    !MCMC steps per parameter to refine after provisional best fit
    real(mcp) :: minimize_refine_temp = 0.01_mcp
    real(mcp) :: minimize_temp_scale_factor = 5._mcp
    real(mcp), parameter :: start_trust_radius = 3._mcp
    logical :: minimize_random_start_pos = .false.

    integer, allocatable  :: rotatable_params_used(:)
    integer :: num_rotatable = 0

    integer :: num_rot = 0
    integer, allocatable :: minimize_indices(:), minimize_indices_used(:)
    real(mcp), allocatable :: used_param_scales(:), inv_cov(:,:)

    real(mcp), allocatable :: paramRot(:,:)
    integer, allocatable :: min_params_rot(:)

    logical :: minimize_uses_MPI = .false.

    public FindBestFit, WriteBestFitParams,Minimize_ReadIni,minimize_uses_MPI

    contains

    subroutine VectToParams(vect, P)
    real(Powell_CO_prec) vect(:)
    Type(ParamSet) P
    integer i,ii

    if (num_rot>0) then
        P%P(minimize_indices(min_params_rot)) = matmul(paramRot,vect(min_params_rot)) &
        + Scales%center(minimize_indices(min_params_rot))
    end if

    do i=1, size(minimize_indices)
        if (.not. any(min_params_rot==i)) then
            ii=minimize_indices(i)
            P%P(ii)=vect(i)*used_param_scales(minimize_indices_used(i)) + Scales%center(ii)
        end if
    end do

    do i=1, size(minimize_indices)
        if (P%P(minimize_indices(i)) < Scales%PMin(minimize_indices(i))-1d-13) &
        call DoAbort(numcat('Minimize:Rotated parameter less than prior boundary:',i))
        if (P%P(minimize_indices(i)) > Scales%PMax(minimize_indices(i))+1d-13) &
        call DoAbort(numcat('Minimize:Rotated parameter greater than prior boundary:',i))
        !Just fix up small rounding error
    end do
    P%P(minimize_indices) = max(Scales%PMin(minimize_indices),P%P(minimize_indices))
    P%P(minimize_indices) = min(Scales%PMax(minimize_indices),P%P(minimize_indices))

    end subroutine VectToparams


    function ffn(n,vect) result(like)
    implicit none
    integer, intent(in) :: n
    real(Powell_CO_prec) :: like
    real(Powell_CO_prec) vect(n)
    Type(ParamSet) P

    P = MinParams
    call VectToParams(vect, P)

    like = GetLogLike(P)

    call AcceptReject(.true.,MinParams%Info, P%Info)
    MinParams = P !want to keep e.g. the Age calculation

    end function ffn

    subroutine Minimize_ReadIni(Ini)
    use IniFile
    Type(TIniFile), intent(in) :: Ini
    character(LEN=Ini_max_string_len) :: diag_params
    integer num, rotparams(num_params), i, rot_params_used(num_params)

    !radius in normalized parameter space to converge
    max_like_radius = Ini_Read_Double_File(Ini,'max_like_radius',max_like_radius)
    max_like_iterations = Ini_Read_Int_File(Ini,'max_like_iterations',max_like_iterations)
    !set points factor above 2 to use a denser sampling of space (may be more robust)
    minimization_points_factor = Ini_Read_Int_File(Ini,'minimization_points_factor',minimization_points_factor)
    !will exit if function difference between iterations less than minimize_loglike_tolerance (even if radius criterion not met)
    minimize_loglike_tolerance = Ini_Read_double_File(Ini,'minimize_loglike_tolerance',minimize_loglike_tolerance)
    minimize_separate_fast = Ini_Read_Logical_File(Ini,'minimize_separate_fast',minimize_separate_fast)
    minimize_mcmc_refine_num = Ini_Read_Int_File(Ini,'minimize_mcmc_refine_num',minimize_mcmc_refine_num)
    if (minimize_mcmc_refine_num>0) then
        minimize_refine_temp = Ini_read_Double_File(Ini,'minimize_refine_temp',minimize_refine_temp)
        minimize_temp_scale_factor = Ini_Read_Double_File(ini,'minimize_temp_scale_factor',minimize_temp_scale_factor)
    end if
    minimize_random_start_pos = Ini_Read_Logical_File(Ini,'minimize_random_start_pos',minimize_random_start_pos)
    minimize_uses_MPI = minimize_random_start_pos

    allocate(used_param_scales(num_params_used))
    allocate(inv_cov(num_params_used, num_params_used))
    inv_cov = Proposer%propose_matrix
    call Matrix_Inverse(inv_cov)
    do i=1,num_params_used
        used_param_scales(i) = min( sqrt( 1/inv_cov(i,i)) , &
        (Scales%PMax(params_used(i))- Scales%PMin(params_used(i)))/start_trust_radius/3.)
    end do

    diag_params = Ini_Read_String_File(Ini,'minimize_diag_params')
    !Like the proposal for MCMC, but potentially a subset of parameters since can only diagonalize
    !parameters that do not have hard boundaries in the high likelihood region
    if (diag_params/='') then
        num= -1
        call ParamNames_ReadIndices(NameMapping,diag_params, rotparams, num)
    else 
        num = 0
    end if

    num_rotatable=0
    do i=1, num_params_used
        if (any(rotparams(1:num)==params_used(i))) then
            num_rotatable=num_rotatable+1
            rot_params_used(num_rotatable) = i
        end if
    end do
    allocate(rotatable_params_used(num_rotatable))
    if (num_rotatable>0) then
        rotatable_params_used = rot_params_used(1:num_rotatable)
    end if

    end subroutine Minimize_ReadIni


    function FindBestFit_indices(num_indices,vect) result(best_like)
    use ParamDef
    integer, intent(in) ::  num_indices
    real(mcp) best_like
    real(Powell_CO_prec) :: vect(num_indices), XL(num_indices), XU(num_indices)
    real(Powell_CO_prec) rhobeg, rhoend
    integer npt,i
    integer rot_params_used(num_indices)

    num_rot=0
    do i=1, num_indices
        if (any(rotatable_params_used(1:num_rotatable)==minimize_indices_used(i))) then
            !Rotatable parameters are assumed not to have any prior boundaries near the optimization region
            num_rot=num_rot+1
            rot_params_used(num_rot) = i
            XL(i) = -1d10
            XU(i) = 1d10
        else
            XL(i) = (Scales%PMin(minimize_indices(i))-Scales%center(minimize_indices(i))) &
            / used_param_scales(minimize_indices_used(i))
            XU(i) = (Scales%PMax(minimize_indices(i))-Scales%center(minimize_indices(i))) &
            / used_param_scales(minimize_indices_used(i))
        end if
    end do

    if (allocated(min_params_rot)) deallocate(min_params_rot)
    if (allocated(paramRot)) deallocate(paramRot)
    allocate(min_params_rot(num_rot))

    if (num_rot > 0) then
        min_params_rot = rot_params_used(1:num_rot)
        allocate(paramRot(num_rot,num_rot))
        paramRot = inv_cov(minimize_indices_used(min_params_rot),minimize_indices_used(min_params_rot))
        call Matrix_Inverse(paramRot)
        do i=1,num_rot
            paramRot(:,i) = paramRot(:,i) / used_param_scales(minimize_indices_used(min_params_rot(i)))
            paramRot(i,:) = paramRot(i,:) / used_param_scales(minimize_indices_used(min_params_rot(i)))
        end do
        call Matrix_Cholesky(paramRot,zeroed=.true.)
        do i=1,num_rot
            paramRot(i,:) = paramRot(i,:) * used_param_scales(minimize_indices_used(min_params_rot(i)))
        end do
    end if

    !Initial and final radius of region required (in normalized units)
    rhobeg =  start_trust_radius   !min(nint(minval(XU-XL)/3),1)
    rhoend = max_like_radius
    if (minimization_points_factor>2) then
        npt = min(minimization_points_factor*num_indices,((num_indices +1)*(num_indices +2))/2)
    else
        npt = 2*num_indices +1 !have had some problems using just this
    end if
    FVAL_Converge_difference = minimize_loglike_tolerance
    if (Feedback>0) then
        print*,'minimizing ',num_indices, ' with rhobeg, rhoend = ',real(rhobeg), real(rhoend)
        print*,'minimization points: ',npt, 'chi2 converge tol:',real(FVAL_Converge_difference)
    end if
    if (.not. BOBYQA (ffn, num_indices ,npt, vect,XL,XU,rhobeg,rhoend,FeedBack+1, max_like_iterations)) then
        !BOBYQA generally uses too few operations to get a Hessian estimate
        !I have asked M Powell, and he indeed recommended calculating the Hessian sepratately afterwards
        best_like = logZero
    else
        best_like = ffn(num_indices,vect)
    end if
    end function FindBestFit_indices


    function FindBestFit(Params,is_best_bestfit) result(best_like)
    use ParamDef
    Type(ParamSet) Params
    logical, intent(out) :: is_best_bestfit
    real(mcp) best_like, last_like
    real(Powell_CO_prec) :: vect(num_params_used), vect_fast(num_fast)
    real(mcp) :: temp, scale, slike, last_best, checklike
    real(mcp), allocatable :: bestfit_loglikes(:)
    integer ierror,i

    vect = 0 !start at zero by definition (from input center values)
    if (minimize_random_start_pos) then
        call SetStartPositions(Params)
    else
        Params%P(1:num_params) = Scales%center(1:num_params)
    end if
    MinParams = Params
    ! scale the params so they are all roughly the same order of magnitude
    !normalized parameter bounds
    if (num_fast>0 .and. num_slow /=0 .and. minimize_separate_fast) then
        if (Feedback>0) print*,'minmizing fast parameters'
        vect_fast=0
        allocate(minimize_indices(num_fast))
        allocate(minimize_indices_used(num_fast))
        minimize_indices_used = Proposer%indices(Proposer%Slow%n+1:Proposer%All%n)
        minimize_indices = params_used(minimize_indices_used)
        best_like = FindBestFit_indices(num_fast,vect_fast)
        if (Feedback>0) print *,'initial fast parameter minimize logLike: ', best_like
        vect(Proposer%indices(Proposer%Slow%n+1:Proposer%All%n)) = vect_fast
        deallocate(minimize_indices,minimize_indices_used)
    end if

    do
        if (Feedback>0) print*,'minmizing all parameters'
        allocate(minimize_indices(num_params_used))
        allocate(minimize_indices_used(num_params_used))
        minimize_indices_used = [1:num_params_used]
        minimize_indices = params_used
        best_like = FindBestFit_indices(num_params_used,vect)
        deallocate(minimize_indices, minimize_indices_used)
        last_like = best_like

        if (num_fast>0 .and. num_slow /=0 .and. minimize_separate_fast) then
            if (Feedback>0) print*,'minmizing fast parameters again'
            vect_fast=vect(Proposer%indices(Proposer%Slow%n+1:Proposer%All%n))
            allocate(minimize_indices(num_fast), minimize_indices_used(num_fast))
            minimize_indices_used = Proposer%indices(Proposer%Slow%n+1:Proposer%All%n)
            minimize_indices = params_used(minimize_indices_used)
            best_like = FindBestFit_indices(num_fast,vect_fast)
            if (Feedback>0) print *,'fast parameter minimize logLike: ', best_like
            vect(Proposer%indices(Proposer%Slow%n+1:Proposer%All%n)) = vect_fast
            deallocate(minimize_indices,minimize_indices_used)
        end if
        !Only finish if sanity check passes
        if (abs(last_like - best_like) < minimize_loglike_tolerance*2 .or. minimize_mcmc_refine_num>0) exit
    end do

    call AcceptReject(.true.,Params%Info, MinParams%Info)
    Params = MinParams

    if (minimize_mcmc_refine_num>0) then
        if (Feedback > 0) print *, 'Refining minimimum using low temp MCMC'
        if (Feedback > 0) print *,'Current logLike: ', best_like
        temp = temperature
        scale = propose_scale
        slike = StartLike
        MCMC_outputs = .false.

        !BOBYQA can stop because of numerical errors, do some MCMC steps to do last bit of numerically noisy convergence
        Temperature = minimize_refine_temp
        do
            if (Feedback > 0) print *,'Minimize MCMC with temp', temperature
            last_best = best_like
            propose_scale = scale*sqrt(temperature)
            StartLike = best_like/temperature
            num_accept=0
            num=0
            call MCMCSample(Params, minimize_mcmc_refine_num * num_params_used)
            if (Feedback > 0) then
                if (MaxLike/=logZero) then
                    print *,MpiRank, 'MCMC MaxLike = ', MaxLike*temperature
                else
                    print *,MpiRank, 'MCMC chain not moved'
                end if
            end if

            if (MaxLike/=logZero .and. MaxLike*temperature < best_like) then
                Params%P  = MaxLikeParams
                checkLike=GetLogLike(Params)*temperature
                print *,MpiRank, 'check like, best_like:', checklike, best_like !avoid recursive IO
                best_like = MaxLike*temperature
            end if
            if (last_best - best_like < minimize_loglike_tolerance &
            .and. temperature < 4*minimize_loglike_tolerance/num_params_used) exit
            if (last_best - best_like < 2*sqrt(num_params_used*temperature)) &
            Temperature = Temperature/ minimize_temp_scale_factor
        end do
        temperature = temp
        propose_scale = scale
        StartLike = slike
        MCMC_outputs = .true.
    end if

    is_best_bestfit=.true.
#ifdef MPI
    if (minimize_uses_MPI) then
        allocate(bestfit_loglikes(MPIchains))
        call MPI_Allgather(best_like, 1, MPI_real_mcp, &
        bestfit_loglikes, 1,  MPI_real_mcp, MPI_COMM_WORLD, ierror)
        if (MpiRank==0 .and. MPIChains>1) then
            print *,'synched bestfits:', bestfit_loglikes
            if (maxval(bestfit_loglikes)-minval(bestfit_loglikes) >1) then
                print *,'WARNING: big spread in log-likes'
            elseif (maxval(bestfit_loglikes)-minval(bestfit_loglikes) >0.2) then
                print *,'WARNING: modest big spread in log-likes'
            end if
        end if
        is_best_bestfit = minval(bestfit_loglikes)==best_like
        deallocate(bestfit_loglikes)
        do i=0, MpiChains-1
            call MPI_Barrier(MPI_COMM_WORLD,ierror)
            if (i==MpiRank .and. .not. is_best_bestfit .and. Feedback>0) then
                print *,MpiRank, 'best_like params:'
                call WriteParamsHumanText(stdout,Params, best_like)
            end if
        end do
        call MPI_Barrier(MPI_COMM_WORLD,ierror)
    end if
#endif

    end function FindBestFit


    subroutine WriteParamsHumanText(aunit, P, like)
    use settings
    use cmbtypes
    use ParamDef
    use DataLikelihoodList
    implicit none
    Type(ParamSet) P
    real(mcp), intent(in), optional :: like
    integer, intent(in) :: aunit
    Type(mc_real_pointer) :: derived
    integer numderived
    integer isused,i

    if (present(like)) then
        write (aunit,*) '-log(Like) = ',like
        write (aunit,*) ' chi-sq    = ',like*2
        write (aunit,*) ''
    end if

    do isused = 0,1
        do i=1, num_params
            if (isused==0 .and. Scales%PWidth(i)/=0 .or. isused==1 .and. Scales%PWidth(i)==0) then
                write(aunit,'(1I5,1E15.7,"   ",1A22)', advance='NO') &
                i, P%P(i), ParamNames_name(NameMapping,i)
                write (aunit,'(a)') trim(NameMapping%label(i))
            end if
        end do
        write (aunit,*) ''
    end do

    if (generic_mcmc) return

    numderived = Parameterization%CalcDerivedParams(P%P,P%Theory, derived)
    do i=1, numderived
        write(aunit,'(1I5,1E15.7,"   ",1A22)', advance='NO') &
        num_params+i, derived%P(i), ParamNames_name(NameMapping,num_params + i )
        write (aunit,'(a)') trim(NameMapping%label(num_params+i))
    end do
    deallocate(derived%P)

    if (present(like)) then
        write(aunit,*) ''
        write(aunit,*) '-log(Like)     chi-sq   data'
        call DataLikelihoods%WriteLikelihoodContribs(aunit, P%likelihoods)
    end if

    end  subroutine WriteParamsHumanText

    subroutine WriteBestFitParams(like, Params, fname)
    real(mcp) like
    Type(ParamSet) Params
    character(LEN=*), intent(in) :: fname

    call CreateTxtFile(fname,tmp_file_unit)
    call WriteParamsHumanText(tmp_file_unit,Params, like)
    close(tmp_file_unit)
    if (Feedback>0) call WriteParamsHumanText(stdout,Params, like)

    end subroutine WriteBestFitParams


    !function BestFitCovmatEstimate(tol) result (M)
    !!Cannot find any way to get usable Hessian matrix from BOBYQA, even with high NPT
    !!This is not used..
    ! use MatrixUtils
    ! real(mcp) M(num_params_used,num_params_used), diag(num_params_used)
    ! real(mcp) tol !tol is the smallest normalized eigenvalue to allow (e..g fraction of input width)
    ! integer i
    !
    !  M = BOBYQA_Hessian
    !  !this may not be invertible, let's make sure all eigenvalues are positive
    !
    !  call  Matrix_Diagonalize(M, diag, num_params_used)
    !       !Does m = U diag U^T, returning U in M
    !  if (Feedback > 0) print *,  'Un-regularized Hessian evalues: ', diag
    !  do i=1,num_params_used
    !    diag(i) = max(diag(i), tol)
    !    M(:,i) = M(:,i) / sqrt(diag(i))
    !  end do
    !  M = matmul(M, transpose(M))
    !
    !  !now go back to un-normalized units
    !  do i=1, num_params_used
    !       M(:,i) = M(:,i) * Scales%PWidth(params_used(i))
    !       M(i,:) = M(i,:) * Scales%PWidth(params_used(i))
    !  end do
    !
    !end function BestFitCovmatEstimate


    end module minimize
