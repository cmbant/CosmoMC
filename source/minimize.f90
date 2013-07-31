
    module minimize
    use ParamDef
    use CalcLike
    use Powell_ConstrainedOptimize
    implicit none
    private

    Type(ParamSet), save :: MinParams
    integer :: minimization_points_factor = 2
    real(mcp) :: minimize_loglike_tolerance = 0.1
    integer, allocatable :: minimize_indices(:)
    public FindBestFit, WriteBestFitParams,minimization_points_factor,minimize_loglike_tolerance

    contains

    subroutine VectToParams(vect, P)
    real(Powell_CO_prec) vect(:)
    Type(ParamSet) P
    integer i,ii

    do i=1, size(minimize_indices)
        ii=minimize_indices(i)
        P%P(ii)=vect(i)*Scales%PWidth(ii) + Scales%center(ii)
    end do

    end subroutine VectToparams

    subroutine ParamsToVect(P,vect)
    real(Powell_CO_prec) vect(:)
    Type(ParamSet) P
    integer i,ii

    do i=1, size(minimize_indices)
        ii=minimize_indices(i)
        vect(i)=(P%P(ii)-Scales%center(ii))/Scales%PWidth(ii)
    end do

    end subroutine ParamsToVect


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

    function FindBestFit_indices(num_indices,vect,sigma_frac_err, max_iterations) result(best_like)
    use ParamDef
    integer, intent(in) :: max_iterations,num_indices
    real(mcp), intent(in) :: sigma_frac_err
    real(mcp) best_like
    real(Powell_CO_prec) :: vect(num_indices), XL(num_indices), XU(num_indices)
    real(Powell_CO_prec) rhobeg, rhoend
    integer npt

    XL = (Scales%PMin(minimize_indices)-Scales%center(minimize_indices)) / Scales%PWidth(minimize_indices)
    XU = (Scales%PMax(minimize_indices)-Scales%center(minimize_indices)) / Scales%PWidth(minimize_indices)

    !Initial and final radius of region required (in normalized units)
    rhobeg = min(nint(minval(XU-XL)/3),1)
    rhoend = sigma_frac_err
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
    if (.not. BOBYQA (ffn, num_indices ,npt, vect,XL,XU,rhobeg,rhoend,FeedBack+1, max_iterations)) then
        !BOBYQA generally uses too few operations to get a Hessian estimate
        !I have asked M Powell, and he indeed recommended calculating the Hessian sepratately afterwards
        best_like = logZero
    else
        best_like = ffn(num_indices,vect)
    end if
    end function FindBestFit_indices


    function FindBestFit(Params,sigma_frac_err, max_iterations) result(best_like)
    use ParamDef
    Type(ParamSet) Params
    integer, intent(in) :: max_iterations
    real(mcp), intent(in) :: sigma_frac_err
    real(mcp) best_like, last_like
    real(Powell_CO_prec) :: vect(num_params_used), vect_fast(num_fast)

    vect = 0 !start at zero by definition (from input center values)

    MinParams = Params
    ! scale the params so they are all roughly the same order of magnitude
    !call ParamsToVect(MinParams, vect)
    !normalized parameter bounds
    if (num_fast>0 .and. num_slow /=0) then
        if (Feedback>0) print*,'minmizing fast parameters'
        vect_fast=0
        allocate(minimize_indices(num_fast))
        minimize_indices = params_used(Proposer%indices(Proposer%Slow%n+1:Proposer%All%n))
        best_like = FindBestFit_indices(num_fast,vect_fast,sigma_frac_err, max_iterations)
        if (Feedback>0) print *,'initial fast parameter minimize logLike: ', best_like
        vect(Proposer%indices(Proposer%Slow%n+1:Proposer%All%n)) = vect_fast
        deallocate(minimize_indices)
    end if

    do
        if (Feedback>0) print*,'minmizing all parameters'
        allocate(minimize_indices(num_params_used))
        minimize_indices = params_used
        best_like = FindBestFit_indices(num_params_used,vect,sigma_frac_err, max_iterations)
        deallocate(minimize_indices)
        last_like = best_like

        if (num_fast>0 .and. num_slow /=0) then
            if (Feedback>0) print*,'minmizing fast parameters again'
            vect_fast=vect(Proposer%indices(Proposer%Slow%n+1:Proposer%All%n))
            allocate(minimize_indices(num_fast))
            minimize_indices = params_used(Proposer%indices(Proposer%Slow%n+1:Proposer%All%n))
            best_like = FindBestFit_indices(num_fast,vect_fast,sigma_frac_err, max_iterations)
            if (Feedback>0) print *,'fast parameter minimize logLike: ', best_like
            vect(Proposer%indices(Proposer%Slow%n+1:Proposer%All%n)) = vect_fast
            deallocate(minimize_indices)
        end if
        !Only finish if sanity check passes
        if (abs(last_like - best_like) < minimize_loglike_tolerance*2) exit
    end do

    call AcceptReject(.true.,Params%Info, MinParams%Info)
    Params = MinParams

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
