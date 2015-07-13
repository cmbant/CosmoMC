
    module minimize
    use CalcLike
    use Powell_ConstrainedOptimize
    use MatrixUtils
    use MonteCarlo
    use ParamPointSet
    use propose
    implicit none
    private

    !TSaveLoadStateObject->TCheckpointable->TLikelihoodUser->TMinimizer
    Type, abstract, extends(TLikelihoodUser) :: TMinimizer
        logical :: uses_MPI = .false.
        Type(ParamSet) :: MinParams
        integer, allocatable :: minimize_indices(:), minimize_indices_used(:)
        real(mcp), allocatable :: used_param_scales(:), inv_cov(:,:)
        integer, allocatable  :: rotatable_params_used(:)
        integer :: num_rotatable = 0
        integer :: num_rot = 0
        real(mcp), allocatable :: origin(:)
        real(mcp), allocatable :: paramRot(:,:)
        integer, allocatable :: min_params_rot(:)
    contains
    procedure :: ffn
    procedure(FindBestFit), deferred :: FindBestFit
    procedure :: VectToParams
    end type

    abstract interface
    function FindBestFit(this,Params,is_best_bestfit) result(best_like)
    import TMinimizer,ParamSet, mcp
    class(TMinimizer) :: this
    Type(ParamSet) Params
    logical, intent(out) :: is_best_bestfit
    real(mcp) best_like
    end function FindBestFit
    end interface

    Type, extends(TBOBYQA) ::TBOBYQAMin
        class(TMinimizer), pointer :: Minimizer
    contains
    procedure :: funkk
    end type

    Type, extends(TMinimizer) :: TPowellMinimizer
        Type(TBOBYQAMin) :: BOBYQA
        real(mcp) :: max_like_radius = 0.01_mcp
        integer :: max_like_iterations  = 10000
        integer :: minimization_points_factor = 2
        real(mcp) :: minimize_loglike_tolerance = 1e-6_mcp
        logical :: minimize_separate_fast = .true.
        integer :: minimize_mcmc_refine_num = 20
        !MCMC steps per parameter to refine after provisional best fit
        real(mcp) :: minimize_refine_temp = 0.01_mcp
        real(mcp) :: minimize_temp_scale_factor = 5._mcp
        real(mcp)  :: start_trust_radius = 3._mcp
        logical :: minimize_random_start_pos = .false.

    contains
    procedure :: ReadParams => TPowellMinimizer_ReadParams
    procedure :: FindBestFit => TPowellMinimizer_FindBestFit
    procedure :: FindBestFit_indices => TPowellMinimizer_FindBestFit_indices
    end type


    public TMinimizer, TPowellMinimizer

    contains

    function funkk(this, n, X)result(like)
    class(TBOBYQAMin) this
    integer, intent(in) :: n
    real(Powell_CO_prec) :: like
    real(Powell_CO_prec) X(n)

    like = this%Minimizer%ffn(n,X)

    end function funkk

    subroutine VectToParams(this,vect, P)
    class(TMinimizer) :: this
    real(Powell_CO_prec) vect(:)
    Type(ParamSet) P
    integer i,ii

    associate (minimize_indices => this%minimize_indices)
        if (this%num_rot>0) then
            P%P(minimize_indices(this%min_params_rot)) = matmul(this%paramRot,vect(this%min_params_rot)) &
                + this%origin(this%minimize_indices(this%min_params_rot))
        end if

        do i=1, size(minimize_indices)
            if (.not. any(this%min_params_rot==i)) then
                ii=minimize_indices(i)
                P%P(ii)=vect(i)*this%used_param_scales(this%minimize_indices_used(i)) + this%origin(ii)
            end if
        end do

        do i=1, size(minimize_indices)
            if (P%P(minimize_indices(i)) < BaseParams%PMin(minimize_indices(i))-1d-13) &
                call DoAbort(numcat('Minimize:Rotated parameter less than prior boundary:',i))
            if (P%P(minimize_indices(i)) > BaseParams%PMax(minimize_indices(i))+1d-13) &
                call DoAbort(numcat('Minimize:Rotated parameter greater than prior boundary:',i))
            !Just fix up small rounding error
        end do
        P%P(minimize_indices) = max(BaseParams%PMin(minimize_indices),P%P(minimize_indices))
        P%P(minimize_indices) = min(BaseParams%PMax(minimize_indices),P%P(minimize_indices))
    end associate

    end subroutine VectToparams


    function ffn(this, n,vect) result(like)
    class(TMinimizer) :: this
    integer, intent(in) :: n
    real(Powell_CO_prec) :: like
    real(Powell_CO_prec) vect(n)
    Type(ParamSet) P

    P = this%MinParams
    call this%VectToParams(vect, P)

    like = this%LikeCalculator%GetLogLike(P)
    if (like == logZero) then
        print *, 'Warning: Minimizer does not currently properly support non-boundary LogZero rejections'
!        call MpiStop('Minimizer does not currently support non-boundary LogZero rejections')
        call P%Clear(keep=this%MinParams)
    else
        call this%MinParams%Clear(keep=P)
        this%MinParams = P !want to keep e.g. the Age calculation
    end if

    end function ffn

    subroutine TPowellMinimizer_ReadParams(this,Ini)
    class(TPowellMinimizer) :: this
    class(TSettingIni) :: Ini
    character(LEN=:), pointer :: diag_params
    integer num, rotparams(num_params), i, rot_params_used(num_params)

    !radius in normalized parameter space to converge
    call Ini%Read('max_like_radius',this%max_like_radius)
    call Ini%Read('max_like_iterations',this%max_like_iterations)
    !set points factor above 2 to use a denser sampling of space (may be more robust)
    call Ini%Read('minimization_points_factor',this%minimization_points_factor)
    !will exit if function difference between iterations less than minimize_loglike_tolerance (even if radius criterion not met)
    call Ini%Read('minimize_loglike_tolerance',this%minimize_loglike_tolerance, min=0.d0)
    call Ini%Read('minimize_separate_fast',this%minimize_separate_fast)
    call Ini%Read('minimize_mcmc_refine_num',this%minimize_mcmc_refine_num)
    if (this%minimize_mcmc_refine_num>0) then
        call Ini%Read('minimize_refine_temp',this%minimize_refine_temp)
        call Ini%Read('minimize_temp_scale_factor',this%minimize_temp_scale_factor)
    end if
    call Ini%Read('minimize_random_start_pos',this%minimize_random_start_pos)
    this%uses_MPI = this%minimize_random_start_pos

    allocate(this%used_param_scales(num_params_used))
    allocate(this%inv_cov(num_params_used, num_params_used))
    this%inv_cov = BaseParams%covariance_estimate
    call Matrix_Inverse(this%inv_cov)
    do i=1,num_params_used
        this%used_param_scales(i) = min( sqrt( 1/this%inv_cov(i,i)) , &
            (BaseParams%PMax(params_used(i))- BaseParams%PMin(params_used(i)))/this%start_trust_radius/3.)
    end do

    diag_params => Ini%Read_String('minimize_diag_params')
    !Like the proposal for MCMC, but potentially a subset of parameters since can only diagonalize
    !parameters that do not have hard boundaries in the high likelihood region
    if (diag_params/='') then
        num= -1
        call BaseParams%NameMapping%ReadIndices(diag_params, rotparams, num)
    else
        num = 0
    end if

    this%num_rotatable=0
    do i=1, num_params_used
        if (any(rotparams(1:num)==params_used(i))) then
            this%num_rotatable=this%num_rotatable+1
            rot_params_used(this%num_rotatable) = i
        end if
    end do
    allocate(this%rotatable_params_used(this%num_rotatable))
    if (this%num_rotatable>0) then
        this%rotatable_params_used = rot_params_used(1:this%num_rotatable)
    end if

    end subroutine TPowellMinimizer_ReadParams


    function TPowellMinimizer_FindBestFit_indices(this,num_indices,vect) result(best_like)
    class(TPowellMinimizer), target :: this
    integer, intent(in) ::  num_indices
    real(mcp) best_like
    real(Powell_CO_prec) :: vect(num_indices), XL(num_indices), XU(num_indices)
    real(Powell_CO_prec) rhobeg, rhoend
    integer npt,i
    integer rot_params_used(num_indices)

    this%BOBYQA%Minimizer => this
    this%num_rot=0
    do i=1, num_indices
        if (any(this%rotatable_params_used(1:this%num_rotatable)==this%minimize_indices_used(i))) then
            !Rotatable parameters are assumed not to have any prior boundaries near the optimization region
            this%num_rot=this%num_rot+1
            rot_params_used(this%num_rot) = i
            XL(i) = -1d10
            XU(i) = 1d10
        else
            XL(i) = (BaseParams%PMin(this%minimize_indices(i))-this%origin(this%minimize_indices(i))) &
                / this%used_param_scales(this%minimize_indices_used(i))
            XU(i) = (BaseParams%PMax(this%minimize_indices(i))-this%origin(this%minimize_indices(i))) &
                / this%used_param_scales(this%minimize_indices_used(i))
        end if
    end do

    if (allocated(this%min_params_rot)) deallocate(this%min_params_rot)
    if (allocated(this%paramRot)) deallocate(this%paramRot)
    allocate(this%min_params_rot(this%num_rot))

    if (this%num_rot > 0) then
        this%min_params_rot = rot_params_used(1:this%num_rot)
        allocate(this%paramRot(this%num_rot,this%num_rot))
        this%paramRot = this%inv_cov(this%minimize_indices_used(this%min_params_rot),&
            & this%minimize_indices_used(this%min_params_rot))
        call Matrix_Inverse(this%paramRot)
        associate(scales => this%used_param_scales(this%minimize_indices_used(this%min_params_rot)))
            do i=1,this%num_rot
                this%paramRot(:,i) = this%paramRot(:,i) / scales(i)
                this%paramRot(i,:) = this%paramRot(i,:) / scales(i)
            end do
            call Matrix_Cholesky(this%paramRot,zeroed=.true.)
            do i=1,this%num_rot
                this%paramRot(i,:) = this%paramRot(i,:) * scales(i)
            end do
        end associate
    end if

    !Initial and final radius of region required (in normalized units)
    rhobeg =  this%start_trust_radius   !min(nint(minval(XU-XL)/3),1)
    rhoend = this%max_like_radius
    if (this%minimization_points_factor>2) then
        npt = min(this%minimization_points_factor*num_indices,((num_indices +1)*(num_indices +2))/2)
    else
        npt = 2*num_indices +1 !have had some problems using just this
    end if
    This%Bobyqa%FVAL_Converge_difference = this%minimize_loglike_tolerance
    if (Feedback>0) then
        print*,'minimizing ',num_indices, ' with rhobeg, rhoend = ',real(rhobeg), real(rhoend)
        print*,'minimization points: ',npt, 'chi2 converge tol:',real(this%Bobyqa%FVAL_Converge_difference)
    end if
    if (.not. this%Bobyqa%BOBYQA (num_indices ,npt, vect,XL,XU,rhobeg,rhoend,FeedBack+1, this%max_like_iterations)) then
        !BOBYQA generally uses too few operations to get a Hessian estimate
        !I have asked M Powell, and he indeed recommended calculating the Hessian sepratately afterwards
        best_like = logZero
    else
        best_like = this%ffn(num_indices,vect)
    end if
    end function TPowellMinimizer_FindBestFit_indices


    function TPowellMinimizer_FindBestFit(this,Params,is_best_bestfit) result(best_like)
    class(TPowellMinimizer) :: this
    Type(ParamSet) Params, MCParams
    logical, intent(out) :: is_best_bestfit
    real(mcp) best_like, last_like
    real(Powell_CO_prec) :: vect(num_params_used), vect_fast(BaseParams%num_fast)
    real(mcp) :: temperature, scale, last_best, checklike
    integer i
    Type(BlockedProposer) :: Proposer
    Class(TLikeCalculator), allocatable :: LikeCalcMCMC
    class(TMetropolisSampler), allocatable :: MCMC
    real(mcp) StartLike
#ifdef MPI
    real(mcp), allocatable :: bestfit_loglikes(:)
    integer ierror
#endif

    vect = 0 !start at zero by definition (from input center values)
    if (this%minimize_random_start_pos) then
        call BaseParams%SetStartPositions(Params)
    else
        Params%P(1:num_params) = BaseParams%center(1:num_params)
    end if
    if (.not. allocated(this%origin)) allocate(this%Origin(num_params))
    this%origin = Params%P(1:num_params)
    this%MinParams = Params
    ! scale the params so they are all roughly the same order of magnitude
    !normalized parameter bounds
    if (BaseParams%num_fast>0 .and. BaseParams%num_slow /=0 .and. this%minimize_separate_fast) then
        if (Feedback>0) print*,'minmizing fast parameters'
        vect_fast=0
        call Proposer%Init(BaseParams%param_blocks, slow_block_max= slow_tp_max, oversample_fast=1)
        allocate(this%minimize_indices_used(BaseParams%num_fast), source=Proposer%indices(Proposer%Slow%n+1:Proposer%All%n))
        allocate(this%minimize_indices(BaseParams%num_fast))
        this%minimize_indices=params_used(this%minimize_indices_used)
        best_like = this%FindBestFit_indices(BaseParams%num_fast,vect_fast)
        if (Feedback>0) print *,'initial fast parameter minimize logLike: ', best_like
        vect(Proposer%indices(Proposer%Slow%n+1:Proposer%All%n)) = vect_fast
        deallocate(this%minimize_indices,this%minimize_indices_used)
    end if

    do
        if (Feedback>0) print*,'minmizing all parameters'
        allocate(this%minimize_indices, source=params_used)
        allocate(this%minimize_indices_used, source=[(I, I=1, num_params_used)])

        best_like = this%FindBestFit_indices(num_params_used,vect)
        deallocate(this%minimize_indices, this%minimize_indices_used)
        last_like = best_like

        if (BaseParams%num_fast>0 .and. BaseParams%num_slow /=0 .and. this%minimize_separate_fast) then
            if (Feedback>0) print*,'minmizing fast parameters again'
            vect_fast=vect(Proposer%indices(Proposer%Slow%n+1:Proposer%All%n))
            allocate(this%minimize_indices(BaseParams%num_fast), this%minimize_indices_used(BaseParams%num_fast))
            this%minimize_indices_used = Proposer%indices(Proposer%Slow%n+1:Proposer%All%n)
            this%minimize_indices = params_used(this%minimize_indices_used)
            best_like = this%FindBestFit_indices(BaseParams%num_fast,vect_fast)
            if (Feedback>0) print *,'fast parameter minimize logLike: ', best_like
            vect(Proposer%indices(Proposer%Slow%n+1:Proposer%All%n)) = vect_fast
            deallocate(this%minimize_indices,this%minimize_indices_used)
        end if
        !Only finish if sanity check passes
        if (abs(last_like - best_like) < this%minimize_loglike_tolerance*2 .or. &
            this%minimize_mcmc_refine_num>0 .and. abs(last_like - best_like) < max(this%minimize_loglike_tolerance,0.5_mcp)) exit
    end do

    call Params%Clear(keep=this%MinParams)
    Params = this%MinParams

    if (this%minimize_mcmc_refine_num>0) then
        if (Feedback > 0) print *,MpiRank,'Refining minimimum using low temp MCMC'
        if (Feedback > 0) print *,MpiRank, 'Current logLike: ', best_like

        scale = 2
        temperature = this%minimize_refine_temp
        allocate(LikeCalcMCMC, source=this%LikeCalculator)
        !BOBYQA can stop because of numerical errors, do some MCMC steps to do last bit of numerically noisy convergence
        do
            MCParams = Params

            allocate(TMetropolisSampler::MCMC)

            if (Feedback > 0) print *,MpiRank,'Minimize MCMC with temp', temperature
            last_best = best_like
            StartLike = best_like/temperature
            LikeCalcMCMC%temperature = temperature

            call MCMC%InitWithPropose(LikecalcMCMC,null(), propose_scale=scale*sqrt(temperature))
            call MCMC%SetCovariance(BaseParams%covariance_estimate)
            call MCMC%SampleFrom(MCParams, StartLike, this%minimize_mcmc_refine_num * num_params_used)

            if (Feedback > 0) then
                if (MCMC%MaxLike/=logZero) then
                    print *,MpiRank, 'MCMC MaxLike = ', MCMC%MaxLike*temperature
                else
                    print *,MpiRank, 'MCMC chain not moved'
                end if
            end if

            if (MCMC%MaxLike/=logZero .and. MCMC%MaxLike*temperature < best_like) then
                Params%P  = MCMC%MaxLikeParams
                checkLike=LikeCalcMCMC%GetLogLike(Params)*temperature
                !same as MCMC%MaxLike*temperature but want to get everything computed
                if (Feedback>0) print *,MpiRank, 'check likes, best_like:', &
                    & real([checklike, MCMC%MaxLike*temperature, best_like]) !this
                best_like = MCMC%MaxLike*temperature
            end if
            call MCParams%Clear(keep=Params)
            deallocate(MCMC)
            if (last_best - best_like < this%minimize_loglike_tolerance &
                .and. temperature < 4*this%minimize_loglike_tolerance/num_params_used) exit
            if (last_best - best_like < 2*sqrt(num_params_used*temperature)) &
                Temperature = Temperature/ this%minimize_temp_scale_factor
        end do
        deallocate(LikeCalcMCMC)
    end if

    is_best_bestfit=.true.
#ifdef MPI
    if (this%uses_MPI) then
        allocate(bestfit_loglikes(MPIchains))
        call MPI_Allgather(best_like, 1, MPI_real_mcp, &
            bestfit_loglikes, 1,  MPI_real_mcp, MPI_COMM_WORLD, ierror)
        if (MpiRank==0 .and. MPIChains>1) then
            print *,'synched bestfits:', bestfit_loglikes
            if (maxval(bestfit_loglikes)-minval(bestfit_loglikes) >1) then
                print *,'WARNING: big spread in log-likes'
            elseif (maxval(bestfit_loglikes)-minval(bestfit_loglikes) >0.2) then
                print *,'WARNING: modest spread in log-likes'
            end if
        end if
        is_best_bestfit = minval(bestfit_loglikes)==best_like
        deallocate(bestfit_loglikes)
        do i=0, MpiChains-1
            call MPI_Barrier(MPI_COMM_WORLD,ierror)
            if (i==MpiRank .and. .not. is_best_bestfit .and. Feedback>0) then
                print *,MpiRank, 'best_like params:'
                call this%LikeCalculator%WriteParamsHumanText(stdout,Params, best_like)
            end if
        end do
        call MPI_Barrier(MPI_COMM_WORLD,ierror)
    end if
#endif

    end function TPowellMinimizer_FindBestFit


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
    !       M(:,i) = M(:,i) * BaseParams%PWidth(params_used(i))
    !       M(i,:) = M(i,:) * BaseParams%PWidth(params_used(i))
    !  end do
    !
    !end function BestFitCovmatEstimate


    end module minimize
