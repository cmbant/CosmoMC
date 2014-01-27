    module CalcLike
    use DataLikelihoodList
    implicit none

    real(mcp) :: Temperature  = 1

    Type BaseLikeCalculator
    contains
    procedure :: AddLike
    procedure :: GetLogLikeBounds
    procedure :: GetLogPriors
    end type BaseLikeCalculator

    type, extends(BaseLikeCalculator) :: GenericLikeCalculator
    contains
    procedure :: GetLogLike => Generic_GetLogLike
    procedure, nopass :: GenericLikelihoodFunction
    end type GenericLikeCalculator

    type, extends(BaseLikeCalculator) :: TheoryLikeCalculator
        class(TTheoryParams), allocatable :: TheoryParams
        class(TCalculationAtParamPoint), pointer :: Params
        logical changeMask(max_num_params)
        logical SlowChanged, SemiSlowchanged
    contains
    procedure :: GetLogLikePost => TheoryLike_GetLogLikePost
    procedure :: SetTheoryParams => TheoryLike_SetTheoryParams
    procedure :: CheckProirCuts => TheoryLike_CheckPriorCuts
    procedure :: CalculateRequiredTheoryChanges =>TheoryLike_CalculateRequiredTheoryChanges
    procedure :: GetTheoryForLike=>TheoryLike_GetTheoryForLike
    procedure :: GetLogLikeWithTheorySet => TheoryLike_LogLikeWithTheorySet
    procedure :: TestLikelihoodFunction
    end type TheoryLikeCalculator

    contains

    subroutine AddLike(this, CurrentLike, LikeToAdd)
    class(BaseLikeCalculator) :: this
    real(mcp), intent(in) :: LikeToAdd
    real(mcp) CurrentLike
    if (CurrentLike/=LogZero) then
        if (LikeToAdd == logZero) then
            CurrentLike = LogZero
        else
            CurrentLike = CurrentLike + LikeToAdd/Temperature
        end if
    end if
    end subroutine AddLike

    function GetLogLikeBounds(this,Params)
    class(BaseLikeCalculator) :: this
    class(TCalculationAtParamPoint):: Params
    real(mcp) :: GetLogLikeBounds

    if (any(Params%P(:num_params) > Scales%PMax(:num_params)) .or. &
    any(Params%P(:num_params) < Scales%PMin(:num_params))) then
        GetLogLikeBounds = logZero
    else
        GetLogLikeBounds=0
    end if

    end function GetLogLikeBounds

    function GetLogPriors(this, P) result(logLike)
    class(BaseLikeCalculator) :: this
    integer i
    real(mcp), intent(in) :: P(num_params)
    real(mcp) logLike

    logLike=0
    do i=1,num_params
        if (Scales%PWidth(i)/=0 .and. GaussPriors%std(i)/=0) then
            logLike = logLike + ((P(i)-GaussPriors%mean(i))/GaussPriors%std(i))**2
        end if
    end do
    logLike=logLike/2

    end function GetLogPriors

    function GenericLikelihoodFunction(Params)
    class(TCalculationAtParamPoint)  Params
    real(mcp) :: GenericLikelihoodFunction

    !Used when you want to plug in your own CMB-independent likelihood function:
    !set generic_mcmc=.true. in settings.f90, then write function here returning -Ln(Likelihood)
    !Parameter array is Params%P, so e.g. 2D unit Gaussian would be
    !GenericLikelihoodFunction = (Params%P(1)**2+Params%P(2)**2)/2
    GenericLikelihoodFunction = LogZero
    call MpiStop('GenericLikelihoodFunction: need to write this function!')

    end function GenericLikelihoodFunction


    function Generic_GetLogLike(this, Params) result(GetLogLike)!Get -Ln(Likelihood) for chains
    class(GenericLikeCalculator) :: this
    class(TCalculationAtParamPoint) :: Params
    real(mcp) GetLogLike

    GetLogLike = this%GetLogLikeBounds(Params)
    if (GetLogLike==LogZero) return
    call this%AddLike(GetLogLike,this%GenericLikelihoodFunction(Params))
    if (GetLogLike==LogZero) return
    call this%AddLike(GetLogLike,this%getLogPriors(Params%P))

    end function Generic_GetLogLike


    function TestLikelihoodFunction(this,Params)
    class(TheoryLikeCalculator) :: this
    class(TCalculationAtParamPoint) Params
    real(mcp) :: TestLikelihoodFunction
    real(mcp), allocatable, save :: covInv(:,:)
    real(mcp) X(num_params_used)

    if (.not. allocated(covInv)) then
        allocate(covInv(num_params_used,num_params_used))
        covInv = test_cov_matrix
        call Matrix_Inverse(covInv)
    end if
    X = Params%P(params_used) - scales%Center(params_used)
    TestLikelihoodFunction= dot_product(X, matmul(covInv, X))/2

    end function TestLikelihoodFunction

    subroutine TheoryLike_SetTheoryParams(this, Params)
    class(TheoryLikeCalculator) :: this
    class(TCalculationAtParamPoint), target :: Params

    if (.not. allocated(this%TheoryParams)) call Parameterization%NewTheoryParams(this%TheoryParams)
    call Parameterization%ParamArrayToTheoryParams(Params%P,this%TheoryParams)
    this%Params => Params

    end subroutine TheoryLike_SetTheoryParams


    function TheoryLike_GetLogLike(this, Params) result(GetLogLike)!Get -Ln(Likelihood) for chains
    class(TheoryLikeCalculator) :: this
    class(TCalculationAtParamPoint) :: Params
    real(mcp) GetLogLike

    GetLogLike = this%GetLogLikeBounds(Params)
    if (GetLogLike==LogZero) return

    if (test_likelihood) then
        call this%AddLike(GetLogLike,this%TestLikelihoodFunction(Params))
        call this%AddLike(GetLogLike,this%GetLogPriors(Params%P))
    else
        call this%SetTheoryParams(Params)
        call this%AddLike(GetLogLike, Parameterization%NonBaseParameterPriors(this%TheoryParams))
        if (GetLogLike == logZero) return
        if (.not. Params%Info%validInfo) then
            this%changeMask(1:num_params) = .true.
        else
            this%changeMask(1:num_params) = Params%Info%lastParamArray(1:num_params)/=Params%P(1:num_params)
        end if

        if (this%CalculateRequiredTheoryChanges()) then
            call this%AddLike(GetLogLike, this%GetLogLikeWithTheorySet())
        else
            GetLogLike = logZero
        end if

        if (GetLogLike/=logZero) Params%Info%lastParamArray(1:num_params) = Params%P(1:num_params)
    end if

    if (Feedback>2 .and. GetLogLike/=LogZero) &
    call DataLikelihoods%WriteLikelihoodContribs(stdout, Params%likelihoods)

    end function TheoryLike_GetLogLike

    function TheoryLike_CheckPriorCuts(this, Params) result(checkPriorCuts)
    class(TheoryLikeCalculator) :: this
    real(mcp)  CheckPriorCuts
    class(TCalculationAtParamPoint) :: Params

    CheckPriorCuts = this%GetLogLikeBounds(Params)
    if (CheckPriorCuts==LogZero) return

    call this%SetTheoryParams(Params)
    CheckPriorCuts = Parameterization%NonBaseParameterPriors(this%TheoryParams)

    end function TheoryLike_CheckPriorCuts


    function TheoryLike_GetLogLikePost(this,Params, do_like) result(GetLogLikePost)
    !for importance sampling where theory may be pre-stored
    class(TheoryLikeCalculator) :: this
    real(mcp)  GetLogLikePost
    class(TCalculationAtParamPoint) :: Params
    logical, optional, intent(in) :: do_like(DataLikelihoods%count)

    GetLogLikePost = this%GetLogLikeBounds(Params)
    if (GetLogLikePost==LogZero) return

    this%SlowChanged = .false.
    this%ChangeMask = .true.

    call this%SetTheoryParams(Params)
    call this%AddLike(GetLogLikePost,Parameterization%NonBaseParameterPriors(this%TheoryParams))
    if (GetLogLikePost == logZero) return

    call this%AddLike(GetLogLikePost, this%GetLogLikeWithTheorySet(do_like))

    end function TheoryLike_GetLogLikePost

    
    function TheoryLike_LogLikeWithTheorySet(this, likelihood_mask) result(logLike)
    class(TheoryLikeCalculator) :: this
    logical, intent(in), optional :: likelihood_mask(DataLikelihoods%count)
    real(mcp) logLike
    real(mcp) itemLike
    class(DataLikelihood), pointer :: like => null()
    integer i
    logical :: do_like(DataLikelihoods%count)

    if (present(likelihood_Mask)) then
        do_like = likelihood_mask
    else
        do_like = .true.
    end if
    logLike = logZero
    call GetTheoryForLike(like) !chance to initalize
    do i= 1, DataLikelihoods%count
        if (do_like(i)) then
            like => DataLikelihoods%Item(i)
            if (any(like%dependent_params(1:num_params) .and. this%changeMask(1:num_params) )) then
                call this%GetTheoryForLike(like)
                itemLike = like%LogLike(this%TheoryParams, this%Params%Theory, this%Params%P(like%nuisance_indices))

                if (itemLike == logZero) return
                this%Params%Likelihoods(i) = itemLike
            end if
        end if
    end do
    logLike = sum(this%Params%likelihoods(1:DataLikelihoods%Count))
    logLike = logLike + this%GetLogPriors(this%Params%P)

    end function TheoryLike_LogLikeWithTheorySet

    
    logical function TheoryLike_CalculateRequiredTheoryChanges(this) 
    class(TheoryLikeCalculator) :: this

    !Set this%Params theory entries (this%Params%Theory) as required for current this%changeMask and this%TheoryParams
    this%Params%Info%validInfo = .true.
    TheoryLike_CalculateRequiredTheoryChanges = .true.
    
    end function TheoryLike_CalculateRequiredTheoryChanges

    subroutine TheoryLike_GetTheoryForLike(this,Like)
    class(TheoryLikeCalculator) :: this
    class(DataLikelihood), pointer :: like

    !If needed, likelihood specific calculation/initalization
    end subroutine TheoryLike_GetTheoryForLike


    end module CalcLike
