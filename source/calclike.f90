    module CalcLike
    use DataLikelihoodList
    use GeneralTypes
    use BaseParameters
    implicit none

    real(mcp) :: Temperature  = 1

    Type, extends(TConfigClass) :: TLikeCalculator
        logical :: test_likelihood= .false.
        real(mcp), allocatable :: test_cov_matrix(:,:)
    contains
    procedure :: AddLike
    procedure :: GetLogLikeBounds
    procedure :: GetLogPriors
    procedure :: GetLogLike
    procedure :: GetLogLikeMain
    procedure :: ReadParams => TLikeCalculator_ReadParams
    procedure :: TestLikelihoodFunction => TLikeCalculator_TestLikelihoodFunction
    end type TLikeCalculator

    type, extends(TLikeCalculator) :: GenericLikeCalculator
    contains
    procedure :: GetLogLikeMain => Generic_GetLogLikeMain
    end type GenericLikeCalculator

    type, extends(TLikeCalculator) :: TheoryLikeCalculator
        class(TTheoryParams), allocatable :: TheoryParams
        class(TCalculationAtParamPoint), pointer :: Params
        logical changeMask(max_num_params)
        logical SlowChanged, SemiSlowchanged
    contains
    procedure :: GetLogLikeMain => TheoryLike_GetLogLikeMain
    procedure :: GetLogLikePost => TheoryLike_GetLogLikePost
    procedure :: SetTheoryParams => TheoryLike_SetTheoryParams
    procedure :: CheckProirCuts => TheoryLike_CheckPriorCuts
    procedure :: CalculateRequiredTheoryChanges =>TheoryLike_CalculateRequiredTheoryChanges
    procedure :: GetTheoryForLike=>TheoryLike_GetTheoryForLike
    procedure :: GetLogLikeWithTheorySet => TheoryLike_LogLikeWithTheorySet
    end type TheoryLikeCalculator

    contains

    subroutine AddLike(this, CurrentLike, LikeToAdd)
    class(TLikeCalculator) :: this
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
    class(TLikeCalculator) :: this
    class(TCalculationAtParamPoint):: Params
    real(mcp) :: GetLogLikeBounds

    if (any(Params%P(:num_params) > BaseParams%PMax(:num_params)) .or. &
    any(Params%P(:num_params) < BaseParams%PMin(:num_params))) then
        GetLogLikeBounds = logZero
    else
        GetLogLikeBounds=0
    end if

    end function GetLogLikeBounds

    function GetLogPriors(this, P) result(logLike)
    class(TLikeCalculator) :: this
    integer i
    real(mcp), intent(in) :: P(num_params)
    real(mcp) logLike

    logLike=0
    do i=1,num_params
        if (BaseParams%varying(i) .and. BaseParams%GaussPriors%std(i)/=0) then
            logLike = logLike + ((P(i)-BaseParams%GaussPriors%mean(i))/BaseParams%GaussPriors%std(i))**2
        end if
    end do
    logLike=logLike/2

    end function GetLogPriors

    function GetLogLike(this, Params) !Get -Ln(Likelihood) for chains
    class(TLikeCalculator) :: this
    class(TCalculationAtParamPoint) :: Params
    real(mcp) GetLogLike

    GetLogLike = this%GetLogLikeBounds(Params)
    if (GetLogLike==LogZero) return
    if (this%test_likelihood) then
        call this%AddLike(GetLogLike,this%TestLikelihoodFunction(Params))
    else
        call this%AddLike(GetLogLike,this%GetLogLikeMain(Params))
    end if
    if (GetLogLike==LogZero) return
    call this%AddLike(GetLogLike,this%getLogPriors(Params%P))

    end function GetLogLike

    function GetLogLikeMain(this, Params) result(LogLike)!Get -Ln(Likelihood) for chains
    class(TLikeCalculator) :: this
    class(TCalculationAtParamPoint) :: Params
    real(mcp) :: LogLike

    LogLike = LogZero !error - must be overridden

    end function GetLogLikeMain
    
    subroutine TLikeCalculator_ReadParams(this, Ini, initial_propose_matrix)
    class(TLikeCalculator) :: this
    class(TIniFile) :: Ini
    real(mcp) :: initial_propose_matrix(:,:)
    character(LEN=Ini_Max_String_Len) :: covMatrix

    this%test_likelihood = Ini%Read_Logical('test_likelihood', .false.)
    if (this%test_likelihood) then
        print *,'** Using test Gaussian likelihood from covariance + hard priors **'
        covMatrix = trim(Ini%Read_String('test_covariance'))
        allocate(this%test_cov_matrix(num_params_used, num_params_used))
        if (covMatrix=='') then
            this%test_cov_matrix = initial_propose_matrix
        else
            call BaseParams%ReadSetCovMatrix(covMatrix, this%test_cov_matrix)
        end if
    end if

    end subroutine TLikeCalculator_ReadParams

    function TLikeCalculator_TestLikelihoodFunction(this,Params) result(LogLike)
    class(TLikeCalculator) :: this
    class(TCalculationAtParamPoint) Params
    real(mcp) :: LogLike
    real(mcp), allocatable, save :: covInv(:,:)
    real(mcp) X(num_params_used)

    if (.not. allocated(covInv)) then
        allocate(covInv(num_params_used,num_params_used))
        covInv = this%test_cov_matrix
        call Matrix_Inverse(covInv)
    end if
    X = Params%P(params_used) - BaseParams%Center(params_used)
    LogLike = dot_product(X, matmul(covInv, X))/2

    end function TLikeCalculator_TestLikelihoodFunction

    function Generic_GetLogLikeMain(this, Params) result(LogLike)!Get -Ln(Likelihood) for chains
    class(GenericLikeCalculator) :: this
    class(TCalculationAtParamPoint) :: Params
    real(mcp) LogLike

    !Used when you want to plug in your own CMB-independent likelihood function:
    !set generic_mcmc=.true. in settings.f90, then write function here returning -Ln(Likelihood)
    !Parameter array is Params%P, so e.g. 2D unit Gaussian would be
    !GenericLikelihoodFunction = (Params%P(1)**2+Params%P(2)**2)/2
    LogLike = LogZero
    call MpiStop('Generic_GetLogLikeMain: need to write this function!')

    end function Generic_GetLogLikeMain


    subroutine TheoryLike_SetTheoryParams(this, Params)
    class(TheoryLikeCalculator) :: this
    class(TCalculationAtParamPoint), target :: Params

    if (.not. allocated(this%TheoryParams)) call this%Config%NewTheoryParams(this%TheoryParams)
    call this%Config%Parameterization%ParamArrayToTheoryParams(Params%P,this%TheoryParams)
    this%Params => Params

    end subroutine TheoryLike_SetTheoryParams


    function TheoryLike_GetLogLikeMain(this, Params) result(LogLike)!Get -Ln(Likelihood) for chains
    class(TheoryLikeCalculator) :: this
    class(TCalculationAtParamPoint) :: Params
    real(mcp) LogLike

    call this%SetTheoryParams(Params)
    call this%AddLike(LogLike, this%Config%Parameterization%NonBaseParameterPriors(this%TheoryParams))
    if (LogLike == logZero) return
    if (.not. Params%Info%validInfo) then
        this%changeMask(1:num_params) = .true.
    else
        this%changeMask(1:num_params) = Params%Info%lastParamArray(1:num_params)/=Params%P(1:num_params)
    end if
    call this%Params%Theory%AssignNew(this%Params%Theory)
    if (this%CalculateRequiredTheoryChanges()) then
        Params%Info%lastParamArray(1:num_params) = Params%P(1:num_params)
        call this%AddLike(LogLike, this%GetLogLikeWithTheorySet())
    else
        LogLike = logZero
    end if
    if (LogLike==logZero) return

    if (Feedback>2) call DataLikelihoods%WriteLikelihoodContribs(stdout, Params%likelihoods)

    end function TheoryLike_GetLogLikeMain

    function TheoryLike_CheckPriorCuts(this, Params) result(checkPriorCuts)
    class(TheoryLikeCalculator) :: this
    real(mcp)  CheckPriorCuts
    class(TCalculationAtParamPoint) :: Params

    CheckPriorCuts = this%GetLogLikeBounds(Params)
    if (CheckPriorCuts==LogZero) return

    call this%SetTheoryParams(Params)
    CheckPriorCuts = this%Config%Parameterization%NonBaseParameterPriors(this%TheoryParams)

    end function TheoryLike_CheckPriorCuts


    function TheoryLike_GetLogLikePost(this,Params, do_like) result(LogLike)
    !for importance sampling where theory may be pre-stored
    class(TheoryLikeCalculator) :: this
    real(mcp)  LogLike
    class(TCalculationAtParamPoint) :: Params
    logical, optional, intent(in) :: do_like(DataLikelihoods%count)

    LogLike = this%GetLogLikeBounds(Params)
    if (LogLike==LogZero) return

    this%SlowChanged = .false.
    this%ChangeMask = .true.

    call this%SetTheoryParams(Params)
    call this%AddLike(LogLike,this%Config%Parameterization%NonBaseParameterPriors(this%TheoryParams))
    if (LogLike == logZero) return
    call this%AddLike(LogLike, this%GetLogLikeWithTheorySet(do_like))
    if (LogLike == logZero) return
    call this%AddLike(LogLike, this%GetLogPriors(Params%P))

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
