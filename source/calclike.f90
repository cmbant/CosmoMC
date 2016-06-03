    module CalcLike
    use DataLikelihoodList
    use GeneralTypes
    use BaseParameters
    use MatrixUtils
    use MiscUtils
    implicit none
    private

    Type, extends(TConfigClass) :: TLikeCalculator
        real(mcp) :: Temperature =1
        logical :: test_likelihood= .false.
        logical :: timing = .false.
        real(mcp), allocatable :: test_cov_matrix(:,:)
    contains
    procedure :: AddLike
    procedure :: AddLikeTemp
    procedure :: GetLogLikeBounds
    procedure :: GetLogPriors
    procedure :: GetLogLike
    procedure :: GetLogLikeMain
    procedure :: ReadParams => TLikeCalculator_ReadParams
    procedure :: TestLikelihoodFunction => TLikeCalculator_TestLikelihoodFunction
    procedure :: WritePerformanceStats => TLikeCalculator_WritePerformanceStats
    procedure :: WriteParamsHumanText_unit => TLikeCalculator_WriteParamsHumanText
    procedure :: WriteParamsHumanText_file => TLikeCalculator_WriteParamsHumanTextFile
    procedure :: WriteParamPointTextData => TLikeCalculator_WriteParamPointTextData
    generic :: WriteParamsHumanText => WriteParamsHumanText_unit, WriteParamsHumanText_file
    end type TLikeCalculator

    type, extends(TLikeCalculator) :: TGenericLikeCalculator
    contains
    procedure :: GetLogLikeMain => Generic_GetLogLikeMain
    end type TGenericLikeCalculator

    type, extends(TLikeCalculator) :: TTheoryLikeCalculator
        class(TTheoryParams), allocatable :: TheoryParams
        class(TCalculationAtParamPoint), pointer :: Params
        logical changeMask(max_num_params)
        logical SlowChanged, SemiSlowchanged
    contains
    procedure :: GetLogLikeMain => TheoryLike_GetLogLikeMain
    procedure :: GetLogLikePost => TheoryLike_GetLogLikePost
    procedure :: SetTheoryParams => TheoryLike_SetTheoryParams
    procedure :: CheckPriorCuts => TheoryLike_CheckPriorCuts
    procedure :: CalculateRequiredTheoryChanges =>TheoryLike_CalculateRequiredTheoryChanges
    procedure :: GetTheoryForLike=>TheoryLike_GetTheoryForLike
    procedure :: GetTheoryForImportance=>TheoryLike_GetTheoryForImportance
    procedure :: GetLogLikeWithTheorySet => TheoryLike_LogLikeWithTheorySet
    procedure :: UpdateTheoryForLikelihoods => TheoryLike_UpdateTheoryForLikelihoods
    procedure :: SetNewTheoryResults => TheoryLike_SetNewTheoryResults
    procedure :: TestLikelihoodFunction => TheoryLike_TestLikelihoodFunction
    procedure :: WriteParamsHumanText_unit => TTheoryLike_WriteParamsHumanText
    procedure :: WriteParamPointTextData => TTheoryLike_WriteParamPointTextData
    end type TTheoryLikeCalculator

    type, extends(TCheckpointable) :: TLikelihoodUser
        class(TLikeCalculator), pointer :: LikeCalculator => null()
    end type

    type, extends(TCheckpointable) :: TTheoryLikelihoodUser
        class(TTheoryLikeCalculator), pointer :: LikeCalculator => null()
    end type

    public TLikeCalculator, TGenericLikeCalculator, TTheoryLikeCalculator, TLikelihoodUser, TTheoryLikelihoodUser
    contains

    subroutine AddLike(this, CurrentLike, LikeToAdd)
    class(TLikeCalculator) :: this
    real(mcp), intent(in) :: LikeToAdd
    real(mcp) CurrentLike

    if (CurrentLike/=LogZero) then
        if (LikeToAdd == logZero) then
            CurrentLike = LogZero
        else
            CurrentLike = CurrentLike + LikeToAdd
        end if
    end if
    end subroutine AddLike

    subroutine AddLikeTemp(this, CurrentLike, LikeToAdd)
    class(TLikeCalculator) :: this
    real(mcp), intent(in) :: LikeToAdd
    real(mcp) CurrentLike

    if (CurrentLike/=LogZero) then
        if (LikeToAdd == logZero) then
            CurrentLike = LogZero
        else
            CurrentLike = CurrentLike + LikeToAdd/this%Temperature
        end if
    end if
    end subroutine AddLikeTemp


    function GetLogLikeBounds(this,Params)
    class(TLikeCalculator) :: this
    class(TCalculationAtParamPoint):: Params
    real(mcp) :: GetLogLikeBounds

    if (any(Params%P(:num_params) > BaseParams%PMax(:num_params)) .or. &
        & any(Params%P(:num_params) < BaseParams%PMin(:num_params))) then
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
        if ((BaseParams%varying(i) .or. BaseParams%include_fixed_parameter_priors) &
            .and. BaseParams%GaussPriors%std(i)/=0) then
        logLike = logLike + ((P(i)-BaseParams%GaussPriors%mean(i))/BaseParams%GaussPriors%std(i))**2
        end if
    end do

    do i= 1, size(BaseParams%LinearCombinations)
        associate(Comb => BaseParams%LinearCombinations(i))
            if (Comb%std/=0) then
                logLike = logLike + ((dot_product(Comb%Combination,P) -Comb%mean)/Comb%std)**2
            end if
        end associate
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
        call this%AddLikeTemp(GetLogLike,this%TestLikelihoodFunction(Params))
    else
        call this%AddLikeTemp(GetLogLike,this%GetLogLikeMain(Params))
    end if
    if (GetLogLike==LogZero) return
    call this%AddLikeTemp(GetLogLike,this%getLogPriors(Params%P))

    end function GetLogLike

    function GetLogLikeMain(this, Params) result(LogLike)!Get -Ln(Likelihood) for chains
    class(TLikeCalculator) :: this
    class(TCalculationAtParamPoint) :: Params
    real(mcp) :: LogLike

    LogLike = LogZero !error - must be overridden

    end function GetLogLikeMain

    subroutine TLikeCalculator_ReadParams(this, Ini)
    class(TLikeCalculator) :: this
    class(TSettingIni) :: Ini
    character(LEN=:), pointer :: covMatrix

    this%test_likelihood = Ini%Read_Logical('test_likelihood', .false.)
    if (this%test_likelihood) then
        print *,'** Using test Gaussian likelihood from covariance + hard priors **'
        covMatrix => Ini%Read_String('test_covariance')
        if (covMatrix/='') then
            allocate(this%test_cov_matrix(num_params_used, num_params_used))
            call BaseParams%ReadSetCovMatrix(covMatrix, this%test_cov_matrix)
        end if
    end if
    call Ini%Read('temperature',this%Temperature)

    end subroutine TLikeCalculator_ReadParams

    function TLikeCalculator_TestLikelihoodFunction(this,Params) result(LogLike)
    class(TLikeCalculator) :: this
    class(TCalculationAtParamPoint) Params
    real(mcp) :: LogLike
    real(mcp), allocatable, save :: covInv(:,:)
    real(mcp) X(num_params_used)

    if (.not. allocated(covInv)) then
        allocate(covInv(num_params_used,num_params_used))
        if (.not. allocated(this%test_cov_matrix)) then
            covInv = BaseParams%covariance_estimate
        else
            covInv = this%test_cov_matrix
        endif
        call Matrix_Inverse(covInv)
    end if
    X = Params%P(params_used) - BaseParams%Center(params_used)
    LogLike = dot_product(X, matmul(covInv, X))/2

    end function TLikeCalculator_TestLikelihoodFunction

    subroutine TLikeCalculator_WritePerformanceStats(this, unit)
    class(TLikeCalculator) :: this
    integer, intent(in) :: unit

    end subroutine TLikeCalculator_WritePerformanceStats


    subroutine TLikeCalculator_WriteParamsHumanText(this, aunit, P, LogLike, weight)
    class(TLikeCalculator) :: this
    class(TCalculationAtParamPoint) P
    real(mcp), intent(in), optional :: LogLike,weight
    integer, intent(in) :: aunit
    integer isused,i

    if (present(weight)) then
        write (aunit,*) ' weight    = ',weight
    end if

    if (present(LogLike)) then
        write (aunit,*) '-log(Like) = ',LogLike
        write (aunit,*) ' chi-sq    = ',LogLike*2
        write (aunit,*) ''
    end if

    do isused = 0,1
        do i=1, num_params
            if (isused==0 .and. BaseParams%varying(i) .or. isused==1 .and. .not. BaseParams%varying(i)) then
                write(aunit,'(1I5,1E15.7,"   ",1A22)', advance='NO') &
                    i, P%P(i), BaseParams%NameMapping%name(i)
                write (aunit,'(a)') trim(BaseParams%NameMapping%label(i))
            end if
        end do
        write (aunit,*) ''
    end do

    end subroutine TLikeCalculator_WriteParamsHumanText


    subroutine TLikeCalculator_WriteParamsHumanTextFile(this, fname, P, LogLike, weight)
    class(TLikeCalculator) :: this
    character(LEN=*), intent(in) :: fname
    class(TCalculationAtParamPoint) P
    real(mcp), intent(in), optional :: LogLike, weight
    Type(TTextFile) F

    call F%CreateFile(fname)
    call this%WriteParamsHumanText(F%unit, P, LogLike, weight)
    call F%Close()

    end subroutine TLikeCalculator_WriteParamsHumanTextFile


    subroutine TLikeCalculator_WriteParamPointTextData(this, output_root, Params)
    class(TLikeCalculator) :: this
    character(LEN=*), intent(in) :: output_root
    class(TCalculationAtParamPoint) Params
    end subroutine TLikeCalculator_WriteParamPointTextData


    function Generic_GetLogLikeMain(this, Params) result(LogLike)!Get -Ln(Likelihood) for chains
    class(TGenericLikeCalculator) :: this
    class(TCalculationAtParamPoint) :: Params
    real(mcp) LogLike

    !Used when you want to plug in your own CMB-independent likelihood function:
    !Parameter array is Params%P, so e.g. 2D unit Gaussian would be
    LogLike = (Params%P(1)**2+Params%P(2)**2)/2
    !LogLike = LogZero
    !call MpiStop('Generic_GetLogLikeMain: need to write this function!')

    end function Generic_GetLogLikeMain


    subroutine TheoryLike_SetTheoryParams(this, Params)
    class(TTheoryLikeCalculator) :: this
    class(TCalculationAtParamPoint), target :: Params

    if (.not. allocated(this%TheoryParams)) call this%Config%Parameterization%NewTheoryParams(this%TheoryParams)
    call this%Config%Parameterization%ParamArrayToTheoryParams(Params%P,this%TheoryParams)
    this%Params => Params

    end subroutine TheoryLike_SetTheoryParams

    subroutine TheoryLike_SetNewTheoryResults(this, Params)
    class(TTheoryLikeCalculator) :: this
    class(TCalculationAtParamPoint) :: Params

    if (.not. allocated(Params%Theory)) call this%Config%NewTheory(Params%Theory)

    end subroutine TheoryLike_SetNewTheoryResults


    function TheoryLike_GetLogLikeMain(this, Params) result(LogLike)
    !Get -Ln(Likelihood), not accounting for temperature
    class(TTheoryLikeCalculator) :: this
    class(TCalculationAtParamPoint) :: Params
    real(mcp) LogLike

    call this%SetTheoryParams(Params)
    LogLike = this%Config%Parameterization%NonBaseParameterPriors(this%TheoryParams)
    if (LogLike == logZero) return
    if (.not. Params%validInfo) then
        this%changeMask(1:num_params) = .true.
    else
        this%changeMask(1:num_params) = Params%lastParamArray(1:num_params)/=Params%P(1:num_params)
    end if
    call this%SetNewTheoryResults(Params)
    if (this%CalculateRequiredTheoryChanges()) then
        call this%AddLike(LogLike, this%GetLogLikeWithTheorySet())
        Params%lastParamArray(1:num_params) = Params%P(1:num_params)
    else
        LogLike = logZero
    end if
    if (LogLike==logZero) return

    if (Feedback>2) call DataLikelihoods%WriteLikelihoodContribs(stdout, Params%likelihoods)

    end function TheoryLike_GetLogLikeMain

    function TheoryLike_CheckPriorCuts(this, Params) result(checkPriorCuts)
    class(TTheoryLikeCalculator) :: this
    real(mcp)  CheckPriorCuts
    class(TCalculationAtParamPoint) :: Params

    CheckPriorCuts = this%GetLogLikeBounds(Params)
    if (CheckPriorCuts==LogZero) return

    call this%SetTheoryParams(Params)
    CheckPriorCuts = this%Config%Parameterization%NonBaseParameterPriors(this%TheoryParams)

    end function TheoryLike_CheckPriorCuts


    function TheoryLike_GetLogLikePost(this,Params, do_like) result(LogLike)
    !for importance sampling where theory may be pre-stored
    class(TTheoryLikeCalculator) :: this
    real(mcp)  LogLike
    class(TCalculationAtParamPoint) :: Params
    logical, optional, intent(in) :: do_like(DataLikelihoods%count)

    LogLike = this%GetLogLikeBounds(Params)
    if (LogLike==LogZero) return

    this%SlowChanged = .false.
    this%ChangeMask = .true.

    call this%SetTheoryParams(Params)
    call this%AddLikeTemp(LogLike,this%Config%Parameterization%NonBaseParameterPriors(this%TheoryParams))
    if (LogLike == logZero) return
    call this%AddLikeTemp(LogLike, this%GetLogLikeWithTheorySet(do_like))
    if (LogLike == logZero) return
    call this%AddLikeTemp(LogLike, this%GetLogPriors(Params%P))

    end function TheoryLike_GetLogLikePost


    function TheoryLike_LogLikeWithTheorySet(this, likelihood_mask) result(logLike)
    class(TTheoryLikeCalculator) :: this
    logical, intent(in), optional :: likelihood_mask(DataLikelihoods%count)
    real(mcp) logLike
    real(mcp) itemLike
    class(TDataLikelihood), pointer :: like => null()
    integer i
    logical :: do_like(DataLikelihoods%count)
    Type(TTimer) Timer

    if (present(likelihood_Mask)) then
        do_like = likelihood_mask
    else
        do_like = .true.
    end if
    logLike = logZero
    call this%GetTheoryForLike(null()) !chance to initalize
    do i= 1, DataLikelihoods%count
        if (do_like(i)) then
            like => DataLikelihoods%Item(i)
            if (any(like%dependent_params(1:num_params) .and. this%changeMask(1:num_params) )) then
                call this%GetTheoryForLike(like)
                if (this%timing) call Timer%Start
                itemLike = like%GetLogLike(this%TheoryParams, this%Params%Theory, this%Params%P(like%nuisance_indices))
                if (this%timing) call Timer%WriteTime('Time for '//trim(like%name))
                if (itemLike == logZero) return
                this%Params%Likelihoods(i) = itemLike
            end if
        end if
    end do
    logLike = sum(this%Params%likelihoods(1:DataLikelihoods%Count))

    end function TheoryLike_LogLikeWithTheorySet


    logical function TheoryLike_CalculateRequiredTheoryChanges(this)
    class(TTheoryLikeCalculator) :: this

    !Set this%Params theory entries (this%Params%Theory) as required for current this%changeMask and this%TheoryParams
    this%Params%validInfo = .true.
    TheoryLike_CalculateRequiredTheoryChanges = .true.

    end function TheoryLike_CalculateRequiredTheoryChanges

    subroutine TheoryLike_GetTheoryForLike(this,Like)
    class(TTheoryLikeCalculator) :: this
    class(TDataLikelihood), pointer :: like

    !If needed, likelihood specific calculation/initalization; like=null for first initial call
    end subroutine TheoryLike_GetTheoryForLike

    subroutine TheoryLike_GetTheoryForImportance(this,Params, error)
    class(TTheoryLikeCalculator) :: this
    class(TCalculationAtParamPoint), target :: Params
    integer error

    call this%SetTheoryParams(Params)
    call this%Config%Calculator%GetTheoryForImportance(this%TheoryParams, Params%Theory, error)

    end subroutine TheoryLike_GetTheoryForImportance


    subroutine TheoryLike_UpdateTheoryForLikelihoods(this, Params)
    class(TTheoryLikeCalculator) :: this
    class(TCalculationAtParamPoint) :: Params

    end subroutine TheoryLike_UpdateTheoryForLikelihoods

    function TheoryLike_TestLikelihoodFunction(this,Params) result(LogLike)
    class(TTheoryLikeCalculator) :: this
    class(TCalculationAtParamPoint) Params
    real(mcp) :: LogLike

    call this%SetNewTheoryResults(Params) !so derived parameters can be output OK
    LogLike = this%TLikeCalculator%TestLikelihoodFunction(Params)

    end function TheoryLike_TestLikelihoodFunction


    subroutine TTheoryLike_WriteParamsHumanText(this, aunit, P, LogLike, weight)
    class(TTheoryLikeCalculator) :: this
    class(TCalculationAtParamPoint) P
    real(mcp), intent(in), optional :: LogLike, weight
    integer, intent(in) :: aunit
    real(mcp), allocatable :: derived(:)
    integer :: numderived = 0
    integer i

    call this%TLikeCalculator%WriteParamsHumanText(aunit, P, LogLike, weight)

    call this%Config%Parameterization%CalcDerivedParams(P%P,P%Theory, derived)
    call DataLikelihoods%addLikelihoodDerivedParams(P%P, P%Theory, derived, P%Likelihoods, LogLike)
    if (allocated(derived)) numderived = size(derived)
    do i=1, numderived
        write(aunit,'(1I5,1E15.7,"   ",1A22)', advance='NO') &
            num_params+i, derived(i), BaseParams%NameMapping%name(num_params + i )
        write (aunit,'(a)') trim(BaseParams%NameMapping%label(num_params+i))
    end do

    if (present(LogLike)) then
        write(aunit,*) ''
        write(aunit,*) '-log(Like)     chi-sq   data'
        call DataLikelihoods%WriteLikelihoodContribs(aunit, P%likelihoods)
    end if

    end subroutine TTheoryLike_WriteParamsHumanText

    subroutine TTheoryLike_WriteParamPointTextData(this, output_root, Params)
    class(TTheoryLikeCalculator) :: this
    character(LEN=*), intent(in) :: output_root
    class(TCalculationAtParamPoint) Params

    if (allocated(Params%Theory)) then
        call DataLikelihoods%WriteDataForLikelihoods(Params%P, Params%Theory, output_root)
        call Params%Theory%WriteTextData(output_root)
    end if

    end subroutine TTheoryLike_WriteParamPointTextData


    end module CalcLike
