    module GeneralSetup
    use CalcLike
    use MonteCarlo
    use GeneralTypes
    use BaseParameters
    use SampleCollector
    use ImportanceSampling
    use minimize
    use IO
    use ParamPointSet
    implicit none

    integer, parameter :: action_MCMC=0, action_importance=1, action_maxlike=2, action_hessian=3, action_tests=4

    Type :: TSetup
        integer :: action = action_MCMC
        integer :: sampling_method = sampling_metropolis
        integer :: samples
        class(TGeneralConfig), allocatable:: Config
        class(TLikeCalculator), allocatable:: LikeCalculator
        class(TSampleCollector), allocatable :: SampleCollector
        class(TSamplingAlgorithm), allocatable :: SamplingAlgorithm
        class(TImportanceSampler), allocatable :: ImportanceSampler
    contains
    procedure :: Init => TSetup_Init
    procedure :: ReadParams  => TSetup_ReadParams
    procedure :: DoneInitialize => TSetup_DoneInitialize !Called again after other things configured
    procedure :: GetMinimizer => TSetup_GetMinimizer
    procedure :: DoSampling => TSetup_DoSampling
    procedure :: DoTests => TSetup_DoTests
    end type

    class(TSetup), allocatable :: Setup

    contains

    subroutine TSetup_Init(this)
    class(TSetup) :: this

    if (.not. allocated(this%LikeCalculator)) allocate(TGenericLikeCalculator::this%LikeCalculator)
    if (.not. allocated(this%Config)) allocate(TGeneralConfig::this%Config)

    end subroutine TSetup_Init


    subroutine TSetup_ReadParams(this, Ini)
    !Called after Init
    class(TSetup) :: this
    class(TSettingIni) :: Ini

    call this%Config%ReadParams(Ini)
    this%action = Ini%Read_Int('action',action_MCMC,min=0,max=action_tests)
    if (this%action/=action_importance) use_fast_slow = Ini%Read_Logical('use_fast_slow',.true.)

    call this%LikeCalculator%InitWithParams(Ini, this%Config)

    if (this%action==action_MCMC) then
        allocate(TMpiChainCollector::this%SampleCollector)
        call this%SampleCollector%InitWithParams(Ini, this%Config)
        this%sampling_method = Ini%Read_Int('sampling_method',sampling_metropolis)
        if (this%sampling_method == sampling_metropolis) then
            allocate(TMetropolisSampler::this%SamplingAlgorithm)
        else if (this%sampling_method == sampling_fast_dragging) then
            allocate(TFastDraggingSampler::this%SamplingAlgorithm)
            use_fast_slow = .true.
        else
            call MpiStop('Sampling method not currently supported')
        end if
        call this%SamplingAlgorithm%ReadParams(Ini)
        this%samples = Ini%Read_Int('samples')
    else if (this%action == action_importance) then
        allocate(TImportanceSampler::this%ImportanceSampler)
        call this%ImportanceSampler%ReadParams(Ini)
        if (this%ImportanceSampler%redo_theory) call this%Config%Calculator%ReadImportanceParams(Ini)
    end if

    end subroutine TSetup_ReadParams


    subroutine TSetup_DoneInitialize(this)
    !Called after likelihoods, base parameters, etc initialized
    class(TSetup), target :: this

    if (allocated(this%SamplingAlgorithm)) then
        call this%SamplingAlgorithm%Init(this%LikeCalculator, this%SampleCollector)
        select type (Sampler=>this%SamplingAlgorithm)
        class is (TChainSampler)
            if (allocated(this%SampleCollector)) then
                select type (Collector=>this%SampleCollector)
                class is (TMpiChainCollector)
                    Collector%Sampler => Sampler
                    if (BaseParams%covariance_has_new) &
                    & Collector%MPi%MPI_Max_R_ProposeUpdate = Collector%Mpi%MPI_Max_R_ProposeUpdateNew
                end select
            end if
            call Sampler%SetCovariance(BaseParams%covariance_estimate)
        end select
    end if

    if (allocated(this%ImportanceSampler)) then
        call this%ImportanceSampler%Init(this%LikeCalculator)
    end if

    end subroutine TSetup_DoneInitialize

    subroutine TSetup_GetMinimizer(this, Minimizer)
    class(TSetup), target :: this
    class(TMinimizer), allocatable :: Minimizer

    allocate(TPowellMinimizer::Minimizer)
    Minimizer%LikeCalculator=>this%LikeCalculator

    end subroutine TSetup_GetMinimizer

    subroutine TSetup_DoSampling(this, Params)
    class(TSetup), target :: this
    class(ParamSet) :: Params
    real(mcp) mult
    real(mcp) StartLike
    character(LEN=:), allocatable :: fname

    select type (Sampler=>this%SamplingAlgorithm)
    class is (TChainSampler)
        if (new_chains) then
            call BaseParams%SetStartPositions(Params)
            StartLike=LogZero
        else
            call IO_ReadLastChainParams(rootname//'.txt', mult, StartLike, Params%P, params_used)
            if (allocated(this%SampleCollector)) then
                call this%SampleCollector%ReadCheckpoint()
            end if
        end if
        fname = rootname//'.txt'
        call ChainOutFile%CreateOpenFile(fname, append = .not. new_chains)

        if (Feedback > 0 .and. MPIRank==0) write (*,*) 'starting Monte-Carlo'
        call Sampler%SampleFrom(Params, StartLike, this%samples)
        call ChainOutFile%Close()
        class default
        call MpiStop('Sampling not implemented')
    end select
    if (Feedback > 0) write (*,*) 'finished'

    end subroutine TSetup_DoSampling

    subroutine TSetup_DoTests(this, output_root)
    !This runs likelihoods for fixed central values of parameters and outputs the likelihoods
    ! - e.g. for likelihood testing between versions, etc.
    class(TSetup) :: this
    character(LEN=*) :: output_root
    Type(ParamSet) :: Params
    real(mcp) :: logLike
    Type(TTimer) Timer

    call Timer%Start()
    Params%P(:num_params) = BaseParams%Center(:num_params)
    this%LikeCalculator%timing = .true.
    logLike = this%LikeCalculator%GetLogLike(Params)
    write(*,*) '   loglike     chi-sq'
    if (Feedback <=2) call DataLikelihoods%WriteLikelihoodContribs(stdout, Params%likelihoods)
    write(*,*) 'Test likelihoods done, total logLike = '//RealToStr(logLike)//', chi-sq = '//RealToStr(logLike*2)
    write(*,*) 'Likelihood calculation time (seconds)= '//RealToStr(Timer%Time())
    if (allocated(Params%Theory) .and. output_root/='') then
        call DataLikelihoods%WriteDataForLikelihoods(Params%P, Params%Theory, output_root)
        call Params%Theory%WriteTextData(output_root)
    end if
    call DoStop()

    end subroutine TSetup_DoTests

    end module GeneralSetup
