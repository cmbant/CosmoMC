    !Do Metropolis-Hastings and slice sampling algorithms
    !Also directional gridding method for fast-slow parameters
    !Provisional implementation of Wang-Landau and Multicanonical methods

    module MonteCarlo
    use CalcLike
    use RandUtils
    use propose
    use ParamPointSet
    implicit none
    private

    real(mcp), parameter :: eps_1 = 1.00001

    Type, extends(TConfigClass) :: TSampleCollector
    contains
    procedure :: AddNewPoint => TSampleCollector_AddNewPoint
    procedure :: AddNewWeightedPoint => TSampleCollector_AddNewWeightedPoint
    procedure :: ReadCheckpoint => TSampleCollector_ReadCheckpoint
    end Type

    Type, extends(TLikelihoodUser) :: TSamplingAlgorithm
        integer :: num_sample = 0
        real(mcp) :: MaxLike = LogZero
        real(mcp) :: MaxLikeParams(max_num_params)
        class(TSampleCollector), pointer :: SampleCollector => null()
    contains
    procedure :: Init => TSamplingAlgorithm_Init
    procedure :: LogLike => TSamplingAlgorithm_LogLike
    procedure :: SaveState => TSamplingAlgorithm_SaveState
    procedure :: LoadState => TSamplingAlgorithm_LoadState
    end Type TSamplingAlgorithm

    Type, extends(TSamplingAlgorithm) :: TChainSampler
        integer :: oversample_fast =1
        integer :: output_thin =1
        integer :: num_accept = 0
        real(mcp) :: propose_scale = 2.4_mcp
        integer :: burn_in = 2
        class(BlockedProposer), pointer :: Proposer => null()
    contains
    procedure :: MetropolisAccept => TChainSampler_MetropolisAccept
    procedure :: MoveDone => TChainSampler_MoveDone
    procedure :: SampleFrom =>  TChainSampler_SampleFrom
    procedure :: GetNewSample =>  TChainSampler_GetNewSample
    procedure :: SetCovariance => TChainSampler_SetCovariance
    procedure :: SaveState => TChainSampler_SaveState
    procedure :: LoadState => TChainSampler_LoadState
    procedure :: ReadParams => TChainSampler_ReadParams
    procedure :: Init => TChainSampler_Init
    procedure :: InitWithPropose=> TChainSampler_InitWithPropose
    end Type

    Type, extends(TChainSampler) :: TMetropolisSampler
        integer :: num_metropolis_accept=0, num_metropolis=0, last_num=0
    contains
    procedure :: GetNewSample => TMetropolisSampler_GetNewSample
    procedure :: GetNewMetropolisSample => TMetropolisSampler_GetNewSample
    procedure :: FastParameterSample => TMetropolisSampler_FastParameterSample
    end Type

    Type, extends(TMetropolisSampler) :: TFastDraggingSampler
        integer :: num_fast_calls = 0, num_slow_calls = 0, drag_accpt=0
        integer :: num_drag=0
    contains
    procedure :: GetNewSample => TFastDraggingSampler_GetNewSample
    procedure :: ReadParams => TFastDraggingSampler_ReadParams
    end Type

    public TSampleCollector, TSamplingAlgorithm, TChainSampler, TMetropolisSampler, TFastDraggingSampler
    contains

    subroutine TSamplingAlgorithm_Init(this, LikeCalculator, SampleCollector)
    class(TSamplingAlgorithm) :: this
    class(TLikeCalculator), target:: LikeCalculator
    class(TSampleCollector), target, optional :: SampleCollector

    if (present(SampleCollector)) then
        this%SampleCollector => SampleCollector
    else
        this%SampleCollector => null()
    end if
    this%LikeCalculator => LikeCalculator
    this%num_sample = 0
    this%MaxLike = LogZero

    end subroutine TSamplingAlgorithm_Init

    function TSamplingAlgorithm_LogLike(this, Params) result(logLike)
    class(TSamplingAlgorithm) :: this
    class(ParamSet) Params
    real(mcp) :: logLike

    logLike = this%LikeCalculator%GetLogLike(Params)

    end function TSamplingAlgorithm_LogLike

    subroutine TSamplingAlgorithm_LoadState(this,F)
    class(TSamplingAlgorithm) :: this
    class(TFileStream) :: F

    call F%Read(this%num_sample, this%MaxLike)
    call F%Read(this%MaxLikeParams)

    end subroutine TSamplingAlgorithm_LoadState

    subroutine TSamplingAlgorithm_SaveState(this,F)
    class(TSamplingAlgorithm) :: this
    class(TFileStream) :: F

    call F%Write(this%num_sample, this%MaxLike)
    call F%Write(this%MaxLikeParams)

    end subroutine TSamplingAlgorithm_SaveState


    !!!TChainSampler

    function TChainSampler_MetropolisAccept(this, Like, CurLike) result(MetropolisAccept)
    class(TChainSampler) :: this
    real(mcp) Like, CurLike
    logical MetropolisAccept

    if (Like /=LogZero) then
        MetropolisAccept = CurLike > Like
        if (.not. MetropolisAccept) MetropolisAccept = randexp1() > Like - CurLike
    else
        MetropolisAccept = .false.
    end if

    end function TChainSampler_MetropolisAccept


    subroutine TChainSampler_SampleFrom(this,CurParams, CurLike, samples_to_get)
    class(TChainSampler) :: this
    Type(ParamSet) CurParams
    integer samples_to_get
    real(mcp)  CurLike
    real(mcp) :: mult

    mult= 0

    do while (this%num_sample <= samples_to_get*this%Oversample_fast)
        call this%GetNewSample(CurParams, CurLike, mult)

        if (CurLike /= logZero) then
            if (associated(this%SampleCollector) .and. mod(this%num_sample, this%output_thin)==0 .and. mult>0) then
                call this%SampleCollector%AddNewPoint(CurParams%P, CurLike)
                if (mod(this%num_sample,100*this%Oversample_fast)==0) call CheckParamChange
            end if
        else
            if (this%num_sample > 1000) then
                call DoAbort('MCMC.f90: Couldn''t start after 1000 tries - check starting ranges')
            end if
        end if
    end do

    if (Feedback > 0) then
        write(*,*) MPIRank, 'Stopping as have ',samples_to_get ,' samples. '
        call this%LikeCalculator%WritePerformanceStats(stdout)
    end if

    end subroutine TChainSampler_SampleFrom


    subroutine TChainSampler_MoveDone(this,accpt, CurParams, CurLike, mult, thin_fac)
    class(TChainSampler) :: this
    logical, intent(in) :: accpt
    class(TCalculationAtParamPoint), intent(in) :: CurParams
    real(mcp) CurLike
    real(mcp), intent(inout):: mult
    integer, intent(in), optional :: thin_fac

    this%num_sample = this%num_sample + 1
    if (accpt) then
        if (mult>0) then
            if (associated(this%SampleCollector) .and. CurLike /= LogZero .and. this%num_accept> this%burn_in) &
            & call this%SampleCollector%AddNewWeightedPoint(CurParams, CurLike, mult, thin_fac)
            this%num_accept = this%num_accept + 1
        end if
        mult=1
        if (CurLike < this%MaxLike) then
            this%MaxLike = CurLike
            this%MaxLikeParams = CurParams%P
        end if
    else
        mult = mult + 1
    end if

    end subroutine TChainSampler_MoveDone

    subroutine TChainSampler_GetNewSample(this, CurParams, CurLike, mult)
    class(TChainSampler) :: this
    Type(ParamSet) CurParams
    real(mcp) CurLike, mult
    end subroutine TChainSampler_GetNewSample


    subroutine TChainSampler_LoadState(this,F)
    class(TChainSampler) :: this
    class(TFileStream) F

    call this%TSamplingAlgorithm%LoadState(F)
    call F%Read(this%num_accept)
    call this%Proposer%LoadState(F)

    end subroutine TChainSampler_LoadState


    subroutine TChainSampler_SaveState(this,F)
    class(TChainSampler) :: this
    class(TFileStream) F

    call this%TSamplingAlgorithm%SaveState(F)
    call F%Write(this%num_accept)
    call this%Proposer%SaveState(F)

    end subroutine TChainSampler_SaveState

    subroutine TChainSampler_SetCovariance(this, cov)
    class(TChainSampler) :: this
    real(mcp) :: cov(:,:)

    call this%Proposer%SetCovariance(cov)

    end subroutine TChainSampler_SetCovariance


    subroutine TChainSampler_ReadParams(this, Ini)
    class(TChainSampler) :: this
    class(TSettingIni) :: Ini

    if (use_fast_slow) then
        this%oversample_fast = Ini%Read_Int('oversample_fast',1, min=1)
        this%output_thin =  this%oversample_fast
    end if
    call Ini%Read('propose_scale',this%propose_scale, min=0.d0)
    call Ini%Read('burn_in',this%burn_in)

    end subroutine TChainSampler_ReadParams

    subroutine TChainSampler_Init(this, LikeCalculator, SampleCollector)
    class(TChainSampler) :: this
    class(TLikeCalculator), target:: LikeCalculator
    class(TSampleCollector), target, optional :: SampleCollector

    call this%InitWithPropose(LikeCalculator, SampleCollector)

    end subroutine TChainSampler_Init

    subroutine TChainSampler_InitWithPropose(this,LikeCalculator, SampleCollector, propose_scale)
    class(TChainSampler) :: this
    class(TLikeCalculator), target:: LikeCalculator
    class(TSampleCollector), target, optional :: SampleCollector
    real(mcp), intent(in), optional :: propose_scale

    if (present(propose_scale)) this%propose_scale = propose_scale

    call this%TSamplingAlgorithm%Init(LikeCalculator,SampleCollector)
    allocate(BlockedProposer::this%Proposer)
    call this%Proposer%Init(BaseParams%param_blocks, slow_block_max= slow_tp_max, &
    oversample_fast=this%oversample_fast, propose_scale=this%propose_scale)
    this%num_accept=0

    end subroutine TChainSampler_InitWithPropose

    !!! TMetropolisSampler

    subroutine TMetropolisSampler_GetNewSample(this, CurParams, CurLike, mult)
    !Standard metropolis hastings
    class(TMetropolisSampler) :: this
    Type(ParamSet) CurParams, Trial
    real(mcp) CurLike, mult, Like
    logical :: accpt

    Trial = CurParams
    call this%Proposer%GetProposal(Trial%P)

    Like = this%LogLike(Trial)
    if (Feedback > 1) write (*,*) instance, 'Likelihood: ', real(Like), 'Current Like:', real(CurLike)

    if (Like /= logZero) then
        accpt = this%MetropolisAccept(Like, CurLike)
    else
        accpt = .false.
    end if

    this%num_metropolis = this%num_metropolis + 1

    call this%MoveDone(accpt, CurParams, CurLike, mult, this%Oversample_fast)
    call CurParams%AcceptReject(Trial, accpt)

    if (accpt) then
        CurParams = Trial
        CurLike = Like
        this%num_metropolis_accept = this%num_metropolis_accept + 1
        if (Feedback > 1) write (*,*) this%num_metropolis, ' metropolis accept. ratio:', &
        & real(this%num_metropolis_accept)/this%num_metropolis
        if (LogFile%Opened() .and. mod(this%num_metropolis_accept,50*this%Oversample_fast) ==0) then
            write (LogFile%unit,*) 'metropolis rat:',real(this%num_metropolis_accept)/this%num_metropolis,  &
            & ' in ',this%num_metropolis, ', best: ',real(this%MaxLike)
            write (LogFile%unit,*) 'local acceptance ratio:', 50./(this%num_metropolis - this%last_num)
            this%last_num = this%num_metropolis
        end if
    end if

    end subroutine TMetropolisSampler_GetNewSample

    subroutine TMetropolisSampler_FastParameterSample(this, CurParams, CurLike, mult)
    !Metropolis hastings on fast parameters
    class(TMetropolisSampler) :: this
    Type(ParamSet) CurParams, Trial
    real(mcp) CurLike, Like
    real(mcp) mult
    logical :: accpt

    Trial = CurParams
    call this%Proposer%GetProposalFast(Trial%P)

    Like = this%LogLike(Trial)
    if (Like /= logZero) then
        accpt = this%MetropolisAccept(Like, CurLike)
    else
        accpt = .false.
    end if

    call this%MoveDone(accpt, CurParams, CurLike, mult)
    call CurParams%AcceptReject(Trial,accpt)

    if (accpt) then
        CurParams = Trial
        CurLike = Like
    end if

    end subroutine TMetropolisSampler_FastParameterSample


    subroutine TFastDraggingSampler_GetNewSample(this,CurParams, CurLike, mult)
    !Make proposals in fast-marginalized slow parameters
    !'drag' fast parameters using method of Neal
    class(TFastDraggingSampler) :: this
    Type(ParamSet) TrialEnd, TrialStart, CurParams, CurEndParams, CurStartParams
    real(mcp) CurLike, mult
    integer numaccpt
    real(mcp) CurIntLike, IntLike, CurEndLike, CurStartLike, EndLike, StartLike
    real(mcp) CurDragLike, DragLike
    logical :: accpt
    real(mcp)  :: likes_start_sum, likes_end_sum
    integer interp_step, interp_steps
    real(mcp) frac, delta(num_params)

    if (CurLike == LogZero .or. BaseParams%num_fast==0 .or. BaseParams%num_slow ==0) then
        call this%GetNewMetropolisSample(CurParams, CurLike, mult)
        return
    end if

    this%num_drag=this%num_drag+1
    if (mod(this%num_drag, this%oversample_fast)/=0) then
        call this%TMetropolisSampler%FastParameterSample(CurParams, CurLike, mult)
        return
    end if

    if (Feedback > 1) write (*,*) instance, 'Fast dragging, Like: ', CurLike

    call Timer()
    TrialEnd = CurParams
    call this%Proposer%GetProposalSlow(TrialEnd%P)

    CurEndLike = this%LogLike(TrialEnd)
    if (CurEndLike==logZero) then
        call TrialEnd%Clear(keep=CurParams)
        mult = mult + 1
        return
    end if
    CurStartLike = CurLike
    if (Feedback > 1) call Timer('Dragging Slow time')

    this%num_slow_calls = this%num_slow_calls + 1

    likes_end_sum = CurEndLike
    likes_start_sum = CurStartLike

    CurStartParams = CurParams
    CurEndParams = TrialEnd

    interp_steps = max(2,nint(dragging_steps * BaseParams%num_fast) + 1)

    numaccpt = 0
    do interp_step = 1, interp_steps-1
        call this%Proposer%GetProposalFastDelta(delta)
        TrialEnd = CurEndParams
        TrialEnd%P(1:num_params) = TrialEnd%P(1:num_params) + delta
        EndLike = this%LogLike(TrialEnd)
        accpt = EndLike /= logZero

        if (accpt) then
            this%num_fast_calls = this%num_fast_calls + 1
            TrialStart = CurStartParams
            TrialStart%P(1:num_params) = TrialStart%P(1:num_params)  + delta
            StartLike = this%LogLike(TrialStart)
            accpt = StartLike/=logZero

            if (accpt) then
                this%num_fast_calls = this%num_fast_calls + 1
                if (Feedback > 2) print *,'End,start drag: ', interp_step, EndLike, StartLike

                frac = real(interp_step, mcp)/interp_steps
                CurIntLike = CurStartLike*(1-frac) + frac*CurEndLike
                IntLike = StartLike*(1-frac) + frac*EndLike
                accpt = this%MetropolisAccept(IntLike, CurIntLike)
            end if
        end if

        call CurEndParams%AcceptReject(TrialEnd, accpt)
        call CurStartParams%AcceptReject(TrialStart, accpt)

        if (accpt) then
            CurEndParams = TrialEnd
            CurStartParams = TrialStart
            CurEndLike = EndLike
            CurStartLike = StartLike
            CurIntLike = IntLike
            numaccpt = numaccpt + 1
        end if

        likes_start_sum = likes_start_sum + CurStartLike
        likes_end_sum = likes_end_sum + CurEndLike
    end do

    if (Feedback > 1) call Timer('Dragging time')

    if (Feedback > 1) print *,'drag steps accept ratio:', real(numaccpt)/(interp_steps)

    CurDragLike = likes_start_sum/interp_steps !old slow
    DragLike = likes_end_sum/interp_steps !proposed new slow
    if (Feedback > 2) print *,'CurDragLike, DragLike: ', CurDragLike, DragLike

    accpt = this%MetropolisAccept(DragLike, CurDragLike)

    call this%MoveDone(accpt, CurParams, CurLike, mult)
    call CurParams%AcceptReject(CurEndParams, accpt)

    if (accpt) then
        CurParams = CurEndParams
        CurLike = CurEndLike
        this%drag_accpt=this%drag_accpt+1
        if (Feedback > 0 .and. mod(this%drag_accpt,30)==0 .or. Feedback>1) &
        write (*,*) trim(concat('Chain:',MpiRank,' drag accpt:')), real(this%drag_accpt)/(this%num_drag/this%oversample_fast), &
        'fast/slow',real(this%num_fast_calls)/this%num_slow_calls, 'slow:', this%num_slow_calls
    end if

    end subroutine TFastDraggingSampler_GetNewSample


    subroutine TFastDraggingSampler_ReadParams(this, Ini)
    class(TFastDraggingSampler) :: this
    class(TSettingIni) :: Ini

    call this%TMetropolisSampler%ReadParams(Ini)
    this%output_thin = 1

    end subroutine TFastDraggingSampler_ReadParams


    subroutine TSampleCollector_AddNewPoint(this,P,like)
    class(TSampleCollector) :: this
    real(mcp), intent(in) ::like
    real(mcp) P(:)
    !Add each point (duplicated if chain does not move)
    end subroutine TSampleCollector_AddNewPoint

    subroutine TSampleCollector_AddNewWeightedPoint(this, CurParams, CurLike, mult, thin_fac)
    class(TSampleCollector) :: this
    class(TCalculationAtParamPoint), intent(in) :: CurParams
    real(mcp) CurLike
    real(mcp), intent(in):: mult !tbe weight, usually integer for standard chains
    integer, intent(in), optional :: thin_fac
    !Add samples accumulated into weighted sample
    end subroutine TSampleCollector_AddNewWeightedPoint

    subroutine TSampleCollector_ReadCheckpoint(this)
    class(TSampleCollector) :: this
    end subroutine TSampleCollector_ReadCheckpoint


    end module MonteCarlo
