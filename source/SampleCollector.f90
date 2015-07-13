
    module SampleCollector
    use propose
    use BaseParameters
    use likelihood
    use GeneralTypes
    use Samples
    use MonteCarlo
    use IO
    implicit none
    private

    Type TMPIData
        real  :: MPI_R_Stop = 0.05
        integer :: MPI_thin_fac = 1
        !following are numbers of samples / (number of parameters)
        !MPI_Min_Sample_Update*num_fast is min number of chain steps
        !before starting to update covmat and checking convergence
        integer :: MPI_Min_Sample_Update = 200
        integer :: MPI_Sample_update_freq = 40
        logical :: MPI_Check_Limit_Converge = .false.
        !After given R_Stop is reached, can optionally check for limit convergence
        !(which is generally a much more stringent test of enough samples)
        real(mcp) :: MPI_Limit_Converge = 0.025
        real    :: MPI_Limit_Converge_Err = 0.3
        !Tolerated cross-chain error on the limit in units of standard deviation
        integer :: MPI_Limit_Param = 0
        !Which parameter's limits to check. If zero, do them all
        logical :: MPI_LearnPropose = .true.
        !If convergence R is too bad, we don't update proposal matrix if we had one to start with
        !(which should be quite good unless using new parameters or much better data)
        real :: MPI_Max_R_ProposeUpdate = 2., MPI_Max_R_ProposeUpdateNew  = 30.
        real :: MPI_R_StopProposeUpdate = 0.
        real(time_dp) :: MPI_StartTime
        real(mcp), private, allocatable, dimension(:,:) :: MPICovmat
    end type

    Type, extends(TSampleCollector) :: TMpiChainCollector
        integer :: sample_num = 0 !number without storage thinning
        real(mcp) :: StartLike = LogZero
        !bad, unless re-starting in which case it is set to last value in chain
        !zero unless from checkpoint
        logical :: Burn_done = .false.
        logical :: all_burn = .false.
        logical :: DoUpdates = .false.
        logical :: flukecheck = .false.
        integer :: indep_sample = 0
        Type(TBinaryFile) :: OutDataFile
        real(mcp) :: acc = 0, indep_acc=0 !counting numbers of samples added
        !number of iterations between dumping full model data. If zero then skip.
        integer :: checkpoint_burn = 0
        integer :: checkpoint_freq = 100
        integer, allocatable, dimension(:) :: req, buf
        integer, allocatable, dimension(:) :: param_changes
        logical :: Waiting = .false.
        logical :: done_check = .false.

        Type(TMPIData) :: Mpi
        Type(TSampleList) :: Samples
        class(TChainSampler), pointer :: Sampler
    contains
    procedure :: SaveState => TMpiChainCollector_SaveState
    procedure :: ReadState => TMpiChainCollector_ReadState
    procedure :: WriteCheckpoint => TMpiChainCollector_WriteCheckpoint
    procedure :: ReadCheckpoint => TMpiChainCollector_ReadCheckpoint
    procedure :: AddNewPoint => TMpiChainCollector_AddNewPoint
    procedure :: AddNewWeightedPoint => TMpiChainCollector_AddNewWeightedPoint
    procedure :: UpdateCovAndCheckConverge => TMpiChainCollector_UpdateCovAndCheckConverge
    procedure :: ReadParams => TMpiChainCollector_ReadParams
    FINAL :: TMpiChainCollector_Clear
#ifdef MPI
    procedure :: CheckLimitsConverge
#endif
    end Type

    integer, parameter :: chk_id = 3252359

    public TMpiChainCollector, TMPIData
    contains


    subroutine TMpiChainCollector_AddNewWeightedPoint(this, CurParams, CurLike, mult, thin_fac)
    class(TMpiChainCollector) :: this
    class(TCalculationAtParamPoint), intent(in) :: CurParams
    real(mcp) CurLike
    real(mcp), intent(in):: mult !the weight, usually integer for standard chains
    integer, intent(in), optional :: thin_fac
    integer thin
    logical want

    thin=1
    if (present(thin_fac)) thin=thin_fac
    want= this%checkpoint_burn ==0
    if (want) then
        this%acc = this%acc + mult
        this%indep_acc= this%indep_acc + mult
    end if
    if (this%acc >= thin .or. thin==1) then
        if (want) then
            output_lines = output_lines +1
            call CurParams%WriteParams(this%Config,this%acc/thin,CurLike)
        end if
        this%acc = mod(this%acc, real(thin,mcp))
    end if
    if (this%checkpoint_burn/=0) this%checkpoint_burn = this%checkpoint_burn-1
    if (this%indep_sample /= 0 .and. this%indep_acc >= this%indep_sample*thin .and. want) then
        if (.not. this%OutDataFile%Opened()) call this%OutDataFile%CreateOpenFile(rootname//'.data',append=.not. new_chains)
        call CurParams%WriteModel(this%OutDataFile, CurLike, real(floor(this%indep_acc/(this%indep_sample*thin)),mcp))
        this%indep_acc = mod(this%indep_acc, real(this%indep_sample*thin,mcp))
    end if

    end subroutine TMpiChainCollector_AddNewWeightedPoint

    subroutine TMpiChainCollector_ReadParams(this,Ini)
    class(TMpiChainCollector) :: this
    class(TSettingIni) :: Ini

#ifdef MPI
    this%Mpi%MPI_StartTime = MPI_WTime()
    call Ini%Read('MPI_Converge_Stop',this%Mpi%MPI_R_Stop)
    call Ini%Read('MPI_LearnPropose',this%Mpi%MPI_LearnPropose)
    if (this%Mpi%MPI_LearnPropose) then
        call Ini%Read('MPI_R_StopProposeUpdate',this%Mpi%MPI_R_StopProposeUpdate)
        call Ini%Read('MPI_Max_R_ProposeUpdate',this%Mpi%MPI_Max_R_ProposeUpdate)
        call Ini%Read('MPI_Max_R_ProposeUpdateNew',this%Mpi%MPI_Max_R_ProposeUpdateNew)
    end if
    call Ini%Read('MPI_Check_Limit_Converge',this%Mpi%MPI_Check_Limit_Converge)
    if (this%Mpi%MPI_Check_Limit_Converge) then
        call Ini%Read('MPI_Limit_Converge',this%Mpi%MPI_Limit_Converge)
        call Ini%Read('MPI_Limit_Converge_Err',this%Mpi%MPI_Limit_Converge_Err)
        call Ini%Read('MPI_Limit_Param',this%Mpi%MPI_Limit_Param)
    end if
#endif
    this%indep_sample = Ini%Read_Int('indep_sample')

    end subroutine TMpiChainCollector_ReadParams


    subroutine TMpiChainCollector_SaveState(this,F)
    class(TMpiChainCollector) :: this
    class(TFileStream) :: F
    integer :: version=2

    call F%Write(version)
    call F%Write(this%Mpi%MPI_thin_fac, this%Burn_done, this%all_burn,  &
        & this%flukecheck,  this%Mpi%MPI_Min_Sample_Update, this%DoUpdates)
    call F%Write(this%Mpi%MPI_Sample_update_freq)
    call this%Samples%SaveState(F)
    call this%Sampler%SaveState(F)

    end subroutine TMpiChainCollector_SaveState

    subroutine TMpiChainCollector_ReadState(this,F)
    class(TMpiChainCollector) :: this
    class(TFileStream) :: F
    integer version

    if (.not. F%ReadItem(version) .or. version>2) call MpiStop('unknown checkpoint format')
    call F%Read(this%Mpi%MPI_thin_fac, this%Burn_done, this%all_burn, &
        & this%flukecheck, this%Mpi%MPI_Min_Sample_Update, this%DoUpdates)
    if (version>=2) then
        call F%Read(this%Mpi%MPI_Sample_update_freq)
    else
        !bug fix
        this%Mpi%MPI_Sample_update_freq = this%Mpi%MPI_Sample_update_freq*BaseParams%num_slow
    end if
    call this%Samples%LoadState(F)
    this%checkpoint_burn=this%checkpoint_freq/3
    call this%Sampler%LoadState(F)

    end subroutine TMpiChainCollector_ReadState


    subroutine TMpiChainCollector_WriteCheckpoint(this)
    class(TMpiChainCollector) :: this
    Type(TBinaryFile) F

    if (Feedback > 1) write (*,*) instance, 'Writing checkpoint'
    call F%CreateFile(rootname//'.chk_tmp')
    !Use temporary file in case crash/stop during write operation
    call F%Write(chk_id)
    call this%SaveState(F)
    call F%Close()
    call File%Delete(rootname//'.chk')
    call Rename(rootname//'.chk_tmp',rootname//'.chk')

    end subroutine TMpiChainCollector_WriteCheckpoint


    subroutine TMpiChainCollector_ReadCheckpoint(this)
    class(TMpiChainCollector) :: this
    integer :: ID
    Type(TBinaryFile) F

    if (Feedback > 0) write (*,*) instance, 'Reading checkpoint from '//rootname//'.chk'
    if (LogFile%Opened()) call LogFile%Write('Re-starting from checkpoint')
    call F%Open(rootname//'.chk')
    if (.not. F%ReadItem(ID) .or. ID/=chk_id) call DoAbort('invalid checkpoint files')
    call this%ReadState(F)
    call F%Close()

    end subroutine TMpiChainCollector_ReadCheckpoint

    subroutine TMpiChainCollector_Clear(this)
    Type(TMpiChainCollector) :: this

    call this%OutDataFile%Close()

    end subroutine TMpiChainCollector_Clear


    subroutine TMpiChainCollector_UpdateCovAndCheckConverge(this)
    class(TMpiChainCollector) :: this
#ifdef MPI
    integer i,j
    real(mcp), allocatable, dimension(:,:,:) ::MPICovmats
    real(mcp), allocatable, dimension(:,:)   ::MPIMeans, MpiCovmat
    real(mcp), allocatable, dimension(:)     ::MPIMean
    real(mcp) delta(num_params_used)
    integer ierror
    real(mcp) norm, mean(num_params_used), chain_means(MPIChains,num_params_used)
    real(mcp) MeansCov(num_params_used,num_params_used), cov(num_params_used,num_params_used)
    real(mcp) evals(num_params_used), R
    logical :: invertible
    character(LEN=128) logLine

    if (Feedback > 0) write (*,*) 'Chain',instance,' MPI communicating'
    allocate(MPIcovmats(num_params_used,num_params_used,MPIChains))
    allocate(MPIMeans(0:num_params_used,MPIChains))
    allocate(MPIMean(0:num_params_used))
    allocate(MPIcovmat(num_params_used,num_params_used))

    MPICovMat = 0
    MPIMean = 0
    MPImean(0) = this%Samples%Count - this%Samples%Count/2 + 1
    do i = this%Samples%Count/2, this%Samples%Count
        MPiMean(1:num_params_used) = MPiMean(1:num_params_used) + this%Samples%Item(i)
    end do
    MPiMean(1:num_params_used) = MPiMean(1:num_params_used) / MPImean(0)
    do i = this%Samples%Count/2, this%Samples%Count
        delta = this%Samples%Item(i)-MPIMean(1:num_params_used)
        do j = 1, num_params_used
            MPICovmat(:,j) =  MPICovmat(:,j) + delta*(this%Samples%Item(i,j)- MPIMean(j))
        end do
    end do
    MPICovMat = MPICovMat / MPImean(0)

    call MPI_ALLGATHER(MPICovMat,Size(MPICovMat),MPI_real_mcp,MPICovmats,Size(MPICovmat), &
        MPI_real_mcp, MPI_COMM_WORLD,ierror)
    call MPI_ALLGATHER(MPIMean,Size(MPIMean),MPI_real_mcp,MPIMeans,Size(MPIMean), &
        MPI_real_mcp, MPI_COMM_WORLD,ierror)

    if (all(MPIMeans(0,:)> this%Mpi%MPI_Min_Sample_Update/2 + 2)) then
        !check have reasonable number of samples in each)
        norm = sum(MPIMeans(0,:))

        do i=1, num_params_used
            mean(i) = sum(MPIMeans(i,:)*MPIMeans(0,:))/norm
            chain_means(:,i) = MPIMeans(i,:)
        end do

        if (MPIChains > 1) then
            do i=1,num_params_used
                do j=i,num_params_used
                    MPICovMat(i,j) = sum(MPIMeans(0,:)*MPICovMats(i,j,:))/ norm
                    cov(i,j) = sum(MPICovMats(i,j,:))/ MPIChains
                    meanscov(i,j) = sum(MPIMeans(0,:)*&
                        (chain_means(:,i)-mean(i))*(chain_means(:,j)-mean(j)))/norm
                    meanscov(j,i) = meanscov(i,j)
                    cov(j,i) = cov(i,j)
                    MPICovMat(j,i) = MPICovMat(i,j)
                end do
            end do
            meansCov = meansCov * real(MPIChains,mcp)/(MPIChains-1)
            invertible= GelmanRubinEvalues(cov, meanscov, evals, num_params_used)
            if (invertible) then
                R = maxval(evals)
                if (Feedback > 1 .and. MPIRank==0) write (*,*) 'Convergence e-values: ', real(evals)
                if (Feedback > 0 .and. MPIRank==0) then
                    write (*,*) 'Current convergence R-1 = ',real(R), ' chain steps =',this%sample_num
                    call this%Sampler%LikeCalculator%WritePerformanceStats(stdout)
                end if
                if (LogFile%Opened()) then
                    write(LogFile%unit,*) 'Current convergence R-1 = ',real(R), ' chain steps =',this%sample_num
                    if (flush_write) call LogFile%Flush()
                    call this%Sampler%LikeCalculator%WritePerformanceStats(LogFile%unit)
                end if
                if (R < this%Mpi%MPI_R_Stop .and. this%flukecheck) then
                    if (this%Mpi%MPI_Check_Limit_Converge) then
                        !Now check if limits from different chains agree well enough
                        if (this%CheckLImitsConverge(this%Samples, MpiCovmat)) &
                            call ConvergeStatus(.true., R, 'Requested limit convergence achieved')
                    else
                        !If not also checking limits, we are done
                        call ConvergeStatus(.true., R, 'Requested convergence R achieved')
                    end if
                end if
                call ConvergeStatus(.false., R)
                this%flukecheck = R < this%Mpi%MPI_R_Stop
                if (this%Samples%Count > 500000) then
                    !Try not to blow memory by storing too many samples
                    call this%Samples%Thin(2)
                    this%Mpi%MPI_thin_fac = this%Mpi%MPI_thin_fac*2
                end if
            else
                if (Feedback > 0 .and. MPIRank==0) write(*,*) 'Covariance not currently invertible'
                R=1e6
            end if
        end if !MPIChains >1

        if (this%Mpi%MPI_LearnPropose .and. ( MPIChains==1 .or. (BaseParams%covariance_is_diagonal &
            .or. R < this%Mpi%MPI_Max_R_ProposeUpdate) .and. R > this%Mpi%MPI_R_StopProposeUpdate)) then
        !If beginning to converge, update covariance matrix
        if (Feedback > 0 .and. MPIRank==0) write (*,*) 'updating proposal density'

        call this%Sampler%SetCovariance(MPICovMat)

        end if !update propose
    end if !all chains have enough

#endif
    end subroutine TMpiChainCollector_UpdateCovAndCheckConverge

    subroutine TMpiChainCollector_AddNewPoint(this,P,like)
    class(TMpiChainCollector) :: this
    real(mcp), intent(in) ::like
    real(mcp) P(:)
    !Collect thinned samples after a burn-in perdiod
    !Then use second half of the samples to get convergence
    !Use R = worst eigenvalue (variance of chain means)/(mean of chain variances) statistic for
    !convergence test, followed optionally by (variance of limit)/(mean variance) statistic
    !If MPI_LearnPropose then updates proposal density using covariance matrix of last half of chains
#ifdef MPI
    integer i,j
    integer STATUS(MPI_STATUS_SIZE),STATs(MPI_STATUS_SIZE*(MPIChains-1))
    logical flag, ierror

    !Dump checkpoint info
    !Have to be careful if were to dump before burn
    if (checkpoint .and. this%all_burn .and. this%checkpoint_burn==0 .and. &
        (.not. this%done_check .or.  mod(this%sample_num+1, this%checkpoint_freq)==0)) then
    this%done_check=.true.
    call this%WriteCheckpoint()
    end if

    !Do main adding samples functions
    this%sample_num = this%sample_num + 1
    if (mod(this%sample_num, this%Mpi%MPI_thin_fac) /= 0) return

    if (.not. allocated(this%req)) allocate(this%req(MPIChains-1))
    call this%Samples%Add(P(params_used))

    if (.not. this%Burn_done) then
        if (this%Samples%Count > 50 +1) then
            !We're not really after independent samples or all of burn in
            !Make sure all parameters are being explored
            if (.not. allocated(this%param_changes)) then
                allocate(this%param_changes(num_params_used))
                this%param_changes= 0
            end if
            do i=1, num_params_used
                if (this%Samples%Value(this%Samples%Count, i) /= this%Samples%Value(this%Samples%Count-1, i)) this%param_changes(i) =  this%param_changes(i) + 1
            end do
            this%Burn_done = all(this%param_changes > 50 +1)
            if (this%Burn_done) then
                if (Feedback > 0) then
                    write (*,*) trim(concat('Chain',instance, ', MPI done ''burn'', Samples = ',this%sample_num))//', like = ',real(like)
                    write (*,*) 'Time: ', MPI_WTime() - this%Mpi%MPI_StartTime, 'output lines=',output_lines
                    call this%Sampler%LikeCalculator%WritePerformanceStats(stdout)
                end if

                !Here we make something like an MPE_IBARRIER to see if all threads have passed burn in
                !On completion of IRECV all should be OK

                allocate(this%buf(MPIChains-1))

                i = 0
                do j=0, MPIChains-1
                    if (j /= MPIRank) then
                        i=i+1
                        call MPI_ISEND(MPIRank,1,MPI_INTEGER, j,0,MPI_COMM_WORLD,this%req(i),ierror)
                        call MPI_IRECV(this%buf(i),1,MPI_INTEGER, j,0,MPI_COMM_WORLD,this%req(i),ierror)
                    end if
                end do

                this%Mpi%MPI_Min_Sample_Update =  50 + BaseParams%num_slow*4 + BaseParams%num_fast
                if (sampling_method == sampling_fast_dragging) then
                    this%Mpi%MPI_Min_Sample_Update=this%Mpi%MPI_Min_Sample_Update*this%Sampler%oversample_fast
                else
                    this%Mpi%MPI_Min_Sample_Update=this%Mpi%MPI_Min_Sample_Update  + BaseParams%num_fast*4
                end if
                print *,'MPI_Min_Sample_Update',this%Mpi%MPI_Min_Sample_Update,this%Samples%Count
                if (this%Samples%Count>this%Mpi%MPI_Min_Sample_Update) call this%Samples%DeleteRange(1, this%Samples%Count-this%Mpi%MPI_Min_Sample_Update)
                if (sampling_method/= sampling_fast_dragging) then
                    this%Mpi%MPI_Sample_update_freq=this%Mpi%MPI_Sample_update_freq*num_params_used
                else
                    this%Mpi%MPI_Sample_update_freq=this%Mpi%MPI_Sample_update_freq*BaseParams%num_slow
                end if
                this%flukecheck = .false.
                deallocate(this%param_changes)
            end if
        end if

    else
        flag = .false.

        if (.not. this%all_burn) then
            call MPI_TESTALL(MPIChains-1,this%req, this%all_burn, stats, ierror)
            if (this%all_burn) then
                deallocate(this%buf)
                if (Feedback>0) write(*,*) instance, 'all_burn done'
            end if
        end if


        if (.not. this%DoUpdates  .and. this%all_burn .and. this%Samples%Count >= this%MPi%MPI_Min_Sample_Update+1) then
            this%DoUpdates = .true.
            if (Feedback>0) write(*,*) instance, 'DoUpdates'
        end if

        if (this%DoUpdates) then
            if (MPIRank == 0) then
                if (this%Waiting) then
                    call MPI_TESTALL(MPIChains-1,this%req, flag, stats, ierror)
                    this%Waiting = .not. flag
                elseif (mod(this%Samples%Count,max(1,this%Mpi%MPI_Sample_update_freq))==0) then
                    this%Waiting = .true.
                    do j=1, MPIChains-1
                        call MPI_ISSEND(MPIRank,1,MPI_INTEGER, j,0,MPI_COMM_WORLD,this%req(j),ierror)
                    end do
                end if

            else
                !See if notified by root chain that time to do stuff
                call MPI_IPROBE(0,0,MPI_COMM_WORLD,flag, status,ierror)
                if (flag)  then
                    call MPI_RECV(i,1,MPI_INTEGER, 0,0,MPI_COMM_WORLD,status,ierror)
                    !Just get rid of it. Must be neater way to do this...
                end if
            end if

        end if

        if (flag) then
            !update covariances, check for convergence
            call this%UpdateCovAndCheckConverge()
        end if !flag
    end if
#endif

    end subroutine TMpiChainCollector_AddNewPoint


#ifdef MPI

    subroutine ConvergeStatus(isDone, R, Msg)
    logical, intent(in) :: isDone
    character(LEN=*), intent(in), optional :: Msg
    real(mcp), intent(in) :: R
    Type(TTextFile) :: F

    if (MPiRank==0) then
        call F%CreateFile(trim(baseroot)//'.converge_stat')
        call F%Write(R)
        if (isDone) call F%Write('Done')
        call F%Close()
    end if
    if (isDone) call DoStop(Msg)

    end subroutine COnvergeStatus

    function CheckLimitsConverge(this, L, MpiCovmat)
    !Check limits from last half chains agree well enough across chains to be confident of result
    !Slowly explored tails will cause problems (long time till stops)
    class(TMpiChainCollector) :: this
    class(TSampleList), intent(in) :: L
    real(mcp), intent(in) :: MpiCovmat(:,:)
    integer i,j, side, ierror, worsti
    real(mcp), allocatable, dimension(:,:,:) :: Limits
    logical :: CheckLimitsConverge
    integer numCheck
    integer, allocatable, dimension(:) :: params_check
    character(LEN=128) logLine
    real(mcp) MeanLimit, var, LimErr, WorstErr

    if (this%Mpi%MPI_Limit_Param/=0) then
        numCheck = 1
        allocate(params_check(numCheck))
        do j=1,num_params_used
            if (params_used(j) == this%Mpi%MPI_Limit_Param) then
                params_check(1) = j
                exit
            end if
        end do
    else
        numCheck = num_params_used
        allocate(params_check(numCheck))
        params_check = (/ (I, I=1, num_params_used) /)
    end if

    allocate(Limits(2,numCheck,MPIChains))

    do j=1, numCheck
        call L%ConfidVal(params_check(j), this%Mpi%MPI_Limit_Converge, L%Count/2, L%Count, &
            Limits(1,j,instance),Limits(2,j,instance))
    end do
    !Now tell everyone else
    do i=1, MPIChains
        j = i-1
        call MPI_BCAST(Limits(:,:,i),2*numCheck,MPI_real_mcp,j,MPI_COMM_WORLD,ierror)
    end do
    !Take as test statistics the rms deviation from the mean limit in units of the standard deviation
    WorstErr = 0
    do j=1, numCheck
        do side = 1,2
            MeanLimit = Sum(Limits(side,j,:))/MPIChains
            var = sum((Limits(side,j,:) - MeanLimit)**2)/(MPIChains-1)
            LimErr = sqrt(var / MPICovMat(params_check(j),params_check(j)))
            if(LimErr > WorstErr) then
                WorstErr = LimErr
                Worsti = params_check(j)
            end if
            if (Feedback > 0 .and. MPIRank ==0) &
                write(*,'(1A22,1A8,1f8.3)') BaseParams%UsedParamNameOrNumber(params_check(j)),'lim err', LimErr
        end do
    end do
    if (Feedback > 0 .and. MPIRank==0) then
        write (*,'(a)')  'Current limit err = '//trim(RealToStr(WorstErr))// &
            ' for '//trim(BaseParams%UsedParamNameOrNumber(Worsti))//'; samps = '//trim(IntToStr(L%Count*this%Mpi%MPI_thin_fac))
    end if
    if (LogFile%Opened()) then
        write (logLine,'(a)') 'Current limit err = '//trim(RealToStr(WorstErr))// &
            ' for '//trim(BaseParams%UsedParamNameOrNumber(Worsti))//'; samps = '//trim(IntToStr(L%Count*this%Mpi%MPI_thin_fac))
        call LogFile%Write(logLine)
        if (flush_write) call LogFile%Flush()
    end if
    CheckLimitsConverge = (WorstErr < this%Mpi%MPI_Limit_Converge_Err)

    end function CheckLimitsConverge

#endif



    end module SampleCollector
