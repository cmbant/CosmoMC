    !Defines a parameterization, computes likelihoods and defined proposal density
    !Change this file to change these things, and MPI updating
    !The Cls are computed using routines defined in the CMB_Cls module

    module ParamDef
    !module defines a parameterization, and works out power spectra
    use CMB_Cls
    use cmbtypes
    use cmbdata
    use Random
    use settings
    use Propose
    use Samples
    use likelihood
    use IniFile
    use ParamNames
    use IO

    implicit none

    Type :: ParamSet
        real(mcp) :: P(max_num_params)
        real(mcp) likelihoods(max_likelihood_functions)
        Type(TheoryPredictions) :: Theory
        Type(ParamSetInfo) Info
    contains
    procedure :: ReadModel
    procedure :: WriteModel
    end Type

    Type ParamScale
        real(mcp) PMin(max_num_params), PMax(max_num_params), StartWidth(max_num_params), &
        PWidth(max_num_params), center(max_num_params)
    end Type ParamScale

    Type(ParamScale) Scales

    Type ParamGaussPrior
        real(mcp) mean(max_num_params),std(max_num_params) !std=0 for no prior (default)
    end Type ParamGaussPrior
    Type(ParamGaussPrior) GaussPriors

    real(mcp) :: StartLike = LogZero
    !bad, unless re-starting in which case it is set to last value in chain
    integer :: Num = 0, num_accept=0
    !zero unless from checkpoint
    integer :: burn_in = 2 !Minimum to get sensible answers
    integer :: semislow_changes=0, slow_changes=0
    integer :: checkpoint_burn = 0

    logical :: estimate_propose_matrix = .false.
    real(mcp), dimension(:,:), allocatable ::  initial_propose_matrix, test_cov_matrix
    Type(BlockedProposer), save :: Proposer

    logical :: checkpoint = .false.
    integer, parameter :: checkpoint_freq = 100
    character(LEN=500) :: rootname

    real  :: MPI_R_Stop = 0.05

    integer  :: MPI_thin_fac = 1
    !following are numbers of samples / (number of parameters)
    !MPI_Min_Sample_Update*num_fast is min number of chain steps
    !before starting to update covmat and checking convergence
    integer :: MPI_Min_Sample_Update = 200, MPI_Sample_update_freq = 40

    logical :: MPI_Check_Limit_Converge = .false.
    !After given R_Stop is reached, can optionally check for limit convergence
    !(which is generally a much more stringent test of enough samples)
    real(mcp)    :: MPI_Limit_Converge = 0.025
    real    :: MPI_Limit_Converge_Err = 0.3
    !Tolerated cross-chain error on the limit in units of standard deviation
    integer :: MPI_Limit_Param = 0
    !Which parameter's limits to check. If zero, do them all
    logical :: MPI_LearnPropose = .true.

    !If convergence R is too bad, we don't update proposal matrix if we had one to start with
    !(which should be quite good unless using new parameters or much better data)
    real  :: MPI_Max_R_ProposeUpdate = 2., MPI_Max_R_ProposeUpdateNew  = 30.

    real    :: MPI_R_StopProposeUpdate = 0.

    logical    :: MPI_StartSliceSampling = .false.

    real(time_dp) ::  MPI_StartTime

    real(mcp), private, allocatable, dimension(:,:) :: MPICovmat
    logical :: StartCovMat = .false.

    contains

    subroutine DoStop(S, abort)
    character(LEN=*), intent(in), optional :: S
    integer ierror
    logical, intent(in), optional :: abort
    logical wantbort

    if (outfile_handle/=0) call IO_Close(outfile_handle)

    if (present(abort)) then
        wantbort = abort
    else
        wantbort = .false.
    end if

    if (present(S) .and. (wantbort .or. MPIRank==0)) write (*,*) trim(S)
#ifdef MPI
    MPI_StartTime = MPI_WTime() - MPI_StartTime
    if (Feedback > 0 .and. MPIRank==0) then

    write (*,*) 'Total time:', nint(MPI_StartTime), &
    '(',MPI_StartTime/(60*60),' hours)'

    if (slow_proposals/=0) write (*,*) 'Slow proposals: ', slow_proposals

    end if
    ierror =0
    if (wantbort) then
        !Abort all in case other continuing chains want to communicate with us
        !in the case when max number of samples is reached
        call MPI_Abort(MPI_COMM_WORLD,ierror,ierror)
    else
        call mpi_finalize(ierror)
    end if
#endif

#ifdef DECONLY
    pause
#endif
    stop
    end subroutine DoStop


    subroutine ParamError(str,param)
    character(LEN=*), intent(in) :: str
    integer, intent(in) :: param

    call DoAbort(trim(str)//': '//trim(ParamNames_NameOrNumber(NameMapping, param)))

    end subroutine ParamError

    subroutine orderIndices(arr,n)
    integer, intent(in) :: n
    integer arr(:), tmp(n),j, s, i

    tmp=arr(1:n)
    s=n
    do i=1,n
        arr(i) = minval(tmp(1:s))
        j=indexOf(arr(i),tmp,s)
        tmp(j:s-1)= tmp(j+1:s)
        s=s-1
    end do
    end subroutine orderIndices

    subroutine SetStartPositions(Params)
    type (ParamSet) Params
    integer ix,i

    do ix=1,num_params_used
        i = params_used(ix)
        do
            if (Scales%StartWidth(i) < 0) then
                !This case we want half gaussian, width -Scales%StartWidth(i)
                !e.g. for positive definite parameters
                Params%P(i) = Scales%center(i) - abs(Gaussian1())*Scales%StartWidth(i)
            else
                Params%P(i) = Scales%center(i) + Gaussian1()*Scales%StartWidth(i)
            end if
            !Repeat until get acceptable values in range
            if (Params%P(i)>=  Scales%PMin(i) .and. Params%P(i) <= Scales%PMax(i)) exit
        end do
    end do

    end subroutine SetStartPositions

    subroutine InitializeUsedParams(Ini,Params, use_propose_matrix)
    type (ParamSet) Params
    type (TIniFile) :: Ini
    logical :: use_propose_matrix

    num_params = max(num_theory_params, NameMapping%num_MCMC)
    num_data_params = num_params - num_theory_params
    call Initalize(Ini,Params,use_propose_matrix)

    end subroutine InitializeUsedParams

    subroutine Initalize(Ini,Params,use_propose_matrix)
    implicit none
    type (ParamSet) Params
    type (TIniFile) :: Ini
    logical :: use_propose_matrix
    integer i, j, ix
    character(LEN=5000) fname,InLine
    character(LEN=1024) prop_mat
    real(mcp) center
    integer fast_params(num_params)
    integer fast_number, fast_param_index
    integer param_type(num_params)
    integer speed, num_speed
    integer, parameter :: tp_unused = 0, tp_slow=1, tp_semislow=2, tp_semifast=3, tp_fast=4
    Type(int_arr_pointer), allocatable :: param_blocks(:)
    logical :: block_semi_fast =.false., block_fast_likelihood_params=.false.
    integer :: breaks(num_params), num_breaks
    Type(DataLikelihood), pointer :: DataLike
    integer :: oversample_fast= 1
    logical first, ierr

    output_lines = 0

    if (use_fast_slow) then
        if (Ini_HasKey_File(Ini,'fast_parameters')) then
            !List of parmeter names to treat as fast
            fast_number = num_params
            call ParamNames_ReadIndices(NameMapping,Ini_Read_String_File(Ini,'fast_parameters'), fast_params, fast_number)
        else
            fast_param_index = max(index_data,DataLikelihoods%first_fast_param)
            !all parameters at and above fast_param_index
            fast_param_index= Ini_Read_Int_File(Ini,'fast_param_index',fast_param_index)
            fast_number = 0
            do i=fast_param_index, num_params
                fast_number = fast_number+1
                fast_params(fast_number) = i
            end do
        end if
        block_semi_fast = Ini_Read_Logical_File(Ini,'block_semi_fast',.true.)
        block_fast_likelihood_params = Ini_Read_logical_File(Ini,'block_fast_likelihood_params',.true.)
        oversample_fast = Ini_Read_Int('oversample_fast',1)
        if (oversample_fast<1) call DoAbort('oversample_fast must be >= 1')
        if (sampling_method /=sampling_fast_dragging) then
            MPI_Thin_fac = MPI_Thin_fac*Proposer%Oversample_fast
        end if
    else
        fast_number = 0
    end if

    Ini_fail_on_not_found = .false.
    if (use_propose_matrix) then
        prop_mat = trim(Ini_Read_String_File(Ini,'propose_matrix'))
        if (prop_mat /= '' .and. prop_mat(1:1) /= '/') prop_mat = concat(LocalDir,prop_mat)
    else
        prop_mat=''
    end if

    fname = Ini_Read_String_File(Ini,'continue_from')
    if (fname /= '') call DoAbort('continue_from replaced by checkpoint')

    param_type = tp_unused

    GaussPriors%std=0 !no priors by default
    do i=1,num_params
        InLine =  ParamNames_ReadIniForParam(NameMapping,Ini,'param',i)
        if (InLine=='') call ParamError('parameter ranges not found',i)
        if (TxtNumberColumns(InLine)==1) then
            !One number means just fix the parameter
            read(InLine, *, IOSTAT = ierr) center
            if (ierr/=0) call ParamError('fixed parameter value not valid',i)
            Scales%PMin(i)= center
            Scales%PMax(i)= center
            Scales%StartWidth(i)=0
            Scales%PWidth(i)=0
        else
            read(InLine, *, IOSTAT = ierr) center, Scales%PMin(i), Scales%PMax(i), Scales%StartWidth(i), Scales%PWidth(i)
            if (ierr/=0) call ParamError('Error reading param details',i) 
        end if
        if (Scales%PWidth(i)/=0) then
            InLine =  ParamNames_ReadIniForParam(NameMapping,Ini,'prior',i)
            if (InLine/='') read(InLine, *, IOSTAT = ierr) GaussPriors%mean(i), GaussPriors%std(i)
            if (ierr/=0) call DoAbort('Error reading prior mean and stdd dev: '//trim(InLIne))
        else
            scales%PMin(i) = center
            Scales%PMax(i) = center
        end if

        Scales%center(i) = center
        Params%P(i) = center
        if (Scales%PMax(i) < Scales%PMin(i)) call ParamError('You have param Max < Min',i)
        if (Scales%center(i) < Scales%PMin(i)) call ParamError('You have param center < Min',i)
        if (Scales%center(i) > Scales%PMax(i)) call ParamError('You have param center > Max',i)
        if (Scales%PWidth(i) /= 0) then !to get sizes for allocation arrays
            if (use_fast_slow .and. any(i==fast_params(1:fast_number))) then
                if (i >= index_data .or. .not. block_semi_fast) then
                    param_type(i) = tp_fast
                else
                    param_type(i) = tp_semifast
                end if
            else
                if (use_fast_slow .and. index_semislow >=0 .and. i >= index_semislow .and. block_semi_fast) then
                    param_type(i) = tp_semislow
                else
                    param_type(i) = tp_slow
                end if
            end if
        end if
    end do

    allocate(params_used( count(param_type /= tp_unused)))
    num_params_used=0
    do i=1,num_params
        if (param_type(i)/=tp_unused) then
            num_params_used=num_params_used+1
            params_used(num_params_used)=i
        end if
    end do

    num_breaks=0
    if (block_fast_likelihood_params) then
        !put parameters for different likelihoods in separate blocks,
        !so not randomly mix them and hence don't all need to be recomputed
        first=.true.
        do i=1,DataLikelihoods%Count
            DataLike=>DataLikelihoods%Item(i)
            do j=1, num_params_used-1
                if (param_type(params_used(j))==tp_fast .and. params_used(j) >= DataLike%new_param_block_start &
                .and. params_used(j) < DataLike%new_param_block_start + DataLike%new_params) then
                    if (first) then
                        first = .false.
                    else
                        num_breaks = num_breaks+1
                        breaks(num_breaks)=j
                    end if
                    exit
                end if
            end do
        end do
    end if
    num_breaks = num_breaks + 1
    breaks(num_breaks) = num_params_used
    if (Feedback>0 .and. MpiRank==0) then
        write(*,*) 'Fast divided into ',num_breaks,' blocks'
        if (num_breaks>1) write(*,*) 'Block breaks at: ',breaks(1:num_breaks-1)
    end if

    call orderIndices(breaks, num_breaks)

    allocate(param_blocks(tp_semifast+num_breaks))
    do speed= tp_slow, tp_semifast
        allocate(param_blocks(speed)%P(count(param_type == speed)))
        num_speed=0
        do i=1,num_params_used
            if (param_type(params_used(i))==speed) then
                num_speed=num_speed+1
                param_blocks(speed)%P(num_speed) = i
            end if
        end do
    end do
    ix=1
    do speed = tp_fast, tp_fast+num_breaks-1
        j = breaks(speed-tp_fast+1)
        allocate(param_blocks(speed)%P(count(param_type(params_used(ix:j)) == tp_fast)))
        num_speed=0
        do i=ix,j
            if (param_type(params_used(i))==tp_fast) then
                num_speed=num_speed+1
                param_blocks(speed)%P(num_speed) = i
            end if
        end do
        ix=j+1
    end do

    call Proposer%Init(param_blocks, slow_block_max= 2, oversample_fast=oversample_fast)
    num_slow = Proposer%Slow%n
    num_fast = Proposer%Fast%n

    if (Feedback > 0 .and. MpiRank==0) then
        write(*,'(" Varying ",1I2," parameters (",1I2," slow (",1I2," semi-slow), ",1I2," fast (",1I2," semi-fast))")') &
        num_params_used,num_slow, size(param_blocks(tp_semislow)%P), num_fast,size(param_blocks(tp_semifast)%P)
    end if

    allocate(initial_propose_matrix(num_params_used, num_params_used))
    StartCovMat  = prop_mat/=''
    call ReadSetCovMatrix(prop_mat, initial_propose_matrix)
    call Proposer%SetCovariance(initial_propose_matrix)

    test_likelihood = Ini_read_Logical_file(Ini,'test_likelihood', .false.)
    if (test_likelihood) then
        print *,'** Using test Gaussian likelihood from covariance + hard priors **'
        prop_mat = trim(Ini_Read_String_File(Ini,'test_covariance'))
        allocate(test_cov_matrix(num_params_used, num_params_used))
        if (prop_mat=='') then
            test_cov_matrix = initial_propose_matrix
        else
            call ReadSetCovMatrix(prop_mat, test_cov_matrix)
        end if
    end if

    end subroutine Initalize

    subroutine ReadSetCovMatrix(prop_mat, matrix)
    character(LEN=*), intent(in) :: prop_mat
    real(mcp) :: matrix(num_params_used, num_params_used)
    real(mcp) pmat(num_params,num_params)
    integer i

    if (prop_mat/='') then
        call IO_ReadProposeMatrix(pmat, prop_mat)

        !If generated with constrained parameters, assume diagonal in those parameters
        do i=1,num_params
            if (pmat(i,i) ==0 .and. Scales%PWidth(i)/=0) then
                pmat(i,i) = Scales%PWidth(i)**2
                MPI_Max_R_ProposeUpdate = MPI_Max_R_ProposeUpdateNew
            end if
            !Enforce new constraints (should really be fixing the inverse...)
            if (Scales%PWidth(i)==0) then
                pmat(i,:) = 0
                pmat(:,i) = 0
            end if
        end do

        matrix = pmat(params_used, params_used)
    else
        matrix = 0
        do i=1,num_params_used
            matrix(i,i) = Scales%PWidth(params_used(i))**2
        end do
    end if
    end subroutine ReadSetCovMatrix

    subroutine WriteIndepSample(P, like, mult)
    Type(ParamSet) P
    real(mcp) like, mult
    if (indepfile_handle ==0) return
    call P%WriteModel(indepfile_handle, like, mult)
    end subroutine WriteIndepSample


    subroutine AddMPIParams(P,like, checkpoint_start)
    real(mcp), intent(in) ::like
    real(mcp) P(:)
    logical, intent(in), optional :: checkpoint_start
    !Collect thinned samples after a burn-in perdiod
    !Then use second half of the samples to get convergence
    !Use R = worst eigenvalue (variance of chain means)/(mean of chain variances) statistic for
    !convergence test, followed optionally by (variance of limit)/(mean variance) statistic
    !If MPI_LearnPropose then updates proposal density using covariance matrix of last half of chains
#ifdef MPI

    integer, save :: sample_num = 0
    integer i,j,k
    logical, save :: Burn_done = .false.
    integer index, STATUS(MPI_STATUS_SIZE),STATs(MPI_STATUS_SIZE*(MPIChains-1))
    logical flag

    real(mcp), allocatable, dimension(:), save ::MPIMean
    real(mcp), allocatable, dimension(:,:,:) ::MPICovmats
    real(mcp), allocatable, dimension(:,:)   ::MPIMeans
    real(mcp) norm, mean(num_params_used), chain_means(MPIChains,num_params_used)
    real(mcp) delta(num_params_used)
    integer, parameter :: chk_id = 3252358
    integer ID
    real(mcp) MeansCov(num_params_used,num_params_used), cov(num_params_used,num_params_used)
    real(mcp) evals(num_params_used), last_P, R
    integer ierror, dummy
    logical DoCheckpoint, invertible

    integer, allocatable, dimension(:), save :: req, buf
    integer,  allocatable, dimension(:), save :: param_changes
    logical, save :: flukecheck, Waiting = .false.

    Type(TSampleList), save :: S
    integer, save :: slice_fac = 1
    logical, save :: all_burn = .false., done_check = .false., DoUpdates = .false.
    character(LEN=128) logLine


    !Read in checkpoing stuff at restart
    if (checkpoint .and. present(checkpoint_start)) then
        if (Feedback > 0) write (*,*) instance, 'Reading checkpoint from '//trim(rootname)//'.chk'
        call OpenFile(trim(rootname)//'.chk',tmp_file_unit,'unformatted')
        read (tmp_file_unit) ID
        if (ID/=chk_id) call DoAbort('invalid checkpoint files')
        read(tmp_file_unit) num, sample_num, num_accept, MPI_thin_fac, dummy, Burn_done, all_burn,sampling_method, &
        slice_fac, S%Count, flukecheck, StartCovMat, MPI_Min_Sample_Update, DoUpdates
        read(tmp_file_unit) cov
        call Proposer%SetCovariance(cov)
        call S%ReadBinary(tmp_file_unit)
        close(tmp_file_unit)
        allocate(req(MPIChains-1))
        allocate(MPIcovmat(num_params_used,num_params_used))
        allocate(MPIMean(0:num_params_used))
        checkpoint_burn=checkpoint_freq/3
        return
    end if

    !Dump checkpoint info
    !Have to be careful if were to dump before burn
    if (checkpoint .and. all_burn .and. checkpoint_burn==0 .and. &
    (.not. done_check .or.  mod(sample_num+1, checkpoint_freq*Proposer%Oversample_fast)==0)) then
        done_check=.true.
        if (Feedback > 1) write (*,*) instance, 'Writing checkpoint'
        call CreateFile(trim(rootname)//'.chk_tmp',tmp_file_unit,'unformatted')
        !Use temporary file in case crash/stop during write operation
        write (tmp_file_unit) chk_id
        write(tmp_file_unit) num, sample_num, num_accept, MPI_thin_fac, S%Count, Burn_done, all_burn, sampling_method, &
        slice_fac, S%Count, flukecheck, StartCovMat, MPI_Min_Sample_Update, DoUpdates
        write(tmp_file_unit) Proposer%propose_matrix
        call S%SaveBinary(tmp_file_unit)
        close(tmp_file_unit)
        call DeleteFile(trim(rootname)//'.chk')
        call Rename(trim(rootname)//'.chk_tmp',trim(rootname)//'.chk')
    end if

    !Do main adding samples functions
    sample_num = sample_num + 1
    if (mod(sample_num, MPI_thin_fac) /= 0) return

    if (S%Count == 0 .and. .not. Burn_done)then
        allocate(param_changes(num_params_used))
        param_changes= 0
        if (MPI_StartSliceSampling) then
            sampling_method = sampling_slice
        end if

        if (sampling_method == sampling_slice) then
            if (Feedback > 0  .and. .not. Burn_done) write (*,*) 'Starting with slice sampling'
            slice_fac = 4
        end if
    end if

    call S%Add(P(params_used))

    if (.not. Burn_done) then
        if (S%Count > 50 +1) then
            !We're not really after independent samples or all of burn in
            !Make sure all parameters are being explored
            do i=1, num_params_used
                if (S%Value(S%Count, i) /= S%Value(S%Count-1, i)) param_changes(i) =  param_changes(i) + 1
            end do
            Burn_done = all(param_changes > 50/slice_fac+2)
            if (Burn_done) then
                if (Feedback > 0) then
                    write (*,*) trim(concat('Chain',instance, ', MPI done ''burn'', Samples = ',sample_num))//', like = ',real(like)
                    write (*,*) 'Time: ', MPI_WTime() - MPI_StartTime, 'output lines=',output_lines
                    if (use_fast_slow) write(*,*) 'slow changes', slow_changes, 'semi-slow changes', semislow_changes
                end if

                !Here we make something like an MPE_IBARRIER to see if all threads have passed burn in
                !On completion of IRECV all should be OK

                allocate(req(MPIChains-1), buf(MPIChains-1))

                i = 0
                do j=0, MPIChains-1
                    if (j /= MPIRank) then
                        i=i+1
                        call MPI_ISEND(MPIRank,1,MPI_INTEGER, j,0,MPI_COMM_WORLD,req(i),ierror)
                        call MPI_IRECV(buf(i),1,MPI_INTEGER, j,0,MPI_COMM_WORLD,req(i),ierror)
                    end if
                end do

                MPI_Min_Sample_Update =  50 + num_slow*4 + num_fast
                if (sampling_method == sampling_fast_dragging) then
                    MPI_Min_Sample_Update=MPI_Min_Sample_Update*Proposer%oversample_fast
                else
                    MPI_Min_Sample_Update=MPI_Min_Sample_Update  + num_fast*4
                end if
                print *,'MPI_Min_Sample_Update',MPI_Min_Sample_Update,S%Count
                if (S%Count>MPI_Min_Sample_Update) call S%DeleteRange(1, S%Count-MPI_Min_Sample_Update)
                if (sampling_method/= sampling_fast_dragging) then
                    MPI_Sample_update_freq=MPI_Sample_update_freq*num_params_used
                else
                    MPI_Sample_update_freq=MPI_Sample_update_freq*num_slow*Proposer%oversample_fast
                end if
                flukecheck = .false.
                deallocate(param_changes)
                allocate(MPIcovmat(num_params_used,num_params_used))
                allocate(MPIMean(0:num_params_used))
            end if
        end if

    else
        flag = .false.

        if (.not. all_burn) then
            call MPI_TESTALL(MPIChains-1,req, all_burn, stats, ierror)
            if (all_burn) then
                deallocate(buf)
                if (Feedback>0) write(*,*) instance, 'all_burn done'
            end if
        end if


        if (.not. DoUpdates  .and. all_burn .and. S%Count >= MPI_Min_Sample_Update+1) then
            DoUpdates = .true.
            if (Feedback>0) write(*,*) instance, 'DoUpdates'
        end if

        if (DoUpdates) then
            if (MPIRank == 0) then
                if (Waiting) then
                    call MPI_TESTALL(MPIChains-1,req, flag, stats, ierror)
                    Waiting = .not. flag
                elseif (mod(S%Count,max(1,MPI_Sample_update_freq/(slice_fac)))==0) then
                    Waiting = .true.
                    do j=1, MPIChains-1
                        call MPI_ISSEND(MPIRank,1,MPI_INTEGER, j,0,MPI_COMM_WORLD,req(j),ierror)
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
            if (Feedback > 0) write (*,*) 'Chain',instance,' MPI communicating'
            allocate(MPIcovmats(num_params_used,num_params_used,MPIChains))
            allocate(MPIMeans(0:num_params_used,MPIChains))
            MPICovMat = 0
            MPIMean = 0
            MPImean(0) = S%Count - S%Count/2 + 1
            do i = S%Count/2, S%Count
                MPiMean(1:num_params_used) = MPiMean(1:num_params_used) + S%Item(i)
            end do
            MPiMean(1:num_params_used) = MPiMean(1:num_params_used) / MPImean(0)
            do i = S%Count/2, S%Count
                delta = S%Item(i)-MPIMean(1:num_params_used)
                do j = 1, num_params_used
                    MPICovmat(:,j) =  MPICovmat(:,j) + delta*(S%Item(i,j)- MPIMean(j))
                end do
            end do
            MPICovMat = MPICovMat / MPImean(0)

            !Replace this with ALLGATHER
            !MPICovMats(:,:,instance) = MPICovMat
            !MPIMeans(:,instance) = MPIMean
            !do i=1, MPIChains
            !    j = i-1
            !    call MPI_BCAST(MPICovMats(:,:,i),Size(MPICovMat),MPI_real_mcp,j,MPI_COMM_WORLD,ierror)
            !    call MPI_BCAST(MPIMeans(:,i),Size(MPIMean),MPI_real_mcp,j,MPI_COMM_WORLD,ierror)
            !end do

            call MPI_ALLGATHER(MPICovMat,Size(MPICovMat),MPI_real_mcp,MPICovmats,Size(MPICovmat), &
            MPI_real_mcp, MPI_COMM_WORLD,ierror)
            call MPI_ALLGATHER(MPIMean,Size(MPIMean),MPI_real_mcp,MPIMeans,Size(MPIMean), &
            MPI_real_mcp, MPI_COMM_WORLD,ierror)

            if (all(MPIMeans(0,:)> MPI_Min_Sample_Update/2 + 2)) then
                !check have reasonable number of samples in each)
                norm = sum(MPIMeans(0,:))

                do i=1, num_params_used
                    mean(i) = sum(MPIMeans(i,:)*MPIMeans(0,:))/norm
                    chain_means(:,i) = MPIMeans(i,:)
                end do

                if (MPIChains > 1) then
                    do i=1,num_params_used
                        do j=i,num_params_used
                            cov(i,j) = sum(MPIMeans(0,:)*MPICovMats(i,j,:))/ norm
                            meanscov(i,j) = sum(MPIMeans(0,:)*&
                            (chain_means(:,i)-mean(i))*(chain_means(:,j)-mean(j)))/norm
                            meanscov(j,i) = meanscov(i,j)
                            cov(j,i) = cov(i,j)
                        end do
                    end do
                    MPICovMat = Cov !+ meansCov !Estimate global covariance for proposal density
                    meansCov = meansCov * real(MPIChains,mcp)/(MPIChains-1)
                    invertible= GelmanRubinEvalues(cov, meanscov, evals, num_params_used)
                    if (invertible) then
                        R = maxval(evals)
                        if (Feedback > 1 .and. MPIRank==0) write (*,*) 'Convergence e-values: ', real(evals)
                        if (Feedback > 0 .and. MPIRank==0) then
                            write (*,*) 'Current convergence R-1 = ',real(R), ' chain steps =',sample_num
                            if (use_fast_slow) write(*,*) 'slow changes', slow_changes, 'power changes', semislow_changes
                        end if
                        if (logfile_unit/=0) then
                            write(logLine,*) 'Current convergence R-1 = ',real(R), ' chain steps =',sample_num
                            if (use_fast_slow) write(logLine,*) 'slow changes', slow_changes, 'power changes', semislow_changes
                            call IO_WriteLog(logfile_unit,logLine)
                        end if
                        if (R < MPI_R_Stop .and. flukecheck) then
                            if (MPI_Check_Limit_Converge) then
                                !Now check if limits from different chains agree well enough
                                if (CheckLImitsConverge(S)) call ConvergeStatus(.true., R, 'Requested limit convergence achieved')
                            else
                                !If not also checking limits, we are done
                                call ConvergeStatus(.true., R, 'Requested convergence R achieved')
                            end if
                        end if
                        call ConvergeStatus(.false., R)
                        flukecheck = R < MPI_R_Stop
                        if (S%Count > 500000) then
                            !Try not to blow memory by storing too many samples
                            call S%Thin(2)
                            MPI_thin_fac = MPI_thin_fac*2
                        end if
                    else
                        if (Feedback > 0 .and. MPIRank==0) write(*,*) 'Covariance not currently invertible'
                        R=1e6
                    end if

                end if !MPIChains >1

                if (MPI_LearnPropose .and. ( MPIChains==1 .or. (.not. StartCovMat &
                .or. R < MPI_Max_R_ProposeUpdate) .and. R > MPI_R_StopProposeUpdate)) then
                    !If beginning to converge, update covariance matrix
                    if (Feedback > 0 .and. MPIRank==0) write (*,*) 'updating proposal density'

                    if (MPI_StartSliceSampling .and. sampling_method == sampling_slice) then
                        if (Feedback > 0 .and. MPIRank==0) write (*,*) 'Switching from slicing to Metropolis'
                        sampling_method = sampling_metropolis
                    end if

                    call Proposer%SetCovariance(MPICovMat)

                end if !update propose
            end if !all chains have enough

            deallocate(MPICovMats, MPIMeans)
        end if !flag
    end if
#endif

    end subroutine AddMPIParams

#ifdef MPI

    subroutine ConvergeStatus(isDone, R, Msg)
    logical, intent(in) :: isDone
    character(LEN=*), intent(in), optional :: Msg
    real(mcp), intent(in) :: R

    if (MPiRank==0) then
        call CreateTxtFile(trim(baseroot)//'.converge_stat',tmp_file_unit)
        write(tmp_file_unit,*) real(R)
        if (isDone) write(tmp_file_unit,'(a)') 'Done'
        call CloseFile(tmp_file_unit)
    end if
    if (isDone) call DoStop(Msg)

    end subroutine COnvergeStatus

    function CheckLimitsConverge(L)
    !Check limits from last half chains agree well enough across chains to be confident of result
    !Slowly explored tails will cause problems (long time till stops)
    Type(TSampleList), intent(in) :: L
    integer i,j, side, ierror, worsti
    real(mcp), allocatable, dimension(:,:,:) :: Limits
    real(mcp) MeanLimit, var, LimErr, WorstErr
    logical :: CheckLimitsConverge
    logical :: AllOK
    integer numCheck
    integer, allocatable, dimension(:) :: params_check
    character(LEN=128) logLine


    if (MPI_Limit_Param/=0) then
        numCheck = 1
        allocate(params_check(numCheck))
        do j=1,num_params_used
            if (params_used(j) ==MPI_Limit_Param) then
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
        call L%ConfidVal(params_check(j), MPI_Limit_Converge, L%Count/2, L%Count, &
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
            write(*,'(1A22,1A8,1f8.3)') UsedParamNameOrNumber(params_check(j)),'lim err', LimErr
        end do
    end do
    if (Feedback > 0 .and. MPIRank==0) then
        write (*,'(a)')  'Current limit err = '//trim(RealToStr(WorstErr))// &
        ' for '//trim(UsedParamNameOrNumber(Worsti))//'; samps = '//trim(IntToStr(L%Count*MPI_thin_fac))
    end if
    if (logfile_unit/=0) then
        write (logLine,'(a)') 'Current limit err = '//trim(RealToStr(WorstErr))// &
        ' for '//trim(UsedParamNameOrNumber(Worsti))//'; samps = '//trim(IntToStr(L%Count*MPI_thin_fac))
        call IO_WriteLog(logfile_unit,logLine)
    end if
    CheckLimitsConverge = (WorstErr < MPI_Limit_Converge_Err)
    deallocate(Limits)
    deallocate(params_check)
    end function CheckLimitsConverge

#endif

    function UsedParamNameOrNumber(i) result(name)
    character(len=ParamNames_maxlen)  :: name
    integer, intent(in) :: i

    name = ParamNames_NameOrNumber(NameMapping, params_used(i))

    end function UsedParamNameOrNumber


    subroutine WriteCovMat(fname, matrix)
    use IO
    integer i
    character(LEN=*), intent(in) :: fname
    character(LEN=4096) outline
    real(mcp), intent(in) :: matrix(:,:)

    if (NameMapping%nnames/=0) then
        outline=''
        do i=1, num_params_used
            outline = trim(outline)//' '//trim(UsedParamNameOrNumber(i))
        end do
        call IO_WriteProposeMatrix(matrix ,fname, outline)
    else
        call Matrix_write(fname,matrix,forcetable=.true.)
    end if
    end subroutine WriteCovMat


    subroutine WriteModel(Params, i, like, mult)
    Class(ParamSet) :: Params
    integer i
    real(mcp), intent(in) :: mult, like
    integer j , len, unused
    logical, save :: first = .true.
    Type(DataLikelihood), pointer :: DataLike

    if (first .and. new_chains) then
        first = .false.
        if (mcp==kind(1.0)) then
            j=3
        else
            j=4
        end if
        write(i) j, num_params_used
        if (.not. any (NameMapping%Name=='')) then
            write(i) .true.
            do j=1,num_params_used
                len=len_trim(NameMapping%name(params_used(j)))
                write(i) len
                write(i) NameMapping%name(params_used(j))(1:len)
            end do
        else
            write(i) .false.
        end if
        write(i) DataLikelihoods%Count
        do j=1, DataLikelihoods%Count
            DataLike => DataLikelihoods%Item(j)
            len = len_trim(dataLIke%name)
            write(i) len
            write(i) dataLIke%name(1:len)
        end do
        unused=0
        write(i) unused
    end if

    write(i) mult, like
    !  write(i) num_matter_power
    write(i) Params%Likelihoods(1:DataLikelihoods%Count)
    write(i) Params%P(params_used)

    call Params%Theory%WriteTheory(i)

    if (flush_write) call FlushFile(i)

    end subroutine WriteModel

    subroutine ReadModel(Params,  i, has_likes, mult, like, error)
    Class (ParamSet) :: Params
    integer, intent(in) :: i
    integer, intent(out) :: error
    real(mcp), intent(out) :: mult, like
    logical, intent(out) :: has_likes(:)
    real(mcp), allocatable, save :: likes(:)
    integer j, k, np, len, unused
    character(LEN=80) :: name
    logical, save :: first = .true.
    integer, save :: numlikes, tmp(1)
    integer, allocatable, save :: like_indices(:)
    Type(DataLikelihood), pointer :: DataLike
    character(LEN=ParamNames_maxlen) ::  pname
    integer, allocatable, save ::  current_param_indices(:)
    logical :: has_names

    error = 0
    if (first) then
        first = .false.
        read(i,end=100,err=100) j, np
        if (j/=3 .and. mcp==kind(1.0) .or. j/=4 .and. mcp/=kind(1.0)) &
        call MpiStop('ReadModel: wrong file format (old cosmomc version?)')
        if (np/=num_params_used) call MpiStop('ReadModel: number of parameters changed')
        read(i) has_names
        allocate(current_param_indices(num_params_used))
        current_param_indices=-1
        if (has_names) then
            do j=1,num_params_used
                read(i) len
                pname=''
                read(i) pname(1:len)
                current_param_indices(j) = ParamNames_index(NameMapping, pname)
            end do
            if (any(current_param_indices==-1)) call MpiStop('ReadModel: parameters in .data files could not be matched')
        else
            current_param_indices = params_used
        end if
        read(i) numlikes
        allocate(likes(numlikes))
        allocate(like_indices(numlikes))
        like_indices=0
        do j=1, numlikes
            read(i) len
            name=''
            read(i) name(1:len)
            do k=1, DataLikelihoods%Count
                DataLike => DataLikelihoods%Item(k)
                if (DataLike%name==name) then
                    like_indices(j)=k
                    exit
                end if
            end do
        end do
        do j=1, DataLikelihoods%Count
            has_likes(j) = any(like_indices==j)
        end do
        read(i) unused
        if (unused/=0) call MpiStop('ReadModel: Don''t know what extra info is')
        if (unused>0) read(i) tmp(1:unused)
    end if

    Params%Likelihoods=0
    read(i,end = 100, err=100) mult, like
    read(i) likes(1:numlikes)
    do j=1,numlikes
        if (like_indices(j)/=0) Params%Likelihoods(like_indices(j)) = likes(j)
    end do
    Params%P= Scales%center
    read(i) Params%P(current_param_indices)
    call Params%Theory%ReadTheory(i)

    return
100 error = 1

    end subroutine ReadModel

    subroutine WriteParams(P, mult, like)
    Type(ParamSet) P
    real(mcp), intent(in) :: mult, like
    real(mcp), allocatable :: output_array(:)
    Type(mc_real_pointer) :: derived
    integer numderived

    if (outfile_handle ==0) return

    if (generic_mcmc) then
        call IO_OutputChainRow(outfile_handle, mult, like, P%P(params_used), num_params_used)
    else
        numderived = Parameterization%CalcDerivedParams(P%P,P%Theory, derived)

        allocate(output_array(num_params_used + numderived))
        output_array(1:num_params_used) =  P%P(params_used)
        output_array(num_params_used+1:num_params_used+numderived) =  derived%P
        deallocate(derived%P)

        call IO_OutputChainRow(outfile_handle, mult, like, output_array)
        deallocate(output_array)
    end if

    end  subroutine WriteParams


    subroutine WriteParamsAndDat(P, mult, like)
    Type(ParamSet) P
    real(mcp), intent(in) :: mult, like
    real(mcp),allocatable :: output_array(:)
    Type(mc_real_pointer) :: derived
    integer numderived

    if (outfile_handle ==0) return

    numderived =Parameterization%CalcDerivedParams(P%P,P%Theory, derived)

    allocate(output_array(num_params_used + numderived + P%Theory%num_k ))
    output_array(1:num_params_used) =  P%P(params_used)
    output_array(num_params_used+1:num_params_used+numderived) =  derived%P
    deallocate(derived%P)

    output_array(num_params_used+numderived+1:num_params_used+numderived+P%Theory%num_k) = &
    P%Theory%matter_power(:,1)

    call IO_OutputChainRow(outfile_handle, mult, like, output_array)
    deallocate(output_array)

    end  subroutine WriteParamsAndDat

    subroutine OutputParamRanges(Names, fname)
    use ParamNames
    Type(TParamNames) :: Names
    character(len=*), intent(in) :: fname
    integer unit,i

    unit = new_file_unit()
    call CreateTxtFile(fname,unit)
    do i=1, Names%num_MCMC
        write(unit,'(1A22,2E17.7)') ParamNames_NameOrNumber(Names,i), Scales%PMin(i),Scales%PMax(i)
    end do
    call CloseFile(unit)

    end subroutine OutputParamRanges

    end module ParamDef
