    program CosmoMC
    ! This is a driving routine that illustrates the use of the program.
    use settings
    use ImportanceSampling
    use CalcLike
    !    use EstCovmatModule
    use minimize
    use CosmologyConfig
    use GeneralSetup
    use DataLikelihoodList
    use RandUtils
    implicit none

    character(LEN=:), allocatable :: LogFileName,  numstr, fname, rootdir
    character(LEN=:), allocatable :: InputFile, test_output_root, test_used_params
    Type(TSettingIni) :: Ini
    integer  i, status
    Type(ParamSet) Params !, EstParams
    real(mcp) bestfit_loglike
    logical want_minimize
    logical is_best_bestfit
    logical :: start_at_bestfit = .false.
    Class(TDataLikelihood), pointer :: Like
    character(LEN=:), allocatable :: prop_mat
    class(TMinimizer), allocatable :: Minimizer
    logical :: estimate_propose_matrix = .false.
    real(mcp), allocatable :: test_paramvals(:)
    real(mcp) :: test_check_compare = logZero
    
#ifdef MPI
    integer ierror

    call mpi_init(ierror)
    if (ierror/=MPI_SUCCESS) stop 'MPI fail: rank'
    MPIRun_start_time = MPI_WTime()
#endif

    instance = 0
#ifndef MPI
    InputFile = GetParam(1)
    if (InputFile == '') call DoAbort('No parameter input file')
#endif
    numstr = GetParam(2)
    if (numstr /= '') then
        read(numstr,*) instance
        rand_inst = instance
    end if
    call InitializeGlobalSettingDefaults()

#ifdef MPI

    if (instance /= 0) call DoAbort('With MPI should not have second parameter')
    call MpiStat(MpiRank, Mpichains)

    instance = MPIrank +1 !start at 1 for chains
    rand_inst = instance

    if (instance == 1) then
        print *, 'Number of MPI processes:',mpichains
        InputFile=GetParam(1)
    end if
    call MpiShareString(InputFile,0)

#endif
    call Ini%Open(InputFile)

    if (Ini%HasKey('local_dir')) LocalDir=Ini%ReadFileName('local_dir')
    if (Ini%HasKey('data_dir')) DataDir=Ini%ReadFileName('data_dir')

    if (Ini%HasKey('custom_params')) then
        fname = Ini%ReadFileName('custom_params')
        if (fname/='') then
            call CustomParams%Open(fname)
        end if
    end if

    call Ini%Read('accuracy_level',AccuracyLevel)
    call Ini%Read('stop_on_error',stop_on_error)

    baseroot = Ini%ReadFileName('file_root', NotFoundFail = .true.)
    if(instance<=1) then
        write(*,*) 'file_root:'//baseroot
    end if
    if (Ini%HasKey('root_dir')) then
        !Begin JD modifications for output of filename in output file
        rootdir = Ini%ReadFileName('root_dir')
        baseroot = trim(rootdir)//baseroot
    end if

    rootname = baseroot

    FileChangeIniAll = rootname//'.read'

    if (instance /= 0) then
        rootname = rootname//'_'//IntToStr(Instance)
    end if

    FeedBack = Ini%Read_Int('feedback',0)

    call Ini%Read('generic_mcmc',generic_mcmc)
    if (generic_mcmc) then
        allocate(TSetup::Setup)
    else
        allocate(TCosmologySetup::Setup)
    end if
    call Setup%Init()
    call Setup%ReadParams(Ini)

    if (Setup%action==action_MCMC) then
        checkpoint = Ini%Read_Logical('checkpoint',.false.)
        flush_write =  Ini%Read_Logical('flush_write',checkpoint)
        start_at_bestfit = Ini%Read_logical('start_at_bestfit',.false.)
    else if (Setup%action==action_tests) then
        test_output_root = Ini%Read_String('test_output_root')
        test_used_params = Ini%Read_String('test_used_params') !If not specified, just take from central value
        test_check_compare = Ini%Read_Double('test_check_compare',logZero)
    end if

    new_chains = .true.

    if (checkpoint) then
#ifdef MPI
        new_chains = .not. File%Exists(rootname //'.chk')
#else
        stop 'Checkpointing without MPI not supported'
        !new_chains = FileSize(rootname //'.txt')<=0
#endif
    end if

    FileChangeIni = rootname//'.read'

    if (Setup%action == action_MCMC .and. .not. Ini%Read_Logical('no_log', .false.)) then
        LogFileName = rootname//'.log'
        call LogFile%CreateOpenFile(LogFileName,append=.not. new_chains)
    end if

    i = Ini%Read_Int('rand_seed',-1)
    if (i/=-1) then
        call InitRandom(i)
    else
        call InitRandom()
    end if

    if (.not. generic_mcmc) then
        call SetDataLikelihoods(Ini)
        call DataLikelihoods%CheckAllConflicts
        call Setup%Config%InitForLikelihoods
    end if

    num_threads = Ini%Read_Int('num_threads',0)
    !$ if (num_threads /=0) call OMP_SET_NUM_THREADS(num_threads)

    estimate_propose_matrix = Ini%Read_Logical('estimate_propose_matrix',.false.)
    if (estimate_propose_matrix) then
        if (Ini%Read_String('propose_matrix') /= '') &
            call DoAbort('Cannot have estimate_propose_matrix and propose_matrix')
    end if
    want_minimize = Setup%action == action_maxlike .or. Setup%action==action_Hessian &
        .or. Setup%action == action_MCMC .and. estimate_propose_matrix .or. &
        start_at_bestfit .and. new_chains

    call Setup%Config%SetTheoryParameterization(Ini, BaseParams%NameMapping,'theta')

    call DataLikelihoods%AddNuisanceParameters(BaseParams%NameMapping)
    call DataLikelihoods%AddOutputLikelihoodParams(BaseParams%NameMapping)

    call BaseParams%InitializeUsedParams(Ini)
    call BaseParams%SetFastSlowParams(Ini, use_fast_slow)
    if (BaseParams%block_semi_fast .and. CosmoSettings%use_nonlinear_lensing &
        & .and. BaseParams%num_slow/=BaseParams%num_semi_slow) &
        call MpiStop('use_nonlinear_lensing does not work rigorously with block_semi_fast')

    Params%P(1:num_params) = BaseParams%Center(1:num_params)
    if (Setup%action /= action_importance) then
        prop_mat = Ini%Read_String('propose_matrix',.false.)
        if (prop_mat /= '' .and. prop_mat(1:min(len(prop_mat),1)) /= '/') prop_mat = concat(LocalDir,prop_mat)
    else
        prop_mat=''
    end if
    call BaseParams%SetCovmat(prop_mat)

    if (want_minimize) then
        call Setup%GetMinimizer(Minimizer)
        call Minimizer%ReadParams(Ini)
    end if
    if (MpiRank==0) then
        call Ini%ReadValues%Add('CosmoMC_Version',CosmoMC_Version)
        do i=1, DataLikelihoods%Count
            like => DataLikelihoods%Item(i)
            if (like%version/='') call Ini%ReadValues%Add(concat('Compiled_data_',like%name),like%version)
        end do
        if (associated(Setup%Config%Calculator)) call Setup%Config%Calculator%VersionTraceOutput(Ini%ReadValues)
        if (Setup%action==action_importance) then
            call Ini%SaveReadValues(trim(Setup%ImportanceSampler%redo_outroot) //'.inputparams')
        else if (Setup%action==action_maxlike .or. Setup%action==action_Hessian) then
            call Ini%SaveReadValues(baseroot //'.minimum.inputparams')
        else
            call Ini%SaveReadValues(baseroot //'.inputparams')
        end if
    end if

    call Ini%Close()

    call Setup%DoneInitialize()


    if (MpiRank==0 .and. Setup%action==action_MCMC .and. BaseParams%NameMapping%nnames/=0) then
        call BaseParams%OutputParamNames(baseroot, params_used, add_derived = .true.)
        call BaseParams%OutputParamRanges(baseroot)
    end if

    if (MpiRank==0 .and. Setup%action==action_MCMC) then
        call DataLikelihoods%OutputDescription(baseroot)
    end if

    call SetIdlePriority !If running on Windows

    if (allocated(Minimizer)) then
        !New Powell 2009 minimization, AL Sept 2012, update Sept 2013
        if (Setup%action /= action_MCMC .and. MPIchains>1 .and. .not. Minimizer%uses_MPI) &
            & call DoAbort('Mimization only uses one MPI thread, use -np 1 or compile without MPI (don''t waste CPUs!)')
        if (MpiRank==0) write(*,*) 'finding best fit point...'
        if (minimizer%uses_MPI .or. MpiRank==0) then
            bestfit_loglike = Minimizer%FindBestFit(Params,is_best_bestfit)
            if (is_best_bestfit) then
                if (bestfit_loglike==logZero) write(*,*) MpiRank,'WARNING: FindBestFit did not converge'
                if (Feedback >0) write(*,*) 'Best-fit results: '
                call Minimizer%LikeCalculator%WriteParamsHumanText(baseroot//'.minimum', Params, bestfit_loglike)
                if (Feedback>0) call Minimizer%LikeCalculator%WriteParamsHumanText(stdout,Params, bestfit_loglike)
                call Minimizer%LikeCalculator%WriteParamPointTextData(baseroot//'.minimum', Params)
                if (Setup%action==action_maxlike) call DoStop('Wrote the minimum to file '//baseroot//'.minimum')
            else
                if (Setup%action==action_maxlike) call DoStop()
            end if
        end if
#ifdef MPI
        if (.not. Minimizer%uses_MPI) then
            CALL MPI_Bcast(Params%P, size(Params%P), MPI_real_mcp, 0, MPI_COMM_WORLD, ierror)
        end if
#endif
        BaseParams%center(1:num_params) = Params%P(1:num_params)
    end if

    if (estimate_propose_matrix .and. Setup%action == action_MCMC .or. Setup%action==action_Hessian) then
        call MpiStop('hessian evaluation disabled for now')
        !        ! slb5aug04 with AL updates
        !        if (MpiRank==0) then
        !            EstParams = Params
        !            write (*,*) 'Estimating propose matrix from Hessian at bfp...'
        !            initial_propose_matrix=EstCovmat(EstParams,4._mcp,status)
        !            ! By default the grid used to estimate the covariance matrix has spacings
        !            ! such that deltaloglike ~ 4 for each parameter.
        !            call EstParams%Clear(keep=Params)
        !            if (status==0) call DoAbort('estimate_propose_matrix: estimating propose matrix failed')
        !            if (Feedback>0) write (*,*) 'Estimated covariance matrix:'
        !            call WriteCovMat(baseroot //'.hessian.covmat', initial_propose_matrix)
        !            write(*,*) 'Wrote the local inv Hessian to file ',baseroot//'.hessian.covmat'
        !            if (action==action_Hessian) call DoStop
        !        end if
        !#ifdef MPI
        !        CALL MPI_Bcast(initial_propose_matrix, size(initial_propose_matrix), MPI_real_mcp, 0, MPI_COMM_WORLD, ierror)
        !#endif
        !
        !        call Proposer%SetCovariance(initial_propose_matrix)
    end if

    if (Setup%action == action_MCMC) then
        call Setup%DoSampling(Params)
    else if (Setup%action==action_importance) then
        if (Feedback > 0 .and. MPIRank==0) write (*,*) 'starting post processing'
        call Setup%ImportanceSampler%ImportanceSample(rootname)
        call DoStop('Postprocesing done',.false.)
    else if (Setup%action == action_tests) then
        allocate(test_paramvals, source = BaseParams%Center(:num_params))
        if (test_used_params/='') then
            read(test_used_params,*, iostat=status) test_paramvals(params_used)
            if (status/=0) call MpiStop('Error reading test_used_params array')
        end if
        call Setup%DoTests(test_output_root,test_paramvals, test_check_compare)
    else
        call DoAbort('undefined action')
    end if
    call LogFile%Close()

    call DoStop('Total requested samples obtained',.true.)

    contains


    subroutine SetIdlePriority
#ifdef RUNIDLE
    USE DFWIN
    Integer dwPriority
    Integer CheckPriority

    dwPriority = 64 ! idle priority
    CheckPriority = SetPriorityClass(GetCurrentProcess(), dwPriority)
#endif
    end subroutine SetIdlePriority

    end program CosmoMC
