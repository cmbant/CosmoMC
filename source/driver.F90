    program SolveCosmology
    ! This is a driving routine that illustrates the use of the program.
    use settings
    use MonteCarlo
    use ImportanceSampling
    use CalcLike
    !    use EstCovmatModule
    use minimize
    use CosmologyConfig
    implicit none

    character(LEN=:), allocatable :: InputFile, LogFile,  numstr, fname, rootdir
    Type(TSettingIni) :: Ini
    logical bad
    integer  i, numtoget, action
    Type(ParamSet) Params, EstParams
    integer file_unit, status
    real(mcp) bestfit_loglike
    integer, parameter :: action_MCMC=0, action_importance=1, action_maxlike=2, &
    action_Hessian=3
    integer unit
    logical want_minimize
    logical :: start_at_bestfit = .false.
    Class(DataLikelihood), pointer :: Like
    real(mcp) mult
    logical is_best_bestfit
    character(LEN=:), allocatable :: prop_mat
    class(TImportanceSampler), pointer :: ImportanceSampler
#ifdef MPI
    double precision intime
    integer ierror

    call mpi_init(ierror)
    if (ierror/=MPI_SUCCESS) stop 'MPI fail: rank'
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
    call mpi_comm_rank(mpi_comm_world,MPIrank,ierror)

    instance = MPIrank +1 !start at 1 for chains
    write (numstr,*) instance
    rand_inst = instance
    if (ierror/=MPI_SUCCESS) call DoAbort('MPI fail')

    call mpi_comm_size(mpi_comm_world,MPIchains,ierror)

    if (instance == 1) then
        print *, 'Number of MPI processes:',mpichains
        InputFile=GetParam(1)
    end if


    CALL MPI_Bcast(InputFile, LEN(InputFile), MPI_CHARACTER, 0, MPI_COMM_WORLD, ierror)

#endif
    call Ini%Open(InputFile)


    if (Ini%HasKey('local_dir')) LocalDir=Ini%ReadFileName('local_dir')
    if (Ini%HasKey('data_dir')) DataDir=Ini%ReadFileName('data_dir')

    if (Ini%HasKey('custom_params')) then
        fname = Ini%ReadFileName('custom_params')
        if (fname/='') then
            call CustomParams%Open(fname, bad)
            if (bad) call DoAbort('Error reading custom_params parameter file')
        end if
    end if

    action = Ini%Read_Int('action',action_MCMC)

    generic_mcmc = Ini%Read_Logical('generic_mcmc',generic_mcmc)

    AccuracyLevel = Ini%Read_Real('accuracy_level',1.)

    if (action==action_MCMC) then
        checkpoint = Ini%Read_Logical('checkpoint',.false.)
        if (checkpoint) flush_write = .true.
        start_at_bestfit= Ini%read_logical('start_at_bestfit',.false.)
        propose_scale = Ini%Read_Real('propose_scale',2.4)


    end if


    stop_on_error = Ini%Read_logical('stop_on_error',stop_on_error)

    baseroot = Ini%ReadFileName(Ini,'file_root', NotFoundFail = .true.)
    if(instance<=1) then
        write(*,*) 'file_root:'//trim(baseroot)
    end if
    if (Ini%HasKey('root_dir')) then
        !Begin JD modifications for output of filename in output file
        rootdir = Ini%ReadFileName(Ini,'root_dir')
        baseroot = trim(rootdir)//trim(baseroot)
    end if

    rootname = trim(baseroot)

    FileChangeIniAll = trim(rootname)//'.read'

    if (instance /= 0) then
        rootname = trim(rootname)//'_'//trim(adjustl(numstr))
    end if

    new_chains = .true.

    if (checkpoint) then
#ifdef MPI
        new_chains = .not. IO_Exists(trim(rootname) //'.chk')
#else
        new_chains = .not. IO_Exists(trim(rootname) //'.txt')
        if(.not. new_chains) new_chains = IO_Size(trim(rootname) //'.txt')<=0
#endif
    end if


    if (action == action_importance) then
        allocate(TImportanceSampler::ImportanceSampler)
        call ImportanceSampler%ReadPostParams(baseroot)
    end if
    FeedBack = Ini%Read_Int('feedback',0)
    FileChangeIni = trim(rootname)//'.read'

    if (action == action_MCMC) then
        LogFile = trim(rootname)//'.log'

        if (LogFile /= '') then
            logfile_unit = IO_OutputOpenForWrite(LogFile, append=.not. new_chains, isLogFile = .true.)
        else
            logfile_unit = 0
        end if

        indep_sample = Ini%Read_Int('indep_sample')
        if (indep_sample /=0) then
            fname = trim(rootname)//'.data'
            indepfile_handle = IO_DataOpenForWrite(fname, append = .not. new_chains)
        end if

        burn_in = Ini%Read_Int('burn_in',0)
        sampling_method = Ini%Read_Int('sampling_method',sampling_metropolis)
        if (sampling_method > 7 .or. sampling_method<1) call DoAbort('Unknown sampling method')
        if (sampling_method==sampling_slowgrid) directional_grid_steps = Ini%Read_Int('directional_grid_steps',20)
        if (sampling_method==sampling_fast_dragging) then
            dragging_steps = Ini%Read_Real('dragging_steps',2.)
            use_fast_slow = .true.
        else
            use_fast_slow = Ini%read_Logical('use_fast_slow',.true.)
        end if
    else
        if (action == action_maxlike) use_fast_slow = Ini%read_Logical('use_fast_slow',.true.)
    end if

    numstr = Ini%Read_String('rand_seed')
    if (numstr /= '') then
        read(numstr,*) i
        call InitRandom(i)
    else
        call InitRandom()
    end if

    if (.not. generic_mcmc) then

    call SetDataLikelihoods(Ini)
    call DataLikelihoods%CheckAllConflicts
    if(use_LSS) call Initialize_PKSettings()
    end if


    num_threads = Ini%Read_Int('num_threads',0)
    !$ if (num_threads /=0) call OMP_SET_NUM_THREADS(num_threads)

    estimate_propose_matrix = Ini%Read_Logical('estimate_propose_matrix',.false.)
    if (estimate_propose_matrix) then
        if (Ini%Read_String('propose_matrix') /= '') &
        call DoAbort('Cannot have estimate_propose_matrix and propose_matrix')
    end if
    want_minimize = action == action_maxlike .or. action==action_Hessian &
    .or. action == action_MCMC .and. estimate_propose_matrix .or. &
    start_at_bestfit .and. new_chains

    numtoget = Ini%Read_Int('samples')

    call SetTheoryParameterization(Ini, NameMapping)
    call DataLikelihoods%AddNuisanceParameters(NameMapping)
    if (.not. generic_mcmc) call CosmoTheory_ReadParams(Ini)

    call BaseParams%InitializeUsedParams(Ini,Params)
    call BaseParams%SetFastSlowParams(Ini, use_fast_slow)
    call SampleCollector%ReadParams(Ini,action==action_MCMC)


    if (action /= action_importance) then
        prop_mat = Ini%Read_String('propose_matrix',.false.)
        if (prop_mat /= '' .and. prop_mat(1:1) /= '/') prop_mat = concat(LocalDir,prop_mat)
    else
        prop_mat=''
    end if

    call BaseParams%SetCovmat(prop_mat)

    if (BaseParams%covariance_has_new) this%Mpi%MPI_Max_R_ProposeUpdate = this%Mpi%MPI_Max_R_ProposeUpdateNew
    if (sampling_method /=sampling_fast_dragging) then
        this%Mpi%MPI_Thin_fac = this%Mpi%MPI_Thin_fac*Sampler%Oversample_fast
    end if


    call Sampler%SetCovariance(BaseParms%covariance_estimate)

    call LikeCalculator%ReadParams(Ini)


    if (want_minimize) call Minimize_ReadIni(Ini)

    if (MpiRank==0) then
        call Ini%ReadValues%Add('CosmoMC_Version',CosmoMC_Version)
        do i=1, DataLikelihoods%Count
            like => DataLikelihoods%Item(i)
            if (like%version/='')  call Ini%ReadValues%Add( &
            concat('Compiled_data_',like%name),like%version)
        end do
        call TheoryCalculator%VersionTraceOutput(Ini%ReadValues)
        if (action==action_importance) then
            call Ini%SaveReadValues(trim(PostParams%redo_outroot) //'.inputparams')
        else if (action==action_maxlike .or. action==action_Hessian) then
            call Ini%SaveReadValues(trim(baseroot) //'.minimum.inputparams')
        else
            call Ini%SaveReadValues(trim(baseroot) //'.inputparams')
        end if
    end if

    call Ini%Close()

    if (MpiRank==0 .and. action==action_MCMC .and. NameMapping%nnames/=0) then
        call IO_OutputParamNames(NameMapping,trim(baseroot), params_used, add_derived = .true.)
        call BaseParams%OutputParamRanges(NameMapping, trim(baseroot)//'.ranges')
    end if

    call SetIdlePriority !If running on Windows

    if (want_minimize) then
        !New Powell 2009 minimization, AL Sept 2012, update Sept 2013
        if (action /= action_MCMC .and. MPIchains>1 .and. .not. minimize_uses_MPI) call DoAbort( &
        'Mimization only uses one MPI thread, use -np 1 or compile without MPI (don''t waste CPUs!)')
        if (MpiRank==0) write(*,*) 'finding best fit point...'
        if (minimize_uses_MPI .or. MpiRank==0) then
            bestfit_loglike = Minimize%FindBestFit(Params,is_best_bestfit)
            if (is_best_bestfit) then
                if (bestfit_loglike==logZero) write(*,*) MpiRank,'WARNING: FindBestFit did not converge'
                if (Feedback >0) write(*,*) 'Best-fit results: '
                call WriteBestFitParams(bestfit_loglike,Params, trim(baseroot)//'.minimum')
                call DataLikelihoods%WriteDataForLikelihoods(Params%P, Params%Theory, trim(baseroot))
                if (use_CMB) call Params%Theory%WriteBestFitData(trim(baseroot))
                if (action==action_maxlike) call DoStop('Wrote the minimum to file '//trim(baseroot)//'.minimum')
            else
                if (action==action_maxlike) call DoStop()
            end if
        end if
#ifdef MPI
        if (.not. minimize_uses_MPI) then
            CALL MPI_Bcast(Params%P, size(Params%P), MPI_real_mcp, 0, MPI_COMM_WORLD, ierror)
        end if
#endif
        Scales%center(1:num_params) = Params%P(1:num_params)
    end if

    if (estimate_propose_matrix .and. action == action_MCMC .or. action==action_Hessian) then
        ! slb5aug04 with AL updates
        if (MpiRank==0) then
            EstParams = Params
            write (*,*) 'Estimating propose matrix from Hessian at bfp...'
            initial_propose_matrix=EstCovmat(EstParams,4._mcp,status)
            ! By default the grid used to estimate the covariance matrix has spacings
            ! such that deltaloglike ~ 4 for each parameter.
            call EstParams%Clear(keep=Params)
            if (status==0) call DoAbort('estimate_propose_matrix: estimating propose matrix failed')
            if (Feedback>0) write (*,*) 'Estimated covariance matrix:'
            call WriteCovMat(trim(baseroot) //'.hessian.covmat', initial_propose_matrix)
            write(*,*) 'Wrote the local inv Hessian to file ',trim(baseroot)//'.hessian.covmat'
            if (action==action_Hessian) call DoStop
        end if
#ifdef MPI
        CALL MPI_Bcast(initial_propose_matrix, size(initial_propose_matrix), MPI_real_mcp, 0, MPI_COMM_WORLD, ierror)
#endif

        call Proposer%SetCovariance(initial_propose_matrix)
    end if

    if (action == action_MCMC) then
        fname = trim(rootname)//'.txt'
        if (new_chains) then
            call SetStartPositions(Params)
        else
            call IO_ReadLastChainParams(fname, mult, StartLike, Params%P, params_used)
            call AddMPIParams(Params%P, StartLike, .true.)
        end if
        outfile_handle = IO_OutputOpenForWrite(fname, append = .not. new_chains)

        if (Feedback > 0 .and. MPIRank==0) write (*,*) 'starting Monte-Carlo'
        call MCMC%Sample_From(Params, numtoget)

        if (Feedback > 0) write (*,*) 'finished'

        if (logfile_unit /=0) call IO_Close(logfile_unit, isLogFile=  .true.)
        if (indepfile_handle /=0) call IO_DataCloseWrite(indepfile_handle)

        close(outfile_handle)
    else if (action==action_importance) then
        if (Feedback > 0 .and. MPIRank==0) write (*,*) 'starting post processing'
        call postprocess(rootname)
        call DoStop('Postprocesing done',.false.)
    else
        call DoAbort('undefined action')
    end if

    call DoStop('Total requested samples obtained',.true.)

    end program SolveCosmology