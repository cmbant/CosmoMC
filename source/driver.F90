    program SolveCosmology
    ! This is a driving routine that illustrates the use of the program.

    use IniFile
    use MonteCarlo
    use ParamDef
    use settings
    use cmbdata
    use posthoc
    use CalcLike
    use EstCovmatModule
    use minimize
    use mpk
    use MatrixUtils
    use IO
    use ParamNames
    use GaugeInterface, only : Eqns_name
    use DefineParameterization

    implicit none

    character(LEN=Ini_max_string_len) InputFile, LogFile

    logical bad
    integer  i, numtoget, action
    character(LEN=Ini_max_string_len)  numstr, fname
    Type(ParamSet) Params, EstParams
    integer file_unit, status
    real(mcp) bestfit_loglike
    real(mcp) max_like_radius
    integer max_like_iterations
    integer, parameter :: action_MCMC=0, action_importance=1, action_maxlike=2, &
    action_Hessian=3
    integer unit
    logical want_minimize
    logical :: start_at_bestfit = .false.
    Class(DataLikelihood), pointer :: Like
    real(mcp) mult
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
    call IO_Ini_Load(DefIni,InputFile, bad)

    if (bad) call DoAbort('Error opening parameter file')

    Ini_fail_on_not_found = .false.

    call TNameValueList_Init(CustomParams%L)
    call TNameValueList_Init(CustomParams%ReadValues)

    if (Ini_HasKey('local_dir')) LocalDir=ReadIniFileName(DefIni,'local_dir')
    if (Ini_HasKey('data_dir')) DataDir=ReadIniFileName(DefIni,'data_dir')

    if (Ini_HasKey('custom_params')) then
        fname = ReadIniFileName(DefIni,'custom_params')
        if (fname/='') then
            file_unit = new_file_unit()
            call  Ini_Open_File(CustomParams,  fname, file_unit, bad)
            call ClearFileUnit(file_unit)
            if (bad) call DoAbort('Error reading custom_params parameter file')
        end if
    end if

    action = Ini_Read_Int('action',action_MCMC)

    propose_scale = Ini_Read_Real('propose_scale',2.4)

    HighAccuracyDefault = Ini_Read_Logical('high_accuracy_default',.false.)
    AccuracyLevel = Ini_Read_Real('accuracy_level',1.)

    if (Ini_HasKey('highL_unlensed_cl_template')) then
        highL_unlensed_cl_template=  ReadIniFilename(DefIni,'highL_unlensed_cl_template')
    else
        highL_unlensed_cl_template = concat(LocalDir,'camb/',highL_unlensed_cl_template)
    end if

    if (action==action_MCMC) then
        checkpoint = Ini_Read_Logical('checkpoint',.false.)
        if (checkpoint) flush_write = .true.
        start_at_bestfit= Ini_read_logical('start_at_bestfit',.false.)
    end if


#ifdef MPI
    MPI_StartTime = MPI_WTime()
    if (action==action_MCMC) then
        MPI_R_Stop = Ini_Read_Real('MPI_Converge_Stop',MPI_R_Stop)
        MPI_LearnPropose = Ini_Read_Logical('MPI_LearnPropose',.true.)
        if (MPI_LearnPropose) then
            MPI_R_StopProposeUpdate=Ini_Read_Real('MPI_R_StopProposeUpdate',0.)
            MPI_Max_R_ProposeUpdate=Ini_Read_Real('MPI_Max_R_ProposeUpdate',2.)
            MPI_Max_R_ProposeUpdateNew=Ini_Read_Real('MPI_Max_R_ProposeUpdateNew',30.)
        end if
        MPI_Check_Limit_Converge = Ini_Read_Logical('MPI_Check_Limit_Converge',.false.)
        MPI_StartSliceSampling = Ini_Read_Logical('MPI_StartSliceSampling',.false.)
        if (MPI_Check_Limit_Converge) then
            MPI_Limit_Converge = Ini_Read_Real('MPI_Limit_Converge',0.025)
            MPI_Limit_Converge_Err = Ini_Read_Real('MPI_Limit_Converge_Err',0.3)
            MPI_Limit_Param = Ini_Read_Int('MPI_Limit_Param',0)
        end if
    end if

#endif

    stop_on_error = Ini_Read_logical('stop_on_error',stop_on_error)

    Ini_fail_on_not_found = .true.

    baseroot = ReadIniFileName(DefIni,'file_root')

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
#endif
    end if


    if (action == action_importance) call ReadPostParams(baseroot)

    FeedBack = Ini_Read_Int('feedback',0)
    FileChangeIni = trim(rootname)//'.read'

    if (action == action_MCMC) then
        LogFile = trim(rootname)//'.log'

        if (LogFile /= '') then
            logfile_unit = IO_OutputOpenForWrite(LogFile, append=.not. new_chains, isLogFile = .true.)
        else
            logfile_unit = 0
        end if

        indep_sample = Ini_Read_Int('indep_sample')
        if (indep_sample /=0) then
            fname = trim(rootname)//'.data'
            indepfile_handle = IO_DataOpenForWrite(fname, append = .not. new_chains)
            !          call CreateOpenFile(fname,indepfile_unit,'unformatted',.not. new_chains)
        end if

        Ini_fail_on_not_found = .false.
        burn_in = Ini_Read_Int('burn_in',0)
        sampling_method = Ini_Read_Int('sampling_method',sampling_metropolis)
        if (sampling_method > 7 .or. sampling_method<1) call DoAbort('Unknown sampling method')
        if (sampling_method==sampling_slowgrid) directional_grid_steps = Ini_Read_Int('directional_grid_steps',20)
        if (sampling_method==sampling_fast_dragging) then
            dragging_steps = Ini_Read_Real('dragging_steps',2.)
            use_fast_slow = .true.
        else
            use_fast_slow = Ini_read_Logical('use_fast_slow',.true.)
        end if
    else
        if (action == action_maxlike) use_fast_slow = Ini_read_Logical('use_fast_slow',.true.)
        Ini_fail_on_not_found = .false.
    end if

    numstr = Ini_Read_String('rand_seed')
    if (numstr /= '') then
        read(numstr,*) i
        call InitRandom(i)
    else
        call InitRandom()
    end if

    use_nonlinear = Ini_Read_Logical('nonlinear_pk',.false.)
    pivot_k = Ini_Read_Real('pivot_k',0.05)
    inflation_consistency = Ini_read_Logical('inflation_consistency',.false.)
    bbn_consistency = Ini_Read_Logical('bbn_consistency',.true.)
    num_massive_neutrinos = Ini_read_int('num_massive_neutrinos',-1)

    call SetDataLikelihoods(DefIni)

    Temperature = Ini_Read_Real('temperature',1.)

    num_threads = Ini_Read_Int('num_threads',0)
    !$ if (num_threads /=0) call OMP_SET_NUM_THREADS(num_threads)

    estimate_propose_matrix = Ini_Read_Logical('estimate_propose_matrix',.false.)
    if (estimate_propose_matrix) then
        if (Ini_Read_String('propose_matrix') /= '') &
        call DoAbort('Cannot have estimate_propose_matrix and propose_matrix')
    end if
    want_minimize = action == action_maxlike .or. action==action_Hessian &
    .or. action == action_MCMC .and. estimate_propose_matrix .or. &
    start_at_bestfit .and. new_chains

    if (want_minimize) then
        !radius in normalized parameter space to converge
        max_like_radius = Ini_Read_Real('max_like_radius',0.01)
        max_like_iterations = Ini_Read_Int('max_like_iterations',6000)
        !set points factor above 2 to use a denser sampling of space (may be more robust)
        minimization_points_factor = Ini_Read_Int('minimization_points_factor',minimization_points_factor)
        !will exit if function difference between iterations less than minimize_loglike_tolerance (even if radius criterion not met)
        minimize_loglike_tolerance = Ini_Read_double('minimize_loglike_tolerance',minimize_loglike_tolerance)
    end if

    Ini_fail_on_not_found = .true.

    !    'This version of CosmoMC does not currently have mpk implemented'
    !    nummpksets = Ini_Read_Int('mpk_numdatasets',0)
    !    if (Use_mpk) then
    !do i= 1, nummpksets
    !    mpk_filename(i) = ReadIniFileName(DefIni,numcat('mpk_dataset',i))
    !    call ReadMpkDataset(mpk_filename(i))
    !end do
    !if (Feedback>1) write(*,*) 'read mpk datasets'
    !    end if

    numtoget = Ini_Read_Int('samples')

    call SetTheoryParameterization(DefIni, NameMapping)
    call DataLikelihoods%AddNuisanceParameters(NameMapping)
    call CMB_Initialize(Params%Info)
    call InitializeUsedParams(DefIni,Params)

    if (MpiRank==0) then
        do i=1, DataLikelihoods%Count
            like => DataLikelihoods%Item(i)
            if (like%version/='')  call TNameValueList_Add(DefIni%ReadValues, &
            concat('Compiled_data_',like%name),like%version)
        end do
        call TNameValueList_Add(DefIni%ReadValues, 'Compiled_CAMB_version', version)
        call TNameValueList_Add(DefIni%ReadValues, 'Compiled_Recombination', Recombination_Name)
        call TNameValueList_Add(DefIni%ReadValues, 'Compiled_Equations', Eqns_name)
        call TNameValueList_Add(DefIni%ReadValues, 'Compiled_Reionization', Reionization_Name)
        call TNameValueList_Add(DefIni%ReadValues, 'Compiled_InitialPower', Power_Name)
        unit = new_file_unit()
        if (action==action_importance) then
            call Ini_SaveReadValues(trim(PostParams%redo_outroot) //'.inputparams',unit)
        else if (action==action_maxlike .or. action==action_Hessian) then
            call Ini_SaveReadValues(trim(baseroot) //'.minimum.inputparams',unit)
        else
            call Ini_SaveReadValues(trim(baseroot) //'.inputparams',unit)
        end if
        call ClearFileUnit(unit)
    end if

    call Ini_Close

    if (MpiRank==0 .and. action==action_MCMC .and. NameMapping%nnames/=0) then
        call IO_OutputParamNames(NameMapping,trim(baseroot), params_used, add_derived = .true.)
        call OutputParamRanges(NameMapping, trim(baseroot)//'.ranges')
    end if

    call SetIdlePriority !If running on Windows

    if (want_minimize) then
        !New Powell 2009 minimization, AL Sept 2012
        if (action /= action_MCMC .and. MPIchains>1) call DoAbort( &
        'Mimization only uses one MPI thread, use -np 1 or compile without MPI (don''t waste CPUs!)')
        if (MpiRank==0) then
            write(*,*) 'finding best fit point...'
            Params%P(1:num_params) = Scales%center(1:num_params)
            bestfit_loglike = FindBestFit(Params,max_like_radius,max_like_iterations)
            if (bestfit_loglike==logZero) write(*,*) 'WARNING: FindBestFit did not converge'
            if (Feedback >0) write(*,*) 'Best-fit results: '
            call WriteBestFitParams(bestfit_loglike,Params, trim(baseroot)//'.minimum')
            if (use_CMB) call Params%Theory%WriteBestFitData(trim(baseroot))
            if (action==action_maxlike) call DoStop('Wrote the minimum to file '//trim(baseroot)//'.minimum')
        end if
#ifdef MPI
        CALL MPI_Bcast(Params%P, size(Params%P), MPI_real_mcp, 0, MPI_COMM_WORLD, ierror)
#endif
        Scales%center(1:num_params) = Params%P(1:num_params)
    end if

    if (estimate_propose_matrix .and. action == action_MCMC .or. action==action_Hessian) then
        ! slb5aug04 with AL updates
        if (MpiRank==0) then
            EstParams = Params
            write (*,*) 'Estimating propose matrix from Hessian at bfp...'
            propose_matrix=EstCovmat(EstParams,4._mcp,status)
            ! By default the grid used to estimate the covariance matrix has spacings
            ! such that deltaloglike ~ 4 for each parameter.          
            call AcceptReject(.true., EstParams%Info, Params%Info)
            if (status==0) call DoAbort('estimate_propose_matrix: estimating propose matrix failed')
            if (Feedback>0) write (*,*) 'Estimated covariance matrix:'
            call WriteCovMat(trim(baseroot) //'.hessian.covmat', propose_matrix)
            write(*,*) 'Wrote the local inv Hessian to file ',trim(baseroot)//'.hessian.covmat'
            if (action==action_Hessian) call DoStop
        end if
#ifdef MPI
        CALL MPI_Bcast(propose_matrix, size(propose_matrix), MPI_real_mcp, 0, MPI_COMM_WORLD, ierror)
#endif
        call Proposer%SetCovariance(propose_matrix)
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
        call MCMCSample(Params, numtoget)

        if (Feedback > 0) write (*,*) 'finished'

        if (logfile_unit /=0) call IO_Close(logfile_unit, isLogFile=  .true.)
        if (indepfile_handle /=0) call IO_DataCloseWrite(indepfile_handle)

        call IO_Close(outfile_handle)

    else if (action==action_importance) then
        if (Feedback > 0 .and. MPIRank==0) write (*,*) 'starting post processing'
        call postprocess(rootname)
        call DoStop('Postprocesing done',.false.)

    else
        call DoAbort('undefined action')
    end if

    call DoStop('Total requested samples obtained',.true.)

    end program SolveCosmology