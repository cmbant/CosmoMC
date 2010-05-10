

program SolveCosmology
        use IniFile
        use MonteCarlo
        use ParamDef
        use settings
        use cmbdata
        use posthoc
        use WeakLen
        use CalcLike
        use EstCovmatModule
        use ConjGradModule
        use mpk
        use MatrixUtils
#ifdef WMAP_PARAMS
        use WMAP_OPTIONS
#endif
        implicit none
      
        
! This is a driving routine that illustrates the use of the program.
        character(LEN=Ini_max_string_len) InputFile, LogFile


        logical bad, est_bfp_before_covmat
        integer numsets, nummpksets, i, numtoget, action
        character(LEN=Ini_max_string_len) baseroot, filename(100), &
         mpk_filename(100),  SZTemplate(100), numstr, fname
        real SZscale(100)
        Type(ParamSet) Params, EstParams
        integer num_points
        integer status          
        integer file_unit
        real delta_loglike
        

#ifdef MPI
        double precision intime
        integer ierror

 call mpi_init(ierror)
 if (ierror/=MPI_SUCCESS) stop 'MPI fail: rank'
#endif


#ifndef MPI 
        InputFile = GetParam(1)
        if (InputFile == '') call DoStop('No parameter input file')
#endif
        numstr = GetParam(2)
        if (numstr /= '') then
         read(numstr,*) instance
         rand_inst = instance   
        else
         instance = 0
        end if

#ifdef MPI 

        if (instance /= 0) call DoStop('With MPI should not have second parameter')
        call mpi_comm_rank(mpi_comm_world,MPIrank,ierror)


        instance = MPIrank +1 !start at 1 for chains                                             
        write (numstr,*) instance
        rand_inst = instance
        if (ierror/=MPI_SUCCESS) call DoStop('MPI fail')

        call mpi_comm_size(mpi_comm_world,MPIchains,ierror)

        if (instance == 1) then
          print *, 'Number of MPI processes:',mpichains
          InputFile=GetParam(1)
        end if


        CALL MPI_Bcast(InputFile, LEN(InputFile), MPI_CHARACTER, 0, MPI_COMM_WORLD, ierror) 

#endif
        file_unit = new_file_unit()
        call Ini_Open(InputFile, file_unit, bad, .false.)
        if (bad) call DoStop('Error opening parameter file')
        Ini_fail_on_not_found = .false.

        propose_scale = Ini_Read_Real('propose_scale',2.4)
        AccuracyLevel = Ini_Read_Real('accuracy_level',1.)

        checkpoint = Ini_Read_Logical('checkpoint',.false.)
        if (checkpoint) flush_write = .true.
        
#ifdef WMAP_PARAMS
        use_TT_beam_ptsrc = Ini_read_Logical('use_TT_beam_ptsrc')
        use_TE = Ini_read_Logical('use_TE')
        use_TT = Ini_Read_Logical('use_TT')
        print *, 'WMAP beam TE TT', use_TT_beam_ptsrc, use_TE, use_TT
#endif

#ifdef MPI 
        
        
        MPI_StartTime = MPI_WTime()
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
      
#endif

        Ini_fail_on_not_found = .true.
        
        baseroot = Ini_Read_String('file_root')
        
        rootname = trim(baseroot)

        FileChangeIniAll = trim(rootname)//'.read'

        if (instance /= 0) then
            rootname = trim(rootname)//'_'//trim(adjustl(numstr))
        end if

        new_chains = .true. 

        action = Ini_Read_Int('action',0)

        if (checkpoint) then
         if (action /= 0)  call DoStop('checkpoint only with action =0')
#ifdef MPI 
         new_chains = .not. FileExists(trim(rootname) //'.chk')
#else
         new_chains = .not. FileExists(trim(rootname) //'.txt')
#endif
        end if


        if (action ==1) call ReadPostParams

        FeedBack = Ini_Read_Int('feedback',0)
        FileChangeIni = trim(rootname)//'.read'

        if (action /= 1) then

        LogFile = trim(rootname)//'.log'

        if (LogFile /= '') then
         logfile_unit = 49
         call CreateOpenTxtFile(LogFile,logfile_unit,.not. new_chains)
        else
         logfile_unit = 0
        end if

        outfile_unit = 48
        fname = trim(rootname)//'.txt'
        if (new_chains) call CreateTxtFile(fname,outfile_unit)

        indep_sample = Ini_Read_Int('indep_sample')
        if (indep_sample /=0) then
          indepfile_unit = 47
          fname = trim(rootname)//'.data' 
          call CreateOpenFile(fname,indepfile_unit,'unformatted',.not. new_chains)
          !open(unit=indepfile_unit,file=fname,form='unformatted',status='replace')
        end if
 
        Ini_fail_on_not_found = .false.
        burn_in = Ini_Read_Int('burn_in',0)     
        sampling_method = Ini_Read_Int('sampling_method',sampling_metropolis)
        if (sampling_method > 6 .or. sampling_method<1) call DoStop('Unknown sampling method')
        if (sampling_method==4) directional_grid_steps = Ini_Read_Int('directional_grid_steps',20)
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

        w_is_w = Ini_Read_Logical ('w_is_w',.true.)     
        oversample_fast = Ini_Read_Int('oversample_fast',1)
        use_fast_slow = Ini_read_Logical('use_fast_slow',.true.)
 
        if (Ini_Read_Logical('cmb_hyperparameters', .false.)) &
            call DoStop( 'Hyperparameters not supported any more')

        if (Ini_Read_String('use_2dF') /= '') stop 'use_2dF now replaced with use_mpk'
        Use_Clusters = Ini_Read_Logical('use_clusters',.false.)
        Use_mpk = Ini_Read_Logical('use_mpk',.false.) ! matter power spectrum, incl 2dF
        Use_HST = Ini_Read_Logical('use_HST',.true.)
        Use_BBN = Ini_Read_Logical('use_BBN',.false.)
        Use_Age_Tophat_Prior= Ini_Read_Logical('use_Age_Tophat_Prior',.true.)
        Use_SN = Ini_Read_Logical('use_SN',.false.)
        Use_CMB = Ini_Read_Logical('use_CMB',.true.)
        Use_WeakLen = Ini_Read_Logical('use_WeakLen',.false.)
        Use_min_zre = Ini_Read_Double('use_min_zre',0.d0) 
        Use_Lya = Ini_Read_logical('use_lya',.false.)
        if (Use_Lya .and. use_nonlinear) &
             call DoStop('Lya.f90 assumes LINEAR power spectrum input')


        !flag to force getting sigma8 even if not using LSS data 
        use_LSS = Ini_Read_Logical('get_sigma8',.false.)
        ! use_LSS = Use_2dF .or. Use_Clusters .or. Use_WeakLen
        use_LSS = Use_LSS .or. Use_mpk .or. Use_Clusters .or. Use_WeakLen .or. Use_Lya

        Temperature = Ini_Read_Real('temperature',1.)
        
        num_threads = Ini_Read_Int('num_threads',0)
        !$ if (num_threads /=0) call OMP_SET_NUM_THREADS(num_threads)


        estimate_propose_matrix = Ini_Read_Logical('estimate_propose_matrix',.false.)
        if (estimate_propose_matrix) then
         if (Ini_Read_String('propose_matrix') /= '') &
           call DoStop('Cannot have estimate_propose_matrix and propose_matrix')
        end if

        delta_loglike = Ini_Read_Real('delta_loglike',2.)
        est_bfp_before_covmat = Ini_Read_Logical('est_bfp_before_covmat',.true.) ! for testing

        Ini_fail_on_not_found = .true.

        numsets = Ini_Read_Int('cmb_numdatasets')
        num_points = 0
        nuisance_params_used = 0
        if (Use_CMB) then
         do i= 1, numsets
          filename(i) = Ini_Read_String(numcat('cmb_dataset',i)) 
          call ReadDataset(filename(i))
          num_points = num_points + datasets(i)%num_points
          SZTemplate(i) = Ini_Read_String(numcat('cmb_dataset_SZ',i), .false.) 
          if (SZTemplate(i)/='') then
           SZScale(i) = Ini_read_Real(numcat('cmb_dataset_SZ_scale',i),1.0)
           call ReadSZTemplate(datasets(i), SZTemplate(i),SZScale(i))
          end if
          nuisance_params_used = nuisance_params_used + datasets(i)%nuisance_parameters
         end do
         if (Feedback > 1) write (*,*) 'read datasets'
        end if

        call Initialize(Params)

        Ini_fail_on_not_found = .true.
        
        nummpksets = Ini_Read_Int('mpk_numdatasets',0)
        do i= 1, nummpksets
         mpk_filename(i) = Ini_Read_String(numcat('mpk_dataset',i)) 
        end do

        numtoget = Ini_Read_Int('samples')

        call Ini_Close
        call ClearFileUnit(file_unit)

        call SetIdlePriority !If running on Windows


        if (Use_mpk) then
         do i=1,nummpksets
          call ReadMpkDataset(mpk_filename(i))
         end do
        if (Feedback>1) write(*,*) 'read mpk datasets'
        end if

        if (estimate_propose_matrix .or. action == 2) then
         ! slb5aug04  
            EstParams = Params
            EstParams%P(params_used) = Scales%center(params_used)
            if (est_bfp_before_covmat) then
              if (Feedback > 0) write(*,*) 'Finding max-like point' 
              call conjgrad_wrapper(EstParams,delta_loglike,status)  
              if (MPIRank == 0) then
                if (Feedback>0) write (*,*) 'Best fit parameters values:'
                call CreateTxtFile(trim(baseroot)//'.minimum',tmp_file_unit)
                write (tmp_file_unit,*) 'Best fit -log(Like) found: ', Bestfit_loglike
                write (tmp_file_unit,*) '' 
                do i=1, num_params_used
                  if (Feedback>0) write (*,*) params_used(i), ' : ', EstParams%P(params_used(i))
                  write (tmp_file_unit,*) params_used(i), ' : ', EstParams%P(params_used(i))
                end do
                close(tmp_file_unit)
              end if
            end if
              if (action == 2) then
                  if (Feedback>0) then
                    write(*,*) 'Have estimated the minimum, now exiting since action=2'
                    write(*,*) 'Wrote the minimum to file ',trim(baseroot)//'.minimum'
                  end if
                call DoStop
              end if
              if (Feedback>0) write (*,*) 'Now estimating propose matrix from Hessian'
              allocate(propose_matrix(num_params_used, num_params_used))
              propose_matrix=EstCovmat(EstParams,4.,status)
              ! By default the grid used to estimate the covariance matrix has spacings
              ! such that deltaloglike ~ 4 for each parameter.               
              call AcceptReject(.true., EstParams%Info, Params%Info)
              has_propose_matrix = status > 0
              if (MPIRank == 0) then ! write out the covariance matrix
                if (Feedback>0) write (*,*) 'Estimated covariance matrix:'
                call Matrix_write(trim(baseroot) //'.local_invhessian',propose_matrix,forcetable=.true.)
                if (Feedback>0) write(*,*) 'Wrote the local inv Hessian to file ',trim(baseroot)//'.local_hessian'
              end if
              if (has_propose_matrix) then
                   call SetProposeMatrix
              else
                   write (*,*) 'estimating propose matrix failed'
                   deallocate(propose_matrix)  
              end if
        end if
    
        if (action == 0) then
         if (Feedback > 0) write (*,*) 'starting Monte-Carlo'
         call MCMCSample(Params, numtoget)
 
         if (Feedback > 0) write (*,*) 'finished'

         if (logfile_unit /=0) close(logfile_unit)
         if (indepfile_unit /=0) close(indepfile_unit)
        
         close(outfile_unit)

        else if (action==1) then
          if (Feedback > 0) write (*,*) 'starting post processing'
          call postprocess(rootname)
        else
         call DoStop('undefined action')
        end if

        call DoStop

end program SolveCosmology

