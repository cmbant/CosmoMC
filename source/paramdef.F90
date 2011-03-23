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
 implicit none

 Type ParamSet
    real P(num_params)
    Type(ParamSetInfo) Info
 end Type

 Type ParamScale
    real PMin(num_params), PMax(num_params), PWidth(num_params), center(num_params)
 end Type ParamScale

 Type(ParamScale) Scales

 Real :: StartLike = LogZero
   !bad, unless re-starting in which case it is set to last value in chain
 integer :: Num = 0
   !zero unless from checkpoint
 integer :: burn_in = 2 !Minimum to get sensible answers

 logical :: has_propose_matrix = .false., estimate_propose_matrix = .false.
 real, dimension(:,:), allocatable :: propose_matrix, propose_matrix_fast 
 real, dimension(:),   allocatable :: sigmas, propose_diag, propose_diag_fast

 integer, dimension(:), allocatable :: slow_evecs

 logical :: checkpoint = .false.
 logical :: new_chains = .true.
 integer, parameter :: checkpoint_freq = 100
 character(LEN=500) :: rootname
 real pmat(num_params,num_params) !actual covariance matrix

 real    :: MPI_R_Stop = 0.05

 integer  :: MPI_thin_fac = 3
 !following are numbers of samples / (number of parameters)
 !MPI_Min_Sample_Update*num_fast is min number of chain steps
 !before starting to update covmat and checking convergence
 integer :: MPI_Min_Sample_Update = 200, MPI_Sample_update_freq = 50


 logical :: MPI_Check_Limit_Converge = .false.
!After given R_Stop is reached, can optionally check for limit convergence
!(which is generally a much more stringent test of enough samples)
 real    :: MPI_Limit_Converge = 0.025
 real    :: MPI_Limit_Converge_Err = 0.3 
   !Tolerated cross-chain error on the limit in units of standard deviation
 integer :: MPI_Limit_Param = 0
   !Which parameter's limits to check. If zero, do them all
 logical :: MPI_LearnPropose = .true.
 
!If convergence R is too bad, we don't update proposal matrix if we had one to start with
!(which should be quite good unless using new parameters or much better data)
 real    :: MPI_Max_R_ProposeUpdate = 2., MPI_Max_R_ProposeUpdateNew  = 30.

 real    :: MPI_R_StopProposeUpdate = 0.

 logical    :: MPI_StartSliceSampling = .false.

 double precision    :: MPI_StartTime

 real, private, allocatable, dimension(:,:) :: MPICovmat


contains

subroutine DoStop(S)
 character(LEN=*), intent(in), optional :: S
 integer ierror
        if (present(S)) write (*,*) trim(S)
#ifdef MPI 
        MPI_StartTime = MPI_WTime() - MPI_StartTime 
        if (Feedback > 0 .and. MPIRank==0) then

           write (*,*) 'Total time:', nint(MPI_StartTime), &
                                   '(',MPI_StartTime/(60*60),' hours)'

           write (*,*) 'Slow proposals: ', slow_proposals

        end if
        call DoAbort('') !Abort all in case other continuing chains want to communicate with us
        !call mpi_finalize(ierror)
#endif

#ifdef DECONLY
       pause
#endif
       stop
end subroutine DoStop


subroutine DoAbort(S)

 character(LEN=*), intent(in), optional :: S
#ifdef MPI 
 integer ierror
#endif
        if (present(S)) write (*,*) trim(S)
#ifdef MPI 
        call MPI_Abort(MPI_COMM_WORLD,ierror,ierror)
#endif

#ifdef DECONLY
       pause
#endif
       stop
end subroutine DoAbort

subroutine Initialize(Ini,Params)
        use IniFile
        use ParamNames
        use settings
        use IO

        implicit none
        type (ParamSet) Params
        type (TIniFile) :: Ini
        integer i
        character(LEN=5000) fname,InLine
        character(LEN=1024) prop_mat
        real center, wid, mult, like
        
        output_lines = 0
        
        call SetParamNames(NameMapping)
        
        Ini_fail_on_not_found = .false.

        call CMB_Initialize(Params%Info)

        prop_mat = trim(Ini_Read_String_File(Ini,'propose_matrix'))
        has_propose_matrix = prop_mat /= ''
        prop_mat = concat(LocalDir,prop_mat)
         
        fname = Ini_Read_String_File(Ini,'continue_from')
        if (fname /= '') call DoStop('continue_from replaced by checkpoint')

       if (.not. new_chains) then
            fname = trim(rootname) //'.txt'
            call IO_ReadLastChainParams(fname, mult, like, Params%P, num_params)
            StartLike = like
            burn_in = 4
            outfile_handle = IO_OutputOpenForWrite(fname, append = .true.)
         end if
    
        num_params_used = 0
        num_fast = 0
        num_slow = 0
        do i=1,num_params
   
           if (i>=index_nuisance) then
            !All unit Gaussians
            center=0
            Scales%PMin(i)=-10
            Scales%PMax(i)=10
            wid=1
            if (i-index_nuisance < nuisance_params_used) then
              Scales%PWidth(i)=1
            else
              Scales%PWidth(i)=0
            end if  
           else 
            InLine =  ParamNames_ReadIniForParam(NameMapping,DefIni,'param',i)
            if (InLine=='') then
             if (NameMapping%nnames/=0) then
              call MpiStop('parameter ranges not found, param '//trim(IntToStr(i))// &
                             ' - '//trim(NameMapping%Name(i)))
             else
              call MpiStop('parameter ranges not found, param '//trim(IntToStr(i)))
             end if
             !Ini_Read_String(numcat('param',i), .true.)
            end if
            read(InLine, *, err = 100) center, Scales%PMin(i), Scales%PMax(i), wid, Scales%PWidth(i)
           end if  
           Scales%center(i) = center
           if (Scales%PMax(i) < Scales%PMin(i)) call DoStop('You have param Max < Min')
           if (Scales%PWidth(i) /= 0) then
               num_params_used = num_params_used + 1
              if (i > num_hard .and. use_fast_slow) then
                num_fast = num_fast + 1
              else
                num_slow = num_slow +1
              end if
           end if
           if (new_chains) then
           do
            if (wid < 0) then
              !This case we want half gaussian, width -wid
              !e.g. for positive definite parameters
              Params%P(i) = center - abs(Gaussian1())*wid
            else
              Params%P(i) = center + Gaussian1()*wid
            end if
            !Repeat until get acceptable values in range
            if (Params%P(i)>=  Scales%PMin(i) .and. Params%P(i) <= Scales%PMax(i)) exit
      
           end do
           end if
        end do
  
        if (Feedback > 0 ) write(*,'(" Varying ",1I2," parameters (",1I2," fast)")') &
               num_params_used,num_fast
        
        allocate(params_used(num_params_used))
        allocate(fast_params_used(num_fast))

        num_params_used = 0
        num_fast = 0
        do i=1,num_params
           if (Scales%PWidth(i) /= 0) then
              num_params_used = num_params_used + 1
              params_used(num_params_used) = i
              if (i > num_hard .and. use_fast_slow) then
                num_fast = num_fast + 1
                fast_params_used(num_fast) = i          
              end if
           end if
        end do
   
        if (has_propose_matrix) then
           
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

           allocate(propose_matrix(num_params_used, num_params_used))
           propose_matrix = pmat(params_used, params_used)

           call SetProposeMatrix
          
        end if

        if (.not. new_chains) call AddMPIParams(Params%P,like, .true.)    
   
        return
100     call DoStop('Error reading param details: '//trim(InLIne))

end subroutine Initialize



subroutine SetProposeMatrix
  use MAtrixUtils
  real proj_len(num_params_used)
  real U(num_params_used,num_params_used)
  real vecs(num_slow,num_params_used)
  integer i, ii,j


   pmat(params_used, params_used) = propose_matrix

!std devs of base parameters

   if (.not. allocated(sigmas)) allocate(sigmas(num_params_used))
   do i = 1, num_params_used
     sigmas(i) = sqrt(propose_matrix(i,i))
   end do

   do i = 1, num_params_used
     propose_matrix(i,:) = propose_matrix(i,:) / sigmas(i)
     propose_matrix(:,i) = propose_matrix(:,i) / sigmas(i)
   end do



   if (num_fast /= 0) then

      !Get the conditional covariance by projecting the inverse covariances of fast parameters
            if (.not. allocated(propose_matrix_fast)) then  
              allocate(propose_matrix_fast(num_fast, num_fast))
              allocate(propose_diag_fast(num_fast))
            end if

            U = propose_matrix
            call Matrix_Inverse(U)
            propose_matrix_fast = U(num_slow+1:num_params_used, num_slow+1:num_params_used)
            call Matrix_Diagonalize(propose_matrix_fast,propose_diag_fast, num_fast)
            !propose_matrix^-1 = U D U^T, returning U in propose_matrix
            propose_matrix_fast = transpose(propose_matrix_fast)

            if (any(propose_diag_fast <= 0)) &
              call DoStop('Fast proposal matrix has negative or zero eigenvalues')


            propose_diag_fast = 1/sqrt(propose_diag_fast)

   end if



   if (num_slow /= 0) then

    if (.not. allocated(propose_diag)) allocate(propose_diag(num_params_used))

    call Matrix_Diagonalize(propose_matrix, propose_diag, num_params_used)
      !propose_matrix = U D U^T, returning U in propose_matrix


    if (any(propose_diag <= 0)) &
        call DoStop('Proposal matrix has negative or zero eigenvalues')
    propose_diag = sqrt(max(1e-12,propose_diag))


!Get projected lengths 

    do i = 1, num_params_used
          vecs(:,i) = propose_diag(i)*propose_matrix(1:num_slow,i)
          proj_len(i) = sum(vecs(:,i)**2)  
    end do



     if (.not. allocated(slow_evecs)) allocate(slow_evecs(num_slow))

!keep evectors with longest projected lengths in slow dimensions, orthogonal to previous longest       

       do i = 1, num_slow
         j = MaxIndex(proj_len, num_params_used)
         slow_evecs(i) = j
         do ii= 1, num_params_used
           if (proj_len(ii) /= 0. .and. ii/=j) then
            vecs(:,ii) = vecs(:,ii) - sum(vecs(:,j)*vecs(:,ii))*vecs(:,j)/proj_len(j)
                  !Take out projection onto jth eigendirection
            proj_len(ii) = sum(vecs(:,ii)**2)
           end if
         end do

         proj_len(j) = 0.
       end do
   end if

end subroutine SetProposeMatrix



 
  subroutine WriteCMBParams(CMB,Theory,mult,like, with_data)     
     Type (CosmoTheory) Theory 
     real, intent(in) :: mult, like
     logical, intent(in), optional :: with_data
     Type(ParamSet) P
     Type(CMBParams) CMB


     call CMBParamsToParams(CMB,P%P)

     P%Info%Theory = Theory
     if (present(with_data)) then
      if (with_data) then
      call WriteParamsAndDat(P, mult,like)
      else
      call WriteParams(P, mult,like)
      end if
     else
      call WriteParams(P, mult,like)
     end if

  end subroutine WriteCMBParams
   
  subroutine WriteIndepSample(P, like)
    Type(ParamSet) P
    real like
    Type(CMBParams) C
    if (indepfile_handle ==0) return
    call ParamsToCMBParams(P%P,C)
    call WriteModel(indepfile_handle, C,P%Info%Theory,like)
   end subroutine WriteIndepSample


   subroutine AddMPIParams(P,like, checkpoint_start)     
     real, intent(in) ::like
     real P(:)
     logical, intent(in), optional :: checkpoint_start
!Collect thinned samples after a burn-in perdiod
!Then use second half of the samples to get convergence
!Use R = worst eigenvalue (variance of chain means)/(mean of chain variances) statistic for
!convergence test, followed optionally by (variance of limit)/(mean variance) statistic
!If MPI_LearnPropose then update proposal density using covariance matrix of last half of chains 
#ifdef MPI

     integer, save :: sample_num = 0
     integer i,j,k 
     logical, save :: Burn_done = .false.
     integer, save :: npoints = 0
     integer index, STATUS(MPI_STATUS_SIZE),STATs(MPI_STATUS_SIZE*(MPIChains-1))
     logical flag

 
     real, allocatable, dimension(:), save ::MPIMean
     real, allocatable, dimension(:,:,:) ::MPICovmats
     real, allocatable, dimension(:,:)   ::MPIMeans
     real norm, mean(num_params_used), chain_means(MPIChains,num_params_used)
     integer, parameter :: chk_id = 3252356
     integer ID
     real MeansCov(num_params_used,num_params_used), cov(num_params_used,num_params_used)
     real sc, evals(num_params_used), last_P, R
     integer ierror
     logical DoCheckpoint

     integer, allocatable, dimension(:), save :: req, buf 
     integer,  allocatable, dimension(:), save :: param_changes
     logical, save :: StartCovMat, flukecheck, Waiting = .false.

     Type(TList_RealArr), save :: S
     integer, save :: slice_fac = 1
     logical, save :: all_burn = .false., done_check = .false., DoUpdates = .false.
     character(LEN=128) logLine


!Read in checkpoing stuff at restart
     if (checkpoint .and. present(checkpoint_start)) then
      if (Feedback > 0) write (*,*) instance, 'Reading checkpoint from '//trim(rootname)//'.chk'
      call OpenFile(trim(rootname)//'.chk',tmp_file_unit,'unformatted')
      read (tmp_file_unit) ID
      if (ID/=chk_id) call DoStop('invalid checkpoint files')
      read(tmp_file_unit) num, sample_num, MPI_thin_fac, npoints, Burn_done, all_burn,sampling_method, &
            slice_fac, has_propose_matrix, S%Count, flukecheck, StartCovMat, MPI_Min_Sample_Update, DoUpdates
      if (has_propose_matrix) then
        if (.not. allocated(propose_matrix)) then
              allocate(propose_matrix(num_params_used, num_params_used))
        end if
        read(tmp_file_unit) propose_matrix  
        call SetProposeMatrix
      end if
      call TList_RealArr_Init(S)
      call TList_RealArr_ReadBinary(S,tmp_file_unit)
      close(tmp_file_unit)
      allocate(req(MPIChains-1))      
      allocate(MPIcovmat(num_params_used,num_params_used))
      allocate(MPIMean(0:num_params_used))
      return  
     end if

!Dump checkpoint info
!Have to be careful if were to dump before burn
     if (checkpoint .and. all_burn .and. (.not. done_check .or. &
            mod(sample_num+1, checkpoint_freq)==0)) then
      done_check=.true.
      if (Feedback > 1) write (*,*) instance, 'Writing checkpoint'
      call CreateFile(trim(rootname)//'.chk',tmp_file_unit,'unformatted')
      write (tmp_file_unit) chk_id
      write(tmp_file_unit) num, sample_num, MPI_thin_fac, npoints, Burn_done, all_burn, sampling_method, &
            slice_fac, has_propose_matrix, S%Count, flukecheck, StartCovMat, MPI_Min_Sample_Update, DoUpdates
      if (has_propose_matrix) write(tmp_file_unit) pmat(params_used, params_used) 
      call TList_RealArr_SaveBinary(S,tmp_file_unit)
      close(tmp_file_unit)      
     end if


!Do main adding samples functions
     sample_num = sample_num + 1
     if (mod(sample_num, MPI_thin_fac) == 0) return 


     if (npoints == 0 .and. .not. Burn_done)then
        allocate(param_changes(num_params_used))
        param_changes= 0
        call TList_RealArr_Init(S)
        if (MPI_StartSliceSampling) then
          sampling_method = sampling_slice
        end if 

        if (sampling_method == sampling_slice) then
          if (Feedback > 0  .and. .not. Burn_done) write (*,*) 'Starting with slice sampling'
          slice_fac = 4
        end if     
     end if
 

     call TList_RealArr_Add(S, P(params_used))
     npoints = npoints + 1       


     if (.not. Burn_done) then
         
       if (npoints > 200/MPI_thin_fac +1) then
        !We're not really after independent samples or all of burn in
        !Make sure all parameters are being explored
             do i=1, num_params_used             
              if (S%Items(S%Count)%P(i) /= S%Items(S%Count-1)%P(i)) &
                   param_changes(i) =  param_changes(i) + 1
             end do
             Burn_done = all(param_changes > 100/MPI_thin_fac/slice_fac+2)
           if (Burn_done) then
               if (Feedback > 0) then
                  write (*,*) 'Chain',instance, &
                  ' MPI done ''burn'', like = ',like, 'Samples = ',sample_num       
                  write (*,*) 'Time: ', MPI_WTime() - MPI_StartTime, 'output lines=',output_lines 
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


               call TList_RealArr_Clear(S)         
               MPI_Min_Sample_Update = &
                  max((MPI_Min_Sample_Update*max(1,num_fast))/(MPI_thin_fac*slice_fac), npoints)
               npoints = 0
               StartCovMat = has_propose_matrix
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


        if (.not. DoUpdates  .and. all_burn .and. npoints >= MPI_Min_Sample_Update+50/MPI_thin_fac/slice_fac+1) then
             DoUpdates = .true.
             if (Feedback>0) write(*,*) instance, 'DoUpdates'
        end if

        if (DoUpdates) then
            if (MPIRank == 0) then
 
              if (Waiting) then  
                call MPI_TESTALL(MPIChains-1,req, flag, stats, ierror)
                Waiting = .not. flag

              elseif (mod(npoints,max(1,(MPI_Sample_update_freq*num_params_used)/(slice_fac*MPI_thin_fac)))==0) then

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
              MPiMean(1:num_params_used) = MPiMean(1:num_params_used) + S%Items(i)%P     
             end do
             MPiMean(1:num_params_used) = MPiMean(1:num_params_used) / MPImean(0)
             do i = S%Count/2, S%Count
               do j = 1, num_params_used
                 MPICovmat(:,j) =  MPICovmat(:,j) + &
                 (S%Items(i)%P(:)-MPIMean(1:num_params_used))*(S%Items(i)%P(j)- MPIMean(j))
              end do
             end do
             MPICovMat = MPICovMat / MPImean(0)
 
            MPICovMats(:,:,instance) = MPICovMat
            MPIMeans(:,instance) = MPIMean

            do i=1, MPIChains  
               j = i-1
               call MPI_BCAST(MPICovMats(:,:,i),Size(MPICovMat),MPI_REAL,j,MPI_COMM_WORLD,ierror)
               call MPI_BCAST(MPIMeans(:,i),Size(MPIMean),MPI_REAL,j,MPI_COMM_WORLD,ierror)
             end do

             

  ! These should be better but don't work
  !          call MPI_ALLGATHER(MPICovMat,Size(MPICovMat),MPI_REAL,MPICovmats,Size(MPICovmats), &
  !               MPI_REAL, MPI_COMM_WORLD,ierror)
  !          call MPI_ALLGATHER(MPIMean,Size(MPIMean),MPI_REAL,MPIMeans,Size(MPIMeans), &
  !               MPI_REAL, MPI_COMM_WORLD,ierror)
 
            if (all(MPIMeans(0,:)> MPI_Min_Sample_Update/2 + 2)) then
               !check have reasonable number of samples in each)   
                  norm = sum(MPIMeans(0,:))

                  do i=1, num_params_used
                   mean(i) = sum(MPIMeans(i,:)*MPIMeans(0,:))/norm
                   chain_means(:,i) = MPIMeans(i,:) 
                  end do 
                  
                  if (MPIChains > 1) then
                  !Do generalized Gelman & Rubin to assess convergence
                      do i=1,num_params_used
                        do j=i,num_params_used
                         cov(i,j) = sum(MPIMeans(0,:)*MPICovMats(i,j,:))/ norm
                         meanscov(i,j) = sum(MPIMeans(0,:)*&
                          (chain_means(:,i)-mean(i))*(chain_means(:,j)-mean(j)))/norm 
                         meanscov(j,i) = meanscov(i,j)
                         cov(j,i) = cov(i,j)
                        end do
                      end do
                      MPICovMat = Cov + meansCov !Estimate global covariance for proposal density   
                      meansCov = meansCov * real(MPIChains)/(MPIChains-1)
                      do j=1,num_params_used
                            sc = sqrt(cov(j,j))
                            cov(j,:) = cov(j,:) / sc
                            cov(:,j) = cov(:,j) / sc
                            meanscov(j,:) = meanscov(j,:) /sc
                            meanscov(:,j) = meanscov(:,j) /sc                 
                      end do
 
                     call Matrix_Diagonalize(meanscov, evals, num_params_used)
                     cov = matmul(matmul(transpose(meanscov),cov),meanscov)
                     R = 0
                     do j=1,num_params_used
                       R = max(R,evals(j)/max(1e-12,cov(j,j)))
                     end do
                     if (Feedback > 0 .and. MPIRank==0) &
                      write (*,*) 'Current convergence R-1 = ',R, 'chain steps =',sample_num
                     if (logfile_unit/=0) then
                      write(logLine,*) 'Current convergence R-1 = ',R, 'chain steps =',sample_num
                      call IO_WriteLog(logfile_unit,logLine)
                     end if 
                     if (R < MPI_R_Stop .and. flukecheck) then
                       if (MPI_Check_Limit_Converge) then
                        !Now check if limits from different chains agree well enough
                        call CheckLImitsConverge(S)
                       else
                        !If not also checking limits, we are done
                        call DoStop
                       end if
                     end if
                     flukecheck = R < MPI_R_Stop
                     if (S%Count > 100000) then
                        !Try not to blow memory by storing too many samples
                          call TList_RealArr_Thin(S, 2)
                          MPI_thin_fac = MPI_thin_fac*2 
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
              

               if (.not. allocated(propose_matrix)) then
                 allocate(propose_matrix(num_params_used, num_params_used))
               end if

               propose_matrix = MPICovMat
        

               call SetProposeMatrix

               has_propose_matrix = .true.

             end if !update propose
            end if !all chains have enough

            deallocate(MPICovMats, MPIMeans)
           end if !flag 
            
          
    end if

#endif
     

 end subroutine AddMPIParams

#ifdef MPI

 subroutine CheckLimitsConverge(L)
  !Check limits from last half chains agree well enough across chains to be confident of result
  !Slowly explored tails will cause problems (long time till stops)
  Type(TList_RealArr), intent(in) :: L
  integer i,j, side, ierror, worsti
  real, allocatable, dimension(:,:,:) :: Limits
  real MeanLimit, var, LimErr, WorstErr
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
      call TList_RealArr_ConfidVal(L, params_check(j), MPI_Limit_Converge, L%Count/2, L%Count, &
                Limits(1,j,instance),Limits(2,j,instance))
     end do
      !Now tell everyone else
     do i=1, MPIChains  
         j = i-1
         call MPI_BCAST(Limits(:,:,i),2*numCheck,MPI_REAL,j,MPI_COMM_WORLD,ierror)
     end do
!Take as test statistics the rms deviation from the mean limit in units of the standard deviation
     WorstErr = 0.
     do j=1, numCheck
      do side = 1,2
       MeanLimit = Sum(Limits(side,j,:))/MPIChains
       var = sum((Limits(side,j,:) - MeanLimit)**2)/(MPIChains-1)
       LimErr = sqrt(var / MPICovMat(params_check(j),params_check(j))) 
       if(LimErr > WorstErr) then
        WorstErr = LimErr
        Worsti = params_check(j)
       end if
       if (MPIRank ==0) print *, 'param ',params_check(j), 'lim err', LimErr
      end do
     end do
   if (Feedback > 0 .and. MPIRank==0) then
     write (*,*) 'Current worst limit error = ', WorstErr
     write (*,*) 'for parameter ',Worsti, 'samps = ',L%Count*MPI_thin_fac   
   end if
   if (logfile_unit/=0) then
      write (logLine,*) 'Current limit err = ',WorstErr, ' param ',Worsti, 'samps = ',L%Count*MPI_thin_fac
      call IO_WriteLog(logfile_unit,logLine)
   end if
   if (WorstErr < MPI_Limit_Converge_Err) call DoStop
     
   deallocate(Limits)       
   deallocate(params_check)
 end subroutine CheckLimitsConverge

#endif

end module ParamDef



 