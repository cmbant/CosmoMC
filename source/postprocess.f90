! This module does post-processing of .data files. For example to importance sample new data,
! correct approximate theory (eg. add in lensing), or to compute missing theory (e.g. matter power).


!March 2006: added redo_add parameter to add new data rather than re-computing all likelihoods
!Tensors now always recomputed if redo_cls=T

module posthoc
  use settings
  use cmbtypes
  use CMB_Cls
  use cmbdata
  use CalcLike
  implicit none

  Type TPostParams
     logical  redo_like, redo_theory, redo_cls, redo_Pk
     integer redo_skip, redo_thin
     character(LEN=120) :: redo_datafile, redo_outroot
     real redo_likeoffset 
     real redo_temperature
     logical redo_change_like_only 

       !This last one is for comparing goodness of fit
       !After importance sampling, you can recompute the likelihoods without the new data, but
       !keeping the weights from the importance sampling, and thereby asses whether the mean
       !likelihood wrt the original distribution of the parameter space after importance sampling 
       !is similar to that after, in which case the datasets intersect in a region of high likelihood

     logical redo_add 
       !if just want to add new datasets rather than re-computing the entire likelihood

     logical redo_from_text
       !Redo from text files if .data files not available

  end Type TPostParams

  Type(TPostParams) :: PostParams

  logical :: txt_theory = .false. !True to put P_k in output chains

contains

  subroutine ReadPostParams(baseroot)
    character(LEN=*), intent(in):: baseroot


    Ini_Fail_On_Not_Found = .false.
    PostParams%redo_like = Ini_Read_Logical('redo_likelihoods')
    PostParams%redo_theory = Ini_read_Logical('redo_theory')
    PostParams%redo_cls= Ini_read_Logical('redo_cls')
    PostParams%redo_Pk= Ini_read_Logical('redo_pk')
    if (Ini_read_Logical('redo_lensed', .false.)) then

      stop 'redo_lensed obsolete. Set CMB_lensing = T and redo_cls'

    end if

    PostParams%redo_skip = Ini_Read_Int('redo_skip',100)
    PostParams%redo_thin = max(1,Ini_Read_Int('redo_thin',1))
    PostParams%redo_datafile = Ini_Read_String('redo_datafile')
    PostParams%redo_outroot = Ini_Read_String('redo_outroot')
    PostParams%redo_likeoffset = Ini_Read_Double('redo_likeoffset',0.d0)
    PostParams%redo_temperature = Ini_Read_Double('redo_temp',1.d0)
    PostParams%redo_change_like_only = Ini_Read_Logical('redo_change_like_only',.false.)
    PostParams%redo_add = Ini_Read_Logical('redo_add',.false.)
    PostParams%redo_from_text = Ini_Read_Logical('redo_from_text',.false.)

    txt_theory = Ini_Read_Logical('txt_theory',.false.)

    if (PostParams%redo_outroot == '') then
      PostParams%redo_outroot = trim(ExtractFilePath(baseroot))//'post_' &
               // trim(ExtractFileName(baseroot))
    end if

    
  end subroutine ReadPostParams

   subroutine postprocess(InputFile)
        use IO
        character(LEN=*), intent(INOUT):: InputFile
        Type(CMBParams) LastCMB,CMB, newCMB
        Type(CosmoTheory) Theory, CorrectTheory
        Type(ParamSetInfo) Info
        real truelike,mult,like
        real weight_min, weight_max, mult_sum, mult_ratio, mult_max,weight

        real max_like, max_truelike
        integer error,num, debug
        character (LEN=120) :: post_root
        real Params(num_params)
        integer infile_handle, outdata_handle

        flush_write = .false.
        weight_min= 1e30
        weight_max = -1e30
        mult_sum = 0
        mult_ratio = 0
        mult_max = -1e30
        max_like = 1e30

        max_truelike =1e30

        debug = 0

        infile_handle = 0
        Temperature = PostParams%redo_temperature 

        if (Feedback>0 .and. PostParams%redo_change_like_only) &
              write (*,*) 'Warning: only changing likelihoods not weights'

        if (PostParams%redo_datafile /= '') InputFile = PostParams%redo_datafile

        if (PostParams%redo_from_text) then

         infile_handle = IO_OpenChainForRead(trim(InputFile)//'.txt')
         
         if (.not. PostParams%redo_theory) write (*,*) '**You probably want to set redo_theory**'
!         if (.not. PostParams%redo_cls .and. Use_CMB) write (*,*) '**You probably want to set redo_cls**'
!         if (.not. PostParams%redo_pk .and. Use_mpk) write (*,*) '**You probably want to set redo_pk**'

         if (PostParams%redo_thin>1) write (*,*) 'redo_thin only OK with redo_from_text if input weights are 1'

        else
         infile_handle = IO_OpenDataForRead(trim(InputFile)//'.data')
        end if

        post_root = PostParams%redo_outroot 
    
        if (MpiRank==0 .and. NameMapping%nnames/=0) &
            call IO_OutputParamNames(NameMapping,trim(post_root))
        
        if (instance /= 0) post_root = numcat(trim(post_root)//'_',instance)

        if (Feedback > 0) then
          if (PostParams%redo_from_text) then
           write (*,*) 'reading from: ' //  trim(InputFile)//'.txt'
          else
           write (*,*) 'reading from: ' //  trim(InputFile)//'.data'
          end if
           write (*,*) 'writing to: ' // trim(post_root)//'.*'
        end if

        write (*,*) 'Using temperature: ', Temperature

        outfile_handle = IO_OutputOpenForWrite(trim(post_root)//'.txt')
        
        outdata_handle = IO_DataOpenForWrite(trim(post_root)//'.data')
     
        num = 0
        do

        if (PostParams%redo_from_text) then
          error = 0
          if (.not. IO_ReadChainRow(infile_handle, mult, like, Params, num_params)) exit
          call ParamsToCMBParams(Params, CMB)

        else
          call ReadModel(infile_handle,CMB,Theory,mult,like, error)
        end if

        if (error ==1) then
          if (num==0) call MpiStop('Error reading data file.')

          exit
        end if
        num=num+1
        if (num>PostParams%redo_skip .and. mod(num,PostParams%redo_thin) == 0) then
     
        LastCMB = CMB

        newCMB = CMB  
        if (PostParams%redo_theory) then
    
         call GetClsInfo(newCMB, CorrectTheory, error, PostParams%redo_cls, PostParams%redo_pk)
    
         if (PostParams%redo_cls) then
            Theory%cl = CorrectTheory%cl
            Theory%cl_tensor = CorrectTheory%cl_tensor
         end if
  
         if (PostParams%redo_pk) then
            Theory%sigma_8 = CorrectTheory%sigma_8
            Theory%Matter_Power = CorrectTheory%Matter_Power
         end if

         Theory%derived_parameters = CorrectTheory%derived_parameters
         Theory%numderived = CorrectTheory%numderived
     
        else
          error = 0
        end if


          CorrectTheory = Theory


         if (error ==0) then
 
           if (PostParams%redo_like) then

              if (Use_LSS .and. CorrectTheory%sigma_8==0) &
                 call MpiStop('Matter power/sigma_8 have not been computed. Use redo_theory and redo_pk.')

              Info%Theory = CorrectTheory
              call CMBParamsToParams(newCMB, Params)
              truelike = GetLogLikePost(newCMB, Params, Info)

              if (truelike == logZero) then
               weight = 0
              else
               if (PostParams%redo_add) truelike = like + truelike

               weight = exp(like-truelike+PostParams%redo_likeoffset)
              end if
           
             if (.not. PostParams%redo_change_like_only)  mult = mult*weight
           else
              truelike = like
              weight = 1
           end if


           max_like = min(max_like,like)

           max_truelike = min(max_truelike,truelike)

           mult_ratio = mult_ratio + weight
           mult_sum = mult_sum + mult
          
           if (mult /= 0) then        
            call WriteCMBParams(newCMB, CorrectTheory, mult, truelike,txt_theory)
            call WriteModel(outdata_handle, newCMB, CorrectTheory,Info%likelihoods,truelike,mult)

           else 

            if (Feedback >1 ) write (*,*) 'Zero weight: new like = ', truelike
           end if
          
           if (Feedback > 1) write (*,*) 'done ', num, 'mult= ', mult,' weight = ', weight
           weight_max = max(weight,weight_max)       
           weight_min = min(weight,weight_min)
           mult_max = max(mult_max,mult)

        end if
        end if
        
        end do

        call IO_Close(infile_handle)
        call IO_Close(outfile_handle)
        call IO_DataCloseWrite(outdata_handle)
      
        num = (num - PostParams%redo_skip) / PostParams%redo_thin
        if (Feedback>0) then 
           write(*,*) 'finished. Processed ',num,' models'
           write (*,*) 'max weight= ',weight_max, ' min weight = ',weight_min
           write (*,*) 'mean mult  = ', mult_sum/num
           write (*,*) 'mean importance weight (approx evidence ratio) = ',mult_ratio/num
           write (*,*) 'effective number of samples =',mult_sum/mult_max
           write (*,*) 'Best redo_likeoffset = ',max_truelike - max_like
         end if
        
        if ((mult_ratio < 1e-6 .or. mult_ratio > 1e8) & 

            .and. .not.PostParams%redo_change_like_only) then

         write (*,*) 'WARNING: use redo_likeoffset to rescale likelihoods'
        end if



        return

        end subroutine postprocess

end module posthoc
