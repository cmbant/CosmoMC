    ! This module does post-processing of .data files. For example to importance sample new data,
    ! correct approximate theory (eg. add in lensing), or to compute missing theory (e.g. matter power).

    module posthoc
    use settings
    use CMB_Cls
    use CalcLike
    use powerspec
    implicit none

    Type TPostParams
        logical  redo_like, redo_theory, redo_cls, redo_pk
        integer redo_thin
        real(mcp) redo_skip
        character(LEN=Ini_max_string_len) :: redo_datafile, redo_outroot
        real(mcp) redo_likeoffset
        real(mcp) redo_temperature
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

        logical redo_no_new_data !true to make no new .data files to save space
        character(LEN=128) :: redo_like_name
    end Type TPostParams

    Type(TPostParams) :: PostParams

    !not supported any more    logical :: txt_theory = .false. !True to put P_k in output chains

    contains

    subroutine ReadPostParams(baseroot)
    character(LEN=*), intent(in):: baseroot

    Ini_Fail_On_Not_Found = .false.
    PostParams%redo_like = Ini%Read_Logical('redo_likelihoods')
    PostParams%redo_theory = Ini%read_Logical('redo_theory')
    PostParams%redo_cls= Ini%read_Logical('redo_cls')
    PostParams%redo_pk= Ini%read_Logical('redo_pk')
    PostParams%redo_skip = Ini%Read_Double('redo_skip',100.d0)
    PostParams%redo_thin = max(1,Ini%Read_Int('redo_thin',1))
    PostParams%redo_datafile = Ini%Read_String('redo_datafile')
    PostParams%redo_outroot = Ini%Read_String('redo_outroot')
    PostParams%redo_likeoffset = Ini%Read_Double('redo_likeoffset',0.d0)
    PostParams%redo_temperature = Ini%Read_Double('redo_temp',1.d0)
    PostParams%redo_change_like_only = Ini%Read_Logical('redo_change_like_only',.false.)
    PostParams%redo_add = Ini%Read_Logical('redo_add',.false.)
    PostParams%redo_from_text = Ini%Read_Logical('redo_from_text',.false.)
    PostParams%redo_no_new_data = Ini%Read_Logical('redo_no_new_data',.false.)
    PostParams%redo_like_name = Ini%Read_String('redo_like_name')

    if (PostParams%redo_from_text .and. (PostParams%redo_add .or. PostParams%redo_like_name/='')) &
    call Mpistop('redo_new_likes requires .data files, not from text')

    if (PostParams%redo_from_text  .and. PostParams%redo_skip>0.d0 .and. PostParams%redo_skip<1) &
    call Mpistop('redo_from_text currently requires redo_skip==0 or redo_skip>=1')

    !    txt_theory = Ini%Read_Logical('txt_theory',.false.)

    if (PostParams%redo_outroot == '') then
        PostParams%redo_outroot = trim(ExtractFilePath(baseroot))//'post_' &
        // trim(ExtractFileName(baseroot))
    end if

    end subroutine ReadPostParams

    subroutine postprocess(InputFile)
    use IO
    USE IFPOSIX
    character(LEN=*), intent(INOUT):: InputFile
    Type(TheoryPredictions) newTheory
    real(mcp) truelike,mult,like
    real(mcp) weight_min, weight_max, mult_sum, mult_ratio, mult_max,weight
    real(mcp) max_like, max_truelike
    integer error,num, debug
    character (LEN=Ini_max_string_len) :: post_root
    integer i, infile_handle
    integer :: outdata_handle=-1
    Type (ParamSet) :: Params
    logical :: has_likes(DataLikelihoods%Count)
    Type(DataLikelihood), pointer :: DataLike
    logical :: first = .false., has_chain = .true.
    integer last_file_loc,file_loc, file_size
    integer :: at_beginning=0, ierror, num_used
    integer :: numz, index_error

    flush_write = .false.
    weight_min= 1e30_mcp
    weight_max = -1e30_mcp
    mult_sum = 0
    mult_ratio = 0
    mult_max = -1e30_mcp
    max_like = 1e30_mcp

    max_truelike =1e30_mcp

    debug = 0

    infile_handle = 0
    Temperature = PostParams%redo_temperature

    if (Feedback>0 .and. PostParams%redo_change_like_only) &
    write (*,*) 'Warning: only changing likelihoods not weights'

    if (PostParams%redo_datafile /= '') InputFile = PostParams%redo_datafile

    if (PostParams%redo_from_text) then
        infile_handle = IO_OpenChainForRead(trim(InputFile)//'.txt')
        if (.not. PostParams%redo_theory) write (*,*) '**You probably want to set redo_theory**'
        if (PostParams%redo_thin>1) write (*,*) 'redo_thin only OK with redo_from_text if input weights are 1'
    else
        if (FileExists(trim(InputFile)//'.data')) then
            infile_handle = IO_OpenDataForRead(trim(InputFile)//'.data')
        else
            write(*,*) 'Chain .data files does not exist: ', MpiRank+1
            has_chain =.false.
        end if
    end if

    post_root = PostParams%redo_outroot

    if (MpiRank==0 .and. NameMapping%nnames/=0) then
        call IO_OutputParamNames(NameMapping,trim(post_root),params_used, add_derived=.true.)
        call BaseParams%OutputParamRanges(NameMapping, trim(post_root)//'.ranges')
    end if

    if (has_chain) then
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
        if (.not. PostParams%redo_no_new_data) outdata_handle = IO_DataOpenForWrite(trim(post_root)//'.data')
        num = 0
        num_used = 0

        do
            if (PostParams%redo_from_text) then
                error = 0
                Params%P= Scales%center
                if (.not. IO_ReadChainRow(infile_handle, mult, like, Params%P, params_used)) exit
                num=num+1
            else
                call Params%ReadModel(infile_handle,has_likes, mult,like, error)
                num=num+1
                if (first .and. PostParams%redo_like_name/='') then
                    first=.false.
                    do i=1, DataLikelihoods%Count
                        DataLike => DataLikelihoods%Item(i)
                        if (DataLike%name==PostParams%redo_like_name) then
                            if (.not. has_likes(i)) &
                            call MpiStop('does not currently have like named:'//trim(PostParams%redo_like_name))
                            has_likes(i)=.true.
                            if (any(.not. has_likes)) call MpiStop('not all other likelihoods exist already')
                            has_likes(i)=.false.
                            PostParams%redo_add =.true.
                            exit
                        end if
                    end do
                end if
                if (PostParams%redo_skip>0.d0 .and. PostParams%redo_skip<1) then
                    at_beginning=at_beginning+1
                    if (at_beginning==1) then
                        CALL PXFFTELL (infile_handle,last_file_loc,ierror)
                        cycle
                    elseif (at_beginning==2) then
                        CALL PXFFTELL (infile_handle,file_loc,ierror)
                        inquire(unit=infile_handle, size=file_size)
                        PostParams%redo_skip = file_size/(file_loc-last_file_loc) * PostParams%redo_skip
                        if (Feedback > 0) print *,'skipping ',nint(PostParams%redo_skip), ' models'
                    end if
                end if
            end if

            if (error ==1) then
                if (num==0) call MpiStop('Error reading data file.')
                exit
            end if

            if (num<=PostParams%redo_skip .or. mod(num,PostParams%redo_thin) /= 0) cycle

            num_used=num_used+1

            if (PostParams%redo_like .or. PostParams%redo_add) then
                !Check for new prior before calculating anything
                if (CheckPriorCuts(Params)==logZero) then
                    if (Feedback >1) write(*,*) 'Model outside new prior bounds: skipped'
                    cycle
                end if
            end if

            if (PostParams%redo_theory ) then
                call GetTheoryForImportance(Params%P, newTheory, error, PostParams%redo_cls, PostParams%redo_pk)

                if (PostParams%redo_cls) then
                    Params%Theory%cl = newTheory%cl
                end if

                if (PostParams%redo_pk) then
                    Params%Theory%sigma_8 = newTheory%sigma_8
                    call InitPK(Params%Theory,newTheory%num_k, size(newTheory%redshifts))
                    Params%Theory%num_k = newTheory%num_k
                    Params%Theory%log_kh = newTheory%log_kh
                    Params%Theory%Matter_Power = newTheory%Matter_Power
                    Params%Theory%ddmatter_power = newTheory%ddmatter_power
                    Params%Theory%redshifts = newTheory%redshifts
                    if(use_nonlinear)then
                        Params%Theory%nlMatter_Power = newTheory%nlMatter_Power
                        Params%Theory%ddnlmatter_power = newTheory%ddnlmatter_power
                    end if
                end if

                Params%Theory%derived_parameters = newTheory%derived_parameters
                Params%Theory%numderived = newTheory%numderived
            else
                error = 0
            end if

            if (error ==0) then
                if (PostParams%redo_like .or. PostParams%redo_add) then
                    if (Use_LSS) then
                        if(Params%Theory%sigma_8==0) &
                        call MpiStop('ERROR: Matter power/sigma_8 have not been computed. Use redo_theory and redo_pk')

                        numz = size(Params%Theory%redshifts)
                        if((power_redshifts(num_power_redshifts)-Params%Theory%redshifts(numz))>1.d-3)then
                            write(*,*) 'ERROR: Thes elected datasets call for a higher redshift than has been calculated'
                            write(*,*) '       Use redo_theory and redo_pk'
                            call MpiStop()
                        end if
                        if(num_power_redshifts > numz)then
                            write(*,*) 'ERROR: The selected datasets call for more redshifts than are calculated'
                            write(*,*) '       Use redo_theory and redo_pk'
                            call MpiStop()
                        end if
                        index_error =0
                        call IndexExactRedshifts(Params%Theory%redshifts,index_error)
                        if(index_error>0)then
                            write(*,*) 'ERROR: One of the datasets needs an exact redshift that is not present '
                            write(*,*) '       Use redo_theory and redo_pk'
                            call MpiStop()
                        end if
                    end if
                    if (PostParams%redo_add) then
                        truelike = GetLogLikePost(Params, .not. has_likes)
                    else
                        truelike = GetLogLikePost(Params)
                    end if
                    if (truelike == logZero) then
                        weight = 0
                    else
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
                    !                    if (txt_theory) then
                    !                        call WriteParamsAndDat(Params, mult,like)
                    !                    else
                    call WriteParams(Params, mult,like)
                    !                   end if
                    if (outdata_handle>=0) call Params%WriteModel(outdata_handle, truelike,mult)
                else
                    if (Feedback >1 ) write (*,*) 'Zero weight: new like = ', truelike
                end if

                if (Feedback > 1) write (*,*) num, ' mult= ', real(mult), ' weight = ', real(weight)
                weight_max = max(weight,weight_max)
                weight_min = min(weight,weight_min)
                mult_max = max(mult_max,mult)

            end if

        end do

        call IO_Close(infile_handle)
        call IO_Close(outfile_handle)
        if (outdata_handle >=0) call IO_DataCloseWrite(outdata_handle)

        if (Feedback>0) then
            write(*,*) 'finished. Processed ',num_used,' models'
            write (*,*) 'max weight= ',weight_max, ' min weight = ',weight_min
            write (*,*) 'mean mult  = ', mult_sum/num_used
            write (*,*) 'mean importance weight (approx evidence ratio) = ',mult_ratio/num_used
            write (*,*) 'effective number of samples =',mult_sum/mult_max
            write (*,*) 'Best redo_likeoffset = ',max_truelike - max_like
        end if

        if ((mult_ratio < 1e-6 .or. mult_ratio > 1e8) .and. .not.PostParams%redo_change_like_only) then
            write (*,*) 'WARNING: use redo_likeoffset to rescale likelihoods'
        end if

    end if

    end subroutine postprocess

    end module posthoc