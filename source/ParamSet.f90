
    module ParamPointSet
    use Random
    use settings
    use Propose
    use Samples
    use likelihood
    use BaseParameters
    use MatrixUtils
    implicit none

    Type, extends(TCalculationAtParamPoint) :: ParamSet
    contains
    procedure :: ReadModel
    procedure :: WriteModel
    end Type


    logical :: estimate_propose_matrix = .false.

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


    subroutine WriteIndepSample(P, like, mult)
    Type(ParamSet) P
    real(mcp) like, mult
    if (indepfile_handle ==0) return
    call P%WriteModel(indepfile_handle, like, mult)
    end subroutine WriteIndepSample



    subroutine IO_WriteCovMat(fname, matrix)
    integer i
    character(LEN=*), intent(in) :: fname
    character(LEN=4096) outline
    real(mcp), intent(in) :: matrix(:,:)

    if (NameMapping%nnames/=0) then
        outline=''
        do i=1, num_params_used
            outline = trim(outline)//' '//trim(BaseParams%UsedParamNameOrNumber(i))
        end do
        call IO_WriteProposeMatrix(matrix ,fname, outline)
    else
        call Matrix_write(fname,matrix,forcetable=.true.)
    end if
    end subroutine IO_WriteCovMat

    subroutine WriteModel(Params, i, like, mult)
    Class(ParamSet) :: Params
    integer i
    real(mcp), intent(in) :: mult, like
    integer j , len, unused
    logical, save :: first = .true.
    class(DataLikelihood), pointer :: DataLike

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
    integer j, k, np, len, unused, status
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
        read(i,iostat=status) j, np
        if (status/=0) then
            error=1
            return
        end if
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
                current_param_indices(j) = NameMapping%index(pname)
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
    read(i,iostat=status) mult, like
    if (status/=0) then
        error=1
        return
    end if
    read(i) likes(1:numlikes)
    do j=1,numlikes
        if (like_indices(j)/=0) Params%Likelihoods(like_indices(j)) = likes(j)
    end do
    Params%P= BaseParams%center
    read(i) Params%P(current_param_indices)
    call Params%Theory%ReadTheory(i)

    end subroutine ReadModel


    !subroutine WriteParamsAndDat(P, mult, like)
    !Type(ParamSet) P
    !real(mcp), intent(in) :: mult, like
    !real(mcp),allocatable :: output_array(:)
    !Type(mc_real_pointer) :: derived
    !integer numderived
    !
    !if (outfile_handle ==0) return
    !
    !numderived =Parameterization%CalcDerivedParams(P%P,P%Theory, derived)
    !
    !allocate(output_array(num_params_used + numderived + P%Theory%num_k ))
    !output_array(1:num_params_used) =  P%P(params_used)
    !output_array(num_params_used+1:num_params_used+numderived) =  derived%P
    !deallocate(derived%P)
    !
    !output_array(num_params_used+numderived+1:num_params_used+numderived+P%Theory%num_k) = &
    !P%Theory%matter_power(:,1)
    !
    !call IO_OutputChainRow(outfile_handle, mult, like, output_array)
    !deallocate(output_array)
    !
    !end  subroutine WriteParamsAndDat


    end module ParamPointSet