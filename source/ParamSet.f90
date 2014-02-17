
    module ParamPointSet
    use settings
    use Samples
    use likelihood
    use GeneralTypes
    use BaseParameters
    use MatrixUtils
    implicit none

    Type, extends(TCalculationAtParamPoint) :: ParamSet
        real(mcp), allocatable :: likes(:)
        logical :: first = .true.
        integer :: numlikes
        integer, allocatable  :: like_indices(:)
        integer, allocatable ::  current_param_indices(:)
    contains
    procedure :: ReadModel => ParamSet_ReadModel
    procedure :: WriteModel => ParamSet_WriteModel
    end Type


    logical :: estimate_propose_matrix = .false.

    contains

    subroutine DoStop(S, abort)
    character(LEN=*), intent(in), optional :: S
    integer ierror
    logical, intent(in), optional :: abort
    logical wantbort

    if (outfile_handle/=0) close(outfile_handle)

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


    !subroutine IO_WriteCovMat(fname, matrix)
    !integer unit
    !character(LEN=*), intent(in) :: fname
    !character(LEN=4096) outline
    !real(mcp), intent(in) :: matrix(:,:)
    !
    !if (NameMapping%nnames/=0) then
    !    outline=''
    !    do unit=1, num_params_used
    !        outline = trim(outline)//' '//trim(BaseParams%UsedParamNameOrNumber(unit))
    !    end do
    !    call IO_WriteProposeMatrix(matrix ,fname, outline)
    !else
    !    call Matrix_write(fname,matrix,forcetable=.true.)
    !end if
    !end subroutine IO_WriteCovMat

    subroutine ParamSet_WriteModel(this, unit, like, mult)
    Class(ParamSet) :: this
    integer unit
    real(mcp), intent(in) :: mult, like
    integer j , alen, unused
    logical, save :: first = .true.
    class(DataLikelihood), pointer :: DataLike

    if (first .and. new_chains) then
        first = .false.
        if (mcp==kind(1.0)) then
            j=3
        else
            j=4
        end if
        write(unit) j, num_params_used
        if (.not. any (BaseParams%NameMapping%Name=='')) then
            write(unit) .true.
            do j=1,num_params_used
                alen=len_trim(BaseParams%NameMapping%name(params_used(j)))
                write(unit) alen
                write(unit) BaseParams%NameMapping%name(params_used(j))(1:alen)
            end do
        else
            write(unit) .false.
        end if
        write(unit) DataLikelihoods%Count
        do j=1, DataLikelihoods%Count
            DataLike => DataLikelihoods%Item(j)
            alen = len_trim(dataLIke%name)
            write(unit) alen
            write(unit) dataLIke%name(1:alen)
        end do
        unused=0
        write(unit) unused
    end if

    write(unit) mult, like
    write(unit) this%Likelihoods(1:DataLikelihoods%Count)
    write(unit) this%P(params_used)

    call this%Theory%WriteTheory(unit)

    if (flush_write) call FlushFile(unit)

    end subroutine ParamSet_WriteModel

    subroutine ParamSet_ReadModel(this,  unit, has_likes, mult, like, error)
    Class (ParamSet) :: this
    integer, intent(in) :: unit
    integer, intent(out) :: error
    real(mcp), intent(out) :: mult, like
    logical, intent(out) :: has_likes(:)
    integer j, k, np, alen, unused, status
    character(LEN=:), allocatable :: name
    logical :: has_names
    Type(DataLikelihood), pointer :: DataLike
    integer tmp(1)
    
    error = 0
    if (this%first) then
        this%first = .false.
        read(unit,iostat=status) j, np
        if (status/=0) then
            error=1
            return
        end if
        if (j/=3 .and. mcp==kind(1.0) .or. j/=4 .and. mcp/=kind(1.0)) &
        call MpiStop('ReadModel: wrong file format (old cosmomc version?)')
        if (np/=num_params_used) call MpiStop('ReadModel: number of parameters changed')
        read(unit) has_names
        allocate(this%current_param_indices(num_params_used))
        this%current_param_indices=-1
        if (has_names) then
            do j=1,num_params_used
                read(unit) alen
                allocate(character(alen)::name)
                read(unit) name
                this%current_param_indices(j) = BaseParams%NameMapping%index(name)
                deallocate(name)
            end do
            if (any(this%current_param_indices==-1)) call MpiStop('ReadModel: parameters in .data files could not be matched')
        else
            this%current_param_indices = params_used
        end if
        read(unit) this%numlikes
        allocate(this%likes(this%numlikes))
        allocate(this%like_indices(this%numlikes))
        this%like_indices=0
        do j=1, this%numlikes
            read(unit) alen
            allocate(character(alen)::name)
            read(unit) name
            do k=1, DataLikelihoods%Count
                DataLike => DataLikelihoods%Item(k)
                if (DataLike%name==name) then
                    this%like_indices(j)=k
                    exit
                end if
            end do
            deallocate(name)
        end do
        do j=1, DataLikelihoods%Count
            has_likes(j) = any(this%like_indices==j)
        end do
        read(unit) unused
        if (unused/=0) call MpiStop('ReadModel: Don''t know what extra info is')
        if (unused>0) read(unit) tmp(1:unused)
    end if

    this%Likelihoods=0
    read(unit,iostat=status) mult, like
    if (status/=0) then
        error=1
        return
    end if
    read(unit) this%likes(1:this%numlikes)
    do j=1,this%numlikes
        if (this%like_indices(j)/=0) this%Likelihoods(this%like_indices(j)) = this%likes(j)
    end do
    this%P= BaseParams%center
    read(unit) this%P(this%current_param_indices)
    call this%Theory%ReadTheory(unit)

    end subroutine ParamSet_ReadModel


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