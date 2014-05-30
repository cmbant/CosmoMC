
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
    final :: FreeParams
    end Type

    contains

    subroutine FreeParams(this)
    Type(ParamSet) this
    !This should not be needed? but does solve memory leak in ifort compilation..
    if (allocated(this%Theory)) deallocate(this%Theory)
    end subroutine


    subroutine ParamSet_WriteModel(this, F, like, mult)
    Class(ParamSet) :: this
    class(TFileStream) :: F
    real(mcp), intent(in) :: mult, like
    integer j , unused
    logical isfirst
    class(TDataLikelihood), pointer :: DataLike

    isfirst = F%Position() == 1
    if (isfirst) then
        if (mcp==kind(1.0)) then
            j=3
        else
            j=4
        end if
        write(F%unit) j, num_params_used
        if (.not. any (BaseParams%NameMapping%Name=='')) then
            write(F%unit) .true.
            do j=1,num_params_used
                call F%WriteTrim(BaseParams%NameMapping%name(params_used(j)))
            end do
        else
            write(F%unit) .false.
        end if
        write(F%unit) DataLikelihoods%Count
        do j=1, DataLikelihoods%Count
            DataLike => DataLikelihoods%Item(j)
            call F%WriteTrim(dataLIke%name)
        end do
        unused=0
        write(F%unit) unused
    end if

    write(F%unit) mult, like
    write(F%unit) this%Likelihoods(1:DataLikelihoods%Count)
    write(F%unit) this%P(params_used)

    call this%Theory%WriteTheory(F, isfirst)

    if (flush_write) call F%Flush()

    end subroutine ParamSet_WriteModel

    subroutine ParamSet_ReadModel(this,  F, has_likes, mult, like, error)
    Class (ParamSet) :: this
    class(TFileStream) :: F
    integer, intent(out) :: error
    real(mcp), intent(out) :: mult, like
    logical, intent(out) :: has_likes(:)
    integer j, k, np, unused, status
    character(LEN=:), allocatable :: name
    logical :: has_names
    class(TDataLikelihood), pointer :: DataLike
    integer tmp(1)
    logical isfirst

    error = 0
    isfirst = this%first
    if (isfirst) then
        this%first = .false.
        read(F%unit,iostat=status) j, np
        if (status/=0) then
            error=1
            return
        end if
        if (j/=3 .and. mcp==kind(1.0) .or. j/=4 .and. mcp/=kind(1.0)) &
        call MpiStop('ReadModel: wrong file format (old cosmomc version?)')
        if (np/=num_params_used) call MpiStop('ReadModel: number of parameters changed')
        read(F%unit) has_names
        allocate(this%current_param_indices(num_params_used))
        this%current_param_indices=-1
        if (has_names) then
            do j=1,num_params_used
                name = F%ReadStringItem()
                this%current_param_indices(j) = BaseParams%NameMapping%index(name)
            end do
            if (any(this%current_param_indices==-1)) call MpiStop('ReadModel: parameters in .data files could not be matched')
        else
            this%current_param_indices = params_used
        end if
        read(F%unit) this%numlikes
        allocate(this%likes(this%numlikes))
        allocate(this%like_indices(this%numlikes))
        this%like_indices=0
        do j=1, this%numlikes
            name = F%ReadStringItem()
            do k=1, DataLikelihoods%Count
                DataLike => DataLikelihoods%Item(k)
                if (DataLike%name==name) then
                    this%like_indices(j)=k
                    exit
                end if
            end do
        end do
        do j=1, DataLikelihoods%Count
            has_likes(j) = any(this%like_indices==j)
        end do
        read(F%unit) unused
        if (unused/=0) call MpiStop('ReadModel: Don''t know what extra info is')
        if (unused>0) read(F%unit) tmp(1:unused)
    end if

    this%Likelihoods=0
    read(F%unit,iostat=status) mult, like
    if (status/=0) then
        error=1
        return
    end if
    read(F%unit) this%likes(1:this%numlikes)
    do j=1,this%numlikes
        if (this%like_indices(j)/=0) this%Likelihoods(this%like_indices(j)) = this%likes(j)
    end do
    this%P(1:num_params)= BaseParams%center
    read(F%unit) this%P(this%current_param_indices)
    call this%Theory%ReadTheory(F, isfirst)

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
