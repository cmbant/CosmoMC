    module Sampler
    use propose
    use BaseParameters
    use likelihood
    use GeneralTypes

    Type, extends(TConfigClass) :: TSampler
        logical :: StartCovMat = .false.
        real(mcp), dimension(:,:), allocatable ::  initial_propose_matrix
        class(BlockedProposer), allocatable :: Proposer
    contains
    procedure :: ReadCovmat => TSampler_ReadCovmat
    procedure :: SetProposalParams => TSampler_SetProposalParams
    end Type

    contains


    subroutine orderIndices(arr,n)
    integer, intent(in) :: n
    integer arr(:), tmp(n),j, s, i

    tmp=arr(1:n)
    s=n
    do i=1,n
        arr(i) = minval(tmp(1:s))
        j=indexOf(arr(i),tmp,s)
        tmp(j:s-1)= tmp(j+1:s)
        s=s-1
    end do
    end subroutine orderIndices

    subroutine TSampler_ReadCovmat(this,Ini, use_propose_matrix)
    class(TSampler) :: this
    class(TIniFile) :: Ini
    logical :: use_propose_matrix, hasNew
    character(LEN=Ini_max_string_len) prop_mat

    if (use_propose_matrix) then
        prop_mat = trim(Ini%Read_String_Default('propose_matrix',''))
        if (prop_mat /= '' .and. prop_mat(1:1) /= '/') prop_mat = concat(LocalDir,prop_mat)
    else
        prop_mat=''
    end if

    allocate(this%initial_propose_matrix(num_params_used, num_params_used))
    this%StartCovMat  = prop_mat/=''

    call BaseParams%ReadSetCovMatrix(prop_mat, this%initial_propose_matrix, hasNew) 
    if (hasNew) MPI_Max_R_ProposeUpdate = MPI_Max_R_ProposeUpdateNew

    call this%Proposer%SetCovariance(this%initial_propose_matrix)

    end subroutine TSampler_ReadCovmat


    subroutine TSampler_SetProposalParams(this,Ini)
    class(TSampler) :: this
    class(TIniFile) :: Ini
    integer i, j, ix
    integer fast_params(num_params)
    integer fast_number, fast_param_index
    integer param_type(num_params)
    integer speed, num_speed
    integer, parameter :: tp_unused = 0, tp_slow=1, tp_semislow=2, tp_semifast=3, tp_fast=4
    Type(int_arr_pointer), allocatable :: param_blocks(:)
    logical :: block_semi_fast =.false., block_fast_likelihood_params=.false.
    integer :: breaks(num_params), num_breaks
    class(DataLikelihood), pointer :: DataLike
    integer :: oversample_fast= 1
    logical first
    integer status


    if (use_fast_slow) then
        if (Ini%HasKey('fast_parameters')) then
            !List of parmeter names to treat as fast
            fast_number = num_params
            call NameMapping%ReadIndices(Ini%Read_String('fast_parameters'), fast_params, fast_number)
        else
            fast_param_index = max(index_data,DataLikelihoods%first_fast_param)
            !all parameters at and above fast_param_index
            fast_param_index= Ini%Read_Int('fast_param_index',fast_param_index)
            fast_number = 0
            do i=fast_param_index, num_params
                fast_number = fast_number+1
                fast_params(fast_number) = i
            end do
        end if
        block_semi_fast = Ini%Read_Logical('block_semi_fast',.true.)
        block_fast_likelihood_params = Ini%Read_logical('block_fast_likelihood_params',.true.)
        oversample_fast = Ini%Read_Int('oversample_fast',1)
        if (oversample_fast<1) call DoAbort('oversample_fast must be >= 1')
        if (sampling_method /=sampling_fast_dragging) then
            MPI_Thin_fac = MPI_Thin_fac*this%Proposer%Oversample_fast
        end if
    else
        fast_number = 0
    end if

    param_type = tp_unused
    do i=1,num_params
        if (BaseParams%varying(i)) then !to get sizes for allocation arrays
            if (use_fast_slow .and. any(i==fast_params(1:fast_number))) then
                if (i >= index_data .or. .not. block_semi_fast) then
                    param_type(i) = tp_fast
                else
                    param_type(i) = tp_semifast
                end if
            else
                if (use_fast_slow .and. index_semislow >=0 .and. i >= index_semislow .and. block_semi_fast) then
                    param_type(i) = tp_semislow
                else
                    param_type(i) = tp_slow
                end if
            end if
        end if
    end do

    num_breaks=0
    if (block_fast_likelihood_params) then
        !put parameters for different likelihoods in separate blocks,
        !so not randomly mix them and hence don't all need to be recomputed
        first=.true.
        do i=1,DataLikelihoods%Count
            DataLike=>DataLikelihoods%Item(i)
            do j=1, num_params_used-1
                if (param_type(params_used(j))==tp_fast .and. params_used(j) >= DataLike%new_param_block_start &
                .and. params_used(j) < DataLike%new_param_block_start + DataLike%new_params) then
                    if (first) then
                        first = .false.
                    else
                        num_breaks = num_breaks+1
                        breaks(num_breaks)=j
                    end if
                    exit
                end if
            end do
        end do
    end if
    num_breaks = num_breaks + 1
    breaks(num_breaks) = num_params_used
    if (Feedback>0 .and. MpiRank==0) then
        write(*,*) 'Fast divided into ',num_breaks,' blocks'
        if (num_breaks>1) write(*,*) 'Block breaks at: ',breaks(1:num_breaks-1)
    end if

    call orderIndices(breaks, num_breaks)

    allocate(param_blocks(tp_semifast+num_breaks))
    do speed= tp_slow, tp_semifast
        allocate(param_blocks(speed)%P(count(param_type == speed)))
        num_speed=0
        do i=1,num_params_used
            if (param_type(params_used(i))==speed) then
                num_speed=num_speed+1
                param_blocks(speed)%P(num_speed) = i
            end if
        end do
    end do
    ix=1
    do speed = tp_fast, tp_fast+num_breaks-1
        j = breaks(speed-tp_fast+1)
        allocate(param_blocks(speed)%P(count(param_type(params_used(ix:j)) == tp_fast)))
        num_speed=0
        do i=ix,j
            if (param_type(params_used(i))==tp_fast) then
                num_speed=num_speed+1
                param_blocks(speed)%P(num_speed) = i
            end if
        end do
        ix=j+1
    end do

    call this%Proposer%Init(param_blocks, slow_block_max= 2, oversample_fast=oversample_fast)
    num_slow = this%Proposer%Slow%n
    num_fast = this%Proposer%Fast%n

    if (Feedback > 0 .and. MpiRank==0) then
        write(*,'(" Varying ",1I2," parameters (",1I2," slow (",1I2," semi-slow), ",1I2," fast (",1I2," semi-fast))")') &
        num_params_used,num_slow, size(param_blocks(tp_semislow)%P), num_fast,size(param_blocks(tp_semifast)%P)
    end if

    end subroutine TSampler_SetProposalParams


    end module Sampler