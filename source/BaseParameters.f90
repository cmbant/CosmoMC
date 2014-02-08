    module BaseParameters
    use settings
    use GeneralTypes
    use propose
    implicit none

    private

    integer, parameter :: tp_unused = 0, tp_slow=1, tp_semislow=2, tp_semifast=3, tp_fast=4

    integer, parameter :: slow_tp_max = tp_semislow

    Type ParamGaussPrior
        real(mcp), allocatable :: mean(:),std(:) !std=0 for no prior (default)
    end Type ParamGaussPrior

    Type :: TBaseParameters
        logical :: use_fast_slow = .false.
        integer :: num_fast, num_slow
        real(mcp), allocatable :: PMin(:), PMax(:), StartWidth(:), PWidth(:), center(:)
        logical(mcp), allocatable :: varying(:)
        Type(ParamGaussPrior) :: GaussPriors
        Type(int_arr_pointer), allocatable :: param_blocks(:)
        real(mcp), dimension(:,:), allocatable ::  covariance_estimate
        logical :: covariance_is_diagonal, covariance_has_new
    contains
    procedure :: ParamError => TBaseParameters_ParamError
    procedure :: SetStartPositions => TBaseParameters_SetStartPositions
    procedure :: ReadParams => TBaseParameters_ReadParams
    procedure :: ReadPriors => TBaseParameters_ReadPriors
    procedure :: ReadSetCovMatrix => TBaseParameters_ReadSetCovMatrix
    procedure :: InitializeUsedParams => TBaseParameters_InitializeUsedParams
    procedure :: OutputParamRanges => TBaseParameters_OutputParamRanges 
    procedure :: SetFastSlowParams => TBaseParameters_SetFastSlowParams
    procedure :: UsedParamNameOrNumber => TBaseParameters_UsedParamNameOrNumber
    procedure :: SetCovmat => TBaseParameters_SetCovmat
    end Type

    Type(TBaseParameters) :: BaseParams

    public TBaseParameters, BaseParams, slow_tp_max

    contains

    subroutine TBaseParameters_InitializeUsedParams(this, Ini, Params)
    class(TBaseParameters) :: this
    class(TIniFile) :: Ini
    class(TCalculationAtParamPoint) Params

    num_params = max(num_theory_params, NameMapping%num_MCMC)
    num_data_params = num_params - num_theory_params
    output_lines = 0

    call this%ReadParams(Ini)
    call this%ReadPriors(Ini)

    Params%P(1:num_params) = this%Center(1:num_params)

    end subroutine TBaseParameters_InitializeUsedParams


    subroutine TBaseParameters_ParamError(this, str,param)
    class(TBaseParameters) :: this
    character(LEN=*), intent(in) :: str
    integer, intent(in) :: param

    call DoAbort(trim(str)//': '//trim(NameMapping%NameOrNumber(param)))

    end subroutine TBaseParameters_ParamError


    subroutine TBaseParameters_SetStartPositions(this,Params)
    class(TBaseParameters) :: this
    class(TCalculationAtParamPoint) Params
    integer ix,i

    do ix=1,num_params_used
        i = params_used(ix)
        do
            if (this%StartWidth(i) < 0) then
                !This case we want half gaussian, width -BaseParameters%StartWidth(i)
                !e.g. for positive definite parameters
                Params%P(i) = this%center(i) - abs(Gaussian1())*this%StartWidth(i)
            else
                Params%P(i) = this%center(i) + Gaussian1()*this%StartWidth(i)
            end if
            !Repeat until get acceptable values in range
            if (Params%P(i)>=  this%PMin(i) .and. Params%P(i) <= this%PMax(i)) exit
        end do
    end do

    end subroutine TBaseParameters_SetStartPositions

    subroutine TBaseParameters_ReadParams(this,Ini)
    class(TBaseParameters) :: this
    class(TIniFile) :: Ini
    integer i, j, status
    character(LEN=:), allocatable :: InLine
    real(mcp) center

    allocate(this%PMin(num_params))
    allocate(this%PMax(num_params))
    allocate(this%StartWidth(num_params))
    allocate(this%PWidth(num_params))
    allocate(this%varying(num_params))

    do i=1,num_params
        InLine =  NameMapping%ReadIniForParam(Ini,'param',i)
        if (InLine=='') call this%ParamError('parameter ranges not found',i)
        if (TxtNumberColumns(InLine)==1) then
            !One number means just fix the parameter
            read(InLine, *, IOSTAT = status) center
            if (status/=0) call this%ParamError('fixed parameter value not valid',i)
            this%PMin(i)= center
            this%PMax(i)= center
            this%StartWidth(i)=0
            this%PWidth(i)=0
        else
            read(InLine, *, iostat=status) center, this%PMin(i), this%PMax(i), this%StartWidth(i), this%PWidth(i)
            if (status/=0) call DoAbort('Error reading param details: '//trim(InLIne))
        end if

        this%varying(i) = this%PWidth(i)/=0

        if (.not. this%varying(i)) then
            this%PMin(i) = center
            this%PMax(i) = center
        end if

        this%center(i) = center
        if (this%PMax(i) < this%PMin(i)) call this%ParamError('You have param Max < Min',i)
        if (this%center(i) < this%PMin(i)) call this%ParamError('You have param center < Min',i)
        if (this%center(i) > this%PMax(i)) call this%ParamError('You have param center > Max',i)
    end do

    num_params_used=count(this%varying)
    allocate(params_used(num_params_used))
    j=0
    do i=1,num_params
        if (this%varying(i)) then
            j=j+1
            params_used(j)=i
        end if
    end do

    end subroutine TBaseParameters_ReadParams

    subroutine TBaseParameters_ReadPriors(this,Ini)
    class(TBaseParameters) :: this
    class(TIniFile) :: Ini
    integer i, status
    character(LEN=:), allocatable :: InLine

    allocate(this%GaussPriors%std(num_params))
    allocate(this%GaussPriors%mean(num_params))
    this%GaussPriors%std=0 !no priors by default
    do i=1,num_params
        if (this%varying(i)) then
            InLine =  NameMapping%ReadIniForParam(Ini,'prior',i)
            if (InLine/='') read(InLine, *, iostat=status) this%GaussPriors%mean(i), this%GaussPriors%std(i)
            if (status/=0) call DoAbort('Error reading prior mean and stdd dev: '//trim(InLIne))
        end if
    end do

    end subroutine TBaseParameters_ReadPriors


    subroutine TBaseParameters_SetCovmat(this, prop_mat)
    class(TBaseParameters) :: this
    character(LEN=*), intent(in) :: prop_mat
    logical :: hasNew

    allocate(this%covariance_estimate(num_params_used, num_params_used))

    call this%ReadSetCovMatrix(prop_mat, this%covariance_estimate, hasNew)

    this%covariance_is_diagonal = prop_mat==''
    this%covariance_has_new = hasNew

    end subroutine TBaseParameters_SetCovmat

    subroutine TBaseParameters_ReadSetCovMatrix(this, prop_mat, matrix, NewParams)
    class(TBaseParameters) :: this
    character(LEN=*), intent(in) :: prop_mat
    real(mcp) :: matrix(num_params_used, num_params_used)
    logical, optional :: NewParams
    logical :: HasNewParams
    real(mcp) pmat(num_params,num_params)
    integer i

    if (prop_mat/='') then
        call IO_ReadProposeMatrix(pmat, prop_mat)
        HasNewParams = .false.
        !If generated with constrained parameters, assume diagonal in those parameters
        do i=1,num_params
            if (pmat(i,i) ==0 .and. this%varying(i)) then
                pmat(i,i) = this%PWidth(i)**2
                HasNewParams = .true.
            end if
            !Enforce new constraints (should really be fixing the inverse...)
            if (.not. this%varying(i)) then
                pmat(i,:) = 0
                pmat(:,i) = 0
            end if
        end do

        matrix = pmat(params_used, params_used)
    else
        matrix = 0
        HasNewParams = .true.
        do i=1,num_params_used
            matrix(i,i) = this%PWidth(params_used(i))**2
        end do
    end if
    if (present(NewParams)) NewParams = HasNewParams

    end subroutine TBaseParameters_ReadSetCovMatrix

    subroutine TBaseParameters_OutputParamRanges(this,Names, fname)
    class(TBaseParameters) :: this
    class(TParamNames) :: Names
    character(len=*), intent(in) :: fname
    integer unit,i

    unit = CreateNewTxtFile(fname)
    do i=1, Names%num_MCMC
        write(unit,'(1A22,2E17.7)') Names%NameOrNumber(i), this%PMin(i),this%PMax(i)
    end do
    close(unit)

    end subroutine TBaseParameters_OutputParamRanges


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


    subroutine TBaseParameters_SetFastSlowParams(this,Ini, use_fast_slow)
    class(TBaseParameters) :: this
    class(TIniFile) :: Ini
    logical, intent(in) :: use_fast_slow
    integer i, j, ix
    integer fast_params(num_params)
    integer fast_number, fast_param_index
    integer param_type(num_params)
    integer speed, num_speed
    logical :: block_semi_fast =.false., block_fast_likelihood_params=.false.
    integer :: breaks(num_params), num_breaks
    class(DataLikelihood), pointer :: DataLike
    logical first
    integer status

    this%use_fast_slow= use_fast_slow
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

    this%num_fast=0
    this%num_slow=0
    allocate(this%param_blocks(tp_semifast+num_breaks))
    do speed= tp_slow, tp_semifast
        allocate(this%param_blocks(speed)%P(count(param_type == speed)))
        num_speed=0
        do i=1,num_params_used
            if (param_type(params_used(i))==speed) then
                if (speed <=slow_tp_max) then
                    this%num_slow=this%num_slow+1
                else
                    this%num_fast = this%num_fast+1
                end if
                num_speed=num_speed+1
                this%param_blocks(speed)%P(num_speed) = i
            end if
        end do
    end do
    ix=1
    do speed = tp_fast, tp_fast+num_breaks-1
        j = breaks(speed-tp_fast+1)
        allocate(this%param_blocks(speed)%P(count(param_type(params_used(ix:j)) == tp_fast)))
        num_speed=0
        do i=ix,j
            if (param_type(params_used(i))==tp_fast) then
                num_speed=num_speed+1
                this%param_blocks(speed)%P(num_speed) = i
            end if
        end do
        ix=j+1
    end do

    if (Feedback > 0 .and. MpiRank==0) then
        if (use_fast_slow) then
            write(*,'(1I3," parameters (",1I2," slow (",1I2," semi-slow), ",1I2," fast (",1I2," semi-fast))")') &
            num_params_used,this%num_slow, size(this%param_blocks(tp_semislow)%P), this%num_fast,size(this%param_blocks(tp_semifast+1)%P)
        else
            write(*,'(1I3," parameters")') num_params_used
        end if
    end if

    end subroutine TBaseParameters_SetFastSlowParams


    function TBaseParameters_UsedParamNameOrNumber(this,i) result(name)
    class(TBaseParameters) :: this
    character(len=:), allocatable  :: name
    integer, intent(in) :: i

    name = NameMapping%NameOrNumber(params_used(i))

    end function TBaseParameters_UsedParamNameOrNumber

    end module BaseParameters

