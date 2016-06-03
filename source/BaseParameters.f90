    module BaseParameters
    use settings
    use GeneralTypes
    use likelihood
    use ObjectLists
    use IO
    use RandUtils
    implicit none
    private

    integer, parameter :: tp_unused = 0, tp_slow=1, tp_semislow=2, tp_semifast=3, tp_fast=4

    integer, parameter :: slow_tp_max = tp_semislow

    Type ParamGaussPrior
        real(mcp), allocatable :: mean(:),std(:) !std=0 for no prior (default)
    end Type ParamGaussPrior

    Type TLinearCombination
        real(mcp) :: mean =0._mcp
        real(mcp) :: std = 0._mcp
        character(LEN=:), allocatable :: Name
        real(mcp), allocatable :: Combination(:)
    end Type TLinearCombination

    Type :: TBaseParameters
        logical :: use_fast_slow = .false.
        integer :: num_fast, num_slow
        integer :: num_semi_fast, num_semi_slow
        logical :: include_fixed_parameter_priors = .false.
        real(mcp), allocatable :: PMin(:), PMax(:), StartWidth(:), PWidth(:), center(:)
        logical(mcp), allocatable :: varying(:)
        logical :: block_semi_fast = .true.
        logical :: block_fast_likelihood_params = .true.
        Type(ParamGaussPrior) :: GaussPriors
        Type(TLinearCombination), allocatable :: LinearCombinations(:)
        Type(int_arr), allocatable :: param_blocks(:)
        real(mcp), dimension(:,:), allocatable ::  covariance_estimate
        logical :: covariance_is_diagonal, covariance_has_new
        Type(TParamNames) :: NameMapping
    contains
    procedure :: ParamError => TBaseParameters_ParamError
    procedure :: SetStartPositions => TBaseParameters_SetStartPositions
    procedure :: ReadParams => TBaseParameters_ReadParams
    procedure :: ReadPriors => TBaseParameters_ReadPriors
    procedure :: ReadSetCovMatrix => TBaseParameters_ReadSetCovMatrix
    procedure :: InitializeUsedParams => TBaseParameters_InitializeUsedParams
    procedure :: OutputParamRanges => TBaseParameters_OutputParamRanges
    procedure :: OutputParamNames => TBaseParameters_OutputParamNames
    procedure :: SetFastSlowParams => TBaseParameters_SetFastSlowParams
    procedure :: UsedParamNameOrNumber => TBaseParameters_UsedParamNameOrNumber
    procedure :: SetCovmat => TBaseParameters_SetCovmat
    end Type

    Type(TBaseParameters), save :: BaseParams

    public TBaseParameters, BaseParams, slow_tp_max

    contains

    subroutine TBaseParameters_InitializeUsedParams(this, Ini)
    class(TBaseParameters) :: this
    class(TSettingIni) :: Ini

    num_params = max(num_theory_params, this%NameMapping%num_MCMC)
    num_data_params = num_params - num_theory_params
    output_lines = 0

    call this%ReadParams(Ini)
    call this%ReadPriors(Ini)

    end subroutine TBaseParameters_InitializeUsedParams


    subroutine TBaseParameters_ParamError(this, str,param)
    class(TBaseParameters) :: this
    character(LEN=*), intent(in) :: str
    integer, intent(in) :: param

    call DoAbort(trim(str)//': '//trim(this%NameMapping%NameOrNumber(param)))

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
    class(TSettingIni) :: Ini
    integer i, j, status
    character(LEN=:), allocatable :: InLine
    real(mcp) center

    allocate(this%center(num_params))
    allocate(this%PMin(num_params))
    allocate(this%PMax(num_params))
    allocate(this%StartWidth(num_params))
    allocate(this%PWidth(num_params))
    allocate(this%varying(num_params))

    do i=1,num_params
        InLine =  this%NameMapping%ReadIniForParam(Ini,'param',i)
        if (InLine=='') call this%ParamError('parameter ranges not found',i)
        if (File%TxtNumberColumns(InLine)==1) then
            !One number means just fix the parameter
            read(InLine, *, IOSTAT = status) center
            if (status/=0) call this%ParamError('fixed parameter value not valid',i)
            this%PMin(i)= center
            this%PMax(i)= center
            this%StartWidth(i)=0
            this%PWidth(i)=0
        else
            read(InLine, *, iostat=status) center, this%PMin(i), this%PMax(i), this%StartWidth(i), this%PWidth(i)
            if (status/=0) call this%ParamError('Error reading param details: '//trim(InLIne),i)
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
    class(TSettingIni) :: Ini
    integer i, status
    character(LEN=:), allocatable :: InLine
    Type(TSettingIni) :: Combs
    integer params(num_params), num_lin

    call Ini%Read('include_fixed_parameter_priors',this%include_fixed_parameter_priors)
    allocate(this%GaussPriors%std(num_params))
    allocate(this%GaussPriors%mean(num_params))
    this%GaussPriors%std=0 !no priors by default
    do i=1,num_params
        if (this%varying(i) .or. this%include_fixed_parameter_priors) then
            InLine =  this%NameMapping%ReadIniForParam(Ini,'prior',i)
            if (InLine/='') then
                read(InLine, *, iostat=status) this%GaussPriors%mean(i), this%GaussPriors%std(i)
                if (status/=0) call this%ParamError('Error reading prior mean and stdd dev: '//trim(InLIne),i)
            end if
        end if
    end do

    call Ini%TagValuesForName('linear_combination', Combs)
    allocate(this%LinearCombinations(Combs%Count))
    do i= 1, Combs%Count
        associate( comb =>this%LinearCombinations(i) )
            Comb%name = Combs%Name(i)
            num_lin = -1
            call this%NameMapping%ReadIndices(Combs%Value(i), params, num_lin)
            allocate(Comb%Combination(num_params),source = 0._mcp)
            InLine = Ini%Read_String(Ini%NamedKey('linear_combination_weights',Comb%name), NotFoundFail=.true.)
            read(InLine,*, iostat=status) Comb%Combination(params(:num_lin))
            if (status/=0) call MpiStop('Error reading linear_combination_weights: '//trim(InLIne))
            InLine = Ini%Read_String(Ini%NamedKey('prior',Comb%name))
            if (InLine/='') then
                read(InLine,*, iostat=status) Comb%mean, Comb%std
                if (status/=0) call MpiStop('Error reading linear_combination_prior: '//trim(InLIne))
            end if
        end associate
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
        call IO_ReadProposeMatrix(this%NameMapping,pmat, prop_mat)
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


    subroutine TBaseParameters_OutputParamNames(this,fname, indices, add_derived)
    class(TBaseParameters) :: this
    character(len=*), intent(in) :: fname
    integer, intent(in), optional :: indices(:)
    logical, intent(in), optional :: add_derived

    call this%NameMapping%WriteFile(trim(fname)//'.paramnames', indices, add_derived)

    end subroutine TBaseParameters_OutputParamNames

    subroutine TBaseParameters_OutputParamRanges(this,fname)
    class(TBaseParameters) :: this
    character(len=*), intent(in) :: fname
    integer i
    Type(TTextFile) :: F

    call F%CreateFile(fname//'.ranges')
    do i=1, this%NameMapping%num_MCMC
        call F%WriteLeftAligned('(1A22)', this%NameMapping%NameOrNumber(i))
        call F%Write([this%PMin(i),this%PMax(i)])
    end do
    call F%Close()

    end subroutine TBaseParameters_OutputParamRanges

    subroutine orderIndices(arr,n)
    use ArrayUtils
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
    class(TSettingIni) :: Ini
    logical, intent(in) :: use_fast_slow
    integer i, j, ix
    integer fast_params(num_params)
    integer fast_number, fast_param_index
    integer param_type(num_params)
    integer speed, num_speed
    integer :: breaks(num_params), num_breaks
    class(TDataLikelihood), pointer :: DataLike
    logical first
    integer status

    this%use_fast_slow= use_fast_slow
    if (use_fast_slow) then
        if (Ini%HasKey('fast_parameters')) then
            !List of parmeter names to treat as fast
            fast_number = num_params
            call this%NameMapping%ReadIndices(Ini%Read_String('fast_parameters'), fast_params, fast_number)
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
        call Ini%Read('block_semi_fast',this%block_semi_fast)
        call Ini%Read('block_fast_likelihood_params',this%block_fast_likelihood_params)
    else
        fast_number = 0
        this%block_semi_fast = .false.
        this%block_fast_likelihood_params = .false.
    end if

    param_type = tp_unused
    do i=1,num_params
        if (BaseParams%varying(i)) then !to get sizes for allocation arrays
            if (use_fast_slow .and. any(i==fast_params(1:fast_number))) then
                if (i >= index_data .or. .not. this%block_semi_fast) then
                    param_type(i) = tp_fast
                else
                    param_type(i) = tp_semifast
                end if
            else
                if (use_fast_slow .and. index_semislow >=0 .and. i >= index_semislow .and. this%block_semi_fast) then
                    param_type(i) = tp_semislow
                else
                    param_type(i) = tp_slow
                end if
            end if
        end if
    end do

    num_breaks=0
    if (this%block_fast_likelihood_params) then
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
                this%num_fast = this%num_fast+1
                num_speed=num_speed+1
                this%param_blocks(speed)%P(num_speed) = i
            end if
        end do
        ix=j+1
    end do

    this%num_semi_slow = size(this%param_blocks(tp_semislow)%P)
    this%num_semi_fast = size(this%param_blocks(tp_semifast)%P)
    if (Feedback > 0 .and. MpiRank==0) then
        if (use_fast_slow) then
            write(*,'(1I3," parameters (",1I2," slow (",1I2," semi-slow), ",1I2," fast (",1I2," semi-fast))")') &
                & num_params_used,this%num_slow, this%num_semi_slow, this%num_fast, this%num_semi_fast
        else
            write(*,'(1I3," parameters")') num_params_used
        end if
    end if

    end subroutine TBaseParameters_SetFastSlowParams


    function TBaseParameters_UsedParamNameOrNumber(this,i) result(name)
    class(TBaseParameters) :: this
    character(len=:), allocatable  :: name
    integer, intent(in) :: i

    name = this%NameMapping%NameOrNumber(params_used(i))

    end function TBaseParameters_UsedParamNameOrNumber

    end module BaseParameters
