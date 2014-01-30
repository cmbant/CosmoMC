    module BaseParameters
    use settings
    use GeneralTypes
    use propose
    implicit none

    private

    Type ParamGaussPrior
        real(mcp), allocatable :: mean(:),std(:) !std=0 for no prior (default)
    end Type ParamGaussPrior

    Type :: TBaseParameters
        real(mcp), allocatable :: PMin(:), PMax(:), StartWidth(:), PWidth(:), center(:)
        logical(mcp), allocatable :: varying(:)
        Type(ParamGaussPrior) :: GaussPriors
    contains
    procedure :: ParamError => TBaseParameters_ParamError
    procedure :: SetStartPositions => TBaseParameters_SetStartPositions
    procedure :: ReadParams => TBaseParameters_ReadParams
    procedure :: ReadPriors => TBaseParameters_ReadPriors
    procedure :: ReadSetCovMatrix => TBaseParameters_ReadSetCovMatrix
    procedure :: InitializeUsedParams => TBaseParameters_InitializeUsedParams
    procedure :: OutputParamRanges => TBaseParameters_OutputParamRanges 
    end Type

    Type(TBaseParameters) :: BaseParams

    public TBaseParameters, BaseParams
    
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
    character(Ini_max_string_len) :: InLine
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
    character(Ini_max_string_len) :: InLine

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


    end module BaseParameters

