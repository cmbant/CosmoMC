    module GeneralTypes
    use settings

    Type TTheoryParams
        real(mcp) :: BaseParams(max_theory_params) = 0._mcp
    end Type TTheoryParams

    Type TTheoryIntermediateCache
        !Cache of intermediate steps in the theory calculation that may be (partly) reused for other points
        real(mcp) lastParamArray(max_num_params)
        logical :: validInfo = .false.
    contains
    procedure :: Clear => TTheoryIntermediateCache_Clear
    procedure :: Initialize => TTheoryIntermediateCache_Initialize
    end Type TTheoryIntermediateCache

    Type TCalculationAtParamPoint
        !Parameter values, calculated theory and intermediates and likelihoods at a particular point in parameter space
        real(mcp) :: P(max_num_params)
        real(mcp) :: likelihoods(max_likelihood_functions)
        class(TTheoryPredictions), pointer :: Theory => null()
        class(TTheoryIntermediateCache), pointer :: Info => null()
    contains
    procedure :: Clear => TCalculationAtParamPoint_Clear
    end Type TCalculationAtParamPoint

    !Each descendant of the base TConfigClass below has a Config element which point to a TGeneralConfig instance taht determines 
    !which related class implementations are actually used
    Type TGeneralConfig
        class(TParameterization), pointer :: Parameterization => null()
        class(TTheoryCalculator), pointer :: Calculator => null()
    contains
    procedure :: SetTheoryParameterization
    procedure :: SetParameterizationName
    end Type

    Type :: TConfigClass
        class(TGeneralConfig), pointer :: Config
    contains
    procedure :: InitConfig  => TConfigClass_InitConfig
    procedure :: InitWithParams  => TConfigClass_InitWithParams
    generic:: Init => InitWithParams,InitConfig
    end Type 

    Type, extends(TConfigClass) :: TTheoryPredictions
        !Actual computed theory predictions used by likelihoods
        !Config%Calculator can in some cases be used to provide theory functions
    contains
    procedure :: Clear => TTheoryPredictions_Clear
    procedure :: AssignNew => TTheoryPredictions_AssignNew
    end Type TTheoryPredictions

    Type, extends(TConfigClass) :: TTheoryCalculator
        character(LEN=128) :: calcName = 'TheoryCalculator'
    contains
    procedure :: AcceptReject => TTheoryCalculator_AcceptReject
    procedure :: Error => TTheoryCalculator_Error
    procedure :: ErrorNotImplemented => TTheoryCalculator_ErrorNotImplemented
    procedure :: VersionTraceOutput => TTheoryCalculator_VersionTraceOutput
    contains

    end Type TTheoryCalculator

    type, extends(TConfigClass)  :: TParameterization
    contains
    procedure :: ParamArrayToTheoryParams
    procedure :: CalcDerivedParams=> TParameterization_CalcDerivedParams
    procedure :: Initialize => TParameterization_Initialize
    procedure :: NonBaseParameterPriors => TParameterization_NonBaseParameterPriors
    procedure :: NewTheoryParams =>TParameterization_NewTheoryParams
    end type TParameterization

    Type, extends(TParameterization) :: GenericParameterization
    contains
    end type GenericParameterization

    contains

    subroutine ParamArrayToTheoryParams(this, Params, CMB)
    class(TParameterization) :: this
    real(mcp) Params(:)
    Class(TTheoryParams), target :: CMB

    CMB%BaseParams(1:num_theory_params) = Params(1:num_theory_params)

    end subroutine ParamArrayToTheoryParams

    subroutine TParameterization_Initialize(this, Ini, Names, DefaultName, Config)
    class(TParameterization) :: this
    class(TIniFile) :: Ini
    class(TParamNames) :: Names
    character(LEN=*), intent(in) :: DefaultName
    class(TGeneralConfig), target :: Config
    character(LEN=Ini_max_string_len) :: ParamNamesFile = ''

    call this%TConfigClass%Init(Config)

    ParamNamesFile = ReadIniFileName(Ini,'ParamNamesFile', NotFoundFail=.false.)
    if (ParamNamesFile =='' .and. DefaultName/='') ParamNamesFile= trim(LocalDir)//trim(DefaultName)

    if (ParamNamesFile /='') then
        call ParamNames_init(Names, ParamNamesFile)
        num_theory_params= Names%num_MCMC
    else
        Names%nnames=0
        num_theory_params= Ini_Read_Int_File(Ini, 'num_theory_params')
    end if
    if (num_theory_params> max_theory_params) call MpiStop('see settings.f90: num_theory_params> max_theory_params')
    index_data =  num_theory_params+1

    end subroutine TParameterization_Initialize

    subroutine TParameterization_NewTheoryParams(this,TheoryParams)
    class(TParameterization) :: this
    class(TTheoryParams), allocatable :: TheoryParams

    allocate(TTheoryParams::TheoryParams)

    end subroutine TParameterization_NewTheoryParams

    function TParameterization_NonBaseParameterPriors(this,CMB)
    class(TParameterization) :: this
    class(TTheoryParams) :: CMB
    real(mcp):: TParameterization_NonBaseParameterPriors

    TParameterization_NonBaseParameterPriors = 0

    end function TParameterization_NonBaseParameterPriors

    function TParameterization_CalcDerivedParams(this, P, Theory, derived) result (num_derived)
    class(TParameterization) :: this
    Type(mc_real_pointer) :: derived
    class(TTheoryPredictions) :: Theory
    real(mcp) :: P(:)
    integer num_derived

    num_derived = 0
    allocate(Derived%P(num_derived))

    end function TParameterization_CalcDerivedParams

    !!!TTheoryPredictions

    subroutine TTheoryPredictions_Clear(this)
    class(TTheoryPredictions) :: this

    end subroutine TTheoryPredictions_Clear

    subroutine TTheoryPredictions_AssignNew(this, NewTheory)
    class(TTheoryPredictions) :: this
    class(TTheoryPredictions), pointer :: NewTheory

    allocate(NewTheory, source = this)

    end subroutine TTheoryPredictions_AssignNew

    !!! TCalculationAtParamPoint

    subroutine TCalculationAtParamPoint_Clear(this, Keep)
    class(TCalculationAtParamPoint) :: this, Keep

    if (associated(this%Info) .and. .not. associated(this%Info, Keep%Info)) then
        call this%Info%Clear()
        deallocate(this%Info)
        nullify(This%Info)
    end if
    if (associated(this%Theory).and. .not. associated(this%Theory, Keep%Theory)) then
        call this%Theory%Clear()
        deallocate(this%Theory)
        nullify(This%Theory)
    end if
    end subroutine TCalculationAtParamPoint_Clear


    !!! TTheoryIntermediateCache

    subroutine TTheoryIntermediateCache_Clear(Info)
    class(TTheoryIntermediateCache) Info

    call Info%Initialize()

    end subroutine TTheoryIntermediateCache_Clear

    subroutine TTheoryIntermediateCache_Initialize(Info)
    class(TTheoryIntermediateCache) Info

    Info%validInfo = .false.

    end subroutine TTheoryIntermediateCache_Initialize


    subroutine TTheoryCalculator_VersionTraceOutput(this, ReadValues)
    class(TTheoryCalculator) :: this
    class(TNameValueList) :: ReadValues

    !Store for the record any useful info about version etc.

    end subroutine TTheoryCalculator_VersionTraceOutput


    subroutine TTheoryCalculator_AcceptReject(this, accpt, CurParams, Trial)
    !Handle freeing of memory of internal info: if accpt then clear CurParams, otherwise clear Trial
    class(TTheoryCalculator) :: this
    logical, intent(in) :: accpt
    class(TCalculationAtParamPoint), pointer :: CurParams, Trial

    if (.not. associated(CurParams%Info, Trial%Info)) then
        !If they point to same memory don't need to free anything
        if (accpt) then
            call CurParams%Clear(keep = Trial)
        else
            call Trial%Clear(keep = CurParams)
        end if
    end if

    end subroutine TTheoryCalculator_AcceptReject

    subroutine TTheoryCalculator_ErrorNotImplemented(this,S)
    class(TTheoryCalculator) :: this
    character(LEN=*), intent(in) :: S

    call MpiStop(trim(this%CalcName)//': Not implemented: '//trim(S))

    end subroutine TTheoryCalculator_ErrorNotImplemented


    subroutine TTheoryCalculator_Error(this,S)
    class(TTheoryCalculator) :: this
    character(LEN=*), intent(in) :: S

    call MpiStop(trim(this%CalcName)//': Error: '//trim(S))

    end subroutine TTheoryCalculator_Error


    !!TConfigClass
    subroutine TConfigClass_InitConfig(this, Config)
    class(TConfigClass) :: this
    class(TGeneralConfig), target :: config

    this%Config => config

    end subroutine TConfigClass_InitConfig


    subroutine TConfigClass_InitWithParams(this, Ini, Config)
    class(TConfigClass) :: this
    class(TIniFile) :: Ini
    class(TGeneralConfig), target :: Config

    call this%InitConfig(Config)

    end subroutine TConfigClass_InitWithParams


    !!TGeneralConfig
    subroutine SetTheoryParameterization(this, Ini, Names, defaultParam)
    class(TGeneralConfig) :: this
    class(TIniFile) :: Ini
    class(TParamNames) :: Names
    character(LEN=*), intent(in) :: defaultParam
    character(LEN=Ini_max_string_len) :: paramtxt

    paramtxt = Ini_Read_String_Default_File(Ini,'parameterization', defaultParam)
    if (.not. this%SetParameterizationName(paramtxt,Ini, Names)) then
        call MpiStop('GeneralConfig: unknown parameterization :'//trim(paramtxt))
    end if

    end subroutine SetTheoryParameterization

    function SetParameterizationName(this, nametag, Ini, Names) result(OK)
    class(TGeneralConfig) :: this
    character(LEN=*), intent(in) :: nametag
    class(TIniFile) :: Ini
    class(TParamNames) :: Names
    character(LEN=Ini_max_string_len) :: paramtxt
    logical OK
    Type(GenericParameterization), pointer :: GenParam

    OK = .true.

    if (nametag =='generic') then
        allocate(GenParam)
        this%Parameterization => GenParam
        call GenParam%Initialize(Ini,Names,'',this)
    else 
        OK = .false.
    end if

    end function SetParameterizationName


    end module GeneralTypes
