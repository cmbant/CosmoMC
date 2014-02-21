
    module GeneralTypes
    use settings
    use likelihood
    use ObjectLists
    use IO
    implicit none
    private

    type int_arr_pointer
        integer, dimension(:), pointer :: p
    end type int_arr_pointer

    type, extends(TSaveLoadStateObject) :: TCheckpointable
    contains
    procedure :: ReadParams => TCheckpointable_ReadParams
    end Type

    type TObjectWithParams
    contains
    procedure :: ReadParams => TObjectWithParams_ReadParams
    end Type

    Type TTheoryParams
        real(mcp) :: BaseParams(max_theory_params) = 0._mcp
    end Type TTheoryParams

    Type TTheoryIntermediateCache
        !Cache of intermediate steps in the theory calculation that may be (partly) reused for other points
    contains
    procedure :: Clear => TTheoryIntermediateCache_Clear
    end Type TTheoryIntermediateCache

    Type TCalculationAtParamPoint
        !Parameter values, calculated theory and intermediates and likelihoods at a particular point in parameter space
        real(mcp) :: P(max_num_params)
        real(mcp) :: likelihoods(max_likelihood_functions)
        real(mcp) :: lastParamArray(max_num_params)
        logical :: validInfo = .false.
        class(TTheoryPredictions), pointer :: Theory => null()
        class(TTheoryIntermediateCache), pointer :: Info => null()
    contains
    procedure :: Clear => TCalculationAtParamPoint_Clear
    procedure :: WriteParams => TCalculationAtParamPoint_WriteParams
    procedure :: AcceptReject => TCalculationAtParamPoint_AcceptReject
    procedure :: ReadModel
    procedure :: WriteModel
    end Type TCalculationAtParamPoint

    !Each descendant of the base TConfigClass below has a Config element which point to a TGeneralConfig instance taht determines
    !which related class implementations are actually used
    Type, extends(TObjectWithParams) :: TGeneralConfig
        class(TParameterization), pointer :: Parameterization => null()
        class(TTheoryCalculator), pointer :: Calculator => null()
    contains
    procedure :: SetTheoryParameterization => TGeneralConfig_SetTheoryParameterization
    procedure :: SetParameterizationName => TGeneralConfig_SetParameterizationName
    procedure :: NewTheory => TGeneralConfig_NewTheory
    procedure :: InitForLikelihoods => TGeneralConfig_InitForLikelihoods
    end Type

    Type, extends(TCheckpointable) :: TConfigClass
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
    procedure :: AssignNew => TTheoryPredictions_AssignNew
    procedure :: Clear => TTheoryPredictions_Clear
    procedure :: WriteTheory
    procedure :: ReadTheory
    procedure :: WriteBestFitData
    end Type TTheoryPredictions

    abstract interface

    subroutine WriteTheory(T, unit)
    import TTheoryPredictions
    class(TTheoryPredictions) T
    integer, intent(in) :: unit
    end subroutine WriteTheory

    subroutine ReadTheory(T, unit)
    import TTheoryPredictions
    class(TTheoryPredictions) T
    integer, intent(in) :: unit
    end subroutine ReadTheory

    subroutine WriteBestFitData(Theory,fnameroot)
    import TTheoryPredictions
    class(TTheoryPredictions) Theory
    character(LEN=*), intent(in) :: fnameroot
    end subroutine WriteBestFitData
    end interface


    Type, extends(TConfigClass) :: TTheoryCalculator
        character(LEN=128) :: calcName = 'TheoryCalculator'
    contains
    procedure :: InitForLikelihoods => TTheoryCalculator_InitForLikelihoods
    procedure :: Error => TTheoryCalculator_Error
    procedure :: ErrorNotImplemented => TTheoryCalculator_ErrorNotImplemented
    procedure :: VersionTraceOutput => TTheoryCalculator_VersionTraceOutput
    procedure :: ReadImportanceParams => TTheoryCalculator_ReadImportanceParams
    procedure :: GetTheoryForImportance => TTheoryCalculator_GetTheoryForImportance
    end Type TTheoryCalculator

    type, extends(TConfigClass)  :: TParameterization
    contains
    procedure :: ParamArrayToTheoryParams
    procedure :: CalcDerivedParams=> TParameterization_CalcDerivedParams
    procedure :: Initialize => TParameterization_Initialize
    procedure :: NonBaseParameterPriors => TParameterization_NonBaseParameterPriors
    procedure :: NewTheoryParams => TParameterization_NewTheoryParams
    end type TParameterization

    Type, extends(TParameterization) :: GenericParameterization
    contains
    end type GenericParameterization


    abstract interface
    subroutine WriteModel(this, unit, like, mult)
    import TCalculationAtParamPoint, mcp
    Class(TCalculationAtParamPoint) :: this
    integer unit
    real(mcp), intent(in) :: mult, like
    end subroutine WriteModel

    subroutine ReadModel(this,  unit, has_likes, mult, like, error)
    import TCalculationAtParamPoint, mcp
    Class(TCalculationAtParamPoint) :: this
    integer, intent(in) :: unit
    integer, intent(out) :: error
    real(mcp), intent(out) :: mult, like
    logical, intent(out) :: has_likes(:)
    end subroutine ReadModel

    end interface

    public int_arr_pointer,TCheckpointable, TTheoryParams, TTheoryIntermediateCache, TCalculationAtParamPoint, TGeneralConfig, &
    & TConfigClass, TTheoryPredictions, TTheoryCalculator, TParameterization, GenericParameterization

    contains


    subroutine ParamArrayToTheoryParams(this, Params, CMB)
    class(TParameterization) :: this
    real(mcp) Params(:)
    Class(TTheoryParams), target :: CMB

    CMB%BaseParams(1:num_theory_params) = Params(1:num_theory_params)

    end subroutine ParamArrayToTheoryParams

    subroutine TParameterization_Initialize(this, Ini, Names, DefaultName, Config)
    class(TParameterization) :: this
    class(TSettingIni) :: Ini
    class(TParamNames) :: Names
    character(LEN=*), intent(in) :: DefaultName
    class(TGeneralConfig), target :: Config
    character(LEN=:), allocatable :: ParamNamesFile

    call this%TConfigClass%Init(Config)

    ParamNamesFile = Ini%ReadFileName('ParamNamesFile', NotFoundFail=.false.)
    if (ParamNamesFile =='' .and. DefaultName/='') ParamNamesFile= trim(LocalDir)//trim(DefaultName)

    if (ParamNamesFile /='') then
        call Names%init(ParamNamesFile)
        num_theory_params= Names%num_MCMC
    else
        Names%nnames=0
        num_theory_params= Ini%Read_Int('num_theory_params')
    end if
    if (num_theory_params> max_theory_params) call MpiStop('see settings.f90: num_theory_params> max_theory_params')
    index_data =  num_theory_params+1

    end subroutine TParameterization_Initialize


    function TParameterization_NonBaseParameterPriors(this,CMB)
    class(TParameterization) :: this
    class(TTheoryParams) :: CMB
    real(mcp):: TParameterization_NonBaseParameterPriors

    TParameterization_NonBaseParameterPriors = 0

    end function TParameterization_NonBaseParameterPriors

    function TParameterization_CalcDerivedParams(this, P, Theory, derived) result (num_derived)
    class(TParameterization) :: this
    Type(mc_real_pointer) :: derived
    class(TTheoryPredictions), pointer :: Theory !can be null for simple cases (e.g. generic)
    real(mcp) :: P(:)
    integer num_derived

    num_derived = 0
    allocate(Derived%P(num_derived))

    end function TParameterization_CalcDerivedParams

    subroutine TParameterization_NewTheoryParams(this,TheoryParams)
    class(TParameterization) :: this
    class(TTheoryParams), allocatable :: TheoryParams

    allocate(TTheoryParams::TheoryParams)

    end subroutine TParameterization_NewTheoryParams


    !!!TTheoryPredictions

    subroutine TTheoryPredictions_AssignNew(this, NewTheory)
    class(TTheoryPredictions) :: this
    class(TTheoryPredictions), pointer :: NewTheory

    allocate(NewTheory, source = this)

    end subroutine TTheoryPredictions_AssignNew

    subroutine TTheoryPredictions_Clear(this)
    class(TTheoryPredictions) :: this
    end subroutine TTheoryPredictions_Clear

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


    subroutine TCalculationAtParamPoint_AcceptReject(this, Trial, accpt)
    !Handle freeing of memory of internal info: if accpt then clear this, otherwise clear Trial
    class(TCalculationAtParamPoint) :: this
    class(TCalculationAtParamPoint) :: Trial
    logical, intent(in) :: accpt

    if (.not. associated(this%Info, Trial%Info)) then
        !If they point to same memory don't need to free anything
        if (accpt) then
            call this%Clear(keep = Trial)
        else
            call Trial%Clear(keep = this)
        end if
    end if

    end subroutine TCalculationAtParamPoint_AcceptReject


    subroutine TCalculationAtParamPoint_WriteParams(this, Config, mult, like)
    class(TCalculationAtParamPoint) this
    class(TGeneralConfig) :: Config
    real(mcp), intent(in) :: mult, like
    real(mcp), allocatable :: output_array(:)
    Type(mc_real_pointer) :: derived
    integer :: numderived = 0

    if (outfile_handle ==0) return

    numderived = Config%Parameterization%CalcDerivedParams(this%P, this%Theory, derived)
    numderived = DataLikelihoods%addLikelihoodDerivedParams(this%P, this%Theory, derived)

    if (numderived==0) then
        call IO_OutputChainRow(outfile_handle, mult, like, this%P(params_used), num_params_used)
    else
        allocate(output_array(num_params_used + numderived))
        output_array(1:num_params_used) =  this%P(params_used)
        output_array(num_params_used+1:num_params_used+numderived) =  derived%P
        call IO_OutputChainRow(outfile_handle, mult, like, output_array)
        deallocate(output_array)
    end if
    deallocate(derived%P)

    end subroutine TCalculationAtParamPoint_WriteParams


    !!! TTheoryIntermediateCache

    subroutine TTheoryIntermediateCache_Clear(Info)
    class(TTheoryIntermediateCache) Info

    end subroutine TTheoryIntermediateCache_Clear

    !!! TTheoryCalculator

    subroutine TTheoryCalculator_InitForLikelihoods(this)
    class(TTheoryCalculator) :: this

    !Called after likelihoods etc loaded
    end subroutine TTheoryCalculator_InitForLikelihoods


    subroutine TTheoryCalculator_ReadImportanceParams(this, Ini)
    class(TTheoryCalculator) :: this
    class(TSettingIni) :: Ini

    end subroutine TTheoryCalculator_ReadImportanceParams

    subroutine TTheoryCalculator_GetTheoryForImportance(this, CMB, Theory, error)
    class(TTheoryCalculator) :: this
    class(*) :: CMB
    class(*) :: Theory
    integer error

    error=0
    !calculate power spectra from scratch (for importance sampling)
    !Theory may already be calculated, so only fill in missing bits (DoCls, DoPk) + derived
    call this%ErrorNotImplemented('GetTheoryForImportance')

    end subroutine TTheoryCalculator_GetTheoryForImportance


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
    class(TSettingIni) :: Ini
    class(TGeneralConfig), target :: Config

    call this%InitConfig(Config)
    call this%ReadParams(Ini)

    end subroutine TConfigClass_InitWithParams


    !!TGeneralConfig
    subroutine TGeneralConfig_SetTheoryParameterization(this, Ini, Names, defaultParam)
    class(TGeneralConfig) :: this
    class(TSettingIni) :: Ini
    class(TParamNames) :: Names
    character(LEN=*), intent(in) :: defaultParam
    character(LEN=:), allocatable :: paramtxt

    paramtxt = Ini%Read_String_Default('parameterization', defaultParam)
    if (.not. this%SetParameterizationName(paramtxt,Ini, Names)) then
        call MpiStop('GeneralConfig: unknown parameterization :'//trim(paramtxt))
    end if

    end subroutine TGeneralConfig_SetTheoryParameterization

    function TGeneralConfig_SetParameterizationName(this, nametag, Ini, Names) result(OK)
    class(TGeneralConfig) :: this
    character(LEN=*), intent(in) :: nametag
    class(TSettingIni) :: Ini
    class(TParamNames) :: Names
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

    end function TGeneralConfig_SetParameterizationName

    function TGeneralConfig_NewTheory(this) result(Theory)
    class(TGeneralConfig) :: this
    class(TTheoryPredictions), pointer :: Theory

    allocate(TTheoryPredictions::Theory)
    call Theory%Init(this)

    end function TGeneralConfig_NewTheory

    subroutine TGeneralConfig_InitForLikelihoods(this)
    class(TGeneralConfig) :: this
    class(TDataLikelihood), pointer :: DataLike
    integer i

    if (associated(this%Calculator)) then
        call this%Calculator%InitForLikelihoods()
        do i=1,DataLikelihoods%Count
            DataLike=>DataLikelihoods%Item(i)
            call DataLike%InitWithCalculator(this%Calculator)
        end do
    end if

    end subroutine TGeneralConfig_InitForLikelihoods

    !!!TCheckpointable

    subroutine TCheckpointable_ReadParams(this, Ini)
    class(TCheckpointable) :: this
    class(TSettingIni) :: Ini
    end subroutine TCheckpointable_ReadParams

    !!!TObjectWithParams

    subroutine TObjectWithParams_ReadParams(this, Ini)
    class(TObjectWithParams) :: this
    class(TSettingIni) :: Ini
    end subroutine TObjectWithParams_ReadParams

    end module GeneralTypes