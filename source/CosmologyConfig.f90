    module CosmologyConfig
    use GeneralTypes
    use CalcLike
    use CosmologyTypes
    use CosmoTheory
    use Calculator_Cosmology
    use CosmologyParameterizations
    use CalcLike_Cosmology
    use SampleCollector
    use GeneralSetup
    use Likelihood_Cosmology
    implicit none
    private

    Type, extends(TGeneralConfig) :: TCosmologyConfig
    contains
    procedure :: SetParameterizationName => TCosmologyConfig_SetParameterizationName
    procedure :: NewTheory => TCosmologyConfig_NewTheory
    procedure :: InitForLikelihoods => TCosmologyConfig_InitForLikelihoods
    procedure :: ReadParams => TCosmologyConfig_ReadParams
    end Type

    type, extends(TSetup) :: TCosmologySetup
    contains
    procedure :: Init => TCosmologySetup_Init
    end type

    public TCosmologyConfig, TCosmologySetup
    contains

    subroutine TCosmologyConfig_ReadParams(this, Ini)
    use Calculator_CAMB
#ifdef PICO
    use Calculator_PICO
#endif
    class(TCosmologyConfig) :: this
    class(TSettingIni) :: Ini
    character(LEN=:), allocatable :: CalcName

    CalcName = Ini%Read_String_Default('cosmology_calculator', 'CAMB')
    if (calcName=='CAMB') then
        allocate(CAMB_Calculator::this%Calculator)
    else if (calcName=='PICO') then
#ifdef PICO
        allocate(PICO_Calculator::this%Calculator)
#else
        call MpiStop('Compile with -DPICO to use PICO calculator (e.g. make cosmomc PICO=/path/to/pic/pypico)')
#endif
    else
        call MpiStop('Calculator not supported: ' //CalcName)
    end if
    call this%Calculator%InitWithParams(Ini,this)

    call CosmoSettings%ReadParams(Ini)

    end subroutine TCosmologyConfig_ReadParams

    function TCosmologyConfig_SetParameterizationName(this, nametag, Ini, Names) result(OK)
    class(TCosmologyConfig) :: this
    character(LEN=*), intent(in) :: nametag
    class(TSettingIni) :: Ini
    class(TParamNames) :: Names
    logical OK
    Type(ThetaParameterization), pointer :: CMBParameterization
    Type(BackgroundParameterization), pointer :: BackgroundParam

    OK = .true.
    if (nametag =='background') then
        allocate(BackgroundParam)
        this%Parameterization => BackgroundParam
        call BackgroundParam%InitWithSetNames(Ini,Names,this)
    else if (nametag=='theta') then
        allocate(CMBParameterization)
        this%Parameterization => CMBParameterization
        call CMBParameterization%InitWithSetNames(Ini,Names,this)
    else
        OK =  this%TGeneralConfig%SetParameterizationName(nametag,Ini,Names)
    end if

    end function TCosmologyConfig_SetParameterizationName


    subroutine TCosmologyConfig_NewTheory(this, Theory)
    class(TCosmologyConfig) :: this
    class(TTheoryPredictions), allocatable :: Theory

    allocate(TCosmoTheoryPredictions::Theory)
    call Theory%Init(this)
    select type(Theory)
    class is (TCosmoTheoryPredictions)
        call Theory%AllocateForSettings(CosmoSettings)
    end select

    end subroutine TCosmologyConfig_NewTheory

    subroutine TCosmologyConfig_InitForLikelihoods(this)
    class(TCosmologyConfig) :: this

    call CosmoSettings%InitForLikelihoods()
    call this%TGeneralConfig%InitForLikelihoods()

    end subroutine TCosmologyConfig_InitForLikelihoods


    !!!TCosmologySetup

    subroutine TCosmologySetup_Init(this)
    class(TCosmologySetup) :: this

    allocate(TCosmologyConfig::this%Config)
    allocate(TCosmoLikeCalculator::this%LikeCalculator)
    call this%TSetup%Init()

    end subroutine TCosmologySetup_Init

    end module CosmologyConfig
