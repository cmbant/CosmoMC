    module CosmologyConfig
    use GeneralTypes
    use CalcLike
    use Calculator_Cosmology
    use cmbtypes
    use CosmologyParameterizations
    implicit none

    Type, extends(TGeneralConfig) :: TCosmologyConfig
    contains
    procedure :: SetParameterizationName => TCosmologyConfig_SetParameterizationName
    procedure :: NewTheory => TCosmologyConfig_NewTheory
    end Type

    contains

    subroutine TCosmologyConfig_ReadParams(this, Ini)
    use Calculator_CAMB
    class(TCosmologyConfig) :: this
    class(TSettingIni) :: Ini
    character(LEN=:), allocatable :: CalcName

    CalcName = Ini%Read_String_Default('cosmology_calculator', 'CAMB')
    if (calcName=='CAMB') then
        allocate(CAMB_Calculator::this%Calculator)
    else
        call MpiStop('only CAMB currently supported')
    end if
    call this%Calculator%Init(this)

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

    function TCosmologyConfig_NewTheory(this) result(Theory)
    class(TCosmologyConfig) :: this
    class(TTheoryPredictions), pointer :: Theory

    allocate(CosmoTheoryPredictions::Theory)
    call Theory%Init(this)

    end function TCosmologyConfig_NewTheory


    end module CosmologyConfig
