    module CosmologyConfig
    use GeneralTypes
    use Calculator_Cosmology
    use cmbtypes
    use CosmologyParameterizations

    Type, extends(TGeneralConfig) :: TCosmologyConfig
    contains
    procedure :: SetParameterizationName => Cosmo_SetParameterizationName

    end Type

    contains

    function Cosmo_SetParameterizationName(this, nametag, Ini, Names) result(OK)
    class(TCosmologyConfig) :: this
    character(LEN=*), intent(in) :: nametag
    class(TIniFile) :: Ini
    class(TParamNames) :: Names
    character(LEN=Ini_max_string_len) :: paramtxt
    logical OK
    Type(ThetaParameterization), pointer :: CMBParameterization
    Type(BackgroundParameterization), pointer :: BackgroundParam

    OK = .true.
    if (paramtxt =='background') then
        allocate(BackgroundParam)
        this%Parameterization => BackgroundParam
        call BackgroundParam%InitWithSetNames(Ini,Names,this)
    else if (paramtxt=='theta') then
        allocate(CMBParameterization)
        this%Parameterization => CMBParameterization
        call CMBParameterization%InitWithSetNames(Ini,Names,this)
    else
        OK =  this%TGeneralConfig%SetParameterizationName(nametag,Ini,Names)
    end if

    end function Cosmo_SetParameterizationName


    end module CosmologyConfig
