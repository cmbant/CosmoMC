    module GeneralTypes
    use settings

    Type TTheoryParams
        real(mcp) :: BaseParams(max_theory_params) = 0._mcp
    end Type TTheoryParams

    Type TTheoryPredictions
    end Type TTheoryPredictions

    type :: TParameterization
    contains
    procedure :: ParamArrayToTheoryParams
    procedure :: CalcDerivedParams=> TParameterization_CalcDerivedParams
    procedure :: Init => TParameterization_init
    procedure :: NonBaseParameterPriors => TParameterization_NonBaseParameterPriors
    end type TParameterization

    contains

    subroutine ParamArrayToTheoryParams(this, Params, CMB)
    class(TParameterization) :: this
    real(mcp) Params(:)
    Class(TTheoryParams), target :: CMB

    CMB%BaseParams(1:num_theory_params) = Params(1:num_theory_params)

    end subroutine ParamArrayToTheoryParams

    subroutine TParameterization_Init(this, Ini, Names, DefaultName)
    Class(TParameterization) :: this
    Type(TIniFile) :: Ini
    Type(TParamNames) :: Names
    character(LEN=*), intent(in) :: DefaultName
    character(LEN=Ini_max_string_len) :: ParamNamesFile = ''

    ParamNamesFile = ReadIniFileName(Ini,'ParamNamesFile', NotFoundFail=.false.)

    if (ParamNamesFile /='') then
        call ParamNames_init(Names, ParamNamesFile)
    else
        if (generic_mcmc) then
            Names%nnames=0
            if (Feedback>0) write (*,*) 'Change the Parameterization type implementation to use Parameter names'
        else
            call ParamNames_init(Names, trim(LocalDir)//trim(DefaultName))
        end if
    end if
    end subroutine TParameterization_Init

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


    end module GeneralTypes
