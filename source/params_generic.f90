    module DefineParameterization
    use GeneralTypes
    implicit none

    Type, extends(TParameterization) :: GenericParameterization
    contains
    procedure :: CalcDerivedParams => Gen_CalcDerivedParams
    procedure :: Initialize => Gen_Initialize
    end type GenericParameterization

    contains

    subroutine SetTheoryParameterization(Ini, Names)
    Type(TIniFile) :: Ini
    Type(TParamNames) :: Names
    character(LEN=Ini_max_string_len) :: paramtxt
    Type(GenericParameterization), pointer :: GenParam

    paramtxt = Ini_Read_String_Default_File(Ini,'parameterization', 'generic')
    if (paramtxt =='generic') then
        allocate(GenParam)
        Parameterization => GenParam
        call GenParam%Initialize(Ini,Names)
    else
        call MpiStop('params_generic: unknown parameterization :'//trim(paramtxt))
    end if

    end subroutine SetTheoryParameterization

    subroutine Generic_Initialize(this, Ini, Names)
    class(TParameterization) :: this
    Type(TIniFile) :: Ini
    Type(TParamNames) :: Names

    Names%nnames=0
    call this%TParameterization%Init(Ini,Names,'')
    if (Names%nnames==0) then
        num_theory_params= Ini_Read_Int_File(Ini, 'num_theory_params')
        if (Feedback>0) write (*,*) 'Change the params_generic type implementation to use parameter names'
    else
        !unnamed parameters (just numered)
        num_theory_params= Names%num_MCMC
    end if

    end subroutine Generic_Initialize

    function Gen_CalcDerivedParams(this, P, Theory, derived) result (num_derived)
    class(GenericParameterization) :: this
    Type(mc_real_pointer) :: derived
    class(TTheoryPredictions) :: Theory
    real(mcp) :: P(:)
    integer num_derived

    !If you want to output any derived parameters for each chain position, define them here
    num_derived = 0
    allocate(Derived%P(num_derived))

    end function Gen_CalcDerivedParams

    end module DefineParameterization

