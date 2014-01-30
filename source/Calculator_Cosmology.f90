    module Calculator_Cosmology
    use cmbtypes
    use settings
    use IniObjects
    use GeneralTypes
    implicit none
    private

    Type, extends(TTheoryCalculator) :: BaseCosmologyCalculator
    contains
    procedure :: BAO_D_v
    procedure :: CMBToTheta
    procedure :: GetNewBackgroundData
    procedure :: GetNewPowerData
    procedure :: GetNewTransferData
    procedure :: GetTheoryForImportance
    procedure :: GetZreFromTau
    procedure :: GetOpticalDepth
    procedure :: SetBackgroundTheoryData
    procedure :: SetParamsForBackground
    procedure :: ReadParams => CosmoCalc_ReadParams
    end Type BaseCosmologyCalculator

    public BaseCosmologyCalculator
    contains

    subroutine GetNewBackgroundData(this, CMB,Theory,error)
    class(BaseCosmologyCalculator) :: this
    class(CMBParams), intent(in) :: CMB
    class(CosmoTheoryPredictions) Theory
    integer error

    call this%SetParamsForBackground(CMB)
    select type (Param => this%Config%Parameterization)
    class is (TCosmologyParameterization)
        if (.not. Param%late_time_only) then
            call this%SetBackgroundTheoryData(CMB, Theory, error)
        end if
        class default
        call MpiStop('BaseCosmologyCalculator: Must have CosmologyParameterization')
    end select

    end subroutine GetNewBackgroundData

    subroutine SetParamsForBackground(this,CMB)
    class(BaseCosmologyCalculator) :: this
    class(CMBParams) CMB

    call this%ErrorNotImplemented('SetParamsForBackground')

    end subroutine SetParamsForBackground

    subroutine SetBackgroundTheoryData(this, CMB,Theory,error)
    class(BaseCosmologyCalculator) :: this
    class(CMBParams), intent(in) :: CMB
    class(CosmoTheoryPredictions) Theory
    integer error

    !calculate thermal history, e.g. z_drag etc.,
    !save derived parameters in Theory%derived_parameters, Theory%numderived
    call this%ErrorNotImplemented('SetBackgroundTheoryData')

    end subroutine SetBackgroundTheoryData

    subroutine GetNewTransferData(this, CMB, Info,Theory,error)
    class(BaseCosmologyCalculator) :: this
    class(CMBParams), intent(in) :: CMB
    class(TTheoryIntermediateCache), pointer :: Info
    class(CosmoTheoryPredictions) :: Theory
    integer error

    !gets transfer functions in Info, and also set any derived parameters
    call this%ErrorNotImplemented('GetNewTransferData')

    end subroutine GetNewTransferData

    subroutine GetNewPowerData(this, CMB, Info, Theory, error)
    class(BaseCosmologyCalculator) :: this
    class(CMBParams) :: CMB
    class(TTheoryIntermediateCache), pointer :: Info
    class(CosmoTheoryPredictions)  :: Theory
    integer error

    !calculate Theory power spectra from transfer functions in Info
    call this%ErrorNotImplemented('GetNewPowerData')

    end subroutine GetNewPowerData


    subroutine GetTheoryForImportance(this, CMB, Theory, error, DoCls, DoPk)
    class(BaseCosmologyCalculator) :: this
    class(CMBParams) :: CMB
    class(CosmoTheoryPredictions) :: Theory
    integer error
    logical, intent(in) :: DoCls, DoPk

    error=0
    !calculate power spectra from scratch (for importance sampling)
    !Theory may already be calculated, so only fill in missing bits (DoCls, DoPk) + derived
    call this%ErrorNotImplemented('GetTheoryForImportance')

    end subroutine GetTheoryForImportance

    function GetOpticalDepth(this,CMB)
    class(BaseCosmologyCalculator) :: this
    type(CMBParams) CMB
    real(mcp) GetOpticalDepth

    call this%ErrorNotImplemented('GetOpticalDepth')
    GetOpticalDepth = 0

    end  function GetOpticalDepth

    function GetZreFromTau(this, CMB, tau)
    class(BaseCosmologyCalculator) :: this
    class(CMBParams) CMB
    real(mcp), intent(in) :: tau
    real(mcp) GetZreFromTau

    call this%ErrorNotImplemented('GetZreFromTau')
    GetZreFromTau=0

    end function GetZreFromTau

    real(mcp) function CMBToTheta(this, CMB)
    class(BaseCosmologyCalculator) :: this
    class(CMBParams) CMB

    call this%ErrorNotImplemented('CMBToTheta')
    CMBToTheta=0

    end function CMBToTheta

    real(mcp) function BAO_D_v(this, z)
    class(BaseCosmologyCalculator) :: this
    real(mcp), intent(IN) :: z

    call this%ErrorNotImplemented('BAO_D_v')
    BAO_D_v = 0

    end function BAO_D_v

    subroutine CosmoCalc_ReadParams(this, Info, Ini)
    class(BaseCosmologyCalculator) :: this
    class(TTheoryIntermediateCache) Info
    class(TIniFile) :: Ini

    call Info%Initialize()

    end subroutine CosmoCalc_ReadParams



    end module Calculator_Cosmology