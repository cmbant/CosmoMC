    module Calculator_Cosmology
    use CosmologyTypes
    use CosmoTheory
    use settings
    use GeneralTypes
    implicit none
    private

    real(mcp), parameter :: neutrino_mass_fac= 94.07_mcp !conversion factor for thermal with Neff=3 TCMB-2.7255
    !93.014 for 3.046
    real(mcp), parameter :: standard_neutrino_neff = 3.046_mcp

    Type TCosmologyImportanceOptions
        logical :: redo_cls, redo_pk
    end Type TCosmologyImportanceOptions

    Type, extends(TTheoryCalculator) :: TCosmologyCalculator
        Type(TCosmologyImportanceOptions) :: ImportanceOptions
    contains
    procedure :: BAO_D_v
    procedure :: Hofz
    procedure :: Hofz_Hunit
    procedure :: AngularDiameterDistance
    procedure :: ComovingRadialDistance
    procedure :: AngularDiameterDistance2
    procedure :: LuminosityDistance
    procedure :: CMBToTheta
    procedure :: GetNewBackgroundData
    procedure :: GetNewPowerData
    procedure :: GetNewTransferData
    procedure :: GetTheoryForImportance
    procedure :: GetCosmoTheoryForImportance
    procedure :: GetZreFromTau
    procedure :: GetOpticalDepth
    procedure :: SetBackgroundTheoryData
    procedure :: SetParamsForBackground
    procedure :: ReadImportanceParams => TCosmologyCalculator_ReadImportanceParams
    end Type TCosmologyCalculator

    public TCosmologyCalculator, TCosmologyImportanceOptions, neutrino_mass_fac, standard_neutrino_neff
    contains

    subroutine TCosmologyCalculator_ReadImportanceParams(this, Ini)
    class(TCosmologyCalculator) :: this
    class(TSettingIni) :: Ini

    this%ImportanceOptions%redo_cls= Ini%Read_Logical('redo_cls')
    this%ImportanceOptions%redo_pk= Ini%Read_Logical('redo_pk')

    end subroutine TCosmologyCalculator_ReadImportanceParams

    subroutine GetNewBackgroundData(this, CMB,Theory,error)
    class(TCosmologyCalculator) :: this
    class(CMBParams) :: CMB
    class(TCosmoTheoryPredictions) Theory
    integer error

    call this%SetParamsForBackground(CMB)
    select type (Param => this%Config%Parameterization)
    class is (TCosmologyParameterization)
        if (.not. Param%late_time_only) then
            call this%SetBackgroundTheoryData(CMB, Theory, error)
        end if
        class default
        call MpiStop('TCosmologyCalculator: Must have CosmologyParameterization')
    end select

    end subroutine GetNewBackgroundData

    subroutine SetParamsForBackground(this,CMB)
    class(TCosmologyCalculator) :: this
    class(CMBParams) CMB

    call this%ErrorNotImplemented('SetParamsForBackground')

    end subroutine SetParamsForBackground

    subroutine SetBackgroundTheoryData(this, CMB,Theory,error)
    class(TCosmologyCalculator) :: this
    class(CMBParams) :: CMB
    class(TCosmoTheoryPredictions) Theory
    integer error

    !calculate thermal history, e.g. z_drag etc.,
    !save derived parameters in Theory%derived_parameters, Theory%numderived
    call this%ErrorNotImplemented('SetBackgroundTheoryData')

    end subroutine SetBackgroundTheoryData

    subroutine GetNewTransferData(this, CMB, Info,Theory,error)
    class(TCosmologyCalculator) :: this
    class(CMBParams) :: CMB
    class(TTheoryIntermediateCache), pointer :: Info
    class(TCosmoTheoryPredictions) :: Theory
    integer error

    !gets transfer functions in Info, and also set any derived parameters
    call this%ErrorNotImplemented('GetNewTransferData')

    end subroutine GetNewTransferData

    subroutine GetNewPowerData(this, CMB, Info, Theory, error)
    class(TCosmologyCalculator) :: this
    class(CMBParams) :: CMB
    class(TTheoryIntermediateCache), pointer :: Info
    class(TCosmoTheoryPredictions)  :: Theory
    integer error

    !calculate Theory power spectra from transfer functions in Info
    call this%ErrorNotImplemented('GetNewPowerData')

    end subroutine GetNewPowerData

    subroutine GetTheoryForImportance(this, CMB, Theory, error)
    class(TCosmologyCalculator) :: this
    class(TTheoryParams) :: CMB
    class(TTheoryPredictions) :: Theory
    integer error

    select type (CMB)
    class is (CMBParams)
        select type (Theory)
        class is (TCosmoTheoryPredictions)
            call this%GetCosmoTheoryForImportance(CMB, Theory, error)
        end select
    end select

    end subroutine GetTheoryForImportance

    subroutine GetCosmoTheoryForImportance(this, CMB, Theory, error)
    class(TCosmologyCalculator) :: this
    class(CMBParams) :: CMB
    class(TCosmoTheoryPredictions) :: Theory
    integer error

    error=0
    !calculate power spectra from scratch (for importance sampling)
    !Theory may already be calculated, so only fill in missing bits (DoCls, DoPk) + derived
    call this%ErrorNotImplemented('GetTheoryForImportance')

    end subroutine GetCosmoTheoryForImportance

    function GetOpticalDepth(this,CMB)
    class(TCosmologyCalculator) :: this
    class(CMBParams) CMB
    real(mcp) GetOpticalDepth

    call this%ErrorNotImplemented('GetOpticalDepth')
    GetOpticalDepth = 0

    end  function GetOpticalDepth

    function GetZreFromTau(this, CMB, tau)
    class(TCosmologyCalculator) :: this
    class(CMBParams) CMB
    real(mcp), intent(in) :: tau
    real(mcp) GetZreFromTau

    call this%ErrorNotImplemented('GetZreFromTau')
    GetZreFromTau=0

    end function GetZreFromTau

    real(mcp) function CMBToTheta(this, CMB)
    class(TCosmologyCalculator) :: this
    class(CMBParams) CMB

    call this%ErrorNotImplemented('CMBToTheta')
    CMBToTheta=0

    end function CMBToTheta

    real(mcp) function BAO_D_v(this, z)
    class(TCosmologyCalculator) :: this
    real(mcp), intent(IN) :: z

    call this%ErrorNotImplemented('BAO_D_v')
    BAO_D_v = 0

    end function BAO_D_v

    real(mcp) function AngularDiameterDistance(this, z)
    class(TCosmologyCalculator) :: this
    real(mcp), intent(IN) :: z

    call this%ErrorNotImplemented('AngularDiameterDistance')
    AngularDiameterDistance = 0

    end function AngularDiameterDistance

    real(mcp) function ComovingRadialDistance(this, z)
    class(TCosmologyCalculator) :: this
    real(mcp), intent(IN) :: z

    call this%ErrorNotImplemented('ComovingRadialDistance')
    ComovingRadialDistance = 0

    end function ComovingRadialDistance

    real(mcp) function AngularDiameterDistance2(this, z1, z2)
    class(TCosmologyCalculator) :: this
    real(mcp), intent(IN) :: z1, z2

    call this%ErrorNotImplemented('AngularDiameterDistance2')
    AngularDiameterDistance2 = 0

    end function AngularDiameterDistance2

    real(mcp) function LuminosityDistance(this, z)
    class(TCosmologyCalculator) :: this
    real(mcp), intent(IN) :: z

    call this%ErrorNotImplemented('LuminosityDistance')
    LuminosityDistance = 0

    end function LuminosityDistance

    real(mcp) function Hofz(this, z)
    class(TCosmologyCalculator) :: this
    real(mcp), intent(IN) :: z

    call this%ErrorNotImplemented('Hofz')
    Hofz = 0

    end function Hofz

    real(mcp) function Hofz_Hunit(this, z)
    class(TCosmologyCalculator) :: this
    real(mcp), intent(IN) :: z

    Hofz_Hunit = const_c*this%Hofz(z)/1.d3

    end function Hofz_Hunit

    end module Calculator_Cosmology
