    module GeneralConfig
    use GeneralTypes
    use IniObjects
    implicit none

    Type TGeneralConfig
        class(TParameterization), pointer :: Parameterization => null()
        class(TTheoryCalculator), pointer :: TheoryCalculator => null()
    contains
    procedure :: SetTheoryParameterization
    procedure :: SetParameterizationName
    end Type

    contains

    end module GeneralConfig
