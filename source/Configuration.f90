    module Configuration
    use CalcLike
    use MonteCarlo
    use GeneralTypes

    Type TConfiguration
        class(TGeneralConfig), allocatable :: Config
        class(TLikeCalculator), allocatable :: Calculator
        class(TSampleCollector), allocatable :: SampleCollector
        class(TSamplingAlgorithm), allocatable :: TSamplingAlgorithm
    end type

    contains

    end module Configuration