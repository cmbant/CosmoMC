    module Likelihood_Cosmology
    use cmbtypes
    use Calculator_Cosmology
    use CosmoTheory
    implicit none
    private

    Type, extends(TCosmologyRequirementsLikelihood) :: TCosmologyLikelihood
        !This just type cases the arguments of the generic GetLogLlike function to the specific cosmology types
        !for use by derived classes (which override LogLike etc, not GetLogLike)
    contains
    procedure :: LogLikeDataParams !same as above when extra theory info not needed
    procedure :: LogLikeTheory !same as above when extra theory and nuisance info not needed
    procedure :: LogLike
    procedure :: GetLogLike
    end type

    Type, extends(TCosmologyLikelihood) :: TCosmoCalcLikelihood
        !Just convenience class for likelihoods which use a CosmologyCalculator
        class(TCosmologyCalculator), pointer :: Calculator => null()
    contains
    procedure :: InitConfig => TCosmoCalcLikelihood_InitConfig
    end type

    public TCosmoCalcLikelihood, TCosmologyLikelihood
    contains


    function GetLogLike(like, Params, Theory, DataParams) result(LogLike)
    class(TCosmologyLikelihood) :: like
    class(TTheoryParams) :: Params
    class(TTheoryPredictions) :: Theory
    real(mcp) :: DataParams(:)
    real(mcp) LogLike

    select type(Params)
    class is (CMBParams)
        select type(Theory)
        class is (TCosmoTheoryPredictions)
            logLike = like%LogLike(Params,Theory,DataParams)
            class default 
            call MpiStop('TCosmologyLikelihood requires TCosmoTheoryPredictions')
        end select
        class default 
        call MpiStop('TCosmologyLikelihood requires CMBParams')
    end select

    end function GetLogLike


    function logLikeTheory(like,CMB) result(LogLike)
    !For likelihoods that don't need Theory or DataParams
    class(TCosmologyLikelihood) :: like
    class(CMBParams) :: CMB
    real(mcp) LogLike

    LogLike= logZero
    stop 'logLikeTheory or logLike should not be overridden'

    end function logLikeTheory

    function LogLikeDataParams(like, CMB, DataParams) result(LogLike)
    class(TCosmologyLikelihood) :: like
    class(CMBParams) :: CMB
    real(mcp) :: DataParams(:)
    real(mcp) LogLike

    LogLike = like%logLikeTheory(CMB)
    end function LogLikeDataParams

    function LogLike(like, CMB, Theory, DataParams)
    class(TCosmologyLikelihood) :: like
    class(CMBParams) :: CMB
    class(TCosmoTheoryPredictions), target :: Theory
    real(mcp) :: DataParams(:)
    real(mcp) LogLike

    logLike = like%LogLikeDataParams(CMB, DataParams)
    end function LogLike


    subroutine TCosmoCalcLikelihood_InitConfig(this, Config)
    class(TCosmoCalcLikelihood) :: this
    class(TGeneralConfig), target :: Config

    call this%TCosmologyLikelihood%InitConfig(Config)

    select type (Calc => this%Calculator)
    class is (TCosmologyCalculator)
        this%Calculator => Calc
        class default 
        call MpiStop('TCosmoCalcLikelihood requires TCosmologyCalculator')
    end select

    end subroutine TCosmoCalcLikelihood_InitConfig

    end module Likelihood_Cosmology