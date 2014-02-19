    module Likelihood_Cosmology
    use cmbtypes
    use Calculator_Cosmology
    implicit none
    private
    
    Type, extends(TCosmologyLikelihood) :: TCosmoCalcLikelihood
        class(TCosmologyCalculator), pointer :: Calculator => null()
    end type

    public TCosmoCalcLikelihood
    contains

    end module Likelihood_Cosmology