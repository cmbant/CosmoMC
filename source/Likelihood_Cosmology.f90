    module Likelihood_Cosmology
    use cmbtypes
    use Calculator_Cosmology
    implicit none
    private

    Type, extends(TCosmologyLikelihood) :: TCosmoCalcLikelihood
        class(TCosmologyCalculator), pointer :: Calculator => null()
    contains
    procedure :: InitWithCalculator => TCosmoCalcLikelihood_InitWithCalculator
    end type

    public TCosmoCalcLikelihood
    contains

    subroutine TCosmoCalcLikelihood_InitWithCalculator(this, Calc)
    class(TCosmoCalcLikelihood) :: this
    class(*), target :: Calc

    select type (Calc)
    class is (TCosmologyCalculator)
        this%Calculator => Calc
        class default 
            call MpiStop('TCosmoCalcLikelihood requires TCosmologyCalculator')
    end select

    end subroutine TCosmoCalcLikelihood_InitWithCalculator

    end module Likelihood_Cosmology