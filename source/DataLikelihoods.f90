    module DataLikelihoodList
    use likelihood
    use settings
    implicit none

    contains

    subroutine SetDataLikelihoods(Ini)
    use HST
    use snovae
    use cmbdata
    use bao
    use mpk
    use wigglez
    class(TSettingIni), intent(in) :: Ini

    get_sigma8 = Ini%Read_Logical('get_sigma8',.false.)

    call CMBDataLikelihoods_Add(DataLikelihoods, Ini)

    call HSTLikelihood_Add(DataLikelihoods, Ini)

    call SNLikelihood_Add(DataLikelihoods, Ini)

    call MPKLikelihood_Add(DataLikelihoods, Ini)
    
    if(use_mpk) call WiggleZLikelihood_Add(DataLikelihoods, Ini)

    call BAOLikelihood_Add(DataLikelihoods, Ini)

    use_LSS = use_mpk !.or. others

    end subroutine SetDataLikelihoods


    end module DataLikelihoodList