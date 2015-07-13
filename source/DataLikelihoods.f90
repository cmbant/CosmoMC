    module DataLikelihoodList
    use likelihood
    use settings
    use CosmologyTypes
    implicit none

    contains

    subroutine SetDataLikelihoods(Ini)
    use HST
    use snovae
    use CMBLikelihoods
    use bao
    use mpk
    use wigglez
    use szcounts !Anna
    use wl
    use ElementAbundances
    class(TSettingIni), intent(in) :: Ini

    CosmoSettings%get_sigma8 = Ini%Read_Logical('get_sigma8',.false.)

    call CMBLikelihood_Add(DataLikelihoods, Ini)

    call AbundanceLikelihood_Add(DataLikelihoods, Ini)

    call HSTLikelihood_Add(DataLikelihoods, Ini)

    call SNLikelihood_Add(DataLikelihoods, Ini)

    call MPKLikelihood_Add(DataLikelihoods, Ini)

    if (use_mpk) call WiggleZLikelihood_Add(DataLikelihoods, Ini)

    call BAOLikelihood_Add(DataLikelihoods, Ini)

    call SZLikelihood_Add(DataLikelihoods, Ini) !Anna

    call WLLikelihood_Add(DataLikelihoods, Ini)

    end subroutine SetDataLikelihoods


    end module DataLikelihoodList
