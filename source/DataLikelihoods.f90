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
!<<<<<<< HEAD
    use szcounts !Anna
!=======
    use ElementAbundances
!>>>>>>> 1f3653b845a720b95e51f132f9b764ecd21f18fd
    class(TSettingIni), intent(in) :: Ini

    CosmoSettings%get_sigma8 = Ini%Read_Logical('get_sigma8',.false.)

    call CMBLikelihood_Add(DataLikelihoods, Ini)

    call AbundanceLikelihood_Add(DataLikelihoods, Ini)

    call HSTLikelihood_Add(DataLikelihoods, Ini)

    call SNLikelihood_Add(DataLikelihoods, Ini)

    call MPKLikelihood_Add(DataLikelihoods, Ini)

    if (use_mpk) call WiggleZLikelihood_Add(DataLikelihoods, Ini)

    call BAOLikelihood_Add(DataLikelihoods, Ini)

!<<<<<<< HEAD
    call SZLikelihood_Add(DataLikelihoods, Ini) !Anna
    CosmoSettings%use_SZ = Ini%Read_Logical('use_SZ',.false.)

    CosmoSettings%use_LSS = use_mpk

!=======
!>>>>>>> 1f3653b845a720b95e51f132f9b764ecd21f18fd
    end subroutine SetDataLikelihoods


    end module DataLikelihoodList
