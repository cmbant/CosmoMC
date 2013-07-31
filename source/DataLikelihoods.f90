    module DataLikelihoodList
    use likelihood
    use IniFile
    use ParamDef
    implicit none

    contains

    subroutine SetDataLikelihoods(Ini)
    use HST
    use snovae
    use cmbdata
    use bao
    Type(TIniFile), intent(in) :: Ini

    use_LSS = Ini_Read_Logical_File(Ini,'get_sigma8',.false.)

    call CMBDataLikelihoods_Add(DataLikelihoods, Ini)

    call HSTLikelihood_Add(DataLikelihoods, Ini)

    call BAOLikelihood_Add(DataLikelihoods, Ini)

    call SNLikelihood_Add(DataLikelihoods, Ini)

    end subroutine SetDataLikelihoods


    end module DataLikelihoodList