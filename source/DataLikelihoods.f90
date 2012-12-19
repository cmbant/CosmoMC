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

    !        if (Ini_Read_Logical_File(Ini, 'use_mpk',.false.)) &
    !           call DataLikelihoods%Add(MpkLikelihood()%Init(Ini)) 
    !        use_LSS=.true.
    call HSTLikelihood_Add(DataLikelihoods, Ini)

    call BAOLikelihood_Add(DataLikelihoods, Ini)

    call SNLikelihood_Add(DataLikelihoods, Ini)

    if(Ini_Read_Logical_File(Ini,'use_BBN',.false.)) &
    call DoAbort('Use_BBN not supported: use prior[omegabh2]=mean std')

    !        Use_Lya = Ini_Read_logical('use_lya',.false.)

    !        if (Use_Lya .and. use_nonlinear) &
    !             call DoAbort('Lya.f90 assumes LINEAR power spectrum input')


    end subroutine SetDataLikelihoods



    end module DataLikelihoodList