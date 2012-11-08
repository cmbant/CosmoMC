module DataLikelihoodList
    use likelihood
    use IniFile
    use ParamDef
    implicit none

    Type(LikelihoodList), target, save :: DataLikelihoods

    integer :: H0_min = 40, H0_max = 100
    real :: Use_min_zre = 0

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

!        Use_WeakLen = Ini_Read_Logical('use_WeakLen',.false.)
         Use_min_zre = Ini_Read_Double_File(Ini,'use_min_zre',0.d0) 
!        Use_Lya = Ini_Read_logical('use_lya',.false.)
       
!        if (Use_Lya .and. use_nonlinear) &
!             call DoAbort('Lya.f90 assumes LINEAR power spectrum input')

        !flag to force getting sigma8 even if not using LSS data 

 
end subroutine SetDataLikelihoods



end module DataLikelihoodList