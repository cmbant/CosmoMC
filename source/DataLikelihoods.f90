module DataLikelihood
    use likelihood
    implicit none

    Type(LikelihoodList),save :: DataLikelihoods

contains

subroutine SetDataLikelihoods(Ini)
 use IniFile
 use HST
 use cmbdata
 logical useItem
 
        Use_CMB = Ini_Read_Logical_File(DefIni, 'use_CMB',.true.)
        if (Use_CMB) then
           call DataLikelihoods%Add(CMBDataLikelihood()%Init(DefIni)) 
        end if

        if (Ini_Read_Logical_File(Ini, 'use_mpk',.false.)) &
           call DataLikelihoods%Add(MpkLikelihood()%Init(DefIni)) 

        if (Ini_Read_Logical_File(Ini, 'use_HST',.false.)) &
           call DataLikelihoods%Add(HSTLikelihood()%Init(DefIni)) 

        if(Ini_Read_Logical_File(Ini,'use_BBN',.false.)) &
          call DoAbort('Use_BBN not supported: use prior[omegabh2]=mean std')

        if (Ini_Read_Logical_File(Ini, 'use_SN',.false.)) &
           call DataLikelihoods%Add(SNLikelihood()%Init(DefIni)) 

        if (Ini_Read_Logical_File(Ini, 'use_BAO',.false.)) &
           call DataLikelihoods%Add(BAOLikelihood()%Init(DefIni)) 
        

        Use_WeakLen = Ini_Read_Logical('use_WeakLen',.false.)
        Use_min_zre = Ini_Read_Double('use_min_zre',0.d0) 
        Use_Lya = Ini_Read_logical('use_lya',.false.)
       
        if (Use_Lya .and. use_nonlinear) &
             call DoAbort('Lya.f90 assumes LINEAR power spectrum input')

        !flag to force getting sigma8 even if not using LSS data 
        use_LSS = Ini_Read_Logical('get_sigma8',.false.)
        ! use_LSS = Use_2dF .or. Use_Clusters .or. Use_WeakLen
        use_LSS = Use_LSS .or. Use_mpk .or. Use_Clusters .or. Use_WeakLen .or. Use_Lya

 
end subroutine SetDataLikelihoods



end module DataLikelihoods