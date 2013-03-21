    module snovae
    implicit none

    contains

    subroutine SNLikelihood_Add(LikeList, Ini)
    use IniFile
    use likelihood
    use SNLS
    use Union2
    class(LikelihoodList) :: LikeList
    Type(TIniFile) :: ini
    integer count

    if (.not. Ini_Read_Logical_File(Ini, 'use_SN',.false.)) return

    count = LikeList%Count
    call SNLSLikelihood_Add(LikeList, Ini)
    call Union2Likelihood_Add(LikeList, Ini)
    if (LikeList%Count > count+1) call MpiStop('SNLikelihood_Add: more than one - datasets not independent')

    end subroutine SNLikelihood_Add

    end module snovae
