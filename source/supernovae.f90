    module snovae
    implicit none

    contains

    subroutine SNLikelihood_Add(LikeList, Ini)
    use SNLS
    use Union2
    use likelihood
    use settings
    use JLA
    class(TLikelihoodList) :: LikeList
    class(TSettingIni) :: ini
    integer count

    if (.not. Ini%Read_Logical('use_SN',.false.)) return

    count = LikeList%Count
    call SNLSLikelihood_Add(LikeList, Ini)
    call JLALikelihood_Add(LikeList, Ini)
    call Union2Likelihood_Add(LikeList, Ini)
    if (LikeList%Count > count+1) call MpiStop('SNLikelihood_Add: more than one - datasets not independent')

    end subroutine SNLikelihood_Add

    end module snovae
