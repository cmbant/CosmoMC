    !HST from http://arxiv.org/abs/0905.0695
    !Thanks to Beth Reid, minor mods by AL, Oct 09

    !Updated by KM to use Riess et al (2011) value of H0 = 73.8 +/- 2.4 km/s/Mpc
    !Riess et al: 1103.2976

    module HST
    use cmbtypes
    use likelihood
    implicit none

    type, extends(CosmologyLikelihood) :: HSTLikelihood
    contains
    procedure :: LogLikeTheory => HST_LnLike
    end type HSTLikelihood

    ! angdistinveffh0 is the inverse of the angular diameter distance at z = 0.04 for H_0 = 74.2
    ! and a fiducial cosmology (omega_k = 0, omega_lambda = 0.7, w = -1); this is proportional to 
    !H_0 but includes the tiny cosmological dependence of the measurement (primarily on w) correctly.
    ! angdistinveffh0err = 3.6 / DL(0.04)
    !real(mcp), parameter :: angdistinvzeffh0 = 6.49405e-3, zeffh0 = 0.04, &
    !                       angdistinvzeffh0errsqr = 9.93e-8

    contains

    subroutine HSTLikelihood_Add(LikeList, Ini)
    class(LikelihoodList) :: LikeList
    Type(TIniFile) :: ini
    Type(HSTLikelihood), pointer :: like

    if (Ini_Read_Logical_File(Ini, 'use_HST',.false.)) then
        allocate(like)
        like%LikelihoodType = 'Hubble'
        like%name='HST'
        like%needs_background_functions = .true.
        call LikeList%Add(like)
    end if

    end subroutine HSTLikelihood_Add

    real(mcp) function HST_LnLike(like, CMB)
    use CAMB, only : AngularDiameterDistance  !!physical angular diam distance also in Mpc no h units
    use constants
    Class(HSTLikelihood) :: like
    Class(CMBParams) CMB
    real(mcp) :: theoryval

    real(mcp), parameter :: angdistinvzeffh0 = 6.45904e-3, zeffh0 = 0.04, &
    angdistinvzeffh0errsqr = 4.412e-8

    theoryval = 1.0/AngularDiameterDistance(real(zeffh0,dl))
    HST_LnLike = (theoryval - angdistinvzeffh0)**2/(2*angdistinvzeffh0errsqr)

    end function  HST_LnLike

    end module HST
