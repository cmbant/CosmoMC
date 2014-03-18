    !HST from http://arxiv.org/abs/0905.0695
    !Thanks to Beth Reid, minor mods by AL, Oct 09

    !Updated by KM to use Riess et al (2011) value of H0 = 73.8 +/- 2.4 km/s/Mpc
    !Riess et al: 1103.2976

    module HST
    use CosmologyTypes
    use Likelihood_Cosmology
    implicit none
    private

    type, extends(TCosmoCalcLikelihood) :: HSTLikelihood
    contains
    procedure :: LogLikeTheory => HST_LnLike
    end type HSTLikelihood

    ! angdistinveffh0 is the inverse of the angular diameter distance at z = 0.04 for H_0 = 74.2
    ! and a fiducial cosmology (omega_k = 0, omega_lambda = 0.7, w = -1); this is proportional to
    !H_0 but includes the tiny cosmological dependence of the measurement (primarily on w) correctly.
    ! angdistinveffh0err = 3.6 / DL(0.04)
    !real(mcp), parameter :: angdistinvzeffh0 = 6.49405e-3, zeffh0 = 0.04, &
    !                       angdistinvzeffh0errsqr = 9.93e-8

    public HSTLikelihood, HSTLikelihood_Add
    contains

    subroutine HSTLikelihood_Add(LikeList, Ini)
    class(TLikelihoodList) :: LikeList
    class(TSettingIni) :: ini
    Type(HSTLikelihood), pointer :: this

    if (Ini%Read_Logical('use_HST',.false.)) then
        allocate(this)
        this%LikelihoodType = 'Hubble'
        this%name='HST'
        this%needs_background_functions = .true.
        call LikeList%Add(this)
    end if

    end subroutine HSTLikelihood_Add

    real(mcp) function HST_LnLike(this, CMB)
    Class(HSTLikelihood) :: this
    Class(CMBParams) CMB
    real(mcp) :: theoryval

    real(mcp), parameter :: angdistinvzeffh0 = 6.45904e-3, zeffh0 = 0.04, &
    angdistinvzeffh0errsqr = 4.412e-8

    theoryval = 1.0/this%Calculator%AngularDiameterDistance(zeffh0)
    HST_LnLike = (theoryval - angdistinvzeffh0)**2/(2*angdistinvzeffh0errsqr)

    end function  HST_LnLike

    end module HST
