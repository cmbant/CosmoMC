    ! Hubble measurement constraint, based on corresponding angular diamter distance value and error

    module HST
    use CosmologyTypes
    use Likelihood_Cosmology
    implicit none
    private

    type, extends(TCosmoCalcLikelihood) :: HSTLikelihood
        real(mcp) :: angconversion = 11425.8d0
        !angconversion converts inverse of the angular diameter distance at z = zeff to H0
        !for the fiducial cosmology (omega_k = 0, omega_lambda = 0.7, w = -1)
        !likelihood is in terms of inverse of the angular diameter distance, so includes the tiny cosmological
        !dependence of the measurement (primarily on w) correctly.
        real(mcp) :: zeff = 0.04d0
        real(mcp) :: H0, H0_err
    contains
    procedure :: LogLikeTheory => HST_LnLike
    end type HSTLikelihood


    public HSTLikelihood, HSTLikelihood_Add
    contains

    subroutine HSTLikelihood_Add(LikeList, Ini)
    class(TLikelihoodList) :: LikeList
    class(TSettingIni) :: ini
    Type(HSTLikelihood), pointer :: this

    if (Ini%Read_Logical('use_HST',.false.)) then
        allocate(this)
        this%LikelihoodType = 'Hubble'
        this%name= Ini%Read_String('Hubble_name')
        this%H0 = Ini%Read_Double('Hubble_H0')
        this%H0_err = Ini%Read_Double('Hubble_H0_err')
        call Ini%Read('Hubble_zeff',this%zeff)
        call Ini%Read('Hubble_angconversion',this%angconversion)
        this%needs_background_functions = .true.
        call LikeList%Add(this)
    end if

    end subroutine HSTLikelihood_Add

    real(mcp) function HST_LnLike(this, CMB)
    Class(HSTLikelihood) :: this
    Class(CMBParams) CMB
    real(mcp) :: theoryval

    theoryval = this%angconversion/this%Calculator%AngularDiameterDistance(this%zeff)
    HST_LnLike = (theoryval - this%H0)**2/(2*this%H0_err**2)

    end function  HST_LnLike

    end module HST
