    ! Element abundance measurements compared to BBN predictions for given omegabh2 and Neff

    module ElementAbundances
    use CosmologyTypes
    use Likelihood_Cosmology
    use Calculator_Cosmology
    use bbn
    implicit none
    private

    type, extends(TCosmologyLikelihood) :: AbundanceLikelihood
        character(LEN=:), allocatable :: measurement
        class(TInterpGrid2D), pointer :: BBN_value, BBN_error
        real(mcp) mean, error
    contains
    procedure :: ReadIni => Abundance_ReadIni
    procedure :: LogLikeTheory => Abundance_LnLike
    end type AbundanceLikelihood


    public AbundanceLikelihood, AbundanceLikelihood_Add
    contains

    subroutine AbundanceLikelihood_Add(LikeList, Ini)
    class(TLikelihoodList) :: LikeList
    class(TSettingIni) :: ini
    Type(AbundanceLikelihood), pointer :: like
    Type(TSettingIni) :: DataSets
    integer i

    call Ini%TagValuesForName('abundance_dataset', DataSets, filename=.true.)

    do i= 1, DataSets%Count
        allocate(AbundanceLikelihood::like)
        call like%ReadDatasetFile(Datasets%Value(i))
        call like%ReadParams(Ini)
        like%speed = 10
        like%needs_background_functions = .true.
        like%needs_powerspectra = .false.
        like%tag = Datasets%Name(i)
        like%LikelihoodType = 'Abund'
        call LikeList%Add(like)
    end do
    if (Feedback > 1 .and. DataSets%Count>0 ) write (*,*) 'read abundance data sets'

    end subroutine AbundanceLikelihood_Add

    subroutine Abundance_ReadIni(this, Ini)
    class(AbundanceLikelihood) :: this
    class(TSettingIni) :: Ini

    this%measurement = Ini%Read_String('measurement', NotFoundFail=.true.)
    if (this%measurement=='Yp') then
        !Use nucleon fraction here, not mass fraction
        this%BBN_value => BBN_YpBBN
        this%BBN_error => BBN_YpBBN_err
    else if (this%measurement=='D/H') then
        this%BBN_value => BBN_DH
        this%BBN_error => BBN_DH_err
    else
        call MpiStop('Un-recognised measurement name: '//this%measurement)
    end if
    this%mean = Ini%Read_Double('mean')
    this%error = Ini%Read_Double('error')
    call this%TCosmologyLikelihood%ReadIni(Ini)

    end subroutine Abundance_ReadIni

    real(mcp) function Abundance_LnLike(this, CMB)
    Class(AbundanceLikelihood) :: this
    Class(CMBParams) CMB
    real(mcp) :: theoryval, theoryerr

    theoryval = this%BBN_value%Value(CMB%ombh2,CMB%nnu - standard_neutrino_neff)
    theoryerr = this%BBN_error%Value(CMB%ombh2,CMB%nnu - standard_neutrino_neff)

    Abundance_LnLike = (theoryval - this%mean)**2/(2*(this%error**2+theoryerr**2))

    end function Abundance_LnLike

    end module ElementAbundances

