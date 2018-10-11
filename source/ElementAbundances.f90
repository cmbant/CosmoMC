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
        class(TBBNPredictions), pointer :: BBN_value, BBN_error
        logical :: non_bbn_yhe = .false.
        real(mcp) mean, error
        real(mcp) :: theory_bias_offset = 0._mcp
        real(mcp) :: theory_effective_error = 0._mcp
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
    Type(TSettingIni) :: DataSets, OverrideSettings
    integer i

    call Ini%TagValuesForName('abundance_dataset', DataSets, filename=.true.)

    do i= 1, DataSets%Count
        call Ini%SettingValuesForTagName('abundance_dataset',DataSets%Name(i),OverrideSettings)
        allocate(AbundanceLikelihood::like)
        call like%ReadDatasetFile(Datasets%Value(i), OverrideSettings)
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
    class(TBBNPredictions), pointer :: tmp
    character(LEN=:), allocatable :: theory_table

    this%measurement = Ini%Read_String('measurement', NotFoundFail=.true.)
    if (this%measurement=='Yp') then
        !Use nucleon fraction here, not mass fraction
        if (CosmoSettings%bbn_consistency) then
            this%BBN_value => BBN_YpBBN
            this%BBN_error => BBN_YpBBN_err
        else
            this%non_bbn_yhe = .true.
        end if
    else if (this%measurement=='D/H') then
        if (.not. CosmoSettings%bbn_consistency) &
              call MpiStop('D/H abundance measurement requires BBN consistency')
        this%BBN_value => BBN_DH
        this%BBN_error => BBN_DH_err
    else
        call MpiStop('Un-recognised measurement name: '//this%measurement)
    end if
    this%mean = Ini%Read_Double('mean')
    this%error = Ini%Read_Double('error')
    if (.not. this%non_bbn_yhe) then
        this%theory_bias_offset = Ini%Read_Double('theory_bias_offset',0.d0)
        this%theory_effective_error = Ini%Read_Double('theory_effective_error',0.d0)
        theory_table = Ini%Read_String('theory_table')
        if (theory_table/='' .and. theory_table /=this%BBN_value%BBN_data_file) then
            tmp => this%BBN_value
            allocate(this%BBN_value, source = tmp)
            this%BBN_value%BBN_data_file = theory_table
            if (this%theory_effective_error == 0.d0) then
                tmp => this%BBN_error
                allocate(this%BBN_error, source = tmp)
                this%BBN_value%BBN_data_file = theory_table
            end if
        end if
    end if
    call this%TCosmologyLikelihood%ReadIni(Ini)

    end subroutine Abundance_ReadIni

    real(mcp) function Abundance_LnLike(this, CMB)
    Class(AbundanceLikelihood) :: this
    Class(CMBParams) CMB
    real(mcp) :: theoryval, theoryerr

    if (this%non_bbn_yhe) then
        Abundance_LnLike = (GetYPBBN(CMB%Yhe) - this%mean)**2/this%error**2/2
    else
        theoryval = this%BBN_value%Value(CMB%ombh2,CMB%nnu - standard_neutrino_neff)

        if (this%theory_effective_error>0) then
            theoryerr = this%theory_effective_error
        else
            theoryerr = this%BBN_error%Value(CMB%ombh2,CMB%nnu - standard_neutrino_neff)
        end if

        Abundance_LnLike = (theoryval  + this%theory_bias_offset - this%mean)**2 &
            /(2*(this%error**2+theoryerr**2))

    end if

    end function Abundance_LnLike

    end module ElementAbundances

