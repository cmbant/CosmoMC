    module CalcLike_Cosmology
    use CalcLike
    use DataLikelihoodList
    use cmbtypes
    use Calculator_Cosmology
    use ParamPointSet
    implicit none

    private

    type, extends(TTheoryLikeCalculator) :: TCosmoLikeCalculator
        integer :: semislow_changes=0, slow_changes=0
        class(TCosmologyCalculator), pointer :: CosmoCalc
    contains
    procedure :: CalulateRequiredTheoryChanges =>Cosmo_CalculateRequiredTheoryChanges
    procedure :: GetTheoryForLike=>Cosmo_GetTheoryForLike
    procedure :: WritePerformanceStats => Cosmo_WritePerformanceStats
    procedure :: InitConfig => TCosmoLikeCalculator_InitConfig
    procedure :: UpdateTheoryForLikelihoods => Cosmo_UpdateTheoryForLikelihoods
    end Type TCosmoLikeCalculator

    public TCosmoLikeCalculator
    contains

    subroutine TCosmoLikeCalculator_InitConfig(this, Config)
    class(TCosmoLikeCalculator) :: this
    class(TGeneralConfig), target :: config

    this%Config => config
    select type (p=>config%Calculator)
    class is (TCosmologyCalculator)
        this%CosmoCalc => p
    end select
    end subroutine TCosmoLikeCalculator_InitConfig

    subroutine Cosmo_GetTheoryForLike(this,Like)
    class(TCosmoLikeCalculator) :: this
    class(DataLikelihood), pointer :: like
    logical, save :: backgroundSet 

    if (.not. associated(like)) then
        !initialize
        backgroundSet = this%slowChanged
        return
    end if

    if (any(like%dependent_params(1:num_hard)) .and. .not. backgroundSet) then
        select type (CMB=>this%TheoryParams)
        class is (CMBParams)
            call this%CosmoCalc%SetParamsForBackground(CMB)
        end select
        backgroundSet = .true.
    end if

    end subroutine Cosmo_GetTheoryForLike


    logical function Cosmo_CalculateRequiredTheoryChanges(this)
    class(TCosmoLikeCalculator) :: this
    integer error

    this%SlowChanged = any(this%changeMask(1:num_hard))
    this%SemiSlowchanged = any(this%changeMask(index_initpower:index_initpower+num_initpower-1))
    error=0

    select type (Theory=>this%Params%Theory)
    class is (TCosmoTheoryPredictions)
        select type (CMB=>this%TheoryParams)
        class is (CMBParams)
            if (Use_CMB .or. Use_LSS .or. get_sigma8) then
                if (this%SlowChanged) then
                    this%slow_changes = this%slow_changes + 1
                    this%Params%validInfo = .false.
                    call this%CosmoCalc%GetNewTransferData(CMB, this%Params%Info,Theory, error)
                end if
                if ((this%SlowChanged .or. this%SemiSlowchanged) .and. error==0) then
                    if (.not. this%SlowChanged) this%semislow_changes = this%semislow_changes  + 1
                    call this%CosmoCalc%GetNewPowerData(CMB, this%Params%Info, Theory,error)
                end if
            else
                if (this%SlowChanged) call this%CosmoCalc%GetNewBackgroundData(CMB, Theory, error)
            end if
        end select
    end select

    this%Params%validInfo = .true.
    Cosmo_CalculateRequiredTheoryChanges = error==0

    end function Cosmo_CalculateRequiredTheoryChanges

    subroutine Cosmo_WritePerformanceStats(this, unit)
    class(TCosmoLikeCalculator) :: this
    integer, intent(in) :: unit

    write(unit,*) 'slow changes', this%slow_changes, 'power changes', this%semislow_changes

    end subroutine Cosmo_WritePerformanceStats

    subroutine Cosmo_UpdateTheoryForLikelihoods(this, Params)
    class(TCosmoLikeCalculator) :: this
    class(TCalculationAtParamPoint) :: Params
    integer :: numz, index_error

    select type (Theory=>Params%Theory)
    class is (TCosmoTheoryPredictions)
        if (Use_LSS) then
            if(Theory%sigma_8==0) &
            call MpiStop('ERROR: Matter power/sigma_8 have not been computed. Use redo_theory and redo_pk')

            numz = size(Theory%redshifts)
            if((power_redshifts(num_power_redshifts)-Theory%redshifts(numz))>1.d-3)then
                write(*,*) 'ERROR: Thes elected datasets call for a higher redshift than has been calculated'
                write(*,*) '       Use redo_theory and redo_pk'
                call MpiStop()
            end if
            if(num_power_redshifts > numz)then
                write(*,*) 'ERROR: The selected datasets call for more redshifts than are calculated'
                write(*,*) '       Use redo_theory and redo_pk'
                call MpiStop()
            end if
            index_error =0
            call IndexExactRedshifts(Theory%redshifts,index_error)
            if(index_error>0)then
                write(*,*) 'ERROR: One of the datasets needs an exact redshift that is not present '
                write(*,*) '       Use redo_theory and redo_pk'
                call MpiStop()
            end if
        end if
    end select

    end subroutine Cosmo_UpdateTheoryForLikelihoods


    end module CalcLike_Cosmology

