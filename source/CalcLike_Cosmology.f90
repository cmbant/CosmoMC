    module CalcLike_Cosmology
    use CalcLike
    use DataLikelihoodList
    use CosmologyTypes
    use CosmoTheory
    use Calculator_Cosmology
    use ParamPointSet
    implicit none

    private

    type, extends(TTheoryLikeCalculator) :: TCosmoLikeCalculator
        integer :: semislow_changes=0, slow_changes=0
        class(TCosmologyCalculator), pointer :: CosmoCalc
    contains
    procedure :: CalculateRequiredTheoryChanges =>Cosmo_CalculateRequiredTheoryChanges
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
    class(TDataLikelihood), pointer :: like
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
    Type(TTimer) :: Timer

    this%SlowChanged = any(this%changeMask(1:num_hard))
    this%SemiSlowchanged = any(this%changeMask(index_initpower:index_initpower+num_initpower-1))
    error=0
    if (this%timing) call Timer%Start

    select type (Theory=>this%Params%Theory)
    class is (TCosmoTheoryPredictions)
        select type (CMB=>this%TheoryParams)
        class is (CMBParams)
            if (CosmoSettings%Use_CMB .or. CosmoSettings%Use_LSS .or. CosmoSettings%get_sigma8) then
                if (this%SlowChanged .or. this%SemiSlowchanged .and. .not. BaseParams%block_semi_fast) then
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

    if (this%timing) call Timer%WriteTime('Time for theory')

    this%Params%validInfo = .true.
    Cosmo_CalculateRequiredTheoryChanges = error==0

    end function Cosmo_CalculateRequiredTheoryChanges

    subroutine Cosmo_WritePerformanceStats(this, unit)
    class(TCosmoLikeCalculator) :: this
    integer, intent(in) :: unit

    write(unit,*) 'slow changes', this%slow_changes, ' power changes', this%semislow_changes

    end subroutine Cosmo_WritePerformanceStats

    subroutine Cosmo_UpdateTheoryForLikelihoods(this, Params)
    class(TCosmoLikeCalculator) :: this
    class(TCalculationAtParamPoint) :: Params
    integer :: index_error

    select type (Theory=>Params%Theory)
    class is (TCosmoTheoryPredictions)
        if (CosmoSettings%Use_LSS) then
            if(Theory%sigma_8==0) &
                call MpiStop('ERROR: Matter power/sigma_8 have not been computed. Use redo_theory and redo_pk')
            if (allocated(Theory%MPK)) then
                if((CosmoSettings%power_redshifts(CosmoSettings%num_power_redshifts)-Theory%MPK%y(Theory%MPK%ny))>1.d-3)then
                    write(*,*) 'ERROR: Thes elected datasets call for a higher redshift than has been calculated'
                    write(*,*) '       Use redo_theory and redo_pk'
                    call MpiStop()
                end if
                if(CosmoSettings%num_power_redshifts > Theory%MPK%ny)then
                    write(*,*) 'ERROR: The selected datasets call for more redshifts than are calculated'
                    write(*,*) '       Use redo_theory and redo_pk'
                    call MpiStop()
                end if
            end if
            index_error =0
            !AL for the moment only allow reading with matching mpk
            !call IndexExactRedshifts(Theory%MPK%y,index_error)
            !if(index_error>0)then
            !    write(*,*) 'ERROR: One of the datasets needs an exact redshift that is not present '
            !    write(*,*) '       Use redo_theory and redo_pk'
            !    call MpiStop()
            !end if
            if(CosmoSettings%use_nonlinear .and. .not. allocated(Theory%NL_MPK))then
                write(*,*) 'ERROR: One of the datasets wants a non-linear MPK which is not present '
                write(*,*) '       Use redo_theory and redo_pk or turn off non-linear'
                call MpiStop()
            end if
        end if
    end select

    end subroutine Cosmo_UpdateTheoryForLikelihoods


    end module CalcLike_Cosmology
