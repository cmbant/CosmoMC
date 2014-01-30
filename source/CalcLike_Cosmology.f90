    module CalcLike_Cosmology
    use CalcLike
    use DataLikelihoodList
    use cmbtypes
    use Calculator_Cosmology
    implicit none

    type, extends(TheoryLikeCalculator) :: CosmoLikeCalculator
        class(BaseCosmologyCalculator), pointer :: CosmoCalc
    contains
    procedure :: CalulateRequiredTheoryChanges =>Cosmo_CalculateRequiredTheoryChanges
    procedure :: GetTheoryForLike=>Cosmo_GetTheoryForLike
    end Type CosmoLikeCalculator

    contains

    subroutine Cosmo_GetTheoryForLike(this,Like)
    class(CosmoLikeCalculator) :: this
    class(DataLikelihood), pointer :: like
    logical, save :: backgroundSet 

    if (.not. associated(like)) then
        !initialize
        backgroundSet = this%slowChanged
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
    class(CosmoLikeCalculator) :: this
    integer error

    this%SlowChanged = any(this%changeMask(1:num_hard))
    this%SemiSlowchanged = any(this%changeMask(index_initpower:index_initpower+num_initpower-1))
    error=0

    select type (Theory=>this%Params%Theory)
    class is (CosmoTheoryPredictions)
        select type (CMB=>this%TheoryParams)
        class is (CMBParams)
            if (Use_CMB .or. Use_LSS .or. get_sigma8) then
                if (this%SlowChanged) then
                    slow_changes = slow_changes + 1
                    call this%CosmoCalc%GetNewTransferData(CMB, this%Params%Info,Theory, error)
                end if
                if ((this%SlowChanged .or. this%SemiSlowchanged) .and. error==0) then
                    if (.not. this%SlowChanged) semislow_changes = semislow_changes  + 1
                    call this%CosmoCalc%GetNewPowerData(CMB, this%Params%Info, Theory,error)
                end if
            else
                if (this%SlowChanged) call this%CosmoCalc%GetNewBackgroundData(CMB, Theory, error)
            end if
        end select
    end select

    this%Params%Info%validInfo = .true.
    Cosmo_CalculateRequiredTheoryChanges = error==0

    end function Cosmo_CalculateRequiredTheoryChanges

    end module CalcLike_Cosmology

