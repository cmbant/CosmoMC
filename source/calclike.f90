    module CalcLike
    use cmbtypes
    use DataLikelihoodList
    implicit none

    real(mcp) :: Temperature  = 1

    Type LikeCalculator
        type(ParamSet), pointer :: Params
        type (CMBParams), pointer :: CMB
        logical changeMask(max_num_params)
        logical SlowChanged, PowerChanged
    end Type LikeCalculator

    contains

    subroutine AddLike(CurrentLike, LikeToAdd)
    real(mcp), intent(in) :: LikeToAdd
    real(mcp) CurrentLike
    if (CurrentLike/=LogZero) then
        if (LikeToAdd == logZero) then
            CurrentLike = LogZero
        else
            CurrentLike = CurrentLike + LikeToAdd/Temperature
        end if
    end if
    end subroutine AddLike

    function GenericLikelihoodFunction(Params)
    type(ParamSet)  Params
    real(mcp) :: GenericLikelihoodFunction
    real(mcp), allocatable, save :: covInv(:,:)
    real(mcp) X(num_params_used)

    if (test_likelihood) then
        if (.not. allocated(covInv)) then
            allocate(covInv(num_params_used,num_params_used))
            covInv = test_cov_matrix
            call Matrix_Inverse(covInv)
        end if
        X = Params%P(params_used) - scales%Center(params_used)
        GenericLikelihoodFunction= dot_product(X, matmul(covInv, X))/2
    else
        !Used when you want to plug in your own CMB-independent likelihood function:
        !set generic_mcmc=.true. in settings.f90, then write function here returning -Ln(Likelihood)
        !Parameter array is Params%P
        GenericLikelihoodFunction = LogZero
        call MpiStop('GenericLikelihoodFunction: need to write this function!')
    end if

    end function GenericLikelihoodFunction

    function GetLogLikeBounds(Params)
    type(ParamSet):: Params
    real(mcp) :: GetLogLikeBounds

    if (any(Params%P(:num_params) > Scales%PMax(:num_params)) .or. &
    any(Params%P(:num_params) < Scales%PMin(:num_params))) then
        GetLogLikeBounds = logZero
    else
        GetLogLikeBounds=0
    end if

    end function GetLogLikeBounds

    function GetLogLike(Params) !Get -Ln(Likelihood) for chains
    type(ParamSet), target :: Params
    Type (CMBParams), target :: CMB
    real(mcp) GetLogLike
    Type(LikeCalculator) :: Calc

    GetLogLike = GetLogLikeBounds(Params)
    if (GetLogLike==LogZero) return

    if (generic_mcmc .or. test_likelihood) then
        call AddLike(GetLogLike,GenericLikelihoodFunction(Params))
        call AddLike(GetLogLike,getLogPriors(Params%P))
    else
        Calc%Params => Params
        Calc%CMB=>CMB
        call Parameterization%ParamArrayToTheoryParams(Params%P,CMB)
        call AddLike(GetLogLike, Parameterization%NonBaseParameterPriors(CMB))
        if (GetLogLike == logZero) return
        if (.not. Params%Info%validInfo) then
            Calc%changeMask(1:num_params) = .true.
        else
            Calc%changeMask(1:num_params) = Params%Info%lastParamArray(1:num_params)/=Params%P(1:num_params)
        end if

        if (CalculateRequiredTheoryChanges(Calc)) then
            call AddLike(GetLogLike, GetLogLikeWithTheorySet(Calc))
        else
            GetLogLike = logZero
        end if

        if (GetLogLike/=logZero) Params%Info%lastParamArray(1:num_params) = Params%P(1:num_params)
    end if

    if (Feedback>2 .and. GetLogLike/=LogZero) &
    call DataLikelihoods%WriteLikelihoodContribs(stdout, Params%likelihoods)

    end function GetLogLike

    function CheckPriorCuts(Params)
    real(mcp)  CheckPriorCuts
    Type(ParamSet), target :: Params
    Type (CMBParams), target :: CMB

    CheckPriorCuts = GetLogLikeBounds(Params)
    if (CheckPriorCuts==LogZero) return

    call Parameterization%ParamArrayToTheoryParams(Params%P,CMB)
    CheckPriorCuts = Parameterization%NonBaseParameterPriors(CMB)

    end function CheckPriorCuts

    function GetLogLikePost(Params, do_like)
    !for importance sampling where theory may be pre-stored
    real(mcp)  GetLogLikePost
    Type(ParamSet), target :: Params
    Type (CMBParams), target :: CMB
    Type(LikeCalculator) :: Calc
    logical, optional, intent(in) :: do_like(DataLikelihoods%count)

    GetLogLikePost = GetLogLikeBounds(Params)
    if (GetLogLikePost==LogZero) return

    Calc%slowChanged = .false.
    Calc%ChangeMask = .true.

    call Parameterization%ParamArrayToTheoryParams(Params%P,CMB)
    call AddLike(GetLogLikePost,Parameterization%NonBaseParameterPriors(CMB))
    if (GetLogLikePost == logZero) return

    Calc%Params => Params
    Calc%CMB=>CMB
    call AddLike(GetLogLikePost, GetLogLikeWithTheorySet(Calc, do_like))

    end function GetLogLikePost

    function getLogPriors(P) result(logLike)
    integer i
    real(mcp), intent(in) :: P(num_params)
    real(mcp) logLike

    logLike=0
    do i=1,num_params
        if (Scales%PWidth(i)/=0 .and. GaussPriors%std(i)/=0) then
            logLike = logLike + ((P(i)-GaussPriors%mean(i))/GaussPriors%std(i))**2
        end if
    end do
    logLike=logLike/2

    end function getLogPriors

    logical function CalculateRequiredTheoryChanges(Calc)
    Type(LikeCalculator) :: Calc
    integer error

    Calc%SlowChanged = any(Calc%changeMask(1:num_hard))
    Calc%PowerChanged = any(Calc%changeMask(index_initpower:index_initpower+num_initpower-1))
    error=0
    if (Use_CMB .or. Use_LSS) then
        if (Calc%SlowChanged) then
            slow_changes = slow_changes + 1
            call GetNewTransferData(Calc%CMB, Calc%Params%Info,Calc%Params%Theory, error)
        end if
        if ((Calc%SlowChanged .or. Calc%PowerChanged) .and. error==0) then
            if (.not. Calc%SlowChanged) semislow_changes = semislow_changes  + 1
            call GetNewPowerData(Calc%CMB, Calc%Params%Info, Calc%Params%Theory,error)
        end if
    else
        if (Calc%SlowChanged) call GetNewBackgroundData(Calc%CMB, Calc%Params%Theory, error)
    end if
    Calc%Params%Info%validInfo = .true.
    CalculateRequiredTheoryChanges = error==0

    end function CalculateRequiredTheoryChanges


    function GetLogLikeWithTheorySet(Calc, likelihood_mask) result(logLike)
    real(mcp) logLike
    Type(LikeCalculator) :: Calc
    logical, intent(in), optional :: likelihood_mask(DataLikelihoods%count)
    real(mcp) itemLike
    Class(DataLikelihood), pointer :: like
    integer i
    logical backgroundSet
    logical :: do_like(DataLikelihoods%count)

    if (present(likelihood_Mask)) then
        do_like = likelihood_mask
    else
        do_like = .true.
    end if
    backgroundSet = Calc%slowChanged
    logLike = logZero
    do i= 1, DataLikelihoods%count
        if (do_like(i)) then
            like => DataLikelihoods%Item(i)
            if (any(like%dependent_params(1:num_params) .and. Calc%changeMask(1:num_params) )) then
                if (any(like%dependent_params(1:num_hard)) .and. .not. backgroundSet) then
                    call SetTheoryForBackground(Calc%CMB)
                    backgroundSet = .true.
                end if
                itemLike = like%LogLike(Calc%CMB, Calc%Params%Theory, Calc%Params%P(like%nuisance_indices))

                if (itemLike == logZero) return
                Calc%Params%Likelihoods(i) = itemLike
            end if
        end if
    end do
    logLike = sum(Calc%Params%likelihoods(1:DataLikelihoods%Count))
    logLike = logLike + getLogPriors(Calc%Params%P)

    end function GetLogLikeWithTheorySet


    end module CalcLike
