    module likelihood
    !DataLikelihood is an instance of data of a particular kind, CMB, BOA, etc.
    !Can be multiple of same type
    use AMLUtils
    use settings
    use IniFile
    use ObjectLists, only: TObjectList
    use ParamNames
    implicit none

    integer, parameter :: max_likelihood_functions = 50

    type :: DataLikelihood
        integer :: speed = 0  !negative for slow likelihoods, larger positive for faster
        character(LEN=80) :: name = ''
        character(LEN=80) :: LikelihoodType= ''
        character(LEN=80) :: version = ''
        Type(TParamNames) :: nuisance_params
        !Internally calculated
        logical :: dependent_params(max_num_params) = .false.
        integer, allocatable :: nuisance_indices(:)
        integer :: new_param_block_start, new_params
    contains
    procedure :: LogLike
    procedure :: LogLikeTheory !same as above when extra info not needed
    procedure :: loadParamNames
    end type DataLikelihood

    !This is the global list of likelihoods we will use
    Type, extends(TObjectList) :: LikelihoodList
    integer :: first_fast_param =0
    contains
    procedure :: Item => LikelihoodItem
    procedure :: WriteLikelihoodContribs
    procedure :: AddNuisanceParameters
    procedure :: Compare => CompareLikes
    end type LikelihoodList

    Type(LikelihoodList), target, save :: DataLikelihoods

    contains

    function LikelihoodItem(L, i) result(P)
    Class(LikelihoodList) :: L
    integer, intent(in) :: i
    Class(DataLikelihood), pointer :: P

    select type (like => L%Items(i)%P)
    class is (DataLikelihood)
        P => like
        class default
        stop 'List contains non-DataLikelihood item'
    end select

    end function LikelihoodItem

    subroutine WriteLikelihoodContribs(L, aunit, likelihoods)
    Class(LikelihoodList) :: L
    integer, intent(in) :: aunit
    real(mcp), intent(in) :: likelihoods(*)
    integer i
    Class(DataLikelihood), pointer :: LikeItem

    do i=1,L%Count
        LikeItem =>  L%Item(i)
        write (aunit,'(2f11.3)',advance='NO') likelihoods(i),likelihoods(i)*2
        write(aunit,'(a)',advance='NO') '   '//trim(LikeItem%LikelihoodType)//': '//trim(LikeItem%name)
        if (LikeItem%Version/='') write(aunit,'(a)',advance='NO') ' '//trim(LikeItem%Version)
        write(aunit,'(a)') ''
    end do

    end subroutine WriteLikelihoodContribs

    integer function CompareLikes(this, R1, R2) result(comp)
    Class(LikelihoodList) :: this
    class(*) R1,R2

    select type (RR1 => R1)
    class is (DataLikelihood)
        select type (RR2 => R2)
        class is (DataLikelihood)
            comp = RR1%speed - RR2%speed
            return
        end select
    end select

    end function CompareLikes



    subroutine AddNuisanceParameters(L, Names)
    use ParamNames
    Class(LikelihoodList) :: L
    Type(TParamNames) :: Names
    Type(DataLikelihood), pointer :: DataLike
    integer i,j

    call DataLikelihoods%Sort
    L%first_fast_param=0
    do i=1,DataLikelihoods%Count
        DataLike=>DataLikelihoods%Item(i)
        if (Feedback>0 .and. MPIrank==0) print *,'adding parameters for: '//trim(DataLIke%name)
        DataLike%new_param_block_start = Names%num_MCMC +1
        if (DataLike%nuisance_params%num_derived>0) call MpiStop('No support for likelihood derived params yet')
        call ParamNames_Add(Names, DataLike%nuisance_params)
        if (Names%num_MCMC > max_num_params) call MpiStop('increase max_data_params in settings.f90')
        DataLike%new_params = Names%num_MCMC - DataLike%new_param_block_start + 1
        allocate(DataLike%nuisance_indices(DataLike%nuisance_params%num_MCMC))
        if (DataLike%nuisance_params%num_MCMC/=0) then
            do j=1, DataLike%nuisance_params%num_MCMC
                DataLike%nuisance_indices(j) = ParamNames_index(Names,DataLike%nuisance_params%name(j))
            end do
            if (any(DataLike%nuisance_indices==-1)) call MpiStop('AddNuisanceParameters: unmatched data param')
            DataLike%dependent_params(DataLike%nuisance_indices) = .true.
            if (Feedback>1 .and. MPIrank==0) print *,trim(DataLike%name)//' data param indices:', DataLike%nuisance_indices
            if (L%first_fast_param==0 .and. DataLike%speed >=0 .and. &
                    DataLike%new_params>0) L%first_fast_param = DataLike%new_param_block_start
        end if
    end do

    end subroutine AddNuisanceParameters

    function logLikeTheory(like, CMB)
    !For likelihoods that don't need Theory or DataParams
    class(DataLikelihood) :: like
    class(*) :: CMB
    real(mcp) LogLikeTheory

    logLikeTheory= logZero
    stop 'logLikeTheory or logLike should not be overridden'
    end function

    function LogLike(like, CMB, Theory, DataParams)
    class(DataLikelihood) :: like
    class(*) :: CMB
    class(*) :: Theory
    real(mcp) :: DataParams(:)
    real(mcp) LogLike

     logLike = like%logLikeTheory(CMB)
    end function

    subroutine loadParamNames(like, fname)
    class(DataLikelihood) :: like
    character(LEN=*), intent(in) :: fname

    call ParamNames_init(like%nuisance_params, fname)

    end subroutine loadParamNames


    end module likelihood
