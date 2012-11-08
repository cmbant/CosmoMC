    module likelihood
    !DataLikelihood is an instance of data of a particular kind, CMB, BOA, etc.
    !Can be multiple of same type
    use AMLUtils
    use settings
    use cmbtypes
    use IniFile
    use ObjectLists, only: TObjectList
    implicit none

    integer, parameter :: max_likelihood_functions = 50

    type :: DataLikelihood
        logical :: dependent_params(num_params) = .false.
        logical :: needs_background_functions = .true.
!not implemented yet..
!        logical :: needs_linear_pk = .false.
!        integer :: needs_cl_lmax = 0
        character(LEN=80) :: name = ''
        character(LEN=80) :: LikelihoodType= ''
    contains
    procedure :: LogLike
    end type DataLikelihood

 !This is the global list of likelihoods we will use
    Type, extends(TObjectList) :: LikelihoodList 
    contains
    procedure :: Item => LikelihoodItem
    end type LikelihoodList

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

    function LogLike(like, CMB, Theory)
    class(DataLikelihood) :: like
    Type(CMBParams) :: CMB
    Type(CosmoTheory) :: Theory
    real LogLike
    
    LogLike= logZero
    stop 'likeliood Init should not be called'
    end function



    end module likelihood