    module likelihood
    !DataLikelihood is data of a particular kind, CMB, BOA, etc.
    !Can contain set of DataItem, where DataItem stores the separate likelihoods and names
    use AMLUtils
    use settings
    use cmbtypes
    use IniFile
    use ObjectLists, only: TObjectList
    implicit none

    integer, parameter :: max_likelihood_functions = 50

    Type DataItem
        character(Len=80) :: name
        real :: LogLike = 0.
    end Type DataItem
    
    Type, extends(TObjectList) :: DataItemList
    contains
    procedure :: Item => DataItemListItem
    end type DataItemList

    type :: DataLikelihood
        Type(DataItemList) :: datasets = DataItemList()
        logical :: dependent_params(num_params) = .false.
        real :: TotalLogLike = 0.
        logical :: needs_background_functions = .true.
        logical :: needs_linear_pk = .false.
        integer :: needs_cl_lmax = 0
        character(LEN=80) :: LikelihoodName = ''
    contains
    procedure :: Init 
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

    function DataItemListItem(L, i) result(P)
    Class(DataItemList) :: L
    integer, intent(in) :: i
    Class(DataItem), pointer :: P

    select type (Item => L%Items(i)%P)
    Class is (DataItem)
        P => Item 
        class default
        stop 'List contains non-DataItem item'
    end select

    end function DataItemListItem
    
    
    subroutine Init(like, ini)
    use IniFile
    class(DataLikelihood) :: like
    Type(TIniFile) :: ini
    end subroutine 

    function LogLike(like, CMB, Theory)
    class(DataLikelihood) :: like
    Type(CMBParams) :: CMB
    Type(CosmoTheory) :: Theory
    real LogLike
    
    LogLike= logZero
    stop 'likeliood Init should not be called'
    end function



    end module likelihood