    module likelihood
    !DataLikelihood is an instance of data of a particular kind, CMB, BOA, etc.
    !Can be multiple of same type
    use AMLUtils
    use settings
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
        character(LEN=80) :: version = ''
    contains
    procedure :: LogLike
    end type DataLikelihood

 !This is the global list of likelihoods we will use
    Type, extends(TObjectList) :: LikelihoodList 
    contains
    procedure :: Item => LikelihoodItem
    procedure :: WriteLikelihoodContribs
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
        write(aunit,'(a)') '   '//trim(LikeItem%LikelihoodType)//': '//trim(LikeItem%name)
     end do

    end subroutine WriteLikelihoodContribs
    
    function LogLike(like, CMB, Theory)
    class(DataLikelihood) :: like
    class(*) :: CMB
    class(*) :: Theory
    real(mcp) LogLike
    
    LogLike= logZero
    stop 'likeliood Init should not be called'
    end function



    end module likelihood