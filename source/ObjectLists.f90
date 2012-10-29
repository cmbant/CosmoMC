
    module ObjectLists
    !Implement lists of arbitrary objects
    !AL Oct 2012
    implicit none

    type Object_pointer
        class(*), pointer :: p
    end type Object_pointer

    Type TObjectList
        integer :: Count =0
        integer :: Delta = 64
        integer :: Capacity = 0
        Type(Object_pointer), dimension(:), allocatable :: Items
    contains
    procedure :: Init
    procedure :: Add
    procedure :: DeleteItem
    procedure :: SetCapacity
    FINAL :: Clear
    end Type TObjectList

    contains

    subroutine Init(L)
    Class(TObjectList) :: L

    L%Count = 0
    L%Capacity = 0
    L%Delta = 32

    end subroutine Init


    subroutine Clear(L)
    Type(TObjectList) :: L

    if (allocated(L%Items)) deallocate (L%Items)
    L%Count = 0
    L%Capacity = 0

    end subroutine Clear


    subroutine Add(L, C)
    Class(TObjectList) :: L
    class(*), intent(in), target :: C

    if (L%Count == L%Capacity) call SetCapacity(L, L%Capacity + L%Delta)
    L%Count = L%Count + 1
    L%Items(L%Count)%P => C

    end subroutine Add

    subroutine SetCapacity(L, C)
    Class(TObjectList) :: L
    integer C
    Type(Object_pointer), dimension(:), allocatable :: TmpItems

    if (L%Count > 0) then
        if (C < L%Count) stop 'SetCapacity: smaller than Count'
        allocate(TmpItems(L%Count))
        TmpItems = L%Items(1:L%Count)
        deallocate(L%Items)
        allocate(L%Items(C))
        L%Items(1:L%Count) = TmpItems
        deallocate(TmpItems)
    else
        allocate(L%Items(C))
    end if  
    L%Capacity = C
    end subroutine SetCapacity

    subroutine DeleteItem(L, i)
    Class(TObjectList) :: L
    integer, intent(in) :: i
    !
    if (L%Count > 1) L%Items(i:L%Count-1) = L%Items(i+1:L%Count)
    L%Count = L%Count -1

    end subroutine DeleteItem

    end module ObjectLists