
    module ObjectLists
    !Implement lists of arbitrary objects
    !AL Oct 2012
    implicit none

    type Object_pointer
        class(*), pointer :: p => null()
        class(*), pointer :: tag => null()
    end type Object_pointer

    Type TObjectList
        integer :: Count =0
        integer :: Delta = 64
        integer :: Capacity = 0
        logical :: OwnsObjects = .true.
        Type(Object_pointer), dimension(:), allocatable :: Items
    contains
    procedure :: Init
    procedure :: Add
    procedure :: DeleteItem
    procedure :: FreeItem
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
    integer i
    
    if (allocated(L%Items)) then
        if (L%OwnsObjects) then
         do i=1,L%count
          call L%FreeItem(i)
         end do
        end if
        deallocate (L%Items)
    end if
    L%Count = 0
    L%Capacity = 0

    end subroutine Clear

    subroutine Add(L, C, copy)
    Class(TObjectList) :: L
    class(*), intent(in), target :: C
    logical, intent(in), optional :: copy

    if (L%Count == L%Capacity) call L%SetCapacity(L%Capacity + L%Delta)
    L%Count = L%Count + 1
    if (present(copy)) then
     if (L%OwnsObjects .and. copy) then
          allocate(L%Items(L%Count)%P, source=C)
     elseif (copy) then 
         stop 'ObjectLists: Cannot add copy to un-owned list'
     else 
         L%Items(L%Count)%P=>C
     end if
    else
     L%Items(L%Count)%P=>C
    end if

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

    subroutine FreeItem(L, i)
    Class(TObjectList) :: L
    integer, intent(in) :: i

    if (associated(L%Items(i)%P)) deallocate(L%Items(i)%P)
    if (associated(L%Items(i)%tag)) deallocate(L%Items(i)%tag)
    L%Items(i)%P=> null()
    L%Items(i)%tag=> null()

    end subroutine FreeItem

    subroutine DeleteItem(L, i)
    Class(TObjectList) :: L
    integer, intent(in) :: i
    
    if (L%OwnsObjects) call L%FreeItem(i)
    if (L%Count > 1) L%Items(i:L%Count-1) = L%Items(i+1:L%Count)
    L%Count = L%Count -1

    end subroutine DeleteItem

    end module ObjectLists