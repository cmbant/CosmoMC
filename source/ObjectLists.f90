
    module ObjectLists
    !Implement lists of arbitrary objects
    !AL Oct 2012
    implicit none

    !type ReadWriteObject
    !!Ancestor for list items that can be read and written
    !contains
    ! procedure, abstract :: WriteBinary
    ! procedure, abstract :: ReadBinary
    !end type ReadWriteObject

    type Object_pointer
        class(*), pointer :: p => null()
    end type Object_pointer

    type Object_array_pointer
        class(*), pointer :: p(:) => null()
    end type Object_array_pointer

#ifdef SINGLE
    integer, parameter :: list_prec = Kind(1.0)
#else
    integer, parameter :: list_prec = Kind(1.d0)
#endif

    Type TObjectList
        integer :: Count =0
        integer :: Delta = 32
        integer :: DeltaScale = 10 
        !expanding expanding Items, expand array by Delta + Count/DeltaScale
        integer :: Capacity = 0
        logical :: OwnsObjects = .true.
        Type(Object_pointer), allocatable :: Items(:)
    contains
    procedure :: AddArray
    procedure :: AddItem
    procedure :: AddCopy
    procedure :: AssignPointers
    procedure :: DeleteItem
    procedure :: FreeItem
    procedure :: SetCapacity
    procedure :: ArrayItem
    procedure :: ArrayItemIndex
    procedure :: SaveBinary
    procedure :: ReadBinary
    procedure :: QuickSortArr
    procedure :: Thin
    procedure :: SortArr
    procedure :: Swap
    procedure :: Compare => CompareReal
    procedure :: Clear
    procedure :: DeltaSize
    FINAL :: finalize
    generic :: Add => AddItem, AddArray
    end Type TObjectList

    Type, extends(TObjectList):: TRealArrayList
    contains
    procedure :: Value
    procedure :: RealItem
    generic :: Item => Value, RealItem
    end Type TRealArrayList


    contains

    subroutine Clear(L)
    Class(TObjectList) :: L
    integer i

    if (allocated(L%Items)) then
        do i=1,L%count
            call L%FreeItem(i)
        end do
        deallocate (L%Items)
    end if
    L%Count = 0
    L%Capacity = 0

    end subroutine Clear

    subroutine finalize(L)
    Type(TObjectList) :: L
    call L%Clear()
    end subroutine finalize


    subroutine AddItem(L, C)
    Class(TObjectList) :: L
    class(*), intent(in), target :: C

    if (L%Count == L%Capacity) call L%SetCapacity(L%Capacity + L%DeltaSize())
    L%Count = L%Count + 1
    L%Items(L%Count)%P=>C

    end subroutine AddItem

    subroutine AddCopy(L, C)
    Class(TObjectList) :: L
    class(*), intent(in) :: C
    class(*), pointer :: P

    if (L%OwnsObjects) then
        allocate(P, source=C)
    else
        stop 'ObjectLists: Cannot add copy to un-owned list'
    end if
    call L%AddItem(P)

    end subroutine AddCopy

    subroutine AssignPointers(L, L2, ixmin, ixmax)
    Class(TObjectList) :: L, L2
    integer, intent(in), optional :: ixmin, ixmax
    integer i1,i2

    call L%Clear()
    i1=1
    i2=L2%Count
    if (present(ixmin)) i1=ixmin
    if (present(ixmax)) i2=ixmax
    call L%SetCapacity(i2-i1+1)
    L%Items = L2%Items(i1:i2)
    L%Count = i2-i1+1
    L%OwnsObjects = .false.

    end subroutine AssignPointers

    integer function DeltaSize(L) 
    Class(TObjectList) :: L

    DeltaSize= L%Delta + L%Count/L%DeltaScale

    end  function

    subroutine SetCapacity(L, C)
    Class(TObjectList) :: L
    integer C
    Type(Object_pointer), dimension(:), allocatable :: TmpItems

    if (L%Count > 0) then
        if (C < L%Count) stop 'SetCapacity: smaller than Count'
        allocate(TmpItems(L%Count), source=L%Items(1:L%Count))
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
    logical want_Dealloc
    if (associated(L%Items(i)%P)) then
        want_Dealloc =  L%OwnsObjects
        select type (point => L%Items(i)%P)
        class is (object_array_pointer)
            if (L%OwnsObjects .and. associated(Point%P)) deallocate(Point%P)
            want_Dealloc = .true.
        end select
        if (want_Dealloc) deallocate(L%Items(i)%P)
    end if

    L%Items(i)%P=> null()
    !    if (associated(L%Items(i)%tag)) deallocate(L%Items(i)%tag)
    !    L%Items(i)%tag=> null()

    end subroutine FreeItem

    subroutine DeleteItem(L, i)
    Class(TObjectList) :: L
    integer, intent(in) :: i

    call L%FreeItem(i)
    if (L%Count > 1) L%Items(i:L%Count-1) = L%Items(i+1:L%Count)
    L%Items(L%Count)%P => null()
    L%Count = L%Count -1

    end subroutine DeleteItem

    subroutine AddArray(L, P)
    Class (TObjectList) :: L
    class(*), target, intent(in) :: P(:)
    !    Type(object_array_pointer), pointer :: Pt
    class(*), pointer :: Pt

    allocate(object_array_pointer::Pt)
    call L%AddItem(Pt)
    select type (Pt)
    class is (object_array_pointer)
        if (L%ownsObjects) then
            allocate(Pt%P(1:SIZE(P)), source= P)
        else
            Pt%P => P
        end if
    end select
    end subroutine AddArray

    function ArrayItem(L, i) result(P)
    Class(TObjectList) :: L
    integer, intent(in) :: i
    Class(*), pointer :: P(:)

    select type (Point=> L%Items(i)%P)
    class is (object_array_pointer)
        P => Point%P
        class default
        stop 'TObjectList: item is not array item'
    end select

    end function ArrayItem

    function ArrayItemIndex(L, i, j) result(P)
    Class(TObjectList) :: L
    integer, intent(in) :: i, j
    Class(*), pointer :: P
    Class(*), pointer :: Arr(:)

    Arr => L%ArrayItem(i)
    P => Arr(j)

    end function ArrayItemIndex


    subroutine SaveBinary(L,fid)
    Class(TObjectList) :: L
    integer, intent(in) :: fid
    integer i,k
    class(*), pointer :: P(:)

    write (fid) L%Count
    do i=1,L%Count
        select type (Item=> L%Items(i)%P)
        class is (object_array_pointer)
            P => L%ArrayItem(i)
            select type (Point=> P)
            Type is (real)
                k=1
                write(fid) size(P),k 
                write(fid) Point
            Type is (double precision)
                k=2
                write(fid) size(P),k 
                write(fid) Point
            Type is (integer)
                k=3
                write(fid) size(P),k 
                write(fid) Point
            Type is (logical)
                k=4
                write(fid) size(P),k 
                write(fid) Point
                class default
                stop 'TObjectList: Unknown type to save'
            end select
            class default
            stop 'TObjectList: not implemented non-array save'
        end select
    end do

    end subroutine SaveBinary

    subroutine ReadBinary(L,fid)
    Class(TObjectList) :: L
    integer, intent(in) :: fid
    integer num,i,sz, k
    real, pointer :: ArrR(:)
    double precision, pointer :: ArrD(:)
    integer, pointer :: ArrI(:)
    logical, pointer :: ArrL(:)

    call L%Clear()
    L%OwnsObjects = .false.
    read (fid) num
    call L%SetCapacity(num)
    do i=1,num
        read(fid) sz, k
        if (k==1) then
            allocate(ArrR(sz))
            read(fid) ArrR
            call L%AddArray(ArrR)
        else if (k==2) then
            allocate(ArrD(sz))
            read(fid) ArrD
            call L%AddArray(ArrD)
        else if (k==3) then
            allocate(ArrI(sz))
            read(fid) ArrI
            call L%AddArray(ArrI)
        else if (k==4) then
            allocate(ArrL(sz))
            read(fid) ArrL
            call L%AddArray(ArrL)
        end if
    end do
    L%Count = num
    L%OwnsObjects = .true.

    end subroutine ReadBinary

    subroutine Thin(L, i)
    Class(TObjectList):: L
    integer, intent(in) :: i
    integer newCount, j
    Type(Object_pointer), dimension(:), pointer :: TmpItems

    if (L%Count > 1) then
        newCount = (L%Count-1)/i+1
        allocate(TmpItems(newCount))
        TmpItems= L%Items(1:L%Count:i)
        if (L%OwnsObjects) then
            do j=1,L%count
                if (mod(j-1,i)/=0) call L%FreeItem(j)
            end do
        end if
        deallocate(L%Items)
        L%Capacity = newCount
        allocate(L%Items(L%Capacity), source = TmpItems)
        L%Count = newCount
        deallocate(TmpItems)
    end if

    end subroutine Thin

    subroutine Swap(L, i, j)
    Class(TObjectList) :: L
    integer, intent(in) :: i, j
    type(Object_pointer) :: temp

    temp = L%Items(i)
    L%Items(i) = L%Items(j)
    L%Items(j) = temp

    end subroutine Swap

    recursive subroutine QuickSortArr(this, Lin, R, index)
    !Sorts an array of pointers by the value of the index'th entry
    Class(TObjectList) :: this
    integer, intent(in) :: Lin, R, index
    integer I, J, L
    class(*), pointer :: P

    L = Lin
    do

    I = L
    J = R
    P => this%ArrayItemIndex((L + R)/2, index)

    do
        do while (this%Compare(this%ArrayItemIndex(I, Index),P) <  0)
            I = I + 1
        end do

        do while (this%Compare(this%ArrayItemIndex(J,Index), P) > 0)
            J = J - 1
        end do

        if (I <= J) then
            call this%Swap(I,J)
            I = I + 1
            J = J - 1
        end if
        if (I > J) exit

    end do
    if (L < J) call this%QuickSortArr(L, J, index)
    L = I
    if (I >= R) exit
    end do

    end subroutine QuickSortArr

    subroutine SortArr(L, index)
    Class(TObjectList) :: L
    integer, intent(in) :: index

    call L%QuickSortArr(1, L%Count, index)

    end subroutine SortArr

    ! List of arrays of reals

    function RealItem(L, i) result(P)
    Class(TRealArrayList) :: L
    integer, intent(in) :: i
    real(list_prec), pointer :: P(:)
    class(*), pointer :: Item(:)

    Item => L%ArrayItem(i)
    select type (pt=>Item)
    type is (real(kind=list_prec))
        P=> pt
        class default
        stop 'TRealArrayList: object of wrong type'
    end select

    end function RealItem

    function Value(L, i, j) result(P)
    Class(TRealArrayList) :: L
    integer, intent(in) :: i, j
    real(list_prec) :: P
    class(*), pointer :: C

    C => L%ArrayItemIndex(i,j)
    select type (Arr=> C)
    Type is (real(list_prec))
        P = Arr
    end select

    end function Value

    integer function CompareReal(this, R1, R2) result(comp)
    Class(TObjectList) :: this
    class(*) R1,R2 
    real(list_prec) R

    select type (RR1 => R1)
    type is (real(list_prec))
        select type (RR2 => R2)
        type is (real(list_prec))
            R = RR1-RR2
            if (R< 0) then
                comp =-1
            elseif (R>0) then
                comp = 1
            else 
                comp = 0
            end if
            return
        end select
        class default
        stop 'TObjectList: Compare not defined for this type' 
    end select

    end function CompareReal

    end module ObjectLists

