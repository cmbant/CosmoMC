    module classes
    implicit none
    Type AData
    end Type

    Type A
        integer :: aa =1
    contains
    procedure :: Work
    end type

    Type, extends(AData) :: BData
    end Type

    Type, extends(A) :: B
    contains
    procedure :: Work => Work2
    end type

    contains

    subroutine Work(this, D)
    class(A) :: this
    class(*), target :: D

    associate (F => A())
        print *,F.aa
    end associate
    end subroutine

    subroutine Work2(this, D)
    class(B) :: this
    class(*), target :: D
    class(BData), pointer :: B

    B=>D	

    end subroutine


    end module classes


    module tests
    use ObjectLists
    implicit none

    type Object_pointer
        class(*), pointer :: p => null()
        class(*), pointer :: Object => null()
    end type Object_pointer

    contains

    subroutine RunTests

    call test_TRealArrayList
    call test_TRealLIst
    call test_TStringLIst
    end subroutine RunTests

    subroutine test_TRealLIst
    Type(TRealList) T
    integer i

    print *, 'test_TRealLIst'
    call T%Add(0.5d0)
    call T%Add(-3.d0)
    call T%Add(12.d0)
    do i=1, 40
        call T%Add(i*1.d0)
    end do
    call T%Add(-5.d0)

    call T%Sort()
    if (T%Item(1) == -5.d0 .and. T%item(3) ==0.5d0) then
        print *,'OK'
    else
        print *,'error'
        print *, T%AsArray()
    end if

    end subroutine

    subroutine test_TRealArrayList
    Type(TRealArrayList) T

    print *, 'test_TRealArrayList'
    call T%Add([0.5d0])
    call T%Add([-3.d0, 1.d0])
    if (all(T%Item(1) == [0.5d0]) .and. T%item(2,2) ==1.d0) then
        print *,'OK'
    else
        print *,'error'
    end if

    end subroutine test_TRealArrayList


    subroutine ReadWrite(T)
    class(TObjectList) T
    integer unit

    open (newunit=unit,file='test.bin',status='replace', form='unformatted')
    call T%SaveBinary(unit)
    close(unit)
    open (newunit=unit,file='test.bin',status='old', form='unformatted')
    call T%Clear()
    call T%ReadBinary(unit)
    close(unit, status = 'DELETE')

    end subroutine

    subroutine test_TStringLIst
    Type(TStringList) :: T
    integer i 


    print *, 'test_TStringLIst'
    call T%Add('here')
    call T%Add('there')
    call T%Add('alpha')
    call T%sort()

    do i=1,2
        if (T%Item(1) /= 'alpha' .or. T%Item(3) /='there') then
            print *, 'Error'
            call T%WriteItems()
        else
            print *, 'OK'
        end if
        if (i==2) exit
        call ReadWrite(T)
    end do
    call T%Clear()
    call T%Add('here','there')
    call T%Add('hot','alpha')
    if (T%ValueOf('here')/='there') then
        print *,'value Of error'
    else
        print *,'OK'
    end if

    end subroutine test_TStringLIst


    subroutine doadd(P, C)
    Type(Object_pointer) P
    class(*), target :: C

    P%P=> C

    end subroutine

    subroutine ifortbug !ifort 14, inconistent freeing - reported cf. 
    ! http://software.intel.com/en-us/forums/topic/390944
    character(LEN=:), pointer :: St
    Type(Object_pointer) P
    Type(Object_pointer), pointer :: PP
    class(*), pointer :: P2

    allocate(Object_pointer::P2)
    call doadd(P,P2)
    deallocate(P%P) !No problem

    allocate(PP)
    call doadd(P,PP)
    deallocate(P%P) !forrtl: severe (173): A pointer passed to DEALLOCATE points to an object that cannot be deallocated

    allocate(character(5)::St)
    call doadd(P,St)
    deallocate(P%P) !forrtl: severe (173): A pointer passed to DEALLOCATE points to an object that cannot be deallocated

    end subroutine ifortbug

    end module tests

    program tester
    use tests
    call RunTests
    end program
