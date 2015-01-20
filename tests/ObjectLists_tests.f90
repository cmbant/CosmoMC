    module ListTests
    use ObjectLists
    implicit none
    character(LEN=0), target :: Empty_String = ''

    type Object_pointer
        class(*), pointer :: p => null()
        class(*), pointer :: Object => null()
    end type Object_pointer

    contains


    function RunListTests() result (fails)
    integer fails
    fails = 0
    call test_TRealArrayList(fails)
    call test_TRealLIst(fails)
    call test_TStringLIst(fails)
    end function RunListTests

    subroutine test_TRealLIst(fails)
    Type(TRealList) T
    integer i, fails

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
        fails= fails + 1
        print *,'error'
        print *, T%AsArray()
    end if

    end subroutine

    subroutine test_TRealArrayList(fails)
    Type(TRealArrayList) T
    integer fails

    print *, 'test_TRealArrayList'
    call T%Add([0.5d0])
    call T%Add([-3.d0, 1.d0])
    if (all(T%Item(1) == [0.5d0]) .and. T%item(2,2) ==1.d0) then
        print *,'OK'
    else
        fails= fails+1
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

    subroutine test_TStringLIst(fails)
    Type(TStringList) :: T
    integer i, fails
    
    print *, 'test_TStringLIst'
    call T%Add('here')
    call T%Add('there')
    call T%Add('alpha')
    call T%sort()

    do i=1,2
        if (T%Item(1) /= 'alpha' .or. T%Item(3) /='there') then
            print *, 'Error'
            fails = fails + 1
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
         fails = fails + 1
    else
        print *,'OK'
    end if

    end subroutine test_TStringLIst

    end module ListTests
