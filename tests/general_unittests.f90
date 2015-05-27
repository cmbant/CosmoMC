
    program tester
    use ListTests
    use InterpolationTests
    use IniTests
    implicit none
    
    integer :: fails =0

    fails = fails + RunListTests()
    fails = fails + RunInterpolationTests()
    fails = fails + RunIniTests()

    print *, 'Total fails: ', fails
    
    contains
    
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


    end program
