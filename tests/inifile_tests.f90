    module IniTests
    use IniObjects
    use StringUtils
    implicit none
    character(LEN=0), target :: Empty_String = ''

    contains


    function RunIniTests() result (fails)
    integer fails
    fails = 0

    call test_Read(fails)
    end function RunIniTests

    subroutine test_read(fails)
    integer fails
    Type(TIniFile) :: Ini
    character(LEN=:), allocatable :: S
    character (LEN=200), dimension(1) :: Lines
    
    Lines(1) = 'parameter = test$(PATH)/mypath$(PATH)'
    
    call Ini%Open_FromLines(Lines,1)
    S = Ini%Read_String('parameter')
    if (S /= 'test'//GetEnvironmentVariable('PATH')//'/mypath'//GetEnvironmentVariable('PATH')) then
        fails = fails + 1
        print *, 'error'
    else
        print *, 'OK Ini path'
    end if
        
    end subroutine test_read
    
    end module IniTests
