    module InterpolationTests
    use Interpolation
    use FileUtils
    implicit None
    integer, parameter, private :: sp= KIND(1.d0)

    contains

    function func(x)
    real(sp) func,x

    func = (x/47.2_sp)**3  + (x/5.5_sp)**2 +x*3.5_sp+4.4579_sp
    end function

    function funcderiv(x)
    real(sp) funcderiv,x

    funcderiv = 3*(x/47.2_sp)**2/47.2_sp  + 2*(x/5.5_sp)/5.5_sp +3.5_sp
    end function

    function RunInterpolationTests() result(fails)
    real(sp), allocatable :: x(:),f(:)
    integer i
    integer fails
    Type(TCubicSpline) :: Irreg, IrregLoad
    Type(TRegularCubicSpline) :: Reg, RegLoad
    Type(TLogRegularCubicSpline) :: Reglog
    real(sp) xx
    real(sP), parameter:: testval = 13.2623_sp
    real(sP), parameter:: testarrayval(3) = [7.3_sp, 9._sp,34.34643_sp]
    real(sp) :: outarray(3), funcarray(3)
    Type(TBinaryFile) :: FB
    
    fails = 0
    allocate(x(100),f(100))
    do i=1, 100
        x(i) = 0.5367_sp*real(i,sp) + 0.3_sp
        f(i) = func(x(i))
    end do
    call Irreg%Init(x,f)
    call Reg%Init(x(1),x(size(x)),100,values=f)
    call RegLog%init(x(1),x(size(x)),100)
    do i=1, 100
        xx= exp(log(x(1)) + RegLog%delta_x*(i-1))
        RegLog%F(i) = func(xx)
    end do

    if (all(abs([Reg%Value(testval), Irreg%Value(testval), RegLog%Value(testval)]-func(testval))<1e-5)) then
        print *,'Value OK'
    else
        fails = fails+1
        print *, 'error'
        print *, Reg%Value(testval), Irreg%Value(testval),RegLog%Value(testval), func(testval)
    end if

    if (all(abs([Reg%Value(x(1)), Irreg%Value(x(1)), RegLog%Value(x(1))]-func(x(1)))<1e-5)) then
        print *,'bottom end value OK'
    else
        fails = fails+1
        print *, 'end error'
        print *, Reg%Value(testval), Irreg%Value(testval),RegLog%Value(testval), func(testval)
    end if

    if (all(abs([Reg%Value(x(100)), Irreg%Value(x(100)), RegLog%Value(x(100))]-func(x(100)))<1e-5)) then
        print *,'top end value OK'
    else
        fails = fails+1
        print *, 'end error'
        print *, Reg%Value(testval), Irreg%Value(testval),RegLog%Value(testval), func(testval)
    end if

    if (all(abs([Reg%Derivative(testval), Irreg%Derivative(testval), RegLog%Derivative(testval)]-funcderiv(testval))<1e-5)) then
        print *,'Derivative OK'
    else
        fails = fails+1
        print *, 'derivative error'
        print *, Reg%Derivative(testval), Irreg%Derivative(testval),RegLog%Derivative(testval), funcderiv(testval)
    end if
    
    do i=1,3
        funcarray(i) = func(testarrayval(i))
    end do
    call Reg%Array(testarrayval, outarray)
    if (all(abs(outarray-funcarray)<1e-7)) then
        print *, 'array value OK'
    else
        fails = fails+1
        print *, 'array error'
        print *, outarray
        print *, funcarray
    end if

    call FB%CreateFile('test.bin')
    call Irreg%SaveState(FB)
    call FB%Close()
    call FB%Open('test.bin')
    call IrregLoad%loadState(FB)
    call FB%Close()
    if (abs( IrregLoad%Value(testval) -Irreg%Value(testval))<1e-8) then
        print *,'Load OK'
    else
        fails = fails+1
        print *, 'Load error'
    end if
    call FB%CreateFile('test.bin')
    call Reg%SaveState(FB)
    call FB%Close()
    call FB%Open('test.bin')
    call RegLoad%loadState(FB)
    call FB%Close()
    if (abs( RegLoad%Value(testval) -Reg%Value(testval))<1e-8) then
        print *,'Load OK'
    else
        fails = fails+1
        print *, 'Load error'
    end if
    call File%Delete('test.bin')
    
    end function RunInterpolationTests

    end module InterpolationTests