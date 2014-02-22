    module settings
    use RandUtils
    use FileUtils
    use StringUtils
    use MpiUtils
    use IniObjects
    use ParamNames
#ifdef f2003
    use, intrinsic :: iso_fortran_env, only : input_unit, output_unit,error_unit
#else
#define input_unit  5
#define output_unit 6
#define error_unit 0
#endif
    implicit none

#ifdef SINGLE
    integer, parameter :: mcp= KIND(1.0)
#else
    integer, parameter :: mcp= KIND(1.d0)
#endif
#ifdef MPI
#ifdef SINGLE
    integer, parameter :: MPI_real_mcp = MPI_REAL
#else
    integer, parameter :: MPI_real_mcp = MPI_DOUBLE_PRECISION
#endif
#endif
    integer,parameter :: time_dp = KIND(1.d0)


    double precision, parameter :: pi=3.14159265358979323846264338328d0, &
    twopi=2*pi, fourpi=4*pi
    double precision, parameter :: root2 = 1.41421356237309504880168872421d0, sqrt2 = root2
    double precision, parameter :: log2 = 0.693147180559945309417232121458d0

    real, parameter :: pi_r = 3.141592653, twopi_r = 2*pi_r, fourpi_r = twopi_r*2

    real(mcp), parameter :: const_c = 2.99792458e8_mcp

    logical :: use_fast_slow = .false.

    character(LEN=*), parameter :: CosmoMC_Version = 'Feb2014_oop'

    Type, extends(TIniFile) :: TSettingIni
    contains
    procedure :: FailStop => TSettingIni_FailStop
    procedure :: ReadFilename => TSettingIni_ReadFilename
    end type

    real(mcp) :: AccuracyLevel = 1.
    !Set to >1 to use theory calculation etc on higher accuracy settings.
    !Does not affect MCMC (except making it all slower)

    logical :: flush_write = .true.

    logical :: new_chains = .true.

    integer, parameter :: max_likelihood_functions = 50

    integer, parameter :: max_data_params = 100
    integer, parameter :: max_theory_params = 30
    integer, parameter :: max_num_params = max_theory_params + max_data_params

    !Set to false if using a slow likelihood function so no there's point is treating
    !'fast' parameters differently (in fact, doing so will make performance worse)

    integer, parameter :: sampling_metropolis = 1, sampling_slice = 2, sampling_fastslice =3, &
    sampling_slowgrid = 4,  sampling_multicanonical = 5,  sampling_wang_landau = 6, &
    sampling_fast_dragging = 7

    integer :: sampling_method = sampling_metropolis

    !For fast dragging method, baseline number of intermediate drag steps
    real(mcp) :: dragging_steps = 3._mcp

    !The rest are set up automatically
    logical  ::  generic_mcmc= .false.
    !set to true to not call CAMB, etc.
    !write GenericLikelihoodFunction in calclike.f90

    character(LEN=:), allocatable :: DataDir, LocalDir

    Type(TSettingIni), save :: CustomParams

    logical :: stop_on_error = .true. !whether to stop with error, or continue ignoring point

    integer :: num_theory_params, index_data, index_semislow=-1 !set later depending on datasets and theory parameterization

    integer, dimension(:), allocatable :: params_used
    integer num_params, num_params_used, num_data_params

    integer :: num_threads = 0
    integer :: instance = 0
    integer :: MPIchains = 1, MPIrank = 0

    logical :: checkpoint = .false.

    integer :: logfile_unit  = 0
    integer :: outfile_handle = 0
    integer :: output_lines = 0

    integer :: Feedback = 0

    real(mcp), parameter :: logZero = 1e30_mcp
    character (LEN =1024) :: FileChangeIni = '', FileChangeIniAll = ''
    character(LEN=:), allocatable :: baseroot, rootname

    integer, parameter :: stdout = output_unit

    type mc_real_pointer
        real(mcp), dimension(:), pointer :: p => null()
    end type mc_real_pointer


    contains

    subroutine InitializeGlobalSettingDefaults

    DataDir='data/'
    LocalDir='./'

    end subroutine InitializeGlobalSettingDefaults

    function ReplaceDirs(S, repdir) result (filename)
    character(LEN=*), intent(in) :: S, repdir
    character(LEN=:), allocatable :: filename

    filename=S
    call StringReplace('%DATASETDIR%',repdir,filename)
    call StringReplace('%LOCALDIR%',LocalDir,filename)

    end function ReplaceDirs


    subroutine DoStop(S, abort)
    character(LEN=*), intent(in), optional :: S
    integer ierror
    logical, intent(in), optional :: abort
    logical wantbort
    real MPI_StartTime

    if (outfile_handle/=0) close(outfile_handle)

    if (present(abort)) then
        wantbort = abort
    else
        wantbort = .false.
    end if

    if (present(S) .and. (wantbort .or. MPIRank==0)) write (*,*) trim(S)
#ifdef MPI
    MPI_StartTime = MPI_WTime() - MPI_StartTime
    if (Feedback > 0 .and. MPIRank==0) then
        write (*,*) 'Total time:', nint(MPI_StartTime), &
        '(',MPI_StartTime/(60*60),' hours)'
    end if
    ierror =0
    if (wantbort) then
        !Abort all in case other continuing chains want to communicate with us
        !in the case when max number of samples is reached
        call MPI_Abort(MPI_COMM_WORLD,ierror,ierror)
    else
        call mpi_finalize(ierror)
    end if
#endif

#ifdef DECONLY
    pause
#endif
    stop
    end subroutine DoStop

    subroutine TSettingIni_FailStop(L)
    class(TSettingIni) :: L

    call MpiStop()

    end subroutine TSettingIni_FailStop

    function TSettingIni_ReadFilename(Ini,key, ADir, NotFoundFail) result (OutName)
    class(TSettingIni) :: Ini
    character(LEN=*), intent(in) :: Key
    character(LEN=*), optional, intent(in) :: ADir
    character(LEN=:), allocatable :: filename, repdir
    character(LEN=:), allocatable :: OutName
    logical, optional :: NotFoundFail
    integer i

    if (present(NotFoundFail)) then
        filename = Ini%Read_String(key, NotFoundFail)
    else
        filename = Ini%Read_String(key)
    end if
    if (present(ADir)) then
        repdir=ADir
    else
        repdir=DataDir
    end if
    filename= ReplaceDirs(filename, repdir)

    do i=1, CustomParams%Count
        call StringReplace('%'//CustomParams%Items(i)%P%Name//'%',&
        trim(ReplaceDirs(CustomParams%Items(i)%P%Value, repdir)) ,filename)
    end do

    OutName = trim(filename)

    end function TSettingIni_ReadFilename

    subroutine CheckParamChangeF(F)
    character(LEN=*), intent(in) ::  F
    logical bad, doexit
    Type(TSettingIni) :: Ini

    if (F /= '') then
        call Ini%Open(F, bad, .false.)
        if (bad) return
        doexit = (Ini%Read_Int('exit',0) == 1)
        FeedBack = Ini%Read_Int('feedback',Feedback)
        num_threads = Ini%Read_Int('num_threads',num_threads)
        call Ini%Close()
        if (F== FileChangeIni) call DeleteFile(FileChangeini)
        if (doexit) call MpiStop('exit requested')
    end if

    end subroutine CheckParamChangeF

    subroutine CheckParamChange

    call CheckParamChangeF(FileChangeIni)
    if (FileChangeIni/=FileChangeIniAll) call CheckParamChangeF(FileChangeIniAll)

    end subroutine CheckParamChange

    subroutine ReadVector(aname, vec, n)
    character(LEN=*), intent(IN) :: aname
    integer, intent(in) :: n
    real(mcp), intent(out) :: vec(n)
    integer j, unit, status

    if (Feedback > 0) write(*,*) 'reading: '//trim(aname)

    unit = OpenNewTxtFile(aname)
    do j=1,n
        read (unit,*, iostat=status) vec(j)
        if (status/=0) then 
            write (*,*) 'vector file '//trim(aname)//' is the wrong size'
            stop
        end if
    end do
    close(unit)

    end subroutine ReadVector

    subroutine WriteVector(aname, vec, n)
    character(LEN=*), intent(IN) :: aname
    integer, intent(in) :: n
    real(mcp), intent(in) :: vec(n)
    integer j, unit

    unit = CreateNewTxtFile(aname)
    do j=1,n
        write (unit,'(1E15.6)') vec(j)
    end do
    close(unit)

    end subroutine WriteVector



    subroutine ReadMatrix(aname, mat, m,n)
    character(LEN=*), intent(IN) :: aname
    integer, intent(in) :: m,n
    real(mcp), intent(out) :: mat(m,n)
    integer j,k, unit
    real(mcp) tmp

    if (Feedback > 0) write(*,*) 'reading: '//trim(aname)
    unit = OpenNewTxtFile(aname)

    do j=1,m
        read (unit,*, end = 200, err=100) mat(j,1:n)
    end do
    goto 120

100 rewind(unit)  !Try other possible format
    do j=1,m
        do k=1,n
            read (unit,*, end = 200) mat(j,k)
        end do
    end do

120 read (unit,*, err = 150, end =150) tmp
    goto 200

150 close(unit)
    return

200 write (*,*) 'matrix file '//trim(aname)//' is the wrong size'
    stop

    end subroutine ReadMatrix


    function TimerTime()
    real(mcp) time
    real(time_dp) :: TimerTime
#ifdef MPI
    TimerTime = MPI_WTime()
#else
    call cpu_time(time)
    TimerTime=  time
#endif
    end function TimerTime

    subroutine Timer(Msg, start)
    character(LEN=*), intent(in), optional :: Msg
    real(time_dp), save :: timer_start
    real(time_dp), optional :: start
    real(time_dp) T

    if (present(start)) then
        T=start
    else
        T=timer_start
    end if

    if (present(Msg)) then
        write (*,*) trim(Msg)//': ', TimerTime() - T
    end if
    if (.not. present(start)) timer_start= TimerTime()

    end subroutine Timer

    subroutine DoAbort(S)
    character(LEN=*), intent(in), optional :: S
#ifdef MPI
    integer ierror
#endif
    if (present(S)) write (*,*) trim(S)
#ifdef MPI
    call MPI_Abort(MPI_COMM_WORLD,ierror,ierror)
#endif

#ifdef DECONLY
    pause
#endif
    stop
    end subroutine DoAbort


    end module settings