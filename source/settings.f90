    module settings
    use AMLutils
    use Random
    use IniFile
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
    real(mcp) :: AccuracyLevel = 1.
    !Set to >1 to use CAMB etc on higher accuracy settings.
    !Does not affect MCMC (except making it all slower)

    logical :: test_likelihood= .false.

    logical :: new_chains = .true.

    integer, parameter :: max_data_params = 100
    integer, parameter :: max_theory_params = 30
    integer, parameter :: max_num_params = max_theory_params + max_data_params

    logical :: use_fast_slow = .false.
    !Set to false if using a slow likelihood function so no there's point is treating
    !'fast' parameters differently (in fact, doing so will make performance worse)

    Type(TParamNames), save :: NameMapping

    integer, parameter :: sampling_metropolis = 1, sampling_slice = 2, sampling_fastslice =3, &
    sampling_slowgrid = 4,  sampling_multicanonical = 5,  sampling_wang_landau = 6, &
    sampling_fast_dragging = 7

    integer :: sampling_method = sampling_metropolis

    !scale of the proposal distribution in units of the posterior standard deviation
    real(mcp)    :: propose_scale  = 2.4_mcp

    !For fast dragging method, baseline number of intermediate drag steps
    real(mcp) :: dragging_steps = 3._mcp

    !The rest are set up automatically
    logical, parameter ::  generic_mcmc= .false.
    !set to true to not call CAMB, etc.
    !write GenericLikelihoodFunction in calclike.f90

    character(LEN=Ini_max_string_len) :: DataDir='data/'
    character(LEN=Ini_max_string_len) :: LocalDir='./'

    Type(TIniFile) :: CustomParams

    logical :: stop_on_error = .true. !whether to stop with error, or continue ignoring point

    integer :: num_theory_params, index_data, index_semislow=-1 !set later depending on datasets and theory parameterization

    integer, dimension(:), allocatable :: params_used
    integer num_params, num_params_used, num_data_params, num_fast, num_slow

    integer :: num_threads = 0
    integer :: instance = 0
    integer :: MPIchains = 1, MPIrank = 0

    logical :: Use_LSS = .false.
    logical :: Use_CMB = .false.

    integer :: logfile_unit  = 0
    integer :: outfile_handle = 0
    integer :: indepfile_handle = 0
    integer :: slow_proposals = 0
    integer :: output_lines = 0

    real(mcp), parameter :: logZero = 1e30_mcp
    character (LEN =1024) :: FileChangeIni = '', FileChangeIniAll = ''
    character(LEN=Ini_max_string_len) baseroot

    integer, parameter :: stdout = output_unit

    type mc_real_pointer
        real(mcp), dimension(:), pointer :: p
    end type mc_real_pointer


    contains

    function ReplaceDirs(S, repdir) result (filename)
    character(LEN=*), intent(in) :: S, repdir
    character(LEN=Ini_max_string_len) :: filename

    filename=S
    call StringReplace('%DATASETDIR%',repdir,filename)
    call StringReplace('%LOCALDIR%',LocalDir,filename)

    end function ReplaceDirs

    function ReadIniFileName(Ini,key, ADir, NotFoundFail) result (filename)
    Type(TIniFile) :: Ini
    character(LEN=*), intent(in) :: Key
    character(LEN=*), optional, intent(in) :: ADir
    character(LEN=Ini_max_string_len) :: filename, repdir
    logical, optional :: NotFoundFail
    integer i

    if (present(NotFoundFail)) then
        filename = Ini_Read_String_File(Ini, key, NotFoundFail)
    else
        filename = Ini_Read_String_File(Ini, key)
    end if
    if (present(ADir)) then
        repdir=ADir
    else
        repdir=DataDir
    end if
    filename= ReplaceDirs(filename, repdir)

    do i=1, CustomParams%L%Count
        call StringReplace('%'//trim(CustomParams%L%Items(i)%P%Name)//'%',&
        trim(ReplaceDirs(CustomParams%L%Items(i)%P%Value, repdir)) ,filename)
    end do

    end function ReadIniFileName

    subroutine CheckParamChangeF(F)
    character(LEN=*), intent(in) ::  F
    logical bad, doexit

    if (F /= '') then

    call Ini_Open(F, tmp_file_unit, bad, .false.)
    if (bad) return
    Ini_fail_on_not_found = .false.
    doexit = (Ini_Read_Int('exit',0) == 1)
    FeedBack = Ini_Read_Int('feedback',Feedback)
    num_threads = Ini_Read_Int('num_threads',num_threads)
    call Ini_Close
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
    integer j

    if (Feedback > 0) write(*,*) 'reading: '//trim(aname)

    call OpenTxtFile(aname, tmp_file_unit)

    do j=1,n
        read (tmp_file_unit,*, end = 200) vec(j)
    end do


    close(tmp_file_unit)
    return

200 write (*,*) 'vector file '//trim(aname)//' is the wrong size'
    stop

    end subroutine ReadVector

    subroutine WriteVector(aname, vec, n)
    character(LEN=*), intent(IN) :: aname
    integer, intent(in) :: n
    real(mcp), intent(in) :: vec(n)
    integer j

    call CreateTxtFile(aname, tmp_file_unit)

    do j=1,n
        write (tmp_file_unit,'(1E15.6)') vec(j)
    end do

    close(tmp_file_unit)

    end subroutine WriteVector



    subroutine ReadMatrix(aname, mat, m,n)
    character(LEN=*), intent(IN) :: aname
    integer, intent(in) :: m,n
    real(mcp), intent(out) :: mat(m,n)
    integer j,k
    real(mcp) tmp

    if (Feedback > 0) write(*,*) 'reading: '//trim(aname)
    call OpenTxtFile(aname, tmp_file_unit)

    do j=1,m
        read (tmp_file_unit,*, end = 200, err=100) mat(j,1:n)
    end do
    goto 120

100 rewind(tmp_file_unit)  !Try other possible format
    do j=1,m
        do k=1,n
            read (tmp_file_unit,*, end = 200) mat(j,k)
        end do
    end do

120 read (tmp_file_unit,*, err = 150, end =150) tmp
    goto 200

150 close(tmp_file_unit)
    return

200 write (*,*) 'matrix file '//trim(aname)//' is the wrong size'
    stop

    end subroutine ReadMatrix



    end module settings

