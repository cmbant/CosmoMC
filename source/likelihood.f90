    module likelihood
    !DataLikelihood is an instance of data of a particular kind, CMB, BOA, etc.
    !Can be multiple of same type
    use AMLUtils
    use settings
    use IniFile
    use ObjectLists, only: TObjectList
    use ParamNames
    implicit none

    integer, parameter :: max_likelihood_functions = 50
    integer, parameter :: LikeNameLen = 80

    type :: DataLikelihood
        integer :: speed = 0  !negative for slow likelihoods, larger positive for faster
        character(LEN=LikeNameLen) :: name = ''
        character(LEN=LikeNameLen) :: LikelihoodType= ''
        character(LEN=LikeNameLen) :: version = ''
        Type(TParamNames) :: nuisance_params
        !Internally calculated
        logical :: dependent_params(max_num_params) = .false.
        integer, allocatable :: nuisance_indices(:)
        integer, allocatable :: derived_indices(:)
        integer :: new_param_block_start, new_params
    contains
    procedure :: LogLike
    procedure :: LogLikeTheory !same as above when extra info not needed
    procedure :: loadParamNames
    procedure :: checkConflicts
    procedure :: derivedParameters  !derived parameters calculated by likelihood function
    procedure  :: WriteLikelihoodData
    end type DataLikelihood

    type, extends(DataLikelihood) :: DatasetFileLikelihood
        !likelihood that reads from a text file description, e.g. .dataset file
        !conflict names read from file and matched by type and name
        integer :: num_conflicts = 0
        character(LEN=LikeNameLen), pointer, dimension(:) :: conflict_type
        character(LEN=LikeNameLen), pointer, dimension(:) :: conflict_name
        class(DatasetFileLikelihood), pointer :: CommonData => null()
    contains
    procedure :: ReadDatasetFile !open file, read standard things, call ReadIni
    procedure :: ReadIni  !read custom settings
    procedure :: checkConflicts => Dataset_CheckConflicts
    end type DatasetFileLikelihood

    !This is the global list of likelihoods we will use
    Type, extends(TObjectList) :: LikelihoodList
        integer :: first_fast_param =0
        integer :: num_derived_parameters = 0
    contains
    procedure :: Item => LikelihoodItem
    procedure :: WriteLikelihoodContribs
    procedure :: AddNuisanceParameters
    procedure :: Compare => CompareLikes
    procedure :: checkAllConflicts
    procedure :: WriteDataForLikelihoods
    procedure :: addLikelihoodDerivedParams
    end type LikelihoodList

    Type(LikelihoodList), target, save :: DataLikelihoods

    contains

    function LikelihoodItem(L, i) result(P)
    Class(LikelihoodList) :: L
    integer, intent(in) :: i
    Class(DataLikelihood), pointer :: P

    select type (like => L%Items(i)%P)
    class is (DataLikelihood)
        P => like
        class default
        stop 'List contains non-DataLikelihood item'
    end select

    end function LikelihoodItem

    subroutine WriteLikelihoodContribs(L, aunit, likelihoods)
    Class(LikelihoodList) :: L
    integer, intent(in) :: aunit
    real(mcp), intent(in) :: likelihoods(*)
    integer i
    Class(DataLikelihood), pointer :: LikeItem

    do i=1,L%Count
        LikeItem =>  L%Item(i)
        write (aunit,'(2f11.3)',advance='NO') likelihoods(i),likelihoods(i)*2
        write(aunit,'(a)',advance='NO') '   '//trim(LikeItem%LikelihoodType)//': '//trim(LikeItem%name)
        if (LikeItem%Version/='') write(aunit,'(a)',advance='NO') ' '//trim(LikeItem%Version)
        write(aunit,'(a)') ''
    end do

    end subroutine WriteLikelihoodContribs

    subroutine WriteDataForLikelihoods(L, P, Theory, fileroot)
    Class(LikelihoodList) :: L
    real(mcp), intent(in) :: P(:)
    character(LEN=*), intent(in) :: fileroot
    class(*), intent(in) :: Theory
    integer i
    Class(DataLikelihood), pointer :: LikeItem

    do i=1,L%Count
        LikeItem => L%Item(i)
        call LikeItem%WriteLikelihoodData(Theory,P(LikeItem%nuisance_indices),fileroot)
    end do

    end subroutine WriteDataForLikelihoods

    integer function CompareLikes(this, R1, R2) result(comp)
    Class(LikelihoodList) :: this
    class(*) R1,R2

    select type (RR1 => R1)
    class is (DataLikelihood)
        select type (RR2 => R2)
        class is (DataLikelihood)
            comp = RR1%speed - RR2%speed
            return
        end select
    end select

    end function CompareLikes



    subroutine AddNuisanceParameters(L, Names)
    use ParamNames
    Class(LikelihoodList) :: L
    Type(TParamNames) :: Names
    Type(TParamNames), pointer :: NewNames
    Class(DataLikelihood), pointer :: DataLike
    integer i,j, baseDerived

    call L%Sort
    L%first_fast_param=0
    baseDerived = Names%num_derived
    do i=1,L%Count
        DataLike=>L%Item(i)
        NewNames => DataLike%nuisance_params
        if (Feedback>0 .and. MPIrank==0) print *,'adding parameters for: '//trim(DataLIke%name)
        DataLike%new_param_block_start = Names%num_MCMC +1
        !        if (DataLike%nuisance_params%num_derived>0) call MpiStop('No support for likelihood derived params yet')
        call ParamNames_Add(Names, NewNames)
        if (Names%num_MCMC > max_num_params) call MpiStop('increase max_data_params in settings.f90')
        DataLike%new_params = Names%num_MCMC - DataLike%new_param_block_start + 1
        allocate(DataLike%nuisance_indices(NewNames%num_MCMC))
        if (NewNames%num_MCMC/=0) then
            do j=1, NewNames%num_MCMC
                DataLike%nuisance_indices(j) = ParamNames_index(Names,NewNames%name(j))
            end do
            if (any(DataLike%nuisance_indices==-1)) call MpiStop('AddNuisanceParameters: unmatched data param')
            DataLike%dependent_params(DataLike%nuisance_indices) = .true.
            if (Feedback>1 .and. MPIrank==0) print *,trim(DataLike%name)//' data param indices:', DataLike%nuisance_indices
            if (L%first_fast_param==0 .and. DataLike%speed >=0 .and. &
            DataLike%new_params>0) L%first_fast_param = DataLike%new_param_block_start
        end if
    end do
    do i=1,L%Count
        !Add likelihood-derived parameters, after full set numbering has been dermined above
        DataLike=>L%Item(i)
        NewNames => DataLike%nuisance_params
        if (NewNames%num_derived>0) allocate(DataLike%derived_indices(NewNames%num_derived))
        do j=1, NewNames%num_derived
            DataLike%derived_indices(j) = ParamNames_index(Names, NewNames%name(j+NewNames%num_MCMC)) - Names%num_MCMC
        end do
        if (Feedback>1 .and. MPIrank==0) print *,trim(DataLike%name)//' derived param indices:', DataLike%derived_indices
        if (any(DataLike%derived_indices<=0)) &
        call MpiStop('AddNuisanceParameters: unmatched derived param')
    end do
    L%num_derived_parameters = Names%num_derived - baseDerived

    end subroutine AddNuisanceParameters

    subroutine checkAllConflicts(L)
    Class(LikelihoodList) :: L
    Class(DataLikelihood), pointer :: DataLike
    integer i

    do i=1,L%Count
        DataLike=>L%Item(i)
        if (.not. DataLike%checkConflicts(L)) &
        call MpiStop('Likelihood conflict reported by '//trim(DataLike%Name))
    end do

    end subroutine checkAllConflicts

    function addLikelihoodDerivedParams(L, P, Theory, derived) result(num_derived)
    Class(LikelihoodList) :: L
    Type(mc_real_pointer) :: derived
    class(*) :: Theory
    real(mcp) :: P(:)
    real(mcp), pointer :: allDerived(:)
    integer num_derived
    Class(DataLikelihood), pointer :: DataLike
    integer i

    num_derived = L%num_derived_parameters + size(derived%P)
    if (L%num_derived_parameters==0) return

    allocate(allDerived(num_derived))
    allDerived(1:size(derived%P)) = derived%P
    deallocate(derived%P)
    derived%P => allDerived

    do i=1,L%Count
        DataLike=>L%Item(i)
        if (allocated(DataLike%derived_indices)) then
            Derived%P(DataLike%derived_indices) = DataLike%derivedParameters(Theory, P(DataLike%nuisance_indices))
        end if
    end do

    end function addLikelihoodDerivedParams

    function derivedParameters(like, Theory, DataParams) result(derived)
    class(DataLikelihood) :: like
    class(*) :: Theory
    real(mcp) :: derived(like%nuisance_params%num_derived)
    real(mcp) :: DataParams(:)
    !Calculate any derived parameters internal to the likelihood that should be output
    !Number matches derived names defined in nuisance_params .paramnames file
    derived=0
    end function derivedParameters

    subroutine WriteLikelihoodData(like,Theory,DataParams, root)
    class(DataLikelihood) :: like
    class(*) :: Theory
    real(mcp), intent(in) :: DataParams(:)
    character(LEN=*), intent(in) :: root
    !Write out any derived data that might be useful for the likelihood (e.g. foreground model)
    end subroutine WriteLikelihoodData

    function logLikeTheory(like, CMB)
    !For likelihoods that don't need Theory or DataParams
    class(DataLikelihood) :: like
    class(*) :: CMB
    real(mcp) LogLikeTheory

    logLikeTheory= logZero
    stop 'logLikeTheory or logLike should not be overridden'
    end function

    function LogLike(like, CMB, Theory, DataParams)
    class(DataLikelihood) :: like
    class(*) :: CMB
    class(*) :: Theory
    real(mcp) :: DataParams(:)
    real(mcp) LogLike

    logLike = like%logLikeTheory(CMB)
    end function


    subroutine loadParamNames(like, fname)
    class(DataLikelihood) :: like
    character(LEN=*), intent(in) :: fname

    call ParamNames_init(like%nuisance_params, fname)

    end subroutine loadParamNames

    function checkConflicts(like, full_list) result(OK)
    !if for some reasons various likelihoods cannot be used at once
    !check here for conflicts after full list of likelihoods has been read in
    class(DataLikelihood) :: like
    class(LikelihoodList) :: full_list
    logical :: OK

    OK=.true.

    end function checkConflicts

    !!!!!! DatasetFileLikelihood

    subroutine ReadIni(like, Ini)
    class(DatasetFileLikelihood) :: like
    Type(TIniFile) :: ini

    end subroutine ReadIni

    subroutine ReadDatasetFile(like, fname)
    class(DatasetFileLikelihood) :: like
    character(LEN=*), intent(in) :: fname
    logical bad
    integer file_unit
    Type(TIniFile) :: ini
    integer i_conflict

    file_unit = new_file_unit()
    call Ini_Open_File(Ini, fname, file_unit, bad, .false.)
    if (bad) then
        call MpiStop('Error opening dataset file '//trim(fname))
    end if

    like%name = Ini_Read_String_File(Ini,'name')

    like%num_conflicts = Ini_Read_Int_File(Ini,'num_conflicts',0)
    allocate(like%conflict_name(like%num_conflicts))
    allocate(like%conflict_type(like%num_conflicts))
    do i_conflict=1,like%num_conflicts
        like%conflict_type(i_conflict) = Ini_Read_String_File(Ini,numcat('type_conflict',i_conflict))
        like%conflict_name(i_conflict) = Ini_Read_String_File(Ini,numcat('name_conflict',i_conflict))
    end do

    call like%ReadIni(Ini)

    call Ini_Close_File(Ini)
    call ClearFileUnit(file_unit)

    end subroutine ReadDatasetFile

    recursive function Dataset_CheckConflicts(like, full_list) result(OK)
    !if for some reasons various likelihoods cannot be used at once
    !check here for conflicts after full list of likelihoods has been read in
    class(DatasetFileLikelihood) :: like
    class(LikelihoodList) :: full_list
    logical :: OK
    class(DataLikelihood), pointer :: like_other
    integer i, i_conflict

    OK=.true.
    do i_conflict=1, like%num_conflicts
        do i= 1, full_list%count
            like_other => full_list%Item(i)
            if (like_other%LikelihoodType==trim(like%conflict_type(i_conflict)) &
            .and. like_other%name==trim(like%conflict_name(i_conflict))) then
                write(*,*) 'ERROR: Cannot use '//trim(like%LikelihoodType)//' dataset: '//trim(like%name)//&
                ' and '//trim(like_other%LikelihoodType)//' dataset: '&
                //trim(like_other%name)//' at the same time.'
                OK = .false.
            end if
        end do
    end do
    if (associated(like%CommonData)) then
        OK = like%CommonData%checkConflicts(full_list) .and. OK
    end if

    end function Dataset_CheckConflicts

    end module likelihood
