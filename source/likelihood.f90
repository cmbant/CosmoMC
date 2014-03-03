    module likelihood
    use settings
    use GeneralTypes
    implicit none
    private

    !TDataLikelihood is an instance of data of a particular kind, CMB, BOA, etc.
    !Can be multiple of same type
    !TDataLikelihood is defined in GeneralTypes, along with TLikelihoodList
    !TDatasetFileLikelihood extends this with some standard methods for reading settings
    !and checking for conflicts

    type, extends(TDataLikelihood) :: TDatasetFileLikelihood
        !likelihood that reads from a text file description, e.g. .dataset file
        !conflict names read from file and matched by type and name
        integer :: num_conflicts = 0
        character(LEN=LikeNameLen), allocatable, dimension(:) :: conflict_type
        character(LEN=LikeNameLen), allocatable, dimension(:) :: conflict_name
        class(TDatasetFileLikelihood), pointer :: CommonData => null()
    contains
    procedure :: ReadDatasetFile !open file, read standard things, call ReadIni
    procedure :: ReadIni  !read custom settings
    procedure :: checkConflicts => Dataset_CheckConflicts
    end type TDatasetFileLikelihood

    public TDatasetFileLikelihood, TLikelihoodList, TDataLikelihood
    contains

    !!!!!! TDatasetFileLikelihood

    subroutine ReadIni(like, Ini)
    class(TDatasetFileLikelihood) :: like
    class(TSettingIni) :: ini

    end subroutine ReadIni

    subroutine ReadDatasetFile(like, fname)
    class(TDatasetFileLikelihood) :: like
    character(LEN=*), intent(in) :: fname
    logical bad
    Type(TSettingIni) :: ini
    integer i_conflict

    call Ini%Open(fname,  bad, .false.)
    if (bad) then
        call MpiStop('Error opening dataset file '//trim(fname))
    end if

    like%name = Ini%Read_String('name')

    like%num_conflicts = Ini%Read_Int('num_conflicts',0)
    allocate(like%conflict_name(like%num_conflicts))
    allocate(like%conflict_type(like%num_conflicts))
    do i_conflict=1,like%num_conflicts
        like%conflict_type(i_conflict) = Ini%Read_String(numcat('type_conflict',i_conflict))
        like%conflict_name(i_conflict) = Ini%Read_String(numcat('name_conflict',i_conflict))
    end do

    call like%ReadIni(Ini)

    call Ini%Close()

    end subroutine ReadDatasetFile

    recursive function Dataset_CheckConflicts(like, full_list) result(OK)
    !if for some reasons various likelihoods cannot be used at once
    !check here for conflicts after full list of likelihoods has been read in
    class(TDatasetFileLikelihood) :: like
    class(TLikelihoodList) :: full_list
    logical :: OK
    class(TDataLikelihood), pointer :: like_other
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
