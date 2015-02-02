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

    subroutine ReadIni(this, Ini)
    class(TDatasetFileLikelihood) :: this
    class(TSettingIni) :: ini

    end subroutine ReadIni

    subroutine ReadDatasetFile(this, fname, OverrideSettings)
    class(TDatasetFileLikelihood) :: this
    character(LEN=*), intent(in) :: fname
    class(TSettingIni), optional, intent(in) :: OverrideSettings
    logical bad
    Type(TSettingIni) :: ini
    integer i_conflict

    call Ini%Open(fname,  bad, .false.)
    if (bad) then
        call MpiStop('Error opening dataset file '//trim(fname))
    end if
    if (present(OverrideSettings)) call Ini%Override(OverrideSettings)

    this%name = Ini%Read_String('name')
    if (this%name=='') this%name = File%ExtractName(fname, no_ext=.true.)

    this%num_conflicts = Ini%Read_Int('num_conflicts',0)
    allocate(this%conflict_name(this%num_conflicts))
    allocate(this%conflict_type(this%num_conflicts))
    do i_conflict=1,this%num_conflicts
        this%conflict_type(i_conflict) = Ini%Read_String(numcat('type_conflict',i_conflict))
        this%conflict_name(i_conflict) = Ini%Read_String(numcat('name_conflict',i_conflict))
    end do

    call this%ReadIni(Ini)

    call Ini%Close()

    end subroutine ReadDatasetFile

    recursive function Dataset_CheckConflicts(this, full_list) result(OK)
    !if for some reasons various likelihoods cannot be used at once
    !check here for conflicts after full list of likelihoods has been read in
    class(TDatasetFileLikelihood) :: this
    class(TLikelihoodList) :: full_list
    logical :: OK
    class(TDataLikelihood), pointer :: like_other
    integer i, i_conflict

    OK=.true.
    do i_conflict=1, this%num_conflicts
        do i= 1, full_list%count
            like_other => full_list%Item(i)
            if (like_other%LikelihoodType==trim(this%conflict_type(i_conflict)) &
            .and. like_other%name==trim(this%conflict_name(i_conflict))) then
                write(*,*) 'ERROR: Cannot use '//trim(this%LikelihoodType)//' dataset: '//trim(this%name)//&
                ' and '//trim(like_other%LikelihoodType)//' dataset: '&
                //trim(like_other%name)//' at the same time.'
                OK = .false.
            end if
        end do
    end do
    if (associated(this%CommonData)) then
        OK = this%CommonData%checkConflicts(full_list) .and. OK
    end if

    end function Dataset_CheckConflicts

    end module likelihood
