
    module GeneralTypes
    use ObjectLists
    use IO
    implicit none
    private

    type int_arr
        integer, dimension(:), allocatable :: p
    end type int_arr

    type, extends(TSaveLoadStateObject) :: TCheckpointable
    contains
    procedure :: ReadParams => TCheckpointable_ReadParams
    end Type

    type TObjectWithParams
    contains
    procedure :: ReadParams => TObjectWithParams_ReadParams
    end Type

    Type TTheoryParams
        real(mcp) :: BaseParams(max_theory_params) = 0._mcp
    end Type TTheoryParams

    Type TTheoryIntermediateCache
        !Cache of intermediate steps in the theory calculation that may be (partly) reused for other points
    contains
    procedure :: Clear => TTheoryIntermediateCache_Clear
    end Type TTheoryIntermediateCache

    !Each descendant of the base TConfigClass below has a Config element which point to a TGeneralConfig instance taht determines
    !which related class implementations are actually used
    Type, extends(TObjectWithParams) :: TGeneralConfig
        class(TParameterization), pointer :: Parameterization => null()
        class(TTheoryCalculator), pointer :: Calculator => null()
    contains
    procedure :: SetTheoryParameterization => TGeneralConfig_SetTheoryParameterization
    procedure :: SetParameterizationName => TGeneralConfig_SetParameterizationName
    procedure :: NewTheory => TGeneralConfig_NewTheory
    procedure :: InitForLikelihoods => TGeneralConfig_InitForLikelihoods
    end Type

    Type, extends(TCheckpointable) :: TConfigClass
        class(TGeneralConfig), pointer :: Config
    contains
    procedure :: InitConfig  => TConfigClass_InitConfig
    procedure :: InitWithParams  => TConfigClass_InitWithParams
    generic:: Init => InitWithParams,InitConfig
    end Type

    Type, extends(TConfigClass) :: TTheoryPredictions
        !Actual computed theory predictions used by likelihoods
        !Config%Calculator can in some cases be used to provide theory functions
    contains
    procedure :: Clear => TTheoryPredictions_Clear
    procedure :: WriteTheory
    procedure :: ReadTheory
    procedure :: WriteTextData
    end Type TTheoryPredictions

    Type, extends(TConfigClass) :: TTheoryCalculator
        character(LEN=128) :: calcName = 'TheoryCalculator'
    contains
    procedure :: InitForLikelihoods => TTheoryCalculator_InitForLikelihoods
    procedure :: Error => TTheoryCalculator_Error
    procedure :: ErrorNotImplemented => TTheoryCalculator_ErrorNotImplemented
    procedure :: VersionTraceOutput => TTheoryCalculator_VersionTraceOutput
    procedure :: ReadImportanceParams => TTheoryCalculator_ReadImportanceParams
    procedure :: GetTheoryForImportance => TTheoryCalculator_GetTheoryForImportance
    end Type TTheoryCalculator

    type, extends(TConfigClass)  :: TParameterization
    contains
    procedure :: ParamArrayToTheoryParams
    procedure :: CalcDerivedParams=> TParameterization_CalcDerivedParams
    procedure :: Initialize => TParameterization_Initialize
    procedure :: NonBaseParameterPriors => TParameterization_NonBaseParameterPriors
    procedure :: NewTheoryParams => TParameterization_NewTheoryParams
    end type TParameterization

    Type, extends(TParameterization) :: GenericParameterization
    end type GenericParameterization

    Type TCalculationAtParamPoint
        !Parameter values, calculated theory and intermediates and likelihoods at a particular point in parameter space
        real(mcp) :: P(max_num_params)
        real(mcp) :: likelihoods(max_likelihood_functions)
        real(mcp) :: lastParamArray(max_num_params)
        logical :: validInfo = .false.
        class(TTheoryPredictions), allocatable :: Theory
        class(TTheoryIntermediateCache), pointer :: Info => null()
    contains
    procedure :: Clear => TCalculationAtParamPoint_Clear
    procedure :: WriteParams => TCalculationAtParamPoint_WriteParams
    procedure :: AcceptReject => TCalculationAtParamPoint_AcceptReject
    procedure :: ReadModel
    procedure :: WriteModel
    end Type TCalculationAtParamPoint

    integer, parameter :: LikeNameLen = 80
    integer, parameter :: LikeTagLen = 20

    !Would like to have likelihood stuff in separate file, but fortran module dependency rules make that very difficult.
    type, extends(TConfigClass) :: TDataLikelihood
        integer :: speed = 0  !negative for slow likelihoods, larger positive for faster
        character(LEN=LikeNameLen) :: name = ''
        character(LEN=:), allocatable :: tag !Short name used to tag output likelihoods
        character(LEN=LikeNameLen) :: LikelihoodType= ''
        character(LEN=LikeNameLen) :: version = ''
        Type(TParamNames) :: nuisance_params
        !Internally calculated
        integer :: original_index = 0
        logical :: dependent_params(max_num_params) = .false.
        integer, allocatable :: nuisance_indices(:)
        integer, allocatable :: derived_indices(:)
        integer :: new_param_block_start, new_params
    contains
    procedure :: InitConfig => TDataLikelihood_InitConfig
    procedure :: GetLogLike => TDataLikelihood_GetLogLike
    procedure :: WriteLikelihoodData => TDataLikelihood_WriteLikelihoodData
    procedure :: derivedParameters  => TDataLikelihood_derivedParameters
    procedure :: loadParamNames => TDataLikelihood_loadParamNames
    procedure :: checkConflicts => TDataLikelihood_checkConflicts
    procedure :: GetTag => TDataLikelihood_GetTag
    end type TDataLikelihood

    !This is the global list of likelihoods we will use
    Type, extends(TObjectList) :: TLikelihoodList
        integer :: first_fast_param =0
        integer :: num_derived_parameters = 0
        Type(TIntegerArrayList) :: LikelihoodTypeIndices
        integer, allocatable :: original_order(:)
    contains
    procedure :: Item => LikelihoodItem
    procedure :: WriteLikelihoodContribs
    procedure :: AddNuisanceParameters
    procedure :: AddOutputLikelihoodParams
    procedure :: Compare => CompareLikes
    procedure :: checkAllConflicts
    procedure :: WriteDataForLikelihoods
    procedure :: addLikelihoodDerivedParams
    procedure :: OutputDescription
    end type TLikelihoodList

    Type(TLikelihoodList), target, save :: DataLikelihoods

    public int_arr,TCheckpointable, TTheoryParams, TTheoryIntermediateCache, TCalculationAtParamPoint, TGeneralConfig, &
        & TConfigClass, TTheoryPredictions, TTheoryCalculator, TParameterization, GenericParameterization, &
        & TDataLikelihood, TLikelihoodList, DataLikelihoods, LikeNameLen

    contains


    subroutine ParamArrayToTheoryParams(this, Params, CMB)
    class(TParameterization) :: this
    real(mcp) Params(:)
    Class(TTheoryParams), target :: CMB

    CMB%BaseParams(1:num_theory_params) = Params(1:num_theory_params)

    end subroutine ParamArrayToTheoryParams

    subroutine TParameterization_Initialize(this, Ini, Names, DefaultName, Config)
    class(TParameterization) :: this
    class(TSettingIni) :: Ini
    class(TParamNames) :: Names
    character(LEN=*), intent(in) :: DefaultName
    class(TGeneralConfig), target :: Config
    character(LEN=:), allocatable :: ParamNamesFile

    call this%TConfigClass%Init(Config)

    ParamNamesFile = Ini%ReadFileName('ParamNamesFile', NotFoundFail=.false.)
    if (ParamNamesFile =='' .and. DefaultName/='') ParamNamesFile= trim(LocalDir)//trim(DefaultName)

    if (ParamNamesFile /='') then
        call Names%init(ParamNamesFile)
        num_theory_params= Names%num_MCMC
    else
        num_theory_params= Ini%Read_Int('num_theory_params')
        call Names%SetUnnamed(num_theory_params)
    end if
    if (num_theory_params> max_theory_params) call MpiStop('see settings.f90: num_theory_params> max_theory_params')
    index_data =  num_theory_params+1

    end subroutine TParameterization_Initialize


    function TParameterization_NonBaseParameterPriors(this,CMB)
    class(TParameterization) :: this
    class(TTheoryParams) :: CMB
    real(mcp):: TParameterization_NonBaseParameterPriors

    TParameterization_NonBaseParameterPriors = 0

    end function TParameterization_NonBaseParameterPriors

    subroutine TParameterization_CalcDerivedParams(this, P, Theory, derived)
    class(TParameterization) :: this
    real(mcp), allocatable :: derived(:)
    class(TTheoryPredictions), allocatable :: Theory !can be unallocated for simple cases (e.g. generic)
    real(mcp) :: P(:)


    end subroutine TParameterization_CalcDerivedParams

    subroutine TParameterization_NewTheoryParams(this,TheoryParams)
    class(TParameterization) :: this
    class(TTheoryParams), allocatable :: TheoryParams

    allocate(TTheoryParams::TheoryParams)

    end subroutine TParameterization_NewTheoryParams


    !!!TTheoryPredictions

    subroutine TTheoryPredictions_Clear(this)
    class(TTheoryPredictions) :: this
    end subroutine TTheoryPredictions_Clear

    !!! TCalculationAtParamPoint

    subroutine TCalculationAtParamPoint_Clear(this, Keep)
    class(TCalculationAtParamPoint) :: this, Keep

    if (associated(this%Info) .and. .not. associated(this%Info, Keep%Info)) then
        call this%Info%Clear()
        deallocate(this%Info)
        nullify(This%Info)
    end if
    if (allocated(this%Theory)) then
        call this%Theory%Clear()
        deallocate(this%Theory)
    end if
    end subroutine TCalculationAtParamPoint_Clear


    subroutine TCalculationAtParamPoint_AcceptReject(this, Trial, accpt)
    !Handle freeing of memory of internal info: if accpt then clear this, otherwise clear Trial
    class(TCalculationAtParamPoint) :: this
    class(TCalculationAtParamPoint) :: Trial
    logical, intent(in) :: accpt

    if (accpt) then
        call this%Clear(keep = Trial)
    else
        call Trial%Clear(keep = this)
    end if

    end subroutine TCalculationAtParamPoint_AcceptReject

    subroutine TCalculationAtParamPoint_WriteParams(this, Config, mult, like)
    class(TCalculationAtParamPoint) this
    class(TGeneralConfig) :: Config
    real(mcp), intent(in) :: mult, like
    real(mcp), allocatable :: output_array(:), derived(:)
    integer :: numderived = 0

    if (ChainOutFile%unit==0) return

    call Config%Parameterization%CalcDerivedParams(this%P, this%Theory, derived)
    call DataLikelihoods%addLikelihoodDerivedParams(this%P, this%Theory, derived, this%Likelihoods, like)

    if (allocated(derived)) numderived = size(derived)

    allocate(output_array(num_params_used+numderived))
    output_array(:num_params_used) = this%P(params_used)
    if (numderived>0) output_array(num_params_used+1:) = derived

    call IO_OutputChainRow(ChainOutFile, mult, like, output_array)

    end subroutine TCalculationAtParamPoint_WriteParams

    subroutine WriteModel(this, F, like, mult)
    Class(TCalculationAtParamPoint) :: this
    class(TFileStream) :: F
    real(mcp), intent(in) :: mult, like
    end subroutine WriteModel

    subroutine  ReadModel(this,  F, has_likes, mult, like, error)
    Class(TCalculationAtParamPoint) :: this
    class(TFileStream) :: F
    integer, intent(out) :: error
    real(mcp), intent(out) :: mult, like
    logical, intent(out) :: has_likes(:)

    mult=0
    like=0
    error=0
    has_likes=.false.
    call MpiStop('ReadModel not implemented')

    end subroutine ReadModel


    !!! TTheoryIntermediateCache

    subroutine TTheoryIntermediateCache_Clear(Info)
    class(TTheoryIntermediateCache) Info

    end subroutine TTheoryIntermediateCache_Clear

    !!! TTheoryCalculator

    subroutine TTheoryCalculator_InitForLikelihoods(this)
    class(TTheoryCalculator) :: this

    !Called after likelihoods etc loaded
    end subroutine TTheoryCalculator_InitForLikelihoods


    subroutine TTheoryCalculator_ReadImportanceParams(this, Ini)
    class(TTheoryCalculator) :: this
    class(TSettingIni) :: Ini

    end subroutine TTheoryCalculator_ReadImportanceParams

    subroutine TTheoryCalculator_GetTheoryForImportance(this, CMB, Theory, error)
    class(TTheoryCalculator) :: this
    class(TTheoryParams) :: CMB
    class(TTheoryPredictions) :: Theory
    integer error

    error=0
    !calculate power spectra from scratch (for importance sampling)
    !Theory may already be calculated, so only fill in missing bits (DoCls, DoPk) + derived
    call this%ErrorNotImplemented('GetTheoryForImportance')

    end subroutine TTheoryCalculator_GetTheoryForImportance


    subroutine TTheoryCalculator_VersionTraceOutput(this, ReadValues)
    class(TTheoryCalculator) :: this
    class(TNameValueList) :: ReadValues

    !Store for the record any useful info about version etc.

    end subroutine TTheoryCalculator_VersionTraceOutput


    subroutine TTheoryCalculator_ErrorNotImplemented(this,S)
    class(TTheoryCalculator) :: this
    character(LEN=*), intent(in) :: S

    call this%Error('Not implemented: '//trim(S))

    end subroutine TTheoryCalculator_ErrorNotImplemented


    subroutine TTheoryCalculator_Error(this,S)
    class(TTheoryCalculator) :: this
    character(LEN=*), intent(in) :: S

    call MpiStop(trim(this%CalcName)//': Error: '//trim(S))

    end subroutine TTheoryCalculator_Error


    !!TConfigClass
    subroutine TConfigClass_InitConfig(this, Config)
    class(TConfigClass) :: this
    class(TGeneralConfig), target :: config

    this%Config => config

    end subroutine TConfigClass_InitConfig


    subroutine TConfigClass_InitWithParams(this, Ini, Config)
    class(TConfigClass) :: this
    class(TSettingIni) :: Ini
    class(TGeneralConfig), target :: Config

    call this%InitConfig(Config)
    call this%ReadParams(Ini)

    end subroutine TConfigClass_InitWithParams

    subroutine WriteTheory(this, F, first)
    class(TTheoryPredictions) this
    class(TFileStream) :: F
    logical, intent(in) :: first
    end subroutine WriteTheory

    subroutine ReadTheory(this, F, first)
    class(TTheoryPredictions) this
    class(TFileStream) :: F
    logical, intent(in) :: first
    end subroutine ReadTheory

    subroutine WriteTextData(this,fnameroot)
    class(TTheoryPredictions) this
    character(LEN=*), intent(in) :: fnameroot
    end subroutine WriteTextData


    !!TGeneralConfig
    subroutine TGeneralConfig_SetTheoryParameterization(this, Ini, Names, defaultParam)
    class(TGeneralConfig) :: this
    class(TSettingIni) :: Ini
    class(TParamNames) :: Names
    character(LEN=*), intent(in) :: defaultParam
    character(LEN=:), allocatable :: paramtxt

    paramtxt = Ini%Read_String_Default('parameterization', defaultParam)
    if (.not. this%SetParameterizationName(paramtxt,Ini, Names)) then
        call MpiStop('GeneralConfig: unknown parameterization :'//trim(paramtxt))
    end if

    end subroutine TGeneralConfig_SetTheoryParameterization

    function TGeneralConfig_SetParameterizationName(this, nametag, Ini, Names) result(OK)
    class(TGeneralConfig) :: this
    character(LEN=*), intent(in) :: nametag
    class(TSettingIni) :: Ini
    class(TParamNames) :: Names
    logical OK
    Type(GenericParameterization), pointer :: GenParam

    OK = .true.

    if (nametag =='generic') then
        allocate(GenParam)
        this%Parameterization => GenParam
        call GenParam%Initialize(Ini,Names,'',this)
    else
        OK = .false.
    end if

    end function TGeneralConfig_SetParameterizationName

    subroutine TGeneralConfig_NewTheory(this, Theory)
    class(TGeneralConfig) :: this
    class(TTheoryPredictions), allocatable :: Theory

    allocate(TTheoryPredictions::Theory)
    call Theory%Init(this)

    end subroutine TGeneralConfig_NewTheory

    subroutine TGeneralConfig_InitForLikelihoods(this)
    class(TGeneralConfig) :: this
    class(TDataLikelihood), pointer :: DataLike
    integer i

    if (associated(this%Calculator)) call this%Calculator%InitForLikelihoods()
    do i=1,DataLikelihoods%Count
        DataLike=>DataLikelihoods%Item(i)
        call DataLike%InitConfig(this)
    end do

    end subroutine TGeneralConfig_InitForLikelihoods

    !!!TCheckpointable

    subroutine TCheckpointable_ReadParams(this, Ini)
    class(TCheckpointable) :: this
    class(TSettingIni) :: Ini
    end subroutine TCheckpointable_ReadParams

    !!!TObjectWithParams

    subroutine TObjectWithParams_ReadParams(this, Ini)
    class(TObjectWithParams) :: this
    class(TSettingIni) :: Ini
    end subroutine TObjectWithParams_ReadParams


    !!!TDataLikelihood

    function TDataLikelihood_GetLogLike(this, Params, Theory, DataParams) result(LogLike)
    class(TDataLikelihood) :: this
    class(TTheoryParams) :: Params
    class(TTheoryPredictions) :: Theory
    real(mcp) :: DataParams(:)
    real(mcp) LogLike

    stop 'GetLogLike should not be overridden'
    logLike = LogZero

    end function TDataLikelihood_GetLogLike


    subroutine TDataLikelihood_WriteLikelihoodData(this,Theory,DataParams, root)
    class(TDataLikelihood) :: this
    class(TTheoryPredictions) :: Theory
    real(mcp), intent(in) :: DataParams(:)
    character(LEN=*), intent(in) :: root
    !Write out any derived data that might be useful for the likelihood (e.g. foreground model)
    end subroutine TDataLikelihood_WriteLikelihoodData


    subroutine TDataLikelihood_InitConfig(this, Config)
    class(TDataLikelihood) :: this
    class(TGeneralConfig), target :: config

    call this%TConfigClass%InitConfig(Config)

    end subroutine TDataLikelihood_InitConfig


    function TDataLikelihood_derivedParameters(this, Theory, DataParams) result(derived)
    class(TDataLikelihood) :: this
    class(TTheoryPredictions) :: Theory
    real(mcp) :: derived(this%nuisance_params%num_derived)
    real(mcp) :: DataParams(:)
    !Calculate any derived parameters internal to the likelihood that should be output
    !Number matches derived names defined in nuisance_params .paramnames file
    derived=0
    end function TDataLikelihood_derivedParameters

    subroutine TDataLikelihood_loadParamNames(this, fname)
    class(TDataLikelihood) :: this
    character(LEN=*), intent(in) :: fname

    call this%nuisance_params%init(fname)

    end subroutine TDataLikelihood_loadParamNames

    function TDataLikelihood_checkConflicts(this, full_list) result(OK)
    !if for some reasons various likelihoods cannot be used at once
    !check here for conflicts after full list of likelihoods has been read in
    class(TDataLikelihood) :: this
    class(TLikelihoodList) :: full_list
    logical :: OK

    OK=.true.

    end function TDataLikelihood_checkConflicts

    function TDataLikelihood_GetTag(this) result(tag)
    class(TDataLikelihood) :: this
    character(LEN=:), allocatable :: tag

    if (allocated(this%tag)) then
        tag = this%tag
    else
        tag = trim(this%Name)
    end if

    end function TDataLikelihood_GetTag


    !!!TLikelihoodList

    function LikelihoodItem(L, i) result(P)
    Class(TLikelihoodList) :: L
    integer, intent(in) :: i
    Class(TDataLikelihood), pointer :: P

    select type (like => L%Items(i)%P)
    class is (TDataLikelihood)
        P => like
        class default
        stop 'List contains non-TDataLikelihood item'
    end select

    end function LikelihoodItem

    subroutine WriteLikelihoodContribs(L, aunit, likelihoods)
    Class(TLikelihoodList) :: L
    integer, intent(in) :: aunit
    real(mcp), intent(in) :: likelihoods(*)
    integer i, ix
    Class(TDataLikelihood), pointer :: LikeItem
    character(LEN=:), allocatable :: tagname

    do ix=1,L%Count
        i = L%Original_order(ix)
        LikeItem =>  L%Item(i)
        write (aunit,'(2f11.3)',advance='NO') likelihoods(i),likelihoods(i)*2
        tagname = trim (LikeItem%name)
        if (allocated(LikeItem%tag)) then
            if (LikeItem%tag /= LikeItem%name) tagname = LikeItem%tag //' = '//tagname
        end if
        write(aunit,'(a)',advance='NO') '   '//trim(LikeItem%LikelihoodType)//': '//tagname
        if (LikeItem%Version/='') write(aunit,'(a)',advance='NO') ' '//trim(LikeItem%Version)
        write(aunit,'(a)') ''
    end do

    end subroutine WriteLikelihoodContribs


    subroutine WriteDataForLikelihoods(L, P, Theory, fileroot)
    Class(TLikelihoodList) :: L
    real(mcp), intent(in) :: P(:)
    character(LEN=*), intent(in) :: fileroot
    class(TTheoryPredictions), intent(in) :: Theory
    integer i
    Class(TDataLikelihood), pointer :: LikeItem

    do i=1,L%Count
        LikeItem => L%Item(i)
        call LikeItem%WriteLikelihoodData(Theory,P(LikeItem%nuisance_indices),fileroot)
    end do

    end subroutine WriteDataForLikelihoods


    integer function CompareLikes(this, R1, R2) result(comp)
    Class(TLikelihoodList) :: this
    class(*) R1,R2

    select type (RR1 => R1)
    class is (TDataLikelihood)
        select type (RR2 => R2)
        class is (TDataLikelihood)
            comp = RR1%speed - RR2%speed
            return
        end select
    end select

    end function CompareLikes


    subroutine AddNuisanceParameters(L, Names)
    Class(TLikelihoodList) :: L
    Type(TParamNames) :: Names
    Type(TParamNames), pointer :: NewNames
    Class(TDataLikelihood), pointer :: DataLike
    integer i,j, baseDerived

    do i=1,L%Count
        DataLike=>L%Item(i)
        DataLike%original_index = i
    end do
    allocate(L%Original_order(L%Count))
    call L%Sort
    L%first_fast_param=0
    baseDerived = Names%num_derived
    do i=1,L%Count
        DataLike=>L%Item(i)
        L%Original_order(DataLike%Original_index) = i
        NewNames => DataLike%nuisance_params
        if (Feedback>0 .and. MPIrank==0) print *,'adding parameters for: '//trim(DataLIke%name)
        DataLike%new_param_block_start = Names%num_MCMC +1
        call Names%Add(NewNames)
        if (Names%num_MCMC > max_num_params) call MpiStop('increase max_data_params in settings.f90')
        DataLike%new_params = Names%num_MCMC - DataLike%new_param_block_start + 1
        allocate(DataLike%nuisance_indices(NewNames%num_MCMC))
        if (NewNames%num_MCMC/=0) then
            do j=1, NewNames%num_MCMC
                DataLike%nuisance_indices(j) = Names%index(NewNames%name(j))
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
        if (NewNames%num_derived>0) then
            allocate(DataLike%derived_indices(NewNames%num_derived))
            do j=1, NewNames%num_derived
                DataLike%derived_indices(j) = Names%index(NewNames%name(j+NewNames%num_MCMC)) - Names%num_MCMC
            end do
            if (Feedback>1 .and. MPIrank==0) print *,trim(DataLike%name)//' derived param indices:', DataLike%derived_indices
            if (any(DataLike%derived_indices<=0)) call MpiStop('AddNuisanceParameters: unmatched derived param')
        end if
    end do
    L%num_derived_parameters = Names%num_derived - baseDerived

    end subroutine AddNuisanceParameters

    subroutine AddOutputLikelihoodParams(L, Names)
    Class(TLikelihoodList) :: L
    Type(TParamNames) :: Names, LikeNames
    integer i, j, ix, like_sum_ix
    class(TDataLikelihood), pointer :: Like
    character(LEN=:), pointer :: atype
    character(LEN=:), allocatable :: tag
    integer, allocatable :: counts(:), indices(:)
    Type(TStringList) :: LikelihoodTypes

    if (L%Count==0) return
    call LikeNames%Alloc(L%Count+1)
    allocate(counts(L%Count), source=0)
    do i=1, L%Count
        Like => L%Item(i)
        tag = Like%GetTag()
        LikeNames%name(Like%Original_index) = 'chi2_'//tag
        LikeNames%label(Like%Original_index) = FormatString(trim(chisq_label), StringEscape(tag,'_'))
        LikeNames%is_derived(Like%Original_index) = .true.
        if (Like%LikelihoodType/='') then
            ix = LikelihoodTypes%IndexOf(Like%LikelihoodType)
            if (ix==-1) then
                call LikelihoodTypes%Add(trim(Like%LikelihoodType))
                counts(LikelihoodTypes%Count)=1
            else
                counts(ix) = counts(ix)+1
            end if
        end if
    end do
    !Add a parameter for the prior
    LikeNames%name(L%Count+1) = 'chi2_prior'
    LikeNames%label(L%Count+1) = FormatString(trim(chisq_label), 'prior')
    LikeNames%is_derived(L%Count+1) = .true.

    call Names%Add(LikeNames,check_duplicates=.true.)

    !Add a derived parameters which are sums of all likelihoods of a given type (e.g. CMB, BAO, etc..)
    call LikeNames%Alloc(count(counts(:LikelihoodTypes%Count)>1))
    like_sum_ix = 0
    do i=1, LikelihoodTypes%Count
        if (counts(i)>1) then
            like_sum_ix = like_sum_ix +1
            allocate(indices(counts(i)))
            ix=1
            atype => LikelihoodTypes%Item(i)
            do j=1, L%Count
                Like => L%Item(j)
                if (Like%LikelihoodType == atype) then
                    indices(ix) = j
                    ix = ix +1
                end if
            end do
            call L%LikelihoodTypeIndices%Add(indices)
            LikeNames%name(like_sum_ix) = 'chi2_'//atype
            LikeNames%label(like_sum_ix) = FormatString(trim(chisq_label), StringEscape(trim(atype),'_'))
            LikeNames%is_derived(like_sum_ix) = .true.
            deallocate(indices)
        end if
    end do
    call Names%Add(LikeNames,check_duplicates=.true.)

    end subroutine AddOutputLikelihoodParams

    subroutine checkAllConflicts(L)
    Class(TLikelihoodList) :: L
    Class(TDataLikelihood), pointer :: DataLike
    integer i

    do i=1,L%Count
        DataLike=>L%Item(i)
        if (.not. DataLike%checkConflicts(L)) &
            call MpiStop('Likelihood conflict reported by '//trim(DataLike%Name))
    end do

    end subroutine checkAllConflicts

    subroutine addLikelihoodDerivedParams(L, P, Theory, derived, Likelihoods, logLike)
    class(TLikelihoodList) :: L
    real(mcp), allocatable :: derived(:)
    real(mcp), intent(in), optional :: Likelihoods(:), logLike
    class(TTheoryPredictions), allocatable :: Theory
    real(mcp) :: P(:)
    real(mcp), allocatable :: allDerived(:)
    Class(TDataLikelihood), pointer :: DataLike
    integer i
    integer :: num_in
    integer :: num_derived

    num_in=0
    num_derived=0
    if (allocated(derived)) num_in = size(derived)
    num_derived = L%num_derived_parameters + num_in
    if (L%num_derived_parameters >=0) then
        allocate(allDerived(num_derived))
        if (num_in > 0) allDerived(1:num_in) = derived
        call move_alloc(allDerived, derived)

        do i=1,L%Count
            DataLike=>L%Item(i)
            if (allocated(DataLike%derived_indices)) then
                Derived(DataLike%derived_indices) = DataLike%derivedParameters(Theory, P(DataLike%nuisance_indices))
            end if
        end do
    end if

    if (present(Likelihoods) .and. L%Count>0) then
        if (.not. present(logLike)) call MpiStop('Must have logLike in addLikelihoodDerivedParams')
        allocate(allDerived(num_derived + L%Count + L%LikelihoodTypeIndices%Count +1))
        if (num_derived>0) allDerived(:num_derived) =  derived
        call move_alloc(allDerived, derived)
        !Add the chi2 for each likelihood
        derived(num_derived+1:num_derived+L%Count) =  Likelihoods(L%Original_order)*2
        !Add the chi2 for the prior
        derived(num_derived+L%Count+1)  = 2*(logLike - sum(Likelihoods(:L%Count))) !prior
        !Add the total chi2 for each likelihood type
        do i=1, L%LikelihoodTypeIndices%Count
            derived(num_derived+ L%Count+i +1) = sum(Likelihoods(L%LikelihoodTypeIndices%Item(i)))*2
        end do
    end if

    end subroutine addLikelihoodDerivedParams

    subroutine OutputDescription(L, fname)
    class(TLikelihoodList) :: L
    character(LEN=*), intent(in) :: fname
    Type(TTextFile) :: F
    integer i,ix
    Class(TDataLikelihood), pointer :: DataLike

    call F%CreateFile(fname//'.likelihoods')
    do ix=1,L%Count
        if (allocated(L%original_order)) then
            i = L%Original_order(ix)
        else
            i = ix
        end if
        DataLike => L%Item(i)
        !first entry is always 1 for now (in future maybe 0 for likelihoods not sampled from)
        call F%Write(Join(char(9),'1', DataLike%LikelihoodType, DataLike%GetTag(), &
            DataLike%Name, DataLike%Version, trimmed = .true.))
    end do
    call F%Close()

    end subroutine OutputDescription

    end module GeneralTypes
