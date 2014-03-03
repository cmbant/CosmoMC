    module CosmologyConfig
    use GeneralTypes
    use CalcLike
    use cmbtypes
    use CosmoTheory
    use Calculator_Cosmology
    use CosmologyParameterizations
    use CalcLike_Cosmology
    use SampleCollector
    use GeneralSetup
    use Likelihood_Cosmology
    implicit none
    private

    Type, extends(TGeneralConfig) :: TCosmologyConfig
    contains
    procedure :: SetParameterizationName => TCosmologyConfig_SetParameterizationName
    procedure :: NewTheory => TCosmologyConfig_NewTheory
    procedure :: InitForLikelihoods => TCosmologyConfig_InitForLikelihoods
    procedure :: ReadParams => TCosmologyConfig_ReadParams
    end Type

    type, extends(TSetup) :: TCosmologySetup
    contains
    procedure :: Init => TCosmologySetup_Init
    end type

    public TCosmologyConfig, TCosmologySetup
    contains

    subroutine TCosmologyConfig_ReadParams(this, Ini)
    use Calculator_CAMB
    class(TCosmologyConfig) :: this
    class(TSettingIni) :: Ini
    character(LEN=:), allocatable :: CalcName

    CalcName = Ini%Read_String_Default('cosmology_calculator', 'CAMB')
    if (calcName=='CAMB') then
        allocate(CAMB_Calculator::this%Calculator)
    else
        call MpiStop('only CAMB currently supported')
    end if
    call this%Calculator%InitWithParams(Ini,this)

    compute_tensors = Ini%Read_Logical('compute_tensors',.false.)

    if (num_cls==3 .and. compute_tensors) write (*,*) 'WARNING: computing tensors with num_cls=3 (BB=0)'
    call Ini%Read('CMB_lensing',CMB_lensing)
    call Ini%Read('use_lensing_potential',use_lensing_potential)
    call Ini%Read('use_nonlinear_lensing',use_nonlinear_lensing)

    call Ini%Read('pivot_k',pivot_k)
    call Ini%read('inflation_consistency',inflation_consistency)
    call Ini%Read('bbn_consistency',bbn_consistency)
    num_massive_neutrinos = Ini%read_int('num_massive_neutrinos',-1)

    if (CMB_lensing) num_clsS = num_cls   !Also scalar B in this case
    if (use_lensing_potential .and. num_cls_ext ==0) &
    call MpiStop('num_cls_ext should be > 0 to use_lensing_potential')
    if (use_lensing_potential .and. .not. CMB_lensing) &
    call MpiStop('use_lensing_potential must have CMB_lensing=T')
    lmax_computed_cl = Ini%Read_Int('lmax_computed_cl',lmax)

    end subroutine TCosmologyConfig_ReadParams

    function TCosmologyConfig_SetParameterizationName(this, nametag, Ini, Names) result(OK)
    class(TCosmologyConfig) :: this
    character(LEN=*), intent(in) :: nametag
    class(TSettingIni) :: Ini
    class(TParamNames) :: Names
    logical OK
    Type(ThetaParameterization), pointer :: CMBParameterization
    Type(BackgroundParameterization), pointer :: BackgroundParam

    OK = .true.
    if (nametag =='background') then
        allocate(BackgroundParam)
        this%Parameterization => BackgroundParam
        call BackgroundParam%InitWithSetNames(Ini,Names,this)
    else if (nametag=='theta') then
        allocate(CMBParameterization)
        this%Parameterization => CMBParameterization
        call CMBParameterization%InitWithSetNames(Ini,Names,this)
    else
        OK =  this%TGeneralConfig%SetParameterizationName(nametag,Ini,Names)
    end if

    end function TCosmologyConfig_SetParameterizationName


    subroutine TCosmologyConfig_NewTheory(this, Theory)
    class(TCosmologyConfig) :: this
    class(TTheoryPredictions), allocatable :: Theory

    allocate(TCosmoTheoryPredictions::Theory)
    call Theory%Init(this)

    end subroutine TCosmologyConfig_NewTheory

    subroutine TCosmologyConfig_InitForLikelihoods(this)
    class(TCosmologyConfig) :: this

    if(use_LSS) call Initialize_PKSettings()
    call this%TGeneralConfig%InitForLikelihoods()

    if (Feedback > 0 .and. MPIRank==0) then
        write (*,*) 'Computing tensors:', compute_tensors
        write (*,*) 'Doing CMB lensing:',CMB_lensing
        write (*,*) 'Doing non-linear Pk:', use_nonlinear

        write(*,'(" lmax              = ",1I4)') lmax
        write(*,'(" lmax_computed_cl  = ",1I4)') lmax_computed_cl

        if (compute_tensors) write(*,'(" lmax_tensor    = ",1I4)') lmax_tensor
        write(*,'(" Number of C_ls = ",1I4)') num_cls
    end if

    end subroutine TCosmologyConfig_InitForLikelihoods


    !!!TCosmologySetup

    subroutine TCosmologySetup_Init(this)
    class(TCosmologySetup) :: this

    allocate(TCosmologyConfig::this%Config)
    allocate(TCosmoLikeCalculator::this%LikeCalculator)
    call this%TSetup%Init()

    end subroutine TCosmologySetup_Init

    end module CosmologyConfig
