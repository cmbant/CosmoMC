    module CosmologyConfig
    use GeneralTypes
    use CalcLike
    use cmbtypes
    use Calculator_Cosmology
    use CosmologyParameterizations
    use CalcLike_Cosmology
    use SampleCollector
    use GeneralSetup
    implicit none
    private

    Type, extends(TGeneralConfig) :: TCosmologyConfig
    contains
    procedure :: SetParameterizationName => TCosmologyConfig_SetParameterizationName
    procedure :: NewTheory => TCosmologyConfig_NewTheory
    procedure :: InitForLikelihoods => TCosmologyConfig_InitForLikelihoods
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
    CMB_lensing = Ini%Read_Logical('CMB_lensing',CMB_lensing)
    use_lensing_potential = Ini%Read_logical('use_lensing_potential',use_lensing_potential)
    use_nonlinear_lensing = Ini%Read_logical('use_nonlinear_lensing',use_nonlinear_lensing)

    pivot_k = Ini%Read_Real('pivot_k',0.05)
    inflation_consistency = Ini%read_Logical('inflation_consistency',.false.)
    bbn_consistency = Ini%Read_Logical('bbn_consistency',.true.)
    num_massive_neutrinos = Ini%read_int('num_massive_neutrinos',-1)

    if (CMB_lensing) num_clsS = num_cls   !Also scalar B in this case
    if (use_lensing_potential .and. num_cls_ext ==0) &
    call MpiStop('num_cls_ext should be > 0 to use_lensing_potential')
    if (use_lensing_potential .and. .not. CMB_lensing) &
    call MpiStop('use_lensing_potential must have CMB_lensing=T')
    lmax_computed_cl = Ini%Read_Int('lmax_computed_cl',lmax)

    if (Feedback > 0 .and. MPIRank==0) then
        write (*,*) 'Computing tensors:', compute_tensors
        write (*,*) 'Doing CMB lensing:',CMB_lensing
        write (*,*) 'Doing non-linear Pk:', use_nonlinear

        write(*,'(" lmax              = ",1I4)') lmax
        write(*,'(" lmax_computed_cl  = ",1I4)') lmax_computed_cl

        if (compute_tensors) write(*,'(" lmax_tensor    = ",1I4)') lmax_tensor
        write(*,'(" Number of C_ls = ",1I4)') num_cls
    end if

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

    function TCosmologyConfig_NewTheory(this) result(Theory)
    class(TCosmologyConfig) :: this
    class(TTheoryPredictions), pointer :: Theory

    allocate(TCosmoTheoryPredictions::Theory)
    call Theory%Init(this)

    end function TCosmologyConfig_NewTheory

    subroutine TCosmologyConfig_InitForLikelihoods(this)
    class(TCosmologyConfig) :: this

    if(use_LSS) call Initialize_PKSettings()

    end subroutine TCosmologyConfig_InitForLikelihoods


    !!!TCosmologySetup

    subroutine TCosmologySetup_Init(this)
    class(TCosmologySetup) :: this

    allocate(TCosmologyConfig::this%Config)
    allocate(TCosmoLikeCalculator::this%LikeCalculator)
    call this%TSetup%Init()

    end subroutine TCosmologySetup_Init

    end module CosmologyConfig
