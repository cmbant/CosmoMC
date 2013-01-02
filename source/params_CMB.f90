    !Default parameterization using theta = r_s/D_a instead of H_0, and tau instead of z_re
    !and log(A_s) instead of A_s
    !Less general, but should give better performance
    !
    !The well-determined parameter A_s exp(-2tau) should be found by the covariance matrix
    !parameter 3 is 100*theta, parameter 4 is tau, others same as params_H except A->log(A)
    !Theta is much better constrained than H_0
    !
    !Also a background-only parameterization, e.g. for use with just supernoave etc

    module DefineParameterization
    use cmbtypes
    implicit none

    real(mcp), parameter :: neutrino_mass_fac= 93.04

    Type, extends(TParameterization) :: ThetaParameterization
        real(mcp) :: H0_min = 40, H0_max = 100
        real(mcp) :: use_min_zre = 0._mcp
    contains
    procedure :: ParamsToCMBParams
    procedure :: CMBParamsToParams
    procedure :: NonBaseParameterPriors
    procedure :: CalcDerivedParams
    procedure :: Initialize => ThetaParameterization_Init
    end type ThetaParameterization

    Type, extends(TParameterization) :: BackgroundParameterization
    contains
    procedure :: ParamsToCMBParams => Background_ParamsToCMBParams
    procedure :: CMBParamsToParams => Background_CMBParamsToParams
    procedure :: CalcDerivedParams => Background_CalcDerivedParams
    procedure :: Initialize => BackgroundParameterization_Init
    end type BackgroundParameterization

    contains

    subroutine SetTheoryParameterization(Ini, Names)
    Type(TIniFile) :: Ini
    Type(TParamNames) :: Names
    character(LEN=Ini_max_string_len) :: paramtxt
    Type(ThetaParameterization), pointer :: CMBParameterization
    Type(BackgroundParameterization), pointer :: BackgroundParam

    paramtxt = Ini_Read_String_Default_File(Ini,'parameterization', 'theta')
    if (paramtxt =='background') then
        allocate(BackgroundParam)
        Parameterization => BackgroundParam
        call BackgroundParam%Initialize(Ini,Names)
    else if (paramtxt=='theta') then
        allocate(CMBParameterization)
        Parameterization => CMBParameterization
        call CMBParameterization%Initialize(Ini,Names)
    else
        call MpiStop('params_CMB: unknown parameterization :'//trim(paramtxt))
    end if

    end subroutine SetTheoryParameterization

    subroutine ThetaParameterization_Init(this, Ini, Names)
    Class(ThetaParameterization) :: this
    Type(TIniFile) :: Ini
    Type(TParamNames) :: Names

    call SetTheoryParameterNumbers(14,5)

    this%H0_min = Ini_Read_Double_File(Ini, 'H0_min',this%H0_min)
    this%H0_max = Ini_Read_Double_File(Ini, 'H0_max',this%H0_max)
    this%Use_min_zre = Ini_Read_Double_File(Ini,'use_min_zre',this%use_min_zre)

    call this%Init(Ini,Names, 'params_CMB.paramnames')

    end subroutine ThetaParameterization_Init

    function NonBaseParameterPriors(this,CMB)
    class(ThetaParameterization) :: this
    Type(CMBParams) :: CMB
    real(mcp):: NonBaseParameterPriors

    NonBaseParameterPriors = logZero
    if (CMB%H0 < this%H0_min .or. CMB%H0 > this%H0_max) return
    if (CMB%zre < this%Use_min_zre) return
    NonBaseParameterPriors = 0

    end function NonBaseParameterPriors

    subroutine ParamsToCMBParams(this, Params, CMB)
    use CMB_Cls
    class(ThetaParameterization) :: this
    real(mcp) Params(:)
    integer, parameter :: ncache =2
    real(mcp), save :: LastParams(max_theory_params,ncache) = 0._mcp
    Type(CMBParams) CMB
    Type(CMBParams), save :: LastCMB(ncache)
    real(mcp) DA
    real(mcp)  D_b,D_t,D_try,try_b,try_t, lasttry
    integer, save :: cache=1
    integer i

    do i=1, ncache
        !want to save two slow positions for some fast-slow methods
        if (all(Params(1:num_hard) == Lastparams(1:num_hard, i))) then
            CMB = LastCMB(i)
            call SetFast(Params,CMB)
            return
        end if
    end do
    DA = Params(3)/100
    try_b = this%H0_min
    call SetForH(Params,CMB,try_b, .true.)
    D_b = CMBToTheta(CMB)
    try_t = this%H0_max
    call SetForH(Params,CMB,try_t, .false.)
    D_t = CMBToTheta(CMB)
    if (DA < D_b .or. DA > D_t) then
        if (Feedback>1) write(*,*) instance, 'Out of range finding H0: ', real(Params(3))
        cmb%H0=0 !Reject it
    else
        lasttry = -1
        do
            call SetForH(Params,CMB,(try_b+try_t)/2, .false.)
            D_try = CMBToTheta(CMB)
            if (D_try < DA) then
                try_b = (try_b+try_t)/2
            else
                try_t = (try_b+try_t)/2
            end if
            if (abs(D_try - lasttry)< 1e-7) exit
            lasttry = D_try
        end do

        !!call InitCAMB(CMB,error)
        if (CMB%tau==0._mcp) then
            CMB%zre=0
        else
            CMB%zre = GetZreFromTau(CMB, CMB%tau)
        end if

        LastCMB(cache) = CMB
        LastParams(1:num_hard,cache) = Params(1:num_hard)
        cache = mod(cache,ncache)+1
    end if
    end subroutine ParamsToCMBParams

    subroutine CMBParamsToParams(this,CMB, Params)
    use CMB_Cls
    class(ThetaParameterization) :: this
    real(mcp) Params(:)
    Type(CMBParams) CMB

    Params(1) = CMB%ombh2

    Params(3) = CMBToTheta(CMB)*100
    Params(4) = CMB%tau
    Params(5) = CMB%omk

    if (neutrino_param_mnu) then
        Params(2) = CMB%omch2
        Params(6) = CMB%omnuh2*neutrino_mass_fac
    else
        Params(2) = CMB%omdmh2
        Params(6) = CMB%nufrac
    end if
    Params(7) = CMB%w
    Params(8) = CMB%wa
    Params(9) = CMB%nnu
    if (bbn_consistency) then
        Params(10) = 0
    else
        Params(10) = CMB%YHe
    endif
    Params(11) = CMB%iso_cdm_correlated
    Params(12) = CMB%zre_delta
    Params(13) = CMB%ALens
    Params(14) = CMB%fdm

    Params(index_initpower:index_initpower+num_initpower-1) =CMB%InitPower(1:num_initpower)
    Params(index_initpower + As_index-1 ) = log(CMB%InitPower(As_index))

    end subroutine CMBParamsToParams


    function CalcDerivedParams(this, P, Theory, derived) result (num_derived)
    class(ThetaParameterization) :: this
    Type(mc_real_pointer) :: derived
    class(TheoryPredictions) :: Theory
    real(mcp) :: P(:)
    Type(CMBParams) CMB
    real(mcp) r10
    integer num_derived

    num_derived = 12 +  Theory%numderived

    allocate(Derived%P(num_derived))

    call this%ParamsToCMBParams(P,CMB)

    if (lmax_tensor /= 0 .and. compute_tensors) then
        r10 = Theory%cl_tensor(10,1)/Theory%cl(10,1)
    else
        r10 = 0
    end if

    derived%P(1) = CMB%omv
    derived%P(2) = CMB%omdm+CMB%omb
    derived%P(3) = Theory%Sigma_8
    derived%P(4) = CMB%zre
    derived%P(5) = r10
    derived%P(6) = CMB%H0
    derived%P(7) = Theory%tensor_ratio_02
    derived%P(8) = cl_norm*CMB%InitPower(As_index)*1e9
    derived%P(9) = CMB%omdmh2 + CMB%ombh2
    derived%P(10)= (CMB%omdmh2 + CMB%ombh2)*CMB%h
    derived%P(11)= CMB%Yhe !value actually used, may be set from bbn consistency
    derived%P(12)= derived%P(8)*exp(-2*CMB%tau)  !A e^{-2 tau}

    derived%P(13:num_derived) = Theory%derived_parameters(1: Theory%numderived)

    end function CalcDerivedParams

    subroutine SetFast(Params,CMB)
    real(mcp) Params(num_Params)
    Type(CMBParams) CMB

    CMB%InitPower(1:num_initpower) = Params(index_initpower:index_initpower+num_initpower-1)
    CMB%InitPower(As_index) = exp(CMB%InitPower(As_index))

    end subroutine SetFast

    subroutine SetForH(Params,CMB,H0, firsttime)
    use CMB_Cls
    use bbn
    real(mcp) Params(num_Params)
    logical, intent(in) :: firsttime
    Type(CMBParams) CMB
    real(mcp) h2,H0

    CMB%H0=H0
    if (firsttime) then
        CMB%reserved = 0
        CMB%ombh2 = Params(1)
        CMB%tau = params(4) !tau, set zre later
        CMB%Omk = Params(5)
        if (neutrino_param_mnu) then
            !Params(6) is now mnu, params(2) is omch2
            CMB%omnuh2=Params(6)/neutrino_mass_fac
            if (CMB%omnuh2 > 0 .and. (Params(9) < 3 .or. Params(9)>3.1)) &
            call MpiStop('params_CMB: change for non-standard nnu with massive nu')
            CMB%omch2 = Params(2)
            CMB%omdmh2 = CMB%omch2+ CMB%omnuh2
            CMB%nufrac=CMB%omnuh2/CMB%omdmh2
        else
            CMB%omdmh2 = Params(2)
            CMB%nufrac=Params(6)
            CMB%omnuh2 = CMB%omdmh2*CMB%nufrac
            CMB%omch2 = CMB%omdmh2 - CMB%omnuh2
        end if
        CMB%w = Params(7)
        CMB%wa = Params(8)
        CMB%nnu = Params(9) !3.046

        if (bbn_consistency) then
            CMB%YHe = yp_bbn(CMB%ombh2,CMB%nnu  - 3.046)
        else
            !e.g. set from free parameter..
            CMB%YHe  =Params(10)
        end if

        CMB%iso_cdm_correlated =  Params(11)
        CMB%zre_delta = Params(12)
        CMB%ALens = Params(13)
        CMB%fdm = Params(14)
        call SetFast(Params,CMB)
    end if

    CMB%h = CMB%H0/100
    h2 = CMB%h**2
    CMB%omb = CMB%ombh2/h2
    CMB%omc = CMB%omch2/h2
    CMB%omnu = CMB%omnuh2/h2
    CMB%omdm = CMB%omdmh2/h2
    CMB%omv = 1- CMB%omk - CMB%omb - CMB%omdm

    end subroutine SetForH

    !!! Simple parameterization for background data, e.g. Supernovae only (no thermal history)
    subroutine BackgroundParameterization_Init(this, Ini, Names)
    Class(BackgroundParameterization) :: this
    Type(TIniFile) :: Ini
    Type(TParamNames) :: Names

    call SetTheoryParameterNumbers(7,0)
    this%late_time_only = .true.
    call this%Init(Ini,Names, 'params_background.paramnames')

    end subroutine BackgroundParameterization_Init

    subroutine Background_ParamsToCMBParams(this, Params, CMB)
    class(BackgroundParameterization) :: this
    real(mcp) Params(:)
    Type(CMBParams) CMB
    real(mcp) omegam, h2

    omegam = Params(1)
    CMB%H0 = Params(2)
    CMB%omk = Params(3)
    CMB%omnuh2=Params(4)/neutrino_mass_fac
    CMB%w =    Params(5)
    CMB%wa =    Params(6)
    CMB%nnu =    Params(7)

    CMB%h=CMB%H0/100
    h2 = CMB%h**2
    CMB%Yhe=0.24
    CMB%omnu = CMB%omnuh2/h2
    CMB%omb= omegam - CMB%omnu
    CMB%ombh2 = CMB%omb*h2
    CMB%omc=0
    CMB%omch2 = CMB%omc*h2
    CMB%zre=0
    CMB%tau=0
    CMB%omdmh2 = CMB%omch2+ CMB%omnuh2
    CMB%omdm = CMB%omdmh2/h2
    CMB%omv = 1- CMB%omk - CMB%omb - CMB%omdm
    CMB%nufrac=CMB%omnuh2/CMB%omdmh2
    CMB%reserved=0
    CMB%fdm=0
    CMB%iso_cdm_correlated=0
    CMB%Alens=1

    end subroutine Background_ParamsToCMBParams

    subroutine Background_CMBParamsToParams(this,CMB, Params)
    use CMB_Cls
    class(BackgroundParameterization) :: this
    real(mcp) Params(:)
    Type(CMBParams) CMB

    Params(1) = CMB%omb+CMB%omc+CMB%omnu
    Params(2) = CMB%H0
    Params(3) = CMB%omk
    Params(4) = CMB%omnuh2*neutrino_mass_fac
    Params(5) = CMB%w
    Params(6) = CMB%wa
    Params(7) = CMB%nnu
    end subroutine Background_CMBParamsToParams

    function Background_CalcDerivedParams(this, P, Theory, derived) result (num_derived)
    class(BackgroundParameterization) :: this
    Type(mc_real_pointer) :: derived
    class(TheoryPredictions) :: Theory
    real(mcp) :: P(:)
    Type(CMBParams) CMB
    integer num_derived

    num_derived = 1

    allocate(Derived%P(num_derived))

    call this%ParamsToCMBParams(P,CMB)

    derived%P(1) = CMB%omv

    end function Background_CalcDerivedParams


    end module DefineParameterization