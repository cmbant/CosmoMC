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

    Type, extends(CosmologyParameterization) :: ThetaParameterization
        real(mcp) :: H0_min = 40, H0_max = 100
        real(mcp) :: H0_prior_mean = 0._mcp, H0_prior_std = 0._mcp
        real(mcp) :: use_min_zre = 0._mcp
    contains
    procedure :: ParamArrayToTheoryParams => TP_ParamArrayToTheoryParams
    procedure :: NonBaseParameterPriors => TP_NonBaseParameterPriors
    procedure :: CalcDerivedParams => TP_CalcDerivedParams
    procedure :: Initialize => TP_Init
    end type ThetaParameterization

    Type, extends(CosmologyParameterization) :: BackgroundParameterization
    contains
    procedure :: ParamArrayToTheoryParams => BK_ParamArrayToTheoryParams
    procedure :: CalcDerivedParams => BK_CalcDerivedParams
    procedure :: Initialize => BK_Init
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

    subroutine TP_Init(this, Ini, Names)
    Class(ThetaParameterization) :: this
    Type(TIniFile) :: Ini
    Type(TParamNames) :: Names
    character(LEN=Ini_max_string_len) prior
    call SetTheoryParameterNumbers(14,6)

    this%H0_min = Ini_Read_Double_File(Ini, 'H0_min',this%H0_min)
    this%H0_max = Ini_Read_Double_File(Ini, 'H0_max',this%H0_max)
    this%Use_min_zre = Ini_Read_Double_File(Ini,'use_min_zre',this%use_min_zre)
    prior = Ini_Read_String_File(Ini, 'H0_prior')
    if (prior/='') then
        read(prior,*) this%H0_prior_mean, this%H0_prior_std
    end if

    call this%Init(Ini,Names, 'params_CMB.paramnames')

    end subroutine TP_Init

    function TP_NonBaseParameterPriors(this,CMB)
    class(ThetaParameterization) :: this
    class(TTheoryParams) :: CMB
    real(mcp):: TP_NonBaseParameterPriors

    select type (CMB)
    class is (CMBParams)
        TP_NonBaseParameterPriors = logZero
        if (CMB%H0 < this%H0_min .or. CMB%H0 > this%H0_max) return
        if (CMB%zre < this%Use_min_zre) return
        TP_NonBaseParameterPriors = 0
        if (this%H0_prior_mean/=0._mcp) then
            TP_NonBaseParameterPriors = ((CMB%H0 - this%H0_prior_mean)/this%H0_prior_std)**2/2
        end if
    end select
    end function TP_NonBaseParameterPriors

    subroutine TP_ParamArrayToTheoryParams(this, Params, CMB)
    use CMB_Cls
    class(ThetaParameterization) :: this
    real(mcp) Params(:)
    integer, parameter :: ncache =2
    Class(TTheoryParams), target :: CMB
    Type(CMBParams), save :: LastCMB(ncache)
    real(mcp) DA
    real(mcp)  D_b,D_t,D_try,try_b,try_t, lasttry
    integer, save :: cache=1
    integer i
    Type(CMBParams), pointer :: CP

    select type (CMB)
    class is (CMBParams)
        do i=1, ncache
            !want to save two slow positions for some fast-slow methods
            if (all(Params(1:num_hard) == LastCMB(i)%BaseParams(1:num_hard))) then
                CP => CMB !needed to make next line work for some odd reason CMB=LastCMB(i) does not work
                CP = LastCMB(i)
                call this%CosmologyParameterization%ParamArrayToTheoryParams(Params, CMB)
                call SetFast(Params,CMB)
                return
            end if
        end do
        call this%CosmologyParameterization%ParamArrayToTheoryParams(Params, CMB)

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
            cache = mod(cache,ncache)+1
        end if
    end select
    end subroutine TP_ParamArrayToTheoryParams

    function TP_CalcDerivedParams(this, P, Theory, derived) result (num_derived)
    class(ThetaParameterization) :: this
    Type(mc_real_pointer) :: derived
    class(TTheoryPredictions) :: Theory
    real(mcp) :: P(:)
    Type(CMBParams) CMB
    real(mcp) r10
    integer num_derived

    select type (Theory)
    class is (TheoryPredictions)
        num_derived = 12 +  Theory%numderived
        allocate(Derived%P(num_derived))

        call this%ParamArrayToTheoryParams(P,CMB)

        derived%P(1) = CMB%omv
        derived%P(2) = CMB%omdm+CMB%omb
        derived%P(3) = Theory%Sigma_8
        derived%P(4) = CMB%zre
        derived%P(5) = Theory%tensor_ratio_r10
        derived%P(6) = CMB%H0
        derived%P(7) = Theory%tensor_ratio_02
        derived%P(8) = cl_norm*CMB%InitPower(As_index)*1e9
        derived%P(9) = CMB%omdmh2 + CMB%ombh2
        derived%P(10)= (CMB%omdmh2 + CMB%ombh2)*CMB%h
        derived%P(11)= CMB%Yhe !value actually used, may be set from bbn consistency
        derived%P(12)= derived%P(8)*exp(-2*CMB%tau)  !A e^{-2 tau}

        derived%P(13:num_derived) = Theory%derived_parameters(1: Theory%numderived)
    end select

    end function TP_CalcDerivedParams

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
    subroutine BK_Init(this, Ini, Names)
    Class(BackgroundParameterization) :: this
    Type(TIniFile) :: Ini
    Type(TParamNames) :: Names

    call SetTheoryParameterNumbers(7,0)
    this%late_time_only = .true.
    call this%Init(Ini,Names, 'params_background.paramnames')

    end subroutine BK_Init

    subroutine BK_ParamArrayToTheoryParams(this, Params, CMB)
    class(BackgroundParameterization) :: this
    real(mcp) Params(:)
    class(TTheoryParams), target :: CMB
    real(mcp) omegam, h2

    select type (CMB)
    class is (CMBParams)
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
    end select
    end subroutine BK_ParamArrayToTheoryParams


    function BK_CalcDerivedParams(this, P, Theory, derived) result (num_derived)
    class(BackgroundParameterization) :: this
    Type(mc_real_pointer) :: derived
    class(TTheoryPredictions) :: Theory
    real(mcp) :: P(:)
    Type(CMBParams) CMB
    integer num_derived

    num_derived = 1

    allocate(Derived%P(num_derived))

    call this%ParamArrayToTheoryParams(P,CMB)

    derived%P(1) = CMB%omv

    end function BK_CalcDerivedParams


    end module DefineParameterization