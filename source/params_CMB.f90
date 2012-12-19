    !Parameterization using theta = r_s/D_a instead of H_0, and tau instead of z_re
    !and log(A_s) instead of A_s
    !Less general, but should give better performance
    !
    !The well-determined parameter A_s exp(-2tau) should be found by the covariance matrix
    !parameter 3 is 100*theta, parameter 4 is tau, others same as params_H except A->log(A)
    !Theta is much better constrained than H_0
    !
    !Assumes prior 0.4 < h < 1

    module DefineParameterization
    use cmbtypes

    Type, extends(TParameterization) :: ThetaParameterization
        real(mcp) :: H0_min = 40, H0_max = 100
        real(mcp) :: use_min_zre
    contains  
    procedure :: ParamsToCMBParams
    procedure :: CMBParamsToParams
    procedure :: NonBaseParameterPriors
    procedure :: CalcDerivedParams
    end type ThetaParameterization

    Type(ThetaParameterization), save, target :: CMBParameterization

    contains

    subroutine SetParameterization(Ini, Names)
    Type(TIniFile) :: Ini
    Type(TParamNames) :: Names
    character(LEN=Ini_max_string_len) :: ParamNamesFile = ''

    Parameterization => CMBParameterization

    ParamNamesFile = ReadIniFileName(Ini,'ParamNamesFile', NotFoundFail=.false.)

    CMBParameterization%H0_min = Ini_Read_Double_File(Ini, 'H0_min',CMBParameterization%H0_min)
    CMBParameterization%H0_max = Ini_Read_Double_File(Ini, 'H0_max',CMBParameterization%H0_max)
    CMBParameterization%Use_min_zre = Ini_Read_Double_File(Ini,'use_min_zre',0.d0) 

    if (ParamNamesFile /='') then
        call ParamNames_init(Names, ParamNamesFile)
    else
        if (generic_mcmc) then
            Names%nnames=0
            if (Feedback>0) write (*,*) 'edit SetParamNames in params_CMB.f90 if you want to use named params'
        else
            call ParamNames_init(Names, trim(LocalDir)//'params_CMB.paramnames')
        end if
    end if

    end subroutine SetParameterization


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
    implicit none
    class(ThetaParameterization) :: this
    real(mcp) Params(num_params)
    integer, parameter :: ncache =2
    real(mcp), save :: LastParams(num_theory_params,ncache) = 0._mcp
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
        if (Feedback>1) print *,'Out of range finding H0'
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
        CMB%zre = GetZreFromTau(CMB, CMB%tau)

        LastCMB(cache) = CMB
        LastParams(:,cache) = Params
        cache = mod(cache,ncache)+1
    end if
    end subroutine ParamsToCMBParams

    subroutine CMBParamsToParams(this,CMB, Params)
    use CMB_Cls
    class(ThetaParameterization) :: this
    real(mcp) Params(num_Params)
    Type(CMBParams) CMB

    Params(1) = CMB%ombh2 

    Params(3) = CMBToTheta(CMB)*100
    Params(4) = CMB%tau
    Params(5) = CMB%omk 

    if (neutrino_param_mnu) then
        Params(2) = CMB%omch2
        Params(6) = CMB%omnuh2*93.04
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
    real(mcp) :: P(max_num_params)
    Type(CMBParams) CMB
    real(mcp) r10
    integer num_derived 

    num_derived = 12 +  Theory%numderived

    allocate(Derived%P(num_derived))

    call Parameterization%ParamsToCMBParams(P,CMB)

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
            CMB%omnuh2=Params(6)/93.04
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

    end module DefineParameterization