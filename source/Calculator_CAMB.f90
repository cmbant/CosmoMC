    !Use CAMB
    module Calculator_CAMB
    use cmbtypes
    use CosmoTheory
    use CAMB, only : CAMB_GetResults, CAMB_GetAge, CAMBParams, CAMB_SetDefParams, &
    AccuracyBoost,  Cl_scalar, Cl_tensor, Cl_lensed, outNone, w_lam, wa_ppf,&
    CAMBParams_Set, MT, CAMBdata, NonLinear_Pk, Nonlinear_lens, Reionization_GetOptDepth, CAMB_GetZreFromTau, &
    CAMB_GetTransfers,CAMB_FreeCAMBdata,CAMB_InitCAMBdata, CAMB_TransfersToPowers, Transfer_SetForNonlinearLensing, &
    initial_adiabatic,initial_vector,initial_iso_baryon,initial_iso_CDM, initial_iso_neutrino, initial_iso_neutrino_vel, &
    HighAccuracyDefault, highL_unlensed_cl_template, ThermoDerivedParams, nthermo_derived, BackgroundOutputs, &
    Transfer_SortAndIndexRedshifts, & !JD added for nonlinear lensing of CMB + MPK compatibility
    Recombination_Name, reionization_name, power_name, threadnum, version
    use Errors !CAMB
    use settings
    use likelihood
    use Calculator_Cosmology
    use GeneralTypes
    implicit none
    private

    Type, extends(TTheoryIntermediateCache) :: CAMBTransferCache
        Type (CAMBdata) :: Transfers
    contains
    procedure :: Clear => CAMBTransferCache_Clear
    end Type CAMBTransferCache

    Type, extends(TCosmologyCalculator) :: CAMB_Calculator
        integer :: ncalls = 0
        integer :: nerrors = 0
        logical :: CAMB_timing = .false.
        type(CAMBParams)  CAMBP
        character(LEN=:), allocatable :: highL_theory_cl_template_file
        real(mcp), allocatable :: highL_lensedCL_template(:,:)
    contains
    !New
    procedure :: CMBToCAMB => CAMBCalc_CMBToCAMB
    procedure :: SetDerived => CAMBCalc_SetDerived
    procedure :: SetPowersFromCAMB => CAMBCalc_SetPowersFromCAMB
    procedure :: InitCAMB => CAMBCalc_InitCAMB
    procedure :: InitCAMBParams => CAMBCalc_InitCAMBParams
    procedure :: SetCAMBInitPower => CAMBCalc_SetCAMBInitPower
    procedure :: SetPkFromCAMB => CAMBCalc_SetPkFromCAMB
    procedure :: TransfersOrPowers => CAMBCalc_TransfersOrPowers
    procedure :: GetNLandRatios => CAMBCalc_GetNLandRatios
    !Overridden inherited
    procedure :: ReadParams => CAMBCalc_ReadParams
    procedure :: InitForLikelihoods => CAMBCalc_InitForLikelihoods
    procedure :: BAO_D_v => CAMBCalc_BAO_D_v
    procedure :: AngularDiameterDistance => CAMBCalc_AngularDiameterDistance
    procedure :: Hofz => CAMBCalc_Hofz
    procedure :: CMBToTheta => CAMBCalc_CMBToTheta
    procedure :: GetNewPowerData => CAMBCalc_GetNewPowerData
    procedure :: GetNewTransferData => CAMBCalc_GetNewTransferData
    procedure :: GetCosmoTheoryForImportance => CAMBCalc_GetTheoryForImportance
    procedure :: GetZreFromTau  => CAMBCalc_GetZreFromTau
    procedure :: GetOpticalDepth => CAMBCalc_GetOpticalDepth
    procedure :: SetBackgroundTheoryData => CAMBCalc_SetBackgroundTheoryData
    procedure :: SetParamsForBackground => CAMBCalc_SetParamsForBackground
    procedure :: VersionTraceOutput => CAMBCalc_VersionTraceOutput
    procedure, private :: LoadFiducialHighLTemplate
    end type CAMB_Calculator


    integer, parameter :: ScalClOrder(5) = (/1,3,2,4,5/), TensClOrder(4) = (/1,4,2,3/)
    !Mapping of CAMB CL array ordering to TT , TE, EE, BB, phi, phiT


    public CAMB_Calculator
    contains

    subroutine CAMBCalc_CMBToCAMB(this,CMB,P)
    use LambdaGeneral
    use CAMBmain, only : ALens
    use constants, only : default_nnu
    class(CAMB_Calculator) :: this
    class(CMBParams) CMB
    type(CAMBParams)  P
    real(dl) neff_massive_standard

    P = this%CAMBP
    P%omegab = CMB%omb
    P%omegan = CMB%omnu
    P%omegac = CMB%omc
    P%omegav = CMB%omv
    P%H0 = CMB%H0
    P%Reion%redshift= CMB%zre
    P%Reion%delta_redshift = CMB%zre_delta
    w_lam = CMB%w
    wa_ppf = CMB%wa
    ALens = CMB%ALens
    P%InitialConditionVector(initial_iso_CDM) = &
    sign(sqrt(abs(CMB%iso_cdm_correlated) /(1-abs(CMB%iso_cdm_correlated))),CMB%iso_cdm_correlated)
    P%Num_Nu_Massive = 0
    P%Nu_mass_numbers = 0
    P%Num_Nu_Massless = CMB%nnu
    if (CMB%omnuh2>0) then
        P%Nu_mass_eigenstates=0
        if (CMB%omnuh2>CMB%omnuh2_sterile) then
            neff_massive_standard = num_massive_neutrinos*default_nnu/3
            P%Num_Nu_Massive = num_massive_neutrinos
            P%Nu_mass_eigenstates=P%Nu_mass_eigenstates+1
            if (CMB%nnu > neff_massive_standard) then
                P%Num_Nu_Massless = CMB%nnu - neff_massive_standard
            else
                P%Num_Nu_Massless = 0
                neff_massive_standard=CMB%nnu
            end if
            P%Nu_mass_numbers(P%Nu_mass_eigenstates) = num_massive_neutrinos
            P%Nu_mass_degeneracies(P%Nu_mass_eigenstates) = neff_massive_standard
            P%Nu_mass_fractions(P%Nu_mass_eigenstates) = (CMB%omnuh2-CMB%omnuh2_sterile)/CMB%omnuh2
        else
            neff_massive_standard=0
        end if
        if (CMB%omnuh2_sterile>0) then
            if (CMB%nnu<default_nnu) call MpiStop('nnu < 3.046 with massive sterile')
            P%Num_Nu_Massless = default_nnu - neff_massive_standard
            P%Num_Nu_Massive=P%Num_Nu_Massive+1
            P%Nu_mass_eigenstates=P%Nu_mass_eigenstates+1
            P%Nu_mass_numbers(P%Nu_mass_eigenstates) = 1
            P%Nu_mass_degeneracies(P%Nu_mass_eigenstates) = max(1d-6,CMB%nnu - default_nnu)
            P%Nu_mass_fractions(P%Nu_mass_eigenstates) = CMB%omnuh2_sterile/CMB%omnuh2
        end if
    end if

    P%YHe = CMB%YHe
#ifdef COSMOREC
    if (P%Recomb%fdm/=0._mcp) P%Recomb%runmode = 3
    P%Recomb%fdm = CMB%fdm * 1e-23_mcp
#else
    if (CMB%fdm/=0._mcp) call MpiStop('Compile with CosmoRec to use fdm')
#endif
    call this%SetCAMBInitPower(P,CMB,1)

    end subroutine CAMBCalc_CMBToCAMB

    subroutine CAMBCalc_SetDerived(this,Theory)
    class(CAMB_Calculator) :: this
    Class(TCosmoTheoryPredictions) Theory
    integer noutputs, i

    noutputs = size(BackgroundOutputs%z_outputs)
    Theory%numderived = nthermo_derived + noutputs*3
    if (Theory%numderived > max_derived_parameters) &
    call MpiStop('numderived > max_derived_parameters: increase in cmbtypes.f90')
    Theory%derived_parameters(1:nthermo_derived) = ThermoDerivedParams(1:nthermo_derived)
    do i=1, noutputs
        Theory%derived_parameters(nthermo_derived+(i-1)*3+1) = BackgroundOutputs%rs_by_D_v(i)
        Theory%derived_parameters(nthermo_derived+(i-1)*3+2) = BackgroundOutputs%H(i)
        Theory%derived_parameters(nthermo_derived+(i-1)*3+3) = BackgroundOutputs%DA(i)
    end do
    end subroutine CAMBCalc_SetDerived

    subroutine CAMBCalc_SetParamsForBackground(this,CMB)
    use Camb, only: CAMBParams_Set
    class(CAMB_Calculator) :: this
    class(CMBParams) CMB
    Type(CAMBParams)  P

    !set background dparameters, but don't calculate thermal history
    call this%CMBToCAMB(CMB, P)
    call CAMBParams_Set(P)

    end subroutine CAMBCalc_SetParamsForBackground

    subroutine CAMBCalc_SetBackgroundTheoryData(this, CMB,Theory,error)
    use cambmain, only: initvars
    class(CAMB_Calculator) :: this
    class(CMBParams) CMB
    Class(TCosmoTheoryPredictions) Theory
    integer error

    call InitVars !calculate thermal history, e.g. z_drag etc.
    if (global_error_flag/=0) then
        error=global_error_flag
        return
    end if
    call this%SetDerived(Theory)

    end subroutine CAMBCalc_SetBackgroundTheoryData

    subroutine CAMBCalc_GetNewTransferData(this, CMB,Info,Theory,error)
    use InitialPower
    class(CAMB_Calculator) :: this
    class(CMBParams) CMB
    class(TTheoryIntermediateCache), pointer :: Info
    class(TCosmoTheoryPredictions) :: Theory
    integer error
    type(CAMBParams)  P
    real(mcp) time

    allocate(CAMBTransferCache::Info)

    select type (Info)
    class is (CAMBTransferCache)
        call CAMB_InitCAMBdata(Info%Transfers)
        call this%CMBToCAMB(CMB, P)

        if (Feedback > 1) write (*,*) 'Calling CAMB'
        Threadnum =num_threads

        if (this%CAMB_timing) time = TimerTime()
        call CAMB_GetTransfers(P, Info%Transfers, error)
        if (this%CAMB_timing) call Timer('GetTransfers', time)
        class default
        call MpiStop('CAMB_GetNewTransferData with wrong TTheoryIntermediateCache type')
    end select
    if (error==0) then
        call this%SetDerived(Theory)
    else
        if (stop_on_error) call MpiStop('CAMB error '//trim(global_error_message))
        if (Feedback > 0) write(*,*) 'CAMB returned error '//trim(global_error_message)
        this%nerrors=this%nerrors+1
    end if
    this%ncalls=this%ncalls+1
    if (mod(this%ncalls,100)==0 .and. LogFile%Opened()) then
        write(LogFile%unit,'("CAMB called ",I0," times; ",I0," errors")')this%ncalls, this%nerrors
    end if
    if (Feedback > 1) write (*,*) 'CAMB done'

    end subroutine CAMBCalc_GetNewTransferData

    subroutine CAMBCalc_GetNewPowerData(this, CMB, Info, Theory, error)
    class(CAMB_Calculator) :: this
    class(CMBParams) :: CMB
    class(TTheoryIntermediateCache), pointer :: Info
    class(TCosmoTheoryPredictions) :: Theory
    integer error

    select type (Info)
    class is (CAMBTransferCache)
        call this%SetCAMBInitPower(Info%Transfers%Params,CMB,1)
        call CAMB_TransfersToPowers(Info%Transfers)
        !this sets slow CAMB params correctly from value stored in Transfers
        if (global_error_flag/=0) then
            error=global_error_flag
            return
        end if
        !JD 08/13 added so we dont have to fill Cls unless using CMB
        if(use_CMB)then
            call this%SetPowersFromCAMB(CMB,Theory)

            if (any(Theory%cl(:,1) < 0 )) then
                error = 1
                call MpiStop('CMB_cls_simple: negative C_l (could edit to silent error here)')
                return
            end if
        else
            Theory%cl(:,:)=0
        end if

        !redshifts are in increasing order, so last index is redshift zero
        if (Use_LSS .or. get_sigma8) then
            Theory%sigma_8 = Info%Transfers%MTrans%sigma_8(size(Info%Transfers%MTrans%sigma_8,1),1)
        else
            Theory%sigma_8 = 0
        end if

        if (Use_LSS) then
            call this%SetPkFromCAMB(Info%Transfers%MTrans,Theory)
        end if
    end select

    end subroutine CAMBCalc_GetNewPowerData

    subroutine CAMBCalc_GetTheoryForImportance(this, CMB, Theory, error)
    use ModelParams, only : ThreadNum
    class(CAMB_Calculator) :: this
    class(CMBParams) CMB
    class(TCosmoTheoryPredictions) Theory
    integer error
    logical :: DoCls, DoPk
    type(CAMBParams)  P

    error = 0
    DoCls = this%ImportanceOptions%redo_cls
    DoPk = this%ImportanceOptions%redo_pk

    if (DoCls .or. DoPk) then
        Threadnum =num_threads
        call this%CMBToCAMB(CMB, P)
        P%OnlyTransfers = .false.

        if (DoPk) then
            P%WantTransfer = .true.
            if (.not. DoCls) then
                P%WantTensors = .false.
            end if
        end if
        if (DoCls) then
            !Assume we just want Cls to higher l
            P%WantTensors = compute_tensors
            !!!not OK for non-linear lensing        if (.not. DoPk) P%WantTransfer = .false.
        end if

        if (this%CAMB_timing) call Timer()
        if (Feedback > 1) write (*,*) 'Calling CAMB'
        call CAMB_GetResults(P)
        if (Feedback > 1) write (*,*) 'CAMB Done'
        if (this%CAMB_timing) call Timer('CAMB_GetResults')

        error = global_error_flag !using error optional parameter gives seg faults on SGI
    else
        call this%GetNewBackgroundData(CMB,Theory,error)
    end if
    if (error==0) then
        if (DoCls) call this%SetPowersFromCAMB(CMB,Theory)
        if (DoPK) then
            Theory%sigma_8 = MT%sigma_8(size(MT%sigma_8,1),1)
            call this%SetPkFromCAMB(MT,Theory)
        end if
        call this%SetDerived(Theory)
    end if
    end subroutine CAMBCalc_GetTheoryForImportance

    subroutine CAMBCalc_SetPowersFromCAMB(this,CMB,Theory)
    use constants
    use InitialPower
    class(CAMB_Calculator) :: this
    class(CMBParams) :: CMB
    class(TCosmoTheoryPredictions) Theory
    real(mcp), parameter :: cons =  (COBE_CMBTemp*1e6)**2*2*pi
    real(mcp) nm
    integer l
    real(mcp) highL_norm, lens_recon_scale

    Theory%cl=0
    lens_recon_scale = CMB%InitPower(Aphiphi_index)

    do l = 2, lmax_computed_cl
        nm = cons/(l*(l+1))
        if (CMB_Lensing) then
            Theory%cl(l,1:num_clsS) =  nm*Cl_lensed(l,1, TensClOrder(1:num_clsS))
        else
            Theory%cl(l,1:num_clsS) =  nm*Cl_scalar(l,1, scalClOrder(1:num_clsS))
        end if

        if (num_cls>num_clsS) Theory%cl(l,num_clsS+1:num_cls) = 0

        if (compute_tensors .and. l<=lmax_tensor) then
            Theory%cl(l,1:num_cls) =  Theory%cl(l,1:num_cls) + nm*Cl_tensor(l,1, TensClOrder(1:num_cls))
        end if

        if (num_cls_ext > 0) then
            !CMB lensing potential
            !in camb Cphi is l^4 C_l, we want [l(l+1)]^2Cphi/2pi
            if (.not. CMB_lensing) call MpiStop('Must have lensing on to use lensing potential')
            Theory%cl(l,num_clsS+1) =  Cl_scalar(l,1, scalClOrder(4))*(real(l+1)**2/l**2)/twopi * lens_recon_scale
            if (num_cls_ext>1) then
                !lensing-temp
                if (num_cls_ext>1) call MpiStop('SetTheoryFromCAMB: check defs for num_cls_ext>1')
                Theory%cl(l,num_clsS+2) =   Cl_scalar(l,1, scalClOrder(5))/real(l)**3 * sqrt(lens_recon_scale)
            end if
        end if
    end do

    if (compute_tensors) then
        Theory%tensor_ratio_02 = TensorPower(0.002d0,1)/ScalarPower(0.002d0,1)
        Theory%tensor_ratio_r10 = Cl_tensor(10, 1, 1)/Cl_scalar(10,1, 1)
    else
        Theory%tensor_ratio_02 = 0
        Theory%tensor_ratio_r10 = 0
    end if

    if (lmax_computed_cl/=lmax) then
        !use template for very high L theory tail, scaled not to be discontinuous
        highL_norm = Theory%cl(lmax_computed_cl,1)/this%highL_lensedCL_template(lmax_computed_cl,1)
        do l = lmax_computed_cl+1, lmax
            Theory%cl(l,1:num_ClsS) =  highL_norm*this%highL_lensedCL_template(l,1:num_clsS)
        end do
    end if

    end subroutine CAMBCalc_SetPowersFromCAMB

    subroutine CAMBCalc_SetPkFromCAMB(this,M,Theory)
    use Transfer
    use camb, only : CP
    class(CAMB_Calculator) :: this
    class(TCosmoTheoryPredictions) Theory
    Type(MatterTransferData) M
    real(mcp), allocatable :: k(:), z(:), PK(:,:)
    integer zix,nz,nk
    !For use with for example WL PK's
    real(mcp), allocatable :: NL_Ratios(:,:)

    !Free theory arrays as they may resize between samples
    call Theory%FreePK()
    allocate(Theory%MPK)

    nk=M%num_q_trans
    nz=CP%Transfer%PK_num_redshifts
    allocate(PK(nk,nz))
    allocate(k(nk))
    allocate(z(nz))

    k = log(M%TransferData(Transfer_kh,:,1))
    do zix=1,nz
        z(zix) = CP%Transfer%PK_redshifts(nz-zix+1)
    end do

    call this%TransfersOrPowers(M,PK,transfer_power_var,transfer_power_var)
    PK = Log(PK)
    call Theory%MPK%Init(k,z,PK)

    if(use_nonlinear)then
        call this%GetNLandRatios(M,Theory,NL_Ratios)
    end if

    end subroutine CAMBCalc_SetPkFromCAMB

    subroutine CAMBCalc_TransfersOrPowers(this,M,PK,t1,t2)
    use Transfer
    use camb, only : CP, ScalarPower
    class(CAMB_Calculator) :: this
    Type(MatterTransferData) :: M
    real(mcp), intent(inout):: PK(:,:)
    integer, intent(in) :: t1
    integer, optional, intent(in) :: t2
    real(mcp), allocatable :: temp(:,:)
    real(mcp) h, k
    integer nz, nk, zix, ik

    nk=size(PK,1)
    nz=size(PK,2)

    allocate(temp(nk,nz))

    h = CP%H0/100

    if(present(t2)) then
        do ik=1,nk
            k = M%TransferData(Transfer_kh,ik,1)*h
            temp(ik,:) = M%TransferData(t1,ik,:)*&
            M%TransferData(t2,ik,:)*k*pi*twopi*h**3*scalarPower(k,1)
        end do
    else
        temp = M%TransferData(t1,:,:)
    end if

    do zix=1,nz
        PK(:,zix) = temp(:,nz-zix+1)
    end do

    end subroutine CAMBCalc_TransfersOrPowers

    subroutine CAMBCalc_GetNLandRatios(this,M,Theory,Ratios)
    use Transfer
    class(CAMB_Calculator) :: this
    class(TCosmoTheoryPredictions) Theory
    Type(MatterTransferData) M
    real(mcp), allocatable, intent(out) :: Ratios(:,:)
    Type(MatterPowerData) :: CPK
    real(mcp), allocatable :: PK(:,:)
    integer nk,nz

    CPK%num_k = Theory%MPK%nx
    CPK%num_z = Theory%MPK%ny

    !Allocate Theory arrays
    allocate(Theory%NL_MPK)
    allocate(Ratios(CPK%num_k,CPK%num_z))

    !Allocate Dummy Pointer and fill with Linear MPK
    allocate(PK(CPK%num_k,CPK%num_z))
    PK=Theory%MPK%z

    allocate(CPK%matpower(CPK%num_k,CPK%num_z))
    allocate(CPK%ddmat(CPK%num_k,CPK%num_z))
    allocate(CPK%nonlin_ratio(CPK%num_k,CPK%num_z))
    allocate(CPK%log_kh(CPK%num_k))
    allocate(CPK%redshifts(CPK%num_z))
    CPK%log_kh = Theory%MPK%x
    CPK%redshifts = Theory%MPK%y
    CPK%matpower = PK

    !need splines to get nonlinear ratios
    call MatterPowerdata_getsplines(CPK)
    call NonLinear_GetRatios(CPK)
    Ratios = CPK%nonlin_ratio
    call MatterPowerdata_Free(CPK)

    PK = PK+2*log(Ratios)
    call Theory%NL_MPK%Init(Theory%MPK%x,Theory%MPK%y,PK)

    end subroutine CAMBCalc_GetNLandRatios

    subroutine CAMBCalc_InitCAMB(this,CMB,error, DoReion)
    class(CAMB_Calculator) :: this
    class(CMBParams), intent(in) :: CMB
    logical, optional, intent(in) :: DoReion
    type(CAMBParams)  P
    integer error

    call this%CMBToCAMB(CMB, P)
    call CAMBParams_Set(P,error,PresentDefault(.true.,DoReion))

    end subroutine CAMBCalc_InitCAMB

    function CAMBCalc_GetOpticalDepth(this,CMB) result(GetOpticalDepth)
    class(CAMB_Calculator) :: this
    class(CMBParams) CMB
    real(mcp) GetOpticalDepth
    type(CAMBParams)  P
    integer error

    call this%CMBToCAMB(CMB, P)
    call CAMBParams_Set(P,error)

    if (error/= 0) then
        GetOpticalDepth = -1
    else
        GetOpticalDepth = Reionization_GetOptDepth(P%Reion, P%ReionHist)
    end if
    end function CAMBCalc_GetOpticalDepth

    function CAMBCalc_GetZreFromTau(this,CMB, tau) result(GetZreFromTau)
    class(CAMB_Calculator) :: this
    class(CMBParams) CMB
    real(mcp), intent(in) :: tau
    real(mcp) GetZreFromTau
    type(CAMBParams)  P

    call this%CMBToCAMB(CMB, P)
    GetZreFromTau = CAMB_GetZreFromTau(P,dble(tau))

    end function CAMBCalc_GetZreFromTau

    function CAMBCalc_CMBToTheta(this,CMB) result(CMBToTheta)
    use ModelParams
    class(CAMB_Calculator) :: this
    class(CMBParams) CMB
    real(mcp) CMBToTheta
    integer error

    call this%InitCAMB(CMB,error,.false.)
    CMBToTheta = CosmomcTheta()

    end function CAMBCalc_CMBToTheta


    real(mcp) function CAMBCalc_BAO_D_v(this, z)
    use CAMB, only : BAO_D_v
    class(CAMB_Calculator) :: this
    real(mcp), intent(IN) :: z

    CAMBCalc_BAO_D_v = BAO_D_v(z)

    end function CAMBCalc_BAO_D_v


    real(mcp) function CAMBCalc_AngularDiameterDistance(this, z)
    use CAMB, only : AngularDiameterDistance  !!angular diam distance also in Mpc no h units
    class(CAMB_Calculator) :: this
    real(mcp), intent(IN) :: z

    CAMBCalc_AngularDiameterDistance = AngularDiameterDistance(z)

    end function CAMBCalc_AngularDiameterDistance

    real(mcp) function CAMBCalc_Hofz(this, z)
    use CAMB, only : Hofz  !!angular diam distance also in Mpc no h units
    class(CAMB_Calculator) :: this
    real(mcp), intent(IN) :: z

    CAMBCalc_Hofz = Hofz(z)

    end function CAMBCalc_Hofz


    subroutine CAMBCalc_InitCAMBParams(this,P)
    use lensing
    use ModelParams
    class(CAMB_Calculator) :: this
    type(CAMBParams)  P
    integer zix
    !JD Changed P%Transfer%redshifts and P%Transfer%num_redshifts to
    !P%Transfer%PK_redshifts and P%Transfer%PK_num_redshifts respectively
    !for nonlinear lensing of CMB + LSS compatibility
    Threadnum =num_threads
    w_lam = -1
    wa_ppf = 0._dl
    call CAMB_SetDefParams(P)

    HighAccuracyDefault = .true.
    P%OutputNormalization = outNone

    !JD added to save computation time when only using MPK
    if(.not. use_CMB) lmax_computed_cl = 10

    P%WantScalars = .true.
    P%WantTensors = compute_tensors
    P%WantTransfer = Use_LSS .or. get_sigma8

    P%Max_l=lmax_computed_cl
    P%Max_eta_k=lmax_computed_cl*2

    P%Max_l_tensor=lmax_tensor
    P%Max_eta_k_tensor=lmax_tensor*5./2

    P%Transfer%k_per_logint=0

    if (use_nonlinear) then
        P%NonLinear = NonLinear_pk
        P%Transfer%kmax = max(1.2_mcp,power_kmax)
    else
        P%Transfer%kmax = max(0.8_mcp,power_kmax)
    end if

    if (AccuracyLevel > 1 .or. HighAccuracyDefault) then
        if (USE_LSS .or. get_sigma8) then
            P%Transfer%high_precision=.true.
            P%Transfer%kmax=P%Transfer%kmax + 0.2
        end if
        AccuracyBoost = AccuracyLevel
        lAccuracyBoost = AccuracyLevel
        lSampleBoost = AccuracyLevel
        P%AccurateReionization = .true.
    end if

    if (max_transfer_redshifts < num_power_redshifts) then
        stop 'Need to manually set max_transfer_redshifts larger in CAMB''s modules.f90'
    end if

    if (use_LSS) then
        P%Transfer%PK_num_redshifts = num_power_redshifts
        do zix=1, num_power_redshifts
            !CAMB's ordering is from highest to lowest
            P%Transfer%PK_redshifts(zix) = power_redshifts(num_power_redshifts-zix+1)
        end do
    else
        P%Transfer%PK_num_redshifts = 1
        P%Transfer%PK_redshifts(1) = 0
    end if

    P%Num_Nu_Massive = 3
    P%Num_Nu_Massless = 0.046
    P%InitPower%nn = 1
    P%AccuratePolarization = num_cls/=1
    P%Reion%use_optical_depth = .false.
    P%OnlyTransfers = .true.

    if (CMB_Lensing) then
        P%DoLensing = .true.
        P%Max_l = lmax_computed_cl +100 + 50 !+50 in case accuracyBoost>1 and so odd l spacing
        P%Max_eta_k = P%Max_l*2
    end if

    if (HighAccuracyDefault) then
        P%Max_eta_k=max(min(P%max_l,3000)*2.5_dl*AccuracyLevel,P%Max_eta_k)
        if (CMB_Lensing .and. use_lensing_potential) P%Max_eta_k = max(P%Max_eta_k, 12000*AccuracyLevel)
        if (CMB_Lensing .and. use_nonlinear_lensing) P%Max_eta_k = max(P%Max_eta_k, 11000*AccuracyLevel)
        !k_etamax=18000 give c_phi_phi accurate to sub-percent at L=1000, <4% at L=2000
        !k_etamax=10000 is just < 1% at L<=500
    end if
    !JD 08/13 for nonlinear lensing of CMB + LSS compatibility
    if (CMB_Lensing .and. use_nonlinear_lensing) then
        P%WantTransfer = .true.
        P%NonLinear = NonLinear_lens
        call Transfer_SetForNonlinearLensing(P%Transfer)
        if(use_nonlinear) P%NonLinear = NonLinear_both
    end if
    call Transfer_SortAndIndexRedshifts(P%Transfer)
    !End JD modifications
    lensing_includes_tensors = .false.

    P%Scalar_initial_condition = initial_vector
    P%InitialConditionVector = 0
    P%InitialConditionVector(initial_adiabatic) = -1

    BackgroundOutputs%z_outputs => z_outputs

    end subroutine CAMBCalc_InitCAMBParams

    !Mapping between array of power spectrum parameters and CAMB
    subroutine CAMBCalc_SetCAMBInitPower(this,P,CMB,in)
    class(CAMB_Calculator) :: this
    type(CAMBParams)  P
    class(CMBParams) CMB
    integer, intent(in) :: in

    if (Power_Name == 'power_tilt') then
        P%InitPower%k_0_scalar = pivot_k
        P%InitPower%k_0_tensor = pivot_k

        P%InitPower%ScalarPowerAmp(in) = cl_norm*CMB%InitPower(As_index)
        P%InitPower%rat(in) = CMB%InitPower(amp_ratio_index)

        P%InitPower%an(in) = CMB%InitPower(1)
        P%InitPower%ant(in) = CMB%InitPower(2)
        if (P%InitPower%rat(in)>0 .and. .not. compute_tensors) &
        call MpiStop('computing r>0 but compute_tensors=F')
        P%InitPower%n_run(in) = CMB%InitPower(3)
        if (inflation_consistency) then
            P%InitPower%ant(in) = - CMB%InitPower(amp_ratio_index)/8.
            !note input n_T is ignored, so should be fixed (to anything)
        end if
    else
        stop 'CAMB_Calculator:Wrong initial power spectrum'
    end if

    end subroutine CAMBCalc_SetCAMBInitPower

    subroutine CAMBCalc_ReadParams(this,Ini)
    class(CAMB_Calculator) :: this
    class(TSettingIni) :: Ini

    call this%TCosmologyCalculator%ReadParams(Ini)
    this%calcName ='CAMB'

    this%CAMB_timing = Ini%Read_Logical('CAMB_timing',.false.)

    this%highL_theory_cl_template_file = Ini%ReadFilename('highL_theory_cl_template',DataDir,.true.)

    if (Ini%HasKey('highL_unlensed_cl_template')) then
        highL_unlensed_cl_template=  Ini%ReadFilename('highL_unlensed_cl_template')
    else
        highL_unlensed_cl_template = concat(LocalDir,'camb/',highL_unlensed_cl_template)
    end if

    end subroutine CAMBCalc_ReadParams


    subroutine CAMBCalc_InitForLikelihoods(this)
    !Called later after likelihoods etc loaded
    class(CAMB_Calculator) :: this

    if (lmax_computed_cl /= lmax) then
        if (compute_tensors .and. lmax_tensor > lmax_computed_cl) call MpiStop('lmax_tensor > lmax_computed_cl')
        call this%LoadFiducialHighLTemplate()
    end if

    call this%InitCAMBParams(this%CAMBP)

    if (Feedback > 0 .and. MPIRank==0) then
        write(*,*) 'max_eta_k         = ', this%CAMBP%Max_eta_k
        write(*,*) 'transfer kmax     = ', this%CAMBP%Transfer%kmax
    end if

    this%CAMBP%WantTensors = compute_tensors

    end subroutine CAMBCalc_InitForLikelihoods


    subroutine CAMBCalc_VersionTraceOutput(this, ReadValues)
    use GaugeInterface, only : Eqns_name
    class(CAMB_Calculator) :: this
    class(TNameValueList) :: ReadValues

    !Store for the record any useful info about version etc.
    call ReadValues%Add( 'Compiled_CAMB_version', version)
    call ReadValues%Add('Compiled_Recombination', Recombination_Name)
    call ReadValues%Add('Compiled_Equations', Eqns_name)
    call ReadValues%Add('Compiled_Reionization', Reionization_Name)
    call ReadValues%Add('Compiled_InitialPower', Power_Name)

    end subroutine CAMBCalc_VersionTraceOutput



    subroutine LoadFiducialHighLTemplate(this)
    class(CAMB_Calculator) :: this
    !This should be a lensed scalar CMB power spectrum, e.g. for including at very high L where foregrounds etc. dominate anyway
    integer L,  status
    real(mcp) array(4), nm
    Type(TTextFile) :: F

    allocate(this%highL_lensedCL_template(2:lmax, num_clsS))
    call F%Open(this%highL_theory_cl_template_file)
    do
        read(F%unit,*, iostat=status) L , array
        if (status/=0 .or. L>lmax) exit
        nm = 2*pi/(l*(l+1))
        if (L>=2) this%highL_lensedCL_template(L,1:num_clsS) = nm*array(TensClOrder(1:num_clsS))
    end do
    call F%Close()

    if (this%highL_lensedCL_template(2,1) < 100) &
    call MpiStop('highL_theory_cl_template must be in muK^2')

    if (L<lmax) call MpiStop('highL_theory_cl_template does not go to lmax')
    if (num_cls_ext>0 .and. MpiRank==0) &
    write(*,*) 'warning: zero padding ext cls in LoadFiducialHighLTemplate'

    end subroutine LoadFiducialHighLTemplate


    !!! CAMBTransferCache

    subroutine CAMBTransferCache_Clear(Info)
    class(CAMBTransferCache) Info

    call CAMB_FreeCAMBdata(Info%Transfers)
    call Info%TTheoryIntermediateCache%Clear()

    end subroutine CAMBTransferCache_Clear

    end module Calculator_CAMB
