    !Use CAMB
    module Calculator_CAMB
    use CosmologyTypes
    use CosmoTheory
    use CAMB, only : CAMB_GetResults, CAMB_GetAge, CAMBParams, CAMB_SetDefParams, &
        AccuracyBoost,  Cl_scalar, Cl_tensor, Cl_lensed, outNone, w_lam, wa_ppf,&
        CAMBParams_Set, MT, CAMBdata, NonLinear_Pk, Nonlinear_lens, Reionization_GetOptDepth, CAMB_GetZreFromTau, &
        CAMB_GetTransfers,CAMB_FreeCAMBdata,CAMB_InitCAMBdata, CAMB_TransfersToPowers, Transfer_SetForNonlinearLensing, &
        initial_adiabatic,initial_vector,initial_iso_baryon,initial_iso_CDM, initial_iso_neutrino, initial_iso_neutrino_vel, &
        HighAccuracyDefault, highL_unlensed_cl_template, ThermoDerivedParams, nthermo_derived, BackgroundOutputs, &
        Transfer_SortAndIndexRedshifts,  &
        Recombination_Name, reionization_name, power_name, threadnum, version, tensor_param_rpivot
    use Errors !CAMB
    use settings
    use likelihood
    use Calculator_Cosmology
    use GeneralTypes
    use, intrinsic :: ieee_arithmetic
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
        real(mcp) :: k_eta_max_scalar = -1._mcp
        logical :: accurate_BB =.false.
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
    procedure :: GetNLandRatios => CAMBCalc_GetNLandRatios
    !Overridden inherited
    procedure :: ReadParams => CAMBCalc_ReadParams
    procedure :: InitForLikelihoods => CAMBCalc_InitForLikelihoods
    procedure :: BAO_D_v => CAMBCalc_BAO_D_v
    procedure :: AngularDiameterDistance => CAMBCalc_AngularDiameterDistance
    procedure :: ComovingRadialDistance => CAMBCalc_ComovingRadialDistance
    procedure :: AngularDiameterDistance2 => CAMBCalc_AngularDiameterDistance2
    procedure :: LuminosityDistance => CAMBCalc_LuminosityDistance
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


    !    integer, parameter :: ScalClOrder(5) = (/1,3,2,4,5/), TensClOrder(4) = (/1,4,2,3/)
    !Mapping of CAMB CL array ordering to TT , TE, EE, BB, phi, phiT


    public CAMB_Calculator
    contains

    subroutine CAMBCalc_CMBToCAMB(this,CMB,P)
    use LambdaGeneral
    use camb, only: CAMB_SetNeutrinoHierarchy
    use CAMBmain, only : ALens
    use constants, only : default_nnu,delta_mnu21,delta_mnu31,mnu_min_normal
    use lensing, only : ALens_Fiducial
    use MassiveNu, only : sum_mnu_for_m1
    class(CAMB_Calculator) :: this
    class(CMBParams) CMB
    type(CAMBParams)  P
    real(dl) neff_massive_standard, mnu, m1, m3, normal_frac
    real(dl), external :: Newton_raphson

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
    ALens_Fiducial = CMB%ALensf
    P%InitialConditionVector(initial_iso_CDM) = &
        sign(sqrt(abs(CMB%iso_cdm_correlated) /(1-abs(CMB%iso_cdm_correlated))),CMB%iso_cdm_correlated)
    P%Num_Nu_Massive = 0
    P%Nu_mass_numbers = 0
    P%Num_Nu_Massless = CMB%nnu
    P%share_delta_neff = .false.
    if (CMB%omnuh2>0) then
        call CAMB_SetNeutrinoHierarchy(P, CMB%omnuh2, CMB%omnuh2_sterile, CMB%nnu, &
            CosmoSettings%neutrino_hierarchy, CosmoSettings%num_massive_neutrinos)
    end if

    P%YHe = CMB%YHe
#ifdef COSMOREC
    if (CMB%fdm/=0._mcp) P%Recomb%runmode = 3
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
    Theory%numderived = nthermo_derived + noutputs*4
    if (Theory%numderived > max_derived_parameters) &
        call MpiStop('numderived > max_derived_parameters: increase in CosmologyTypes.f90')
    Theory%derived_parameters(1:nthermo_derived) = ThermoDerivedParams(1:nthermo_derived)
    do i=1, noutputs
        Theory%derived_parameters(nthermo_derived+(i-1)*4+1) = BackgroundOutputs%rs_by_D_v(i)
        Theory%derived_parameters(nthermo_derived+(i-1)*4+2) = BackgroundOutputs%H(i)*const_c/1e3_mcp
        Theory%derived_parameters(nthermo_derived+(i-1)*4+3) = BackgroundOutputs%DA(i)
        Theory%derived_parameters(nthermo_derived+(i-1)*4+4) = (1+BackgroundOutputs%z_outputs(i))* &
            BackgroundOutputs%DA(i) * BackgroundOutputs%H(i) !F_AP parameter
    end do
    end subroutine CAMBCalc_SetDerived

    subroutine CAMBCalc_SetParamsForBackground(this,CMB)
    use Camb, only: CAMBParams_Set
    class(CAMB_Calculator) :: this
    class(CMBParams) CMB
    Type(CAMBParams)  P

    !set background parameters, but don't calculate thermal history
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
        write(LogFile%unit,'("CAMB called ",I0," times; ",I0," errors")') this%ncalls, this%nerrors
    end if
    if (Feedback > 1) write (*,*) 'CAMB done'

    end subroutine CAMBCalc_GetNewTransferData

    subroutine CAMBCalc_GetNewPowerData(this, CMB, Info, Theory, error)
    class(CAMB_Calculator) :: this
    class(CMBParams) :: CMB
    class(TTheoryIntermediateCache), pointer :: Info
    class(TCosmoTheoryPredictions) :: Theory
    integer error,i,j

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
        if(CosmoSettings%use_CMB)then
            call this%SetPowersFromCAMB(CMB,Theory)
            do i=1, min(3,CosmoSettings%num_cls)
                if (CosmoSettings%cl_lmax(i,i)>0) then
                    if (any(Theory%cls(i,i)%Cl(:) < 0 )) then
                        error = 1
                        call MpiStop('Calculator_CAMB: negative C_l (could edit to silent error here)')
                        return
                    end if
                end if
                do j= i, 1, -1
                    if (CosmoSettings%cl_lmax(i,j)>0) then
                        if ( any(ieee_is_nan(Theory%cls(i,j)%Cl))) then
                            error=1
                            write(*,*) 'WARNING: NaN CL?', i, j
                            return
                        end if
                    end if
                end do
            end do
        end if

        !redshifts are in increasing order, so last index is redshift zero
        if (CosmoSettings%Use_LSS .or. CosmoSettings%get_sigma8) then
            Theory%sigma_8 = Info%Transfers%MTrans%sigma_8(size(Info%Transfers%MTrans%sigma_8,1),1)
        else
            Theory%sigma_8 = 0
        end if

        if (CosmoSettings%Use_LSS) then
            call this%SetPkFromCAMB(Info%Transfers%MTrans,Theory,error)
            if (error/=0) then
                write(*,*) 'WARNING: NaN PK?'
                return
            end if
        end if
    end select

    end subroutine CAMBCalc_GetNewPowerData

    subroutine CAMBCalc_GetTheoryForImportance(this, CMB, Theory, error)
    use ModelParams, only : ThreadNum
    class(CAMB_Calculator) :: this
    class(CMBParams) CMB
    class(TCosmoTheoryPredictions) Theory
    class(CAMBTransferCache), allocatable :: Info
    integer error,i,j
    logical :: DoCls, DoPk
    type(CAMBParams)  P

    error = 0
    DoCls = this%ImportanceOptions%redo_cls
    DoPk = this%ImportanceOptions%redo_pk

    if (DoCls .or. DoPk) then
        allocate(Info)
        call CAMB_InitCAMBdata(Info%Transfers)
        call this%CMBToCAMB(CMB, P)
        Threadnum =num_threads

        P%WantCls = DoCls

        if (.not. DoPk .and. .not. (CosmoSettings%CMB_Lensing .and. &
            CosmoSettings%use_nonlinear_lensing)) P%WantTransfer = .false.

        if (this%CAMB_timing) call Timer()
        if (Feedback > 1) write (*,*) 'Calling CAMB'
        call CAMB_GetTransfers(P, Info%Transfers, error)
        if (Feedback > 1) write (*,*) 'CAMB Done'
        if (this%CAMB_timing) call Timer('CAMB_GetTransfers')

        if (error==0) then
            call this%SetCAMBInitPower(Info%Transfers%Params,CMB,1)
            call CAMB_TransfersToPowers(Info%Transfers)
            error=global_error_flag
        end if
    else
        call this%GetNewBackgroundData(CMB,Theory,error)
    end if

    if (DoCls .and. error==0) then
        call this%SetPowersFromCAMB(CMB,Theory)
        if (any(Theory%cls(1,1)%Cl(:) < 0 )) then
            error = 1
            call MpiStop('Calculator_CAMB: negative C_l (could edit to silent error here)')
        end if
        do i=1, min(3,CosmoSettings%num_cls)
            if(error/=0) exit
            do j= i, 1, -1
                if (CosmoSettings%cl_lmax(i,j)>0) then
                    if ( any(ieee_is_nan(Theory%cls(i,j)%Cl))) then
                        error=1
                        write(*,*) 'WARNING: NaN CL?'
                        exit
                    end if
                end if
            end do
        end do
    end if

    if (DoPK .and. error==0) then
        Theory%sigma_8 = Info%Transfers%MTrans%sigma_8(size(Info%Transfers%MTrans%sigma_8,1),1)
        call this%SetPkFromCAMB(Info%Transfers%MTrans,Theory,error)
        if (error/=0) write(*,*) 'WARNING: NaN PK?'
    end if

    if (error==0) call this%SetDerived(Theory)

    if (DoCls .or. DoPk) call Info%Clear()

    end subroutine CAMBCalc_GetTheoryForImportance

    subroutine CAMBCalc_SetPowersFromCAMB(this,CMB,Theory)
    use constants
    use InitialPower
    use ModelData
    class(CAMB_Calculator) :: this
    class(CMBParams) :: CMB
    class(TCosmoTheoryPredictions), target :: Theory
    real(mcp), parameter :: cons =  (COBE_CMBTemp*1e6)**2
    integer l
    real(mcp) :: highL_norm = 0
    real(mcp) lens_recon_scale, rms
    integer i,j, lmx, lmaxCL
    integer, save, allocatable :: indicesS(:,:), indicesT(:,:)

    if (.not. allocated(indicesS)) then
        allocate(indicesS(3,3))
        allocate(indicesT(3,3))
        indicesS=0
        indicesS(1,1) =  C_Temp
        indicesS(2,2) =  C_E
        indicesT = indicesS
        indicesT(2,1) =  CT_Cross
        indicesT(3,3) =  CT_B
        indicesS(2,1) =  C_Cross
    end if

    lens_recon_scale = CMB%InitPower(Aphiphi_index)

    do i=1, min(3,CosmoSettings%num_cls)
        do j= i, 1, -1
            lmaxCL = CosmoSettings%cl_lmax(i,j)
            lmx = min(CosmoSettings%lmax_computed_cl, lmaxCL)
            if (lmx/=0) then
                associate( CL => Theory%Cls(i,j)%CL)
                    if (indicesT(i,j)==0) then
                        CL=0
                    else
                        if (CosmoSettings%CMB_Lensing) then
                            CL(2:lmx) = cons*Cl_lensed(2:lmx,1, indicesT(i,j))
                        else
                            if (indicesS(i,j)/=0) CL(2:lmx) = cons*Cl_Scalar(2:lmx,1, indicesS(i,j))
                        end if
                        if (CosmoSettings%lmax_computed_cl < lmaxCL) then
                            if (highL_norm ==0) & !normally normalize off TT
                            & highL_norm = CL(lmx)/this%highL_lensedCL_template(lmx,indicesT(i,j))
                            CL(lmx+1:lmaxCL) =  highL_norm*this%highL_lensedCL_template(lmx+1:lmaxCL,indicesT(i,j))
                        end if
                        if (CosmoSettings%compute_tensors) then
                            lmx = min(lmx,CosmoSettings%lmax_tensor)
                            CL(2:lmx) =  CL(2:lmx) + cons*Cl_tensor(2:lmx,1, indicesT(i,j))
                        end if
                    end if
                end associate
            end if
        end do
    end do

    if (CosmoSettings%use_lensing_potential) then
        !CMB lensing potential
        !in camb Cphi is l^4 C_l, we want [l(l+1)]^2Cphi/2pi
        lmx = min(CosmoSettings%lmax_computed_cl, CosmoSettings%cl_lmax(CL_Phi,CL_Phi))
        if (lmx/=0) then
            if (.not. CosmoSettings%CMB_lensing) call MpiStop('Must have lensing on to use lensing potential')
            associate(CL=> Theory%Cls(CL_Phi,CL_Phi)%CL)
                do l=2, lmx
                    CL(L) =  Cl_scalar(L,1, C_Phi)*(real(l+1)**2/l**2)/twopi * lens_recon_scale
                end do
                CL(lmx+1:)=0
            end associate
        end if
        lmx = min(CosmoSettings%lmax_computed_cl, CosmoSettings%cl_lmax(CL_Phi,CL_T))
        if (lmx/=0) then
            !lensing-temp
            do l=2, lmx
                Theory%Cls(CL_phi,CL_T)%CL = Cl_scalar(l,1, C_PhiTemp)/real(l)**3 * sqrt(lens_recon_scale)
            end do
        end if
    end if

    if (CosmoSettings%CMB_Lensing .and. this%CAMBP%max_l>=2000) then
        !Get RMS deflection angle in arcmin
        rms=0
        do L=2, 2000
            rms = rms +  Cl_scalar(L,1, C_Phi)*(real(l+1)**2/l**2)/twopi*(L+0.5_mcp)/(L*(L+1))
        end do
        Theory%Lensing_rms_deflect = sqrt(rms)*180/pi*60
    else
        Theory%Lensing_rms_deflect = 0
    end if

    if (CosmoSettings%compute_tensors) then
        Theory%tensor_ratio_02 = TensorPower(0.002d0,1)/ScalarPower(0.002d0,1)
        Theory%tensor_AT = TensorPower(CosmoSettings%tensor_pivot_k,1)
        Theory%tensor_ratio_BB = TensorPower(0.01d0,1)/ScalarPower(0.01d0,1)
        Theory%tensor_ratio_C10 = Cl_tensor(10, 1, 1)/Cl_scalar(10,1, 1)
    else
        Theory%tensor_ratio_02 = 0
        Theory%tensor_ratio_BB = 0
        Theory%tensor_ratio_C10 = 0
        Theory%tensor_AT = 0
    end if

    end subroutine CAMBCalc_SetPowersFromCAMB

    subroutine CAMBCalc_SetPkFromCAMB(this,M,Theory,error)
    use Transfer
    use camb, only : CP
    class(CAMB_Calculator) :: this
    class(TCosmoTheoryPredictions) Theory
    Type(MatterTransferData) M
    integer :: error
    real(mcp), allocatable :: k(:), z(:), PK(:,:)
    integer zix,nz,nk, nR
    real(mcp), allocatable :: NL_Ratios(:,:)
    real(mcp) :: dR, R, minR
    integer i

    !Free theory arrays as they may resize between samples
    call Theory%FreePK()

    nz=CP%Transfer%PK_num_redshifts

    if (.not. allocated(Theory%growth_z)) allocate(Theory%growth_z, Theory%sigma8_z)
    call Theory%growth_z%InitForSize(nz)
    call Theory%sigma8_z%InitForSize(nz)
    allocate(z(nz))

    do zix=1,nz
        z(zix) = CP%Transfer%PK_redshifts(nz-zix+1)
        Theory%sigma8_z%F(zix) = M%sigma_8(nz-zix+1,1)
        Theory%growth_z%F(zix) = M%sigma2_vdelta_8(nz-zix+1,1)/M%sigma_8(nz-zix+1,1)
    end do
    Theory%sigma8_z%X=z
    Theory%growth_z%X=z

    if (CosmoSettings%use_matterpower) then
        nk=M%num_q_trans
        nz=CP%Transfer%PK_num_redshifts
        allocate(PK(nk,nz))
        allocate(k(nk))

        k = log(M%TransferData(Transfer_kh,:,1))

        call Transfer_GetUnsplinedPower(M, PK,transfer_power_var,transfer_power_var)
        PK = Log(PK)
        if (any(ieee_is_nan(PK))) then
            error = 1
            return
        end if
        allocate(Theory%MPK)
        call Theory%MPK%Init(k,z,PK)
    end if


    if (CosmoSettings%use_Weylpower) then
        call Transfer_GetUnsplinedPower(M, PK,transfer_Weyl,transfer_Weyl,hubble_units=.false.)
        PK = Log(PK)
        if (any(ieee_is_nan(PK))) then
            error = 1
            return
        end if
        allocate(Theory%MPK_WEYL)
        call Theory%MPK_WEYL%Init(k,z,PK)
    end if

    if (CosmoSettings%use_SigmaR) then
        !Note R is in k*h units
        dR = log(1.2)/AccuracyLevel
        minR = 1/CP%Transfer%kmax
        nR = nint(log(150/minR)/dR) +1
        if (.not. allocated(Theory%Sigma_R)) allocate(Theory%Sigma_R)
        call Theory%Sigma_R%InitForSize(nR)
        do i=1, nR
            Theory%Sigma_R%X(i) = exp((i-1)*dR)*minR
        end do
        call Transfer_GetSigmaRArray(M, Theory%Sigma_R%X, Theory%Sigma_R%F, &
            var1 = transfer_nonu,var2=transfer_nonu)
    end if

    if(CosmoSettings%use_nonlinear)then
        call this%GetNLandRatios(M,Theory,NL_Ratios,error)
        if(error/=0) return
    end if

    end subroutine CAMBCalc_SetPkFromCAMB


    subroutine CAMBCalc_GetNLandRatios(this,M,Theory,Ratios,error)
    use Transfer
    class(CAMB_Calculator) :: this
    class(TCosmoTheoryPredictions) Theory
    Type(MatterTransferData) M
    real(mcp), allocatable, intent(out) :: Ratios(:,:)
    Type(MatterPowerData) :: CPK
    real(mcp), allocatable :: PK(:,:)
    integer error,zix,nz

    CPK%num_k = Theory%MPK%nx
    CPK%num_z = Theory%MPK%ny

    !Allocate Theory arrays
    allocate(Theory%NL_MPK)
    allocate(Ratios(CPK%num_k,CPK%num_z))

    !fill PK with Linear MPK
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
    if (any(ieee_is_nan(PK))) then
        error = 1
        return
    end if
    call Theory%NL_MPK%Init(Theory%MPK%x,Theory%MPK%y,PK)

    if (allocated(Theory%MPK_WEYL)) then
        !Assume Weyl scales the same way under non-linear correction
        allocate(Theory%NL_MPK_WEYL)
        PK = Theory%MPK_WEYL%z + 2*log(Ratios)
        call Theory%NL_MPK_WEYL%Init(Theory%MPK%x,Theory%MPK%y,PK)
    end if

    end subroutine CAMBCalc_GetNLandRatios

    subroutine CAMBCalc_InitCAMB(this,CMB,error, DoReion)
    class(CAMB_Calculator) :: this
    class(CMBParams), intent(in) :: CMB
    logical, optional, intent(in) :: DoReion
    logical WantReion
    type(CAMBParams)  P
    integer error

    if (present(DoReion)) then
        WantReion = DoReion
    else
        WantReion = .true.
    end if

    call this%CMBToCAMB(CMB, P)
    call CAMBParams_Set(P,error,WantReion)

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

    real(mcp) function CAMBCalc_ComovingRadialDistance(this, z)
    use CAMB, only : ComovingRadialDistance  !!comoving radial distance also in Mpc no h units
    class(CAMB_Calculator) :: this
    real(mcp), intent(IN) :: z

    CAMBCalc_ComovingRadialDistance = ComovingRadialDistance(z)

    end function CAMBCalc_ComovingRadialDistance

    real(mcp) function CAMBCalc_AngularDiameterDistance2(this, z1, z2)
    use CAMB, only : AngularDiameterDistance2  !!angular diam distance also in Mpc no h units
    class(CAMB_Calculator) :: this
    real(mcp), intent(IN) :: z1, z2

    CAMBCalc_AngularDiameterDistance2 = AngularDiameterDistance2(z1, z2)

    end function CAMBCalc_AngularDiameterDistance2

    real(mcp) function CAMBCalc_LuminosityDistance(this, z)
    use CAMB, only : LuminosityDistance  !! distance also in Mpc no h units
    class(CAMB_Calculator) :: this
    real(mcp), intent(IN) :: z

    CAMBCalc_LuminosityDistance = LuminosityDistance(z)

    end function CAMBCalc_LuminosityDistance

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

    !JD Modified to save computation time when only using MPK
    if(CosmoSettings%use_CMB) then
        P%WantScalars = .true.
        P%WantTensors = CosmoSettings%compute_tensors
        P%Max_l=CosmoSettings%lmax_computed_cl
        P%Max_eta_k=CosmoSettings%lmax_computed_cl*2
        P%Max_l_tensor=CosmoSettings%lmax_tensor
        P%Max_eta_k_tensor=CosmoSettings%lmax_tensor*5./2
    else
        P%WantCls = .false.
    end if

    P%WantTransfer = CosmoSettings%Use_LSS .or. CosmoSettings%get_sigma8
    P%Transfer%k_per_logint=0

    if (CosmoSettings%use_nonlinear) then
        P%NonLinear = NonLinear_pk
        P%Transfer%kmax = max(1.2_mcp,CosmoSettings%power_kmax)
    else
        P%Transfer%kmax = max(0.8_mcp,CosmoSettings%power_kmax)
    end if

    if (AccuracyLevel > 1 .or. HighAccuracyDefault) then
        if (CosmoSettings%Use_LSS .or. CosmoSettings%get_sigma8) then
            P%Transfer%high_precision=.true.
            P%Transfer%kmax=P%Transfer%kmax + 0.2
        end if
        AccuracyBoost = AccuracyLevel
        lAccuracyBoost = AccuracyLevel
        lSampleBoost = AccuracyLevel
        P%AccurateReionization = .true.
    end if

    P%AccurateBB = this%accurate_BB

    if (max_transfer_redshifts < CosmoSettings%num_power_redshifts) then
        stop 'Need to manually set max_transfer_redshifts larger in CAMB''s modules.f90'
    end if

    if (CosmoSettings%num_power_redshifts>1) then
        P%Transfer%PK_num_redshifts = CosmoSettings%num_power_redshifts
        do zix=1, CosmoSettings%num_power_redshifts
            !CAMB's ordering is from highest to lowest
            P%Transfer%PK_redshifts(zix) = CosmoSettings%power_redshifts(CosmoSettings%num_power_redshifts-zix+1)
        end do
    else
        P%Transfer%PK_num_redshifts = 1
        P%Transfer%PK_redshifts(1) = 0
    end if

    P%Num_Nu_Massive = 3
    P%Num_Nu_Massless = 0.046
    P%InitPower%nn = 1
    P%AccuratePolarization = CosmoSettings%num_cls/=1
    P%Reion%use_optical_depth = .false.
    P%OnlyTransfers = .true.

    if (CosmoSettings%CMB_Lensing) then
        P%DoLensing = .true.
        P%Max_l = CosmoSettings%lmax_computed_cl +100 + 50 !+50 in case accuracyBoost>1 and so odd l spacing
        P%Max_eta_k = P%Max_l*2
    end if

    if (HighAccuracyDefault) then
        P%Max_eta_k=max(min(P%max_l,3000)*2.5_dl*AccuracyLevel,P%Max_eta_k)
        if (CosmoSettings%CMB_Lensing .and. (CosmoSettings%use_lensing_potential .or. CosmoSettings%use_nonlinear_lensing)) &
            P%Max_eta_k = max(P%Max_eta_k, 14000*AccuracyLevel)
        !k_etamax=18000 give c_phi_phi accurate to sub-percent at L=1000, <4% at L=2000
        !k_etamax=10000 is just < 1% at L<=500
    end if
    if (this%k_eta_max_scalar>0) then
        P%Max_eta_k = this%k_eta_max_scalar
    end if
    !JD 08/13 for nonlinear lensing of CMB + LSS compatibility
    if (CosmoSettings%CMB_Lensing .and. CosmoSettings%use_nonlinear_lensing) then
        P%WantTransfer = .true.
        P%NonLinear = NonLinear_lens
        call Transfer_SetForNonlinearLensing(P%Transfer)
        if(CosmoSettings%use_nonlinear) P%NonLinear = NonLinear_both
    end if
    call Transfer_SortAndIndexRedshifts(P%Transfer)
    !End JD modifications
    lensing_includes_tensors = .false.

    P%Scalar_initial_condition = initial_vector
    P%InitialConditionVector = 0
    P%InitialConditionVector(initial_adiabatic) = -1

    BackgroundOutputs%z_outputs => CosmoSettings%z_outputs

    end subroutine CAMBCalc_InitCAMBParams

    !Mapping between array of power spectrum parameters and CAMB
    subroutine CAMBCalc_SetCAMBInitPower(this,P,CMB,ix)
    class(CAMB_Calculator) :: this
    type(CAMBParams)  P
    class(CMBParams) CMB
    integer, intent(in) :: ix

    if (Power_Name == 'power_tilt') then
        P%InitPower%k_0_scalar = CosmoSettings%pivot_k
        P%InitPower%k_0_tensor = CosmoSettings%tensor_pivot_k
        if (P%InitPower%k_0_tensor/=P%InitPower%k_0_scalar) P%InitPower%tensor_parameterization = tensor_param_rpivot
        P%InitPower%ScalarPowerAmp(ix) = cl_norm*CMB%InitPower(As_index)
        P%InitPower%rat(ix) = CMB%InitPower(amp_ratio_index)

        P%InitPower%an(ix) = CMB%InitPower(ns_index)

        if (P%InitPower%rat(ix)>0 .and. .not. CosmoSettings%compute_tensors) &
            call MpiStop('computing r>0 but compute_tensors=F')
        P%InitPower%n_run(ix) = CMB%InitPower(nrun_index)
        P%InitPower%n_runrun(ix) = CMB%InitPower(nrunrun_index)

        if (CosmoSettings%inflation_consistency) then
            if (CMB%InitPower(nt_index)/=0 .or. CMB%InitPower(ntrun_index)/=0) &
                & call MpiStop('Error: inflation_consistency but n_t not set to zero')
            ! first-order consistency relation
            !P%InitPower%ant(ix) = - CMB%InitPower(amp_ratio_index)/8
            !next order consistency relation
            P%InitPower%ant(ix) = - CMB%InitPower(amp_ratio_index)/8*(2-CMB%InitPower(ns_index) - CMB%InitPower(amp_ratio_index)/8)
            P%InitPower%nt_run(ix) = CMB%InitPower(amp_ratio_index)/8* &
                & (CMB%InitPower(amp_ratio_index)/8 + CMB%InitPower(ns_index) - 1)
            !note input n_T, nt run is ignored, so should be fixed
        else
            P%InitPower%ant(ix) = CMB%InitPower(nt_index)
            P%InitPower%nt_run(ix) = CMB%InitPower(ntrun_index)
        end if
    else
        stop 'CAMB_Calculator:Wrong initial power spectrum'
    end if

    end subroutine CAMBCalc_SetCAMBInitPower

    subroutine CAMBCalc_ReadParams(this,Ini)
    use NonLinear
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

    this%k_eta_max_scalar = Ini%Read_Double('k_eta_max_scalar',-1._mcp)
    this%accurate_BB = Ini%Read_Logical('accurate_BB',.false.)

    halofit_version = Ini%Read_Int('halofit_version',halofit_default)

    end subroutine CAMBCalc_ReadParams


    subroutine CAMBCalc_InitForLikelihoods(this)
    !Called later after likelihoods etc loaded
    class(CAMB_Calculator) :: this

    if (CosmoSettings%use_CMB .and. CosmoSettings%lmax_computed_cl /= CosmoSettings%lmax) then
        if (CosmoSettings%compute_tensors .and. CosmoSettings%lmax_tensor > CosmoSettings%lmax_computed_cl) &
            & call MpiStop('lmax_tensor > lmax_computed_cl')
        call this%LoadFiducialHighLTemplate()
    end if

    call this%InitCAMBParams(this%CAMBP)

    if (Feedback > 0 .and. MPIRank==0) then
        if(CosmoSettings%use_CMB) write(*,*) 'max_eta_k         = ', real(this%CAMBP%Max_eta_k)
        write(*,*) 'transfer kmax     = ', real(this%CAMBP%Transfer%kmax)
    end if

    this%CAMBP%WantTensors = CosmoSettings%compute_tensors

    end subroutine CAMBCalc_InitForLikelihoods


    subroutine CAMBCalc_VersionTraceOutput(this, ReadValues)
    use GaugeInterface, only : Eqns_name
    class(CAMB_Calculator) :: this
    class(TNameValueList) :: ReadValues

    !Store for the record any useful info about version etc.
    call ReadValues%Add('Compiled_CAMB_version', version)
    call ReadValues%Add('Compiled_Recombination', Recombination_Name)
    call ReadValues%Add('Compiled_Equations', Eqns_name)
    call ReadValues%Add('Compiled_Reionization', Reionization_Name)
    call ReadValues%Add('Compiled_InitialPower', Power_Name)

    end subroutine CAMBCalc_VersionTraceOutput



    subroutine LoadFiducialHighLTemplate(this)
    class(CAMB_Calculator) :: this
    !This should be a lensed scalar CMB power spectrum, e.g. for including at very high L where foregrounds etc. dominate anyway
    integer L,  status
    real(mcp) array(4)
    Type(TTextFile) :: F

    allocate(this%highL_lensedCL_template(2:CosmoSettings%lmax, 4))
    call F%Open(this%highL_theory_cl_template_file)
    do
        read(F%unit,*, iostat=status) L , array
        if (status/=0 .or. L>CosmoSettings%lmax) exit
        if (L>=2) this%highL_lensedCL_template(L,:) = array
    end do
    call F%Close()

    if (this%highL_lensedCL_template(2,1) < 100) &
        call MpiStop('highL_theory_cl_template must be in muK^2')

    if (L<CosmoSettings%lmax) call MpiStop('highL_theory_cl_template does not go to lmax')

    end subroutine LoadFiducialHighLTemplate


    !!! CAMBTransferCache

    subroutine CAMBTransferCache_Clear(Info)
    class(CAMBTransferCache) Info

    call CAMB_FreeCAMBdata(Info%Transfers)
    call Info%TTheoryIntermediateCache%Clear()

    end subroutine CAMBTransferCache_Clear

    end module Calculator_CAMB
