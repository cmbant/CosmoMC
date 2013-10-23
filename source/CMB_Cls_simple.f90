    !Use CAMB
    module CMB_Cls
    use cmbtypes
    use CAMB, only : CAMB_GetResults, CAMB_GetAge, CAMBParams, CAMB_SetDefParams, &
    AccuracyBoost,  Cl_scalar, Cl_tensor, Cl_lensed, outNone, w_lam, wa_ppf,&
    CAMBParams_Set, MT, CAMBdata, NonLinear_Pk, Nonlinear_lens, Reionization_GetOptDepth, CAMB_GetZreFromTau, &
    CAMB_GetTransfers,CAMB_FreeCAMBdata,CAMB_InitCAMBdata, CAMB_TransfersToPowers, Transfer_SetForNonlinearLensing, &
    initial_adiabatic,initial_vector,initial_iso_baryon,initial_iso_CDM, initial_iso_neutrino, initial_iso_neutrino_vel, &
    HighAccuracyDefault, highL_unlensed_cl_template, ThermoDerivedParams, nthermo_derived, BackgroundOutputs, &
    Transfer_SortAndIndexRedshifts  !JD added for nonlinear lensing of CMB + MPK compatibility
    use Errors !CAMB
    use settings
    use IO
    use likelihood
    use powerspec
    implicit none

    logical :: CMB_lensing = .false.
    logical :: use_lensing_potential = .false.
    !logical :: use_nonlinear = .false.  !JD 08/13 moved to settings, needed  WiggleZ module
    logical :: use_nonlinear_lensing = .false.
    real(mcp) :: lens_recon_scale = 1._mcp

    Type ParamSetInfo
        Type (CAMBdata) :: Transfers
        real(mcp) lastParamArray(max_num_params)
        logical :: validInfo = .false.
    end Type ParamSetInfo

    integer :: lmax_computed_cl = lmax !value used for CAMB

    integer, parameter :: ScalClOrder(5) = (/1,3,2,4,5/), TensClOrder(4) = (/1,4,2,3/)
    !Mapping of CAMB CL array ordering to TT , TE, EE, BB, phi, phiT
    integer :: ncalls = 0
    integer :: nerrors = 0
    type(CAMBParams)  CAMBP

    real(mcp), allocatable :: highL_lensedCL_template(:,:)

    contains
    subroutine CMBToCAMB(CMB,P)
    use LambdaGeneral
    use CAMBmain, only : ALens
    use constants, only : default_nnu
    type(CMBParams) CMB
    type(CAMBParams)  P
    real(dl) neff_massive_standard

    P = CAMBP
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
    call SetCAMBInitPower(P,CMB,1)

    end subroutine CMBToCAMB


    subroutine SetTheoryForBackground(CMB)
    use Camb, only: CAMBParams_Set
    type(CMBParams) CMB
    type(CAMBParams)  P
    !set background dparameters, but don't calculate thermal history
    call CMBToCAMB(CMB, P)
    call CAMBParams_Set(P)

    end subroutine SetTheoryForBackground

    subroutine GetNewBackgroundData(CMB,Theory,error)
    use cambmain, only: initvars
    type(CMBParams) CMB
    integer error
    Class(TheoryPredictions) Theory

    call SetTheoryForBackground(CMB)
    select type (Parameterization)
    class is (CosmologyParameterization)
        if (.not. Parameterization%late_time_only) then
            call InitVars !calculate thermal history, e.g. z_drag etc.
            if (global_error_flag/=0) then
                error=global_error_flag
                return
            end if
            call SetDerived(Theory)
        end if
        class default
        call MpiStop('CMB_Cls_Simple: Must have CosmologyParameterization')
    end select

    end subroutine GetNewBackgroundData

    subroutine SetDerived(Theory)
    Class(TheoryPredictions) Theory
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
    end subroutine SetDerived

    subroutine GetNewTransferData(CMB,Info,Theory,error)
    use ModelParams, only : ThreadNum
    use InitialPower
    type(CMBParams) CMB
    integer error
    Type(ParamSetInfo) Info
    Type(TheoryPredictions) Theory
    type(CAMBParams)  P
    character(LEN=128) :: LogLine

    call CAMB_InitCAMBdata(Info%Transfers)
    call CMBToCAMB(CMB, P)

    if (Feedback > 1) write (*,*) 'Calling CAMB'
    Threadnum =num_threads

    call CAMB_GetTransfers(P, Info%Transfers, error)
    if (error==0) then
        call SetDerived(Theory)
    else
        if (stop_on_error) call MpiStop('CAMB error '//trim(global_error_message))
        if (Feedback > 0) write(*,*) 'CAMB returned error '//trim(global_error_message)
        nerrors=nerrors+1
    end if
    ncalls=ncalls+1
    if (mod(ncalls,100)==0 .and. logfile_unit/=0) then
        write (logLine,*) 'CAMB called ',ncalls, ' times; ', nerrors,' errors'
        call IO_WriteLog(logfile_unit,logLine)
    end if
    if (Feedback > 1) write (*,*) 'CAMB done'

    end subroutine GetNewTransferData

    subroutine GetNewPowerData(CMB, Info, Theory, error)
    use ModelParams, only : ThreadNum
    type(CMBParams) CMB
    integer error
    Type(ParamSetInfo) Info
    Type(TheoryPredictions) Theory

    call SetCAMBInitPower(Info%Transfers%Params,CMB,1)
    call CAMB_TransfersToPowers(Info%Transfers)
    !this sets slow CAMB params correctly from value stored in Transfers
    if (global_error_flag/=0) then
        error=global_error_flag
        return
    end if
    !JD 08/13 added so we dont have to fill Cls unless using CMB
    if(use_CMB)then 
        call SetPowersFromCAMB(Theory)

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
        call SetPkFromCAMB(Theory,Info%Transfers%MTrans)
    end if

    end subroutine GetNewPowerData

    subroutine GetTheoryForImportance(Params, Theory, error, DoCls, DoPk)
    use ModelParams, only : ThreadNum
    type(CMBParams) CMB
    Type(TheoryPredictions) Theory
    integer error
    logical, intent(in) :: DoCls, DoPk
    real(mcp):: Params(:)
    type(CAMBParams)  P

    call Parameterization%ParamArrayToTheoryParams(Params,CMB)
    error = 0

    if (DoPk .or. DoCls) then
        Threadnum =num_threads
        call CMBToCAMB(CMB, P)
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

        call CAMB_GetResults(P)

        error = global_error_flag !using error optional parameter gives seg faults on SGI
    else
        call GetNewBackgroundData(CMB,Theory,error)
    end if
    if (error==0) then
        if (DoCls) call SetPowersFromCAMB(Theory)
        if (DoPK) then
            Theory%sigma_8 = MT%sigma_8(size(MT%sigma_8,1),1)
            call SetPkFromCAMB(Theory,MT)
        end if
        call SetDerived(Theory)
    end if
    end subroutine GetTheoryForImportance

    subroutine SetPowersFromCAMB(Theory)
    use constants
    use InitialPower
    Type(TheoryPredictions) Theory
    real(mcp), parameter :: cons =  (COBE_CMBTemp*1e6)**2*2*pi
    real(mcp) nm
    integer l
    real(mcp) highL_norm

    Theory%cl=0
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
        highL_norm = Theory%cl(lmax_computed_cl,1)/highL_lensedCL_template(lmax_computed_cl,1)
        do l = lmax_computed_cl+1, lmax
            Theory%cl(l,1:num_ClsS) =  highL_norm*highL_lensedCL_template(l,1:num_clsS)
        end do
    end if

    end subroutine SetPowersFromCAMB


    subroutine SetPkFromCAMB(Theory,M)
    use camb, only : MatterTransferData
    Type(TheoryPredictions) Theory
    Type(MatterTransferData) M

    call Theory_GetMatterPowerData(M,Theory,1)

    end subroutine SetPkFromCAMB

    subroutine InitCAMB(CMB,error, DoReion)
    type(CMBParams), intent(in) :: CMB
    logical, optional, intent(in) :: DoReion
    logical WantReion
    type(CAMBParams)  P
    integer error

    if (present(DoReion)) then
        WantReion = DoReion
    else
        WantReion = .true.
    end if

    call CMBToCAMB(CMB, P)
    call CAMBParams_Set(P,error,WantReion)

    end subroutine InitCAMB

    function GetOpticalDepth(CMB)
    type(CMBParams) CMB
    real(mcp) GetOpticalDepth
    type(CAMBParams)  P
    integer error

    call CMBToCAMB(CMB, P)
    call CAMBParams_Set(P,error)

    if (error/= 0) then
        GetOpticalDepth = -1
    else
        GetOpticalDepth = Reionization_GetOptDepth(P%Reion, P%ReionHist)
    end if
    end  function GetOpticalDepth

    function GetZreFromTau(CMB, tau)
    type(CMBParams) CMB
    real(mcp), intent(in) :: tau
    real(mcp) GetZreFromTau
    type(CAMBParams)  P

    call CMBToCAMB(CMB, P)
    GetZreFromTau = CAMB_GetZreFromTau(P,dble(tau))

    end  function GetZreFromTau

    function CMBToTheta(CMB)
    use ModelParams
    implicit none
    Type(CMBParams) CMB
    real(mcp) CMBToTheta
    integer error

    call InitCAMB(CMB,error,.false.)
    CMBToTheta = CosmomcTheta()

    end function CMBToTheta




    subroutine InitCAMBParams(P)
    use lensing
    use ModelParams
    use mpk
    type(CAMBParams)  P
    integer zix
    !JD Changed P%Transfer%redshifts and P%Transfer%num_redshifts to 
    !P%Transfer%PK_redshifts and P%Transfer%PK_num_redshifts respectively
    !for nonlinear lensing of CMB + LSS compatibility
    Threadnum =num_threads
    w_lam = -1
    wa_ppf = 0._dl
    call CAMB_SetDefParams(P)

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
        P%Transfer%kmax = max(1.2,power_kmax)
    else
        P%Transfer%kmax = max(0.8,power_kmax)
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

    end subroutine InitCAMBParams

    !Mapping between array of power spectrum parameters and CAMB
    subroutine SetCAMBInitPower(P,CMB,in)
    use camb
    use settings
    use cmbtypes
    implicit none
    type(CAMBParams)  P
    Type(CMBParams) CMB

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
        lens_recon_scale = CMB%InitPower(Aphiphi_index)
    else
        stop 'params_CMB:Wrong initial power spectrum'
    end if

    end subroutine SetCAMBInitPower


    subroutine LoadFiducialHighLTemplate
    !This should be a lensed scalar CMB power spectrum, e.g. for including at very high L where foregrounds etc. dominate anyway
    integer L
    real(mcp) array(4), nm
    character(LEN=Ini_max_string_len) :: fname

    fname = ReadIniFilename(DefIni,'highL_theory_cl_template',DataDir,.true.)
    allocate(highL_lensedCL_template(2:lmax, num_clsS))
    call OpenTxtFile(fname,tmp_file_unit)
    do
        read(tmp_file_unit,*, end=500) L , array
        if (L>lmax) exit
        nm = 2*pi/(l*(l+1))
        if (L>=2) highL_lensedCL_template(L,1:num_clsS) = nm*array(TensClOrder(1:num_clsS))
    end do
500 close(tmp_file_unit)

    if (highL_lensedCL_template(2,1) < 100) &
    call MpiStop('highL_theory_cl_template must be in muK^2')

    if (L<lmax) call MpiStop('highL_theory_cl_template does not go to lmax')
    if (num_cls_ext>0 .and. MpiRank==0) &
    write(*,*) 'warning: zero padding ext cls in LoadFiducialHighLTemplate'

    end subroutine LoadFiducialHighLTemplate

    subroutine ParamSetInfo_Clear(Info)
    Type(ParamSetInfo) Info

    call CAMB_InitCAMBdata(Info%Transfers)
    Info%validInfo = .false.

    end subroutine

    subroutine CMB_Initialize(Info)
    Type(ParamSetInfo) Info
    type(CAMBParams)  P
    compute_tensors = Ini_Read_Logical('compute_tensors',.false.)
    if (num_cls==3 .and. compute_tensors) write (*,*) 'WARNING: computing tensors with num_cls=3 (BB=0)'
    CMB_lensing = Ini_Read_Logical('CMB_lensing',CMB_lensing)
    use_lensing_potential = Ini_Read_logical('use_lensing_potential',use_lensing_potential)
    use_nonlinear_lensing = Ini_Read_logical('use_nonlinear_lensing',use_nonlinear_lensing)
    if (CMB_lensing) num_clsS = num_cls   !Also scalar B in this case
    if (use_lensing_potential .and. num_cls_ext ==0) &
    call MpiStop('num_cls_ext should be > 0 to use_lensing_potential')
    if (use_lensing_potential .and. .not. CMB_lensing) &
    call MpiStop('use_lensing_potential must have CMB_lensing=T')
    lmax_computed_cl = Ini_Read_Int('lmax_computed_cl',lmax)
    if (lmax_computed_cl /= lmax) then
        if (lmax_tensor > lmax_computed_cl) call MpiStop('lmax_tensor > lmax_computed_cl')
        call LoadFiducialHighLTemplate
    end if

    call InitCAMBParams(P)

    if (Feedback > 0 .and. MPIRank==0) then
        write (*,*) 'Computing tensors:', compute_tensors
        write (*,*) 'Doing CMB lensing:',CMB_lensing
        write (*,*) 'Doing non-linear Pk:',P%NonLinear

        write(*,'(" lmax              = ",1I4)') lmax
        write(*,'(" lmax_computed_cl  = ",1I4)') lmax_computed_cl
        write(*,*) 'max_eta_k         = ', P%Max_eta_k
        write(*,*) 'transfer kmax     = ', P%Transfer%kmax

        if (compute_tensors) write(*,'(" lmax_tensor    = ",1I4)') lmax_tensor
        write(*,'(" Number of C_ls = ",1I4)') num_cls
    end if

    call ParamSetInfo_Clear(Info)

    P%WantTensors = compute_tensors
    CAMBP = P

    end subroutine CMB_Initialize


    subroutine AcceptReject(accpt, CurParams, Trial)
    logical, intent(in) :: accpt
    Type(ParamSetInfo) CurParams, Trial

    if (.not. associated(CurParams%Transfers%ClTransScal%Delta_p_l_k,&
    Trial%Transfers%ClTransScal%Delta_p_l_k)) then
        !If they point to same memory don't need to free anything
        if (accpt) then
            call CAMB_FreeCAMBdata(CurParams%Transfers)
        else
            call CAMB_FreeCAMBdata(Trial%Transfers)
        end if
    end if

    end subroutine AcceptReject

    end module CMB_Cls

