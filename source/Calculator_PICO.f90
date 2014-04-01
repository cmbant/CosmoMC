    !Use CAMB
    module Calculator_PICO
    use CosmologyTypes
    use CosmoTheory
    use Calculator_CAMB
    use CAMB, only : CAMB_GetResults, CAMB_GetAge, CAMBParams, CAMB_SetDefParams, &
    AccuracyBoost,  Cl_scalar, Cl_tensor, Cl_lensed, outNone, w_lam, wa_ppf,&
    CAMBParams_Set, MT, CAMBdata, NonLinear_Pk, Nonlinear_lens, Reionization_GetOptDepth, CAMB_GetZreFromTau, &
    CAMB_GetTransfers,CAMB_FreeCAMBdata,CAMB_InitCAMBdata, CAMB_TransfersToPowers, Transfer_SetForNonlinearLensing, &
    initial_adiabatic,initial_vector,initial_iso_baryon,initial_iso_CDM, initial_iso_neutrino, initial_iso_neutrino_vel, &
    HighAccuracyDefault, highL_unlensed_cl_template, ThermoDerivedParams, nthermo_derived, BackgroundOutputs, &
    Transfer_SortAndIndexRedshifts, & !JD added for nonlinear lensing of CMB + MPK compatibility
    Recombination_Name, reionization_name, power_name, threadnum, version, lmin
    use Errors !CAMB
    use settings
    use likelihood
    use Calculator_Cosmology
    use GeneralTypes
    implicit none
    private

    Type, extends(CAMB_Calculator) :: PICO_Calculator
    contains
    procedure :: CMBToCAMB => PICO_CMBToCAMB
    procedure :: GetNewTransferData => PICO_GetNewTransferData
    procedure :: GetNewPowerData => PICO_GetNewPowerData
    procedure :: GetCosmoTheoryForImportance => PICO_GetTheoryForImportance
    procedure :: VersionTraceOutput => PICO_VersionTraceOutput
    procedure :: ReadParams => PICO_ReadParams
    end type PICO_Calculator

    public PICO_Calculator
    contains

    subroutine PICO_CMBToCAMB(this,CMB,P)
    class(PICO_Calculator) :: this
    class(CMBParams) CMB
    type(CAMBParams)  P

    call this%CAMB_Calculator%CMBToCamb(CMB,P)
    P%Reion%use_optical_depth = .true.
    P%OnlyTransfers = .true.
    if (.not. CosmoSettings%CMB_Lensing) call MpiStop('PICO assumes lensed CMB')

    end subroutine PICO_CMBToCAMB

    subroutine PICO_GetNewTransferData(this, CMB,Info,Theory,error)
    class(PICO_Calculator) :: this
    class(CMBParams) CMB
    class(TTheoryIntermediateCache), pointer :: Info
    class(TCosmoTheoryPredictions) :: Theory
    integer error
    !Do nothing
    end subroutine PICO_GetNewTransferData

    subroutine PICO_GetOutputArray(Theory,i,j, name)
    class(TCosmoTheoryPredictions) :: Theory
    integer i,j
    character(LEN=*) :: name

    if (CosmoSettings%cl_lmax(i,j)>0) &
    & call fpico_read_output(name,Theory%Cls(i,j)%CL(lmin:),lmin,min(CosmoSettings%lmax_computed_cl,CosmoSettings%cl_lmax(i,j)))

    end subroutine PICO_GetOutputArray

    subroutine PICO_AddTensorOutputArray(Theory,i,j, name)
    class(TCosmoTheoryPredictions) :: Theory
    integer i,j, lmx
    character(LEN=*) :: name
    real(mcp) tmp(lmin:CosmoSettings%lmax_tensor)

    lmx=min(CosmoSettings%lmax_tensor,CosmoSettings%cl_lmax(i,j))
    if (lmx>0) then
        call fpico_read_output(name,tmp(lmin:),lmin,lmx)
        Theory%Cls(i,j)%CL(lmin:lmx) = Theory%Cls(i,j)%CL(lmin:lmx) + tmp(lmin:lmx)
    end if

    end subroutine PICO_AddTensorOutputArray

    subroutine PICO_GetNewPowerData(this, CMB, Info, Theory, error)
    use ModelParams
    use ModelData
    class(PICO_Calculator) :: this
    class(CMBParams) :: CMB
    class(TTheoryIntermediateCache), pointer :: Info
    class(TCosmoTheoryPredictions) :: Theory
    integer error,i,j, lmaxCL,lmx
    type(CAMBParams) P
    logical :: success
    real(mcp) :: silviamassive
    real(mcp) :: highL_norm = 0
    integer, save, allocatable :: indicesT(:,:)

    if (.not. allocated(indicesT)) then
        allocate(indicesT(3,3))
        indicesT=0
        indicesT(1,1) =  C_Temp
        indicesT(2,2) =  C_E
        indicesT(2,1) =  CT_Cross
        indicesT(3,3) =  CT_B
    end if

    call this%CMBToCAMB(CMB, P)
    call CAMBParams_Set(P)
    call this%SetBackgroundTheoryData(CMB,Theory,error)

    call fpico_reset_params()
    call fpico_set_param("ombh2", CMB%ombh2)
    call fpico_set_param("omch2", CMB%omch2)
    call fpico_set_param("omnuh2", CMB%omnuh2)
    call fpico_set_param("omvh2", CMB%omv*CMB%h**2)
    call fpico_set_param("omk", CMB%omk)
    call fpico_set_param("hubble", CMB%H0)
    call fpico_set_param("Alens", CMB%Alens)
    call fpico_set_param("w", CMB%w)
    call fpico_set_param("theta", CosmomcTheta())
    call fpico_set_param("helium_fraction", p%yhe)

    call fpico_set_param("massless_neutrinos", p%Num_Nu_massless)
    !!!Check what's going on with neutrinos
    call fpico_set_param("massive_neutrinos", p%Num_Nu_massless+p%Num_Nu_massive)
    call fpico_set_param("scalar_spectral_index(1)",p%InitPower%an(1))
    call fpico_set_param("tensor_spectral_index(1)",p%InitPower%ant(1))
    call fpico_set_param("scalar_nrun(1)",p%InitPower%n_run(1))
    call fpico_set_param("initial_ratio(1)",p%InitPower%rat(1))
    call fpico_set_param("scalar_amp(1)",p%InitPower%ScalarPowerAmp(1))
    call fpico_set_param("pivot_scalar",p%InitPower%k_0_scalar)
    call fpico_set_param("re_optical_depth",p%Reion%optical_depth)

    call fpico_reset_requested_outputs()
    if (P%WantCls) then
        call fpico_request_output("lensed_TT")
        call fpico_request_output("lensed_TE")
        call fpico_request_output("lensed_EE")
        call fpico_request_output("lensed_BB")
        if (P%WantTensors) then
            call fpico_request_output("tensor_TT")
            call fpico_request_output("tensor_TE")
            call fpico_request_output("tensor_EE")
            call fpico_request_output("tensor_BB")
        end if
        if (P%WantTransfer) then
            call fpico_request_output("k")
            call fpico_request_output("pk")
        end if
    end if

    call fpico_compute_result(success)
    if (.not. success) call mpiStop('PICO failed to get result')
    call PICO_GetOutputArray(Theory,1,1, "lensed_TT")
    call PICO_GetOutputArray(Theory,2,1, "lensed_TE")
    call PICO_GetOutputArray(Theory,2,2, "lensed_EE")
    call PICO_GetOutputArray(Theory,3,3, "lensed_BB")

    do i=1, min(3,CosmoSettings%num_cls)
        do j= i, 1, -1
            lmaxCL = CosmoSettings%cl_lmax(i,j)
            lmx = min(CosmoSettings%lmax_computed_cl, lmaxCL)
            if (lmx/=0) then
                associate( CL => Theory%Cls(i,j)%CL)
                    if (CosmoSettings%lmax_computed_cl < lmaxCL) then
                        if (highL_norm ==0) & !normally normalize off TT
                        & highL_norm = CL(lmx)/this%highL_lensedCL_template(lmx,indicesT(i,j))
                        CL(lmx+1:lmaxCL) =  highL_norm*this%highL_lensedCL_template(lmx+1:lmaxCL,indicesT(i,j))
                    end if
                    if (CosmoSettings%compute_tensors) then
                        call PICO_AddTensorOutputArray(Theory,1,1, "tensor_TT")
                        call PICO_AddTensorOutputArray(Theory,2,1, "tensor_TE")
                        call PICO_AddTensorOutputArray(Theory,2,2, "tensor_EE")
                        call PICO_AddTensorOutputArray(Theory,3,3, "tensor_BB")
                    end if
                    end associate
            end if
        end do
    end do

    !redshifts are in increasing order, so last index is redshift zero
    if (CosmoSettings%Use_LSS .or. CosmoSettings%get_sigma8) then
        call MpiStop('No MPK or sigma8 with PICO')
    end if

    end subroutine PICO_GetNewPowerData

    subroutine PICO_GetTheoryForImportance(this, CMB, Theory, error)
    class(PICO_Calculator) :: this
    class(CMBParams) CMB
    class(TCosmoTheoryPredictions) Theory
    integer :: error

    call MpiStop('PICO: no GetTheoryForImportance yet')
    end subroutine PICO_GetTheoryForImportance

    subroutine PICO_VersionTraceOutput(this, ReadValues)
    use GaugeInterface, only : Eqns_name
    class(PICO_Calculator) :: this
    class(TNameValueList) :: ReadValues

    call ReadValues%Add('Compiled_PICO_version', '3.2') !dynamic version variable?
    call this%CAMB_Calculator%VersionTraceOutput(ReadValues)

    end subroutine PICO_VersionTraceOutput


    subroutine PICO_ReadParams(this,Ini)
    class(PICO_Calculator) :: this
    class(TSettingIni) :: Ini

    call this%CAMB_Calculator%ReadParams(Ini)
    this%calcName ='PICO'

    call fpico_load(Ini%Read_String("pico_datafile"))
    call fpico_set_verbose(Ini%Read_Logical("pico_verbose",.false.))

    end subroutine PICO_ReadParams


    end module Calculator_PICO
