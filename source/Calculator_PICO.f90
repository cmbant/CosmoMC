#ifdef PICO

    module Calculator_PICO
    use CosmologyTypes
    use CosmoTheory
    use Calculator_CAMB
    use CAMB, only : CAMBParams_Set, highL_unlensed_cl_template, ThermoDerivedParams, nthermo_derived, &
        BackgroundOutputs, lmin, CAMBParams
    use Errors !CAMB
    use settings
    use likelihood
    use Calculator_Cosmology
    use GeneralTypes
    use fpico
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
    integer(fpint) lmax, lmn

    if (CosmoSettings%cl_lmax(i,j)>0) then
        lmax = min(CosmoSettings%lmax_computed_cl,CosmoSettings%cl_lmax(i,j))
        lmn = lmin
        call fpico_read_output(name,Theory%Cls(i,j)%CL(lmin:),lmn,lmax)
    end if

    end subroutine PICO_GetOutputArray

    subroutine PICO_GetNewPowerData(this, CMB, Info, Theory, error)
    use ModelParams
    use ModelData
    class(PICO_Calculator) :: this
    class(CMBParams) :: CMB
    class(TTheoryIntermediateCache), pointer :: Info
    class(TCosmoTheoryPredictions) :: Theory
    integer error,i,j, lmaxCL,lmx
    type(CAMBParams) P
    integer(fpint) :: success
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

    if (p%InitPower%n_runrun(1)/=0 .or. p%InitPower%nt_run(1)/=0 .or. p%InitPower%ant(1)/=0) &
        & call MpiStop('PICO: currently unsupported initial power parameter')

    if (CosmoSettings%Use_LSS .or. CosmoSettings%get_sigma8) then
        call MpiStop('PICO: currently no MPK or sigma8')
    end if

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

    call fpico_set_param("massive_neutrinos", CMB%nnu)
    call fpico_set_param("scalar_spectral_index(1)",p%InitPower%an(1))
    call fpico_set_param("scalar_nrun(1)",p%InitPower%n_run(1))
    call fpico_set_param("initial_ratio(1)",p%InitPower%rat(1))
    call fpico_set_param("scalar_amp(1)",p%InitPower%ScalarPowerAmp(1))
    call fpico_set_param("pivot_scalar",p%InitPower%k_0_scalar)
    call fpico_set_param("pivot_tensor",p%InitPower%k_0_tensor)
    call fpico_set_param("re_optical_depth",CMB%tau)
    !    call fpico_set_param("force",1.d0)

    call fpico_reset_requested_outputs()
    if (P%WantCls) then
        call fpico_request_output("cl_TT")
        if (CosmoSettings%num_cls>1) then
            call fpico_request_output("cl_TE")
            call fpico_request_output("cl_EE")
            if (CosmoSettings%num_cls>2) then
                if (CosmoSettings%cl_lmax(3,3)>0) call fpico_request_output("cl_BB")
                if (CosmoSettings%num_cls>3) then
                    !                    call MpiStop('PICO: lensing output currently innaccurate')
                    !                    if (CosmoSettings%cl_lmax(4,4)>0) call fpico_request_output("cl_pp")
                end if
            end if
        end if
    end if
   ! if (P%WantTransfer) then
   !     call fpico_request_output("k")
   !     call fpico_request_output("pk")
   ! end if

    call fpico_compute_result(success)
    if (success == 0) then
        ! for now just reject if out of pico bounds
        !call mpiStop('PICO failed to get result')
        if (Feedback>1) write(*,*) 'PICO out of bounds'
        error = 1
        return
    end if
    error = 0
    if (P%WantCls) then
        call PICO_GetOutputArray(Theory,1,1, "cl_TT")
        if (CosmoSettings%num_cls>1) then
            call PICO_GetOutputArray(Theory,2,1, "cl_TE")
            call PICO_GetOutputArray(Theory,2,2, "cl_EE")
            if (CosmoSettings%num_cls>2) call PICO_GetOutputArray(Theory,3,3, "cl_BB")
            !            if (CosmoSettings%num_cls>3) call PICO_GetOutputArray(Theory,4,4, "cl_pp") need to change units, as in scalCls
        end if
    end if

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
                end associate
            end if
        end do
    end do

    end subroutine PICO_GetNewPowerData

    subroutine PICO_GetTheoryForImportance(this, CMB, Theory, error)
    class(PICO_Calculator) :: this
    class(CMBParams) CMB
    class(TCosmoTheoryPredictions) Theory
    integer :: error

    call this%GetNewPowerData(CMB, null(), Theory, error)
    call this%GetNewBackgroundData(CMB,Theory,error)
    if (error==0) call this%SetDerived(Theory)

    end subroutine PICO_GetTheoryForImportance

    subroutine PICO_VersionTraceOutput(this, ReadValues)
    use GaugeInterface, only : Eqns_name
    class(PICO_Calculator) :: this
    class(TNameValueList) :: ReadValues

    call ReadValues%Add('Compiled_PICO_version', '3.3.0') !dynamic version variable?
    call this%CAMB_Calculator%VersionTraceOutput(ReadValues)

    end subroutine PICO_VersionTraceOutput


    subroutine PICO_ReadParams(this,Ini)
    class(PICO_Calculator) :: this
    class(TSettingIni) :: Ini

    call this%CAMB_Calculator%ReadParams(Ini)
    this%calcName ='PICO'
    
    !$ write(*,*) '**WARNING**: pico may not work when CosmoMC compiled with -openmp (why??)'

    call fpico_init(1_fpint)
    call fpico_load(Ini%Read_String_Default("pico_datafile", EnvDefault=.true.))
    call fpico_set_verbose(int(IfThenElse(Ini%Read_Logical("pico_verbose",.false.),1,0),fpint))

    end subroutine PICO_ReadParams


    end module Calculator_PICO
#endif
