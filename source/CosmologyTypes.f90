    module CosmologyTypes
    use settings
    use likelihood
    use GeneralTypes
    use ObjectLists
    implicit none

    integer, parameter :: derived_age=1, derived_zstar=2, derived_rstar=3, derived_thetastar=4, derived_DAstar = 5, &
        derived_zdrag=6, derived_rdrag=7,derived_kD=8,derived_thetaD=9, derived_zEQ =10, derived_keq =11, &
        derived_thetaEQ=12, derived_theta_rs_EQ = 13 !index in derived parameters array

    integer, parameter :: As_index=1, ns_index =2, nrun_index=3, nrunrun_index=4, amp_ratio_index = 5, &
        & nt_index= 6, ntrun_index = 7, Aphiphi_index = 8, last_power_index = Aphiphi_index

    integer, parameter :: max_inipower_params = 10

    real(mcp), parameter :: cl_norm = 1e-10_mcp !units for As
    integer, parameter :: max_derived_parameters = 30

    integer :: num_hard, num_initpower
    integer :: index_initpower

    character(LEN=4), parameter :: CMB_CL_Fields = 'TEBP'

    integer, parameter :: neutrino_hierarchy_normal = 1, neutrino_hierarchy_inverted = 2, neutrino_hierarchy_degenerate = 3
    character(LEN=Ini_Enumeration_Len), parameter :: neutrino_types(3) = &
        [character(Ini_Enumeration_Len)::'normal','inverted','degenerate']

    Type TCosmoTheoryParams
        logical :: get_sigma8 = .true.
        logical :: Use_LSS = .false.
        logical :: Use_CMB = .false.
        logical :: use_nonlinear = .false.    !JD for WiggleZ MPK

        !l_max. Tensors are not computed unless compute_tensors = T in input file
        !Make these multiples of 50, should be 50 more than you need accurately
        integer :: lmax = 0
        integer :: num_cls = 0
        integer :: lmax_tensor = 600
        !note only lmax_computed_cl is actually calculated
        integer :: lmax_computed_cl = 4500 !Never compute primordial above this, just fit
        integer :: lmin_computed_cl = 2000 !if cl needed, calculate to at least this
        integer :: lmin_store_all_cmb = 0 !>0 if you want output everything even if not used

        !redshifts for output of BAO background parameters
        real(mcp) :: z_outputs(1) = [0.57_mcp]

        logical :: CMB_lensing = .true.
        logical :: use_lensing_potential = .false.
        logical :: use_nonlinear_lensing = .false.
        logical :: compute_tensors = .false.

        !Parameters for calculating/storing the matter power spectrum
        logical :: use_matterpower = .false.
        logical :: use_Weylpower = .false. !power spectrum of Weyl potential for lensing
        logical :: use_sigmaR =.false. !sigma_R, e.g. for clusters
        real(mcp) :: power_kmax = 0.8_mcp
        integer :: num_power_redshifts = 0

        !Only used in params_CMB
        real(mcp) :: pivot_k = 0.05_mcp !Point for defining primordial power spectra
        real(mcp) :: tensor_pivot_k = 0.05_mcp !Point for defining tensor power spectra
        logical :: inflation_consistency = .true. !fix n_T or not

        logical :: bbn_consistency = .true. !JH

        integer :: num_massive_neutrinos = -1 !if neutrino_hierarcy_degenerate, number of massive degenerate eigenstates
        integer :: neutrino_hierarchy = neutrino_hierarchy_normal

    end Type TCosmoTheoryParams

    integer, parameter ::  CL_T = 1, CL_E=2, CL_B=3, CL_Phi=4

    Type, extends(TCosmoTheoryParams):: TCosmoTheorySettings
        !Just add the allocatable components
        integer, allocatable :: cl_lmax(:,:)
        integer, allocatable :: ArraySizes(:)
        !e.g. lmax_cl(1,1) is lmax for TT; zero if CL is not used; order is T, E, B, Phi
        real(mcp), dimension(:), allocatable :: power_redshifts
    contains
    procedure, private :: Initialize_PKSettings
    procedure, private :: Initialize_CMBSettings
    procedure :: InitForLikelihoods => TCosmoTheorySettings_InitForLikelihoods
    procedure :: ReadParams => TCosmoTheorySettings_ReadParams
    end type TCosmoTheorySettings

    type(TCosmoTheorySettings), save, target :: CosmoSettings

    type, extends(TDatasetFileLikelihood) :: TCosmologyRequirementsLikelihood
        logical :: needs_background_functions = .true.
        logical :: needs_powerspectra = .false.

        integer, allocatable :: cl_lmax(:,:)!T, E, B, Phi, with first index higher or equal, e.g. ET not TE

        logical :: needs_nonlinear_pk = .false.
        logical :: needs_exact_z = .false.
        logical :: needs_Weylpower = .false.
        logical :: needs_sigmaR = .false.
        real(mcp) :: kmax = 0.8_mcp
        integer :: num_z = 0
        real(mcp) :: max_z = 0._mcp
        real(mcp), dimension(:), allocatable :: exact_z
        integer, dimension(:), allocatable :: exact_z_index
    contains
    procedure :: InitForSettings => TCosmologyRequirementsLikelihood_InitForSettings
    end type TCosmologyRequirementsLikelihood

    Type, extends(TTheoryParams) :: CMBParams
        real(mcp) InitPower(max_inipower_params)
        !These are fast paramters for the initial power spectrum
        !Now remaining (non-independent) parameters
        real(mcp) omb, omc, omv, omnu, omk, omdm
        real(mcp) ombh2, omch2, omnuh2, omdmh2
        real(mcp) zre, zre_delta, nufrac
        real(mcp) h, H0, tau
        real(mcp) w, wa
        real(mcp) YHe, nnu, iso_cdm_correlated, ALens, Alensf, fdm !fdm is dark matter annihilation, eg,. 0910.3663
        real(mcp) :: omnuh2_sterile = 0._mcp  !note omnhu2 is the sum of this + standard neutrinos
        real(mcp) :: sum_mnu_standard
        real(mcp) reserved(5)
    end Type CMBParams

    Type, extends(TParameterization) :: TCosmologyParameterization
        logical :: late_time_only = .false.
    contains
    procedure :: NewTheoryParams => TCosmologyParameterization_NewTheoryParams
    procedure :: SetTheoryParameterNumbers => TCosmologyParameterization_SetNumbers
    end type

    contains

    subroutine TCosmologyParameterization_SetNumbers(this, slow_num, semi_slow_num)
    class(TCosmologyParameterization) :: this
    integer, intent (in) :: slow_num, semi_slow_num
    integer i
    class(TDataLikelihood), pointer :: DataLike

    !Called after this%init
    num_hard = slow_num
    num_initpower = semi_slow_num
    if (num_hard + num_initpower /= num_theory_params) &
        call MpiStop('SetTheoryParameterNumbers: parameter numbers do not match')
    index_initpower = num_hard+1
    index_semislow = index_initpower
    if (num_initpower> max_inipower_params) call MpiStop('see CosmologyTypes.f90: num_initpower> max_inipower_params')

    do i=1,DataLikelihoods%Count
        DataLike=>DataLikelihoods%Item(i)
        select type (DataLike)
        class is (TCosmologyRequirementsLikelihood)
            if (DataLike%needs_background_functions) DataLike%dependent_params(1:num_hard)=.true.
            if (DataLike%needs_powerspectra) DataLike%dependent_params(1:num_theory_params)=.true.
        end select
    end do

    end subroutine TCosmologyParameterization_SetNumbers

    subroutine TCosmologyParameterization_NewTheoryParams(this,TheoryParams)
    class(TCosmologyParameterization) :: this
    class(TTheoryParams), allocatable :: TheoryParams

    allocate(CMBParams::TheoryParams)

    end subroutine TCosmologyParameterization_NewTheoryParams


    subroutine TCosmoTheorySettings_ReadParams(this, Ini)
    class(TCosmoTheorySettings) this
    class(TSettingIni) :: Ini

    this%compute_tensors = Ini%Read_Logical('compute_tensors',.false.)

    call Ini%Read('CMB_lensing',this%CMB_lensing)

    if (this%CMB_lensing) call Ini%Read('use_nonlinear_lensing',this%use_nonlinear_lensing)
    if (Ini%HasKey('use_lensing_potential')) &
        & write(*,*) 'NOTE: use_lensing_potential now set internally from likelihoods'
    if (Ini%HasKey('use_CMB')) &
        & write(*,*) 'NOTE: use_CMB now set internally from likelihoods'

    call Ini%Read('pivot_k',this%pivot_k)
    this%tensor_pivot_k = this%pivot_k
    call Ini%Read('tensor_pivot_k',this%tensor_pivot_k)

    call Ini%Read('inflation_consistency',this%inflation_consistency)
    call Ini%Read('bbn_consistency',this%bbn_consistency)

    this%neutrino_hierarchy = Ini%Read_Enumeration('neutrino_hierarchy',neutrino_types, neutrino_hierarchy_normal)
    if (this%neutrino_hierarchy == neutrino_hierarchy_degenerate) then
        call Ini%Read('num_massive_neutrinos',this%num_massive_neutrinos)
    else if (Ini%Read_Int('num_massive_neutrinos',0)>0) then
        write(*,*) 'NOTE: num_massive_neutrinos ignored, using specified hierarchy'
    end if
    call Ini%Read('lmax_computed_cl',this%lmax_computed_cl)
    call Ini%Read('lmin_computed_cl',this%lmin_computed_cl)
    call Ini%Read('lmin_store_all_cmb',this%lmin_store_all_cmb)
    call Ini%Read('lmax_tensor',this%lmax_tensor)

    end subroutine TCosmoTheorySettings_ReadParams


    subroutine TCosmoTheorySettings_InitForLikelihoods(this)
    class(TCosmoTheorySettings) this
    class(TDataLikelihood), pointer :: DataLike
    integer i,j

    call this%Initialize_PKSettings()

    call this%Initialize_CMBSettings()

    if (this%use_lensing_potential .and. .not. this%CMB_lensing) &
        & call MpiStop('use_lensing_potential must have CMB_lensing=T')

    if (Feedback > 0 .and. MPIRank==0) then
        write (*,*) 'Doing non-linear Pk:', this%use_nonlinear

        if(this%use_CMB)then
            write (*,*) 'Doing CMB lensing:', this%CMB_lensing
            if (this%CMB_lensing) write (*,*) 'Doing non-linear lensing:', this%use_nonlinear_lensing
            if (allocated(this%cl_lmax)) then
                do i=1, this%num_cls
                    do j= i, 1, -1
                        if (this%cl_lmax(i,j) >0) &
                            write(*,'(" '//CMB_CL_Fields(i:i)//CMB_CL_Fields(j:j)//' lmax = ",(I5))') this%cl_lmax(i,j)
                    end do
                end do
                write(*,'(" lmax_computed_cl  = ",1I5)') this%lmax_computed_cl
                write (*,*) 'Computing tensors:', this%compute_tensors
                if (this%compute_tensors) write(*,'(" lmax_tensor    = ",1I5)') this%lmax_tensor
            end if
        end if
    end if

    do i=1,DataLikelihoods%Count
        DataLike=>DataLikelihoods%Item(i)
        select type (DataLike)
        class is (TCosmologyRequirementsLikelihood)
            call DataLike%InitForSettings(this)
        end select
    end do

    end subroutine TCosmoTheorySettings_InitForLikelihoods


    subroutine Initialize_CMBSettings(this)
    class(TCosmoTheorySettings) this
    class(TDataLikelihood), pointer :: DataLike
    integer i, parse, numcls, a

    numcls=0
    do parse=1,2
        do i=1,DataLikelihoods%Count
            DataLike=>DataLikelihoods%Item(i)
            select type (DataLike)
            class is (TCosmologyRequirementsLikelihood)
                if (allocated(DataLike%cl_lmax)) then
                    if (any(DataLike%cl_lmax>0)) then
                        DataLike%needs_powerspectra = .true.
                    else
                        cycle
                    end if
                    if (parse==1) then
                        numcls = max(numcls,size(DataLike%cl_lmax,2))
                    else
                        if (size(DataLike%cl_lmax,2)/=size(DataLike%cl_lmax,1)) &
                            & call MpiStop('cl_max(i,j) should be square: '//trim(DataLike%Name))
                        do a=1, size(DataLike%cl_lmax,2)
                            if (any(DataLike%cl_lmax(a,a+1:)>0)) &
                                & call MpiStop('cl_max(i,j) should have i>=j: '//trim(DataLike%Name))
                            this%cl_lmax(a,1:a) = max(this%cl_lmax(a,1:a), DataLike%cl_lmax(a,1:a))
                        end do
                    end if
                else
                    if(DataLike%LikelihoodType=='CMB') call MpiStop(DataLike%name//' CMB likelihood seems to have no cl_lmax set')
                end if
            end select
        end do
        if (parse==2) exit
        if (this%lmin_store_all_cmb>0 ) numcls=max(numcls,IfThenElse(this%use_nonlinear_lensing,4,2))
        allocate(this%cl_lmax(numcls,numcls), source=0)
        if (this%lmin_store_all_cmb>0) then
            this%cl_lmax(1,1) = this%lmin_store_all_cmb
            this%cl_lmax(2,1:2) = this%lmin_store_all_cmb
            if (this%use_nonlinear_lensing) then
                !Only force BB and PhiPhi output if likely to be accurate
                this%cl_lmax(3,3) = this%lmin_store_all_cmb
                this%cl_lmax(4,4) = this%lmin_store_all_cmb
            end if
        end if
    end do
    where (this%cl_lmax>0)
        this%cl_lmax = max(this%cl_lmax,this%lmin_computed_cl)
    end where
    this%num_cls = numcls
    this%use_CMB = this%num_cls > 0
    this%use_lensing_potential = numcls>3
    if (this%use_lensing_potential) this%use_lensing_potential= any(this%cl_lmax(CL_Phi,:)/=0)
    this%lmax = maxval(this%cl_lmax)
    this%lmax_computed_cl = min(this%lmax,this%lmax_computed_cl)

    end subroutine Initialize_CMBSettings

    subroutine Initialize_PKSettings(this)
    class(TCosmoTheorySettings) this
    class(TDataLikelihood), pointer :: DataLike
    Type(TRealList) :: exact_z, full_z
    real(mcp) :: dlnz, maxz, zcur
    integer :: i,iz,izprev
    integer :: num_range

    maxz = 0
    dlnz = 30

    call full_z%Add(0.d0)
    this%use_LSS = size(CosmoSettings%z_outputs)>0 .and. this%get_sigma8 !e.g. for growth function

    do i=1,DataLikelihoods%Count
        DataLike=>DataLikelihoods%Item(i)
        select type (DataLike)
        class is (TCosmologyRequirementsLikelihood)
            if (DataLike%needs_powerspectra) then
                if (DataLike%needs_exact_z .or. DataLike%num_z>0 .or. DataLike%needs_sigmaR) then
                    this%Use_LSS = .true.
                else
                    cycle
                end if
                this%power_kmax = max(this%power_kmax,DataLike%kmax)
                this%use_nonlinear = this%use_nonlinear .or. DataLike%needs_nonlinear_pk
                this%use_matterpower = .true.
                this%use_Weylpower = this%use_Weylpower .or. DataLike%needs_Weylpower
                this%use_sigmaR = this%use_sigmaR .or. DataLike%needs_sigmaR
                if(DataLike%needs_exact_z) then
                    call exact_z%AddArrayItems(DataLike%exact_z)
                else
                    num_range = DataLike%num_z
                    maxz = max(maxz,DataLike%max_z)
                    if(num_range >2 .and. maxz>0) then
                        dlnz = min(dlnz,log(DataLike%max_z+1)/(num_range-1))
                    else if(num_range<2 .and. maxz > 0)then
                        write(*,'("ERROR: ",A," dataset: ",A, "wants less than 2 redshifts")')&
                            trim(DataLike%LikelihoodType),trim(DataLike%name)
                        write(*,'("       but wants a maximum redshift of ",F7.2,". A minimum ")')maxz
                        write(*,*)"       of 2 redshifts is required or PowerAtZ will fail."
                        write(*,*)"       Check dataset settings!"
                        call Mpistop()
                    else if(num_range>1 .and. maxz==0.)then
                        write(*,'("ERROR: ",A," dataset: ",A, "wants only ",I0," redshifts")')&
                            trim(DataLike%LikelihoodType),trim(DataLike%name),num_range
                        write(*,*)"       but wants a maximum redshift of 0.0.  Check dataset settings!"
                        call Mpistop()
                    else
                        cycle
                    end if
                end if
            end if
        end select
    end do

    if(.not. this%use_LSS) return

    call exact_z%AddArrayItems(CosmoSettings%z_outputs)

    !Build array of redshifts where the redshift exact value doesn't matter
    if(maxz>0)then
        num_range = ceiling(log(maxz+1)/dlnz)
        dlnz = log(maxz+1)/(num_range)
        do i=1,num_range-1
            zcur = dexp(dlnz*i)-1
            call full_z%Add(zcur)
        end do
        call full_z%Add(maxz)
    end if
    num_range = full_z%Count

    !Sort and remove duplicates of exact_z
    if(exact_z%Count>0)then
        call exact_z%sort()
        call exact_z%RemoveDuplicates()
    end if

    !Add exact redshifts to the full array
    izprev = 1
    do i=1, exact_z%count
        if(exact_z%Item(i)==0.d0) cycle
        iz = nint(log(exact_z%Item(i)+1)/dlnz)+1
        if(iz<=izprev) iz=izprev+1
        zcur = exact_z%Item(i)
        call full_z%Add(zcur)
        if(iz<=num_range) then
            call full_z%Swap(iz,full_z%Count)
            call full_z%DeleteItem(full_z%Count)
        end if
        izprev = iz
    end do

    if(full_z%Item(full_z%Count)< maxz) call full_z%Add(maxz)
    !JD added line below to fix interpolation bug when only z=0 is desired
    !since 2D interpolator requires at least 2 data points in each interpolation direction
    !ie k & z.  (Thanks to Tijmen de Haan for spotting the bug).
    if(full_z%Count<2) call full_z%Add(0.1d0)

    this%num_power_redshifts = full_z%Count
    allocate(this%power_redshifts(this%num_power_redshifts))
    this%power_redshifts= full_z%AsArray()
    call full_z%Clear()
    call exact_z%Clear()

    end subroutine Initialize_PKSettings

    subroutine TCosmologyRequirementsLikelihood_InitForSettings(this, Settings)
    class(TCosmologyRequirementsLikelihood) :: this
    class(TCosmoTheorySettings) :: Settings
    integer :: iz, izprev, numz

    if (this%needs_powerspectra .and. this%needs_exact_z) then
        numz = size(Settings%power_redshifts)
        allocate(this%exact_z_index(this%num_z))
        this%exact_z_index = 0
        do iz=1,this%num_z
            izprev=1
            do while(abs(this%exact_z(iz)-Settings%power_redshifts(izprev))>1.d-4)
                izprev=izprev+1
                if(izprev>numz)then
                    call MpiStop("TCosmologyRequirementsLikelihood_InitForSettings: could not find redshift index")
                end if
            end do
            this%exact_z_index(iz) = izprev
        end do
    end if

    end subroutine TCosmologyRequirementsLikelihood_InitForSettings

    end module CosmologyTypes
