    !Define the data types and read/writes them to disk. Also change l_max here.

    module cmbtypes
    use settings
    use likelihood
    use GeneralTypes
    use ObjectLists
    implicit none

    !Number of CMB Cls, 1 for just temperature, 3 (4) for polarization (with B)
    integer, parameter  :: num_cls  = 4

    integer, parameter  :: num_cls_ext=1
    !number of other C_l
    !e.g. 2 for CMB lensing potential and cross-correlation


    logical :: get_sigma8 = .true.
    logical :: Use_LSS = .false.
    logical :: Use_CMB = .false.
    logical :: use_nonlinear = .false.    !JD for WiggleZ MPK


    !l_max. Tensors are not computed unless compute_tensors = T in input file
    !Make these multiples of 50, should be 50 more than you need accurately
    integer, parameter :: lmax = 6500, lmax_tensor = 400 !note only lmax_computed_cl is actually calculated

    !redshifts for output of BAO_dv background parameters
    real(mcp), target :: z_outputs(1) = [0.57_mcp]

    integer, parameter  :: derived_age=1, derived_zstar=2, derived_rstar=3, derived_thetastar=4,derived_zdrag=5, &
    derived_rdrag=6,derived_kD=7,derived_thetaD=8 , derived_zEQ =9, derived_thetaEQ=10 !index in derived parameters array

    logical :: CMB_lensing = .true.
    logical :: use_lensing_potential = .false.
    logical :: use_nonlinear_lensing = .false.
    integer :: lmax_computed_cl = lmax !value used for CAMB
    logical :: compute_tensors = .false.

    integer, parameter :: As_index=4, amp_ratio_index = 5, Aphiphi_index = 6

    !Parameters for calculating/storing the matter power spectrum
    real(mcp) :: power_kmax = 0.8
    integer :: num_power_redshifts
    real(mcp), dimension(:), allocatable :: power_redshifts

    !Only used in params_CMB
    real(mcp) :: pivot_k = 0.05_mcp !Point for defining primordial power spectra
    logical :: inflation_consistency = .false. !fix n_T or not

    logical :: bbn_consistency = .true. !JH

    integer :: num_massive_neutrinos = 3 !if >0, number of massive degenerate eigenstates

    real(mcp), parameter :: cl_norm = 1e-10_mcp !units for As

    integer, parameter :: max_derived_parameters = 30

    integer, parameter :: num_cls_tot = num_cls + num_cls_ext
    !Number of scalar-only cls
    !if num_cls=4 and CMB_lensing then increased to 4
    integer :: num_clsS=min(num_cls,3)

    integer, parameter :: max_inipower_params = 10
    integer:: num_hard, num_initpower
    integer :: index_initpower


    type, extends(TDatasetFileLikelihood) :: TCosmologyLikelihood
        !Don't have to use extract features of DatasetFileLikelihood
        !not implemented yet..
        !        integer :: needs_cl_lmax = 0
        logical :: needs_background_functions = .true.
        logical :: needs_powerspectra = .false.

        logical :: needs_nonlinear_pk = .false.
        logical :: needs_exact_z = .false.
        integer :: num_z = 0
        real(mcp), dimension(:), allocatable :: exact_z
        integer, dimension(:), allocatable :: exact_z_index
        real(mcp) :: max_z
        real(mcp) :: kmax = 0.8
    end type TCosmologyLikelihood

    Type, extends(TTheoryParams) :: CMBParams
        real(mcp) InitPower(max_inipower_params)
        !These are fast paramters for the initial power spectrum
        !Now remaining (non-independent) parameters
        real(mcp) omb, omc, omv, omnu, omk, omdm
        real(mcp) ombh2, omch2, omnuh2, omdmh2
        real(mcp) zre, zre_delta, nufrac
        real(mcp) h, H0, tau
        real(mcp) w, wa
        real(mcp) YHe, nnu, iso_cdm_correlated, ALens, fdm !fdm is dark matter annihilation, eg,. 0910.3663
        real(mcp) :: omnuh2_sterile = 0._mcp  !note omnhu2 is the sum of this + standard neutrinos
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
    use likelihood
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
    if (num_initpower> max_inipower_params) call MpiStop('see cmbtypes.f90: num_initpower> max_inipower_params')

    do i=1,DataLikelihoods%Count
        DataLike=>DataLikelihoods%Item(i)
        select type (DataLike)
        class is (TCosmologyLikelihood)
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


    subroutine Initialize_PKSettings()
    use likelihood
    class(TDataLikelihood), pointer :: DataLike
    Type(TRealList) :: exact_z, full_z
    real(mcp) :: dlnz, maxz, zcur
    integer :: i,iz,izprev
    integer :: num_range
    maxz = 0.
    dlnz = 30.

    call full_z%Add(0.d0)

    do i=1,DataLikelihoods%Count
        DataLike=>DataLikelihoods%Item(i)
        select type (DataLike)
        class is (TCosmologyLikelihood)
            if (DataLike%needs_powerspectra) then
                power_kmax = max(power_kmax,DataLike%kmax)
                use_nonlinear = use_nonlinear .or. DataLike%needs_nonlinear_pk
                if(DataLike%needs_exact_z) then
                    call exact_z%AddArrayItems(DataLike%exact_z)
                else
                    num_range = DataLike%num_z
                    maxz = max(maxz,DataLike%max_z)
                    if(num_range >2 .and. maxz>0) then
                        dlnz = min(dlnz,log(DataLike%max_z+1)/(num_range-1))
                    else if(num_range<3 .and. maxz > 0)then
                        write(*,'("ERROR: ",A," dataset: ",A, "wants less than 3 redshifts")')&
                        trim(DataLike%LikelihoodType),trim(DataLike%name)
                        write(*,'("       but wants a maximum redshift of ",F7.2,". A minimum ")')maxz
                        write(*,*)"       of 3 redshifts is required or PowerAtZ will fail."
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

    !Build array of redshifts where the redshift exact value doesn't matter
    if(maxz>0)then
        num_range = ceiling(log(maxz+1)/dlnz)
        dlnz = log(maxz+1)/(num_range)
        do i=1,num_range
            zcur = dexp(dlnz*i)-1
            call full_z%Add(zcur)
        end do
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

    num_power_redshifts = full_z%Count
    allocate(power_redshifts(num_power_redshifts))
    power_redshifts= full_z%AsArray()
    call full_z%Clear()
    call exact_z%Clear()
    call IndexExactRedshifts(power_redshifts)

    end subroutine Initialize_PKSettings

    subroutine IndexExactRedshifts(A,error)
    use likelihood
    class(TDataLikelihood), pointer :: DataLike
    real(mcp), dimension(:) :: A
    integer, intent(out), optional :: error
    integer :: i, iz, izprev, numz

    numz = size(A)

    do i=1,DataLikelihoods%Count
        DataLike=>DataLikelihoods%Item(i)
        select type (DataLike)
        class is (TCosmologyLikelihood)
            if (DataLike%needs_powerspectra) then
                if(DataLike%needs_exact_z) then
                    DataLike%exact_z_index = 0
                    do iz=1,DataLike%num_z
                        izprev=1
                        do while(abs(DataLike%exact_z(iz)-A(izprev))>1.d-4)
                            izprev=izprev+1
                            if(izprev>numz)then
                                write(*,*) "ERROR, In IndexExactRedshifts: could not find redshift index"
                                if(present(error))error=1
                                return
                            end if
                        end do
                        DataLike%exact_z_index(iz) = izprev
                    end do
                end if
            end if
        end select
    end do

    end subroutine IndexExactRedshifts

    end module cmbtypes