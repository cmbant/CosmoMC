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
    integer :: num_matter_power

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


    type, extends(DatasetFileLikelihood) :: CosmologyLikelihood
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
    end type CosmologyLikelihood

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

    Type, extends(TTheoryPredictions) :: TCosmoTheoryPredictions
        real(mcp) cl(lmax,num_cls_tot)
        !TT, TE, EE (BB) + other C_l (e.g. lensing)  in that order
        real(mcp) sigma_8
        real(mcp) tensor_ratio_r10, tensor_ratio_02
        integer numderived
        real(mcp) derived_parameters(max_derived_parameters)

        !everything is a function of k/h
        integer   ::  num_k
        real(mcp), dimension(:), allocatable :: log_kh
        !matpower is log(P_k)
        real(mcp), dimension(:,:), allocatable :: matter_power, ddmatter_power
        real(mcp), dimension(:,:), allocatable :: nlmatter_power, ddnlmatter_power
        real(mcp), dimension(:), allocatable :: redshifts
    contains
    procedure :: WriteTheory => TCosmoTheoryPredictions_WriteTheory
    procedure :: ReadTheory => TCosmoTheoryPredictions_ReadTheory
    procedure :: WriteBestFitData => TCosmoTheoryPredictions_WriteBestFitData
    end Type TCosmoTheoryPredictions

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
    class(DataLikelihood), pointer :: DataLike

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
        class is (CosmologyLikelihood)
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
    class(DataLikelihood), pointer :: DataLike
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
        class is (CosmologyLikelihood)
            if (DataLike%needs_powerspectra) then
                power_kmax = max(power_kmax,DataLike%kmax)
                use_nonlinear = use_nonlinear .or. DataLike%needs_nonlinear_pk
                if(DataLike%needs_exact_z) then
                    call exact_z%AddArrayItems(DataLike%exact_z)
                else
                    num_range = DataLike%num_z
                    maxz = max(maxz,DataLike%max_z)
                    if(num_range >1 .and. maxz>0) then
                        dlnz = min(dlnz,log(DataLike%max_z+1)/(num_range-1))
                    else if(num_range==1 .and. maxz > 0)then
                        write(*,'("ERROR: ",A," dataset: ",A, "wants only 1 redshift")')&
                        trim(DataLike%LikelihoodType),trim(DataLike%name)
                        write(*,'("       but wants a maximum redshift of ",F7.2,".  Check dataset settings!")')maxz
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
    class(DataLikelihood), pointer :: DataLike
    real(mcp), dimension(:) :: A
    integer, intent(out), optional :: error
    integer :: i, iz, izprev, numz

    numz = size(A)

    do i=1,DataLikelihoods%Count
        DataLike=>DataLikelihoods%Item(i)
        select type (DataLike)
        class is (CosmologyLikelihood)
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

    subroutine TCosmoTheoryPredictions_WriteTheory(T, unit)
    integer unit
    Class(TCosmoTheoryPredictions) T
    integer, parameter :: varcount = 1
    integer tmp(varcount)
    logical, save :: first = .true.

    if (first .and. new_chains) then
        first = .false.
        write(unit) use_LSS, compute_tensors
        write(unit) lmax, lmax_tensor, num_cls, num_cls_ext
        write(unit) varcount
        if (get_sigma8) then
            tmp(1)=1
        else
            tmp(1)=0
        end if
        write(unit) tmp(1:varcount)
    end if

    write(unit) T%numderived
    write(unit) T%derived_parameters(1:T%numderived)
    write(unit) T%cl(2:lmax,1:num_cls)
    if (num_cls_ext>0) write(unit) T%cl(2:lmax,num_cls+1:num_cls_tot)

    if (compute_tensors) then
        write(unit) T%tensor_ratio_02, T%tensor_ratio_r10
    end if

    if (get_sigma8 .or. use_LSS) write(unit) T%sigma_8

    if (use_LSS) then
        write(unit) T%num_k, size(T%redshifts)
        write(unit) T%log_kh
        write(unit) T%redshifts
        write(unit) T%matter_power
        write(unit) use_nonlinear
        if(use_nonlinear) write(unit) T%nlmatter_power
    end if

    end subroutine TCosmoTheoryPredictions_WriteTheory

    subroutine TCosmoTheoryPredictions_ReadTheory(T, unit)
    Class(TCosmoTheoryPredictions) T
    integer, intent(in) :: unit
    integer unused
    logical, save :: first = .true.
    logical, save :: has_sigma8, has_LSS, has_tensors
    integer, save :: almax, almaxtensor, anumcls, anumclsext, tmp(1)
    logical, save :: planck1_format
    real(mcp) old_matterpower(74,1)
    !JD 10/13 new variables for handling new pk arrays
    integer :: num_z, stat
    logical  has_nonlinear

    if (first) then
        first = .false.
        read(unit) has_LSS, has_tensors
        read(unit) almax, almaxtensor, anumcls, anumclsext
        if (almax > lmax) call MpiStop('ReadTheory: reading file with larger lmax')
        if (anumcls /= num_cls) call MpiStop('ReadTheory: reading file with different Cls')
        if (anumclsext /= num_cls_ext) call MpiStop('ReadTheory: reading file with different ext Cls')
        read(unit) unused
        planck1_format = unused==0
        if (unused>0) read(unit) tmp(1:unused)
        if (.not. planck1_format) has_sigma8 = tmp(1)==1
    end if

    T%cl = 0
    T%derived_parameters=0
    read(unit) T%numderived
    read(unit) T%derived_parameters(1:T%numderived)
    read(unit) T%cl(2:almax,1:anumcls)
    if (anumclsext >0) read(unit) T%cl(2:almax,num_cls+1:num_cls+anumclsext)

    if (has_tensors) then
        read(unit) T%tensor_ratio_02, T%tensor_ratio_r10
    end if

    if (planck1_format) then
        if (has_LSS) then
            read(unit) T%sigma_8, old_matterpower
            !discard unused mpk array
        end if
    else
        if (has_sigma8 .or. has_LSS) read(unit) T%sigma_8
        if (has_LSS) then
            read(unit) T%num_k, num_z
            call InitPK(T,T%num_k,num_z)
            read(unit) T%log_kh
            read(unit) T%redshifts
            read(unit, iostat=stat) T%matter_power
            if (IS_IOSTAT_END(stat) .and. use_nonlinear) then
                deallocate(T%nlmatter_power,T%ddnlmatter_power)
                write(*,*)"ReadTheory:  You want a nonlinear MPK but but your datafile is "
                write(*,*)"in an old format does not include one."
                write(*,*)"Make sure you set redo_pk = T or the program will fail"
            else
                read(unit)has_nonlinear
                if(has_nonlinear .and. use_nonlinear) then
                    read(unit)T%nlmatter_power
                else if(has_nonlinear .and. .not. use_nonlinear) then
                    if(allocated(T%nlmatter_power))deallocate(T%nlmatter_power)
                    if(allocated(T%ddnlmatter_power))deallocate(T%ddnlmatter_power)
                    allocate(T%nlmatter_power(T%num_k,num_z))
                    allocate(T%ddnlmatter_power(T%num_k,num_z))
                    write(*,*)"Your data files have nonlinear power spectra, but you are not using"
                    write(*,*)"nonlinear power spectra.  Be careful that this is what you intended"
                    read(unit)T%nlmatter_power
                else 
                    deallocate(T%nlmatter_power,T%ddnlmatter_power)
                    write(*,*)"ReadTheory:  You want a nonlinear MPK but but your datafile does not include one."
                    write(*,*)"Make sure you set redo_pk = T or the program will fail"
                end if
            end if
            call IOTheory_GetSplines(T)
        end if
    end if

    end subroutine TCosmoTheoryPredictions_ReadTheory


    subroutine ClsFromTheoryData(T, Cls)
    Type(TCosmoTheoryPredictions) T
    real(mcp) Cls(lmax,num_cls_tot)

    Cls(2:lmax,1:num_clsS) =T%cl(2:lmax,1:num_clsS)
    if (num_cls>3 .and. num_ClsS==3) Cls(2:lmax,num_cls)=0
    if (num_cls_ext>0) then
        Cls(2:lmax,num_cls+1:num_cls_tot) =T%cl(2:lmax,num_clsS+1:num_clsS+num_cls_ext)
    end if

    end subroutine ClsFromTheoryData

    subroutine WriteTextCls(aname,T)
    Type(TCosmoTheoryPredictions) T
    character (LEN=*), intent(in) :: aname
    integer l
    real(mcp) Cls(lmax,num_cls_tot), nm
    character(LEN=80) fmt
    integer unit

    unit=CreateNewTxtFile(aname)
    call ClsFromTheoryData(T, Cls)
    fmt = concat('(1I6,',num_cls_tot,'E15.5)')
    do l = 2, lmax
        nm = 2*pi/(l*(l+1))
        if (num_cls_ext > 0) then
            write (unit,fmt) l, cls(l,1:num_cls)/nm, cls(l,num_cls+1:num_cls_tot)
        else
            write (unit,fmt) l, cls(l,:)/nm
        end if
    end do
    close(unit)

    end subroutine WriteTextCls

    subroutine TCosmoTheoryPredictions_WriteBestFitData(Theory,fnameroot)
    class(TCosmoTheoryPredictions) Theory
    character(LEN=*), intent(in) :: fnameroot

    if (use_CMB) call WriteTextCls(fnameroot //'.bestfit_cl', Theory)

    end subroutine TCosmoTheoryPredictions_WriteBestFitData

    subroutine InitPK(Theory, num_k, num_z)
    Type(TCosmoTheoryPredictions) :: Theory
    integer, intent(in) :: num_k, num_z

    if(allocated(Theory%log_kh))deallocate(Theory%log_kh)
    if(allocated(Theory%matter_power))deallocate(Theory%matter_power)
    if(allocated(Theory%ddmatter_power))deallocate(Theory%ddmatter_power)
    if(allocated(Theory%redshifts))deallocate(Theory%redshifts)
    allocate(Theory%log_kh(num_k))
    allocate(Theory%matter_power(num_k,num_z))
    allocate(Theory%ddmatter_power(num_k,num_z))
    allocate(Theory%redshifts(num_z))
    if(use_nonlinear) then
        if(allocated(Theory%nlmatter_power))deallocate(Theory%nlmatter_power)
        if(allocated(Theory%ddnlmatter_power))deallocate(Theory%ddnlmatter_power)
        allocate(Theory%nlmatter_power(num_k,num_z))
        allocate(Theory%ddnlmatter_power(num_k,num_z))
    end if

    end subroutine InitPK

    subroutine IOTheory_GetSplines(Theory)
    use Interpolation, only : spline, SPLINE_DANGLE
    Type(TCosmoTheoryPredictions) Theory
    integer :: zix,num_k,nz

    num_k = Theory%num_k
    nz = size(Theory%redshifts)

    do zix=1, nz
        call spline(Theory%log_kh,Theory%matter_power(:,zix),num_k,SPLINE_DANGLE,&
        SPLINE_DANGLE,Theory%ddmatter_power(:,zix))

        if(use_nonlinear .and. allocated(Theory%nlmatter_power))&
        call spline(Theory%log_kh,Theory%nlmatter_power(:,zix),num_k,SPLINE_DANGLE,&
        SPLINE_DANGLE,Theory%ddnlmatter_power(:,zix))
    end do

    end subroutine IOTheory_GetSplines

    end module cmbtypes