    !Define the data types and read/writes them to disk. Also change l_max here.

    module cmbtypes
    use settings
    use likelihood
    use GeneralTypes
    implicit none

    !Number of CMB Cls, 1 for just temperature, 3 (4) for polarization (with B)
    integer, parameter  :: num_cls  = 4

    integer, parameter  :: num_cls_ext=1
    !number of other C_l
    !e.g. 2 for CMB lensing potential and cross-correlation

    !l_max. Tensors are not computed unless compute_tensors = T in input file
    !Make these multiples of 50, should be 50 more than you need accurately
    integer, parameter :: lmax = 6500, lmax_tensor = 400 !note only lmax_computed_cl is actually calculated

    !redshifts for output of BAO_dv background parameters
    real(mcp), target :: z_outputs(1) = [0.57_mcp]

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
    contains
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

    Type, extends(TTheoryPredictions) :: TheoryPredictions
        real(mcp) cl(lmax,num_cls_tot)
        !TT, TE, EE (BB) + other C_l (e.g. lensing)  in that order
        real(mcp) sigma_8
        real(mcp) tensor_ratio_r10, tensor_ratio_02
        integer numderived
        real(mcp) derived_parameters(max_derived_parameters)

        !everything is a function of k/h
        integer   ::  num_k
        real(mcp), dimension(:), pointer :: log_kh => NULL()
        !matpower is log(P_k)
        real(mcp), dimension(:,:), pointer :: matter_power => NULL(), ddmatter_power => NULL()
        real(mcp), dimension(:,:), pointer :: nlmatter_power => NULL(), ddnlmatter_power => NULL()
        real(mcp), dimension(:), pointer :: redshifts => NULL()

    contains
    procedure :: WriteTheory
    procedure :: ReadTheory
    procedure :: WriteBestFitData
    end Type TheoryPredictions

    Type, extends(TParameterization) :: CosmologyParameterization
        logical :: late_time_only = .false.
    end type

    integer, parameter :: As_index=4, amp_ratio_index = 5, Aphiphi_index = 6
    logical :: compute_tensors = .false.

    contains

    subroutine SetTheoryParameterNumbers(slow_num, semi_slow_num)
    use likelihood
    integer, intent (in) :: slow_num, semi_slow_num
    integer i
    class(DataLikelihood), pointer :: DataLike

    num_hard = slow_num
    num_initpower = semi_slow_num
    num_theory_params= num_hard + num_initpower
    index_initpower = num_hard+1
    index_semislow = index_initpower
    index_data =  num_theory_params+1
    if (num_initpower> max_inipower_params) call MpiStop('see cmbtypes.f90: num_initpower> max_inipower_params')
    if (num_theory_params> max_theory_params) call MpiStop('see settings.f90: num_theory_params> max_theory_params')

    do i=1,DataLikelihoods%Count
        DataLike=>DataLikelihoods%Item(i)
        select type (DataLike)
        class is (CosmologyLikelihood)
            if (DataLike%needs_background_functions) DataLike%dependent_params(1:num_hard)=.true.
            if (DataLike%needs_powerspectra) DataLike%dependent_params(1:num_theory_params)=.true.
        end select
    end do

    end subroutine SetTheoryParameterNumbers

    subroutine Initialize_PKSettings()
    use likelihood
    class(DataLikelihood), pointer :: DataLike
    real(mcp) :: dlnz, maxz
    integer :: i,izexact,izrange,iz,izprev
    real(mcp), dimension(:), allocatable :: exact_z,range_redshifts,tmp
    integer :: num_exact = 0
    integer :: num_range = 0
    maxz = 0.
    dlnz = 30.

    do i=1,DataLikelihoods%Count
        DataLike=>DataLikelihoods%Item(i)
        select type (DataLike)
        class is (CosmologyLikelihood)
            if (DataLike%needs_powerspectra) then
                power_kmax = max(power_kmax,DataLike%kmax)
                use_nonlinear = use_nonlinear .or. DataLike%needs_nonlinear_pk
                if(DataLike%needs_exact_z) then
                    if(num_exact==0)then
                        allocate(exact_z(DataLike%num_z))
                        exact_z = DataLike%exact_z
                    else
                        allocate(tmp(num_exact+DataLike%num_z))
                        tmp = [exact_z,DataLike%exact_z]
                        call move_alloc(tmp,exact_z)
                    end if
                    num_exact = num_exact+DataLike%num_z
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
        allocate(range_redshifts(num_range))
        do i=1,num_range
            range_redshifts(i)=dexp(dlnz*i)-1._mcp
        end do
    end if

    !Sort and remove duplicates of exact_z; add exact redshifts to master array
    if(num_exact>0)then
        call quick_sort(exact_z)
        allocate(tmp(num_exact+1))
        tmp(1) = 0._mcp
        iz = 1
        do i=1, num_exact
            if(exact_z(i)/=tmp(iz)) then
                iz = iz+1
                tmp(iz)=exact_z(i)
            end if
        end do
        num_exact = iz
        deallocate(exact_z)
        allocate(exact_z(num_exact))
        exact_z=tmp(1:num_exact)
        deallocate(tmp)
    end if

    !Combine two redshift arrays
    allocate(tmp(num_exact+num_range))
    tmp(1) = 0._mcp
    tmp(2:num_range+1) = range_redshifts(:)
    i=0
    izprev = 1
    do izexact=2, num_exact
        iz = nint(log(exact_z(izexact)+1)/dlnz)+1
        if(iz<=izprev) iz=izprev+1
        if(iz<=num_range+1) then
            tmp(iz)=exact_z(izexact)
        else
            i=i+1
            tmp(num_range+1+i)=exact_z(izexact)
        end if
        izprev=iz
    end do

    if(tmp(num_range+1+i)< maxz) then
        i = i+1
        tmp(num_range+1+i) = maxz
    end if

    num_power_redshifts = num_range+1+i
    allocate(power_redshifts(num_power_redshifts))
    power_redshifts= tmp(1:num_power_redshifts)

    deallocate(tmp)

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

    subroutine WriteTheory(T, i)
    integer i
    Class(TheoryPredictions) T
    integer, parameter :: varcount = 1
    integer tmp(varcount)
    logical, save :: first = .true.

    if (first .and. new_chains) then
        first = .false.
        write(i) use_LSS, compute_tensors
        write(i) lmax, lmax_tensor, num_cls, num_cls_ext
        write(i) varcount
        if (get_sigma8) then
            tmp(1)=1
        else
            tmp(1)=0
        end if
        write(i) tmp(1:varcount)
    end if

    write(i) T%numderived
    write(i) T%derived_parameters(1:T%numderived)
    write(i) T%cl(2:lmax,1:num_cls)
    if (num_cls_ext>0) write(i) T%cl(2:lmax,num_cls+1:num_cls_tot)

    if (compute_tensors) then
        write(i) T%tensor_ratio_02, T%tensor_ratio_r10
    end if

    if (get_sigma8 .or. use_LSS) write(i) T%sigma_8

    if (use_LSS) then
        write(i) T%num_k, size(T%redshifts)
        write(i) T%log_kh
        write(i) T%redshifts
        write(i) T%matter_power
    end if

    end subroutine WriteTheory

    subroutine ReadTheory(T, i)
    Class(TheoryPredictions) T
    integer, intent(in) :: i
    integer unused
    logical, save :: first = .true.
    logical, save :: has_sigma8, has_LSS, has_tensors
    integer, save :: almax, almaxtensor, anumcls, anumclsext, tmp(1)
    logical, save :: planck1_format
    real(mcp) old_matterpower(74,1)
    !JD 10/13 new variables for handling new pk arrays
    integer :: num_z

    if (first) then
        first = .false.
        read(i) has_LSS, has_tensors
        read(i) almax, almaxtensor, anumcls, anumclsext
        if (almax > lmax) call MpiStop('ReadTheory: reading file with larger lmax')
        if (anumcls /= num_cls) call MpiStop('ReadTheory: reading file with different Cls')
        if (anumclsext /= num_cls_ext) call MpiStop('ReadTheory: reading file with different ext Cls')
        read(i) unused
        planck1_format = unused==0
        if (unused>0) read(i) tmp(1:unused)
        if (.not. planck1_format) has_sigma8 = tmp(1)==1
    end if

    T%cl = 0
    T%derived_parameters=0
    read(i) T%numderived
    read(i) T%derived_parameters(1:T%numderived)
    read(i) T%cl(2:almax,1:anumcls)
    if (anumclsext >0) read(i) T%cl(2:almax,num_cls+1:num_cls+anumclsext)

    if (has_tensors) then
        read(i) T%tensor_ratio_02, T%tensor_ratio_r10
    end if

    if (planck1_format) then
        if (has_LSS) then
            read(i) T%sigma_8, old_matterpower
            !discard unused mpk array
        end if
    else
        if (has_sigma8 .or. has_LSS) read(i) T%sigma_8
        if (has_LSS) then
            call InitPK(T)
            read(i) T%num_k, num_z
            allocate(T%matter_power(T%num_k,num_z))
            allocate(T%ddmatter_power(T%num_k,num_z))
            allocate(T%log_kh(T%num_k))
            allocate(T%redshifts(num_z))
            read(i) T%log_kh
            read(i) T%redshifts
            read(i) T%matter_power
            call IOTheory_GetNLAndSplines(T)
        end if
    end if

    end subroutine ReadTheory

    subroutine ClsFromTheoryData(T, Cls)
    Type(TheoryPredictions) T
    real(mcp) Cls(lmax,num_cls_tot)

    Cls(2:lmax,1:num_clsS) =T%cl(2:lmax,1:num_clsS)
    if (num_cls>3 .and. num_ClsS==3) Cls(2:lmax,num_cls)=0
    if (num_cls_ext>0) then
        Cls(2:lmax,num_cls+1:num_cls_tot) =T%cl(2:lmax,num_clsS+1:num_clsS+num_cls_ext)
    end if

    end subroutine ClsFromTheoryData

    subroutine WriteTextCls(aname,T)
    Type(TheoryPredictions) T
    character (LEN=*), intent(in) :: aname
    integer l
    real(mcp) Cls(lmax,num_cls_tot), nm
    character(LEN=80) fmt

    call CreateTxtFile(aname,tmp_file_unit)
    call ClsFromTheoryData(T, Cls)
    fmt = concat('(1I6,',num_cls_tot,'E15.5)')
    do l = 2, lmax
        nm = 2*pi/(l*(l+1))
        if (num_cls_ext > 0) then
            write (tmp_file_unit,fmt) l, cls(l,1:num_cls)/nm, cls(l,num_cls+1:num_cls_tot)
        else
            write (tmp_file_unit,fmt) l, cls(l,:)/nm
        end if
    end do
    call CloseFile(tmp_file_unit)

    end subroutine WriteTextCls

    subroutine WriteBestFitData(Theory,fnameroot)
    class(TheoryPredictions) Theory
    character(LEN=*), intent(in) :: fnameroot

    call WriteTextCls(fnameroot //'.bestfit_cl', Theory)

    end subroutine WriteBestFitData

    subroutine InitPK(Theory)
    Type(TheoryPredictions) :: Theory
    integer i

    deallocate(Theory%log_kh,stat=i)
    deallocate(Theory%matter_power,stat=i)
    deallocate(Theory%ddmatter_power,stat=i)
    deallocate(Theory%nlmatter_power,stat=i)
    deallocate(Theory%ddnlmatter_power,stat=i)
    deallocate(Theory%redshifts,stat=i)
    call PK_Nullify(Theory)

    end subroutine InitPK

    subroutine PK_Nullify(Theory)
    Type(TheoryPredictions) :: Theory

    nullify(Theory%log_kh)
    nullify(Theory%ddmatter_power, Theory%nlmatter_power)
    nullify(Theory%nlmatter_power, Theory%ddnlmatter_power)
    nullify(Theory%redshifts)

    end subroutine PK_Nullify

    subroutine IOTheory_GetNLAndSplines(Theory)
    use Transfer
    Type(TheoryPredictions) Theory
    Type(MatterPowerData) :: Cosmo_PK
    integer :: nz, zix

    Cosmo_PK%num_k = Theory%num_k
    Cosmo_PK%num_z = 1
    allocate(Cosmo_PK%matpower(Cosmo_PK%num_k,1))
    allocate(Cosmo_PK%ddmat(Cosmo_PK%num_k,1))
    allocate(Cosmo_PK%nonlin_ratio(Cosmo_PK%num_k,1))
    allocate(Cosmo_PK%log_kh(Cosmo_PK%num_k))
    allocate(Cosmo_PK%redshifts(1))
    Cosmo_PK%log_kh = Theory%log_kh

    do zix=1, size(Theory%redshifts)
        Cosmo_PK%redshifts(1) = Theory%redshifts(zix)
        Cosmo_PK%matpower(:,1) = Theory%matter_power(:,zix)
        call MatterPowerdata_getsplines(Cosmo_PK)
        Theory%ddmatter_power(:,zix) = Cosmo_PK%ddmat(:,1)
        if(use_nonlinear) then
            call MatterPowerdata_MakeNonlinear(Cosmo_PK)
            Theory%nlmatter_power(:,zix) = Cosmo_PK%matpower(:,1)
            Theory%ddnlmatter_power(:,zix) = Cosmo_PK%ddmat(:,1)
        end if
    end do

    call MatterPowerdata_Free(Cosmo_PK)

    end subroutine IOTheory_GetNLAndSplines

    recursive subroutine quick_sort(list)
    real(mcp), dimension(:), intent(in out) :: list
    integer :: i, j, n
    real(mcp) :: chosen, temp
    integer, parameter :: max_simple_sort_size = 6

    n = size(list)
    if (n <= max_simple_sort_size) then
        ! Use interchange sort for small lists
        call interchange_sort(list)
    else
        ! Use partition (“quick”) sort chosen = list(n/2)
        i=0
        j=n+1
        do
            ! Scan list from left end
            ! until element >= chosen is found
            do
                i=i+1
                if (list(i) >= chosen) exit
            end do
            ! Scan list from right end
            ! until element <= chosen is found
            do
                j=j-1
                if (list(j) <= chosen) exit
            end do

            if (i < j) then
                ! Swap two out of place elements
                temp = list(i)
                list(i) = list(j)
                list(j) = temp
            else if (i == j) then
                i=i+1
                exit
            else
                exit
            end if
        end do
        if (1 < j) call quick_sort(list(:j))
        if (i < n) call quick_sort(list(i:))
    end if  ! test for small array

    end subroutine quick_sort

    subroutine interchange_sort(list)
    real(mcp), dimension(:), intent(in out) :: list
    integer :: i, j
    real(mcp) :: temp
    do i = 1, size(list) - 1
        do j = i + 1, size(list)
            if (list(i) >  list(j)) then
                temp = list(i)
                list(i) = list(j)
                list(j) = temp
            end if
        end do
    end do
    end subroutine interchange_sort

    end module cmbtypes
