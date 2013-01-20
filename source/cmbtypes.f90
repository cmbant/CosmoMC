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

    !Parameters for calculating/storing the matter power spectrum
    !Note that by default everything is linear

    !Note these are the interpolated/extrapolated values. The k at which matter power is computed up to 
    !by CAMB is set in CMB_Cls_xxx with, e.g. P%Transfer%kmax = 0.6 (which is enough for 2dF)

    !Old mpk settings
#ifdef DR71RG
    !!! BR09: Reid et al 2009 settings for the LRG power spectrum.
    integer, parameter :: num_matter_power = 300 !number of points computed in matter power spectrum
    real(mcp), parameter    :: matter_power_minkh =  0.999e-4_mcp  !minimum value of k/h to store
    real(mcp), parameter    :: matter_power_dlnkh = 0.03_mcp     !log spacing in k/h
    real(mcp), parameter    :: matter_power_maxz = 1._mcp !Not used, but must be non-zero to avoid error when have 4 z steps and use_mpk=F
    integer, parameter :: matter_power_lnzsteps = 4  ! z=0 to get sigma8 (this first entry appears to be coded in some spots in the code!!), plus 3 LRG redshifts.
#else
    integer, parameter :: num_matter_power = 74 !number of points computed in matter power spectrum
    real(mcp), parameter    :: matter_power_minkh =  0.999e-4_mcp  !1e-4 !minimum value of k/h to store
    real(mcp), parameter    :: matter_power_dlnkh = 0.143911568_mcp     !log spacing in k/h
    real(mcp), parameter    :: matter_power_maxz = 0._mcp    !6.0
    integer, parameter :: matter_power_lnzsteps = 1 !20
#endif
    !Only used in params_CMB
    real(mcp) :: pivot_k = 0.05_mcp !Point for defining primordial power spectra
    logical :: inflation_consistency = .false. !fix n_T or not

    logical :: bbn_consistency = .true. !JH

    integer :: num_massive_neutrinos = -1 !if >0, number of massive degenerate eigenstates
    logical :: nonthermal_masive_neutrinos = .false.
    logical :: neutrino_param_mnu = .true. !parameter 6 is sum mnu (false for old behaviour of param(6) is fnu)

    real(mcp), parameter :: cl_norm = 1e-10_mcp !units for As

    integer, parameter :: max_derived_parameters = 20

    integer, parameter :: num_cls_tot = num_cls + num_cls_ext
    !Number of scalar-only cls
    !if num_cls=4 and CMB_lensing then increased to 4 
    integer :: num_clsS=min(num_cls,3) 

    integer, parameter :: max_inipower_params = 10
    integer:: num_hard, num_initpower
    integer :: index_initpower


    type, extends(DataLikelihood) :: CosmologyLikelihood
        !not implemented yet..
        !        logical :: needs_linear_pk = .false.
        !        integer :: needs_cl_lmax = 0
        logical :: needs_background_functions = .true.
        logical :: needs_powerspectra = .false.
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
        real(mcp) reserved(5)
    end Type CMBParams

    Type, extends(TTheoryPredictions) :: TheoryPredictions
        real(mcp) cl(lmax,num_cls_tot)
        !TT, TE, EE (BB) + other C_l (e.g. lensing)  in that order
        real(mcp) sigma_8
        real(mcp) tensor_ratio_r10, tensor_ratio_02
        integer numderived
        real(mcp) derived_parameters(max_derived_parameters)

        real(mcp) matter_power(num_matter_power,matter_power_lnzsteps)
        !second index is redshifts from 0 to matter_power_maxz
        !if custom_redshift_steps = false with equal spacing in
        !log(1+z) and matter_power_lnzsteps points
        !if custom_redshift_steps = true set in mpk.f90 
        ! BR09 additions
        real(mcp) mpk_nw(num_matter_power,matter_power_lnzsteps) !no wiggles fit to matter power spectrum
        real(mcp) mpkrat_nw_nl(num_matter_power,matter_power_lnzsteps) !halofit run on mpk_nw
        real(mcp) finalLRGtheoryPk(num_matter_power)  !! this is the quantity that enters the LRG likelihood calculation
        ! end BR09 additions
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

    subroutine WriteTheory(T, i)
    integer i
    Class(TheoryPredictions) T
    integer unused
    logical, save :: first = .true.

    if (first .and. new_chains) then
        first = .false.
        write(i) use_LSS, compute_tensors
        write(i) lmax, lmax_tensor, num_cls, num_cls_ext
        unused=0
        write(i) unused
    end if

    write(i) T%numderived
    write(i) T%derived_parameters(1:T%numderived) 
    write(i) T%cl(2:lmax,1:num_cls)
    if (num_cls_ext>0) write(i) T%cl(2:lmax,num_cls+1:num_cls_tot)

    if (compute_tensors) then
        write(i) T%tensor_ratio_02, T%tensor_ratio_r10
    end if
    if (use_LSS) then
        write(i) T%sigma_8, T%matter_power
    end if

    end subroutine WriteTheory

    subroutine ReadTheory(T, i)
    Class(TheoryPredictions) T
    integer, intent(in) :: i
    integer unused
    logical, save :: first = .true.
    logical, save :: has_LSS, has_tensors
    integer, save :: almax, almaxtensor, anumcls, anumclsext, tmp(1)

    if (first) then
        first = .false.
        read(i) has_LSS, has_tensors
        read(i) almax, almaxtensor, anumcls, anumclsext
        if (almax > lmax) call MpiStop('ReadTheory: reading file with larger lmax')
        if (anumcls /= num_cls) call MpiStop('ReadTheory: reading file with different Cls')
        if (anumclsext /= num_cls_ext) call MpiStop('ReadTheory: reading file with different ext Cls')
        read(i) unused
        if (unused>0) read(i) tmp(1:unused)
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

    if (has_LSS) then
        read(i) T%sigma_8, T%matter_power
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

    function MatterPowerAt(T,kh)
    !get matter power spectrum today at kh = k/h by interpolation from stored values
    real(mcp), intent(in) :: kh
    Type(TheoryPredictions) T
    real(mcp) MatterPowerAt
    real(mcp) x, d
    integer i

    x = log(kh/matter_power_minkh) / matter_power_dlnkh
    if (x < 0 .or. x >= num_matter_power-1) then
        write (*,*) ' k/h out of bounds in MatterPowerAt (',kh,')'
        call MpiStop('') 
    end if
    i = int(x)
    d = x - i
    MatterPowerAt = exp(log(T%matter_power(i+1,1))*(1-d) &
    + log(T%matter_power(i+2,1))*d)
    !Just do linear interpolation in logs for now..
    !(since we already cublic-spline interpolated to get the stored values)
    !Assume matter_power_lnzsteps is at redshift zero
    end function



    !BR09 this function is just a copy of the one above but with LRG theory put in instead of linear theory
    function LRGPowerAt(T,kh)
    !get LRG matter power spectrum today at kh = k/h by interpolation from stored values
    real(mcp), intent(in) :: kh
    Type(TheoryPredictions) T
    real(mcp) LRGPowerAt
    real(mcp) x, d
    integer i

    x = log(kh/matter_power_minkh) / matter_power_dlnkh
    if (x < 0 .or. x >= num_matter_power-1) then
        write (*,*) ' k/h out of bounds in MatterPowerAt (',kh,')'
        call MpiStop('') 
    end if
    i = int(x)
    d = x - i
    LRGPowerAt = exp(log(T%finalLRGtheoryPk(i+1))*(1-d) + log(T%finalLRGtheoryPk(i+2))*d)
    !Just do linear interpolation in logs for now..
    !(since we already cublic-spline interpolated to get the stored values)
    end function
    !!BRO09 addition end

    function MatterPowerAt_Z(T,kh,z)
    !get matter power spectrum at z at kh = k/h by interpolation from stored values

    real(mcp), intent(in) :: kh
    Type(TheoryPredictions) T
    real(mcp) MatterPowerAt_Z
    real(mcp) x, d, z, y, dz, mup, mdn
    real(mcp) matter_power_dlnz
    integer i, iz

    matter_power_dlnz = log(matter_power_maxz+1) / (matter_power_lnzsteps -1 + 1e-13)
    y = log(1.+ z) / matter_power_dlnz 

    if (z > matter_power_maxz ) then
        write (*,*) ' z out of bounds in MatterPowerAt_Z (',z,')'
        call MpiStop('')
    end if
    x = log(kh/matter_power_minkh) / matter_power_dlnkh
    if (x < 0 .or. x >= num_matter_power-1) then
        write (*,*) ' k/h out of bounds in MatterPowerAt_Z (',kh,')'
        call MpiStop('')
    end if

    iz = int(y*0.99999999)
    dz = y - iz

    i = int(x)
    d = x - i

    mup = log(T%matter_power(i+1,iz+2))*(1-d) + log(T%matter_power(i+2,iz+2))*d
    mdn = log(T%matter_power(i+1,iz+1))*(1-d) + log(T%matter_power(i+2,iz+1))*d

    MatterPowerAt_Z = exp(mdn*(1-dz) + mup*dz)

    end function MatterPowerAt_Z

    end module cmbtypes
