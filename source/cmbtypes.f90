!Define the data types and read/writes them to disk. Also change l_max here.

module cmbtypes
use settings
implicit none


!Number of CMB Cls, 1 for just temperature, 3 (4) for polarization (with B)
  integer, parameter  :: num_cls  = 4

  integer, parameter  :: num_cls_ext=0
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
   logical :: neutrino_param_mnu = .true. !parameter 6 is sum mnu (false for old behaviour of param(6) is fnu)

   integer, parameter :: max_derived_parameters = 20
  
  integer, parameter :: num_cls_tot = num_cls + num_cls_ext
!Number of scalar-only cls
!if num_cls=4 and CMB_lensing then increased to 4 
  integer :: num_clsS=min(num_cls,3) 

  Type CMBParams
     real(mcp) nuisance(1:num_nuisance_params)
      !unit Gaussians for experimental parameters
     real(mcp) data_params(1:num_freq_params)
      !These are fast parameters controling amplitudes, calibrations, etc.
     real(mcp) InitPower(1:num_initpower)
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

  Type TheoryPredictions
     real(mcp) r10
     real(mcp) cl(lmax,num_cls_tot), cl_tensor(lmax_tensor,num_cls) 
      !TT, TE, EE (BB) + other C_l (e.g. lensing)  in that order
     real(mcp) sigma_8, tensor_ratio_02
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
  end Type TheoryPredictions

  integer, parameter :: As_index=4, amp_ratio_index = 5
  logical :: compute_tensors = .false.

contains

   subroutine WriteTheory(T, i)
    integer i
    Class(TheoryPredictions) T
    integer unused
    logical, save :: first = .true.
    
    if (first) then
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
          write(i) T%tensor_ratio_02, T%r10
          write(i) T%cl_tensor(2:lmax_tensor,1:num_cls)
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
        T%cl_tensor = 0
        T%derived_parameters=0
        read(i) T%numderived
        read(i) T%derived_parameters(1:T%numderived)
        read(i) T%cl(2:almax,1:anumcls)
        if (anumclsext >0) read(i) T%cl(2:almax,num_cls+1:num_cls+anumclsext)

        if (has_tensors) then
          read(i) T%tensor_ratio_02, T%r10
          read(i) T%cl_tensor(2:almaxtensor,1:anumcls)
        end if

        if (has_LSS) then
          read(i) T%sigma_8, T%matter_power
        end if

   end subroutine ReadTheory

   subroutine ClsFromTheoryData(T, CMB, Cls)
     Type(TheoryPredictions) T
     Type(CMBParams) CMB
     real(mcp) Cls(lmax,num_cls_tot)
     integer i
     
     Cls(2:lmax,1:num_clsS) =T%cl(2:lmax,1:num_clsS)    !CMB%norm(norm_As)*T%cl(2:lmax,1:num_clsS)
     if (num_cls>3 .and. num_ClsS==3) Cls(2:lmax,num_cls)=0
     if (num_cls_ext>0) then
      Cls(2:lmax,num_cls+1:num_cls_tot) =T%cl(2:lmax,num_clsS+1:num_clsS+num_cls_ext)   
     end if

     i = amp_ratio_index !this convolution is to avoid compile-time bounds-check errors on CMB%norm
     if (CMB%InitPower(amp_ratio_index) /= 0) then
        Cls(2:lmax_tensor,1:num_cls) =  Cls(2:lmax_tensor,1:num_cls)+ T%cl_tensor(2:lmax_tensor,1:num_cls) 
         !CMB%norm(norm_As)*CMB%norm(norm_amp_ratio)*T%cl_tensor(2:lmax_tensor,:)
     end if 

   end subroutine ClsFromTheoryData

   subroutine WriteTextCls(aname,T, CMB)
     Type(TheoryPredictions) T
     Type(CMBParams) CMB
     character (LEN=*), intent(in) :: aname
     integer l
     real(mcp) Cls(lmax,num_cls_tot)

     call ClsFromTheoryData(T,CMB,Cls)
     open(unit = tmp_file_unit, file = aname, form='formatted', status = 'replace')
     do l=2, lmax
        write (tmp_file_unit,*) l, Cls(l,1:num_cls)*l*(l+1)/(2*pi), Cls(l,num_cls+1:num_cls_tot)
     end do
     close(tmp_file_unit)

   end subroutine WriteTextCls

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
