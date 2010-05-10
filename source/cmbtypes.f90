!Define the data types and read/writes them to disk. Also change l_max here.

module cmbtypes
use settings
implicit none


!Number of Cls, 1 for just temperature, 3 (4) for polarization (with B)
  integer, parameter  :: num_cls  = 3

!l_max. Tensors are not computed unless compute_tensors = T in input file
!Make these multiples of 50, should be 50 more than you need accurately
  integer, parameter :: lmax = 2100, lmax_tensor = 400

!Parameters for calculating/storing the matter power spectrum
!Note that by default everything is linear
!Old mpk settings
#ifdef DR71RG
!!! BR09: Reid et al 2009 settings for the LRG power spectrum.
  integer, parameter :: num_matter_power = 300 !number of points computed in matter power spectrum
  real, parameter    :: matter_power_minkh =  0.999e-4  !minimum value of k/h to store
  real, parameter    :: matter_power_dlnkh = 0.03     !log spacing in k/h
  real, parameter    :: matter_power_maxz = 0.
  integer, parameter :: matter_power_lnzsteps = 4  ! z=0 to get sigma8 (this first entry appears to be coded in some spots in the code!!), plus 3 LRG redshifts.
#else
  integer, parameter :: num_matter_power = 74 !number of points computed in matter power spectrum
  real, parameter    :: matter_power_minkh =  0.999e-4  !1e-4 !minimum value of k/h to store
  real, parameter    :: matter_power_dlnkh = 0.143911568     !log spacing in k/h
  real, parameter    :: matter_power_maxz = 0.    !6.0
  integer, parameter :: matter_power_lnzsteps = 1 !20
#endif
!Only used in params_CMB
   real :: pivot_k = 0.05 !Point for defining primordial power spectra
   logical :: inflation_consistency = .false. !fix n_T or not

        !Note these are the interpolated/extrapolated values. The k at which matter power is computed up to 
        !by CAMB is set in CMB_Cls_xxx with, e.g. P%Transfer%kmax = 0.6 (which is enough for 2dF)

!Number of scalar-only cls
!if num_cls=4 and CMB_lensing then increased to 4 
  integer :: num_clsS=min(num_cls,3) 

  integer, parameter :: norm_As=1, norm_amp_ratio=2, norm_freq_ix = 3
  
  Type CMBParams
     real nuisance(1:num_nuisance_params)
      !unit Gaussians for experimental parameters
     real norm(1:num_norm)
      !These are fast parameters controling amplitudes, calibrations, etc.
     real InitPower(1:num_initpower)
      !These are fast paramters for the initial power spectrum
     !Now remaining (non-independent) parameters
     real omb, omc, omv, omnu, omk, omdm
     real ombh2, omch2, omnuh2, omdmh2
     real zre, nufrac
     real h, H0
     real w
     real reserved(5)
  
  end Type CMBParams

  Type CosmoTheory
     real Age, r10
     real SN_loglike, HST_loglike, BAO_loglike, reserved(1)
     real cl(lmax,num_cls), cl_tensor(lmax_tensor,num_cls) !TT, TE, EE and BB  in that order
     real sigma_8
     real matter_power(num_matter_power,matter_power_lnzsteps)
       !second index is redshifts from 0 to matter_power_maxz
       !if custom_redshift_steps = false with equal spacing in
       !log(1+z) and matter_power_lnzsteps points
       !if custom_redshift_steps = true set in mpk.f90 
     ! BR09 additions
     real mpk_nw(num_matter_power,matter_power_lnzsteps) !no wiggles fit to matter power spectrum
     real mpkrat_nw_nl(num_matter_power,matter_power_lnzsteps) !halofit run on mpk_nw
     real finalLRGtheoryPk(num_matter_power)  !! this is the quantity that enters the LRG likelihood calculation
    ! end BR09 additions
  end Type CosmoTheory

  logical, parameter ::  Old_format  = .false.
  logical, parameter :: write_all_Cls = .false. 
   !if false use CAMB's flat interpolation scheme (lossless if models are flat except near lmax when lensed)

contains


   subroutine WriteModel(i,CMB, T, like, mult)
    integer i
    real, intent(in), optional :: mult
    Type(CosmoTheory) T
    real like, amult
    Type(CMBParams) CMB
    integer j 

    if (present(mult)) then
       amult = mult
    else
       amult = 1
    end if
    
    if (Old_format) then

      stop 'old not supported'    
    else

    j = 0 !format ID
    if (write_all_cls) j=1
    write(i) j

    write(i) amult, num_matter_power, lmax, lmax_tensor, num_cls

    write(i) T%SN_loglike, T%HST_loglike, T%reserved

    write(i) like
    write(i) CMB

    write(i) T%Age, T%r10, T%sigma_8, T%matter_power

    if (write_all_cls) then
     write(i) T%cl(2:lmax,1:num_cls)
     write(i) T%cl_tensor(2:lmax_tensor,1:num_cls)
    else

       !Use interpolation scheme CAMB uses for flat models
       !If using significantly non-flat, or increasing interpolation accuracy, save all th cls instead
        write(i) T%cl(2:20,1:num_cls)
        do j=30,90,10 
         write(i) T%cl(j,1:num_cls)
        end do
        do j=110,130, 20
         write(i) T%cl(j,1:num_cls)
        end do
        do j=150,lmax, 50
         write(i) T%cl(j,1:num_cls)
        end do

        if (lmax_tensor /= 0) then
            if (lmax_tensor<150) stop 'lmax_tensor too small'
            write(i) T%cl_tensor(2:20,1:num_cls)
            do j=30,90,10 
             write(i) T%cl_tensor(j,1:num_cls)
            end do
            do j=110,130,20 
             write(i) T%cl_tensor(j,1:num_cls)
            end do
            do j=150,lmax_tensor, 50
             write(i) T%cl_tensor(j,1:num_cls)
            end do
        end if
    end if

    end if 

    if (flush_write) call FlushFile(i)

   end subroutine WriteModel

   
 subroutine ReadModel(i,CMB, T, mult, like, error)
    integer, intent(in) :: i
    integer, intent(out) :: error
    real, intent(out) :: mult
    Type(CosmoTheory) T
    real like
    Type(CMBParams) CMB
    real icl(lmax,1:num_cls)
    integer allcl,j,ind, ix(lmax)
    integer almax,almaxtensor, anumpowers, anumcls
   
    error = 0

    if (old_format) then

       stop 'old not supported'

    else

        read(i,end=100,err=100) allcl

        if (allcl/=0 .and. allcl/=1) stop 'wrong file format'

        read(i,end=100,err=100) mult,anumpowers,almax, almaxtensor, anumcls
        if (almax > lmax) stop 'reading file with larger lmax'
        if (anumcls > num_cls) stop 'reading file with more Cls (polarization)'

        read(i) T%SN_loglike, T%HST_loglike,T%reserved
   
        read(i,end = 100, err=100) like
        read(i) CMB

    
        read(i) T%Age, T%r10, T%sigma_8, T%matter_power(1:anumpowers,1:matter_power_lnzsteps)
        T%cl = 0
        T%cl_tensor = 0
         
        if(allcl==1) then  
         read(i) T%cl(2:almax,1:anumcls)
         read(i) T%cl_tensor(2:almaxtensor,1:anumcls)
        else

            read(i) icl(1:19,1:num_cls)
            ind =1
            do j =2,20 
               ix(ind)=j
               ind=ind+1
            end do
            do j=30,90,10 
             read(i) icl(ind,1:num_cls)
             ix(ind) = j
             ind = ind + 1
            end do
             do j=110,130,20 
             read(i) icl(ind,1:num_cls)
             ix(ind) = j
             ind = ind + 1
            end do
            do j=150,almax, 50
             read(i) icl(ind,1:num_cls)
             ix(ind) = j
             ind = ind+1
            end do
            ind = ind-1
 
           call InterpCls(ix,icl, T%cl, ind, almax, num_Cls)

           if (almaxtensor /= 0) then     
            read(i) icl(1:19,1:num_cls)
            ind =1
            do j =2,20 
               ix(ind)=j
               ind=ind+1
            end do
            do j=30,90,10 
             read(i) icl(ind,1:num_cls)
             ix(ind) = j
             ind = ind + 1
            end do
            do j=110,130,20 
             read(i) icl(ind,1:num_cls)
             ix(ind) = j
             ind = ind + 1
            end do
            do j=150,almaxtensor, 50
             read(i) icl(ind,1:num_cls)
             ix(ind) = j
             ind = ind+1
            end do
            ind = ind-1
            call InterpCls(ix,icl, T%cl_tensor, ind, almaxtensor,num_cls)
           end if
        
        end if 

        return
    100 error = 1

    end if

   end subroutine ReadModel

      subroutine InterpCls(l,iCl, all_Cl, n, almax, ncls)
      integer, intent(in) :: n, almax,ncls
   
      real, intent(in) :: iCl(lmax,1:num_cls)
      integer l(n),p
      real all_Cl(:,:)
   
      integer il,llo,lhi,xi
      real xl(n),  ddCl(n)

      real a0,b0,ho
      real inCl(n)

    
      do p =1, ncls

        do il=1,n
           inCl(il) = iCl(il,p)*l(il)**2
        end do

        xl = l
        call spline_real(xl,inCl,n,ddCl)
     
        llo=1
        do il=2,l(n)
           xi=il
           if ((xi > l(llo+1)).and.(llo < n)) then
              llo=llo+1
           end if
           lhi=llo+1
           ho=l(lhi)-l(llo)
           a0=(xl(lhi)-xi)/ho
           b0=(xi-xl(llo))/ho
          
           all_Cl(il,p) = (a0*inCl(llo)+ b0*inCl(lhi)+((a0**3-a0)* ddCl(llo) &
                   +(b0**3-b0)*ddCl(lhi))*ho**2/6)/il**2
  
        end do

      end do

      all_Cl(l(n)+1:almax,:) = 0
       

      end subroutine InterpCls
  
  
   subroutine ClsFromTheoryData(T, CMB, Cls)
     Type(CosmoTheory) T
     Type(CMBParams) CMB
     real Cls(lmax,num_cls)
     integer i
     
     Cls(2:lmax,1:num_clsS) =T%cl(2:lmax,1:num_clsS)    !CMB%norm(norm_As)*T%cl(2:lmax,1:num_clsS)
     if (num_cls>3 .and. num_ClsS==3) Cls(2:lmax,num_cls)=0

     i = norm_amp_ratio !this convolution is to avoid compile-time bounds-check errors on CMB%norm
     if (CMB%norm(i) /= 0) then
        Cls(2:lmax_tensor,:) =  Cls(2:lmax_tensor,:)+ T%cl_tensor(2:lmax_tensor,:) 
         !CMB%norm(norm_As)*CMB%norm(norm_amp_ratio)*T%cl_tensor(2:lmax_tensor,:)
     end if 

   end subroutine ClsFromTheoryData

   subroutine WriteTextCls(aname,T, CMB)
     Type(CosmoTheory) T
     Type(CMBParams) CMB
     character (LEN=*), intent(in) :: aname
     integer l
     real Cls(lmax,num_cls)

     call ClsFromTheoryData(T,CMB,Cls)
     open(unit = tmp_file_unit, file = aname, form='formatted', status = 'replace')
     do l=2, lmax
        write (tmp_file_unit,*) l, Cls(l,:)*l*(l+1)/(2*pi)
     end do
     close(tmp_file_unit)

   end subroutine WriteTextCls

   function MatterPowerAt(T,kh)
     !get matter power spectrum today at kh = k/h by interpolation from stored values
     real, intent(in) :: kh
     Type(CosmoTheory) T
     real MatterPowerAt
     real x, d
     integer i
   
     x = log(kh/matter_power_minkh) / matter_power_dlnkh
     if (x < 0 .or. x >= num_matter_power-1) then
        write (*,*) ' k/h out of bounds in MatterPowerAt (',kh,')'
        stop 
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
     real, intent(in) :: kh
     Type(CosmoTheory) T
     real LRGPowerAt
     real x, d
     integer i
   
     x = log(kh/matter_power_minkh) / matter_power_dlnkh
     if (x < 0 .or. x >= num_matter_power-1) then
        write (*,*) ' k/h out of bounds in MatterPowerAt (',kh,')'
        stop 
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

     real, intent(in) :: kh
     Type(CosmoTheory) T
     real MatterPowerAt_Z
     real x, d, z, y, dz, mup, mdn
     real matter_power_dlnz
     integer i, iz
   
     matter_power_dlnz = log(matter_power_maxz+1) / (matter_power_lnzsteps -1 + 1e-13)
     y = log(1.+ z) / matter_power_dlnz 

     if (z > matter_power_maxz ) then
        write (*,*) ' z out of bounds in MatterPowerAt_Z (',z,')'
        stop
     end if
     x = log(kh/matter_power_minkh) / matter_power_dlnkh
     if (x < 0 .or. x >= num_matter_power-1) then
        write (*,*) ' k/h out of bounds in MatterPowerAt_Z (',kh,')'
        stop
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
