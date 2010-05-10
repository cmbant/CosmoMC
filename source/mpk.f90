!Module storing observed matter power spectrum datasets, their points and window functions
!and routines for computing the likelihood

!This code is based on that in cmbdata.f90 
!and on Sam Leach's incorporation of Max Tegmark's SDSS code
!
!Originally SLB Sept 2004
!AL April 2006: added covariance matrix support (following 2df 2005)
!LV_06 : incorporation of LRG DR4 from Tegmark et al . astroph/0608632 
!AL: modified LV SDSS to do Q and b^2 or b^2*Q marge internally as for 2df

module mpk
use settings
use cmbtypes
implicit none

 Type mpkdataset
    logical :: use_set
    integer :: num_mpk_points_use ! total number of points used (ie. max-min+1)
    integer :: num_mpk_kbands_use ! total number of kbands used (ie. max-min+1)
    character(LEN=20) :: name
    real, pointer, dimension(:,:) :: N_inv
    real, pointer, dimension(:,:) :: mpk_W, mpk_invcov
    real, pointer, dimension(:) :: mpk_P, mpk_sdev, mpk_k
    logical :: use_scaling !as SDSS_lrgDR3
   !for Q and A see e.g. astro-ph/0501174, astro-ph/0604335
    logical :: Q_marge, Q_flat
    real :: Q_mid, Q_sigma, Ag
   end Type mpkdataset

 integer :: num_mpk_datasets = 0
 Type(mpkdataset) mpkdatasets(10)

 !Note all units are in k/h here
 
  logical :: use_mpk = .false.
  
contains 

  subroutine ReadmpkDataset(gname)   
    character(LEN=*), intent(IN) :: gname
    character(LEN=120) :: kbands_file, measurements_file, windows_file, cov_file
    Type (mpkdataset) :: mset

    integer i,iopb
    real keff,klo,khi,beff
    integer :: num_mpk_points_full ! actual number of bandpowers in the infile
    integer :: num_mpk_kbands_full ! actual number of k positions " in the infile
    integer :: max_mpk_points_use ! in case you don't want the smallest scale modes (eg. sdss)
    integer :: min_mpk_points_use ! in case you don't want the largest scale modes
    integer :: max_mpk_kbands_use ! in case you don't want to calc P(k) on the smallest scales (will truncate P(k) to zero here!)
    integer :: min_mpk_kbands_use ! in case you don't want to calc P(k) on the largest scales (will truncate P(k) to zero here!)
    real, dimension(:,:), allocatable :: mpk_Wfull, mpk_covfull
    real, dimension(:), allocatable :: mpk_kfull, mpk_fiducial

    character(80) :: dummychar
    logical bad
 
    num_mpk_datasets = num_mpk_datasets + 1
    if (num_mpk_datasets > 10) stop 'too many datasets'

    call Ini_Open(gname, 1, bad, .false.)
    if (bad) then
      write (*,*)  'Error opening dataset file '//trim(gname)
      stop
    end if

    mset%name = Ini_Read_String('name') 
    mset%use_set =.true.
    if (Feedback > 0) write (*,*) 'reading: '//trim(mset%name)
    num_mpk_points_full = Ini_Read_Int('num_mpk_points_full',0)
    if (num_mpk_points_full.eq.0) write(*,*) ' ERROR: parameter num_mpk_points_full not set'
    num_mpk_kbands_full = Ini_Read_Int('num_mpk_kbands_full',0)
    if (num_mpk_kbands_full.eq.0) write(*,*) ' ERROR: parameter num_mpk_kbands_full not set'
    min_mpk_points_use = Ini_Read_Int('min_mpk_points_use',1)
    min_mpk_kbands_use = Ini_Read_Int('min_mpk_kbands_use',1)
    max_mpk_points_use = Ini_Read_Int('max_mpk_points_use',num_mpk_points_full)
    max_mpk_kbands_use = Ini_Read_Int('max_mpk_kbands_use',num_mpk_kbands_full)
    mset%num_mpk_points_use = max_mpk_points_use - min_mpk_points_use +1
    mset%num_mpk_kbands_use = max_mpk_kbands_use - min_mpk_kbands_use +1

    allocate(mpk_Wfull(num_mpk_points_full,num_mpk_kbands_full))
    allocate(mpk_kfull(num_mpk_kbands_full))
    allocate(mset%mpk_P(mset%num_mpk_points_use))
    allocate(mset%mpk_sdev(mset%num_mpk_points_use))  ! will need to replace with the covmat
    allocate(mset%mpk_k(mset%num_mpk_kbands_use))
    allocate(mset%mpk_W(mset%num_mpk_points_use,mset%num_mpk_kbands_use))
    allocate(mpk_fiducial(mset%num_mpk_points_use))

    kbands_file  = Ini_Read_String('kbands_file')
    call ReadVector(kbands_file,mpk_kfull,num_mpk_kbands_full)
    mset%mpk_k(1:mset%num_mpk_kbands_use)=mpk_kfull(min_mpk_kbands_use:max_mpk_kbands_use) 
    if (Feedback > 1) then 
       write(*,*) 'reading: ',mset%name,' data'
       write(*,*) 'Using kbands windows between',mset%mpk_k(1),' < k/h < ',mset%mpk_k(mset%num_mpk_kbands_use)      
    endif
    if  (mset%mpk_k(1) < matter_power_minkh) then
       write (*,*) 'WARNING: k_min in '//trim(mset%name)//'less than setting in cmbtypes.f90'
       write (*,*) 'all k<matter_power_minkh will be set to matter_power_minkh' 
    end if

    measurements_file  = Ini_Read_String('measurements_file')
    call OpenTxtFile(measurements_file, tmp_file_unit)
    mset%mpk_P=0.
    read (tmp_file_unit,*) dummychar
    read (tmp_file_unit,*) dummychar
    do i= 1, (min_mpk_points_use-1)
       read (tmp_file_unit,*, iostat=iopb) keff,klo,khi,beff,beff,beff
    end do
    if (Feedback > 1 .and. min_mpk_points_use>1) write(*,*) 'Not using bands with keff=  ',keff,' or below'
    do i =1, mset%num_mpk_points_use
       read (tmp_file_unit,*, iostat=iopb) keff,klo,khi,mset%mpk_P(i),mset%mpk_sdev(i),mpk_fiducial(i)
    end do
    close(tmp_file_unit) 
    if (Feedback > 1) write(*,*) 'bands truncated at keff=  ',keff
    
    windows_file  = Ini_Read_String('windows_file')
    if (windows_file.eq.'') write(*,*) 'ERROR: mpk windows_file not specified'
    call ReadMatrix(windows_file,mpk_Wfull,num_mpk_points_full,num_mpk_kbands_full)
    mset%mpk_W(1:mset%num_mpk_points_use,1:mset%num_mpk_kbands_use)= &
       mpk_Wfull(min_mpk_points_use:max_mpk_points_use,min_mpk_kbands_use:max_mpk_kbands_use)
    

    cov_file  = Ini_Read_String('cov_file')
    if (cov_file /= '') then
     allocate(mpk_covfull(num_mpk_points_full,num_mpk_points_full))
     call ReadMatrix(cov_file,mpk_covfull,num_mpk_points_full,num_mpk_points_full)
     allocate(mset%mpk_invcov(mset%num_mpk_points_use,mset%num_mpk_points_use))
     mset%mpk_invcov=  mpk_covfull(min_mpk_points_use:max_mpk_points_use,min_mpk_points_use:max_mpk_points_use)
     call Matrix_Inverse(mset%mpk_invcov)
     deallocate(mpk_covfull)
    else
     nullify(mset%mpk_invcov)
    end if

    mset%use_scaling = Ini_Read_Logical ('use_scaling',.false.)

    mset%Q_marge = Ini_Read_Logical('Q_marge',.false.)
    if (mset%Q_marge) then
     mset%Q_flat = Ini_Read_Logical('Q_flat',.false.)
     if (.not. mset%Q_flat) then
      !gaussian prior on Q
      mset%Q_mid = Ini_Read_Real('Q_mid')
      mset%Q_sigma = Ini_Read_Real('Q_sigma')
     end if
     mset%Ag = Ini_Read_Real('Ag', 1.4)
    end if 
    if (iopb.ne.0) then
       stop 'Error reading mpk file'
    endif
 
   call Ini_Close

   deallocate(mpk_Wfull, mpk_kfull,mpk_fiducial)

   mpkdatasets(num_mpk_datasets) = mset
 
  end subroutine ReadmpkDataset

 
  function LSS_mpklike(Theory,mset,CMB) result(LnLike) ! LV_06 added CMB here
   Type (mpkdataset) :: mset
   Type (CosmoTheory) Theory
   Type(CMBparams) CMB     !LV_06 added for LRGDR4
   real LnLike
   real, dimension(:), allocatable :: mpk_Pth, mpk_k2,mpk_lin,k_scaled !LV_06 added for LRGDR4
   real, dimension(:), allocatable :: w
   real, dimension(:), allocatable :: mpk_WPth, mpk_WPth_k2
   real :: covdat(mset%num_mpk_points_use), covth(mset%num_mpk_points_use),  covth_k2(mset%num_mpk_points_use)
   real :: normV, Ag, Q, minchisq
   real :: a_scl  !LV_06 added for LRGDR4
   integer :: i, iQ
   logical :: do_marge
   integer, parameter :: nQ=6
   real :: tmp, dQ = 0.4
   real chisq(-nQ:nQ)
   real calweights(-nQ:nQ)
   real vec2(2),Mat(2,2)
   
   allocate(mpk_lin(mset%num_mpk_kbands_use) ,mpk_Pth(mset%num_mpk_kbands_use))
   allocate(mpk_WPth(mset%num_mpk_points_use))
   allocate(k_scaled(mset%num_mpk_kbands_use))!LV_06 added for LRGDR4
   allocate(w(mset%num_mpk_points_use))

   chisq = 0

   if (.not. mset%use_set) then
      LnLike = 0
      return
   end if

   ! won't actually want to do this multiple times for multiple galaxy pk data sets?..

   IF(mset%use_scaling) then   
      call compute_scaling_factor(dble(CMB%omk),dble(CMB%omv),dble(CMB%w),a_scl)      
   else
     a_scl = 1
   end if
      

   do i=1, mset%num_mpk_kbands_use 
     !Errors from using matter_power_minkh at lower end should be negligible
         k_scaled(i)=max(matter_power_minkh,a_scl*mset%mpk_k(i))
         mpk_lin(i)=MatterPowerAt(Theory,k_scaled(i))/a_scl**3
   end do
      
      
    do_marge = mset%Q_Marge
    if (do_marge .and. mset%Q_flat) then
        !Marginalize analytically with flat prior on b^2 and b^2*Q
        !as recommended by Max Tegmark for SDSS
          allocate(mpk_k2(mset%num_mpk_kbands_use))
          allocate(mpk_WPth_k2(mset%num_mpk_points_use))
 
          mpk_Pth=mpk_lin/(1+mset%Ag*k_scaled)
          mpk_k2=mpk_Pth*k_scaled**2
          mpk_WPth = matmul(mset%mpk_W,mpk_Pth)
          mpk_WPth_k2 = matmul(mset%mpk_W,mpk_k2)

          if (associated(mset%mpk_invcov)) then
            covdat = matmul(mset%mpk_invcov,mset%mpk_P)
            covth = matmul(mset%mpk_invcov,mpk_WPth)
            covth_k2 = matmul(mset%mpk_invcov,mpk_WPth_k2)
          else
            w=1/(mset%mpk_sdev**2)
            covdat = mset%mpk_P*w
            covth = mpk_WPth*w
            covth_k2 = mpk_WPth_k2*w     
          end if
           
          Mat(1,1) = sum(covth*mpk_WPth)
          Mat(2,2) = sum(covth_k2*mpk_WPth_k2)
          Mat(1,2) = sum(covth*mpk_WPth_k2)
          Mat(2,1) = Mat(1,2)
          LnLike = log( Mat(1,1)*Mat(2,2)-Mat(1,2)**2)
          call inv_mat22(Mat)
          vec2(1) = sum(covdat*mpk_WPth)
          vec2(2) = sum(covdat*mpk_WPth_k2)
          LnLike = (sum(mset%mpk_P*covdat) - sum(vec2*matmul(Mat,vec2)) + LnLike ) /2

          deallocate(mpk_k2,mpk_WPth_k2)      
    else

      if (mset%Q_sigma==0) do_marge = .false.
      
      do iQ=-nQ,nQ
         Q = mset%Q_mid +iQ*mset%Q_sigma*dQ 
         
         if (mset%Q_marge) then
            mpk_Pth=mpk_lin*(1+Q*k_scaled**2)/(1+mset%Ag*k_scaled)
         else 
            mpk_Pth = mpk_lin
         end if
         
         mpk_WPth = matmul(mset%mpk_W,mpk_Pth)
         
         !with analytic marginalization over normalization nuisance (flat prior on b^2)
         !See appendix F of cosmomc paper
         
         if (associated(mset%mpk_invcov)) then
            covdat = matmul(mset%mpk_invcov,mset%mpk_P)
            covth = matmul(mset%mpk_invcov,mpk_WPth)
            normV = sum(mpk_WPth*covth)
            chisq(iQ) = sum(mset%mpk_P*covdat)  - sum(mpk_WPth*covdat)**2/normV  + log(normV)  
            
         else
      
            !with analytic marginalization over normalization nuisance (flat prior on b^2)
            w=1/(mset%mpk_sdev**2)
            normV = sum(mpk_WPth*mpk_WPth*w)
            tmp=sum(mpk_WPth*mset%mpk_P*w)/normV ! avoid subtracting one large number from another
            chisq(iQ) = sum(mset%mpk_P*(mset%mpk_P - mpk_WPth*tmp)*w)  + log(normV)
         end if
         
         if (do_marge) then
            calweights(iQ) = exp(-(iQ*dQ)**2/2)
         else 
            LnLike = chisq(iQ)/2
            exit
         end if
   
      end do
      
      !without analytic marginalization
      !! chisq = sum((mset%mpk_P(:) - mpk_WPth(:))**2*w) ! uncommented for debugging purposes
      
      if (do_marge) then
         minchisq=minval(chisq)
         LnLike = sum(exp(-(chisq-minchisq)/2)*calweights)/sum(calweights)
         if (LnLike == 0) then
            LnLike = LogZero
         else
            LnLike =  -log(LnLike) + minchisq/2
         end if
      end if

   end if !not analytic over Q
      
   if (Feedback>1) write(*,*) 'mpk chi-sq:', LnLike*2
   
   if (LnLike > 1e8) then
      write(*,*) 'Chisq is huge, maybe there is a problem? chisq=',chisq
   end if
   
   deallocate(mpk_Pth,mpk_lin)
   deallocate(mpk_WPth,k_scaled,w)
   
 end function LSS_mpklike


 function LSSLnLike(CMB, Theory)
   Type (CMBParams) CMB
   Type (CosmoTheory) Theory
   real LSSLnLike
   integer i
   real tot(num_mpk_datasets)
  do i=1, num_mpk_datasets
     if (mpkdatasets(i)%name == 'twodf') then
        stop 'twodf no longer supported - use data/2df_2005.dataset'
     else
      tot(i) = LSS_mpklike(Theory,mpkdatasets(i),CMB) !LV_06 added CMB here
     end if
  end do
  LSSLnLike = SUM(tot) 
  
 end function LSSLnLike

 subroutine inv_mat22(M)
    real M(2,2), Minv(2,2), det

    det = M(1,1)*M(2,2)-M(1,2)*M(2,1)
    Minv(1,1)=M(2,2)
    Minv(2,2) = M(1,1)
    Minv(1,2) = - M(2,1)
    Minv(2,1) = - M(1,2)
    M = Minv/det

 end subroutine inv_mat22

!-----------------------------------------------------------------------------
!LV added to include lrg DR4

subroutine compute_scaling_factor(Ok,Ol,w,a)
  ! a = dV for z=0.35 relative to its value for flat Om=0.25 model.
  ! This is the factor by which the P(k) measurement would shift 
  ! sideways relative to what we got for this fiducial flat model.
  ! * a = (a_angular**2 * a_radial)**(1/3)
  ! * a_angular = comoving distance to z=0.35 in Mpc/h relative to its value for flat Om=0.25 model
  !     dA = (c/H)*eta = (2997.92458 Mpc/h)*eta, so we care only about 
  !     eta scaling, not h scaling.
  !     For flat w=-1 models, a ~ (Om/0.25)**(-0.065)
  !     For the LRG mean redshift z=0.35, the power law fit 
  !    dA(z,Om= 0.3253 (Om/0.25)^{-0.065}c H_0^{-1} is quite good within 
  !    our range of interest,
  !     accurate to within about 0.1% for 0.2<Om<0.3.
  ! * a_radial = 1/H(z) relative to its value for flat Om=0.25 model
  implicit none
  real*8 Or, Om, Ok, Ol, w, Ok0, Om0, Ol0, w0, z, eta, eta0, Hrelinv, Hrelinv0, tmp
  real*8 a_radial, a_angular
  real a
  !Or= 0.0000415996*(T_0/2.726)**4 / h**2
  Or= 0! Radiation density totally negligible at  z < 0.35
  Om= 1-Ok-Ol-Or
  z  = 0.35
  Hrelinv= 1/sqrt(Ol*(1+z)**(3*(1+w)) + Ok*(1+z)**2 + Om*(1+z)**3 + Or*(1+z)**4)
!  write(*,*) Ok,Ol,w  
call compute_z_eta(Or,Ok,Ol,w,z,eta)
  tmp = sqrt(abs(Ok))
  if (Ok.lt.-1.d-6) eta = sin(tmp*eta)/tmp
  if (Ok.gt.1d-6)   eta = (exp(tmp*eta)-exp(-tmp*eta))/(2*tmp) ! sinh(tmp*eta)/tmp
  Ok0= 0
  Ol0= 0.75
  w0= -1
  Om0= 1-Ok0-Ol0-Or
  call compute_z_eta(Or,Ok0,Ol0,w0,z,eta0)
  Hrelinv0= 1/sqrt(Ol0*(1+z)**(3*(1+w0)) + Ok0*(1+z)**2 + Om0*(1+z)**3 + Or*(1+z)**4)
  !a_angular = (Om/0.25)**(-0.065) * (-w*Otot)**0.14 ! Approximation based on Taylor expansion
  a_angular = eta/eta0
  a_radial= Hrelinv/Hrelinv0
  a=  (a_angular**2 * a_radial)**(1/3.d0)
  !write(*,'(9f10.5)') Ok,Ol,w,a,a_radial,a_angular,(Om/0.25)**(-0.065) * (-w*(1-Ok))**0.14
  !write(*,'(9f10.5)') Ok,Ol,w,a,a_radial**(2/3.d0),a_angular**(4/3.d0),((Om/0.25)**(-0.065) * (-w*(1-Ok))**0.14)**(4/3.d0)
end subroutine compute_scaling_factor

subroutine eta_demo
  implicit none
  real*8 Or, Ok, Ol, w, h, z, eta
  h  = 0.7
  Ok = 0
  Ol = 0.7
  Or = 0.0000416/h**2 
  w  = -1
  z  = 1090
  call compute_z_eta(Or,Ok,Ol,w,z,eta)
!  print *,'eta.............',eta
!  print *,'dlss in Gpc.....',(2.99792458/h)*eta
end subroutine eta_demo

!INTERFACE
logical function nobigbang2(Ok,Ol,w)
  ! Test if we're in the forbidden zone where the integrand blows up
  ! (where's H^2 < 0 for some a between 0 and 1).
  ! The function f(a) = Omega_m + Omega_k*a + Omega_l*a**(-3*w)
  ! can have at most one local minimum, at (Ok/(3*w*Ol))**(-1/(1+3*w)), 
  ! so simply check if f(a)<0 there or at the endpoints a=0, a=1.
  ! g(0) = Omega_m - Omega_l*a**(-3*w) < 0 if w > 0 & Omega_k > 1
  !                                     or if w = 0 & Omega_l < 1       
  ! g(1) = Omega_m + Omega_k + Omega_l = 1 > 0
  implicit none
  real*8 Ok, Ol, w, Om, tmp, a, epsilon
  integer failure
  failure = 0
  epsilon = 0
  !epsilon = 0.04  ! Numerical integration fails even before H^2 goes negative.
  Om = 1.d0 - Ok - Ol
  if (w*Ol.ne.0) then
     tmp = Ok/(3*w*Ol)
     if ((tmp.gt.0).and.(1+3*w.ne.0)) then ! f'(0)=0 for some a>0
        a = tmp**(-1/(1+3*w))
        if (a.lt.1) then
           if (Om + Ok*a + Ol*a**(-3*w).lt.epsilon) failure = 1
        end if
     end if
  end if
  if ((w.eq.0).and.(Ok.gt.1)) failure = 2
  if ((w.gt.0).and.(Ol.lt.0)) failure = 3
  nobigbang2 = (failure.gt.0)
  if (failure.gt.0) print *,'Big Bang failure mode ',failure
  return
end function nobigbang2
!END INTERFACE

real*8 function eta_integrand(a)
  implicit none
  real*8 Or, Ok, Ox, w
  common/eta/Or, Ok, Ox, w
  real*8 a, Om
  ! eta = int (H0/H)dz = int (H0/H)(1+z)dln(1+z) = int (H0/H)/a dlna = int (H0/H)/a^2 da = 
  ! Integrand = (H0/H)/a^2
  ! (H/H0)**2 = Ox*a**(-3*(1+w)) + Ok/a**2 + Om/a**3 + Or/a**4 
  if (a.eq.0.d0) then 
     eta_integrand = 0.d0
  else
     Om = 1.d0 - Or - Ok - Ox
     eta_integrand = 1.d0/sqrt(Ox*a**(1-3*w) + Ok*a**2 + Om*a + Or)
  end if
  return
end function eta_integrand

subroutine eta_z_integral(Omega_r,Omega_k,Omega_x,w_eos,z,eta)
  ! Computes eta as a function
  ! of the curvature Omega_k, the dark energy density Omega_x
  ! and its equation of state w.
  implicit none
  real*8 Or, Ok, Ox, w
  common/eta/Or, Ok, Ox, w
  real*8 Omega_r, Omega_k,Omega_x,w_eos, z, eta, epsabs, epsrel, amin, amax!, eta_integrand
  Or = Omega_r
  Ok = Omega_k
  Ox = Omega_x
  w  = w_eos
  epsabs  = 0
  epsrel  = 1.d-10
  amin= 1/(1+z)
  amax= 1
  call qromb2(eta_integrand,amin,amax,epsabs,epsrel,eta)
  return
end subroutine eta_z_integral

subroutine compute_z_eta(Or,Ok,Ox,w,z,eta)
  ! Computes the conformal distance eta(z)
  implicit none
  real*8 Or, Ok, Ox, w, z, eta, Ht0
!  logical nobigbang2
  if (nobigbang2(Ok,Ox,w)) then
     print *,'No big bang, so eta undefined if z>zmax.'
     eta = 99 
  else
     call eta_z_integral(Or,Ok,Ox,w,z,eta) 
     ! print *,'Or, Ok, Ox, w, z, H_0 t_0...',Or, Ok, Ox, w, eta
  end if
  return
end subroutine compute_z_eta


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! num rec routines
!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE qromb2(func,a,b,epsabs,epsrel,ss)
! The numerical recipes routine, but modified so that is decides
! it's done when either the relative OR the absolute accuracy has been attained.
! The old version used relative errors only, so it always failed when
! when the integrand was near zero.
! epsabs = epsrel = 1e-6 are canonical choices.
  INTEGER JMAX,JMAXP,K,KM
  REAL*8 a,b,func,ss,epsabs,epsrel
  EXTERNAL func
  PARAMETER (JMAX=20, JMAXP=JMAX+1, K=5, KM=K-1)
                                !    USES polint,trapzd
  INTEGER j
  REAL*8 dss,h(JMAXP),s(JMAXP)
  h(1)=1.d0
  do j=1,JMAX
     call trapzd(func,a,b,s(j),j)
     if (j.ge.K) then
        call polint(h(j-KM),s(j-KM),K,0.d0,ss,dss)
        if (abs(dss).le.epsrel*abs(ss)) return
        if (abs(dss).le.epsabs) return
     endif
     s(j+1)=s(j)
     h(j+1)=0.25d0*h(j)
  ENDDO
  print *,'Too many steps in qromb'
      
  RETURN 
END SUBROUTINE qromb2
  
SUBROUTINE polint(xa,ya,n,x,y,dy) ! From Numerical Recipes
  INTEGER n,NMAX
  REAL*8 dy,x,y,xa(n),ya(n)
  PARAMETER (NMAX=10)
  INTEGER i,m,ns
  REAL*8 den,dif,dift,ho,hp,w,c(NMAX),d(NMAX)
  ns=1
  dif=abs(x-xa(1))
  do  i=1,n
     dift=abs(x-xa(i))
     if (dift.lt.dif) then
        ns=i
        dif=dift
     endif
     c(i)=ya(i)
     d(i)=ya(i)
  enddo
  y=ya(ns)
  ns=ns-1
  do  m=1,n-1
     do  i=1,n-m
        ho=xa(i)-x
        hp=xa(i+m)-x
        w=c(i+1)-d(i)
        den=ho-hp
        if(den.eq.0.) then
           print*, 'failure in polint'
           stop
        endif
        den=w/den
        d(i)=hp*den
        c(i)=ho*den
     enddo
     if (2*ns.lt.n-m)then
        dy=c(ns+1)
     else
        dy=d(ns)
        ns=ns-1
     endif
     y=y+dy
  enddo
  return
END SUBROUTINE polint
      
SUBROUTINE trapzd(func,a,b,s,n) ! From Numerical Recipes
  INTEGER n
  REAL*8 a,b,s,func
  EXTERNAL func
  INTEGER it,j
  REAL*8 del,sum,tnm,x
  if (n.eq.1) then
     s=0.5*(b-a)*(func(a)+func(b))
  else
     it=2**(n-2)
     tnm=it
     del=(b-a)/tnm
     x=a+0.5*del
     sum=0.
     do  j=1,it
        sum=sum+func(x)
        x=x+del
     enddo
     s=0.5*(s+(b-a)*sum/tnm)
  endif
  return
END SUBROUTINE trapzd
            





end module 


