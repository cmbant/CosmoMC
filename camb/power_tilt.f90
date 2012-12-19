!This module provides the initial power spectra, parameterized as an expansion in ln k
!
! ln P = ln A_s + (n_s -1)*ln(k/k_0) + n_{run}/2 * ln(k/k_0)^2 
!
! so if n_run = 0
!
! P = A_s (k/k_0_scalar)^(n_s-1)
!
!for the scalar spectrum, when an(in) is the in'th spectral index. k_0_scalar
!is a pivot scale, fixed here to 0.05/Mpc (change it below as desired).
!
!This module uses the same inputs an(in), ant(in) and rat(in) as CMBFAST, however here
!rat(in) is used to set the ratio of the initial power spectra, so here
!
!** rat(in) is not the Cl quadrupole ratio ***
!
!in general models the quadrupole ratio depends in a complicated way on the ratio of the initial
!power spectra

!The absolute normalization of the Cls is unimportant here, but the relative ratio
!of the tensor and scalar Cls generated with this module will be correct for general models


!The OutputNormalization parameter controls the final output
!Absolute Cls can be obtained by setting OuputNormalization=outNone, otherwise the overall normalization
!of the power spectra doesn't matter

!This version December 2003 - changed default tensor pivot to 0.05 (consistent with CMBFAST 4.5)

     module InitialPower   
     use Precision
     implicit none   

     private
    
      character(LEN=*), parameter :: Power_Name = 'power_tilt_modpk'

      integer, parameter :: nnmax= 5 
       !Maximum possible number of different power spectra to use

      Type InitialPowerParams

       integer nn  !Must have this variable
      !The actual number of power spectra to use
  
      !For the default implementation return power spectra based on spectral indices
       real(dl) an(nnmax) !scalar spectral indices
       real(dl) n_run(nnmax) !running of spectral index 
       real(dl) ant(nnmax) !tensor spectral indices
       real(dl) rat(nnmax) !ratio of scalar to tensor initial power spectrum amplitudes
       real(dl) k_0_scalar, k_0_tensor
       real(dl) ScalarPowerAmp(nnmax)
 
      end Type InitialPowerParams

      real(dl) curv  !Curvature contant, set in InitializePowers     
      
      Type(InitialPowerParams) :: P
  
!Make things visible as neccessary...
 
      public InitialPowerParams, InitialPower_ReadParams, InitializePowers, ScalarPower, TensorPower
      public nnmax,Power_Descript, Power_Name, SetDefPowerParams


    contains
       
  
       subroutine SetDefPowerParams(AP)
        Type (InitialPowerParams) :: AP

         AP%nn     = 1 !number of initial power spectra
         AP%an     = 1 !scalar spectral index
         AP%n_run   = 0 !running of scalar spectral index
         AP%ant    = 0 !tensor spectra index
         AP%rat    = 1
         AP%k_0_scalar = 0.05
         AP%k_0_tensor = 0.05
         AP%ScalarPowerAmp = 1

       end subroutine SetDefPowerParams

       subroutine InitializePowers(AParamSet,acurv)
!MODIFIED P(K)
         
        
 
         use modpkparams, only : use_modpk, findiffdphi, k_pivot, &
              modpk_As, modpk_r, modpk_ns, modpk_nt, modpk_nrun, flag_do_reconstruction,  reconstruction_Nefold_limit
         USE internals, only : phi_ik, pi

         use potential
         use camb_interface, only : pk_bad, pk_initialized
         use access_modpk, ONLY: potinit, evolve, pkspline_n, pkspline_k, &
              pkspline_p, pkspline_p2der, pkspline_pt, pkspline_pt2der, &
              pkspline_kmin, pkspline_kmax
         integer ik
         real(dl) k1, k2, ps0, ps1, ps2, pt0, pt1, pt2, dlnk
!END MODIFIED P(K)
         Type (InitialPowerParams) :: AParamSet
         !Called before computing final Cls in cmbmain.f90
         !Could read spectra from disk here, do other processing, etc.

        real(dl) acurv


        if (AParamSet%nn > nnmax) then
           write (*,*) 'To use ',AParamSet%nn,'power spectra you need to increase'
           write (*,*) 'nnmax in power_tilt.f90, currently ',nnmax
        end if
        P = AParamSet

        curv=acurv         

!Write implementation specific code here...

!MODIFIED P(K)
        if ((use_modpk).and.(.not.pk_initialized)) then
           ! make splines for P(k)

           findiffdphi = epsilon(1.d0)


           CALL potinit()


           if (pk_bad/=0) RETURN


           !write (*,*) 'Computing initial P(k) spline'
           do ik = 1, pkspline_n
              pkspline_k(ik) = exp(pkspline_kmin+real(ik-1)*(pkspline_kmax-pkspline_kmin)/dble(pkspline_n-1))
              call evolve(pkspline_k(ik),pkspline_p(ik),pkspline_pt(ik))
              if (pk_bad/=0) RETURN
           end do

           call splinepow(pkspline_k,pkspline_p,pkspline_n,0.d0,0.d0,pkspline_p2der)
           call splinepow(pkspline_k,pkspline_pt,pkspline_n,0.d0,0.d0,pkspline_pt2der)

! compute amplitude and tilt parameters at k_pivot
           dlnk = 2.0d0
           call evolve(k_pivot,ps1,pt1)
           modpk_As = ps1
           modpk_r = pt1/ps1
           call evolve(k_pivot,ps0,pt0)
           call evolve(k_pivot*exp(-dlnk),ps1,pt1)
           call evolve(k_pivot*exp(dlnk),ps2,pt2)
           modpk_ns = 1.d0+log(ps2/ps1)/dlnk/2.d0
           modpk_nt = log(pt2/pt1)/dlnk/2.d0
           modpk_nrun = log(ps1*ps2/ps0**2)/dlnk**2
           
           
           pk_initialized = .true.

        end if
!END MODIFIED P(K)

       end subroutine InitializePowers
       

      function ScalarPower(k,in)

!MODIFIED P(K)
       USE access_modpk, ONLY: pkspline_n, pkspline_k, &
            pkspline_p, pkspline_p2der, pkspline_kmin, pkspline_kmax
       USE modpkparams, ONLY: use_modpk, flag_do_reconstruction
       USE camb_interface, ONLY : pk_bad
!END MODIFIED P(K)

       !"in" gives the index of the power to return for this k
       !ScalarPower = const for scale invariant spectrum
       !The normalization is defined so that for adiabatic perturbations the gradient of the 3-Ricci 
       !scalar on co-moving hypersurfaces receives power
       ! < |D_a R^{(3)}|^2 > = int dk/k 16 k^6/S^6 (1-3K/k^2)^2 ScalarPower(k) 
       !In other words ScalarPower is the power spectrum of the conserved curvature perturbation given by
       !-chi = \Phi + 2/3*\Omega^{-1} \frac{H^{-1}\Phi' - \Psi}{1+w}
       !(w=p/\rho), so < |\chi(x)|^2 > = \int dk/k ScalarPower(k).
       !Near the end of inflation chi is equal to 3/2 Psi.
       !Here nu^2 = (k^2 + curv)/|curv| 

       !This power spectrum is also used for isocurvature modes where 
       !< |\Delta(x)|^2 > = \int dk/k ScalarPower(k)
       !For the isocurvture velocity mode ScalarPower is the power in the neutrino heat flux.

        real(dl) ScalarPower,k, lnrat
        integer in

!MODIFIED P(K)
        integer ik
        real(dl) pk, a, b

        if (pk_bad/=0) then
           ScalarPower = 0.d0
        else
           
           if (use_modpk) then
              
              if (k.le.pkspline_k(1)) then
                 ScalarPower = pkspline_p(1)
              else
                if (k.ge.pkspline_k(pkspline_n)) then
                    ScalarPower = pkspline_p(pkspline_n)
                else
                    ik = 1+(log(k)-pkspline_kmin)*dble(pkspline_n-1)/(pkspline_kmax-pkspline_kmin)
                    a = log(k/pkspline_k(ik))/log(pkspline_k(ik+1)/pkspline_k(ik))
                    if ((a.lt.0.).or.(a.gt.1.)) then
                       write (*,*) ' Error in spline interpolation in'
                       write (*,*) ' ScalarPower (power_tilt_feature.f90)'
                       stop
                    endif

                    if (.not. flag_do_reconstruction) then
                       b = 1.-a
                       ScalarPower = a*pkspline_p(ik)+b*pkspline_p(ik+1)+ &
                         ((a**3-a)*pkspline_p2der(ik)+(b**3-b)*pkspline_p2der(ik+1))* &
                         (pkspline_k(ik+1)-pkspline_k(ik))**2/6.
                    else 
! Weird SRR models + spline can lead to negative values. 
! Therefore we use just a simple interpolation in log space.
                      ScalarPower = exp(log(pkspline_p(ik))+(log(k)-log(pkspline_k(ik))) & 
                         *(log(pkspline_p(ik+1))-log(pkspline_p(ik)))/(log(pkspline_k(ik+1))-log(pkspline_k(ik))))
                    endif
                endif
              endif
              
           else
              
              lnrat = log(k/P%k_0_scalar)
              ScalarPower=P%ScalarPowerAmp(in)*exp((P%an(in)-1)*lnrat + P%n_run(in)/2*lnrat**2)   
              
           end if
        end if
!END MODIFIED P(K)

      end function ScalarPower

      
      function TensorPower(k,in)

!MODIFIED P(K)
       USE access_modpk, ONLY: pkspline_n, pkspline_k, &
            pkspline_pt, pkspline_pt2der, pkspline_kmin, pkspline_kmax
       USE modpkparams, ONLY: use_modpk
       USE camb_interface, ONLY : pk_bad
!END MODIFIED P(K)
      
       !TensorPower= const for scale invariant spectrum
       !The normalization is defined so that
       ! < h_{ij}(x) h^{ij}(x) > = \sum_nu nu /(nu^2-1) (nu^2-4)/nu^2 TensorPower(k)
       !for a closed model
       ! < h_{ij}(x) h^{ij}(x) > = int d nu /(nu^2+1) (nu^2+4)/nu^2 TensorPower(k)
       !for an open model
       !"in" gives the index of the power spectrum to return 
       !Here nu^2 = (k^2 + 3*curv)/|curv| 


        real(dl) TensorPower,k   
        real(dl), parameter :: PiByTwo=3.14159265d0/2._dl
   
        integer in

!MODIFIED P(K)
        integer ik
        real(dl) pk, a, b

        if (pk_bad/=0) then
           TensorPower = 0.d0
        else

        if (use_modpk) then
           if (k.le.pkspline_k(1)) then
              TensorPower = pkspline_pt(1)
           else
              if (k.ge.pkspline_k(pkspline_n)) then
                 TensorPower = pkspline_pt(pkspline_n)
              else
                 ik = 1+(log(k)-pkspline_kmin)*dble(pkspline_n-1)/(pkspline_kmax-pkspline_kmin)
                 a = log(k/pkspline_k(ik))/log(pkspline_k(ik+1)/pkspline_k(ik))
                 if ((a.lt.0.).or.(a.gt.1.)) then
                    write (*,*) ' Error in spline interpolation in'
                    write (*,*) ' TensorPower (power_tilt_feature.f90)'
                    stop
                 endif
                 b = 1.-a
                 TensorPower = a*pkspline_pt(ik)+b*pkspline_pt(ik+1)+ &
                      ((a**3-a)*pkspline_pt2der(ik)+(b**3-b)*pkspline_pt2der(ik+1))* &
                      (pkspline_k(ik+1)-pkspline_k(ik))**2/6.
              endif
           endif

        else

           TensorPower=P%rat(in)*P%ScalarPowerAmp(in)*exp(P%ant(in)*log(k/P%k_0_tensor))
           if (curv < 0) TensorPower=TensorPower*tanh(PiByTwo*sqrt(-k**2/curv-3)) 

        end if
        end if
!END MODIFIED P(K)
       
      end function TensorPower

      !Get parameters describing parameterisation (for FITS file)
     function Power_Descript(in, Scal, Tens, Keys, Vals)
         character(LEN=8), intent(out) :: Keys(*)
         real(dl), intent(out) :: Vals(*)
         integer, intent(IN) :: in
         logical, intent(IN) :: Scal, Tens
         integer num, Power_Descript
         num=0
         if (Scal) then
         num=num+1
         Keys(num) = 'n_s'
         Vals(num) = P%an(in)
         num=num+1
         Keys(num) = 'n_run'
         Vals(num) = P%n_run(in)
         num=num+1
         Keys(num) = 's_pivot'
         Vals(num) = P%k_0_scalar
         end if
         if (Tens) then
         num=num+1
         Keys(num) = 'n_t'
         Vals(num) = P%ant(in)
         num=num+1
         Keys(num) = 't_pivot'
         Vals(num) = P%k_0_tensor
         if (Scal) then
           num=num+1
           Keys(num) = 'p_ratio'
           Vals(num) = P%rat(in)
         end if
         end if
         Power_Descript = num

       end  function Power_Descript
        
       subroutine InitialPower_ReadParams(InitPower, Ini, WantTensors)
          use IniFile
!MODIFIED P(K)
          use modpkparams, only : use_modpk, max_vparams, vparams, &
               vparams_num, vnderivs, &
               k_pivot, N_pivot, phi_init0, instreheat, &
               potential_choice, slowroll_infl_end, phi_infl_end, &
               reconstruction_Nefold_limit, flag_do_reconstruction, &
               k_min, k_max, modpk_physical_priors, modpk_rho_reheat, &
               modpk_w_primordial_lower, modpk_w_primordial_upper
          use CAMB_interface, only : modpkoutput
!END MODIFIED P(K)
          Type(InitialPowerParams) :: InitPower
          Type(TIniFile) :: Ini
          logical, intent(in) :: WantTensors
          integer i
          
           InitPower%k_0_scalar = Ini_Read_Double_File(Ini,'pivot_scalar',InitPower%k_0_scalar)
           InitPower%k_0_tensor = Ini_Read_Double_File(Ini,'pivot_tensor',InitPower%k_0_tensor) 
           InitPower%nn = Ini_Read_Int('initial_power_num')
           if (InitPower%nn>nnmax) stop 'Too many initial power spectra - increase nnmax in InitialPower'
           InitPower%rat(:) = 1
           do i=1, InitPower%nn

              InitPower%an(i) = Ini_Read_Double_Array_File(Ini,'scalar_spectral_index', i)
              InitPower%n_run(i) = Ini_Read_Double_Array_File(Ini,'scalar_nrun',i,0._dl) 
              
              if (WantTensors) then
                 InitPower%ant(i) = Ini_Read_Double_Array_File(Ini,'tensor_spectral_index',i)
                 InitPower%rat(i) = Ini_Read_Double_Array_File(Ini,'initial_ratio',i)
              end if              

              InitPower%ScalarPowerAmp(i) = Ini_Read_Double_Array_File(Ini,'scalar_amp',i,1._dl) 
              !Always need this as may want to set tensor amplitude even if scalars not computed
           end do

!MODIFIED P(K)
           use_modpk = Ini_Read_Logical('use_modpk',.false.)

           if (use_modpk) then
              modpkoutput=.true.
              potential_choice = Ini_Read_Int('potential_choice',1)
              vnderivs = Ini_Read_Logical('vnderivs',.false.)
              phi_init0 = Ini_Read_Double('phi_init',17.d0)
              slowroll_infl_end = Ini_Read_Logical('slowroll_infl_end',.true.)
              phi_infl_end = Ini_Read_Double('phi_infl_end',0.d0)
              instreheat = Ini_Read_Logical('instreheat',.false.)
              k_pivot = Ini_Read_Double('k_pivot',0.05d0)
              N_pivot = Ini_Read_Double('N_pivot',50.d0)
              vparams_num = Ini_Read_Int('vparams_num',max_vparams)
              if (vparams_num>max_vparams) then
                 write(*,*)
                 write(*,*) 'Error: vparams_num > max_vparams.'
                 write(*,*) 'Increase max_vparams in modpk_modules.f90.'
                 stop
              end if
              do i = 1, vparams_num
                 vparams(i) = Ini_Read_Double_Array('vparams',i,0.d0)
              end do
              modpk_physical_priors = Ini_Read_Logical('modpk_physical_priors',.false.)
              if (modpk_physical_priors) then 
                 modpk_rho_reheat = Ini_Read_Double('modpk_rho_reheat',1.d0) !units are GeV^4
                 modpk_w_primordial_lower = Ini_Read_Double('modpk_w_primordial_lower',-0.333333d0)
                 modpk_w_primordial_upper = Ini_Read_Double('modpk_w_primordial_upper',1.d0)
              endif

              IF (potential_choice .eq. 7.or. potential_choice .eq. 8) THEN
                 flag_do_reconstruction = .true.
                 k_min = Ini_Read_Double('infl_min_k',1.d-5)
                 k_max = Ini_Read_Double('infl_max_k',5.d0)
                 reconstruction_Nefold_limit = Ini_Read_Double('reconstruction_Nefold_limit',20.5d0)
              END IF
              

           end if
!END MODIFIED P(K)


       end  subroutine InitialPower_ReadParams 



!MODIFIED P(K)

! calculates array of second derivatives used by cubic spline
! interpolation. y2 is array of second derivatives, yp1 and ypn are first
! derivatives at end points.


! 
     SUBROUTINE splinepow(x,y,n,yp1,ypn,y2)
      implicit none
      INTEGER, intent(in) :: n
      REAL*8, intent(in) :: x(n), y(n), yp1, ypn
      REAL*8, intent(out) :: y2(n)
      INTEGER i,k
      REAL*8 p,qn,sig,un
      REAL*8, dimension(:), allocatable :: u

       
      Allocate(u(1:n))
      if (yp1.gt..99d30) then
        y2(1)=0.
        u(1)=0.
      else
        y2(1)=-0.5d0
        u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
      
      do i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2. 
   
        y2(i)=(sig-1.)/p
      
         u(i)=(6.*((y(i+1)-y(i))/(x(i+ &
         1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig* &
         u(i-1))/p
      end do
      if (ypn.gt..99d30) then
        qn=0.
        un=0.
      else
        qn=0.5d0
        un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)
      do k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
      end do

      Deallocate(u)
  
!  (C) Copr. 1986-92 Numerical Recipes Software, adapted.
      END SUBROUTINE splinepow

!END MODIFIED P(K)




     end module InitialPower
