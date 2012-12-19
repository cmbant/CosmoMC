PROGRAM driver_modpk
  USE camb_interface
  USE modpkparams
  USE access_modpk, ONLY: potinit, evolve
  USE potential
  USE internals
  USE ode_path
  IMPLICIT NONE
  INTEGER*4 :: i,j,l,ii
  real*8 :: kin, pow, powt, kmin, kmax, dlnk, pmin, pmax
  real*8 :: pspivot, ptpivot
  real*8 :: ps0, ps1, ps2, pt0, pt1, pt2

  modpkoutput=.true.
  
  modpk_physical_priors=.true.
  modpk_rho_reheat = 1.d44
  modpk_w_primordial_lower=-0.333333 
  modpk_w_primordial_upper=1.0

  potential_choice = 10


  IF (potential_choice .eq. 7 .or. potential_choice .eq. 8) THEN
     flag_do_reconstruction = .true.
  ELSE
     flag_do_reconstruction = .false.
  END IF

  vparams_num = 3
  vnderivs=.false.

  slowroll_infl_end =.true.
  instreheat=.false.
  phi_infl_end=0.

  IF (flag_do_reconstruction) THEN
 
     kmin = log(1.d-5)
     kmax=log(5.)
     

     reconstruction_Nefold_limit=20.5
     vparams(1) = -3
     vparams(2) = 0.1
     vparams(3) = -0.02
     vparams(4) = 3.2
     k_pivot=0.05
  ELSE
     kmin = log(5.d-4)
     kmax = log(5.d0)

     vparams(1) = -3.1546d0
     vparams(2) = -3.4d0
     vparams(3) = .5d0

     phi_init0 = 20.
     k_pivot=0.05d0
     N_pivot=57.5d0
  
  END IF


  findiffdphi = epsilon(1.d0)  


  CALL potinit()


  IF(pk_bad==0) THEN    
     
     DO ii=0,1500
        kin=exp(kmin+REAL(ii)*(kmax-kmin)/REAL(15000))
        CALL evolve(kin,pow,powt)
        WRITE(*,*) kin,pow,powt
     END DO
  endif

  STOP
END PROGRAM driver_modpk
