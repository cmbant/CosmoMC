PROGRAM postproc_modpk
  USE camb_interface
  USE modpkparams
  USE access_modpk, ONLY: potinit, evolve
  USE potential
  USE internals
  USE ode_path
  IMPLICIT NONE
  INTEGER*4 :: i,j,l,ii
  REAL*8 :: kin, pow, powt, kmin, kmax, dlnk, pmin, pmax, k002,r002,ps002,pt002,ns002,nrun002
  REAL*8 :: pspivot, ptpivot
  REAL*8 :: ps0, ps1, ps2, pt0, pt1, pt2

  REAL*8 :: infile(83)
  CHARACTER(120) :: filename1,filename2,directory
  CHARACTER(10) :: model
  CHARACTER(10) :: reheat_mode
  INTEGER :: num_models, num_reheat

  modpkoutput=.FALSE.

  DO num_models=1,6 

     SELECT CASE (num_models)
     CASE(1)
        model='m2phi2_'
        potential_choice=1
        vparams_num=1
     CASE(2)
        model='natural_'
        potential_choice=2
        vparams_num=2
     CASE(3)
        model='lphi4_'
        potential_choice=3
        vparams_num=1
     CASE(4)
        model='neq1_'
        potential_choice=4
        vparams_num=1
     CASE(5)
        model='neq2ov3_'
        potential_choice=5
        vparams_num=1
     CASE(6)
        model='hilltop_'
        potential_choice=6
        vparams_num=2
     END SELECT

     DO num_reheat=1,3

        SELECT CASE (num_reheat)
        CASE(1)
           reheat_mode='inst'
           instreheat=.TRUE.
           modpk_physical_priors = .TRUE.
           modpk_rho_reheat = 1.d36
           modpk_w_primordial_lower=-0.333333333333333333333
           modpk_w_primordial_upper=0.333333333333333333333
        CASE(2)
           reheat_mode='reheat2a'
           instreheat=.FALSE.
           modpk_physical_priors = .TRUE.
           modpk_rho_reheat = 1.d36
           modpk_w_primordial_lower=-0.333333
           modpk_w_primordial_upper=0.333333
        CASE(3)
           reheat_mode='reheat3b'
           instreheat=.FALSE.
           modpk_physical_priors = .TRUE.
           modpk_rho_reheat = 1.d12
           modpk_w_primordial_lower=-0.333333
           modpk_w_primordial_upper=1.0
        END SELECT

        filename1='../chains/'//TRIM(model)//TRIM(reheat_mode)//'.txt'
        filename2='../chains/'//TRIM(model)//TRIM(reheat_mode)//'_postproc.txt'

        PRINT*,filename1

        findiffdphi = EPSILON(1.d0)  

        flag_do_reconstruction = .FALSE.
        vnderivs=.FALSE.
        slowroll_infl_end =.TRUE.
        phi_infl_end=0.

        kmin = LOG(5.d-5)
        kmax = LOG(5.d0)
        phi_init0 = 20.
        k_pivot=0.05d0
        k002=0.002d0

        findiffdphi = EPSILON(1.d0)  

        OPEN(10, file=filename1)
        OPEN(12, file=filename2)

        DO i=1,10000

           vparams(1)=0.0
           vparams(2)=0.0
           N_pivot=0.0

           PRINT*,'processing line=',i
           READ(10,*,END=10) (infile(j),j=1,83)

           IF(instreheat) THEN
              vparams(1)=infile(2+19)
              IF(vparams_num.GT.1) vparams(2)=infile(2+20)
           ELSE
              N_pivot=infile(2+19)
              vparams(1)=infile(2+20)
              IF(vparams_num.GT.1) vparams(2)=infile(2+21)
           ENDIF

           CALL potinit()  

           IF(pk_bad==0) THEN    
              ! compute amplitude and tilt parameters at k_pivot
              dlnk = 2.0d0
              CALL evolve(k_pivot,ps1,pt1)
              modpk_As = ps1
              modpk_r = pt1/ps1
              CALL evolve(k_pivot,ps0,pt0)
              CALL evolve(k_pivot*EXP(-dlnk),ps1,pt1)
              CALL evolve(k_pivot*EXP(dlnk),ps2,pt2)
              modpk_ns = 1.d0+LOG(ps2/ps1)/dlnk/2.d0
              modpk_nt = LOG(pt2/pt1)/dlnk/2.d0
              modpk_nrun = LOG(ps1*ps2/ps0**2)/dlnk**2
              ! compute amplitude and tilt parameters at k_star in paper
              call evolve(k002,ps002,pt002)
              r002=pt002/ps002
              CALL evolve(k002,ps0,pt0)
              CALL evolve(k002*EXP(-dlnk),ps1,pt1)
              CALL evolve(k002*EXP(dlnk),ps2,pt2)
              ns002 = 1.d0+LOG(ps2/ps1)/dlnk/2.d0
              nrun002 = LOG(ps1*ps2/ps0**2)/dlnk**2
              
              IF(vparams_num.GT.1) THEN
                 WRITE(12,'(14E27.18)') infile(1:2), ns002, r002, nrun002, LOG(1.d10*ps002), LOG(1.d10*modpk_As), modpk_ns, modpk_r, &
                      modpk_nrun, N_pivot, modpk_w, vparams(1), vparams(2)   
              ELSE
                 WRITE(12,'(13E27.18)') infile(1:2), ns002, r002, nrun002, LOG(1.d10*ps002), LOG(1.d10*modpk_As), modpk_ns, modpk_r, &
                      modpk_nrun, N_pivot, modpk_w, vparams(1)              
              ENDIF

           ENDIF

        END DO

10      CONTINUE

        CLOSE(10)
        CLOSE(12)

     END DO

  END DO

  STOP

END PROGRAM postproc_modpk
