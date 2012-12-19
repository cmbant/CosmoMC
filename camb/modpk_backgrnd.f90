MODULE background_evolution
  USE camb_interface
  USE modpkparams
  USE modpk_odeint
  USE ode_path
  USE powersp
  USE potential, ONLY : pot,getH, getHdot, dVdphi, dHdphi, getEps, getHJH, getHJeps, getHJeta, getHJxi
  USE modpk_utils, ONLY : locate, polint, bderivs, rkqs
  USE internals, ONLY : PI
  IMPLICIT NONE
  PUBLIC :: backgrnd
  PUBLIC :: rec_backgrnd
  PUBLIC :: rec_init_backgrnd
CONTAINS

  SUBROUTINE rec_init_backgrnd
    IMPLICIT NONE

    CHARACTER(16) :: fmt = '(a12,es10.2)'

    REAL*8 :: ah_kmin, dalpha
    REAL*8 :: phi_init_min, phi_init_max ! SR reconstruction
    INTEGER*4 :: i, j, count_trials ! SR reconstruction
    LOGICAL :: phi_init_set ! SR reconstruction
    REAL*8 :: eps ! SR reconstruction

    count_trials=0.
    phi_init_set=.false.
    phi_init_min=0.
    phi_init_max=0.
    

    do while (.not.phi_init_set)

       count_trials=count_trials+1
       if (phi_init_max.eq.0) then
          phi_init=2.*phi_init
       else
          phi_init=0.5*(phi_init_min+phi_init_max)
       endif

!       if (modpkoutput) write(*,*) ' phi_trial, phi_min, phi_max ', phi_init, phi_init_min, phi_init_max

       phidot_sign = -(dHdphi(phi_init))/ABS(dHdphi(phi_init))   !replaced dVdphi by dHdphi
       
       !rolling towards phi=0 ?
       if (phidot_sign.ge.0) then
          phi_init_max=phi_init
       else
          ! how many efolds from phi_init to phi=0
          call rec_trial_phi_init(phi_init)
          ! unnecessary efolds
          if (Nefolds_to_zero.eq.-1) then
!             if (modpkoutput) print*, 'unnecessary efolds'
             phi_init_max=phi_init
          ! inflation ended before reaching phi=0
          else if (Nefolds_to_zero.eq.0) then
!             if (modpkoutput) print*, 'inflation ended prematurely', phi_init
             phi_init_max=phi_init             
          else
             H_star=getHJH(phi_init)
             k_star=k_pivot/exp(Nefolds_to_zero)*H_star/H0

             ! not enough efolds

             if (Nefolds_to_zero>0 .and. k_star .gt. k_min/k_start*0.99) then ! the 0.99 is a safety margin 
 !               if (modpkoutput) print*, 'not enough efolds'
                phi_init_min=phi_init             
             ! phi_init found with enough efolds 
             else if (k_star.le.k_min/k_start) then

                phi_init_set=.true.

                a_init=k_star*Mpc2Mpl/H_star
       
                lna(1:kount)=xp(1:kount)
                DO i=1,kount
                   hubarr(i)=getHJH(phiarr(i))
                   aharr(i)=LOG(a_init*EXP(lna(i))*hubarr(i))
                END DO
                nactual=kount
                
                ! find number of efolds from k_star to k_min
                ah_kmin=LOG(k_min*Mpc2Mpl)
                i=locate(aharr,ah_kmin)
                j=MIN(MAX(i-(4-1)/2,1),nactual+1-4)
                CALL polint(aharr(j:j+4), lna(j:j+4), ah_kmin,  Nefolds_to_kmin, dalpha)

                if (modpkoutput) then
                   print*, 'MODPK: Initial phi found:'
                   print*, 'MODPK: phi_star = ', phi_init
                   print*, 'MODPK: k_star = ', k_star
                   print*, 'MODPK: H_star = ', H_star
                   print*, 'MODPK: H0 = ', H0

                   print*, 'MODPK: Nefolds from k_star to k_pivot',  Nefolds_to_zero
                   print*, 'MODPK: Nefolds from k_star to k_min',  Nefolds_to_kmin

                endif

                
             end if
          end if
       endif

       if (phi_init_min.ne.0.0 .and. abs(1.d0 - phi_init_max/phi_init_min) < 1.d-10) then
           if (modpkoutput) then 
              write (*,*), 'MODPK: phi_min too close to phi_max'
              write (*,*), 'MODPK: could not find initial phi'
              write (*,*), 'MODPK: vparams = ', vparams
           endif
           pk_bad = 13
           return
       end if

          
       if (count_trials > 5000) then
          if (modpkoutput) then
             write (*,*), 'MODPK: could not find initial phi'
             write (*,*), 'MODPK: vparams = ', vparams
          end if
          pk_bad=12
          return
       endif
       
          
    end do
    

!    if (modpkoutput) then
!       write (*,*) 'eps_star:', getHJeps(phi_init)
!       write (*,*) 'eta_star:', getHJeta(phi_init)
!       write (*,*) 'xi_star:' , getHJxi(phi_init)
!       write (*,*) 'H_star:' , H_star
!    endif

    RETURN
  END SUBROUTINE rec_init_backgrnd



  SUBROUTINE rec_trial_phi_init(phi_init)
    IMPLICIT NONE

    REAL*8 :: phi_init, errorN
    CHARACTER(16) :: fmt = '(a12,es10.2)'

    INTEGER*4, PARAMETER :: BNVAR=2
    REAL*8 :: accuracy,h1,hmin,x1,x2
    REAL*8, DIMENSION(BNVAR) :: y 
    REAL*8 :: A1, A2, A3
    REAL*8 :: eps

    !
    !     Set initial conditions
    !
    !     y(1)=phi              dydx(1)=dphi/dx
    !     y(2)=dphi/dx          dydx(1)=d^2phi/dx^2

!    if (modpkoutput) write(*,*) ' try phi_init =', phi_init


    x1=0.0 !starting value
    x2=Nefold_max !ending value
    y(1)=phi_init  !phi(x1)

    if (potential_choice .eq. 7) then
       eps = vparams(1)
    else
       eps = 10.d0**vparams(1)
    end if
    A1 = sqrt(eps/2.) ! vparams(1) is epsilon_star
    A2 = vparams(2)/4. ! vparams(2) is eta_star
    A3 = vparams(3)/(12.*sqrt(2.*eps)) ! vparams(3) is xi_star
!    in the reconstruction case it can be done exactly
    y(2)=(-2*(A1 + 2*A2*y(1) + 3*A3*y(1)**2))/ &
         & (1 + A1*y(1) + A2*y(1)**2 + A3*y(1)**3) 

    !dphi/dalpha(x1) slowroll approx
!    h_init=SQRT(pot(phi_init)/6./M_Pl/M_Pl*(1.+SQRT(1.+2./3.* &
!         &     (M_Pl*dVdphi(phi_init)/pot(phi_init))**2.)))
!        y(2)=-dVdphi(phi_init)/3./h_init/h_init 


    !     
    !     Call the integrator
    !
    ode_trial_phi=.TRUE.
    ode_underflow=.FALSE.
    ode_ps_output=.FALSE.
    ode_infl_end=.TRUE.
    save_steps=.TRUE. 
    pk_bad=0

    IF(getEps(y(1),y(2)) .GT. 1.) THEN
       slowroll_start=.FALSE.
    ELSE
       slowroll_start=.TRUE.
    ENDIF

    h1 = 0.1

    dxsav=1.d-7
    accuracy=1.0d-10
    hmin=0.0 !minimum stepsize
    
    CALL odeint(y,x1,x2,accuracy,h1,hmin,bderivs,rkqs) 

    ode_trial_phi=.FALSE.

              
    IF(ode_underflow) THEN    
       write (*,*) 'ODE_UNDERFLOW!!!'
       pk_bad=4
    ENDIF

    RETURN
  END SUBROUTINE rec_trial_phi_init

  SUBROUTINE rec_backgrnd
    IMPLICIT NONE

    INTEGER*4 :: i,j
    REAL*8 :: alpha_e,dalpha,V_end,dv,dh,ep
    REAL*8 :: a_end, a_end_inst
    REAL*8 ::V_i, ph, alpha_pivot,aa,bb
    REAL*8 :: epsarr(nsteps),vv(nsteps)
    REAL*8 :: Np_last
    CHARACTER(16) :: fmt = '(a12,es10.2)'

    INTEGER*4, PARAMETER :: BNVAR=2
    REAL*8 :: accuracy,h1,hmin,x1,x2
    REAL*8, DIMENSION(BNVAR) :: y 
    REAL*8 :: A1, A2, A3
    REAL*8 :: eps

    LOGICAL :: infl_forever=.false.

    !
    !     Set initial conditions
    !
    !     y(1)=phi              dydx(1)=dphi/dx
    !     y(2)=dphi/dx          dydx(1)=d^2phi/dx^2

!    if (modpkoutput) write(*,fmt) ' phi_init =', phi_init



    x1=0.0 !starting value
    x2=Nefold_max !ending value
    y(1)=phi_init  !phi(x1)

    if (potential_choice .eq. 7) then
       eps = vparams(1)
    else
       eps = 10.d0**vparams(1)
    end if
    A1 = sqrt(eps/2.) ! vparams(1) is epsilon_star
    A2 = vparams(2)/4. ! vparams(2) is eta_star
    A3 = vparams(3)/(12.*sqrt(2.*eps)) ! vparams(3) is xi_star
    ! in the reconstruction case it can be done exactly
    y(2)=(-2*(A1 + 2*A2*y(1) + 3*A3*y(1)**2))/ &
         & (1 + A1*y(1) + A2*y(1)**2 + A3*y(1)**3) 

    !dphi/dalpha(x1) slowroll approx
    !    y(2)=-dVdphi(phi_init)/3./h_init/h_init 


    !     
    !     Call the integrator
    !
    ode_underflow=.FALSE.
    ode_ps_output=.FALSE.
    ode_infl_end=.TRUE.
    save_steps=.TRUE. 
    pk_bad=0

    IF(getEps(y(1),y(2)) .GT. 1.) THEN
       slowroll_start=.FALSE.
    ELSE
       slowroll_start=.TRUE.
    ENDIF

    h1 = 0.1
    dxsav=1.d-8
    accuracy=1.0d-10
    hmin=0.0 !minimum stepsize
    CALL odeint(y,x1,x2,accuracy,h1,hmin,bderivs,rkqs) 


    IF(.NOT. ode_underflow) THEN    
       
       a_init=k_star*Mpc2Mpl/H_star
       
       lna(1:kount)=xp(1:kount)
       phiarr(1:kount)=yp(1,1:kount)
       dphiarr(1:kount)=yp(2,1:kount)
       DO i=1,kount
          vv(i)=pot(phiarr(i))
          hubarr(i)=getHJH(phiarr(i))
          epsarr(i) = getHJeps(phiarr(i))
          aharr(i)=LOG(a_init*EXP(lna(i))*hubarr(i))
       END DO
       nactual=kount
       
       if (lna(kount).lt.reconstruction_Nefold_limit+Nefolds_to_kmin) then
          if (kount.gt.100000) then
             if (abs((yp(1,kount)-yp(1,kount-100000))/yp(1,kount)) .lt. 0.001) then
                infl_forever=.true.
                if (modpkoutput) then
                   PRINT*,'MODPK: Inflating forever.'
                endif
             endif
          endif
          if (infl_forever.eqv..false.) then
               if (modpkoutput) then
                   PRINT*,'MODPK: Did not reach the Nefold_limit', reconstruction_Nefold_limit
                   print*,'MODPK: Inflation ended at ', lna(kount) - Nefolds_to_kmin
                   print*,'MODPK: Nefolds from k_star to k_min: ',  Nefolds_to_kmin
               endif
               pk_bad=9        

          endif
       else if (exp(aharr(kount))/Mpc2Mpl .lt. k_max*k_start) then
          if (modpkoutput) then
             PRINT*,'MODPK: Inflation ended before 100*k_max exits the horizon'
             print*,'100*k_max = ', k_start*k_max
             print*,'MODPK: Inflation ended at ', exp(aharr(kount))/Mpc2Mpl
             print*,'MODPK: Nefolds from k_min: ',  lna(kount) - Nefolds_to_kmin
          endif
          pk_bad=11      
       endif
       

    ELSE
       write (*,*) 'ODE_UNDERFLOW!!!'
       pk_bad=4
       stop
          
    ENDIF



    RETURN
  END SUBROUTINE rec_backgrnd



  SUBROUTINE backgrnd
    INTEGER*4 :: i,j, rescl_count
    REAL*8 :: phi_init_trial
    REAL*8 :: alpha_e,dalpha,V_end,dv,dh,ep,fac_A,fac_B,w_primordial
    REAL*8 :: a_end, a_end_inst
    REAL*8 ::V_i, ph, H_pivot, alpha_pivot,aa,bb
    REAL*8 :: epsarr(nsteps),vv(nsteps)
    REAL*8 :: Np_last
    REAL*8 :: connection_factor = 71.2146d0
    CHARACTER(16) :: fmt = '(a12,es10.2)'
    !
    !     Set initial conditions
    !
    !     y(1)=phi              dydx(1)=dphi/dx
    !     y(2)=dphi/dx          dydx(1)=d^2phi/dx^2

    if (modpkoutput) write(*,fmt) ' phi_init =', phi_init

    phidot_sign = -(dVdphi(phi_init))/ABS(dVdphi(phi_init))

    !if < 0 field rolls from large to small
    !if > 0, the field rolls from small to large

    !Check whether phi_infl_end is in the correct direction from phi_init
    IF(.NOT.(slowroll_infl_end)) THEN
       IF (phidot_sign .GT.0 .AND. phi_init.GE.phi_infl_end) THEN
          PRINT*, 'MODPK: Initial phi is smaller than final phi.' 
          PRINT*, 'MODPK: Please check your initial conditions'
          PRINT*, 'MODPK: QUITTING'
          STOP
       ENDIF
       IF (phidot_sign .LT.0 .AND. phi_init.LE.phi_infl_end) THEN
          PRINT*, 'MODPK: Initial phi is larger than final phi.' 
          PRINT*, 'MODPK: Please check your initial conditions'
          PRINT*, 'MODPK: QUITTING'
          STOP
       ENDIF

    ENDIF

    !Try the background evolution for a user-defined phi_init and rescale if necessary
    rescl_count=0
    rescale_factor = 0.1d0

    DO   ! now checking that we can find the pivot scale inside the inflationary phase

       phi_init_trial=phi_init   ! initial phi_init computed from user-supplied routine (modpk_potential)

       CALL trial_background(phi_init_trial, alpha_e, V_end)   ! solve equations of motion, phi_init_trial replaced if needed. 

       IF ((pk_bad/=0) .OR. phi_init_trial.EQ.phi_init) EXIT   ! if phi_init works, we can stop.

       rescl_count=rescl_count+1  ! keep track of trials. Quit if we make too many attempts (should not get here)
       IF (rescl_count .EQ. 50) THEN
          pk_bad=2
          PRINT*,'MODPK: phi_init rescaling did not work after 50 tries.'
          EXIT
       END IF


       if (modpkoutput) then      ! output
          PRINT*,'MODPK: phi_init was inconsistent. Rescaling:'
          WRITE(*,fmt) ' phi_init =', phi_init_trial
       end if


       phi_init=phi_init_trial   ! replace with guess provided by trial_background, repeat

    END DO

    IF(pk_bad==0) THEN
       !Matching condition 
       V_i=pot(phi_init)
       IF (instreheat) THEN
          a_init=EXP(-connection_factor-alpha_e+LOG(M_Pl**4/(1.5d0 * V_end))/4.d0)
          Np_last = 0.5d0*N_pivot
     
          do while (abs(N_pivot-Np_last)>0.01d0)
             Np_last = N_pivot
             alpha_pivot = alpha_e-N_pivot
             i=locate(lna(1:nactual),alpha_pivot)
             j=MIN(MAX(i-(4-1)/2,1),nactual+1-4)
             CALL polint(lna(j:j+4), hubarr(j:j+4), alpha_pivot, H_pivot, dh)
             CALL polint(lna(j:j+4), phiarr(j:j+4), alpha_pivot, aa, bb)
             N_pivot = -connection_factor + LOG(M_Pl**4/(1.5d0 * V_end))/4.d0 -log(k_pivot*Mpc2Mpl/H_pivot)
          end do
          if (modpkoutput) then
             !WRITE(*,fmt) ' H_pivot =',H_pivot
             WRITE(*,fmt) ' a_end = ',EXP(alpha_e)*a_init
             WRITE(*,fmt) ' phi_pivot =', aa
             WRITE(*,fmt) ' a_pivot =',EXP(alpha_pivot)*a_init
             WRITE(*,fmt) ' N_pivot =',N_pivot
          end if
       ELSE
          alpha_pivot = alpha_e-N_pivot
          i=locate(lna(1:nactual),alpha_pivot)
          j=MIN(MAX(i-(4-1)/2,1),nactual+1-4)
          CALL polint(lna(j:j+4), hubarr(j:j+4), alpha_pivot, H_pivot, dh)
          CALL polint(lna(j:j+4), phiarr(j:j+4), alpha_pivot, aa, bb)
          a_init=k_pivot*Mpc2Mpl/H_pivot/EXP(alpha_pivot)
          if (modpkoutput) then
             !WRITE(*,fmt) ' alpha_piv =',alpha_pivot
             !WRITE(*,fmt) ' H_pivot =',H_pivot
             WRITE(*,fmt) ' a_end =',EXP(alpha_e)*a_init
             WRITE(*,fmt) ' phi_pivot =', aa
             WRITE(*,fmt) ' a_pivot =',EXP(alpha_pivot)*a_init
             WRITE(*,fmt) ' N_pivot =',alpha_e-alpha_pivot
          end if
          a_end=EXP(alpha_e)*a_init

          a_end_inst=EXP(-connection_factor+LOG(M_Pl**4/(1.5d0 * V_end))/4.d0)


          ! If we not using a "physical" prior for the post-inflationary phase, reject models for which
          ! a_end exceeds the value at needed with instant thermalization. (This rejects scenarios
          ! for which the post-inflationary universe has p> rho/3)       
          IF (.not.modpk_physical_priors .and. (a_end .GT. a_end_inst)) THEN 
             PRINT*,'MODPK: inflation ends too late with this N_pivot.'
             pk_bad=3
             RETURN
          ENDIF

          ! With a "physical prior" we reject only models for which the energy density following inflation                                                          
          ! (recalling rho_end = 1.5 V_end) is greater than the energy scale at which we assume thermalization has occurred.
          ! This rejects scenarios which would require a "phantom energy" phase (where rho increases with a(t)) in
          ! post-inflationary universe.
          IF (modpk_physical_priors .AND. ((1.5d0 *V_end)**0.25 .LT. (modpk_rho_reheat**0.25 / 2.436d18))) THEN 
             PRINT*,'MODPK: inflation ends with rho less than assumed thermalization density'
             pk_bad=6
             RETURN
          ENDIF

       END IF
       !if (modpkoutput) WRITE(*,fmt) ' a_init =', a_init

!HVP adding in the physical prior for reheating
       modpk_w = 0.0d0
       if (modpk_physical_priors) then
          if(instreheat) THEN 
             w_primordial = 1.d0/3.d0
          else 
             fac_A = LOG(modpk_rho_reheat/(1.5d0 * V_end))-4.d0*LOG(2.436d18) !modpk_rho_reheat is in GeV^4, convert to M_Pl units
             fac_B = connection_factor+LOG(a_end)-LOG(M_Pl**4/(1.5d0 * V_end))/4.d0
             w_primordial = (fac_A - 12.d0*fac_B)/(3.d0 *fac_A + 12d0*fac_B)
          endif

          IF (w_primordial .LT. modpk_w_primordial_lower .OR. w_primordial .GT. modpk_w_primordial_upper) THEN 
             PRINT*,'MODPK: failed w_primordial physical prior'
             pk_bad = 6
             RETURN
          ENDIF
          modpk_w = w_primordial

       endif
!HVP end physical reheat prior

       DO i=1,nactual
          aharr(i)=LOG(a_init*EXP(lna(i))*hubarr(i))
       END DO
    ENDIF

    RETURN
  END SUBROUTINE backgrnd

  SUBROUTINE trial_background(phi_init_trial, alpha_e, V_end)
    INTEGER*4 :: i,j
    INTEGER*4, PARAMETER :: BNVAR=2
    REAL*8 :: accuracy,h1,hmin,x1,x2
    REAL*8 :: alpha_e,dalpha,V_end,dv,ep, phi_init_trial
    REAL*8 :: ph, H_pivot, alpha_pivot,aa,bb
    REAL*8 :: epsarr(nsteps),vv(nsteps)
    REAL*8, DIMENSION(BNVAR) :: y 
    CHARACTER(16) :: fmt = '(a12,es10.2)'

    h_init=SQRT(pot(phi_init_trial)/6./M_Pl/M_Pl*(1.+SQRT(1.+2./3.* &
         &     (M_Pl*dVdphi(phi_init_trial)/pot(phi_init_trial))**2.)))  ! initial hubble parameter

    x1=0.0 !starting value for scale factor
    x2=Nefold_max !ending value for scale factor (will assume inflation future infinite if we get here)
    y(1)=phi_init_trial  !phi(x1)
    y(2)=-dVdphi(phi_init_trial)/3./h_init/h_init !dphi/dalpha(x1) slowroll approx
    !     
    !     Call the integrator
    !
    ode_underflow=.FALSE.
    ode_ps_output=.FALSE.
    ode_infl_end=.TRUE.
    save_steps=.TRUE. 
    pk_bad=0

    IF(getEps(y(1),y(2)) .GT. 1.) THEN  ! this is "hubble slow roll" epsilon, not (V'/V)^2
       slowroll_start=.FALSE.
    ELSE
       slowroll_start=.TRUE.
    ENDIF

    !guessed start stepsize
    if (potential_choice.eq.6) then
       h1 = 0.001
    else
       h1 = 0.1
    end if

    dxsav=1.d-7
    accuracy=1.0d-11
    hmin=0.0 !minimum stepsize

    CALL odeint(y,x1,x2,accuracy,h1,hmin,bderivs,rkqs)   ! does integration

    IF(.NOT. ode_underflow) THEN    
       lna(1:kount)=xp(1:kount)  ! array with log(scalefactor) values 
       phiarr(1:kount)=yp(1,1:kount)  ! phi values
       dphiarr(1:kount)=yp(2,1:kount) ! phi-dot values 
       DO i=1,kount
          vv(i)=pot(phiarr(i))        ! value of the potential
          hubarr(i)=getH(phiarr(i),dphiarr(i))
          epsarr(i) = getEps(phiarr(i),dphiarr(i)) 
       END DO
       !
       !     Determine the parameters needed for converting k(Mpc^-1) to K
       !      
       nactual=kount

       IF(slowroll_infl_end) THEN  ! we are expecting inflation ends with breakdown of slow roll.
          ep=1.d0                  ! solves for point where epsilon=1 (a''=0)
          i=locate(epsarr(1:kount),ep)
          j=MIN(MAX(i-(4-1)/2,1),nactual+1-4)
          CALL polint(epsarr(j:j+4), lna(j:j+4), ep, alpha_e, dalpha)
          CALL polint(epsarr(j:j+4), vv(j:j+4), ep, V_end, dv)
          CALL polint(epsarr(j:j+4), phiarr(j:j+4), ep, phi_infl_end, bb)
       ELSE
          ep=phi_infl_end          ! we are supplying phi_end by hand
          i=locate(phiarr(1:kount),ep)
          j=MIN(MAX(i-(4-1)/2,1),nactual+1-4)
          CALL polint(phiarr(j:j+4), lna(j:j+4), ep, alpha_e, dalpha)
          CALL polint(phiarr(j:j+4), vv(j:j+4), ep, V_end, dv)
       ENDIF
       if (modpkoutput) WRITE(*,fmt) ' phi_end =', phi_infl_end

       IF (instreheat) THEN
          !Set a plausible pivot (ignoring logarithmic V_k and V_end terms) for 
          !the smallest k likely to be requested by CAMB. 
          N_pivot=-LOG10(1.d-6)+60.
       ENDIF


       IF(alpha_e.LT.(N_pivot+20.)) THEN  ! inflation did not last long enough
          IF ((potential_choice.eq.6).and.(vparams(1)<-2.d0)) THEN ! hillop case + V~phi limit  
             phi_init_trial = phi_init*0.9d0
          ELSE
             phi_init_trial=phi_init+(phi_init-phi_infl_end)*rescale_factor
          ENDIF
          RETURN
       END IF

       IF (alpha_e.GT.N_pivot+55) THEN
          ep=alpha_e-(N_pivot+50)
          i=locate(lna(1:kount),ep)
          j=MIN(MAX(i-(4-1)/2,1),nactual+1-4)
          CALL polint(lna(j:j+4), phiarr(j:j+4), ep, aa, bb)
          phi_init_trial=aa
          RETURN
       END IF

    ELSE
       pk_bad=4
    ENDIF

    RETURN
  END SUBROUTINE TRIAL_BACKGROUND


END MODULE BACKGROUND_EVOLUTION

