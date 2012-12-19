MODULE access_modpk
  USE camb_interface
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: potinit, evolve

! Number of k values for computing spline of P(k).
! Set to lower numbers for smoother potentials.

#ifdef WIGGLY
  INTEGER*4, PARAMETER, PUBLIC :: pkspline_n = 2000
#else
  INTEGER*4, PARAMETER, PUBLIC :: pkspline_n = 500
#endif

  REAL*8, PUBLIC :: pkspline_kmin = log(1.d-5), pkspline_kmax = log(5.d0)
  REAL*8, PUBLIC :: pkspline_k(pkspline_n), pkspline_p(pkspline_n), &
	pkspline_p2der(pkspline_n), pkspline_pt(pkspline_n), &
	pkspline_pt2der(pkspline_n)

CONTAINS

  SUBROUTINE potinit()
    USE modpkparams
    USE powersp
    USE background_evolution, ONLY : backgrnd, rec_backgrnd, rec_init_backgrnd
    USE potential, ONLY : initialphi,getHJH
    USE ode_path
    USE internals, ONLY : PI
    IMPLICIT NONE

    REAL*8 :: k,klo,khi
    REAL*8 :: c,gamma,A_s
    REAL*8 :: ps, pt


    if(potential_choice.eq.10) then ! monodromy
       k_start=2d2
       eval_ps=2d2
    else
       k_start=1d2
       eval_ps=1d2
    endif
 
    !
    !     Solve the background equations
    !
    pk_bad = 0

    phiarr(:) = 0
    aharr(:) = 0
    hubarr(:) = 0

    ! SRR reconstruction
    IF (potential_choice .EQ. 7.or. potential_choice .eq. 8) THEN
       flag_do_reconstruction=.true.

       pkspline_kmin = log(k_min)
       pkspline_kmax = log(k_max)

       gamma = 0.5772156649015d0
       c = -2. + log(2.) + gamma

       A_s = EXP(vparams(4))/1.d10

       IF (potential_choice .eq. 7) then ! only options 7 and 8 for now
          H0 = sqrt(8*PI**2*vparams(1)*A_s)/abs(1-(2*c+1)*vparams(1)+c*vparams(2))
       else If (potential_choice .eq.8) then
          H0 = sqrt(8*PI**2*10.d0**vparams(1)*A_s)/abs(1-(2*c+1)*10.d0**vparams(1)+c*vparams(2))
       end IF

       phi_init=0.1

       CALL rec_init_backgrnd
   
       if (pk_bad==0) then
          CALL rec_backgrnd  
!          if (pk_bad==0) then
!             call evolve(k_min,ps,pt)
!             print *,"MODPK: k_min, P(k_min)", k_min, ps
!             call evolve(k_max,ps,pt)            
!             print *,"MODPK: k_max, P(k_max)", k_max, ps
!          endif        
       endif
   
    ELSE
       phi_init = initialphi(phi_init0)

       CALL backgrnd
       
    ENDIF
    
    RETURN
  END SUBROUTINE potinit

  SUBROUTINE evolve(kin, pow, powt)
    USE modpk_odeint
    USE ode_path
    USE modpkparams
    USE internals
    USE powersp
    USE potential, ONLY: pot,powerspectrum, dVdphi, getH, getHdot
    USE modpk_utils, ONLY : locate, polint, derivs, rkqs
    IMPLICIT NONE
    INTEGER*4, PARAMETER :: NVAR=10
    INTEGER*4 :: i,j
    REAL*8 :: accuracy,h1,hmin,x1,x2 
    REAL*8, DIMENSION(NVAR) :: y 
    REAL*8, INTENT(IN) :: kin
    REAL*8, INTENT(OUT) :: pow, powt
    REAL*8 :: dum,ah,alpha_ik,dalpha,dh,p_ik,delp,dp_ik,deldp
    REAL*8 :: A1,A2,A3
    REAL*8 :: eps
    !
    !     Set initial conditions
    !
    !     y(1)=phi              dydx(1)=dphi/dalpha
    !     y(2)=dphi/dalpha      dydx(2)=d^2phi/dalpha^2
    !     y(3)=u1               dydx(3)=du1/dalpha
    !     y(4)=du1/dalpha       dydx(4)=d^2u1/dalpha^2
    !     y(5)=u2               dydx(5)=du2/dalpha
    !     y(6)=du2/dalpha       dydx(6)=d^2u2/dalpha^2


    k=kin*Mpc2Mpl

    ah=LOG(k/k_start)     
    i= locate(aharr, ah)
    IF(i.eq.0..OR. (aharr(nactual)<ah)) THEN
       PRINT*,'MODPK: The background solution worked, but the k you requested is outside'
       PRINT*,'MODPK: the bounds of the background you solved for. Please reconsider'
       PRINT*,'MODPK: your phi_init and N_pivot combo.'
       PRINT*,'MODPK: k = ', kin
       PRINT*,'MODPK: QUITTING'
       PRINT*, vparams
       STOP
    END IF
    j=MIN(MAX(i-(4-1)/2,1),nactual+1-4)

    CALL polint(aharr(j:j+4), phiarr(j:j+4), ah, p_ik, delp)
    CALL polint(aharr(j:j+4), lna(j:j+4), ah,  alpha_ik, dalpha)
    CALL polint(aharr(j:j+4), hubarr(j:j+4), ah,  h_ik, dh)
    a_ik=EXP(alpha_ik)*a_init
    x1=alpha_ik
    IF(x1.lt.0.) THEN
       PRINT*,'MODPK: The phi_init you specified is too small to give'
       PRINT*,'MODPK: sufficient efolds of inflation. We cannot self-consistently'
       PRINT*,'MODPK: solve this for you. Please adjust phi_init and try again.'
       PRINT*,'MODPK: QUITTING'
       STOP
    END IF

   IF(delp .GT. 0.1 .OR. dalpha .GT. 0.1 .OR. dh .GT. 0.1) THEN
       PRINT*,'MODPK: The interpolation in SUBROUTINE EVOLVE has suspiciously large'
       PRINT*,'MODPK: errors. Your model smells fishy.'
       PRINT*,'MODPK: QUITTING'
       print*, vparams
       STOP
    ENDIF

    y(1)=p_ik             !phi(x1)


    IF(flag_do_reconstruction) THEN
      if(potential_choice .eq. 7) then !only choices 7 and 8 for now
          eps = vparams(1)
       else if (potential_choice .eq. 8) then
          eps = 10.d0**vparams(1)
       end if
       A1 = sqrt(eps/2.) ! vparams(1) is epsilon_star
       A2 = vparams(2)/4. ! vparams(2) is eta_star
       A3 = vparams(3)/(12.*sqrt(2.*eps)) ! vparams(3) is xi_star
       ! in the reconstruction case it can be done exactly
       y(2)=(-2*(A1 + 2*A2*y(1) + 3*A3*y(1)**2))/ &
            & (1 + A1*y(1) + A2*y(1)**2 + A3*y(1)**3) 
    ELSE
       CALL polint(aharr(j:j+4), dphiarr(j:j+4), ah, dp_ik, deldp)
       y(2) = dp_ik
       
       IF(deldp .gt. 0.1) THEN
          PRINT*,'MODPK: The interpolation in SUBROUTINE EVOLVE has suspiciously large'
          PRINT*,'MODPK: errors. Your model smells fishy.'
          PRINT*,'MODPK: QUITTING'
          print*, vparams
          STOP
       ENDIF
        

    ENDIF

    y(3)=1.               !u1(x1)
    y(4)=0.               !du1/dalpha(x1)
    y(5)=0.               !u2(x1)
    y(6)=1.               !du2/dalpha(x2)
    y(7)=1.               !v1(x1)
    y(8)=0.               !dv1/dalpha(x1)
    y(9)=0.               !v2(x1)
    y(10)=1.              !dv2/dalpha(x2)
    !     
    !     Call the integrator 
    !
    ode_underflow=.FALSE.
    ode_ps_output=.TRUE.
    ode_infl_end=.FALSE.
    save_steps=.FALSE. 
    pk_bad=0
    x2=Nefold_max !ending value 
    h1=0.1 !guessed start stepsize
    accuracy=5.0d-6 !has a big impact on the speed of the code
    hmin=0.0 !minimum stepsize
    CALL odeint(y,x1,x2,accuracy,h1,hmin,derivs,rkqs) 
    IF(.NOT. ode_underflow) THEN 
       pow=pow_ik
       powt=powt_ik
    ELSE
       pow=0.
       powt=0.
       pk_bad=1
    ENDIF

    RETURN
  END SUBROUTINE evolve

END MODULE access_modpk
