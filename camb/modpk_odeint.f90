module modpk_odeint
  IMPLICIT NONE
  CONTAINS

SUBROUTINE odeint(ystart,x1,x2,eps,h1,hmin,derivs,rkqs)
  USE ode_path
  USE internals
  USE powersp
  USE modpkparams
  USE potential
  USE modpk_utils, only : reallocate_rv, reallocate_rm

  IMPLICIT NONE
  REAL*8, DIMENSION(:), INTENT(INOUT) :: ystart
  REAL*8, INTENT(IN) :: x1,x2,eps,h1,hmin
  INTERFACE
     SUBROUTINE derivs(x,y,dydx)
       IMPLICIT NONE
       REAL*8, INTENT(IN) :: x
       REAL*8, DIMENSION(:), INTENT(IN) :: y
       REAL*8, DIMENSION(:), INTENT(OUT) :: dydx
     END SUBROUTINE derivs

     SUBROUTINE rkqs(y,dydx,x,htry,eps,yscal,hdid,hnext,derivs)
       IMPLICIT NONE
       REAL*8, DIMENSION(:), INTENT(INOUT) :: y
       REAL*8, DIMENSION(:), INTENT(IN) :: dydx,yscal
       REAL*8, INTENT(INOUT) :: x
       REAL*8, INTENT(IN) :: htry,eps
       REAL*8, INTENT(OUT) :: hdid,hnext
       INTERFACE
          SUBROUTINE derivs(x,y,dydx)
            IMPLICIT NONE
            REAL*8, INTENT(IN) :: x
            REAL*8, DIMENSION(:), INTENT(IN) :: y
            REAL*8, DIMENSION(:), INTENT(OUT) :: dydx
          END SUBROUTINE derivs
       END INTERFACE
     END SUBROUTINE rkqs
  END INTERFACE
  REAL*8, PARAMETER :: TINY=1.0d-30
  INTEGER*4, PARAMETER :: MAXSTP=1000000
  INTEGER*4 :: nstp,i
  REAL*8 :: h,hdid,hnext,x,xsav
  REAL*8, DIMENSION(SIZE(ystart)) :: dydx,y,yscal
  REAL*8 :: z, scalefac
  REAL*8 :: H_pivot

  REAL*8 :: phi_last=1
  REAL*8 :: ah_last=1

  ode_underflow=.FALSE.
  infl_ended=.FALSE.
  x=x1
  h=SIGN(h1,x2-x1)
  nok=0
  nbad=0
  kount=0
  y(:)=ystart(:)
  IF (associated(xp)) DEALLOCATE(xp)
  IF (associated(yp)) DEALLOCATE(yp)
  IF (save_steps) THEN
     xsav=x-2.d0*dxsav
     ALLOCATE(xp(2048))
     ALLOCATE(yp(SIZE(ystart),SIZE(xp)))
  END IF
  DO nstp=1,MAXSTP

!     write (*,*) 'Nefold, phi, V(phi)= ' , x , y(1), pot(y(1))
     
     ! SR reconstruction
     if (flag_do_reconstruction .AND. ode_infl_end) then

        if (ode_trial_phi) then
        
           !inflation ended before reaching phi=0
           if (getHJeps(y(1)) .GT. 1 ) then
              Nefolds_to_zero=0
              return
           endif
           
           !unnecessary efolds
           H_star=getHJH(phi_init)
           k_star=k_pivot/exp(x)*H_star/H0
           if (k_star .lt. 0.1*k_min/k_start) then
              Nefolds_to_zero=-1
              return
           end if


           if (y(1).le.0) then
              ystart(:)=y(:)
              IF (save_steps) CALL save_a_step
              Nefolds_to_zero=x

!          Interpolation to get a more accurate Nefolds_to_zero:
!          Nefolds_to_zero = (xp(kount) - xp(kount-1))/(yp(1,kount)-yp(1,kount-1)) &
!               & *(0.d0-yp(1,kount-1))+ xp(kount-1)

              return
           endif
        
           ! compute background until the number of efolds is reached
        else
           if (x.ge.reconstruction_Nefold_limit+Nefolds_to_kmin) then
!              print *, 'reached the number of specified efolds'

              if (k_star/H_star*EXP(x)*getHJH(y(1)) .ge. k_max*k_start) then
                 ystart(:)=y(:)
                 IF (save_steps) CALL save_a_step
                 return
              endif
           endif
           !        if (getEps(y(1),y(2)) .GT. 1 ) then
           if (getHJeps(y(1)) .GT. 1 ) then
!              PRINT*,'MODPK: The model did not reach the specified number'
!              PRINT*,'MODPK: of e-folds, ended at: ', x
              ystart(:)=y(:)
              IF (save_steps) CALL save_a_step
              return
           endif
        endif
     endif



     CALL derivs(x,y,dydx)
     yscal(:)=ABS(y(:))+ABS(h*dydx(:))+TINY
     IF (save_steps .AND. (ABS(x-xsav) > ABS(dxsav))) &
          CALL save_a_step
     IF ((x+h-x2)*(x+h-x1) > 0.0) h=x2-x
     CALL rkqs(y,dydx,x,h,eps,yscal,hdid,hnext,derivs)
     IF (hdid == h) THEN
        nok=nok+1
     ELSE
        nbad=nbad+1
     END IF


     IF (.NOT.(flag_do_reconstruction)  .AND. (x-x2)*(x2-x1) >= 0.0) THEN
        PRINT*,'MODPK: This could be a model for which inflation does not end.'
        PRINT*,'MODPK: Either adjust phi_init or use slowroll_infl_end for a potential'
        PRINT*,'MODPK: for which inflation does not end by breakdown of slowroll.'
        PRINT*,'MODPK: QUITTING'
        write(*,*) 'vparams: ', (vparams(i),i=1,max_vparams)
        if (.not.instreheat) write(*,*) 'N_pivot: ', N_pivot
        STOP       
        RETURN
     END IF
     
     IF(getEps(y(1),y(2)) .LT. 1 .AND. .NOT.(slowroll_start)) slowroll_start=.true.

     IF(ode_infl_end .AND. .not. flag_do_reconstruction) THEN 
           IF (slowroll_infl_end) THEN
              IF(getEps(y(1),y(2)) .GT. 2 .AND. slowroll_start) THEN 
                 infl_ended=.TRUE.
                 ystart(:)=y(:)
                 IF (save_steps) CALL save_a_step
                 RETURN
              ENDIF
           ELSE
              IF(getEps(y(1),y(2)) .GT. 1 .AND. slowroll_start) THEN
                 print*,getEps(y(1),y(2))
                 PRINT*,'MODPK: You asked for a no-slowroll-breakdown model, but inflation'
                 PRINT*,'MODPK: already ended via slowroll violation before your phi_end was'
                 PRINT*,'MODPK: reached. Please take another look at your inputs.'
                 PRINT*,'MODPK: QUITTING'
                 STOP
              ENDIF
              
              IF (phidot_sign.GT.0..AND.(y(1).GT.(phi_infl_end+0.05))) THEN
                 infl_ended=.TRUE.
                 ystart(:)=y(:)
                 IF (save_steps) CALL save_a_step
                 RETURN
              ENDIF
              IF (phidot_sign.LT.0..AND.(y(1).LT.(phi_infl_end-0.05))) THEN
                 infl_ended=.TRUE.
                 ystart(:)=y(:)
                 IF (save_steps) CALL save_a_step
                 RETURN
              ENDIF
              
           ENDIF
     ENDIF

     IF(ode_ps_output) THEN 
        if (k.ge.a_init*EXP(x)*getH(y(1),y(2))) then
           phi_ik=y(1)
        endif
        
        IF(k.LT.a_init*EXP(x)*getH(y(1),y(2))/eval_ps) THEN !k<aH/eval_ps
           z=a_init*EXP(x)*y(2) ! z=a*phi_dot/H=a*dphi/dalpha
           scalefac=a_init*EXP(x)
           pow_ik=powerspectrum(y(3),y(5),z) !Calculate scalar power at this k
           powt_ik=tensorpower(y(7),y(9),scalefac) !Calculate tensor power at this k
           ystart(:)=y(:)
           IF (save_steps) CALL save_a_step
           RETURN
        END IF
     END IF

     IF (ode_underflow) RETURN
     IF (ABS(hnext) < hmin) THEN
        write(*,*) 'stepsize smaller than minimum in odeint'
        STOP
     END IF
     h=hnext
  END DO

  if (kount.gt.100000) then
     if (abs((yp(1,kount)-yp(1,kount-100000))/yp(1,kount)) .lt. 0.001) then
!        PRINT*,'inflating forever!'
        return
     endif
  endif
     
  
  PRINT*,'too many steps in odeint'
  PRINT*,'vparams = ', vparams
  ode_underflow=.TRUE.
CONTAINS
  SUBROUTINE save_a_step
    kount=kount+1
    IF (kount > SIZE(xp)) THEN
       xp=>reallocate_rv(xp,2*SIZE(xp))
       yp=>reallocate_rm(yp,SIZE(yp,1),SIZE(xp))
    END IF
    xp(kount)=x
    yp(:,kount)=y(:)
    xsav=x
  END SUBROUTINE save_a_step
!  (C) Copr. 1986-92 Numerical Recipes Software, adapted.
END SUBROUTINE odeint

end module modpk_odeint

