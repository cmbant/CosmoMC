MODULE modpk_utils
  IMPLICIT NONE

CONTAINS


  SUBROUTINE bderivs(x,y,yprime)
    USE modpkparams
    USE potential, ONLY: pot,dVdphi,d2Vdphi2,getH,getHdot, getHJeps, getHJH,getHJHdot
    USE camb_interface, ONLY : pk_bad
    REAL*8, INTENT(IN) :: x
    REAL*8, DIMENSION(:), INTENT(IN) :: y
    REAL*8, DIMENSION(:), INTENT(OUT) :: yprime
    REAL*8 :: p,delp,ddp,hubble,dhubble
    REAL*8 :: a,u1,du1,u2,du2
    integer :: i
    !
    !     x=alpha
    !     y(1)=phi              dydx(1)=dphi/dalpha
    !     y(2)=dphi/dalpha      dydx(2)=d^2phi/dalpha^2
    !  
    p=y(1)
    delp=y(2)
    IF(ABS(delp).GT.SQRT(6.0) .AND. .NOT. flag_do_reconstruction ) THEN
       WRITE(*,*) 'MODPK: H is imaginary:',x,p,delp
       !In the case of the hilltop potential, the integrator
       !in a trial step can go here very occasionally because
       !the trial step is too large and it has come too close to V=0.
       !We will stop it going this way, and the code will find the
       !correct epsilon=1 point which has, by definition, to be
       !before this problematic region is reached.
       IF(potential_choice.eq.6) THEN
          yprime(1)=0.0d0
          yprime(2)=0.0d0
          RETURN
       ENDIF
       WRITE(*,*) 'MODPK: QUITTING'
       write(*,*) 'vparams: ', (vparams(i),i=1,max_vparams)
       if (.not.instreheat) write(*,*) 'N_pivot: ', N_pivot
       write (*,*) 'eps(phi):' , getHJeps(p)
       STOP
    END IF
    IF(flag_do_reconstruction) THEN
       hubble=getHJH(p)
       dhubble=getHJHdot(p,delp)
    ELSE
       hubble=getH(p,delp)
       dhubble=getHdot(p,delp)
    ENDIF
    yprime(1)=delp
    yprime(2)=-((3.0+dhubble/hubble)*delp+dVdphi(p)/hubble/hubble)
    RETURN
  END SUBROUTINE bderivs 


  SUBROUTINE derivs(x,y,yprime)
    USE camb_interface, ONLY : pk_bad
    USE modpkparams
    USE internals
    USE potential, ONLY: pot,dVdphi,d2Vdphi2,getH,getHdot, getHJeps, getHJH, getHJHdot
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: x
    REAL*8, DIMENSION(:), INTENT(IN) :: y
    REAL*8, DIMENSION(:), INTENT(OUT) :: yprime
    REAL*8 :: p,delp,ddp,hubble,dhubble
    REAL*8 :: a,u1,du1,u2,du2,v1,dv1,v2,dv2
    integer :: i
    !
    !     x=alpha
    !     y(1)=phi              dydx(1)=dphi/dalpha
    !     y(2)=dphi/dalpha      dydx(2)=d^2phi/dalpha^2
    !     y(3)=u1               dydx(3)=du1/dalpha
    !     y(4)=du1/dalpha       dydx(4)=d^2u1/dalpha^2
    !     y(5)=u2               dydx(5)=du2/dalpha
    !     y(6)=du2/dalpha       dydx(6)=d^2u2/dalpha^2
    !     y(7)=v1               dydx(7)=dv1/dalpha
    !     y(8)=dv1/dalpha       dydx(8)=d^2v1/dalpha^2
    !     y(9)=v2               dydx(9)=dv2/dalpha
    !     y(10)=dv2/dalpha      dydx(10)=d^2v2/dalpha^2
    !
    p=y(1)
    delp=y(2)


    IF(flag_do_reconstruction) THEN
!       print *, 'phi, eps: ', p, getHJeps(p)
       IF ( getHJeps(p) .GT. 1 .AND. pk_bad.ne.99) THEN
          WRITE(*,*) 'MODPK: WARNING'
          WRITE(*,*)  'MODPK: Inflation ended during the mode evolution!'
          WRITE(*,*) 'MODPK: alpha, phi, dphi/dalpha',x,p,delp
          write (*,*) 'MODPK: eps(phi) = ' , getHJeps(p)
!          WRITE(*,*) 'MODPK: QUITTING'
          write(*,*)  'vparams: ', (vparams(i),i=1,max_vparams)
!          STOP
          pk_bad=99
       ENDIF

       hubble=getHJH(p)
       dhubble=getHJHdot(p,delp)
    ELSE
       IF(ABS(delp).GT.SQRT(6.)) THEN
          WRITE(*,*) 'MODPK: H is imaginary:',x,p,delp
          WRITE(*,*) 'MODPK: QUITTING'
          write(*,*) 'vparams: ', (vparams(i),i=1,max_vparams)
          if (.not.instreheat) write(*,*) 'N_pivot: ', N_pivot
          STOP
       END IF

       hubble=getH(p,delp)
       dhubble=getHdot(p,delp)
    ENDIF
    a=EXP(x)
    u1=y(3)
    du1=y(4)
    u2=y(5)
    du2=y(6)
    v1=y(7)
    dv1=y(8)
    v2=y(9)
    dv2=y(10)
    yprime(1)=delp
    yprime(2)=-((3.0+dhubble/hubble)*delp+dVdphi(p)/hubble/hubble)
    ddp=yprime(2)
    yprime(3)=du1
    yprime(4)=-((dhubble/hubble+1.)*du1 + ((k/a_init/a/hubble)**2. &
         &     +(-2.0 + 5.*dhubble/hubble + d2Vdphi2(p)/hubble/hubble &
         &     + 2.*(dhubble/hubble)**2. + 4*(dhubble/hubble)*ddp/delp))*u1)
    yprime(5)=du2
    yprime(6)=-((dhubble/hubble+1.)*du2 + ((k/a_init/a/hubble)**2. &
         &     +(-2.0 + 5.*dhubble/hubble + d2Vdphi2(p)/hubble/hubble &
         &     + 2.*(dhubble/hubble)**2. + 4*(dhubble/hubble)*ddp/delp))*u2)
    yprime(7)=dv1
    yprime(8)=-((dhubble/hubble+1.)*dv1 + ((k/a_init/a/hubble)**2. &
         &     -(dhubble/hubble+2.))*v1)
    yprime(9)=dv2
    yprime(10)=-((dhubble/hubble+1.)*dv2 + ((k/a_init/a/hubble)**2. &
         &     -(dhubble/hubble+2.))*v2)

!    print *, 'x = ', x
!    print *, 'y = ', y(1),y(2)
!    print *, 'yp = ', yprime(1),yprime(2)

    RETURN
  END SUBROUTINE derivs


  SUBROUTINE rkqs(y,dydx,x,htry,eps,yscal,hdid,hnext,derivs)
        USE ode_path
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
	INTEGER*4 :: ndum
	REAL*8 :: errmax,h,htemp,xnew
	REAL*8, DIMENSION(size(y)) :: yerr,ytemp
	REAL*8, PARAMETER :: SAFETY=0.9d0,PGROW=-0.2d0,PSHRNK=-0.25d0,&
		ERRCON=1.89e-4
        if (size(y)==size(dydx) .and. size(dydx)==size(yscal)) then
           ndum = size(y)
        else
           write(*,*) 'Wrong array sizes in rkqs'
           stop
        end if
	h=htry
	do
		call rkck(y,dydx,x,h,ytemp,yerr,derivs)
		errmax=maxval(abs(yerr(:)/yscal(:)))/eps
		if (errmax <= 1.0) exit
		htemp=SAFETY*h*(errmax**PSHRNK)
		h=sign(max(abs(htemp),0.1d0*abs(h)),h)
		xnew=x+h
		if (xnew == x) then 
                   print*, 'stepsize underflow in rkqs'
                   ode_underflow = .true.
                   return
                endif
	end do
	if (errmax > ERRCON) then
		hnext=SAFETY*h*(errmax**PGROW)
	else
		hnext=5.0d0*h
	end if
	hdid=h
	x=x+h
	y(:)=ytemp(:)
!  (C) Copr. 1986-92 Numerical Recipes Software, adapted.
	END SUBROUTINE rkqs

	SUBROUTINE rkck(y,dydx,x,h,yout,yerr,derivs)
	IMPLICIT NONE
	REAL*8, DIMENSION(:), INTENT(IN) :: y,dydx
	REAL*8, INTENT(IN) :: x,h
	REAL*8, DIMENSION(:), INTENT(OUT) :: yout,yerr
	INTERFACE
		SUBROUTINE derivs(x,y,dydx)
		IMPLICIT NONE
		REAL*8, INTENT(IN) :: x
		REAL*8, DIMENSION(:), INTENT(IN) :: y
		REAL*8, DIMENSION(:), INTENT(OUT) :: dydx
		END SUBROUTINE derivs
	END INTERFACE
	INTEGER*4 :: ndum
	REAL*8, DIMENSION(size(y)) :: ak2,ak3,ak4,ak5,ak6,ytemp
	REAL*8, PARAMETER :: A2=0.2d0,A3=0.3d0,A4=0.6d0,A5=1.0d0,&
		A6=0.875d0,B21=0.2d0,B31=3.0d0/40.0d0,B32=9.0d0/40.0d0,&
		B41=0.3d0,B42=-0.9d0,B43=1.2d0,B51=-11.0d0/54.0d0,&
		B52=2.5d0,B53=-70.0d0/27.0d0,B54=35.0d0/27.0d0,&
		B61=1631.0d0/55296.0d0,B62=175.0d0/512.0d0,&
		B63=575.0d0/13824.0d0,B64=44275.0d0/110592.0d0,&
		B65=253.0d0/4096.0d0,C1=37.0d0/378.0d0,&
		C3=250.0d0/621.0d0,C4=125.0d0/594.0d0,&
		C6=512.0d0/1771.0d0,DC1=C1-2825.0d0/27648.0d0,&
		DC3=C3-18575.0d0/48384.0d0,DC4=C4-13525.0d0/55296.0d0,&
		DC5=-277.0d0/14336.0d0,DC6=C6-0.25d0
        if (size(y)==size(dydx) .and. size(dydx)==size(yout) .and. size(yout)==size(yerr)) then
           ndum = size(y)
        else
           write(*,*) 'Wrong array sizes in rkck'
           stop
        end if
	ytemp=y+B21*h*dydx
	call derivs(x+A2*h,ytemp,ak2)
	ytemp=y+h*(B31*dydx+B32*ak2)
	call derivs(x+A3*h,ytemp,ak3)
	ytemp=y+h*(B41*dydx+B42*ak2+B43*ak3)
	call derivs(x+A4*h,ytemp,ak4)
	ytemp=y+h*(B51*dydx+B52*ak2+B53*ak3+B54*ak4)
	call derivs(x+A5*h,ytemp,ak5)
	ytemp=y+h*(B61*dydx+B62*ak2+B63*ak3+B64*ak4+B65*ak5)
	call derivs(x+A6*h,ytemp,ak6)
	yout=y+h*(C1*dydx+C3*ak3+C4*ak4+C6*ak6)
	yerr=h*(DC1*dydx+DC3*ak3+DC4*ak4+DC5*ak5+DC6*ak6)
!  (C) Copr. 1986-92 Numerical Recipes Software, adapted.
      END SUBROUTINE rkck

      FUNCTION locate(xx,x)
	IMPLICIT NONE
	REAL*8, DIMENSION(:), INTENT(IN) :: xx
	REAL*8, INTENT(IN) :: x
	INTEGER*4 :: locate
	INTEGER*4 :: n,jl,jm,ju
	LOGICAL :: ascnd
	n=size(xx)
	ascnd = (xx(n) >= xx(1))
	jl=0
	ju=n+1
	do
		if (ju-jl <= 1) exit
		jm=(ju+jl)/2
		if (ascnd .eqv. (x >= xx(jm))) then
			jl=jm
		else
			ju=jm
		end if
	end do
	if (x == xx(1)) then
		locate=1
	else if (x == xx(n)) then
		locate=n-1
	else
		locate=jl
	end if
!  (C) Copr. 1986-92 Numerical Recipes Software, adapted.
      END FUNCTION locate

      SUBROUTINE polint(xa,ya,x,y,dy)
	IMPLICIT NONE
	REAL*8, DIMENSION(:), INTENT(IN) :: xa,ya
	REAL*8, INTENT(IN) :: x
	REAL*8, INTENT(OUT) :: y,dy
	INTEGER*4 :: m,n,ns
        INTEGER*4, DIMENSION(1) :: imin
	REAL*8, DIMENSION(size(xa)) :: c,d,den,ho,absho
        if (size(xa)==size(ya)) then
           n=size(xa)
        else
           write(*,*) 'Wrong array sizes in polint'
           stop
        end if
	c=ya
	d=ya
	ho=xa-x
        absho=abs(ho)
        imin=minloc(absho(:))
        ns=imin(1)
	y=ya(ns)
	ns=ns-1
	do m=1,n-1
		den(1:n-m)=ho(1:n-m)-ho(1+m:n)
		if (any(den(1:n-m) == 0.0)) then
                     write(*,*) 'polint: calculation failure'
                     stop
                  end if
		den(1:n-m)=(c(2:n-m+1)-d(1:n-m))/den(1:n-m)
		d(1:n-m)=ho(1+m:n)*den(1:n-m)
		c(1:n-m)=ho(1:n-m)*den(1:n-m)
		if (2*ns < n-m) then
			dy=c(ns+1)
		else
			dy=d(ns)
			ns=ns-1
		end if
		y=y+dy
	end do
!  (C) Copr. 1986-92 Numerical Recipes Software, adapted.
      END SUBROUTINE polint

      FUNCTION reallocate_rv(p,n)
        REAL*8, DIMENSION(:), POINTER :: p, reallocate_rv
        INTEGER*4, INTENT(IN) :: n
        INTEGER*4 :: nold,ierr
        allocate(reallocate_rv(n),stat=ierr)
        if (ierr /= 0) then
           write(*,*) 'reallocate_rv: problem in attempt to allocate memory'
           stop
        end if
        if (.not. associated(p)) RETURN
        nold=size(p)
        reallocate_rv(1:min(nold,n))=p(1:min(nold,n))
        deallocate(p)
!  (C) Copr. 1986-92 Numerical Recipes Software, adapted.
      END FUNCTION reallocate_rv

      FUNCTION reallocate_rm(p,n,m)
        REAL*8, DIMENSION(:,:), POINTER :: p, reallocate_rm
        INTEGER*4, INTENT(IN) :: n,m
        INTEGER*4 :: nold,mold,ierr
        allocate(reallocate_rm(n,m),stat=ierr)
        if (ierr /= 0) then
           write(*,*) 'reallocate_rm: problem in attempt to allocate memory'
           stop
        end if
        if (.not. associated(p)) RETURN
        nold=size(p,1)
        mold=size(p,2)
        reallocate_rm(1:min(nold,n),1:min(mold,m))=&
             p(1:min(nold,n),1:min(mold,m))
        deallocate(p)
!  (C) Copr. 1986-92 Numerical Recipes Software, adapted.
      END FUNCTION reallocate_rm

END MODULE modpk_utils
