! Covariance matrix, diagonalization

module utils1
  !$ use omp_lib
  use RandomNS
  implicit none
  
  	integer lwork,liwork
  	data lwork,liwork /1,1/
      logical setBlk
      data setBlk /.false./

contains

!---------------------------------------------------------------------- 
  
 subroutine Diagonalize(a,diag,n,switch)
	!Does m = U diag U^T, returning U in m
	integer n,id
	double precision a(n,n),diag(n)
	logical switch
	integer i,j,ierr,il,iu,m,isuppz(2*n)
	double precision vl,vu,abstol,z(n,n)
	double precision, dimension(:), allocatable :: work
	integer, dimension(:), allocatable :: iwork

	diag=0.
      abstol=0.
      
      id=0
      !$ id=omp_get_thread_num()
      
      if(.not.setBlk .or. switch) then
      	if(id==0) then
      		!find optimal block sizes
      		allocate(work(1))
            	allocate(iwork(1))
            	lwork=-1
            	liwork=-1
		
            	call DSYEVR( 'V', 'A', 'U', n, a, n, vl, vu, il, iu, &
      		abstol, m, diag, Z, n, isuppz, work, lwork, iwork, liwork, ierr )

            	lwork=int(work(1))
            	liwork=iwork(1)
            	deallocate(work,iwork)
            	setBlk=.true.
	else
            	do
                  	!$OMP FLUSH(setBlk)
                  	if(setBlk) exit
		enddo
	end if
      endif 
      	
      allocate(work(lwork),iwork(liwork))
      call DSYEVR( 'V', 'A', 'U', n, a, n, vl, vu, il, iu, &
      abstol, m, diag, Z, n, isuppz, work, lwork, iwork, liwork, ierr )
      a=z
      deallocate(work,iwork)
    
      !check for inf & nan
      do i=1,n
    	if(diag(i)/=diag(i) .or. diag(i)>huge(1d0)) then
		diag=1d0
		a=0d0
		do j=1,n
			a(j,j)=1d0
		enddo
		exit
	endif
      enddo
	  
 end subroutine Diagonalize

!---------------------------------------------------------------------- 

! LogSumExp(x,y)=log(exp(x)+exp(y))
  double precision function LogSumExp(x,y)
    double precision x,y

    if (x.gt.y) then
       LogSumExp=x+log(1+exp(y-x))
    else
       LogSumExp=y+log(1+exp(x-y))
    end if
    return

  end function LogSumExp
  
!----------------------------------------------------------------------

  subroutine calc_covmat(npt,np,p,mean,covmat)
    implicit none
    integer np,npt
    double precision p(:,:)
    double precision covmat(:,:)
    double precision mean(:)
    integer i,j,ip

    covmat=0.
    do i=1,np
       do j=i,np
          do ip=1,npt
             covmat(i,j)=covmat(i,j)+((p(i,ip)-mean(i))*(p(j,ip)-mean(j)))
          end do
          if(j/=i) covmat(j,i)=covmat(i,j)
       end do
       covmat(i,1:np)=covmat(i,1:np)/dble(npt-1)
    end do
   
    return
    
  end subroutine calc_covmat

!----------------------------------------------------------------------

  subroutine calc_covmat_wt(npt,np,p,wt,mean,covmat)
    implicit none
    integer np,npt
    double precision p(:,:),wt(:)
    double precision covmat(:,:)
    double precision mean(:)
    integer i,j,ip

    covmat=0.
    do i=1,np
       do j=i,np
          do ip=1,npt
             covmat(i,j)=covmat(i,j)+((p(i,ip)-mean(i))*(p(j,ip)-mean(j)))*wt(ip)
          end do
          if(j/=i) covmat(j,i)=covmat(i,j)
       end do
    end do
   
    return
    
  end subroutine calc_covmat_wt

!----------------------------------------------------------------------

  !inv_cov=(evec).(inv_eval).Transpose(evec)
  subroutine calc_invcovmat(numdim,evec,eval,invcov)
    implicit none
    integer numdim !num of dimensions & points
    double precision evec(:,:),eval(:) !eigenvectors & eigenvalues
    double precision invcov(:,:)
    integer i,j,k

    invcov=0d0
    do i=1,numdim
       do j=i,numdim
          do k=1,numdim
             invcov(i,j)=invcov(i,j)+(evec(i,k)*evec(j,k)/eval(k))
          end do
          if(i/=j) invcov(j,i)=invcov(i,j)
       end do
    end do
    
    return
    
  end subroutine calc_invcovmat

!----------------------------------------------------------------------

  subroutine ScaleFactor(npt,ndim,pt,mean,inv_cov,kmax)
    implicit none
 
    double precision kmax
    integer npt,ndim
    double precision pt(ndim,npt),ptM(1,ndim),sf(1,1)
    integer i,k
    double precision inv_cov(ndim,ndim)
    double precision mean(ndim)
    double precision temp_p(ndim,npt)
    
    do i=1,ndim
	temp_p(i,1:npt)=pt(i,1:npt)-mean(i)
    enddo
    
    kmax=0d0
    do k=1,npt
    	ptM(1,:)=temp_p(:,k)
    	sf=MatMul(MatMul(ptM,inv_cov),Transpose(ptM))
        if(sf(1,1)>kmax) kmax=sf(1,1)
    enddo
    
    !kmax=sqrt(kmax)
    
  end subroutine ScaleFactor
   
!----------------------------------------------------------------------	
 
 double precision function ptScaleFac(ndim,pt,meanx,invcovx)
 	
	implicit none
      
	integer ndim!dimensionality
      	double precision pt(ndim),point(1,ndim)!point to be checked
      	double precision meanx(ndim),invcovx(ndim,ndim)!cluster attributes
      	double precision kfac(1,1)!k factor of present point
      
      	point(1,:)=pt(:)-meanx(:)
	kfac=MATMUL(MATMUL(point,invcovx),Transpose(point))
	ptScaleFac=kfac(1,1)
      
 end function ptScaleFac
  
!----------------------------------------------------------------------
   
   subroutine getTMatrix(ndim,evec,eval,TMatrix)
    
    implicit none
 
    double precision TMatrix(ndim,ndim),evec(ndim,ndim),eval(ndim)
    integer ndim,np,i,j
    double precision TMat(ndim,ndim)
    
    np=ndim
    TMatrix=evec
    
    do i=1,np
       do j=1,np
          TMatrix(i,j) = evec(i,j) * sqrt(eval(j))
       enddo
    enddo
    TMat=Transpose(TMatrix)
    TMatrix=TMat
    
  end subroutine getTMatrix
  
!----------------------------------------------------------------------

  !return the Mahalanobis distance of a point from the given ellipsoid
  double precision function MahaDis(ndim,pt,mean,inv_cov,kfac)
    implicit none
    
    !input varaibles
    integer ndim !dimensionality
    double precision pt(ndim) !point
    double precision mean(ndim) !centroid of the ellipsoid
    double precision inv_cov(ndim,ndim) !inv covariance matrix of the ellipsoid
    double precision kfac !enlargement factor
    
    MahaDis=ptScaleFac(ndim,pt,mean,inv_cov)/kfac
    
  end function MahaDis
   
!----------------------------------------------------------------------	
  
  double precision function ellVol(ndim,eval,k_fac)
    
    implicit none
    
    integer ndim
    double precision eval(ndim),k_fac,pi
    integer i
    
    pi=4.d0*atan(1.d0)
    
    ellVol=sqrt(product(eval)*(k_fac**dble(ndim)))
    
    if(mod(ndim,2)==0) then
    	do i=2,ndim,2
		ellVol=ellVol*2.d0*pi/dble(i)
	enddo
    else
    	ellVol=ellVol*2.
    	do i=3,ndim,2
		ellVol=ellVol*2.d0*pi/dble(i)
	enddo
    endif
  end function ellVol

!---------------------------------------------------------------------- 
  
  subroutine genPtInSpheroid(np,u,id)
   
    implicit none
    
    integer np,i,id
    double precision u(:), mod, urv
  
    mod=0d0
    do i=1,np
       u(i)=gaussian1NS(id)
       mod=mod+u(i)**2.
    enddo
    urv=ranmarNS(id)
    mod=(urv**dble(1./np))/sqrt(mod)
    u(:)=mod*u(:)
    
  end subroutine genPtInSpheroid

!----------------------------------------------------------------------   
  
  subroutine genPtOnSpheroid(np,u,id)
   
    implicit none
    
    integer np,i,id
    double precision u(:), mod
  
    mod=0.d0
    do i=1,np
       u(i)=gaussian1NS(id)
       mod=mod+u(i)**2.
    enddo
    mod=1.d0/sqrt(mod)
    u(:)=mod*u(:)
    
  end subroutine genPtOnSpheroid

!----------------------------------------------------------------------  

  !calculate the mean, covariance matrix, inv covariance matrix, evalues,
  !evects, det(cov) & enlargement of a given point set
  subroutine CalcEllProp(npt,ndim,pt,mean,covmat,invcov,tMat,evec,eval,detcov, &
  	kfac,eff,vol,pVol,switch)
    implicit none
    !input variables
    integer npt !no. of points
    integer ndim !dimensionality
    double precision pt(ndim,npt) !points
    double precision pVol !min vol this ellipsoid should occupy
    logical switch !initialize the LAPACK eigen analysis routines?
    !output variables
    double precision mean(ndim) !mean
    double precision covmat(ndim,ndim) !covariance matrix
    double precision invcov(ndim,ndim) !inverse covariance matrix
    double precision tMat(ndim,ndim) !transformation matrix
    double precision evec(ndim,ndim) !eigenvectors
    double precision eval(ndim)
    double precision detcov !determinant of the covariance matrix
    double precision kfac !enlargement point factor
    double precision eff !enlargement volume factor
    double precision vol !ellipsoid volume
    !work variables
    integer i,j
    
    !calculate the mean
    do i=1,ndim
    	mean(i)=sum(pt(i,1:npt))/npt
    enddo
    
    if(npt<=1) then
    	kfac=0.d0
	vol=0.d0
	eval=0.d0
	eff=1.d0
	return
    endif
    
    !calculate the covariance matrix
    call calc_covmat(npt,ndim,pt,mean,covmat)
    
    !eigen analysis
    evec=covmat
    call Diagonalize(evec,eval,ndim,switch)
    
    !eigenvalues of covariance matrix can't be zero
    do i=1,ndim-1
	if(eval(i)<=0.d0) eval(1:i)=eval(i+1)/2.
    enddo
	
    !if the no. of points is less than ndim+1 then set the eigenvalues of unconstrained 
    !dimensions equal to the min constrained eigenvalue
    if(npt<ndim+1) then
	eval(1:ndim+1-npt)=eval(ndim+2-npt)
    endif
    
    !calculate the determinant of covariance matrix
    detcov=product(eval)
    
    !calculate the inverse of covariance matrix
    call calc_invcovmat(ndim,evec,eval,invcov)
    
    !calculate the enlargement factor
    call ScaleFactor(npt,ndim,pt,mean,invcov,kfac)
    
    !calculate the ellipsoid volume
    vol=ellVol(ndim,eval,kfac)
    
    !scale the enlargement factor if vol<pVol
    if(pVol>0d0 .and. vol<pVol) then
    	eff=(pVol/vol)**(2.d0/ndim)
	vol=pVol
    else
    	eff=1.d0
    endif
        
    !calculate the transformation matrix
    call getTMatrix(ndim,evec,eval,tMat)
    
  end subroutine CalcEllProp

!---------------------------------------------------------------------- 

  !calculate the inv covariance matrix, evalues & evects of a given point set
  subroutine CalcBEllInfo(npt,ndim,pt,mean,eval,evec,invcov,tMat,kfac,minPt)
    implicit none
    !input variables
    integer npt !no. of points
    integer ndim !dimensionality
    double precision pt(ndim,npt) !points
    integer minPt
    !output variables
    double precision mean(ndim) !centroid
    double precision invcov(ndim,ndim) !inverse covariance matrix
    double precision tMat(ndim,ndim) !transformation matrix
    double precision evec(ndim,ndim) !eigenvectors
    double precision eval(ndim) !eigenvalues
    double precision kfac !multiplicative factor for covmat for bounding ell
    !work variables
    double precision sigma(ndim)
    integer i,j,k,n
    parameter(n=0)
    double precision covmat(ndim,ndim),ptk(ndim,npt),dist(npt)
    double precision maxIndx(n,2)
    integer nOut,ptOut(n),nptk
    logical flag
    
    
    !calculate the mean
    do i=1,ndim
    	mean(i)=sum(pt(i,1:npt))/npt
    	sigma(i)=sum(pt(i,1:npt)**2)/npt
    enddo
    sigma(1:ndim)=sqrt(sigma(1:ndim)-mean(1:ndim)**2)
    
    if(npt==1) then
	eval=0.d0
	kfac=0.d0
	return
    endif
    
    !remove outliers
    if(npt>minPt) then
	k=min(n,npt-minPt)
	if(k>0) then
		maxIndx(:,2)=0d0
    		do i=1,npt
			dist(i)=0d0
			do j=1,ndim
				dist(i)=dist(i)+abs((pt(j,i)-mean(j)))/sigma(j)
			enddo
		
			do j=1,k+1
				if(j==k+1) exit
				if(dist(i)<maxIndx(j,2)) exit
			enddo
			j=j-1
			if(j>0) then
				maxIndx(1:j-1,:) = maxIndx(2:j,:)
				maxIndx(j,1) = dble(i)
				maxIndx(j,2) = dist(i)
			endif
    		enddo
	
		nOut=0
		do i=1,k
			if(maxIndx(i,2)>3d0) then
				nOut=nOut+1
				ptOut(nOut)=int(maxIndx(i,1))
			endif
		enddo
	else
		nOut=0
	endif
	
	if(nOut>0) then
		j=0
		do i=1,npt
			flag=.false.
			do k=1,nOut
				if(i == ptOut(k)) then
					flag=.true.
					exit
				endif
			enddo
			
			if(flag) cycle
			
			ptk(:,j+1)=pt(:,i)
			j=j+1
		enddo
		nptk=npt-nOut
	else
		ptk=pt
		nptk=npt
	endif
    else
	ptk=pt
	nptk=npt
    endif
    
    
    !calculate the covariance matrix
    call calc_covmat(nptk,ndim,ptk,mean,covmat)
    
    !eigen analysis
    evec=covmat
    call Diagonalize(evec,eval,ndim,.false.)
    !eigenvalues of covariance matrix can't be zero
    do i=1,ndim-1
	if(eval(i)<=0.d0) eval(1:i)=eval(i+1)/2.
    enddo
	
    !if the no. of points is less than ndim+1 then set the eigenvalues of unconstrained 
    !dimensions equal to the min constrained eigenvalue
    if(nptk<ndim+1) then
	eval(1:ndim+1-nptk)=eval(ndim+2-nptk)
    endif
    
    !calculate the inverse of covariance matrix
    call calc_invcovmat(ndim,evec,eval,invcov)
    
    !now scale the eigenvalues so that its a bounding ellipsoid with k = kfac
    call ScaleFactor(nptk,ndim,ptk,mean,invcov,kfac)
        
    !calculate the transformation matrix
    call getTMatrix(ndim,evec,eval,tMat)
    
  end subroutine CalcBEllInfo

!----------------------------------------------------------------------

  !generate a point uniformly inside the given ellipsoid
  subroutine genPtInEll(ndim,mean,efac,TMat,id,pt)
  	implicit none
	
	!input variables
	integer ndim !dimensionality
	double precision mean(ndim) !centroid of the given ellipsoid
	double precision efac !enlargement factor of the given ellipsoid
	double precision TMat(ndim,ndim) !transformation matrix of the given ellipsoid
	integer id !processor id (for OpenMP)
	!output variable
	double precision pt(ndim)
	!work variables
	double precision u(1,ndim),pnewM(1,ndim)
	
	call genPtInSpheroid(ndim,u(1,:),id)
	pnewM=MatMul(u,TMat)
	pt(:)=sqrt(efac)*pnewM(1,:)+mean(:)
	
  end subroutine genPtInEll

!----------------------------------------------------------------------

  !generate a point uniformly on the surface of the given ellipsoid
  subroutine genPtOnEll(ndim,mean,efac,TMat,id,pt)
  	implicit none
	
	!input variables
	integer ndim !dimensionality
	double precision mean(ndim) !centroid of the given ellipsoid
	double precision efac !enlargement factor of the given ellipsoid
	double precision TMat(ndim,ndim) !transformation matrix of the given ellipsoid
	integer id !processor id (for OpenMP)
	!output variable
	double precision pt(ndim)
	!work variables
	double precision u(1,ndim),pnewM(1,ndim)
	
	call genPtOnSpheroid(ndim,u(1,:),id)
	pnewM=MatMul(u,TMat)
	pt(:)=sqrt(efac)*pnewM(1,:)+mean(:)
	
  end subroutine genPtOnEll

!----------------------------------------------------------------------
   
   !evolve an ellipsoid from which a point has either been rejected or a point
   !has been inserted
  subroutine evolveEll(a_r,npt,ndim,newpt,pts,mean,eval,invcov,kfac,eff,vol,pVol)
    
	implicit none
    	!input variables
	integer a_r !if 0 then point rejected, 1 then point inserted
	integer npt !no. of points after rejection or before insertion
	integer ndim !dimensionality
	double precision newpt(ndim) !point to be rejected/inserted
	double precision pts(ndim,npt) !point set after rejection or before insertion
	double precision mean(ndim) !centroid of the ellipsoid
	double precision eval(ndim) !eigenvalues of the ellipsoid
	double precision invcov(ndim,ndim) !inverse covariance matrix of the ellipsoid
	double precision pVol !target volume
	
	!input/output variables
	double precision kfac !input (output): point enlargement before (after) rejection/insertion
	double precision eff !input (output): volume enlargement before (after) rejection/insertion
	double precision vol !input (output): volume before (after) rejection/insertion
	
	!work variables
	double precision new_pt(ndim,1),new_kfac
	
	
	!sanity check
	if(a_r < 0 .or. a_r > 1) then
		write(*,*)"a_r in evolveEll cannot be", a_r
		stop
	endif
	
	if(pVol==0.d0 .or. eval(ndim)==0.d0) then
		kfac=0.d0
		eff=1.d0
		vol=0.d0
		return
	endif
	
	
	!point inserted
	if(a_r==1) then
		if(npt==0) then
			write(*,*)"can not insert point in an ellipsoid with no points"
			stop
		else
			!calculate the scale factor
			new_pt(:,1)=newpt(:)
			call ScaleFactor(1,ndim,new_pt,mean,invcov,new_kfac)
			if(new_kfac>kfac) then
				eff=eff*kfac/new_kfac
				kfac=new_kfac
				if(eff<1d0) then
					eff=1d0
					vol=ellVol(ndim,eval,kfac)
				endif
			endif
			
			!if the point has been inserted to a cluster with npt > 0 then kfac remains the same
			!scale eff if target vol is bigger
			if(pVol>vol) then
    				eff=eff*(pVol/(vol))**(2.d0/dble(ndim))
				vol=pVol
			else
				!if target vol is smaller then calculate the point volume & 
				!scale eff if required
				vol=ellVol(ndim,eval,kfac)
			
				!scale the enlargement factor if vol<pVol
    				if(vol<pVol) then
    					eff=max(1d0,(pVol/(vol))**(2.d0/dble(ndim)))
					vol=pVol
    				else
    					eff=1.d0
    				endif
			endif
		endif
	else
		if(npt==0) then
			vol=0.d0
			kfac=0.d0
			eff=1.d0
			return
		endif
		
		!if the point has been rejected then figure out if its a boundary point
		new_pt(:,1)=newpt(:)
		call ScaleFactor(1,ndim,new_pt,mean,invcov,new_kfac)
		
		!if it's a boundary point then find new kfac
		if(new_kfac==kfac) then
			call ScaleFactor(npt,ndim,pts,mean,invcov,kfac)
			if(kfac>new_kfac*eff .and. npt>1) then
				write(*,*)"Problem in evoleell"
				write(*,*)kfac,new_kfac,eff,npt
				stop
			endif
			eff=1d0
			vol=ellVol(ndim,eval,kfac)
		endif
		
		!scale the enlargement factor if vol<pVol
    		if(vol<pVol) then
    			eff=(pVol/(vol))**(2.d0/dble(ndim))
			vol=pVol
    		endif
	endif
    
  end subroutine evolveEll
  
!----------------------------------------------------------------------
   
   !enlarge an ellipsoid because of an additional point being inserted in it
  subroutine enlargeEll(ndim,newpt,mean,eval,invcov,kfac,eff,vol,pVol)
    
	implicit none
    	!input variables
	integer ndim !dimensionality
	double precision newpt(ndim) !point to be inserted
	double precision mean(ndim) !centroid of the ellipsoid
	double precision eval(ndim) !eigenvalues of the ellipsoid
	double precision invcov(ndim,ndim) !inverse covariance matrix of the ellipsoid
	double precision pVol !target volume
	
	!input/output variables
	double precision kfac !input (output): point enlargement before (after) insertion
	double precision eff !input (output): volume enlargement before (after) insertion
	double precision vol !input (output): volume before (after) insertion
	
	!work variables
	double precision new_pt(ndim,1),d1
	
	
	!find new kfac
	new_pt(:,1)=newpt(:)
	call ScaleFactor(1,ndim,new_pt,mean,invcov,d1)
	if(d1>kfac) then
		kfac=d1
		eff=1.d0
		vol=ellVol(ndim,eval,kfac)
	endif
		
	if(pVol>vol) then
    		eff=eff*(pVol/(vol))**(2.d0/dble(ndim))
		vol=pVol
	endif
    
  end subroutine enlargeEll
  
!----------------------------------------------------------------------


  function inprior(np,p)
    
    implicit none
    
    logical inprior
    integer np, i
    double precision p(np)
    
    inprior = .true.
    
    do i=1,np
	if(p(i)<0d0 .or. p(i)>1d0) then
          	inprior = .false.
            	return
	endif
    enddo
  
  end function inprior

!----------------------------------------------------------------------
  
  !Returns the value gamma(ln[(xx)]) for xx > 0.
  FUNCTION gammln(xx) 
  	REAL gammln,xx  
      INTEGER j 
      DOUBLE PRECISION ser,stp,tmp,x,y,cof(6) 
      !Internal arithmetic will be done in double precision, 
      !a nicety that you can omit if  ve- gure accuracy is good enough. 
      SAVE cof,stp 
      
      DATA cof,stp/76.18009172947146d0,-86.50532032941677d0,24.01409824083091d0, &
      -1.231739572450155d0,.1208650973866179d-2,-.5395239384953d-5,2.5066282746310005d0/ 
      
      x=xx 
      y=x 
      tmp=x+5.5d0 
      tmp=(x+0.5d0)*log(tmp)-tmp 
      ser=1.000000000190015d0 
      do j=1,6 
      	y=y+1.d0 
            ser=ser+cof(j)/y 
      enddo
      gammln=tmp+log(stp*ser/x) 
      return 
  END FUNCTION gammln
  
  
  
  
  !USES gcf,gser Returns the incomplete gamma function P(a, x). 
  FUNCTION gammp(a,x) 
  	REAL a,gammp,x 
  	REAL gammcf,gamser,gln 
  
  	if(x.lt.0..or.a.le.0.) then
		write(*,*)  'bad arguments in gammp'
		stop
	endif
  	if(x.lt.a+1.)then 
  		!Use the series representation. 
  		call gser(gamser,a,x,gln) 
      	gammp=gamser
  	else 
  		!Use the continued fraction representation 
      	call gcf(gammcf,a,x,gln) 
      	gammp=1.-gammcf !and take its complement. 
  	endif 
  	return 
  END FUNCTION gammp
  
  
  !USES gcf,gser Returns the incomplete gamma function Q(a,x)=1-P(a,x). 
  FUNCTION gammq(a,x) 
  	REAL a,gammq,x 
      REAL gammcf,gamser,gln 
      
      if(x.lt.0..or.a.le.0.) then
      	write(*,*) 'bad arguments in gammq'
	stop
      endif
      if(x.lt.a+1.)then 
      	!Use the series representation 
      	call gser(gamser,a,x,gln) 
            gammq=1.-gamser !and take its complement. 
      else 
      	!Use the continued fraction representation. 
            call gcf(gammcf,a,x,gln) 
            gammq=gammcf 
	endif 
      return 
  END FUNCTION gammq
  
  
  !USES gammln Returns the incomplete gamma function P(a,x) evaluated by its series 
  !representation as gamser. 
  !Also returns ln (a) as gln. 
  SUBROUTINE gser(gamser,a,x,gln) 
  	INTEGER ITMAX 
      REAL a,gamser,gln,x,EPS 
      PARAMETER (ITMAX=100,EPS=3.e-7) 
      INTEGER n 
      REAL ap,del,sum 
      
      gln=gammln(a) 
      if(x.le.0.)then 
      	if(x.lt.0.) write(*,*)  'x < 0 in gser'  
            gamser=0. 
            return 
      endif 
      ap=a 
      sum=1./a 
      del=sum 
      do n=1,ITMAX 
      	ap=ap+1. 
            del=del*x/ap 
            sum=sum+del 
            if(abs(del).lt.abs(sum)*EPS)goto 91 
      enddo
      write(*,*)  'a too large, ITMAX too small in gser'  
    91 gamser=sum*exp(-x+a*log(x)-gln) 
    	return 
  END SUBROUTINE gser
   
   
  !USES gammln Returns the incomplete gamma function Q(a, x) evaluated by its continued 
  !fraction representation as gammcf. Also returns ln  (a) as gln. 
  !Parameters: ITMAX is the maximum allowed number of iterations; 
  !EPS is the relative accuracy; FPMIN is a number near the smallest representable 
  !floating-point number.  
  SUBROUTINE gcf(gammcf,a,x,gln) 
  	INTEGER ITMAX 
  	REAL a,gammcf,gln,x,EPS,FPMIN 
      PARAMETER (ITMAX=100,EPS=3.e-7,FPMIN=1.e-30) 
      INTEGER i 
      REAL an,b,c,d,del,h 
      
      gln=gammln(a) 
      b=x+1.-a 
      c=1./FPMIN 
      d=1./b 
      h=d 
      do i=1,ITMAX !Iterate to convergence. 
      	an=-i*(i-a) 
            b=b+2. 
            d=an*d+b 
            if(abs(d).lt.FPMIN)d=FPMIN 
            c=b+an/c 
            if(abs(c).lt.FPMIN)c=FPMIN 
            d=1./d 
            del=d*c
            h=h*del 
            if(abs(del-1.).lt.EPS)goto 91 
	enddo
      write(*,*)  'a too large, ITMAX too small in gcf'  
   91 gammcf=exp(-x+a*log(x)-gln)*h !Put factors in front. 
   	return
  END SUBROUTINE gcf
  
  
  !USES gammp Returns the error function erf(x). 
  FUNCTION erf(x) 
  	REAL erf,x
      
      if(x.lt.0.)then 
      	erf=-gammp(.5,x**2) 
      else 
      	erf=gammp(.5,x**2) 
      endif 
      return 
  END FUNCTION erf
  
  
  double precision function stNormalCDF(x)
  	double precision x
     
     	if(x>6.) then
      	stNormalCDF=1.
      else if(x<-6.) then
      	stNormalCDF=0.
	else
      	stNormalCDF=0.5*(1.+erf(real(x)/1.41421356))
      end if
      
  end function stNormalCDF
  
!----------------------------------------------------------------------

  integer function binSearch(array,npt,x)
    implicit none
    integer npt,array(npt) !total no. of points in the array & the array
    integer x !no. whose position has to be found
    integer sp,ep !starting, end & insertion points

    binSearch = 1
    if(npt==0) then
      	return
    elseif(x>=array(1)) then
      	return
    end if
    if(x<=array(npt)) then
    	binSearch=npt+1
      return
    end if
    sp=1
    ep=npt
    do
    	if(sp==ep) return
	if(x>array((ep+sp)/2)) then
		ep=(ep+sp)/2
		binSearch=ep
	elseif(x<array((ep+sp)/2)) then
		sp=(ep+sp)/2+1
		binSearch=sp
	else
		ep=(ep+sp)/2
		binSearch=ep+1
            return
	end if
    end do
    
  end function binSearch

!----------------------------------------------------------------------

  subroutine piksrt(n,n1,arr,arr1)
  	integer n,n1
  	double precision arr(n)
  	double precision arr1(n1,n)
  	integer i,j
  	double precision a
  	double precision a1(n1)
  	do j=2,n
     		a=arr(j)
     		a1(1:n1)=arr1(1:n1,j)
     		do i=j-1,1,-1
        		if(arr(i).le.a) goto 10
        			arr(i+1)=arr(i)
				arr1(1:n1,i+1)=arr1(1:n1,i)
     		end do
     		i=0
10     	arr(i+1)=a
      	arr1(1:n1,i+1)=a1(1:n1)
  	end do
  	return
  end subroutine piksrt

!----------------------------------------------------------------------
  !calculation the multivariate normal function
  double precision function mNormalF(d,x,mu,C)
  	integer d !dimension
      double precision x(d) !data vector
      double precision mu(d) !mean
      double precision C(d,d) !covariance matrix
      double precision a(d) !residual vector (x-mu)
      integer i
      double precision Chisq,sqrD
      double precision TwoPi
      Parameter(TwoPi=6.283185307)
      double precision b(d)
      
      !DCHDC routine variables
      integer INFO
      double precision L(d,d) !cholesky decomposition matrix
      
      !compute the Cholesky decomposition of the covariance matrix 'C'
      !C = L^T.L
      !uses the LAPACK routine DPOTRF
      L=C
      call DPOTRF('L',d,L,d,INFO)
      
	!Compute chi-squared as the dot product b^T b, where
	!Lb=a and a is the vector of residuals
      !(x-mu)^T.C^-1.(x-mu)=b^T.b
      !where L.b=(x-mu) & C=L^T.L
	a(1:d)=x(1:d)-mu(1:d)
	do i=1,d
      	b(i)=(a(i)-sum(L(i,1:i-1)*b(1:i-1)))/L(i,i)
	end do
      Chisq=sum(b(1:d)**2.)
      
      !square root of det(C)
      sqrD=1.
      do i=1,d
      	sqrD=sqrD*L(i,i)
	end do
      
      !calculate the multivariate normal function
      mNormalF=exp(-Chisq/2.)/(sqrD*(TwoPi**(d/2.)))
      
  end function mNormalF

!----------------------------------------------------------------------
  !calculation the cumulative Gaussian mixture probability in 1-D
  double precision function cGaussMix(nClstr,pPt,pMean,pCov,wt)
  	integer nClstr !no. of cluster
      double precision pPt !point
      double precision pMean(nClstr),pCov(nClstr) !projected mean & variance of clusters
      double precision wt(nClstr)
      integer i
      double precision z
      
      cGaussMix=0.
      do i=1,nClstr
      	z=(pPt-pMean(i))/sqrt(pCov(i))
            cGaussMix=cGaussMix+wt(i)*stNormalCDF(z)
	end do
      
  end function cGaussMix

!---------------------------------------------------------------------- 
 
 logical function ptIn1Ell(ndim,pt,meanx,invcovx,kfacx)
 	
      implicit none
      
      integer ndim
      double precision pt(:),point(1,ndim)!point to be checked
      double precision meanx(:),invcovx(:,:),kfacx!cluster attributes
      double precision kfac!k factor of present point
      
      ptIn1Ell=.false.
        
      !if the cluster has vol=0 (i.e. kfacx=0) then point isn't in this cluster
      if(kfacx==0.d0) return
      
      point(1,:)=pt(:)
      call ScaleFactor(1,ndim,point,meanx,invcovx,kfac)
      if(kfac<kfacx .or. abs(kfac-kfacx)<0.0001) then
	ptIn1Ell=.true.
      endif
      
 end function ptIn1Ell
  
!----------------------------------------------------------------------
 
 subroutine ceffEvolve(ndim,vol,eff,npt)
 	
      implicit none
      
      !input variables
      integer ndim
      integer npt !no. of points in the mode
      
      !input/output variables
      double precision vol !ellipsoid volume
      double precision eff !enlargement
      
      !work variables
      double precision d1
      
      if(eff>1d0) then
      	d1=eff
      	eff=max(1d0,eff*exp(-2d0/dble(npt*ndim)))
      	vol=vol*((eff/d1)**(dble(ndim)/2d0))
      endif
      
      
 end subroutine ceffEvolve
  
!----------------------------------------------------------------------

end module utils1
