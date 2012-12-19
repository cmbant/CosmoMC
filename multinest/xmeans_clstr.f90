! X-means, Pelleg, Moore

module xmeans_clstr
  use kmeans_clstr
  use utils1
  implicit none
  
      integer n_dim
      double precision LogTwoPi
      parameter(LogTwoPi=1.83787706641)
      double precision larg
      parameter(larg=huge(1.d0))
      !information about ellipses at each node
      integer numClstr,nCls,ptClstrd,maxClstr !total clusters found yet & total pts in clstrs yet
      double precision, dimension(:,:), allocatable :: p,xclsMean,xclsEval,aux
      double precision, dimension(:), allocatable :: xclsBIC,loglike,xclsVar,xclsVol,xclsKfac,xclsEff,xclsDetcov
      double precision, dimension(:,:,:), allocatable :: xclsInvCov,xclsTMat,xclsCovmat,xclsEvec
      integer, dimension(:), allocatable :: revcPos,xclsPos,ptInClstr
      integer, dimension(:,:), allocatable :: pt_Clstr

 contains
 
  !simple X-means with likelihoods also being rearranged
  subroutine doXmeans(points,npt,np,nClstr,ptClstr,cMean,cCovmat,cInvCov,cEval,cEvec,like,min_pt)
	implicit none
      
      double precision points(:,:),like(:)
      integer npt,np !num of points & dimensionality
      integer nClstr !total clusters found
      integer ptClstr(:),min_pt
      double precision cMean(:,:),cCovmat(:,:,:),cInvCov(:,:,:),cEval(:,:),cEvec(:,:,:)
      integer i,j

	n_dim=np
      maxClstr=npt/min_pt+1
	allocate (p(n_dim,npt),loglike(npt),xclsBIC(maxClstr),revcPos(maxClstr), &
      	xclsPos(maxClstr),ptInClstr(maxClstr))
	allocate (xclsMean(maxClstr,n_dim),xclsCovmat(maxClstr,n_dim,n_dim), &
      	xclsInvCov(maxClstr,n_dim,n_dim),xclsEvec(maxClstr,n_dim,n_dim), &
      	xclsEval(maxClstr,n_dim))
      nCls=1
	numClstr=0
      ptClstrd=0
      ptInClstr=0
      revcPos=0
      
      call xmeans(points,npt,like,min_pt)
      !if(numClstr>1) call postprocess
      j=0
      do i=1,numClstr
      	if(ptInClstr(xclsPos(i))>0) then
            	j=j+1
      		ptClstr(j)=ptInClstr(xclsPos(i))
                  cMean(j,1:n_dim)=xclsMean(i,1:n_dim)
                  cEvec(j,1:n_dim,1:n_dim)=xclsEvec(i,1:n_dim,1:n_dim)
                  cEval(j,1:n_dim)=xclsEval(i,1:n_dim)
                  cCovmat(j,1:n_dim,1:n_dim)=xclsCovmat(i,1:n_dim,1:n_dim)
                  cInvCov(j,1:n_dim,1:n_dim)=xclsInvCov(i,1:n_dim,1:n_dim)
            end if
      end do
      nClstr=j
      points=p
      like=loglike
      deallocate(p,loglike,xclsBIC,revcPos,xclsPos,ptInClstr)
      deallocate(xclsMean,xclsInvCov,xclsCovmat,xclsEvec,xclsEval)
      
  end subroutine doXmeans
  
!----------------------------------------------------------------------
  !sub-clustering X-means with shared points
  !points rearranged with the cluster index of primary points also returned
  subroutine doXmeans5(points,npt,np,nClstr,ptClstr,ad,naux,auxa,min_pt)
      implicit none
      !input/output
      double precision points(:,:) !actual points which are rearranged after clustering
      			 !note: the clusters have shared points too
      double precision auxa(:,:) !auxilliary variables which are rearranged after clustering
      !input
      integer npt,np !num of points & dimensionality
      integer ad !no. of shared points required per cluster
      integer min_pt !min points allowed per cluster
      integer naux !no. of auxilliary variables to be a re-arranged
      !ouput
      integer nClstr !total clusters found
      integer ptClstr(:) !total no. of points in each cluster (including shared points)
      integer cluster(npt) !cluster index of primary points (before re-arrangement)
      
      !work
      integer i,j,maxec,k
      
      if(npt<2*min_pt) then
      	nClstr=1
            ptClstr(1)=npt
            ad=0
            return
	end if

      if(ad<0) ad=0
	n_dim=np
      maxec=npt/min_pt+1
      
	allocate (p(n_dim,npt+maxec*ad),revcPos(maxec),xclsPos(maxec), &
      	ptInClstr(maxec),aux(naux,npt+maxec*ad))
      
      ptInClstr=0
      nCls=1
	numClstr=0
      ptClstrd=0
      revcPos=0
      
      call xmeans5(points(:,1:npt),npt,cluster,ad,naux,auxa,min_pt)
      
      j=0
      k=0
      do i=1,numClstr
      	if(ptInClstr(xclsPos(i))>0) then
            	j=j+1
      		ptClstr(j)=ptInClstr(xclsPos(i))
                  k=k+ptInClstr(xclsPos(i))
            end if
      end do
      
      nClstr=j
      points(1:n_dim,1:k)=p(1:n_dim,1:k)
      auxa(1:naux,1:k)=aux(1:naux,1:k)
      
      deallocate(p,revcPos,xclsPos,ptInClstr,aux)
      
  end subroutine doXmeans5
      
!----------------------------------------------------------------------
  !sub-clustering X-means with 
  !nClstr = no. of clusters found by previous clustering
  !ptClstr(nClstr) = no. of pts in each cluster found by previous clustering
  !nb = every cluster is divided into nb sub-clusters
  !ad = shared pts per sub-cluster
  subroutine doXmeans8(points,like,npt,np,nClstr,ptClstr,ad,min_pt)
	implicit none
      
      double precision points(:,:),like(:)
      integer npt,np !num of points & dimensionality
      integer nClstr !total clusters found
      integer ptClstr(:),min_pt
      integer i,j,maxec,k
      integer ad !no. of shared points per cluster

      if(npt<2*min_pt) then
      	nClstr=1
            ptClstr(1)=npt
            ad=0
            return
	end if
      
      if(ad<0) ad=0

	n_dim=np
      maxec=npt/min_pt+1
      
	allocate (p(n_dim,npt+maxec*ad))
	allocate (loglike(npt+maxec*ad))
	allocate (revcPos(maxec))
	allocate (xclsPos(maxec))
	allocate (ptInClstr(maxec))
      
      ptInClstr=0
      nCls=1
	numClstr=0
      ptClstrd=0
      revcPos=0
      
      call Xmeans8(nClstr,points(:,1:npt),like(1:npt),ptClstr,npt,ad,min_pt)
      
      j=0
      k=0
      do i=1,numClstr
      	if(ptInClstr(xclsPos(i))>0) then
            	j=j+1
      		ptClstr(j)=ptInClstr(xclsPos(i))
                  k=k+ptInClstr(xclsPos(i))
            end if
      end do
      
      nClstr=j
      points(1:n_dim,1:k)=p(1:n_dim,1:k)
      like(1:k)=loglike(1:k)
      
      deallocate(p)
      deallocate(revcPos)
      deallocate(loglike)
      deallocate(xclsPos)
      deallocate(ptInClstr)
      
  end subroutine doXmeans8
      
!----------------------------------------------------------------------
 
  !simple X-means with likelihoods also being rearranged
  subroutine doXmeans6(points,npt,np,nClstr,ptClstr,naux,auxa,min_pt,maxC)
	implicit none
      
      double precision points(:,:),auxa(:,:)
      integer npt,np !num of points & dimensionality
      integer naux !dimension of auxilliary array
      integer nClstr !total clusters found
      integer ptClstr(:),min_pt,maxC
      integer i,j

      if(npt<2*min_pt) then
      	nClstr=1
            ptClstr(1)=npt
            return
	end if

	n_dim=np
      maxClstr=maxC
      allocate (p(n_dim,npt),aux(naux,npt),xclsPos(maxClstr),ptInClstr(maxClstr), &
      	xclsMean(maxClstr,n_dim),xclsVar(maxClstr),xclsBIC(maxClstr),revcPos(maxClstr))
      nCls=1
	numClstr=0
      ptClstrd=0
      ptInClstr=0
      
      call xmeans6(points,npt,naux,auxa,min_pt)
      
      j=0
      do i=1,numClstr
      	if(ptInClstr(xclsPos(i))>0) then
            	j=j+1
      		ptClstr(j)=ptInClstr(xclsPos(i))
            end if
      end do
      
      nClstr=j
      points(1:n_dim,1:npt)=p(1:n_dim,1:npt)
      auxa(1:naux,1:npt)=aux(1:naux,1:npt)
      deallocate(p,aux,xclsPos,ptInClstr,xclsMean,xclsVar,xclsBIC,revcPos)
      
  end subroutine doXmeans6
      
!----------------------------------------------------------------------
 
  !simple X-means
  subroutine doXmeans7(points,npt,np,nClstr,ptClstr,min_pt)
	implicit none
      
      double precision points(:,:)
      integer npt,np !num of points & dimensionality
      integer nClstr !total clusters found
      integer ptClstr(:),min_pt
      integer i,j

      if(npt<2*min_pt) then
      	nClstr=1
            ptClstr(1)=npt
            return
	end if

	n_dim=np
      maxClstr=npt/min_pt+1
      allocate (p(n_dim,npt),xclsPos(maxClstr),ptInClstr(maxClstr), &
      	xclsMean(maxClstr,n_dim),xclsVar(maxClstr),xclsBIC(maxClstr),revcPos(maxClstr))
      nCls=1
	numClstr=0
      ptClstrd=0
      ptInClstr=0
      
      call Xmeans7(points,npt,min_pt)
      
      j=0
      do i=1,numClstr
      	if(ptInClstr(xclsPos(i))>0) then
            	j=j+1
      		ptClstr(j)=ptInClstr(xclsPos(i))
            end if
      end do
      
      nClstr=j
      points(1:n_dim,1:npt)=p(1:n_dim,1:npt)
      deallocate(p,xclsPos,ptInClstr,xclsMean,xclsVar,xclsBIC,revcPos)
      
  end subroutine doXmeans7
      
!----------------------------------------------------------------------
  
  subroutine doGmeans(points,npt,np,nClstr,ptClstr,naux,auxa,min_pt,maxC)
	implicit none
      
      	double precision points(:,:),auxa(:,:)
      	integer npt,naux,np !num of points & dimensionality
      	integer nClstr !total clusters found
      	integer ptClstr(:),min_pt,maxC
      	integer i,j

      	if(npt<2*min_pt) then
      		nClstr=1
            	ptClstr(1)=npt
            	return
	endif

	n_dim=np
      	maxClstr=maxC
	allocate (p(n_dim,npt),aux(naux,npt),xclsPos(maxClstr),ptInClstr(maxClstr))
      	nCls=1
	numClstr=0
      	ptClstrd=0
      	ptInClstr=0
      
      	call Gmeans(points,npt,naux,auxa,min_pt)
      
      	j=0
      	do i=1,numClstr
      		if(ptInClstr(xclsPos(i))>0) then
            		j=j+1
      			ptClstr(j)=ptInClstr(xclsPos(i))
            	endif
      	enddo
      
      	nClstr=j
      	points(1:n_dim,1:npt)=p(1:n_dim,1:npt)
      	auxa(1:naux,1:npt)=aux(1:naux,1:npt)
      
      	deallocate(p,aux,xclsPos,ptInClstr)
      
  end subroutine doGmeans
      
!----------------------------------------------------------------------  
  
  subroutine doGmeans2(points,npt,np,nClstr,ptClstr,min_pt)
	implicit none
      
      double precision points(:,:)
      integer npt,np !num of points & dimensionality
      integer nClstr !total clusters found
      integer ptClstr(:),min_pt
      integer i,j

      if(npt<2*min_pt) then
      	nClstr=1
            ptClstr(1)=npt
            return
	end if

	n_dim=np
      maxClstr=npt/min_pt+1
	allocate (p(n_dim,npt),xclsPos(maxClstr),ptInClstr(maxClstr))
      nCls=1
	numClstr=0
      ptClstrd=0
      ptInClstr=0
      
      call Gmeans2(points,npt,min_pt)

      j=0
      do i=1,numClstr
      	if(ptInClstr(xclsPos(i))>0) then
            	j=j+1
      		ptClstr(j)=ptInClstr(xclsPos(i))
            end if
      end do
      nClstr=j
      points=p
      deallocate(p,xclsPos,ptInClstr)
      
  end subroutine doGmeans2
      
!----------------------------------------------------------------------  

  recursive subroutine Xmeans(pt,npt,like,min_pt)
    	implicit none
 
	integer min_pt,npt !num of points
    	double precision pt(:,:) !points
      double precision like(:) !likelihood values
    	double precision ptk(2,n_dim,npt)!points in clusters
    	double precision likek(2,npt)!loglike, to change order only
      integer nptk(2) !no. of points in the clusters
    	integer cluster(npt) !cluster array having cluster num of each pt
    	integer i,j,ip
      double precision covar(n_dim,n_dim),covark(2,n_dim,n_dim),invcov(n_dim,n_dim),invcovk(2,n_dim,n_dim)
    	double precision evec(n_dim,n_dim),eveck(2,n_dim,n_dim),eval(n_dim),evalk(2,n_dim),detcov,detcovk(2)
      double precision mean(n_dim),meank(2,n_dim)
      double precision BIC1,BIC2 !BIC for 1 & 2 clusters respectively
      integer i1
      logical flag
      
      flag=.false.
      
      do i=1,n_dim
		mean(i)=sum(pt(i,1:npt))/dble(npt)
      enddo
      !calculate cov,inv_cov,evec matrices & eval & det(cov)
      call calc_covmat(npt,n_dim,pt,mean,covar)
      evec=covar
      call diagonalize(evec,eval,n_dim,flag)
      !if any of the eval are 0 or -ve then use the lowest +ve eval
      do j=1,n_dim
		if(eval(j)<=0.) eval(1:j)=eval(j+1)
      enddo
      detcov=product(eval)
      call calc_invcovmat(n_dim,evec,eval,invcov)
      
      !don't cluster 
      !if it produces a cluster with less than 2*(min_pt-1)+1 points
      !since it would cluster into clusters having less than min_pt points
      !or
      !if the number of clusters has reached maxCls
      if(npt<2*min_pt.or.nCls==maxClstr) then
      	BIC2=larg
            BIC1=-larg
            flag=.true.
	end if
      
      if(.not.flag) then
      	!calculate BIC for all the points, no clustering
      	BIC1=calcBIC1(pt,npt,mean,invcov,detcov)
      
      	!breakup the points in 2 clusters
		i1=2
		meank(1,1:n_dim)=mean(1:n_dim)
      	call kmeans3(i1,pt,npt,n_dim,meank(1:i1,1:n_dim),cluster,min_pt)
      	!calculate the no. of points in each cluster & separate them
      	nptk=0
      	ptk=0.
      	do i=1,npt
            	nptk(cluster(i))=nptk(cluster(i))+1
            	ptk(cluster(i),:,nptk(cluster(i)))=pt(:,i)
            	likek(cluster(i),nptk(cluster(i)))=like(i)
		end do
      
      	!don't cluster if either of the clusters have less than min_pt points
      	!since that might be noise
      	do i=1,2
      		if(nptk(i)<min_pt) then
      			BIC2=larg
                        flag=.true.
                        exit
			end if
		end do
      end if
      
      if(.not.flag) then      
     		!calculate covar,inv_cov,evec matrices & eval & det(cov) for the points in each of the clusters
     		do i=1,2
      		!calculate covariance matrix
      		call calc_covmat(nptk(i),n_dim,ptk(i,:,1:nptk(i)),meank(i,:),covark(i,:,:))
			!diagonalize
            	eveck(i,:,:)=covark(i,:,:)
      		call diagonalize(eveck(i,:,:),evalk(i,:),n_dim,flag)
      
      		!if any of the eval are 0 or -ve then use the lowest +ve eval
     			do j=1,n_dim
				if(evalk(i,j)<=0.) evalk(i,1:j)=evalk(i,j+1)
      		enddo
      
      		detcovk(i)=product(evalk(i,:))
       		call calc_invcovmat(n_dim,eveck(i,:,:),evalk(i,:),invcovk(i,:,:))
      	end do

     		BIC2=calcBIC2(ptk(1:2,:,:),nptk(1:2),meank,invcovk,detcovk)
	end if
      
      if(BIC1>BIC2) then
		nCls=nCls+1
      	call xmeans(ptk(1,:,1:nptk(1)),nptk(1),likek(1,1:nptk(1)),min_pt)
      	call xmeans(ptk(2,:,1:nptk(2)),nptk(2),likek(2,1:nptk(2)),min_pt)
      	return
      else
      	p(:,ptClstrd+1:ptClstrd+npt)=pt(:,1:npt)
            loglike(ptClstrd+1:ptClstrd+npt)=like(1:npt)
      	ip=binSearch(ptInClstr,numClstr,npt)
            ptInClstr(ip+1:numClstr+1)=ptInClstr(ip:numClstr)
            ptInClstr(ip)=npt
            do i=1,numClstr
            	if(xclsPos(i)>=ip) xclsPos(i)=xclsPos(i)+1
            end do
            xclsPos(numClstr+1)=ip
            ptClstrd=ptClstrd+npt
            numClstr=numClstr+1
            !save cluster info
            xclsBIC(numClstr)=BIC1
            xclsMean(numClstr,1:n_dim)=mean(1:n_dim)
            xclsCovmat(numClstr,1:n_dim,1:n_dim)=covar(1:n_dim,1:n_dim)
            xclsInvCov(numClstr,1:n_dim,1:n_dim)=invcov(1:n_dim,1:n_dim)
            xclsEvec(numClstr,1:n_dim,1:n_dim)=evec(1:n_dim,1:n_dim)
            xclsEval(numClstr,1:n_dim)=eval(1:n_dim)
      end if
      
  end subroutine Xmeans 
      
!---------------------------------------------------------------------- 

  recursive subroutine Xmeans4(pt,npt,ad,min_pt)
    	implicit none
 
	integer min_pt,npt !num of points
      integer ad !no. of shared points required, returns the no. of shared points found
    	double precision pt(:,:) !points
    	double precision ptk(npt/(n_dim+1),n_dim,npt+ad)!points in clusters
      integer nptk(npt/(n_dim+1)) !no. of points in the clusters
    	integer cluster(npt) !cluster array having cluster num of each pt
    	integer i,j,ip,q
      double precision covar(n_dim,n_dim),covark(2,n_dim,n_dim),invcov(n_dim,n_dim),invcovk(2,n_dim,n_dim)
    	double precision evec(n_dim,n_dim),eveck(2,n_dim,n_dim),eval(n_dim),evalk(2,n_dim),detcov,detcovk(2)
      double precision mean(n_dim),meank(2,n_dim)
      double precision BIC1,BIC2 !BIC for 1 & 2 clusters respectively
      integer i1
      integer cluster2(npt/(n_dim+1),ad)
      logical flag
      
      flag=.false.
      
      do i=1,n_dim
		mean(i)=sum(pt(i,1:npt))/dble(npt)
      enddo
      
      !don't cluster if it produces a cluster with less than 2*(min_pt-1)+1 points
      !since it would cluster into clusters having less than min_pt points
      !don't cluster if the number of clusters has reached maxCls
      if(npt<2*min_pt .or. nCls==maxClstr) then
      	BIC2=larg
            BIC1=-larg
            flag=.true.
      end if
      
      if(.not.flag) then
      	!calculate cov,inv_cov,evec matrices & eval & det(cov)
      	call calc_covmat(npt,n_dim,pt,mean,covar)
      	evec=covar
      	call diagonalize(evec,eval,n_dim,flag)
      	!if any of the eval are 0 or -ve then use the lowest +ve eval
      	do j=1,n_dim
			if(eval(j)<=0.) eval(1:j)=eval(j+1)
      	enddo
      	detcov=product(eval)
      	call calc_invcovmat(n_dim,evec,eval,invcov)
      
      	!calculate BIC for all the points, no clustering
      	BIC1=calcBIC1(pt,npt,mean,invcov,detcov)

      
      	!breakup the points in 2 clusters
		i1=2
		meank(1,1:n_dim)=mean(1:n_dim)
      	call kmeans3(i1,pt,npt,n_dim,meank(1:i1,1:n_dim),cluster,min_pt)
      	!calculate the no. of points in each cluster & separate them
      	nptk=0
      	ptk=0.
      	do i=1,npt
            	nptk(cluster(i))=nptk(cluster(i))+1
            	ptk(cluster(i),:,nptk(cluster(i)))=pt(:,i)
		end do
      
      	!don't cluster if either of the clusters have less than min_pt points
      	!since that might be noise
      	do i=1,2
      		if(nptk(i)<min_pt) then
      			BIC2=larg
            		flag=.true.
                        exit
			end if
		end do
      end if
                  
	if(.not.flag) then
     		!calculate covar,inv_cov,evec matrices & eval & det(cov) for the points in each of the clusters
     		do i=1,2
      		!calculate covariance matrix
            	call calc_covmat(nptk(i),n_dim,ptk(i,:,1:nptk(i)),meank(i,:),covark(i,:,:))
			!diagonalize
            	eveck(i,:,:)=covark(i,:,:)
      		call diagonalize(eveck(i,:,:),evalk(i,:),n_dim,flag)
      
      		!if any of the eval are 0 or -ve then use the lowest +ve eval
     			do j=1,n_dim
				if(evalk(i,j)<=0.) evalk(i,1:j)=evalk(i,j+1)
      		enddo
      
      		detcovk(i)=product(evalk(i,:))
       		call calc_invcovmat(n_dim,eveck(i,:,:),evalk(i,:),invcovk(i,:,:))
            
      	end do
     		BIC2=calcBIC2(ptk(1:2,:,:),nptk(1:2),meank,invcovk,detcovk)
	end if
      
      if(BIC1>BIC2) then
		nCls=nCls+1
            call xmeans4(ptk(1,:,1:nptk(1)),nptk(1),ad,min_pt)
		call xmeans4(ptk(2,:,1:nptk(2)),nptk(2),ad,min_pt)
      	return
      else
      	!i1=min(nb,npt/max(min_pt,min(10,2*min_pt)))
            i1=npt/min_pt+1
      	!if(npt>=2*max(min_pt,min(10,2*min_pt)).and.i1>1) then
      	if(npt>=2*min_pt.and.i1>1) then
            	!break up the cluster into nb clusters
                  call kmeans7(i1,pt,npt,n_dim,cluster,ad,cluster2,min_pt)
			!calculate the no. of points in each cluster & separate them
                  
      		nptk=0
      		ptk=0.
      		do i=1,npt
            		nptk(cluster(i))=nptk(cluster(i))+1
            		ptk(cluster(i),:,nptk(cluster(i)))=pt(:,i)
			end do
                  
                  do i=1,i1
                  	do j=1,ad
                              !if(cluster2(i,j)==0) exit
					nptk(i)=nptk(i)+1
                              ptk(i,:,nptk(i))=pt(:,cluster2(i,j))
				end do
			end do
                  	
                  !arrange the points
                  do i=1,i1
                  	p(:,ptClstrd+1:ptClstrd+nptk(i))=ptk(i,:,1:nptk(i))
      			ip=binSearch(ptInClstr,numClstr,nptk(i))
           	 		ptInClstr(ip+1:numClstr+1)=ptInClstr(ip:numClstr)
            		ptInClstr(ip)=nptk(i)
           			do q=1,numClstr
            			if(xclsPos(q)>=ip) xclsPos(q)=xclsPos(q)+1
            		end do
            		xclsPos(numClstr+1)=ip
            		ptClstrd=ptClstrd+nptk(i)
            		numClstr=numClstr+1
            		!save cluster info
                        do q=1,n_dim
            			xclsMean(numClstr,q)=sum(ptk(i,q,1:nptk(i)))/dble(nptk(i))
                        end do
			end do
		else
      		p(:,ptClstrd+1:ptClstrd+npt)=pt(:,1:npt)
      		ip=binSearch(ptInClstr,numClstr,npt)
            	ptInClstr(ip+1:numClstr+1)=ptInClstr(ip:numClstr)
            	ptInClstr(ip)=npt
            	do i=1,numClstr
            		if(xclsPos(i)>=ip) xclsPos(i)=xclsPos(i)+1
            	end do
            	xclsPos(numClstr+1)=ip
            	ptClstrd=ptClstrd+npt
            	numClstr=numClstr+1
            	!save cluster info
            	xclsInvCov(numClstr,1:n_dim,1:n_dim)=invcov(1:n_dim,1:n_dim)
            	xclsBIC(numClstr)=BIC1
            	xclsMean(numClstr,1:n_dim)=mean(1:n_dim)
		end if
      end if
      
  end subroutine xmeans4 
      
!---------------------------------------------------------------------- 

  recursive subroutine Xmeans5(ptf,tnpt,cluster,ad,naux,auxa,min_pt)
    	implicit none
 
	integer min_pt,npt !num of points
      integer ad !no. of shared points required, returns the no. of shared points found
      integer tnpt
    	double precision ptf(:,:),pt(n_dim,tnpt) !points
    	double precision ptk(tnpt/(n_dim+1),n_dim,tnpt+ad)!points in clusters
      integer nptk(tnpt/(n_dim+1)) !no. of points in the clusters
    	integer cluster(tnpt) !cluster array having cluster num of each pt
    	integer i,j,ip
      integer i1
      integer cluster2(tnpt/(n_dim+1),ad)
      integer naux
      double precision auxa(:,:),auxk(tnpt/(n_dim+1),naux,tnpt+ad)
      
      npt=tnpt
      pt(:,1:npt)=ptf(:,1:npt)
            
      i1=npt/min_pt+1
      	
      if(npt>=2*min_pt.and.i1>1) then
            !break up the cluster into nb clusters
		call kmeans7(i1,pt(:,1:npt),npt,n_dim,cluster(1:npt),ad,cluster2,min_pt)
		!calculate the no. of points in each cluster & separate them
		nptk=0
      	do i=1,npt
            	nptk(cluster(i))=nptk(cluster(i))+1
            	ptk(cluster(i),:,nptk(cluster(i)))=pt(:,i)
            	auxk(cluster(i),:,nptk(cluster(i)))=auxa(:,i)
		end do
                  
		do i=1,i1
			do j=1,ad
				!if(cluster2(i,j)==0) exit
				nptk(i)=nptk(i)+1
				ptk(i,:,nptk(i))=pt(:,cluster2(i,j))
				auxk(i,:,nptk(i))=auxa(:,cluster2(i,j))
			end do
		end do
                  	
		!arrange the points
		do i=1,i1
                  p(:,ptClstrd+1:ptClstrd+nptk(i))=ptk(i,:,1:nptk(i))
                  aux(:,ptClstrd+1:ptClstrd+nptk(i))=auxk(i,:,1:nptk(i))
      		!ip=binSearch(ptInClstr,numClstr,nptk(i))
			ip=numClstr+1
			!ptInClstr(ip+1:numClstr+1)=ptInClstr(ip:numClstr)
            	ptInClstr(ip)=nptk(i)
			!do q=1,numClstr
            		!if(xclsPos(q)>=ip) xclsPos(q)=xclsPos(q)+1
            	!end do
            	xclsPos(numClstr+1)=ip
            	ptClstrd=ptClstrd+nptk(i)
            	numClstr=numClstr+1
		end do
	else
      	p(:,ptClstrd+1:ptClstrd+npt)=pt(:,1:npt)
      	aux(:,ptClstrd+1:ptClstrd+npt)=auxa(:,1:npt)
      	!ip=binSearch(ptInClstr,numClstr,npt)
		ip=numClstr+1
		!ptInClstr(ip+1:numClstr+1)=ptInClstr(ip:numClstr)
            ptInClstr(ip)=npt
            !do i=1,numClstr
			!if(xclsPos(i)>=ip) xclsPos(i)=xclsPos(i)+1
            !end do
            xclsPos(numClstr+1)=ip
            ptClstrd=ptClstrd+npt
            numClstr=numClstr+1
	end if
      
  end subroutine xmeans5 
      
!----------------------------------------------------------------------   

  recursive subroutine Xmeans8(nClsf,ptf,lkf,ptInClsf,tnpt,ad,min_pt)
    	implicit none
 
	integer min_pt,npt !num of points
      integer ad  !no. of shared points required, returns the no. of shared points found
      integer tnpt,nClsf,ptInClsf(:)
    	double precision ptf(:,:),pt(n_dim,tnpt),lkf(:),lk(tnpt) !points & log-like
    	double precision ptk(tnpt/(n_dim+1),n_dim,tnpt+ad)!points in clusters
    	double precision likek(tnpt/(n_dim+1),tnpt+ad)!log-like of points in clusters
      integer nptk(tnpt/(n_dim+1)) !no. of points in the clusters
    	integer cluster(tnpt) !cluster array having cluster num of each pt
    	integer i,j,ip,q,l,m
      integer i1
      integer cluster2(tnpt/(n_dim+1),ad)
      
      m=1
	do l=1,nClsf
      	npt=ptInClsf(l)
            pt(:,1:npt)=ptf(:,m:m+ptInClsf(l))
            lk(1:npt)=lkf(m:m+ptInClsf(l))
            m=m+ptInClsf(l)
            
      	i1=npt/min_pt+1
      	
            if(npt>=2*min_pt.and.i1>1) then
            	!break up the cluster into nb clusters
                  call kmeans7(i1,pt(:,1:npt),npt,n_dim,cluster(1:npt),ad,cluster2,min_pt)
			!calculate the no. of points in each cluster & separate them
                  
      		nptk=0
      		ptk=0.
      		do i=1,npt
            		nptk(cluster(i))=nptk(cluster(i))+1
            		ptk(cluster(i),:,nptk(cluster(i)))=pt(:,i)
            		likek(cluster(i),nptk(cluster(i)))=lk(i)
			end do
                  
                  if(ad>0) then
                  	do i=1,i1
                  		do j=1,ad
                              	!if(cluster2(i,j)==0) exit
                                    nptk(i)=nptk(i)+1
                              	ptk(i,:,nptk(i))=pt(:,cluster2(i,j))
                              	likek(i,nptk(i))=lk(cluster2(i,j))
                              end do
				end do
                  end if
                  	
                  !arrange the points
                  do i=1,i1
                  	p(:,ptClstrd+1:ptClstrd+nptk(i))=ptk(i,:,1:nptk(i))
                  	loglike(ptClstrd+1:ptClstrd+nptk(i))=likek(i,1:nptk(i))
      			ip=binSearch(ptInClstr,numClstr,nptk(i))
           	 		ptInClstr(ip+1:numClstr+1)=ptInClstr(ip:numClstr)
            		ptInClstr(ip)=nptk(i)
           			do q=1,numClstr
            			if(xclsPos(q)>=ip) xclsPos(q)=xclsPos(q)+1
            		end do
            		xclsPos(numClstr+1)=ip
            		ptClstrd=ptClstrd+nptk(i)
            		numClstr=numClstr+1
			end do
		else
      		p(:,ptClstrd+1:ptClstrd+npt)=pt(:,1:npt)
      		loglike(ptClstrd+1:ptClstrd+npt)=lk(1:npt)
      		ip=binSearch(ptInClstr,numClstr,npt)
            	ptInClstr(ip+1:numClstr+1)=ptInClstr(ip:numClstr)
            	ptInClstr(ip)=npt
            	do i=1,numClstr
            		if(xclsPos(i)>=ip) xclsPos(i)=xclsPos(i)+1
            	end do
            	xclsPos(numClstr+1)=ip
            	ptClstrd=ptClstrd+npt
            	numClstr=numClstr+1
		end if
      end do
      
  end subroutine Xmeans8 
      
!----------------------------------------------------------------------   

  recursive subroutine Xmeans6(pt,npt,naux,auxa,min_pt)
    	implicit none
 
	integer min_pt,npt,naux !num of points
    	double precision pt(:,:) !points
      double precision auxa(:,:) !auxilliary array
    	double precision ptk(2,n_dim,npt)!points in clusters
    	double precision auxk(2,naux,npt)!aux, to change order only
      integer nptk(2) !no. of points in the clusters
    	integer cluster(npt) !cluster array having cluster num of each pt
      double precision mean(3,n_dim),var(3)
    	integer i,j,ip
      double precision BIC1,BIC2 !BIC for 1 & 2 clusters respectively
      integer i1
      logical flag
      
      flag=.false.
      
      !don't cluster 
      !if it produces a cluster with less than 2*(min_pt-1)+1 points
      !since it would cluster into clusters having less than min_pt points
      !or
      !if the number of clusters has reached maxCls
      if(npt<2*min_pt.or.nCls==maxClstr) then
      	BIC2=larg
            BIC1=-larg
            flag=.true.
	end if
      
      if(.not.flag) then
      	var=0.
      	!calculate mean & variance
      	do i=1,n_dim
      		mean(3,i)=sum(pt(i,1:npt))/npt
            	var(3)=var(3)+sum(pt(i,1:npt)**2.)/npt-mean(3,i)**2.
     	 	end do
      
      	!calculate BIC for all the points, no clustering
      	BIC1=calcBIC1_iso(npt,var(3))
      
      	!breakup the points in 2 clusters
		i1=2
		mean(1,:)=mean(3,:)
      	call kmeans3(i1,pt(1:n_dim,1:npt),npt,n_dim,mean(1:i1,1:n_dim),cluster,min_pt)
      	!calculate the no. of points in each cluster & separate them
      	nptk=0
      	ptk=0.
      	do i=1,npt
            	nptk(cluster(i))=nptk(cluster(i))+1
            	ptk(cluster(i),:,nptk(cluster(i)))=pt(:,i)
            	auxk(cluster(i),:,nptk(cluster(i)))=auxa(:,i)
		end do
      
      	!don't cluster if either of the clusters have less than min_pt points
      	!since that might be noise
      	do i=1,2
      		if(nptk(i)<min_pt) then
      			BIC2=larg
                        flag=.true.
                        exit
			end if
		end do
	end if
      
      if(.not.flag) then
      	!calculate mean & variance
		var(1:2)=0.
      	do i=1,2
      		do j=1,n_dim
            		var(i)=var(i)+sum(ptk(i,j,1:nptk(i))**2.)/nptk(i)-mean(i,j)**2.
			end do
      	end do
      
      	BIC2=calcBIC2_iso(nptk(1:2),var(1:2))
	end if
      
      if(BIC1>BIC2) then
		nCls=nCls+1
      	call Xmeans6(ptk(1,:,1:nptk(1)),nptk(1),naux,auxk(1,:,1:nptk(1)),min_pt)
      	call Xmeans6(ptk(2,:,1:nptk(2)),nptk(2),naux,auxk(2,:,1:nptk(2)),min_pt)
      	return
      else
      	p(:,ptClstrd+1:ptClstrd+npt)=pt(:,1:npt)
            aux(:,ptClstrd+1:ptClstrd+npt)=auxa(:,1:npt)
      	ip=binSearch(ptInClstr,numClstr,npt)
            ptInClstr(ip+1:numClstr+1)=ptInClstr(ip:numClstr)
            ptInClstr(ip)=npt
            do i=1,numClstr
            	if(xclsPos(i)>=ip) xclsPos(i)=xclsPos(i)+1
            end do
            xclsPos(numClstr+1)=ip
            xclsMean(numClstr+1,1:n_dim)=mean(3,1:n_dim)
            xclsVar(numClstr+1)=var(3)
            xclsBIC(numClstr+1)=BIC1
            ptClstrd=ptClstrd+npt
            numClstr=numClstr+1
      end if
      
  end subroutine Xmeans6
      
!----------------------------------------------------------------------  

  recursive subroutine Xmeans7(pt,npt,min_pt)
    	implicit none
 
	integer min_pt,npt !num of points
    	double precision pt(:,:) !points
    	double precision ptk(2,n_dim,npt)!points in clusters
      integer nptk(2) !no. of points in the clusters
    	integer cluster(npt) !cluster array having cluster num of each pt
      double precision mean(3,n_dim),var(3)
    	integer i,j,ip
      double precision BIC1,BIC2 !BIC for 1 & 2 clusters respectively
      integer i1
      logical flag
      
      flag=.false.
      
      !don't cluster 
      !if it produces a cluster with less than 2*(min_pt-1)+1 points
      !since it would cluster into clusters having less than min_pt points
      !or
      !if the number of clusters has reached maxCls
      if(npt<2*min_pt.or.nCls==maxClstr) then
      	BIC2=larg
            BIC1=-larg
            flag=.true.
	end if
      
      if(.not.flag) then
      	var=0.
      	!calculate mean & variance
      	do i=1,n_dim
      		mean(3,i)=sum(pt(i,1:npt))/npt
            	var(3)=var(3)+sum(pt(i,1:npt)**2.)/npt-mean(3,i)**2.
     	 	end do
      
      	!calculate BIC for all the points, no clustering
      	BIC1=calcBIC1_iso(npt,var(3))
      
      	!breakup the points in 2 clusters
		i1=2
		mean(1,:)=mean(3,:)
      	call kmeans3(i1,pt(1:n_dim,1:npt),npt,n_dim,mean(1:i1,1:n_dim),cluster,min_pt)
      	!calculate the no. of points in each cluster & separate them
      	nptk=0
      	ptk=0.
      	do i=1,npt
            	nptk(cluster(i))=nptk(cluster(i))+1
            	ptk(cluster(i),:,nptk(cluster(i)))=pt(:,i)
		end do
      
      	!don't cluster if either of the clusters have less than min_pt points
      	!since that might be noise
      	do i=1,2
      		if(nptk(i)<min_pt) then
      			BIC2=larg
                        flag=.true.
                        exit
			end if
		end do
	end if
      
      if(.not.flag) then
      	!calculate mean & variance
		var(1:2)=0.
      	do i=1,2
      		do j=1,n_dim
            		var(i)=var(i)+sum(ptk(i,j,1:nptk(i))**2.)/nptk(i)-mean(i,j)**2.
			end do
      	end do
      
      	BIC2=calcBIC2_iso(nptk(1:2),var(1:2))
	end if
      
      if(BIC1>BIC2) then
		nCls=nCls+1
      	call Xmeans7(ptk(1,:,1:nptk(1)),nptk(1),min_pt)
      	call Xmeans7(ptk(2,:,1:nptk(2)),nptk(2),min_pt)
      	return
      else
      	p(:,ptClstrd+1:ptClstrd+npt)=pt(:,1:npt)
      	ip=binSearch(ptInClstr,numClstr,npt)
            ptInClstr(ip+1:numClstr+1)=ptInClstr(ip:numClstr)
            ptInClstr(ip)=npt
            do i=1,numClstr
            	if(xclsPos(i)>=ip) xclsPos(i)=xclsPos(i)+1
            end do
            xclsPos(numClstr+1)=ip
            xclsMean(numClstr+1,1:n_dim)=mean(3,1:n_dim)
            xclsVar(numClstr+1)=var(3)
            xclsBIC(numClstr+1)=BIC1
            ptClstrd=ptClstrd+npt
            numClstr=numClstr+1
      end if
      
  end subroutine Xmeans7
      
!---------------------------------------------------------------------- 
  !sub-clustering X-means with shared points
  !points rearranged with the cluster index of primary points also returned
  subroutine doXmeans9(points,npt,np,nClstr,ptClstr,ad,min_pt)
	implicit none
      !input/output
      double precision points(:,:) !actual points which are rearranged after clustering
      			 !note: the clusters have shared points too
      !input
      integer npt,np !num of points & dimensionality
      integer ad !no. of shared points required per cluster
      integer min_pt !min points allowed per cluster
      !ouput
      integer nClstr !total clusters found
      integer ptClstr(:) !total no. of points in each cluster (including shared points)
      integer cluster(npt) !cluster index of primary points (before re-arrangement)
      
      !work
      integer i,j,maxec,k
      
      if(npt<2*min_pt) then
      	nClstr=1
            ptClstr(1)=npt
            ad=0
            return
	end if

      if(ad<0) ad=0
	n_dim=np
      maxec=npt/min_pt+1
      
	allocate (p(n_dim,npt+maxec*ad))
	allocate (revcPos(maxec))
	allocate (xclsPos(maxec))
	allocate (ptInClstr(maxec))
      
      ptInClstr=0
      nCls=1
	numClstr=0
      ptClstrd=0
      revcPos=0
      
      call Xmeans9(points(:,1:npt),npt,cluster,ad,min_pt)
      
      j=0
      k=0
      do i=1,numClstr
      	if(ptInClstr(xclsPos(i))>0) then
            	j=j+1
      		ptClstr(j)=ptInClstr(xclsPos(i))
                  k=k+ptInClstr(xclsPos(i))
            end if
      end do
      
      nClstr=j
      points(1:n_dim,1:k)=p(1:n_dim,1:k)
      
      deallocate(p)
      deallocate(revcPos)
      deallocate(xclsPos)
      deallocate(ptInClstr)
      
  end subroutine doXmeans9
      
!----------------------------------------------------------------------

  recursive subroutine Xmeans9(ptf,tnpt,cluster,ad,min_pt)
    	implicit none
 
	integer min_pt,npt !num of points
      integer ad !no. of shared points required, returns the no. of shared points found
      integer tnpt
    	double precision ptf(:,:),pt(n_dim,tnpt) !points
    	double precision ptk(tnpt/(n_dim+1),n_dim,tnpt+ad)!points in clusters
      integer nptk(tnpt/(n_dim+1)) !no. of points in the clusters
    	integer cluster(tnpt) !cluster array having cluster num of each pt
    	integer i,j,ip
      integer i1
      integer cluster2(tnpt/(n_dim+1),ad)
      
      npt=tnpt
      pt(:,1:npt)=ptf(:,1:npt)
            
      i1=npt/min_pt+1
      	
      if(npt>=2*min_pt.and.i1>1) then
            !break up the cluster into nb clusters
		call kmeans7(i1,pt(:,1:npt),npt,n_dim,cluster(1:npt),ad,cluster2,min_pt)
		!calculate the no. of points in each cluster & separate them
		nptk=0
      	do i=1,npt
            	nptk(cluster(i))=nptk(cluster(i))+1
            	ptk(cluster(i),:,nptk(cluster(i)))=pt(:,i)
		end do
                  
		do i=1,i1
			do j=1,ad
				!if(cluster2(i,j)==0) exit
				nptk(i)=nptk(i)+1
				ptk(i,:,nptk(i))=pt(:,cluster2(i,j))
			end do
		end do
                  	
		!arrange the points
		do i=1,i1
                  p(:,ptClstrd+1:ptClstrd+nptk(i))=ptk(i,:,1:nptk(i))
      		!ip=binSearch(ptInClstr,numClstr,nptk(i))
			ip=numClstr+1
			!ptInClstr(ip+1:numClstr+1)=ptInClstr(ip:numClstr)
            	ptInClstr(ip)=nptk(i)
			!do q=1,numClstr
            		!if(xclsPos(q)>=ip) xclsPos(q)=xclsPos(q)+1
            	!end do
            	xclsPos(numClstr+1)=ip
            	ptClstrd=ptClstrd+nptk(i)
            	numClstr=numClstr+1
		end do
	else
      	p(:,ptClstrd+1:ptClstrd+npt)=pt(:,1:npt)
      	!ip=binSearch(ptInClstr,numClstr,npt)
		ip=numClstr+1
		!ptInClstr(ip+1:numClstr+1)=ptInClstr(ip:numClstr)
            ptInClstr(ip)=npt
            !do i=1,numClstr
			!if(xclsPos(i)>=ip) xclsPos(i)=xclsPos(i)+1
            !end do
            xclsPos(numClstr+1)=ip
            ptClstrd=ptClstrd+npt
            numClstr=numClstr+1
	end if
      
  end subroutine Xmeans9 
      
!----------------------------------------------------------------------   

  recursive subroutine Gmeans(pt,npt,naux,auxa,min_pt)
    	implicit none
 
	integer min_pt,npt,naux !num of points
    	double precision pt(:,:) !points
      double precision auxa(:,:) !auxilliary points
    	double precision ptk(2,n_dim,npt)!points in clusters
    	double precision auxk(2,naux,npt)!loglike, to change order only
      integer nptk(2) !no. of points in the clusters
    	integer cluster(npt) !cluster array having cluster num of each pt
    	integer i,ip
      integer i1
      double precision alpha !confidence level for Anderson-Darlind test
      parameter(alpha=0.0001)
      double precision mean(n_dim),meank(2,n_dim)
      double precision delMean(n_dim) !difference between the means of the two clusters
      logical flag
      
      flag=.false.
      
      !calculate the mean
      do i=1,n_dim
		mean(i)=sum(pt(i,1:npt))/npt
      enddo
      
      !don't cluster if it produces a cluster with less than 2*(min_pt-1)+1 points
      !since it would cluster into clusters having less than min_pt points
      if(npt<2*min_pt.or.nCls==maxClstr) flag=.true.
      
      if(.not.flag) then
      	!breakup the points in 2 clusters
		i1=2
      	meank(1,1:n_dim)=mean(1:n_dim)
      	call kmeans3(i1,pt,npt,n_dim,meank(1:i1,1:n_dim),cluster,min_pt)
      	!calculate the no. of points in each cluster & separate them
      	nptk=0
      	do i=1,npt
           	 	nptk(cluster(i))=nptk(cluster(i))+1
            	ptk(cluster(i),:,nptk(cluster(i)))=pt(:,i)
            	auxk(cluster(i),:,nptk(cluster(i)))=auxa(:,i)
		end do
      
      	!don't cluster if either of the clusters have less than min_pt points
      	!since that might be noise
      	do i=1,2
      		if(nptk(i)<min_pt) then
                  	flag=.true.
                        exit
			end if
		end do
	end if
      
      if(.not.flag) then
      	!calculate the means of the two clusters & their difference
      	delMean=meank(1,1:n_dim)-meank(2,1:n_dim)
            flag=AndersonDarling(npt,pt,delMean,alpha)
	end if
      
      if(.not.flag) then
		nCls=nCls+1
      	call Gmeans(ptk(1,:,1:nptk(1)),nptk(1),naux,auxk(1,:,1:nptk(1)),min_pt)
      	call Gmeans(ptk(2,:,1:nptk(2)),nptk(2),naux,auxk(2,:,1:nptk(2)),min_pt)
      	return
      else
      	p(:,ptClstrd+1:ptClstrd+npt)=pt(:,1:npt)
            aux(:,ptClstrd+1:ptClstrd+npt)=auxa(:,1:npt)
      	ip=binSearch(ptInClstr,numClstr,npt)
            ptInClstr(ip+1:numClstr+1)=ptInClstr(ip:numClstr)
            ptInClstr(ip)=npt
            do i=1,numClstr
            	if(xclsPos(i)>=ip) xclsPos(i)=xclsPos(i)+1
            end do
            xclsPos(numClstr+1)=ip
            ptClstrd=ptClstrd+npt
            numClstr=numClstr+1
      end if
      
  end subroutine Gmeans 
      
!----------------------------------------------------------------------

  recursive subroutine Gmeans2(pt,npt,min_pt)
    	implicit none
 
	integer min_pt,npt !num of points
    	double precision pt(:,:) !points
    	double precision ptk(2,n_dim,npt)!points in clusters
      integer nptk(2) !no. of points in the clusters
    	integer cluster(npt) !cluster array having cluster num of each pt
    	integer i,ip
      integer i1
      double precision alpha !confidence level for Anderson-Darlind test
      parameter(alpha=0.0001)
      double precision mean(n_dim),meank(2,n_dim)
      double precision delMean(n_dim) !difference between the means of the two clusters
      logical flag
      
      flag=.false.
      
      !calculate the mean
      do i=1,n_dim
		mean(i)=sum(pt(i,1:npt))/npt
      enddo
      
      !don't cluster if it produces a cluster with less than 2*(min_pt-1)+1 points
      !since it would cluster into clusters having less than min_pt points
      if(npt<2*min_pt.or.nCls==maxClstr) flag=.true.
      
      if(.not.flag) then
      	!breakup the points in 2 clusters
		i1=2
      	meank(1,1:n_dim)=mean(1:n_dim)
      	call kmeans3(i1,pt,npt,n_dim,meank(1:i1,1:n_dim),cluster,min_pt)
      	!calculate the no. of points in each cluster & separate them
      	nptk=0
      	do i=1,npt
           	 	nptk(cluster(i))=nptk(cluster(i))+1
            	ptk(cluster(i),:,nptk(cluster(i)))=pt(:,i)
		end do
      
      	!don't cluster if either of the clusters have less than min_pt points
      	!since that might be noise
      	do i=1,2
      		if(nptk(i)<min_pt) then
                  	flag=.true.
                        exit
			end if
		end do
	end if
      
      if(.not.flag) then
      	!calculate the means of the two clusters & their difference
      	delMean=meank(1,1:n_dim)-meank(2,1:n_dim)
            flag=AndersonDarling(npt,pt,delMean,alpha)
	end if
      
      if(.not.flag) then
		nCls=nCls+1
      	call Gmeans2(ptk(1,:,1:nptk(1)),nptk(1),min_pt)
      	call Gmeans2(ptk(2,:,1:nptk(2)),nptk(2),min_pt)
      	return
      else
      	p(:,ptClstrd+1:ptClstrd+npt)=pt(:,1:npt)
      	ip=binSearch(ptInClstr,numClstr,npt)
            ptInClstr(ip+1:numClstr+1)=ptInClstr(ip:numClstr)
            ptInClstr(ip)=npt
            do i=1,numClstr
            	if(xclsPos(i)>=ip) xclsPos(i)=xclsPos(i)+1
            end do
            xclsPos(numClstr+1)=ip
            ptClstrd=ptClstrd+npt
            numClstr=numClstr+1
      end if
      
  end subroutine Gmeans2 
      
!----------------------------------------------------------------------

  !calculate BIC for no clustering in the data
  double precision function calcBIC1(pt,npt,mean,invcov,detcov)
	implicit none
      integer npt !num of points
      double precision pt(:,:)
      double precision detcov !determinant of the covariance matrix
      double precision mean(:),invcov(:,:)
      double precision tpt(n_dim,npt),a,nd,nn
	integer i,j,k
      
      nd=n_dim
      nn=npt
      
      !calculate point-mean
      do i=1,npt
      	do j=1,n_dim
            	tpt(j,i)=pt(j,i)-mean(j)
            end do
      end do
            	
      !a=sum(Transpose(pt(i)-mean).invcov.(pt(i)-mean))
      a=0.
      do i=1,npt
      	do j=1,n_dim
            	a=a+tpt(j,i)*tpt(j,i)*invcov(j,j)
            	do k=j+1,n_dim
            		a=a+2.*tpt(j,i)*tpt(k,i)*invcov(j,k)
                  end do
		end do
	end do
	calcBIC1=nn*nd*LogTwoPi+nn*log(detcov)+a+0.5*nd*(nd+3.)*log(nn)
  	
  end function calcBIC1
      
!----------------------------------------------------------------------

  !calculate BIC for no clustering in the data
  double precision function calcBIC1_iso(npt,variance)
	implicit none
      integer npt !num of points
      double precision variance,var
      double precision nd,nn
      
      nd=n_dim
      nn=npt
      var=variance
      
      !convert into unbiased variance
      var=var*nn/(nn-1.)
      
      calcBIC1_iso=nn*LogTwoPi+nn*nd*log(var)+nn-1.+(nd+1.)*log(nn)
  	
  end function calcBIC1_iso
      
!----------------------------------------------------------------------

  !calculate BIC for clustering in the data
  double precision function calcBIC2(ptk,nptk,meank,invcovk,detcovk)
	implicit none
      integer nptk(2) !num of points
      double precision ptk(2,n_dim,nptk(1)+nptk(2))
      double precision detcovk(2) !determinant of the covariance matrix
      double precision meank(2,n_dim),invcovk(2,n_dim,n_dim)
      double precision tptk(2,n_dim,nptk(1)+nptk(2)),ak(2)
	integer i,j,k,l
      double precision nd,np1,np2,npt
      
      np1=nptk(1)
      np2=nptk(2)
      npt=np1+np2
      nd=n_dim
      
      !calculate point-mean
	do l=1,2
      	do i=1,nptk(l)
      		do j=1,n_dim
            		tptk(l,j,i)=ptk(l,j,i)-meank(l,j)
            	end do
      	end do
      end do
      
      !a=sum(Transpose(pt(i)-mean).invcov.(pt(i)-mean))
      ak=0.
      do l=1,2
      	do i=1,nptk(l)
      		do j=1,n_dim
            		ak(l)=ak(l)+tptk(l,j,i)*tptk(l,j,i)*invcovk(l,j,j)
                        do k=j+1,n_dim
            			ak(l)=ak(l)+2.*tptk(l,j,i)*tptk(l,k,i)*invcovk(l,j,k)
                  	end do
			end do
		end do
	end do
      
	calcBIC2=npt*nd*LogTwoPi+np1*log(detcovk(1))+np2*log(detcovk(2))+ak(1)+ak(2) &
      -2.*(np1*log(np1)+np2*log(np2))+2.*npt*log(npt)+(nd**2.+3.*nd+1.)*log(npt)
  end function calcBIC2
      
!----------------------------------------------------------------------

  !calculate BIC for clustering in the data
  double precision function calcBIC2_iso(nptk,var)
	implicit none
      integer nptk(2) !num of points
      double precision var(2)
      double precision variance
	integer i
      double precision nd,np1,np2,npt
      
      np1=nptk(1)
      np2=nptk(2)
      npt=np1+np2
      nd=n_dim
      
      variance=0.
      do i=1,2
      	variance=variance+var(i)*nptk(i)
	end do
      variance=variance/(npt-2.)
      
      calcBIC2_iso=npt*LogTwoPi+npt*nd*log(variance)+(npt-2.)-2.*np1*log(np1)- &
      2.*np2*log(np2)+(2.+2.*nd+2.*npt)*log(npt)
  end function calcBIC2_iso
      
!----------------------------------------------------------------------

  !calculate BIC for no clustering of the points constituting the given clusters
  double precision function caltBIC(cls1,cls2)
	implicit none
      integer cls1,cls2
      integer n1,n2,n !total num of points
      double precision mean(n_dim),covar(n_dim,n_dim),invcov(n_dim,n_dim),evec(n_dim,n_dim),eval(n_dim),detcov
      integer i,j
      integer ind1,ind2!indices of starting points of the 2 clusters
      double precision, dimension(:,:), allocatable :: ptt
      logical flag
      
      !find total number of points
      n1=ptInClstr(xclsPos(cls1))
      n2=ptInClstr(xclsPos(cls2))
      n=n1+n2
      
      !find the starting indices of points in both the clusters
      ind1=0
      ind2=0
      do i=1,cls1-1
      	ind1=ind1+ptInClstr(xclsPos(i))
	end do
      do i=1,cls2-1
      	ind2=ind2+ptInClstr(xclsPos(i))
	end do
      
      !combine points from both the clusters
      allocate(ptt(n_dim,n))
      ptt(:,1:n1)=p(:,ind1+1:ind1+n1)
      ptt(:,n1+1:n)=p(:,ind2+1:ind2+n2)
      
      !find mean of the merged cluster
      mean(1:n_dim)=(n1*xclsMean(cls1,1:n_dim)+n2*xclsMean(cls2,1:n_dim))/n
      !calculate invcov & det(cov) of the merged cluster
      call calc_covmat(n,n_dim,ptt,mean,covar)
      evec=covar
      flag=.false.
      call diagonalize(evec,eval,n_dim,flag)
      
	!if any of the eval are 0 or -ve then use the lowest +ve eval
	do j=1,n_dim
		if(eval(j)<=0.) eval(1:j)=eval(j+1)
	enddo
            
      detcov=product(eval)
      call calc_invcovmat(n_dim,evec,eval,invcov)
      
      caltBIC=calcBIC1(ptt,n,mean,invcov,detcov)
      
      deallocate(ptt)
      
  end function caltBIC
      
!----------------------------------------------------------------------

  !Anderson Darling test for normal distribution
  logical function AndersonDarling(npt,pt,delMean,alpha)
	implicit none
      integer npt !no. of points
      double precision pt(n_dim,npt) !points
      double precision alpha !significance level
      double precision delMean(n_dim) !difference between the two cluster means
      double precision ppt(npt) !projected & normalized points
      double precision mean,sigma !mean & st.dev. of the projected & normalized points
      double precision A !Anderson Darling statistic
      double precision stn1,stn2
      integer i
      double precision mode
      integer*8 nn,info
      
      mean=0.
      sigma=0.
      !project the points onto delMean & calculate the mean & sigma of the
      !projected points
      mode=sqrt(sum((delMean(1:n_dim)**2)))
      do i=1,npt
      	ppt(i)=sum(pt(1:n_dim,i)*delMean(1:n_dim))/mode
            mean=mean+ppt(i)
            sigma=sigma+ppt(i)**2
      end do
      mean=mean/npt
      sigma=sqrt(sigma/npt-mean**2)*npt/(npt-1.)
      
      !transform the projected points into ones with mean=0 & sigma=1
      ppt(1:npt)=(ppt(1:npt)-mean)/sigma
      
      !sort ppt in ascending order
      !CALL QSORT(ppt,npt,8,COMPARE)
      nn=npt
      call DLASRT('I',nn,ppt,info)
      
      !calculate Anderson Darling statistic
      A=0.
      do i=1,npt
		if(ppt(i)>5.) then
            	stn1=1.
            else if(ppt(i)<-5.) then
            	stn1=0.
		else
            	stn1=stNormalCDF(ppt(i))
		end if
            
		if(ppt(npt+1-i)>5.) then
            	stn2=1.
            else if(ppt(npt+1-i)<-5.) then
            	stn2=0.
		else
            	stn2=stNormalCDF(ppt(npt+1-i))
		end if
            
      	A=A+(2.*i-1.)*(log(stn1)+log(1.-stn2))
      end do
      A=A/(-1.*npt)-npt
      A=A*(1.+4./npt-25./(npt**2))
      
      !inference at alpha=0.0001 
      if(A>1.8692) then
      	AndersonDarling=.false.
      else
      	AndersonDarling=.true.
      end if
      
  end function AndersonDarling
      
!----------------------------------------------------------------------

  INTEGER*2 FUNCTION COMPARE(N1, N2)
	double precision N1
	double precision N2

	IF (N1 .LT. N2) THEN
		COMPARE = -1
	ELSE IF (N1 .EQ. N2) THEN
		COMPARE = 0
	ELSE
		COMPARE = 1
	END IF
      
	RETURN
      
  END FUNCTION COMPARE
      
!----------------------------------------------------------------------

  !learn the Gaussian mixture by expectation maximization
  !some points are constrained to be in specific clusters
  !likelihood constraint
  subroutine GaussMixExpMaxLike(ndim,nCls,npt,ncon,norm,pt,like,wt,wtNorm,locZ,lowlike,switch)
	implicit none
      
      	!input variables
      	integer ndim !dimensionality
      	integer nCls !no. of clusters
      	integer npt !total no. of points
      	integer ncon(nCls) !no. constrained points in each cluster
      	double precision pt(ndim,npt) !points
      	double precision like(2,npt) !log-like & log of the dx
      	logical norm !rescale the points so that all have the same range?
	double precision lowlike(nCls)
      	logical switch !initialize the LAPACK eigenanalysis routines?
      
      	!output variables
      	double precision wt(npt,nCls) !probability weights of points for each cluster
      	double precision wtNorm(npt,nCls) !normalized probability weights of points for each cluster
	double precision locZ(nCls) !local evidence
      
      	!work variables
	integer i,k,count
      	double precision d1
	double precision pts(ndim,npt)
      	double precision mean(nCls,ndim),eval(nCls,ndim),evec(nCls,ndim,ndim),cwt(nCls)
      	double precision old_locZ(nCls),covmat(nCls,ndim,ndim),invcov(nCls,ndim,ndim),detcov(nCls)
      
      	n_dim=ndim
	count=0
      	!rescaling
      	pts=pt
      	if(norm) call rescale(ndim,npt,pts)
      
      	wt=0.d0
      	!set the initial probability weights
      	k=0
      	do i=1,nCls
		wt(k+1:k+ncon(i),i)=1.d0
            
		k=k+ncon(i)
	enddo
      
      	!calculate the initial Gaussian mixture model
	wtNorm=wt
	do i=1,nCls
		!calculate the evidence & normalized weights
	    	call setWt(npt,like,wtNorm(:,i),locZ(i),.true.)
		!now the Gaussian with normalized weights
            	call GaussProp(npt,ndim,pts,wt(:,i),wtNorm(:,i),cwt(i),mean(i,:),covmat(i,:,:), &
			invcov(i,:,:),eval(i,:),evec(i,:,:),detcov(i),switch,.true.)
	enddo
	!normalize cluster prior probabilities
	d1=sum(cwt(1:nCls))
	cwt(1:nCls)=cwt(1:nCls)/d1
      
	k=sum(ncon(1:nCls))
      	!expectation-maximization
      	do
		count=count+1
      		old_locZ=locZ
		!check all the points now
		do i=k+1,npt
            		!calculate the probability of points lying in each cluster
			call normalProbClsGivPt(nCls,ndim,pts(:,i),like(1,i), &
				lowlike(1:nCls),mean(1:nCls,:),invcov(1:nCls,:,:), &
				detcov(1:nCls),cwt(1:nCls),wt(i,1:nCls))
		enddo
            
      		!re-calculate the Gaussian mixture model
		wtNorm=wt
            	do i=1,nCls
			!calculate the evidence & normalized weights
	    		call setWt(npt,like,wtNorm(:,i),locZ(i),.true.)
			!now the Gaussian with normalized weights
            		call GaussProp(npt,ndim,pts,wt(:,i),wtNorm(:,i),cwt(i),mean(i,:),covmat(i,:,:), &
                  	invcov(i,:,:),eval(i,:),evec(i,:,:),detcov(i),.false.,.true.)
		enddo
            	!normalize cluster prior probabilities
            	d1=sum(cwt(1:nCls))
            	cwt(1:nCls)=cwt(1:nCls)/d1
            
            	!check for convergence
            	d1=sqrt(sum((old_locZ(1:nCls)-locZ(1:nCls))**2.))
            	if(d1<0.0001 .or. count==100) then
			!wt=wtNorm
			exit
		endif
	enddo
            	
  end subroutine GaussMixExpMaxLike
      
!----------------------------------------------------------------------
  
  !returns the normalized normal probability for the given point
  double precision function normalProb(d,p,mean,invcov,detcov)
  
  	implicit none
      
      !input variables
      integer d !dimensionality
      double precision p(d) !point
      double precision mean(d) !mean
      double precision invcov(d,d) !inverse covariance matrix
      double precision detcov !determinant of the covariance matrix
      
      !work variables
      integer j,k
      double precision a,tpt(d),pi
      
      pi=4.*atan(1.d0)
      
      !calculate point-mean
      tpt(1:d)=p(1:d)-mean(1:d)
            	
      !a=sum(Transpose(pt-mean).invcov.(pt-mean))
      a=0.
      do j=1,n_dim
            a=a+tpt(j)*tpt(j)*invcov(j,j)
            do k=j+1,n_dim
            	a=a+2.*tpt(j)*tpt(k)*invcov(j,k)
		end do
	end do
	normalProb=exp(-a/2.)/((detcov**0.5)*((2.*pi)**(d/2.)))
  
  end function normalProb
      
!----------------------------------------------------------------------
  
  !returns the normalized normal probability for the points given a cluster
  subroutine setWt(npt,like,wt,locZ,flag)
  
  	implicit none
	
	!input variables
	integer npt !no. of points
	double precision like(2,npt) !log-like & log of dX of points
	logical flag !set the weights
	!input/output variables
	double precision wt(npt) !weights
	double precision locZ !local evidence
	!work variables
	integer j
	
	
	locZ=-huge(1.d0)*epsilon(1.d0) !logZero
	do j=1,npt
		if(wt(j)>0.d0) locZ=logSumExp(locZ,like(1,j)+like(2,j)+log(wt(j)))
	enddo
	
	if(flag) then
		!now calculate the posterior probabilty weights
		do j=1,npt
			if(wt(j)>0.d0) wt(j)=wt(j)*exp(like(1,j)+like(2,j)-locZ)
		enddo
	endif
  
  end subroutine setWt
      
!----------------------------------------------------------------------
  
  !assuming the model is Gaussian Mixuture, return the probability of each cluster given a particular point
  subroutine normalProbClsGivPt(nCls,d,p,like,lowlike,mean,invcov,detcov,cwt,prob)
  
  	implicit none
      
      !input variables
      integer nCls !no. of clusters
      integer d !dimensionality
      double precision p(d) !point
      double precision like !log-like of the point
      double precision lowlike(nCls) !lowest log-like of each cluster
      double precision mean(nCls,d) !mean
      double precision invcov(nCls,d,d) !inverse covariance matrix
      double precision detcov(nCls) !determinant of the covariance matrix
      double precision cwt(nCls) !cluster prior probabilities
      
      !output variables
      double precision prob(nCls)
      
      !work variables
      integer i,j,k
      double precision a(nCls),tpt(nCls,d),pi,d2
      logical flag
      
      pi=4.d0*atan(1.d0)
      
      do i=1,nCls
      	!calculate point-mean
      	tpt(i,1:d)=p(1:d)-mean(i,1:d)
      enddo
            	
      !a=sum(Transpose(pt-mean).invcov.(pt-mean))
      flag=.false.
      a=0.d0
      do i=1,nCls
      	if(lowlike(i)>=like) then
		do j=1,d
			a(i)=a(i)+tpt(i,j)*tpt(i,j)*invcov(i,j,j)
            		do k=j+1,d
            			a(i)=a(i)+2.d0*tpt(i,j)*tpt(i,k)*invcov(i,j,k)
			enddo
		enddo
		prob(i)=log(cwt(i))-log(detcov(i))/2.d0-a(i)/2.d0
		flag=.true.
	else
		prob(i)=-huge(1.d0)*epsilon(1.d0) !logZero
	endif
      enddo
      
      if(.not.flag) then
      	prob(1:nCls)=0.d0
	return
      endif
      
      d2=prob(1)
      do i=2,nCls
      	d2=logSumExp(d2,prob(i))
      enddo
      
      if(prob(i)==-huge(1.d0)*epsilon(1.d0)) then
      	prob(i)=0.d0
      else
      	prob(1:nCls)=exp(prob(1:nCls)-d2)
      endif
  
  end subroutine normalProbClsGivPt
      
!----------------------------------------------------------------------
  
  subroutine GaussProp(n,d,p,wt,wtNorm,cwt,mean,covmat,invcov,eval,evec,detcov,switch,calcMu)
  
	implicit none
      
      	!input variables
      	integer n !no. of points
	integer d !dimensionality
	double precision p(d,n) !points
	double precision wt(n) !probability weights of each point
	double precision wtNorm(n) !normalized (with the evidence) probability weights of each point
	logical switch !initialize the LAPACK eigenanalysis routine?
	logical calcMu !calculate the mean?
      
      	!output variables
      	double precision mean(d) !mean
      	double precision covmat(d,d) !covariance matrix
      	double precision invcov(d,d) !inverse covariance matrix
      	double precision eval(d) !eigen values
      	double precision evec(d,d) !eigen vectors
      	double precision detcov !determinant of the covariance matrix
      	double precision cwt !no. of points in the cluster
      
      	!work variables
      	integer i
  	
      
      	!total no. of points in the cluster
      	cwt=sum(wt(1:n))
      
      	if(calcMu) then
      		mean=0.d0
      		do i=1,n
            		mean(1:d)=mean(1:d)+p(1:d,i)*wtNorm(i)
		enddo
      	endif
      
      	!covariance matrix
	call calc_covmat_wt(n,d,p,wtNorm,mean,covmat)
      	!eigen analysis
      	evec=covmat
      	call diagonalize(evec,eval,d,switch)
    	!eigenvalues of covariance matrix can't be zero
    	do i=1,d
		if(eval(i)<=0.d0) eval(1:i)=eval(i+1)/100.
	enddo
      	
	!determinant of the covariance matrix
	detcov=product(eval)
      
      	!inverse of covariance matrix
      	call calc_invcovmat(d,evec,eval,invcov)
	
  end subroutine GaussProp
      
!----------------------------------------------------------------------
      
  !rescale the points so that all have the same range
  subroutine  rescale(ndim,npt,pt)
	implicit none
      
      !input variables
      integer ndim !dimensionality
      integer npt !total no. of points
      
      !input/output varialbe
      double precision pt(ndim,npt) !points
      
      !work variables
      integer i
      double precision d1,d2
      
      
      do i=1,ndim
      	d2=maxval(pt(i,1:npt))
      	d1=minval(pt(i,1:npt))
	pt(i,1:npt)=(pt(i,1:npt)-d1)/(d2-d1)
      enddo
      
  end subroutine rescale 
   
!----------------------------------------------------------------------
  
  function Dinosaur(points,npt,ndim,nClstr,ptClstr,naux,auxa,min_pt,maxC, &
  	meanx,invcovx,tmatx,evalx,evecx,kfacx,effx,volx,pVol,cVol,neVol,eVolFrac, &
	nIter,dTol,switch,rFlag,cSwitch,nCdim)
	implicit none
      	
	!input variables
	integer npt !total no. of points
	integer ndim !dimensionality
      	double precision points(ndim,npt) !points
	integer naux !dimensionality of auxiliary points
	double precision auxa(naux,npt) !auxiliary points
	integer min_pt !min no. of points per cluster
	integer maxC !max no. of clusters
	double precision pVol !prior volume
	double precision cVol !current vol
	integer neVol
	double precision eVolFrac(2*neVol,2)
	integer nIter
	double precision dTol !volume tolerance
    	logical switch !initialize the LAPACK eigen analysis routines?
	logical cSwitch
	integer nCdim
	logical rFlag !desperate acceptance
	
	!output variables
	integer nClstr !total no. clusters found
	integer ptClstr(maxC) !points in each cluster
	double precision meanx(maxC,ndim) !cluster means
	double precision invcovx(maxC,ndim,ndim) !cluster inv covarianc matrices
	double precision tmatx(maxC,ndim,ndim) !cluster transformation matrices
	double precision evalx(maxC,ndim) !cluster eigenvalues
	double precision evecx(maxC,ndim,ndim) !cluster eigenvectors
	double precision kfacx(maxC) !cluster point enlargement factors
	double precision effx(maxC) !cluster volume enlargement factors
	double precision volx(maxC) !cluster volumes
	logical Dinosaur !successful?
	
	!work variables
	integer i,i2,j,k,n1,n2,n3,j1,k1,m,l,logLloc(1)
	integer, allocatable :: nptx(:)
	logical flag
	double precision, allocatable :: mean(:), covmat(:,:), invcov(:,:), tmat(:,:), evec(:,:), eval(:), p2(:,:), aux2(:,:)
	double precision kfac,eff,detcov,vol
	double precision d1,d2
	logical, allocatable :: updEll(:)
	integer, allocatable :: updPt(:), cluster(:)
	integer nClstrk !total no. clusters found
	integer, allocatable :: ptClstrk(:) !points in each cluster
	double precision, allocatable :: meank(:,:) !cluster means
	double precision, allocatable :: covmatk(:,:,:) !cluster covarianc matrices
	double precision, allocatable :: invcovk(:,:,:) !cluster inv covarianc matrices
	double precision, allocatable :: tmatk(:,:,:) !cluster transformation matrices
	double precision, allocatable :: evalk(:,:) !cluster eigenvalues
	double precision, allocatable :: eveck(:,:,:) !cluster eigenvectors
	double precision, allocatable :: kfack(:) !cluster point enlargement factors
	double precision, allocatable :: effk(:) !cluster volume enlargement factors
	double precision, allocatable :: detcovk(:) !cluster covariance matrix determinants
	double precision, allocatable :: volk(:) !cluster volumes
	double precision, allocatable :: pointsk(:,:), auxk(:,:)
	!merge operation variables
	integer mChk
	double precision, allocatable :: dis(:,:)
	logical, allocatable :: check(:,:)
	logical gflag
	
	
	n_dim=ndim
		
	!sanity checks
	if(min_pt<2) min_pt=2
	if(npt<min_pt) then
		nClstr=1
		ptClstr(1)=npt
		volx(1)=0.d0
		kfacx(1)=0.d0
		effx(1)=1.d0
      		Dinosaur=.true.
            	return
	endif
	
	maxClstr=min(maxC,ceiling(dble(npt)/dble(min_pt)))
	!maxClstr=npt
      	nCls=1
	numClstr=0
      	ptClstrd=0
	
	
	allocate( updEll(maxC), check(maxC,maxC) )
	allocate( mean(ndim), covmat(ndim,ndim), invcov(ndim,ndim), tmat(ndim,ndim), evec(ndim,ndim),eval(ndim), &
	p2(ndim,npt), aux2(naux,npt) )
	allocate( nptx(npt), updPt(maxC), cluster(npt), ptClstrk(maxC) )
	allocate( meank(maxC,ndim), covmatk(maxC,ndim,ndim), invcovk(maxC,ndim,ndim), tmatk(maxC,ndim,ndim), &
	evalk(maxC,ndim), eveck(maxC,ndim,ndim), kfack(maxC), effk(maxC), detcovk(maxC), volk(maxC), pointsk(ndim,npt), &
	auxk(naux,npt), dis(10,2) )
	
	allocate(xclsMean(maxClstr,ndim),xclsInvcov(maxClstr,ndim,ndim),xclsCovmat(maxClstr,ndim,ndim), &
	xclsTMat(maxClstr,ndim,ndim),xclsEval(maxClstr,ndim),xclsEvec(maxClstr,ndim,ndim),xclsKfac(maxClstr), &
	xclsEff(maxClstr),xclsVol(maxClstr),p(ndim,npt),xclsDetcov(maxClstr), &
	aux(naux,npt),ptInClstr(maxClstr))
      	
	check = .false.
	ptInClstr=0
	
	!first calculate the model with 1 cluster
	call CalcEllProp(npt,ndim,points,mean,covmat,invcov,tmat,evec,eval,detcov,kfac,eff,vol,pVol, &
	switch)
	
		
	!now do the clustering
      	call makeDino(points,npt,ndim,naux,auxa,min_pt,mean,covmat,invcov,tmat,evec,eval,detcov,kfac,eff,vol, &
	pVol,cSwitch,nCdim)

	j=0
	do i=1,numClstr
		if(ptInClstr(i)==0) cycle
		cluster(j+1:j+ptInClstr(i))=i
		j=j+ptInClstr(i)
	enddo
	
	flag=.false.
	
	d1=sum(xclsVol(1:numClstr))

	if(postDino(numClstr,p,npt,aux,naux,ndim,cluster,min_pt,ptInClstr(1:numClstr),xclsMean(1:numClstr,:), &
	xclsCovmat(1:numClstr,:,:),xclsInvcov(1:numClstr,:,:),xclsTmat(1:numClstr,:,:),xclsEval(1:numClstr,:), &
	xclsEvec(1:numClstr,:,:),xclsKfac(1:numClstr),xclsEff(1:numClstr),xclsDetcov(1:numClstr), &
	xclsVol(1:numClstr),cVol)) then
		nClstrk=numClstr
		pointsk=p
		auxk=aux
		ptClstrk(1:nClstrk)=ptInClstr(1:nClstrk)
		meank(1:nClstrk,:)=xclsMean(1:nClstrk,:)
		covmatk(1:nClstrk,:,:)=xclsCovmat(1:nClstrk,:,:)
		invcovk(1:nClstrk,:,:)=xclsInvcov(1:nClstrk,:,:)
		tmatk(1:nClstrk,:,:)=xclsTmat(1:nClstrk,:,:)
		evalk(1:nClstrk,:)=xclsEval(1:nClstrk,:)
		eveck(1:nClstrk,:,:)=xclsEvec(1:nClstrk,:,:)
		kfack(1:nClstrk)=xclsKfac(1:nClstrk)
		effk(1:nClstrk)=xclsEff(1:nClstrk)
		detcovk(1:nClstrk)=xclsDetcov(1:nClstrk)
		volk(1:nClstrk)=xclsVol(1:nClstrk)
	
		j=sum(ptClstrk(1:nClstrk-1))
		d1=pVol*ptClstrk(nClstrk)/npt
	
      		nCls=1
		numClstr=0
      		ptClstrd=0
		ptInClstr=0
		
		call makeDino(pointsk(:,j+1:npt),ptClstrk(nClstrk),ndim,naux,auxk(:,j+1:npt),min_pt, &
		meank(nClstrk,:),covmatk(nClstrk,:,:),invcovk(nClstrk,:,:),tmatk(nClstrk,:,:), &
		eveck(nClstrk,:,:),evalk(nClstrk,:),detcovk(nClstrk),kfack(nClstrk),effk(nClstrk), &
		volk(nClstrk),d1,cSwitch,nCdim)
			
		if(numClstr > 1) then
			j=0
			do i=1,numClstr
				cluster(j+1:j+ptInClstr(i))=i
				j=j+ptInClstr(i)
			enddo
			
			d1=pVol*ptClstrk(nClstrk)/npt
		
			flag=postDino(numClstr,p(:,1:ptClstrk(nClstrk)),ptClstrk(nClstrk), &
			aux(:,1:ptClstrk(nClstrk)),naux,ndim,cluster,min_pt,ptInClstr(1:numClstr), &
			xclsMean(1:numClstr,:),xclsCovmat(1:numClstr,:,:),xclsInvcov(1:numClstr,:,:), &
			xclsTmat(1:numClstr,:,:),xclsEval(1:numClstr,:),xclsEvec(1:numClstr,:,:), &
			xclsKfac(1:numClstr),xclsEff(1:numClstr),xclsDetcov(1:numClstr),xclsVol(1:numClstr), &
			d1)
			
			j=sum(ptClstrk(1:nClstrk-1))
			
			pointsk(:,j+1:npt)=p(:,1:ptClstrk(nClstrk))
			auxk(:,j+1:npt)=aux(:,1:ptClstrk(nClstrk))
			ptClstrk(nClstrk:nClstrk+numClstr-1)=ptInClstr(1:numClstr)
			meank(nClstrk:nClstrk+numClstr-1,:)=xclsMean(1:numClstr,:)
			covmatk(nClstrk:nClstrk+numClstr-1,:,:)=xclsCovmat(1:numClstr,:,:)
			invcovk(nClstrk:nClstrk+numClstr-1,:,:)=xclsInvcov(1:numClstr,:,:)
			tmatk(nClstrk:nClstrk+numClstr-1,:,:)=xclsTmat(1:numClstr,:,:)
			evalk(nClstrk:nClstrk+numClstr-1,:)=xclsEval(1:numClstr,:)
			eveck(nClstrk:nClstrk+numClstr-1,:,:)=xclsEvec(1:numClstr,:,:)
			kfack(nClstrk:nClstrk+numClstr-1)=xclsKfac(1:numClstr)
			effk(nClstrk:nClstrk+numClstr-1)=xclsEff(1:numClstr)
			detcovk(nClstrk:nClstrk+numClstr-1)=xclsDetcov(1:numClstr)
			volk(nClstrk:nClstrk+numClstr-1)=xclsVol(1:numClstr)
		else
			flag=.true.
		endif
		
		numClstr=nClstrk+numClstr-1
		p=pointsk
		aux=auxk
		ptInClstr(1:numClstr)=ptClstrk(1:numClstr)
		xclsMean(1:numClstr,:)=meank(1:numClstr,:)
		xclsCovmat(1:numClstr,:,:)=covmatk(1:numClstr,:,:)
		xclsInvcov(1:numClstr,:,:)=invcovk(1:numClstr,:,:)
		xclsTmat(1:numClstr,:,:)=tmatk(1:numClstr,:,:)
		xclsEval(1:numClstr,:)=evalk(1:numClstr,:)
		xclsEvec(1:numClstr,:,:)=eveck(1:numClstr,:,:)
		xclsKfac(1:numClstr)=kfack(1:numClstr)
		xclsEff(1:numClstr)=effk(1:numClstr)
		xclsDetcov(1:numClstr)=detcovk(1:numClstr)
		xclsVol(1:numClstr)=volk(1:numClstr)
		
		logLloc=maxloc(aux(1,1:npt))
		if(logLloc(1)>npt-ptInClstr(numClstr)) then
			d2=0d0
		else
			d2=0.005d0
		endif
		
		if(flag) then
			xclsVol(numClstr)=0.d0
			if(dble(ptInClstr(numClstr))/dble(npt)<d2) flag=.false.
		endif
	endif
	
	if(flag) then
		Dinosaur=.false.
		
		deallocate( updEll, check )
		deallocate( mean, covmat, invcov, tmat, evec, eval, p2, aux2 )
		deallocate( nptx, updPt, cluster, ptClstrk )
		deallocate( meank, covmatk, invcovk, tmatk, evalk, eveck, kfack, effk, detcovk, volk, pointsk, auxk, dis )
		
		deallocate(xclsMean,xclsInvcov,xclsTMat,xclsEval,xclsEvec,xclsKfac,xclsEff,xclsVol,xclsCovmat,xclsDetcov,p,aux,ptInClstr)
		
		!update the vol fractions of the past iterations
		eVolFrac(2:2*neVol,:)=eVolFrac(1:2*neVol-1,:)
		eVolFrac(1,1)=cVol/pVol
		eVolFrac(1,2)=nIter
		
		return
	endif
	
	!sort out the clusters with less than min_pt points
	if(numClstr>1) then
		updPt=0
		updEll=.false.
		nptx(1:numClstr)=ptInClstr(1:numClstr)
		gflag=.false.
		k=0
		do i=1,numClstr
			if(ptInClstr(i)<min_pt .and. xclsVol(i)>0.d0) then
				do j=1,ptInClstr(i)
					d1=1.d99
					do l=1,numClstr
						if(l==i .or. ptInClstr(l)<min_pt .or. xclsVol(i)==0.d0) cycle
						d2=MahaDis(ndim,p(:,k+j),xclsMean(l,:),xclsInvcov(l,:,:),xclsKfac(l)*xclsEff(l))
						d2=d2*xclsVol(l)
						if(d2<d1) then
							d1=d2
							m=l
						endif
					enddo
					!enlarge the ellipsoid
					nptx(m)=nptx(m)+1
					l=sum(ptInClstr(1:m-1))
					d2=pVol*dble(nptx(m))/dble(npt)
					call enlargeEll(ndim,p(:,k+j),xclsMean(m,:),xclsEval(m,:),xclsInvcov(m,:,:),xclsKfac(m),xclsEff(m), &
					xclsVol(m),d2)
					
					updPt(k+j)=m
					updEll(m)=.true.
					gflag=.true.
				enddo
			endif	
			k=k+ptInClstr(i)
		enddo
		
		if(gflag) then
			j=0
			k=0
			n1=0
			do i=1,numClstr
				if(ptInClstr(i)>=min_pt .or. xclsVol(i)==0.d0) then
					n1=n1+1
					p2(:,j+1:j+ptInClstr(i))=p(:,k+1:k+ptInClstr(i))
					aux2(:,j+1:j+ptInClstr(i))=aux(:,k+1:k+ptInClstr(i))
					j=j+ptInClstr(i)
					k=k+ptInClstr(i)
				
					if(updEll(i)) then
						do m=1,npt
							if(updPt(m)==i) then
								p2(:,j+1)=p(:,m)
								aux2(:,j+1)=aux(:,m)
								j=j+1
								ptInClstr(i)=ptInClstr(i)+1
							endif
						enddo
					endif
				else
					k=k+ptInClstr(i)
					ptInClstr(i)=0
					xclsVol(i)=0.d0
				endif
			enddo
			p=p2
			aux=aux2
		else
			n1=numClstr
		endif
		
		if(xclsVol(numClstr)==0.d0) then
			i2=numClstr-1
		else
			i2=numClstr
		endif
		
		n1=0
		do i=1,i2
			if(ptInClstr(i)>0) n1=n1+1
		enddo
		
		j1=0
		n3=0
		do i=1,i2
			if(ptInClstr(i)==0) cycle
			n3=n3+1
				
			do
				mChk=min(10,n1-1)
				dis=huge(1.d0)
				do j=1,i2
					if(ptInClstr(j)==0 .or. j==i) cycle
					
					d1=sum((xclsMean(i,:)-xclsMean(j,:))**2)
					
					k1=0
					do k=mChk,1,-1
						if(d1>dis(k,1)) then
							exit
						else
							k1=k
						endif
					enddo
					if(k1/=0) then
						dis(k1+1:mChk,:)=dis(k1:mChk-1,:)
						dis(k1,1)=d1
						dis(k1,2)=dble(j)
					endif
				enddo
					
				!now do the merge check
				flag=.false.
				p2(:,1:ptInClstr(i))=p(:,j1+1:j1+ptInClstr(i))
				aux2(:,1:ptInClstr(i))=aux(:,j1+1:j1+ptInClstr(i))
				do j=1,mChk
					k=int(dis(j,2)) !ellipsoid to check
						
					!already checked
					if(check(i,k)) cycle
					
					check(i,k)=.true.
					check(k,i)=.true.
					
					k1=sum(ptInClstr(1:k-1))
					n2=ptInClstr(i)+ptInClstr(k)
					p2(:,ptInClstr(i)+1:n2)=p(:,k1+1:k1+ptInClstr(k))
					aux2(:,ptInClstr(i)+1:n2)=aux(:,k1+1:k1+ptInClstr(k))
					d1=pVol*dble(n2)/dble(npt)
						
					call CalcEllProp(n2,ndim,p2(:,1:n2),mean,covmat,invcov,tmat,evec,eval,detcov,kfac,eff,vol, &
					d1,switch)
					
					!success?
					if(vol/1.1<xclsVol(i)+xclsVol(k)) then
						xclsMean(i,:)=mean(:)
						xclsInvCov(i,:,:)=invcov(:,:)
						xclsCovmat(i,:,:)=covmat(:,:)
						xclsTMat(i,:,:)=tMat(:,:)
						xclsEval(i,:)=eval(:)
						xclsEvec(i,:,:)=evec(:,:)
						xclsKfac(i)=kfac
						xclsEff(i)=eff
						xclsVol(i)=vol
						xclsVol(k)=0.d0
						xclsDetcov(i)=detcov
						
						!re-arrangement
						if(k<i) then
							p(:,k1+1:j1-ptInClstr(k))=p(:,k1+ptInClstr(k)+1:j1)
							aux(:,k1+1:j1-ptInClstr(k))=aux(:,k1+ptInClstr(k)+1:j1)
							j1=j1-ptInClstr(k)
						else
							p(:,j1+n2+1:k1+ptInClstr(k))=p(:,j1+ptInClstr(i)+1:k1)
							aux(:,j1+n2+1:k1+ptInClstr(k))=aux(:,j1+ptInClstr(i)+1:k1)
						endif
						p(:,j1+1:j1+n2)=p2(:,1:n2)
						aux(:,j1+1:j1+n2)=aux2(:,1:n2)
							
						ptInClstr(i)=n2
						ptInClstr(k)=0
						n1=n1-1
						check(i,:)=.false.
						check(:,i)=.false.
							
						flag=.true.
					endif

					if(flag) exit
				enddo
					
				if(.not.flag) exit
			enddo
			j1=j1+ptInClstr(i)
		enddo
	endif
		
	!check if this was a successful partition
	flag=.false.
	!if the actual vol is an improvement then success
	if(.not.flag) then
		if(rFlag .and. maxClstr==numClstr) then
			flag=.true.
		else
			!if(sum(xclsVol(1:numClstr))<100.d0 .and. &
			!sum(xclsVol(1:numClstr))/cVol-1d0<dTol) flag=.true.
			if(sum(xclsVol(1:numClstr))/cVol-1d0<dTol) flag=.true.
		endif
	endif
	
	!update the vol fractions of the past iterations
	eVolFrac(2:2*neVol,:)=eVolFrac(1:2*neVol-1,:)
	eVolFrac(1,1)=sum(xclsVol(1:numClstr))/pVol
	eVolFrac(1,2)=nIter
	
	if(.not.flag) then
		Dinosaur=.false.
	else
		Dinosaur=.true.
      
      		j=0
      		do i=1,numClstr
      			if(ptInClstr(i)>0) then
            			j=j+1
      				ptClstr(j)=ptInClstr(i)
				meanx(j,:)=xclsMean(i,:)
				invcovx(j,:,:)=xclsInvcov(i,:,:)
				tmatx(j,:,:)=xclsTMat(i,:,:)
				evalx(j,:)=xclsEval(i,:)
				evecx(j,:,:)=xclsEvec(i,:,:)
				kfacx(j)=xclsKfac(i)
				effx(j)=xclsEff(i)
				volx(j)=xclsVol(i)
            		endif
      		enddo
		
		nClstr=j
      		points(1:ndim,1:npt)=p(1:ndim,1:npt)
	      	auxa(1:naux,1:npt)=aux(1:naux,1:npt)
	endif
      	
	
	deallocate( updEll, check )
	deallocate( mean, covmat, invcov, tmat, evec, eval, p2, aux2 )
	deallocate( nptx, updPt, cluster, ptClstrk )
	deallocate( meank, covmatk, invcovk, tmatk, evalk, eveck, kfack, effk, detcovk, volk, pointsk, auxk, dis )
	
      	deallocate(xclsMean,xclsInvcov,xclsTMat,xclsEval,xclsEvec,xclsKfac,xclsEff, &
	xclsVol,xclsCovmat,xclsDetcov,p,aux,ptInClstr)
      
  end function Dinosaur
      
!---------------------------------------------------------------------- 

  recursive subroutine makeDino(pt,npt,ndim,naux,auxa,min_pt,mean,covmat,invcov, &
  	tmat,evec,eval,detcov,kfac,eff,vol,pVol,cSwitch,nCdim)
    	implicit none
	
	!input variables
	integer npt !no. of points to analyze
	double precision pVol !total prior volume
	integer ndim !dimensionality
      	double precision pt(ndim,npt) !points
	integer naux !dimensionality of auxiliary points
	double precision auxa(naux,npt) !auxiliary points
	integer min_pt !min no. of points per cluster
	double precision mean(ndim) !overall mean
	double precision covmat(ndim,ndim) !overall covariance matrix
	double precision invcov(ndim,ndim) !overall inv covariance matrix
	double precision tmat(ndim,ndim) !overall transformation matrix
	double precision evec(ndim,ndim) !overall eigenvector
	double precision eval(ndim) !overall eigenvalues
	double precision detcov !overall covariance matrix determinant
	double precision kfac !overall cluster point enlargement factor
	double precision eff !overall cluster volume enlargement factor
	double precision vol !overall volume
	logical cSwitch
	integer nCdim
	
	!work variables
	integer nk
	parameter(nk=2)
	integer i,k,d,ip
	integer, allocatable :: cluster(:), nptk(:)
	logical flag
	double precision, allocatable :: meank(:,:), covmatk(:,:,:), invcovk(:,:,:), tmatk(:,:,:)
	double precision, allocatable :: eveck(:,:,:), evalk(:,:), detcovk(:),kfack(:), effk(:),volk(:)
	double precision, allocatable :: ptk(:,:,:), auxk(:,:,:)
	double precision d1
	
	if(npt==0) return
	
	allocate( cluster(npt), nptk(nk) )
	allocate( meank(nk,ndim), covmatk(nk,ndim,ndim), invcovk(nk,ndim,ndim), tmatk(nk,ndim,ndim), eveck(nk,ndim,ndim), &
	evalk(nk,ndim), detcovk(nk), kfack(nk), effk(nk), volk(nk), ptk(nk,ndim,npt), auxk(nk,naux,npt) )
		
	if(abs((vol-pVol)/pVol)<0.01 .or. npt<2*min_pt .or. nCls>=maxClstr) then
		if(abs((vol-pVol)/pVol)<0.01) d=0
		flag=.true.
	else
		flag=.false.
	endif
      
      	if(.not.flag) then
		!assign all points to the same cluster & calculate its properties
		cluster=1
		meank(1,:)=mean(:)
		covmatk(1,:,:)=covmat(:,:)
		invcovk(1,:,:)=invcov(:,:)
		tmatk(1,:,:)=tmat(:,:)
		evalk(1,:)=eval(:)
		eveck(1,:,:)=evec(:,:)
		kfack(1)=kfac
		effk(1)=eff
		detcovk(1)=detcov
		volk(1)=vol
		nptk(1)=npt
		nptk(2:nk)=0
		volk(2:nk)=0.d0
		
      		!breakup the points in 2 or 3 clusters
		k=2
		do
			if(k>nk) then
				flag=.true.
				exit
			endif
			
			nptk(k)=0
			
      			d=Dmeans(k,pt,npt,auxa(naux,:),ndim,cluster,min_pt,nptk(1:k),meank(1:k,:),covmatk(1:k,:,:), &
			invcovk(1:k,:,:),tmatk(1:k,:,:),evalk(1:k,:),eveck(1:k,:,:),kfack(1:k),effk(1:k), &
			detcovk(1:k),volk(1:k),pVol,vol,cSwitch,nCdim)
			
			if(d==2) then
				!couldn't partition at all, exit the loop
				flag=.true.
				exit
			else
				!found a partition
      				!calculate the no. of points in each partition & separate them
      				nptk=0
      				do i=1,npt
           	 			nptk(cluster(i))=nptk(cluster(i))+1
            				ptk(cluster(i),:,nptk(cluster(i)))=pt(:,i)
            				auxk(cluster(i),:,nptk(cluster(i)))=auxa(:,i)
				enddo
				
				!was it a better partition?
				if(d==0 .or. vol>2.*pVol) then
					!yes, then exit
					flag=.false.
					exit
				else
					!no, then partition further
					k=k+1
				endif
			endif
		enddo
	endif
      
      	if(.not.flag) then
		nCls=nCls+k-1
		do i=1,k
			d1=pVol*dble(nptk(i))/dble(npt)
			call makeDino(ptk(i,:,1:nptk(i)),nptk(i),ndim,naux,auxk(i,:,1:nptk(i)),min_pt, &
			meank(i,:),covmatk(i,:,:),invcovk(i,:,:),tmatk(i,:,:),eveck(i,:,:),evalk(i,:), &
			detcovk(i),kfack(i),effk(i),volk(i),d1,cSwitch,nCdim)
		enddo
      	else
      		p(:,ptClstrd+1:ptClstrd+npt)=pt(:,1:npt)
            	aux(:,ptClstrd+1:ptClstrd+npt)=auxa(:,1:npt)
      		ip=numClstr+1
            	ptInClstr(ip)=npt
		xclsMean(ip,:)=mean(:)
		xclsInvcov(ip,:,:)=invcov(:,:)
		xclsCovmat(ip,:,:)=covmat(:,:)
		xclsTMat(ip,:,:)=tmat(:,:)
		xclsEval(ip,:)=eval(:)
		xclsEvec(ip,:,:)=evec(:,:)
		xclsKfac(ip)=kfac
		xclsEff(ip)=eff
		xclsVol(ip)=vol
		xclsDetcov(ip)=detcov
            	ptClstrd=ptClstrd+npt
            	numClstr=numClstr+1
      	endif
	
	deallocate( cluster, nptk )
	deallocate( meank, covmatk, invcovk, tmatk, eveck, evalk, detcovk, kfack, effk, volk, ptk, auxk )
      
  end subroutine makeDino 
      
!---------------------------------------------------------------------- 
  
  logical function kDinosaur(points,npt,ndim,nClstr,ptClstr,naux,auxa,min_pt,maxC, &
  	meanx,invcovx,tmatx,evalx,evecx,kfacx,effx,volx,pVol,cVol,neVol,eVolFrac, &
	nIter,dTol,switch,rFlag)
	implicit none
      	
	!input variables
	integer npt !total no. of points
	integer ndim !dimensionality
      	double precision points(ndim,npt) !points
	integer naux !dimensionality of auxiliary points
	double precision auxa(naux,npt) !auxiliary points
	integer min_pt !min no. of points per cluster
	integer maxC !max no. of clusters
	double precision pVol !prior volume
	double precision cVol !current vol
	integer neVol
	double precision eVolFrac(2*neVol,2)
	integer nIter
	double precision dTol !volume tolerance
    	logical switch !initialize the LAPACK eigen analysis routines?
	logical rFlag !desperate acceptance
	!output variables
	integer nClstr !total no. clusters found
	integer ptClstr(maxC) !points in each cluster
	double precision meanx(maxC,ndim) !cluster means
	double precision invcovx(maxC,ndim,ndim) !cluster inv covarianc matrices
	double precision tmatx(maxC,ndim,ndim) !cluster tranformation matrices
	double precision evalx(maxC,ndim) !cluster eigenvalues
	double precision evecx(maxC,ndim,ndim) !cluster eigenvectors
	double precision kfacx(maxC) !cluster point enlargement factors
	double precision effx(maxC) !cluster volume enlargement factors
	double precision volx(maxC) !cluster volumes
	!work variables
	integer i,j,k,nChkd
	logical flag
	integer cluster(npt),nptk(maxC)
	double precision meank(maxC,ndim),covmatk(maxC,ndim,ndim),invcovk(maxC,ndim,ndim),tmatk(maxC,ndim,ndim)
	double precision evalk(maxC,ndim),eveck(maxC,ndim,ndim),kfack(maxC),effk(maxC),detcovk(maxC),volk(maxC)
	double precision pts(ndim,npt),auxk(naux,npt)
	
	kDinosaur = .false.
	n_dim=ndim
	
	!sanity checks
	if(min_pt<ndim+1) min_pt=ndim+1
	if(npt<2*min_pt .or. maxC==1) then
      		nClstr=1
		ptClstr(1)=npt
		
		call CalcEllProp(npt,ndim,points,meanx(1,:),covmatk(1,:,:),invcovx(1,:,:),tmatx(1,:,:), &
		evecx(1,:,:),evalx(1,:),detcovk(1),kfacx(1),effx(1),volx(1),0d0,switch)
		return
	endif
	
	maxClstr=min(maxC,npt/min_pt+1)
	k=maxClstr
	
	call sKmeans(k,points,npt,ndim,cluster,min_pt)
	
	nChkd=0
	nptk=0
	do i=1,k
		do j=1,npt
			if(cluster(j)==i) then
				nptk(i)=nptk(i)+1
				pts(:,nChkd+nptk(i))=points(:,j)
				auxk(:,nChkd+nptk(i))=auxa(:,j)
			endif
		enddo
		
		call CalcEllProp(nptk(i),ndim,pts(:,nChkd+1:nChkd+nptk(i)),meank(i,:),covmatk(i,:,:), &
		invcovk(i,:,:),tmatk(i,:,:),eveck(i,:,:),evalk(i,:),detcovk(i),kfack(i),effk(i),volk(i), &
		0d0,switch)
		
		nChkd=nChkd+nptk(i)
	enddo
	
	!check if this was a successful partition
	flag=.false.
	!if the actual vol is an improvement then success
	if(.not.flag) then
		if(rFlag .and. maxClstr==k) then
			flag=.true.
		else
			if(sum(volk(1:k))<100.d0 .and. sum(volk(1:k))/cVol-1d0<dTol) flag=.true.
		endif
		!write(*,*)"xxx",sum(volk(1:k)),sum(volk(1:k)),cVol,dTol,flag
	endif
	
	!update the vol fractions of the past iterations
	eVolFrac(2:2*neVol,:)=eVolFrac(1:2*neVol-1,:)
	eVolFrac(1,1)=sum(volk(1:k))/pVol
	eVolFrac(1,2)=nIter
	
	if(.not.flag) then
		kDinosaur=.false.
	else
		kDinosaur=.true.
	
		ptClstr=nptk
		nClstr=k
      		points(1:ndim,1:npt)=pts(1:ndim,1:npt)
		auxa(1:naux,1:npt)=auxk(1:naux,1:npt)
		meanx(1:nClstr,:)=meank(1:nClstr,:)
		invcovx(1:nClstr,:,:)=invcovk(1:nClstr,:,:)
		tmatx(1:nClstr,:,:)=tmatk(1:nClstr,:,:)
		evalx(1:nClstr,:)=evalk(1:nClstr,:)
		evecx(1:nClstr,:,:)=eveck(1:nClstr,:,:)
		kfacx(1:nClstr)=kfack(1:nClstr)
		effx(1:nClstr)=effk(1:nClstr)
		volx(1:nClstr)=volk(1:nClstr)
	endif
      
  end function kDinosaur
  
!  ---------------------------------------------------------------------------------

  !Dinosaur clustering
  function postDino(k,pt,npt,auxa,naux,ndim,cluster,min_pt,nptk,meank,covmatk,invcovk, &
  tmatk,evalk,eveck,kfack,effk,detcovk,volk,cVol)
	implicit none
    
    	!input variables
    	integer k !no. of clusters required
    	integer npt !no. of points
    	integer ndim !dimensionality
    	double precision pt(ndim,npt) !points
	integer naux !no. of aux parameters
	double precision auxa(naux,npt)
    	integer min_pt !min no. of points allowed in a cluster
    	double precision cVol !prior volume
    
    	!input/output variables
    	integer nptk(k) !no. of points in each of k clusters
    	integer cluster(npt) !cluster membership of each point for k clusters
    	double precision meank(k,ndim) !input: means of k clusters
    	double precision covmatk(k,ndim,ndim)
    	double precision invcovk(k,ndim,ndim)
    	double precision tmatk(k,ndim,ndim)
    	double precision evalk(k,ndim)
    	double precision eveck(k,ndim,ndim)
    	double precision kfack(k)
    	double precision effk(k)
    	double precision detcovk(k)
    	double precision volk(k)
    	
    	!output variables
    	logical postDino !F if everything OK, T if something more needs to be done
    
    	!work variables
    	integer i,j,i1,i2,indx(1),count
    	double precision d1
    	double precision, allocatable :: h(:,:), mdis(:,:), ptk(:,:), auxk(:,:)
    	logical doCal(k),flag
    	!ellipsoid properties
    	integer, allocatable :: nptx(:), clusterx(:)
    	double precision, allocatable :: meanx(:,:), covmatx(:,:,:), invcovx(:,:,:), tmatx(:,:,:), evalx(:,:)
    	double precision, allocatable :: evecx(:,:,:), kfacx(:), effx(:), detcovx(:), volx(:), fVal(:)
    	double precision fTol
	
	if(k==1) then
		postDino=.false.
		return
	endif
	
	allocate( fVal(k) )
	
	fTol=20.d0
	
	fVal(1:k)=volk(1:k)/(cVol*nptk(1:k)/npt)

	!calculate the distance measure
	postDino=.false.
    	do i=1,k
		if(fVal(i)>fTol) then
			postDino=.true.
			exit
		endif
	enddo
	
	if(.not.postDino) then
		deallocate(fVal)
		return
	else
		allocate( h(k,npt), mdis(k,npt), ptk(ndim,npt), auxk(naux,npt) )
		allocate( nptx(k), clusterx(npt) )
		allocate( meanx(k,ndim), covmatx(k,ndim,ndim), invcovx(k,ndim,ndim), tmatx(k,ndim,ndim), evalx(k,ndim), &
		evecx(k,ndim,ndim), kfacx(k), effx(k), detcovx(k), volx(k) )
	endif
	
	!initialization
	nptx=nptk
	clusterx=cluster
	meanx=meank
	covmatx=covmatk
	invcovx=invcovk
	tmatx=tmatk
	evalx=evalk
	evecx=eveck
	kfacx=kfack
	effx=effk
	detcovx=detcovk
	volx=volk
	count=0
	doCal=.false.
	
    	do i=1,k	
		if(fVal(i)>fTol) then
			h(i,:)=1.d99
		else
			do j=1,npt
				if(fVal(clusterx(j))<=fTol) cycle
				mdis(i,j)=MahaDis(ndim,pt(:,j),meanx(i,:),invcovx(i,:,:),kfacx(i)*effx(i))
				if(mdis(i,j)>1.d0) then
					h(i,j)=1.d99
				else
					h(i,j)=volx(i)*mdis(i,j)/nptx(i)
					h(i,j)=mdis(i,j)
				endif
			enddo
		endif
    	enddo
    
    	!now assign point j to the ellipsoid i such that h(i,j) is min
    	flag=.false.
    	do j=1,npt
		if(fVal(clusterx(j))<=fTol) cycle
		
		indx=minloc(h(1:k,j))
		if(h(indx(1),j)/=1.d99) then
			doCal(indx(1))=.true.
			flag=.true.
			nptx(clusterx(j))=nptx(clusterx(j))-1
			if(nptx(clusterx(j))==0) then
				doCal(clusterx(j))=.false.
				fVal(clusterx(j))=0.d0
			else
				doCal(clusterx(j))=.true.
			endif
			nptx(indx(1))=nptx(indx(1))+1
			clusterx(j)=indx(1)
		endif
    	enddo
		
	if(flag) then
		!point reassigned, recalculate the ellipsoid properties
		do i=1,k
			if(.not.doCal(i)) cycle
			i1=0
			do j=1,npt
				if(clusterx(j)==i) then
					i1=i1+1
					ptk(:,i1)=pt(:,j)
					auxk(:,i1)=auxa(:,j)
				endif
			enddo
			
			d1=cVol*dble(i1)/dble(npt)
			call CalcEllProp(i1,ndim,ptk(:,1:i1),meanx(i,:),covmatx(i,:,:), &
				invcovx(i,:,:),tmatx(i,:,:),evecx(i,:,:),evalx(i,:),detcovx(i),kfacx(i), &
				effx(i),volx(i),d1,.false.)
					
			fVal(i)=volx(i)/(cVol*nptx(i)/npt)
		enddo
		
		nptk=nptx
    		cluster=clusterx
    		meank=meanx
    		covmatk=covmatx
	    	invcovk=invcovx
    		tmatk=tmatx
    		evalk=evalx
    		eveck=evecx
    		kfack=kfacx
	    	effk=effx
    		detcovk=detcovx
    		volk=volx
	endif
	
	!rearrange
	i1=0
	i2=0
	
	postDino=.false.
	do i=1,k
		if(fVal(i)<=fTol .and. nptk(i)>=min_pt) then
			i2=i2+1
			nptx(i2)=nptk(i)
    			meanx(i2,:)=meank(i,:)
    			covmatx(i2,:,:)=covmatk(i,:,:)
	   		invcovx(i2,:,:)=invcovk(i,:,:)
    			tmatx(i2,:,:)=tmatk(i,:,:)
    			evalx(i2,:)=evalk(i,:)
    			evecx(i2,:,:)=eveck(i,:,:)
    			kfacx(i2)=kfack(i)
	    		effx(i2)=effk(i)
    			detcovx(i2)=detcovk(i)
    			volx(i2)=volk(i)
		else
			postDino=.true.
			cycle
		endif		
			
		do j=1,npt
			if(cluster(j)==i) then
				i1=i1+1
				clusterx(i1)=i
				ptk(:,i1)=pt(:,j)
				auxk(:,i1)=auxa(:,j)
			endif
    		enddo
	enddo
		
	if(postDino) then
		nptx(i2+1)=0
		do i=1,k
			if(fVal(i)>fTol .or. nptk(i)<min_pt) then
				nptx(i2+1)=nptx(i2+1)+nptk(i)
				
				do j=1,npt
					if(cluster(j)==i) then
						i1=i1+1
						clusterx(i1)=i
						ptk(:,i1)=pt(:,j)
						auxk(:,i1)=auxa(:,j)
					endif
				enddo
			endif
		enddo
		k=i2+1
		
		d1=cVol*dble(nptx(k))/dble(npt)
		call CalcEllProp(nptx(k),ndim,ptk(:,npt-nptx(k)+1:npt),meanx(k,:),covmatx(k,:,:), &
			invcovx(k,:,:),tmatx(k,:,:),evecx(k,:,:),evalx(k,:),detcovx(k),kfacx(k), &
			effx(k),volx(k),d1,.false.)
	endif
	
	nptk=nptx
    	cluster=clusterx
    	meank=meanx
    	covmatk=covmatx
	invcovk=invcovx
    	tmatk=tmatx
    	evalk=evalx
    	eveck=evecx
    	kfack=kfacx
	effk=effx
    	detcovk=detcovx
    	volk=volx
	pt=ptk
	auxa=auxk
	
	deallocate( h, mdis, ptk, auxk )
	deallocate( nptx, clusterx )
	deallocate( meanx, covmatx, invcovx, tmatx, evalx, evecx, kfacx, effx, detcovx, volx, fVal )
	
	
  end function postDino  
 
      
!---------------------------------------------------------------------- 

end module xmeans_clstr
