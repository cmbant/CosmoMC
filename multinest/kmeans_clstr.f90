! soft k-means clustering with k=2
! D.J.C. MacKay, Information Theory, Inference & Learning Algorithms, 2003, p.304
! Aug 2006

module kmeans_clstr
  use RandomNS
  use utils1
  implicit none
      

contains
  
!----------------------------------------------------------------------
  
  !soft k-means corresponding to the model of axis aligned Gaussians (see
  !MacKay, 2003)
  subroutine kmeans(k,pt,npt,numdim,cluster,min_pt)
    implicit none
 
    integer k,numdim,npt
    double precision pt(numdim,npt)
    double precision sigsq(k,numdim),means(k,numdim),wt(k),r(npt,k)
    integer cluster(npt)
    double precision sumr,totR(k),temp,covmat(numdim,numdim),mean(numdim)
    integer i,j,indx(1),x,nochg,count,min_pt
    logical clstrd
    double precision urv

    if(k>npt/min_pt+1) k=npt/min_pt+1

    clstrd=.false.
    nochg=0
    count=0
    
    do i=1,numdim
    	mean(i)=sum(pt(i,1:npt))/dble(npt)
    end do
    
    call calc_covmat(npt,numdim,pt,mean,covmat)
    
    do i=1,k
    	wt(i)=ranmarns(0)
      urv=ranmarns(0)
      x=int(dble(npt)*urv)+1
    	means(i,:)=pt(:,x)
    end do

!    means(1,1)=0
!    means(1,2)=0
!    means(2,1)=1.
!    means(2,2)=1.

    temp=huge(1.d0)
    do j=1,numdim
	       if(covmat(j,j)<temp) temp=covmat(j,j)
    end do
    sigsq=temp
    temp=sum(wt)
    wt=wt/temp
    
    do
    	count=count+1
      clstrd=.true.
    	do i=1,npt
		do j=1,k
			temp=1.
			do x=1,numdim
				temp=temp*2.506628275*sqrt(sigsq(j,x))
			end do
			r(i,j)=wt(j)*exp(-sum(((means(j,:)-pt(:,i))**2.)/(2.*sigsq(j,:))))/temp
		end do
		sumr=sum(r(i,:))
		r(i,:)=r(i,:)/sumr
		indx=maxloc(r(i,:))
            if(cluster(i)/=indx(1)) clstrd=.false.
		cluster(i)=indx(1)
	end do
      
	if(clstrd) then
		nochg=nochg+1
		if(nochg==2) exit
	else
		nochg=0
	end if
	
	do i=1,k
		totR(i)=sum(r(:,i))
		
		! update means & sigsq
		do j=1,numdim
			means(i,j)=sum(r(:,i)*pt(j,:))/totR(i)
			sigsq(i,j)=sum(r(:,i)*((pt(j,:)-means(i,j))**2.))/totR(i)
		end do
	end do
		
	!update weights
	temp=sum(totR(:))
	wt(:)=totR(:)/temp	
		
    end do
  
	
  end subroutine kmeans  
 
!----------------------------------------------------------------------
 
  !soft k-means corresponding to the model of spherical Gaussians (see MacKay 2003)
  subroutine kmeans2(k,pt,npt,numdim,cluster,min_pt)
    implicit none
 
    integer k,numdim,npt
    double precision pt(numdim,npt)
    double precision sigsq(k),means(k,numdim),wt(k),r(npt,k)
    integer cluster(npt),old_cluster(npt)
    double precision sumr,totR(k),temp,covmat(numdim,numdim),mean(numdim)
    integer i,j,indx(1),x,nochg,min_pt
    logical clstrd
    double precision urv
    
    if(k>npt/min_pt+1) k=npt/min_pt+1
    clstrd=.false.
    nochg=0
    
    do i=1,numdim
    	mean(i)=sum(pt(i,1:npt))/dble(npt)
    end do
    
    call calc_covmat(npt,numdim,pt,mean,covmat)
    
    do i=1,k
    	wt(i)=ranmarns(0)
      urv=ranmarns(0)
      x=int(dble(npt)*urv)+1
    	means(i,:)=pt(:,x)
    end do
    
    temp=huge(1.d0)
    do j=1,numdim
	       if(covmat(j,j)<temp) temp=covmat(j,j)
    end do 
    sigsq=temp      
    temp=sum(wt)
    wt=wt/temp
    
    old_cluster=0
    
    do
    	do i=1,npt
		do j=1,k
			temp=(2.506628275*sqrt(sigsq(j)))**dble(numdim)
			r(i,j)=wt(j)*exp(-sum((means(j,:)-pt(:,i))**2.)/(2.*sigsq(j)))/temp
		end do
		sumr=sum(r(i,:))
		r(i,:)=r(i,:)/sumr
		indx=maxloc(r(i,:))
		cluster(i)=indx(1)
	end do
	
	clstrd=.true.
	do i=1,npt
		if(old_cluster(i).ne.cluster(i)) then
			clstrd=.false.
			exit
		end if
	end do
	if(clstrd) then
		nochg=nochg+1
		if(nochg==2) exit
	else
		nochg=0
	end if
	old_cluster=cluster
	
	do i=1,k
		totR(i)=sum(r(:,i))
		sigsq(i)=0.
		
		! update means & sigsq
		do j=1,numdim
			means(i,j)=sum(r(:,i)*pt(j,:))/totR(i)
			sigsq(i)=sigsq(i)+sum(r(:,i)*((pt(j,:)-means(i,j))**2.))
		end do
		sigsq(i)=sigsq(i)/(numdim*totR(i))
	end do
		
	!update weights
	temp=sum(totR(:))
	wt(:)=totR(:)/temp	
		
    end do
  
	
  end subroutine kmeans2     
 
!----------------------------------------------------------------------
 
  !simple k-means
  subroutine kmeans3(k,pt,npt,numdim,means,cluster,min_pt)
	implicit none
 
    	integer k,numdim,npt,i1,i2,i3
    	double precision pt(numdim,npt)
    	double precision means(k,numdim),dis(min_pt,2)
    	integer cluster(npt),old_cluster(npt),r(npt,k),totR(k)
    	double precision temp,dist
    	integer i,j,x,nochg,scrap(k),min_pt
    	logical clstrd,flag
    	double precision urv,d1

    	if(k>npt/min_pt+1) k=npt/min_pt+1
	
	if(k==1) then
		do i=1,numdim
			means(1,i)=sum(pt(i,1:npt))/dble(npt)
		enddo
		cluster=1
		return
	endif
	
    	clstrd=.false.
    	nochg=0
    
!    call ranmarns(urv)
!    x=int(dble(npt)*urv)+1
!    means(1,:)=pt(:,x)
!    temp=-1.
!    do i=2,k
!    	lmean(:)=sum(means(1:i-1,:))/dble(i-1)
!    	do j=1,npt
!      	dist=sum((lmean(:)-pt(:,j))**2.)
!            if(dist>temp) then
!			temp=dist
!                  x=j
!            end if
!	end do
!	means(i,:)=pt(:,x)
!    end do

!    lmean(:)=means(1,:)
!    do i=1,k/2+1
!    	do
!      	flag=.false.
!    		call ranmarns(urv)
!    		x=int(dble(npt)*urv)+1
!      	do j=1,i-1
!      		if(x==scrap(j)) then
!            		flag=.true.
!                  	exit
!			end if
!		end do
!            if(.not.flag) exit
!	end do
!      scrap(i)=x
!    	means(i*2-1,1:numdim)=pt(1:numdim,x)
!      if(i*2>k) exit
!    	means(i*2,1:numdim)=2.*lmean(1:numdim)-pt(1:numdim,x)
!    end do

    	!choose random points as starting positions
    	do i=1,k
    		do
      			flag=.false.
            		urv=ranmarns(0)
    			x=int(dble(npt)*urv)+1
      			do j=1,i-1
      				if(x==scrap(j)) then
                  			flag=.true.
                        		exit
				endif
			enddo
            		if(flag) then
            			cycle
            		else
            			scrap(i)=x
                  		exit
			endif
		enddo
      		means(i,:)=pt(:,x)
    	enddo
    
    
    	old_cluster=0
    
    	do
    		do i=1,npt
    			temp=huge(1.d0)
			do j=1,k
				dist=sum((means(j,:)-pt(:,i))**2.)
                  		if(dist<temp) then
                  			temp=dist
                        		x=j
                  		endif
			enddo
            		r(i,:)=0
            		r(i,x)=1
			cluster(i)=x
		enddo
	
		clstrd=.true.
		do i=1,npt
			if(old_cluster(i)/=cluster(i)) then
				clstrd=.false.
				exit
			endif
		enddo
      
		if(clstrd) then
			!check if all the clusters have more than min_pt points
			do i=1,k
				if(totR(i)<min_pt) then
					dis=1.d99
					i1=min_pt-totR(i)
					do j=1,npt
						if(cluster(j)/=i .and. totR(cluster(j))>min_pt) then
							d1=sum((means(i,:)-pt(:,j))**2.)
							i3=0
							do i2=i1,1,-1
								if(d1<dis(i2,1)) then
									i3=i2
								else	
									exit
								endif
							enddo
							if(i3/=0) then
								dis(i3+1:i1,:)=dis(i3:i1-1,:)
								dis(i3,1)=d1
								dis(i3,2)=dble(j)
							endif
						endif
					enddo
					do j=1,i1
						i3=int(dis(j,2))
						i2=cluster(i3)
						cluster(i3)=i
						totR(i)=totR(i)+1
						totR(i2)=totR(i2)-1
						r(i3,:)=0
            					r(i3,i)=1
					enddo
				endif
			enddo
			do i=1,k
				!update means
				do j=1,numdim
					means(i,j)=sum(r(:,i)*pt(j,:))/totR(i)
				enddo
			enddo
			exit
		endif
      
		old_cluster=cluster
	
		do i=1,k
			totR(i)=sum(r(:,i))
            		if(totR(i)==0) cycle
		
			!update means
			do j=1,numdim
				means(i,j)=sum(r(:,i)*pt(j,:))/totR(i)
			enddo
		enddo
    	enddo
  
	
  end subroutine kmeans3
 
!---------------------------------------------------------------------- 

  !Incremental K-means (Pham, Dimov, Nguyen, 2005)
  subroutine kmeans4(k,pt,npt,numdim,cluster,min_pt)
    implicit none
 
    integer k,numdim,npt
    double precision pt(:,:)
    double precision means(k,numdim)
    integer cluster(:)
    integer totR(k)
    integer i,j,x,m,min_pt
    logical clstrd
    double precision urv
    double precision distortion(k)!distortion of each cluster
    integer indx(1)
    
    if(k>npt/min_pt+1) k=npt/min_pt+1
    distortion=0.
    totR=0
    
    !choose a random point to be the cluster center
    cluster=1
    totR(1)=npt
    do j=1,numdim
    	means(1,j)=sum(pt(j,:))/dble(npt)
    end do
    
    if(k==1) return
    
    do m=2,k
    	clstrd=.false.
      
    	!pick the cluster with the max distortion
      indx=maxloc(distortion(1:m-1))
            
      !pick a random point from this cluster as new cluster center
      j=0
      urv=ranmarns(0)
    	x=int(dble(totR(indx(1)))*urv)+1
      do i=1,npt
      	if(cluster(i)==indx(1)) j=j+1
            if(j==x) exit
      end do
      means(cluster(i),:)=(means(cluster(i),:)*dble(totR(cluster(i)))-pt(:,i))/dble(totR(cluster(i))-1)
      means(m,:)=pt(:,i)
      totR(cluster(i))=totR(cluster(i))-1
      totR(m)=1
      cluster(i)=m
            
      do
      	clstrd=.true.
            
            !check whether each point is closest to its own cluster or to the newly formed
            !one & assign accordingly
            do i=1,npt
            	if(sum((pt(1:numdim,i)-means(cluster(i),1:numdim))**2.)> &
                  sum((pt(1:numdim,i)-means(m,1:numdim))**2.))then
                  	!point moved so not finished yet
                        clstrd=.false.
                              
                        !update cluster data
                        totR(cluster(i))=totR(cluster(i))-1
                        cluster(i)=m
                        totR(m)=totR(m)+1
                  end if
            end do
                        
		!update means & distortion of each cluster
            means=0.
            distortion=0.
            do i=1,npt
                  means(cluster(i),1:numdim)=means(cluster(i),1:numdim)+pt(1:numdim,i)
                  distortion(cluster(i))=distortion(cluster(i))+sum(pt(1:numdim,i)**2.)
            end do
            do j=1,numdim
                  means(1:m,j)=means(1:m,j)/dble(totR(1:m))
                  distortion(1:m)=distortion(1:m)-dble(totR(1:m))*(means(1:m,j)**2.)
            end do
                  
            if(clstrd) exit
	end do
    end do
  
	
  end subroutine kmeans4  
 
!----------------------------------------------------------------------

  !2-means with k=2 & clusters placed in their expected positions as the starting points
  subroutine kmeans6(pt,npt,numdim,cluster)
    implicit none
 
    integer numdim,npt
    double precision pt(numdim,npt)
    double precision means(2,numdim)
    integer cluster(npt)
    integer totR(2)
    integer i,j,x
    logical clstrd
    double precision covmat(numdim,numdim),evec(numdim,numdim),eval(numdim)
    
    !calculate the mean of the data
    do j=1,numdim
    	means(1,j)=sum(pt(j,1:npt))/dble(npt)
    end do
    
    !do principal component analysis
    call calc_covmat(npt,numdim,pt,means(1,:),covmat)
    evec=covmat
    call diagonalize(evec,eval,numdim,.false.)
    
    !place the cluster centers in their ex[ected locations
    means(2,:)=means(1,:)+evec(:,numdim)*sqrt(2.*eval(numdim)/3.1416)
    means(1,:)=means(1,:)-evec(:,numdim)*sqrt(2.*eval(numdim)/3.1416)
    
    cluster=1
    totR(1)=npt
    totR(2)=0

    do
      clstrd=.true.
            
	!check whether each point is closest to its own cluster or to the newly formed
	!one & assign accordingly
	do i=1,npt
      	if(cluster(i)==1) then
            	x=2
            else
            	x=1
            end if
            if(sum((pt(1:numdim,i)-means(cluster(i),1:numdim))**2.)> &
            sum((pt(1:numdim,i)-means(x,1:numdim))**2.))then
                  !point moved so not finished yet
                  clstrd=.false.
                              
                  !update cluster data
                  totR(cluster(i))=totR(cluster(i))-1
                  cluster(i)=x
                  totR(x)=totR(x)+1
		end if
	end do
                        
	!update means of each cluster
      means=0.
	do i=1,npt
		means(cluster(i),1:numdim)=means(cluster(i),1:numdim)+pt(1:numdim,i)
	end do
	do j=1,numdim
		means(1:2,j)=means(1:2,j)/dble(totR(1:2))
	end do
                  
	if(clstrd) exit
    end do
  	
  end subroutine kmeans6  
 
!----------------------------------------------------------------------
  !Incremental K-means with Unknown K 
  !(Pham, Dimov, Nguyen, 2005, "Incremental K-means Algorithm")
  !(Pham, Dimov, Nguyen, 2005, "Selection of K in K-means Clustering")
  subroutine kmeans5(k,pt,npt,numdim,cluster,min_pt)
    implicit none
 
    integer k,numdim,npt
    double precision pt(numdim,npt)
    double precision means(2,npt/(numdim+1),numdim),meanst(npt/(numdim+1),numdim)
    double precision mean(npt/(numdim+1),numdim)
    integer cls(2,npt),clst(npt),cluster(npt)
    integer totR(2,npt/(numdim+1)),totRt(npt/(numdim+1))!total points in each cluster
    integer i,j,x,m,min_pt
    logical clstrd
    !distortion of each cluster & total distortion of the recent two iterations
    double precision distortion(2,npt/(numdim+1)),totDistor(2),distor(npt/(numdim+1))
    integer indx(1)
    double precision f(npt/(numdim+1)),fmin,a(npt/(numdim+1))!evaluation function
    double precision S(npt/(numdim+1))!total distortion
    integer count!counter for evaluation function
    
    if(k>npt/min_pt+1) k=npt/min_pt+1
    distortion=0.
    totR=0
    f=0.
    a=0.
    S=0.
    count=0
    
    !choose a random point to be the cluster center
    cluster=1
    cls(2,:)=1
    totR(2,1)=npt
    do j=1,numdim
    	means(2,1,j)=sum(pt(j,:))/dble(npt)
      S(1)=S(1)+sum((pt(j,1:npt)-means(2,1,j))**2.)
    end do
    mean(:,:)=means(2,:,:)
    
    f(1)=1.
    fmin=1.
    k=1
    
    a(2)=1.-3./dble(4*numdim)
    m=1
    do 
    	m=m+1
    	clstrd=.false.
      
    	!pick the cluster with the max distortion
      indx=maxloc(distortion(2,1:m-1))
      !break this cluster
      clst(:)=cls(2,:)
      totRt(:)=totR(2,:)
      meanst(:,:)=means(2,:,:) 
      x=m-1 
      j=indx(1)
      call kmeansDis(x,pt,npt,numdim,clst,j,totRt,distor,meanst)
      cls(1,:)=clst(:)
      totR(1,:)=totRt(:)
      means(1,:,:)=meanst(:,:)
      distortion(1,:)=distor(:)
      totDistor(1)=sum(distor(1:m))
      
      do i=1,m-1
      	if(i==indx(1)) cycle
            if(distortion(2,i)>1.5*(totDistor(2)-totDistor(1))) then
            	clst(:)=cls(2,:)
      		totRt(:)=totR(2,:)
     	 		meanst(:,:)=means(2,:,:)
                  x=m-1
                  call kmeansDis(x,pt,npt,numdim,clst,i,totRt,distor,meanst)
                  if(sum(distor(1:m))<totDistor(1)) then
                  	cls(1,:)=clst(:)
      			totR(1,:)=totRt(:)
      			means(1,:,:)=meanst(:,:)
      			distortion(1,:)=distor(:)
      			totDistor(1)=sum(distor(1:m))
                  end if
            end if
      end do
      
            
      !calculate the total distortion
      S(m)=0.
      do i=1,npt
            S(m)=S(m)+sum((pt(1:numdim,i)-means(1,cls(1,i),1:numdim))**2.)
      end do
      
	!calculate the evaluation function
      f(m)=S(m)/(S(m-1)*a(m))
      a(m+1)=a(m)+(1.-a(m))/6.
      !write(*,*)m,fmin,f(m),f(m-1),a(m)
      
      if(f(m)<fmin.and.f(m)<0.85) then
      	fmin=f(m)
            cluster(:)=cls(1,:)
            mean(:,:)=means(1,:,:)
            k=m
            count=0
      else
      	count=count+1
      end if
      
      !if 5 consecutive f values are greater than 0.85 then exit
      if(count==5) exit
      
      means(2,:,:)=means(1,:,:)
      cls(2,:)=cls(1,:)
      totR(2,:)=totR(1,:)
      distortion(2,:)=distortion(1,:)
      totDistor(2)=totDistor(1)
    end do
    
    !meanst=mean
    !do i=1,k
      !do j=1,numdim
    		!call Rescale_nest(meanst(i,j),j)
      !end do
    !end do
    !write(*,*)"means",meanst(1:k,:)
    !write(*,*)"fmin=",fmin
	
  end subroutine kmeans5  
 
!---------------------------------------------------------------------- 

  !Incremental K-means (Pham, Dimov, Nguyen, 2005) with if npt in any cluster is
  !less than numdim+1 than (numdim+1-npt) closest points borrowed from its closest
  !neighbour cluster 
  subroutine kmeans7(k,pt,npt,numdim,cluster,ad,cluster2,min_p)
    implicit none
 
    integer k,numdim,npt
    double precision pt(:,:)
    double precision means(npt/min_p+10,numdim)
    integer cluster(:)
    integer totR(npt/min_p+10),min_p
    integer i,j,x,m,q,r,k0
    logical clstrd
    double precision urv
    !double precision distortion(k)!distortion of each cluster
    integer indx(1),indx2(1)
    double precision adis(npt)!for keeping a record of the point's distances
    double precision dis1,dis2,dis(npt,2)
    integer ad,nc,cluster2(:,:),cc(npt),cpt(npt)
    
    !distortion=0.
    means=0.
    totR=0
    
    !choose a random point to be the cluster center
    cluster=1
    totR(1)=npt
    do j=1,numdim
    	means(1,j)=sum(pt(j,1:npt))/dble(npt)
    end do
    
    k=min(npt/min_p,k+10)
    if(k<=1) then
    	k=1
      ad=0
      return
    end if
    k0=k

    do m=2,k
    	clstrd=.false.
      
      !pick the cluster with most points
      indx=maxloc(totR(1:m-1))
      
      !if all the clusters have less than npt<2(numdim+1) points then exit
      if(totR(indx(1))<2*min_p) then
      	k=m-1
            exit
      end if
            
      !pick a random point from this cluster as new cluster center
      j=0
      urv=ranmarns(0)
    	x=int(dble(totR(indx(1)))*urv)+1
      do i=1,npt
      	if(cluster(i)==indx(1)) j=j+1
            if(j==x) exit
      end do
      means(cluster(i),:)=(means(cluster(i),:)*dble(totR(cluster(i)))-pt(:,i))/dble(totR(cluster(i))-1)
	means(m,:)=pt(:,i)
      totR(cluster(i))=totR(cluster(i))-1
      totR(m)=1
      cluster(i)=m
      
      do
      	clstrd=.true.
            
            !check whether each point is closest to its own cluster or to the newly formed
            !one & assign accordingly
            do i=1,npt
                  	dis1=sum((pt(1:numdim,i)-means(cluster(i),1:numdim))**2)
                        dis2=sum((pt(1:numdim,i)-means(m,1:numdim))**2)
                        adis(i)=dis1-dis2
            end do
            
            do i=1,numdim
            	means(1:m,i)=means(1:m,i)*dble(totR(1:m))
            end do
            do
            	indx2=maxloc(adis)
                  if(adis(indx2(1))<=0.) then
                  	exit
                  else if(totR(cluster(indx2(1)))<=min_p) then
                  	adis(indx2(1))=-1.
                  else
				!point moved so not finished yet
                        clstrd=.false.
                        !update cluster data
                        totR(cluster(indx2(1)))=totR(cluster(indx2(1)))-1
                       	means(cluster(indx2(1)),1:numdim)=means(cluster(indx2(1)),1:numdim)-pt(1:numdim,indx2(1))
                       	cluster(indx2(1))=m
                        totR(m)=totR(m)+1
                        means(m,1:numdim)=means(m,1:numdim)+pt(1:numdim,indx2(1))
                        adis(indx2(1))=-1.
                  end if
		end do
    
		!update means
            do j=1,numdim
                  means(1:m,j)=means(1:m,j)/dble(totR(1:m))
            end do
                  	
                  
            !if num of points in the newly formed cluster are less than
            !(numdim+1) then borrow closest points from cluster from which the new cluster
            !was formed
            if(clstrd.and.(totR(m)<min_p.or.totR(indx(1))<min_p)) then
            	if(totR(m)<min_p) then
                  	x=m
                        q=indx(1)
                  else
                  	x=indx(1)
                        q=m
                  end if
                  
            	adis=1.d99
            	do i=1,npt
                  	if(cluster(i)==q) then
                        	adis(i)=sum((pt(1:numdim,i)-means(x,1:numdim))**2)
                        end if
                  end do
                  
                  means(q,1:numdim)=means(q,1:numdim)*dble(totR(q))
                  means(x,1:numdim)=means(x,1:numdim)*dble(totR(x))
                  
                  do i=1,min_p-totR(x)
                  	indx2=minloc(adis)
                        adis(indx2(1))=1.d99
                        cluster(indx2(1))=x
                        means(q,1:numdim)=means(q,1:numdim)-pt(1:numdim,indx2(1))
                        means(x,1:numdim)=means(x,1:numdim)+pt(1:numdim,indx2(1))
                  end do
    
                  totR(q)=totR(q)-min_p+totR(x)
                  totR(x)=min_p
                  means(q,1:numdim)=means(q,1:numdim)/dble(totR(q))
                  means(x,1:numdim)=means(x,1:numdim)/dble(totR(x))
            end if
                  
            if(clstrd) exit
	end do
    end do
    
    !find shared points
    ad=min(ad,npt-1)
    nc=min(k-1,ad)
    cluster2=0
    if(ad==0) return
    do i=1,k
      dis(1:nc,1)=1.d99
    	do j=1,k
      	if(i==j) cycle
         	dis1=sum((means(i,1:numdim)-means(j,1:numdim))**2)
            do r=1,nc
            	if(dis1<dis(r,1)) then
                  	dis(r+1:nc,1)=dis(r:nc-1,1)
                  	dis(r+1:nc,2)=dis(r:nc-1,2)
                        dis(r,1)=dis1
                        dis(r,2)=dble(j)
                        exit
                  end if
		end do
	end do
      
      cc(1:nc)=int(dis(1:nc,2))
      
      !make a list of point indices in cluster i
      q=0
      do j=1,npt
      	if(cluster(j)==i) then
                  q=q+1
            	cpt(q)=j
		end if
            if(q==totR(i)) exit
	end do
      
      !adis=1.d99
      dis(1:ad,1)=1.d99
      do j=1,npt
		do r=1,nc
                  if(cluster(j)==cc(r)) then
                  	do q=1,totR(i)
                        	dis1=sum((pt(1:numdim,cpt(q))-pt(1:numdim,j))**2)
                        	do k0=1,ad
                              	if(dis(k0,1)>dis1) then
            					dis(k0+1:ad,1)=dis(k0:ad-1,1)
                  				dis(k0+1:ad,2)=dis(k0:ad-1,2)
                        			dis(k0,1)=dis1
                        			dis(k0,2)=dble(j)
                                          exit
            				end if
					end do
                        end do
                        exit
                  end if
            end do
	end do
      
      cluster2(i,1:ad)=int(dis(1:ad,2))
    end do
	
  end subroutine kmeans7  
 
!----------------------------------------------------------------------

  !Incremental K-means (Pham, Dimov, Nguyen, 2005) with if npt in any cluster is
  !less than numdim+1 than (numdim+1-npt) closest points borrowed from its closest
  !neighbour cluster 
  subroutine sKmeans(k,pt,npt,numdim,cluster,min_p)
	implicit none
 
    	integer k,numdim,npt
    	double precision pt(:,:)
    	double precision means(npt/min_p+10,numdim)
    	integer cluster(npt)
    	integer totR(npt/min_p+10),min_p
    	integer i,j,x,m,q,k0
    	logical clstrd
    	double precision urv
    	integer indx(1),indx2(1)
    	double precision adis(npt)!for keeping a record of the point's distances
    	double precision dis1,dis2
    
	means=0.
    	totR=0
    
    	!choose a random point to be the cluster center
    	cluster=1
    	totR(1)=npt
    	do j=1,numdim
    		means(1,j)=sum(pt(j,1:npt))/dble(npt)
    	end do
    
    	k=min(npt/min_p,k+10)
    	if(k<=1) then
    		k=1
      	return
    	end if
    	k0=k

    	do m=2,k
    		clstrd=.false.
      
      		!pick the cluster with most points
      		indx=maxloc(totR(1:m-1))
      
      		!if all the clusters have less than npt<2(numdim+1) points then exit
      		if(totR(indx(1))<2*min_p) then
      			k=m-1
            		exit
      		endif
            
      		!pick a random point from this cluster as new cluster center
      		j=0
      		urv=ranmarns(0)
    		x=int(dble(totR(indx(1)))*urv)+1
      		do i=1,npt
      			if(cluster(i)==indx(1)) j=j+1
            		if(j==x) exit
      		enddo
      		means(cluster(i),:)=(means(cluster(i),:)*dble(totR(cluster(i)))-pt(:,i))/dble(totR(cluster(i))-1)
		means(m,:)=pt(:,i)
      		totR(cluster(i))=totR(cluster(i))-1
      		totR(m)=1
      		cluster(i)=m
      
      		do
      			clstrd=.true.
            
            		!check whether each point is closest to its own cluster or to the newly formed
            		!one & assign accordingly
            		do i=1,npt
                  		dis1=sum((pt(1:numdim,i)-means(cluster(i),1:numdim))**2)
                        	dis2=sum((pt(1:numdim,i)-means(m,1:numdim))**2)
                        	adis(i)=dis1-dis2
            		enddo
            
			do i=1,numdim
            			means(1:m,i)=means(1:m,i)*dble(totR(1:m))
            		enddo
            		do
            			indx2=maxloc(adis)
                  		if(adis(indx2(1))<=0.) then
                  			exit
                  		elseif(totR(cluster(indx2(1)))<=min_p) then
                  			adis(indx2(1))=-1.
                  		else
					!point moved so not finished yet
                        		clstrd=.false.
                        		!update cluster data
                        		totR(cluster(indx2(1)))=totR(cluster(indx2(1)))-1
                       			means(cluster(indx2(1)),1:numdim)=means(cluster(indx2(1)),1:numdim)-pt(1:numdim,indx2(1))
                       			cluster(indx2(1))=m
                        		totR(m)=totR(m)+1
                        		means(m,1:numdim)=means(m,1:numdim)+pt(1:numdim,indx2(1))
                        		adis(indx2(1))=-1.
                  		endif
			enddo
    
			!update means
            		do j=1,numdim
                  		means(1:m,j)=means(1:m,j)/dble(totR(1:m))
            		enddo
                  	
                  
            		!if num of points in the newly formed cluster are less than
            		!(numdim+1) then borrow closest points from cluster from which the new cluster
            		!was formed
            		if(clstrd.and.(totR(m)<min_p.or.totR(indx(1))<min_p)) then
            			if(totR(m)<min_p) then
                  			x=m
                        		q=indx(1)
                  		else
                  			x=indx(1)
                        		q=m
                  		endif
                  
            			adis=1.d99
            			do i=1,npt
                  			if(cluster(i)==q) then
                        			adis(i)=sum((pt(1:numdim,i)-means(x,1:numdim))**2)
                        		endif
                  		enddo
                  
                  		means(q,1:numdim)=means(q,1:numdim)*dble(totR(q))
                  		means(x,1:numdim)=means(x,1:numdim)*dble(totR(x))
                  
                  		do i=1,min_p-totR(x)
                  			indx2=minloc(adis)
                        		adis(indx2(1))=1.d99
                        		cluster(indx2(1))=x
                        		means(q,1:numdim)=means(q,1:numdim)-pt(1:numdim,indx2(1))
                        		means(x,1:numdim)=means(x,1:numdim)+pt(1:numdim,indx2(1))
                  		enddo
    
                  		totR(q)=totR(q)-min_p+totR(x)
                  		totR(x)=min_p
                  		means(q,1:numdim)=means(q,1:numdim)/dble(totR(q))
                  		means(x,1:numdim)=means(x,1:numdim)/dble(totR(x))
            		endif
                  
            		if(clstrd) exit
		enddo
    	enddo
	
  end subroutine sKmeans
 
!----------------------------------------------------------------------
  
  !soft k-means corresponding to the model of arbitrary Gaussians
  subroutine kmeans8(k,pt,npt,numdim,means,Covmat,cluster,min_pt)
    implicit none
 
    integer k,numdim,npt
    double precision pt(numdim,npt)
    double precision means(k,numdim),wt(k),r(npt,k)
    integer cluster(npt),min_pt
    double precision sumr,totR(k),temp,covmat(k,numdim,numdim),mean(numdim)
    integer i,j,indx(1),x,nochg,count,q
    logical clstrd
    double precision urv
	
    if(k>npt/min_pt+1) k=npt/min_pt+1
    clstrd=.false.
    nochg=0
    count=0
    
    do i=1,numdim
    	mean(i)=sum(pt(i,1:npt))/dble(npt)
    end do
    
    call calc_covmat(npt,numdim,pt,mean,covmat(1,1:numdim,1:numdim))
    do i=2,k
    	covmat(i,:,:)=covmat(1,:,:)
    end do
    
    do i=1,k
      wt(i)=ranmarns(0)
      urv=ranmarns(0)
      x=int(dble(npt)*urv)+1
    	means(i,:)=pt(:,x)
    end do
    temp=sum(wt)
    wt=wt/temp
    
    do
    	count=count+1
      clstrd=.true.
    	do i=1,npt
		do j=1,k
			r(i,j)=wt(j)*mNormalF(numdim,pt(:,i),means(j,:),Covmat(j,:,:))
		end do
		sumr=sum(r(i,:))
		r(i,:)=r(i,:)/sumr
		indx=maxloc(r(i,:))
            if(cluster(i)/=indx(1)) clstrd=.false.
		cluster(i)=indx(1)
	end do
      
	if(clstrd) then
		nochg=nochg+1
		if(nochg==2) exit
	else
		nochg=0
	end if
	
	do i=1,k
		totR(i)=sum(r(:,i))
		
		!update means
		do j=1,numdim
			means(i,j)=sum(r(:,i)*pt(j,:))/totR(i)
		end do
            !update the covariance matrix
            do q=1,numdim
            	do j=q,numdim
            		Covmat(i,q,j)=sum(r(:,i)*pt(q,:)*pt(j,:))/totR(i)-means(i,q)*means(i,j)
			end do
                  if(q/=j) Covmat(i,j,q)=Covmat(i,q,j)
		end do	
	end do
		
	!update weights
	temp=sum(totR(:))
	wt(:)=totR(:)/temp	
    end do
	
  end subroutine kmeans8  
 
!----------------------------------------------------------------------

  !Dinosaur clustering
  function Dmeans(k,pt,npt,like,ndim,cluster,min_pt,nptk,meank,covmatk,invcovk,tmatk,evalk,eveck,kfack,effk, &
  detcovk,volk,pVol,fVol,cSwitch,nCdim)
    implicit none
    
    !input variables
    integer k !no. of clusters required
    integer npt !no. of points
    double precision like(npt) !scaled log-like
    integer ndim !dimensionality
    double precision pt(ndim,npt) !points
    integer min_pt !min no. of points allowed in a cluster
    double precision pVol !prior volume
    logical cSwitch
    integer nCdim
    
    !input/output variables
    integer nptk(k) !input: no. of points in each of k-1 clusters, output: no. of points in each of k clusters
    integer cluster(npt) !input: cluster membership of each point for k-1 clusters, 
    		    !output: cluster membership of each point for k clusters
    double precision meank(k,ndim) !input: means of k-1 clusters, output: means of k clusters
    double precision covmatk(k,ndim,ndim) !input: covmat of k-1 clusters, output: covmat of k clusters
    double precision invcovk(k,ndim,ndim) !input: invcov of k-1 clusters, output: invcov of k clusters
    double precision tmatk(k,ndim,ndim) !input: tmat of k-1 clusters, output: tmat of k clusters
    double precision evalk(k,ndim) !input: eval of k-1 clusters, output: eval of k clusters
    double precision eveck(k,ndim,ndim) !input: evec of k-1 clusters, output: evec of k clusters
    double precision kfack(k) !input: kfac of k-1 clusters, output: kfac of k clusters
    double precision effk(k) !input: eff of k-1 clusters, output: eff of k clusters
    double precision detcovk(k) !input: detcov of k-1 clusters, output: detcov of k clusters
    double precision volk(k) !input: vol of k-1 clusters, output: vol of k clusters
    double precision fVol !volume of the father ellipsoid
    
    !output variables
    integer Dmeans !0 if found a better partition, 1 if found a partition but not better, 2 if couldn't partition
    
    !work variables
    integer i,j,x,i1,i2,indx(1),count,gcount,cls(2)
    double precision, allocatable :: h(:,:),mdis(:,:),ptk(:,:),mu_tmp(:,:)
    double precision d1,d2
    logical flag
    logical, allocatable :: doCal(:), broken(:)
    !ellipsoid properties
    integer, allocatable :: nptx(:), clusterx(:),cluster2(:)
    double precision, allocatable :: meanx(:,:),covmatx(:,:,:),invcovx(:,:,:),tmatx(:,:,:),evalx(:,:)
    double precision, allocatable :: evecx(:,:,:),kfacx(:),effx(:),detcovx(:),volx(:)
    
    
    !sanity check
    if(npt<min_pt*k) then
    	Dmeans=2
	return
    endif
    
    
    allocate( h(k,npt), mdis(k,npt), ptk(ndim+1,npt), mu_tmp(2,ndim+1), meanx(k,ndim), covmatx(k,ndim,ndim), &
    invcovx(k,ndim,ndim), tmatx(k,ndim,ndim), evalx(k,ndim), evecx(k,ndim,ndim), kfacx(k),effx(k), detcovx(k), volx(k) )
    allocate( doCal(k), broken(k) )
    allocate( nptx(k), clusterx(npt), cluster2(npt) )
    
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
    gcount=1
    broken=.false.
    
    do
    	flag=.false.
    	!first find the ellipsoid with most fractional wastage
    	d2=-1.d99
    	do i=1,k-1
    		if(broken(i)) cycle
    		d1=(volx(i)-pVol*nptx(i)/npt)/(pVol*nptx(i)/npt)
		if(d1>d2) then
			flag=.true.
			d2=d1
			j=i
		endif
    	enddo
    
    	!if none of the ellipsoids can be split then return
    	if(.not.flag) exit
	
    	broken(j)=.true.
	
	!split ellipsoid j
	i1=0
	do i=1,npt
		if(clusterx(i)==j) then
			i1=i1+1
			if(cSwitch) then
				ptk(1:nCdim,i1)=pt(1:nCdim,i)
				ptk(nCdim+1,i1)=like(i)
			else
				ptk(1:ndim,i1)=pt(1:ndim,i)
				ptk(ndim+1,i1)=like(i)
			endif	
		endif
	enddo
	
	i2=2
!	if(cSwitch) then
!		call Kmeans3(i2,ptk(1:nCdim+1,1:i1),i1,nCdim+1,mu_tmp(:,1:nCdim+1),cluster2(1:i1),min_pt)
!	else
!		call Kmeans3(i2,ptk(1:ndim+1,1:i1),i1,ndim+1,mu_tmp(:,1:ndim+1),cluster2(1:i1),min_pt)
!	endif
	if(cSwitch) then
		call Kmeans3(i2,ptk(1:nCdim,1:i1),i1,nCdim,mu_tmp(:,1:nCdim),cluster2(1:i1),min_pt)
	else
		call Kmeans3(i2,ptk(1:ndim,1:i1),i1,ndim,mu_tmp(:,1:ndim),cluster2(1:i1),min_pt)
	endif
	
	!separate out the points
	nptx(j)=0
	nptx(k)=0
	i2=1
	do i=1,i1
		do x=i2,npt
			if(clusterx(x)==j) then
				if(cluster2(i)==1) then
					clusterx(x)=j
					nptx(j)=nptx(j)+1
				else
					clusterx(x)=k
					nptx(k)=nptx(k)+1
				endif
				i2=x+1
				exit
			endif
		enddo
	enddo
	
	cls(1)=j
	cls(2)=k
	do i2=1,2
		x=0
		!separate out the points
		do i=1,npt
			if(clusterx(i)==cls(i2)) then
				x=x+1
				ptk(1:ndim,x)=pt(1:ndim,i)
			endif
		enddo
		!min volume this ellipsoid should occupy
		d1=pVol*x/npt
		!calculate the ellipsoid properties
		call CalcEllProp(x,ndim,ptk(1:ndim,1:x),meanx(cls(i2),:),covmatx(cls(i2),:,:), &
		invcovx(cls(i2),:,:),tmatx(cls(i2),:,:),evecx(cls(i2),:,:),evalx(cls(i2),:),detcovx(cls(i2)),kfacx(cls(i2)), &
		effx(cls(i2)),volx(cls(i2)),d1,.false.)
	enddo
	
	!find if it's a better partition
!	if(nptx(j)<min_pt .or. nptx(k)<min_pt) then
!		count=0
!		Dmeans=2
!	else
		count=1
		if(gcount==1) then
			!save the best partition
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
		if(sum(volx(1:k))<fVol) then
			Dmeans=0
		else
			Dmeans=1
		endif
!	endif
    
    	doCal=.true.
		
	!calculate the distance measure
!    	do i=1,k
!		do j=1,npt
!			mdis(i,j)=MahaDis(ndim,pt(:,j),meanx(i,:),invcovx(i,:,:),kfacx(i)*effx(i))
!		enddo
!    	enddo
    
    	do
    		gcount=gcount+1
		!calculate the distance measure
    		do i=1,k
    			if(.not.doCal(i)) cycle
			do j=1,npt
				mdis(i,j)=MahaDis(ndim,pt(:,j),meanx(i,:),invcovx(i,:,:),kfacx(i)*effx(i))
				h(i,j)=volx(i)*mdis(i,j)/nptx(i)
			enddo
    		enddo
    		
!		do x=1,npt
!			i=clusterx(x) !cluster in which point x lies
!			
!			if(doCal(i)) mdis(x,i)=MahaDis(ndim,pt(:,x),meanx(i,:),invcovx(i,:,:),kfacx(i)*effx(i))
!			d2=kfacx(i)*effx(i)
!			d1=(((nptx(i)/(nptx(i)-1.d0))**3)*(1.d0-mdis(x,i)/(nptx(i)-1.d0))-1.d0)*sqrt(detcovx(i)*(d2**(ndim)))/ &
!    				(((nptx(i)/(nptx(i)-1.d0))**1.5d0)*sqrt(1.d0-mdis(x,i)/(nptx(i)-1.d0))+1.d0)
!			
!			do j=1,k
!				if(j==i) then
!					h(j,x)=0.d0
!				else
!					if(doCal(j)) mdis(x,j)=MahaDis(ndim,pt(:,x),meanx(j,:),invcovx(j,:,:),kfacx(j)*effx(j))
!					d2=kfacx(j)*effx(j)
!					h(j,x)=(((nptx(j)/(nptx(j)+1.d0))**3)*(1.d0+mdis(x,j)/(nptx(j)+1.d0))-1.d0)*sqrt(detcovx(j)*(d2**(ndim)))/ &
!    						(((nptx(j)/(nptx(j)+1.d0))**1.5d0)*(sqrt(1.d0+mdis(x,j)/(nptx(j)+1.d0)))+1.d0)
!					h(j,x)=(h(j,x)+d1)!/dble(nptx(i)+nptx(j))
!				endif
!			enddo
!    		enddo
    
    		!now assign point j to the ellipsoid i such that h(i,j) is min
    		flag=.false.
    		do j=1,npt
    			indx=minloc(h(1:k,j))
			if(clusterx(j)/=indx(1)) then
				doCal(clusterx(j))=.true.
				doCal(indx(1))=.true.
				flag=.true.
				nptx(clusterx(j))=nptx(clusterx(j))-1
				nptx(indx(1))=nptx(indx(1))+1
				clusterx(j)=indx(1)
			endif
    		enddo
		
		!if a cluster has less than min_pt points then kill it
		if(flag) then
			do j=1,k
				if(nptx(j)<min_pt) then
					do i1=1,npt
						if(clusterx(i1)==j) then
							h(j,i1)=1.d99
							indx=minloc(h(1:k,i1))
							doCal(indx(1))=.true.
							nptx(j)=nptx(j)-1
							nptx(indx(1))=nptx(indx(1))+1
							clusterx(i1)=indx(1)
						endif
					enddo
					doCal(j)=.false.
					mdis(j,1:npt)=1.d99
					volx(j)=0.d0
					kfacx(j)=0.d0
					evalx(j,:)=0.d0
				endif
			enddo
			exit
		endif

		if(.not.flag) then
			!no point reassigned
			!check if clustering resulted in the new cluster having more than min_pt points
			if(nptx(k)>0 .and. nptx(k)<npt) then
				!check the boundary point of each ellipsoid
				do i=1,k
					if(nptx(i)==min_pt) cycle
			
					!find the boundary point
					d1=0.d0
					do j=1,npt
						if(clusterx(j)==i) then
							if(mdis(i,j)>d1) then
								d1=mdis(i,j)
								i1=j
							endif
						endif
					enddo
			
					!now calculate delF for all the other ellipsoids
					do j=1,k
						if(i==j .or. nptx(j)==0) cycle
						d1=kfacx(i)*effx(i)
						d2=kfacx(j)*effx(j)
						if(delF(ndim,nptx(i),nptx(j),mdis(i,i1),mdis(j,i1),detcovx(i), &
						detcovx(j),d1,d2)<0.) then
							doCal(i)=.true.
							doCal(j)=.true.
							flag=.true.
							nptx(i)=nptx(i)-1
							nptx(j)=nptx(j)+1
							clusterx(i1)=j
							exit
						endif
					enddo
					if(flag) exit
				enddo
			endif
		endif
		
		if(flag) then
			!point reassigned, recalculate the ellipsoid properties
			count=count+1
			!if no better partition found in the past 5 iterations then return
			if(count==20) exit
			if(nptx(k)>0 .and. nptx(k)<npt) then
				do i=1,k
					if(.not.doCal(i)) cycle
					i1=0
					do j=1,npt
						if(clusterx(j)==i) then
							i1=i1+1
							ptk(1:ndim,i1)=pt(1:ndim,j)
						endif
					enddo
					d1=pVol*dble(i1)/dble(npt)
					call CalcEllProp(i1,ndim,ptk(1:ndim,1:i1),meanx(i,:),covmatx(i,:,:), &
						invcovx(i,:,:),tmatx(i,:,:),evecx(i,:,:),evalx(i,:),detcovx(i),kfacx(i), &
						effx(i),volx(i),d1,.false.)
				enddo
				!check if this resulted in a better partition
				if(sum(volx(1:k))<sum(volk(1:k))) then
					!reset the counter
					count=0
					!save the best partition
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
				
					if(sum(volx(1:k))<fVol) then
						Dmeans=0
						if(abs(sum(volx(1:k))-pVol)/pVol<0.01) return
					else
						Dmeans=1
					endif
				endif
			endif
		else
			exit
		endif
	enddo
	if(Dmeans==0) exit
    enddo
    
    deallocate( h, mdis, ptk, mu_tmp, meanx, covmatx, invcovx, tmatx, evalx, evecx, kfacx, effx, detcovx, volx )
    deallocate( doCal, broken )
    deallocate( nptx, clusterx, cluster2 )
    
	
  end function Dmeans  
 
!----------------------------------------------------------------------
  !calculate the variation in weighted average of e-tightness functions of 2 
  !components for re-assigning a point from component 1 to 2
  !Choi, Wang & Kim, Eurographics 2007, vol. 26, No. 3
  function delF(ndim,n1,n2,mdis1,mdis2,detcov1,detcov2,kfac1,kfac2)
    implicit none
    
    !input variables
    integer ndim !dinemsionality
    integer n1,n2 !no. of points in each ellipsoid
    double precision mdis1,mdis2 !Mahalanobis distance of the point from both ellipsoids
    double precision detcov1,detcov2 !determinant of the covariance matrices of the ellipsoids
    double precision kfac1,kfac2 !overall enlargement factors of the ellipsoids
    
    !output variables
    double precision delF
    
    
    delF=(((n1/(n1-1.))**3)*(1.-mdis1/(n1-1.))-1.)*sqrt(detcov1*(kfac1**(ndim)))/ &
    	(((n1/(n1-1.))**1.5)*sqrt(1.-mdis1/(n1-1.))+1.)
    delF=delF+(((n2/(n2+1.))**3)*(1.+mdis2/(n2+1.))-1.)*sqrt(detcov2*(kfac2**(ndim)))/ &
    	(((n2/(n2+1.))**1.5)*(sqrt(1.+mdis2/(n2+1.)))+1.)
	
    delF=delF/(dble(n1+n2))
	
  end function delF  
 
!----------------------------------------------------------------------
  !Incremental K-means with a given cluster 'm' to be split, returns means & distortions 
  !(Pham, Dimov, Nguyen, 2005, "Incremental K-means Algorithm")
  !(Pham, Dimov, Nguyen, 2005, "Selection of K in K-means Clustering")
  subroutine kmeansDis(k,pt,npt,numdim,cluster,m,totR,distortion,means)
    implicit none
 
    integer k,numdim,npt
    double precision pt(numdim,npt)
    double precision means(:,:)
    integer cluster(npt)
    integer totR(:)!total points in each cluster
    integer m!cluster to be split
    integer i,j,x
    logical clstrd
    double precision urv
    double precision distortion(:)!distortion of each cluster
    
    if(k==0) then
    	cluster=1
    	do j=1,numdim
    		means(1,j)=sum(pt(j,:))/dble(npt)
    	end do
     	return
    end if
    
    clstrd=.false.
            
    !pick a random point from the cluster to be split as new cluster center
    j=0
    urv=ranmarns(0)
    x=int(dble(totR(m))*urv)+1
    do i=1,npt
      if(cluster(i)==m) j=j+1
	if(j==x) exit
    end do
    means(k+1,1:numdim)=pt(1:numdim,i)
            
    do
      clstrd=.true.
            
      !check whether each point is closest to its own cluster or to the newly formed
      !one & assign accordingly
      do i=1,npt
		if(sum((pt(1:numdim,i)-means(cluster(i),1:numdim))**2)> &
            sum((pt(1:numdim,i)-means(k+1,1:numdim))**2))then
			!point moved so not finished yet
			clstrd=.false.
                         
			!update cluster data
                  totR(cluster(i))=totR(cluster(i))-1
			cluster(i)=k+1
			totR(k+1)=totR(k+1)+1
		end if
	end do
                        
	!update means & distortion of each cluster
	means=0.
	distortion=0.
	do i=1,npt
		means(cluster(i),1:numdim)=means(cluster(i),1:numdim)+pt(1:numdim,i)
		distortion(cluster(i))=distortion(cluster(i))+sum(pt(1:numdim,i)**2)
	end do
	do j=1,numdim
		means(1:k+1,j)=means(1:k+1,j)/dble(totR(1:k+1))
		distortion(1:k+1)=distortion(1:k+1)-dble(totR(1:k+1))*(means(1:k+1,j)**2)
	end do
                  
	if(clstrd) exit
    end do
	
  end subroutine kmeansDis
 
!----------------------------------------------------------------------  
  
  subroutine cal_means(numdim,npt,p,mean)
    implicit none
    integer numdim,npt
    double precision p(numdim,npt)
    double precision mean(numdim)
    integer i,j

    mean=0.
    do i=1,numdim
       do j=1,npt
          mean(i)=mean(i)+p(i,j)
       end do
       mean(i)=mean(i)/dble(npt)
    end do
    
    return
    
  end subroutine cal_means

!---------------------------------------------------------------------- 

end module kmeans_clstr
