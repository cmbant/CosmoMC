! posterior samples code

module posterior
  use utils1
  
  implicit none
      
  double precision, dimension(:,:,:), allocatable :: evdatp
  integer, dimension(:), allocatable :: nbranchp,nPtPerNode,ncon,nSamp,clstrdNode
  double precision, dimension(:,:), allocatable :: branchp
  double precision, dimension(:,:), allocatable :: pts,pts2,consP,unconsP,pwt,pNwt
  logical, dimension(:), allocatable :: check
  integer nClst,nUncon
  double precision, dimension(:,:), allocatable :: stMu,stSigma

contains 
  
!------------------------------------------------------------------------
      
 subroutine pos_samp(Ztol,nIter,root,nLpt,ndim,nCdim,nPar,multimodal,outfile,globZ,globinfo,ic_n,ic_Z,ic_info,ic_reme, &
 ic_vnow,ic_npt,ic_nBrnch,ic_brnch,phyP,l,evDataAll,dumper,context)
 !subroutine pos_samp
  	implicit none
  	
	double precision Ztol !null evidence
	integer nIter !globff (total no. replacements)
	character(LEN=100)root !base root
	integer nLpt !no. of live points
	integer ndim !dimensionality
	integer nCdim !no. of parameters to cluster on
	integer nPar !total no. of parameters to save
	logical multimodal
	logical outfile !write output files?
	integer i,j,k,i1,ios,ic_n,m,indx
	integer ic_npt(ic_n),ic_nptloc(ic_n),ic_nBrnch(ic_n)
	double precision ic_Z(ic_n),ic_info(ic_n),ic_vnow(ic_n),ic_brnch(ic_n,ic_n),phyP(nPar,nLpt),l(nLpt),evDataAll(:)
	logical ic_reme(ic_n)
  	character(len=100) evfile,livefile,postfile,resumefile
      	character(len=100) sepFile,statsFile,postfile4,strictSepFile,summaryFile
  	character(len=32) fmt,fmt2
      	logical l1
      	double precision d1,d2,urv
      	double precision ltmp(nPar+2),ic_zloc(ic_n),ic_infoloc(ic_n),gzloc,ginfoloc
	integer phyID(nLpt)
	integer context
      
      	!posterior info
      	double precision lognpost,globZ,globInfo,gZold,maxWt
      	integer npost !no. of equally weighted posterior samples
      	double precision, dimension(:,:), allocatable :: wt,temp !probability weight of each posterior sample
      	double precision, dimension(:), allocatable :: ic_zold,llike !local evidence
	
	! parameters for dumper
	double precision, pointer :: physLive(:,:), posterior(:,:), paramConstr(:)
	double precision maxLogLike, logZ, logZerr
	integer nSamples
	
	INTERFACE
		!the user dumper function
    		subroutine dumper(nSamples, nlive, nPar, physLive, posterior, paramConstr, maxLogLike, logZ, logZerr, context_pass)
			integer nSamples, nlive, nPar, context_pass
			double precision, pointer :: physLive(:,:), posterior(:,:), paramConstr(:)
			double precision maxLogLike, logZ, logZerr
		end subroutine dumper
	end INTERFACE
	
	
	ic_zloc=ic_z
	ic_infoloc=ic_info
	gzloc=globz
	ginfoloc=globinfo
	ic_nptloc=ic_npt
	
	
	!Ztol=-1.d99
      	!file names
      	postfile=TRIM(root)//'.txt'
      	statsFile=TRIM(root)//'stats.dat'
      	postfile4=TRIM(root)//'post_equal_weights.dat'
      	strictSepFile=TRIM(root)//'post_separate_strict.dat'
      	sepFile=TRIM(root)//'post_separate.dat'
      	evfile=TRIM(root)//'ev.dat'
      	livefile = TRIM(root)//'phys_live.points'
      	resumefile = TRIM(root)//'resume.dat'
      	summaryFile = TRIM(root)//'summary.txt'
	
      	allocate(branchp(0:ic_n,ic_n),evdatp(ic_n,nIter,nPar+2),wt(ic_n,nIter))
      	allocate(nbranchp(0:ic_n),nPtPerNode(ic_n))
      	allocate(pts2(ndim+1,nIter),pts(nPar+2,nIter),consP(nIter,nPar+2),unconsP(nIter,nPar+2),pwt(nIter,ic_n), &
	pNwt(nIter,ic_n))
      	allocate(ncon(ic_n),nSamp(ic_n),clstrdNode(ic_n))
	allocate(check(ic_n),ic_zold(ic_n),temp(ic_n,3))
	allocate(stMu(ic_n,nPar),stSigma(ic_n,nPar),llike(ic_n))
	
	nPtPerNode=0
	clstrdNode=0
      	nbranchp=0
	check=.false.
	nSamp=0
	
	nbranchp(1:ic_n)=ic_nbrnch(1:ic_n)
	branchp(1:ic_n,:)=ic_brnch(1:ic_n,:)
	
	
	!set the membership of live points to their coresponding nodes
	k=0
	do i=1,ic_n
		if(ic_nptloc(i)==0) cycle
		phyID(k+1:k+ic_nptloc(i))=i
		k=k+ic_nptloc(i)
	enddo
	
	
	!add in the contribution of the remaining live points to the evidence 
	i1=1
	do i=1,nLpt
		ltmp(1:nPar)=phyP(1:nPar,i)
		ltmp(nPar+1)=l(i)
		i1=phyID(i)
	
		d1=ltmp(nPar+1)+log(ic_vnow(i1)/dble(ic_nptloc(i1)))
		
		!global evidence & info
		gZold=gzloc
		gzloc=LogSumExp(gzloc,d1)
!		ginfoloc=exp(d1-gzloc)*ltmp(nPar+1)+exp(gZold-gzloc)*(ginfoloc+gZold)-gzloc
		ginfoloc=ginfoloc*exp(gzold-gzloc)+exp(d1-gzloc)*ltmp(nPar+1)
		
		!local evidence & info
!		gZold=ic_Z(i1)
		ic_zold(i1)=ic_zloc(i1)
		ic_zloc(i1)=LogSumExp(ic_zloc(i1),d1)
!		ic_infoloc(i1)=exp(d1-ic_zloc(i1))*ltmp(nPar+1)+exp(gZold-ic_zloc(i1))*(ic_infoloc(i1)+gZold)-ic_zloc(i1)
		ic_infoloc(i1)=ic_infoloc(i1)*exp(ic_zold(i1)-ic_zloc(i1))+exp(d1-ic_zloc(i1))*ltmp(nPar+1)
	enddo
	
	ginfoloc=ginfoloc-gzloc
	ic_infoloc(1:ic_n)=ic_infoloc(1:ic_n)-ic_zloc(1:ic_n)

      	!make the top level branch
      	nbranchp(0)=1
      	branchp(0,1)=1.d0
        
      	lognpost=0.d0
	if(outfile) then
      		!read the ev.dat file & calculate the probability weights
      		open(unit=55,file=evfile,status='old')
	endif
	i=0
    	do
		if(outfile) then
    			read(55,*,IOSTAT=ios) ltmp(1:nPar+2),i1
		
			!end of file?
        		if(ios<0) then
				close(55)
				exit
			endif
		else
			i=i+1
			if(i*(nPar+3)>size(evDataAll)) exit
			ltmp(1:nPar+2)=evDataAll((i-1)*(nPar+3)+1:i*(nPar+3)-1)
			i1=int(evDataAll(i*(nPar+3)))
		endif
		
		nPtPerNode(i1)=nPtPerNode(i1)+1
		evdatp(i1,nPtPerNode(i1),1:nPar+2)=ltmp(1:nPar+2)
		
		!probability weight  	
		wt(i1,nPtPerNode(i1))=exp(evdatp(i1,nPtPerNode(i1),nPar+1)+evdatp(i1,nPtPerNode(i1),nPar+2)-gzloc)
		if(wt(i1,nPtPerNode(i1))>0.d0) lognpost=lognpost-wt(i1,nPtPerNode(i1))*log(wt(i1,nPtPerNode(i1)))
    	enddo
	
	allocate(posterior(nIter, nPar+2), physLive(nLpt, nPar+1), paramConstr(4*nPar))
	paramConstr(1:2*nPar) = 0d0
	maxLogLike = -huge(1d0)
	logZ = gzloc
	nSamples = nIter
	maxWt = 0d0
	m = 0
      	
 	!read in final remaining points & calculate their probability weights
     	do i=1,nLpt
		ltmp(1:nPar)=phyP(1:nPar,i)
		ltmp(nPar+1)=l(i)
		i1=phyID(i)
		
		physLive(i,1:nPar+1) = ltmp(1:nPar+1)
		if( ltmp(nPar+1) > maxLogLike ) then
			maxLogLike = ltmp(nPar+1)
			indx = i
		endif
		nPtPerNode(i1)=nPtPerNode(i1)+1
            	evdatp(i1,nPtPerNode(i1),1:nPar+1)=ltmp(1:nPar+1)
            	evdatp(i1,nPtPerNode(i1),nPar+2)=log(ic_vnow(i1)/dble(ic_nptloc(i1)))
		wt(i1,nPtPerNode(i1))=exp(evdatp(i1,nPtPerNode(i1),nPar+1)+evdatp(i1,nPtPerNode(i1),nPar+2)-gzloc)
		if(wt(i1,nPtPerNode(i1))>0.d0) lognpost=lognpost-wt(i1,nPtPerNode(i1))*log(wt(i1,nPtPerNode(i1)))
	enddo
	! global maxlike parameters
	paramConstr(nPar*2+1:nPar*3) = physLive(indx,1:nPar)
      	
      	!no. of equally weighted posterior samples
      	npost=nint(exp(lognpost))
	
      	!write the global posterior files
	if(outfile) then
      		open(55,file=postfile,form='formatted',status='replace')
      		open(56,file=postfile4,form='formatted',status='replace')
      		write(fmt,'(a,i4,a)')  '(',nPar+2,'E28.18)'
      		write(fmt2,'(a,i4,a)')  '(',nPar+1,'E28.18)'
	endif
      	do i=1,ic_n
      		do j=1,nPtPerNode(i)
			m = m + 1
			posterior(m, 1:nPar+1) = evdatp(i, j, 1:nPar+1)
			posterior(m, nPar+2) = wt(i,j)
			if( wt(i,j) > maxWt ) then
				indx = m
				maxWt = wt(i,j)
			endif
			! global paramater means
			paramConstr(1:nPar) = paramConstr(1:nPar) + evdatp(i, j, 1:nPar) * wt(i,j)
			! global paramater standard deviations
			paramConstr(nPar+1:2*nPar) = paramConstr(nPar+1:2*nPar) + (evdatp(i, j, 1:nPar)**2.0) * wt(i,j)
			
            		if(wt(i,j)>1.d-99) then
            			if(outfile) write(55,fmt) wt(i,j),-2.d0*evdatp(i,j,nPar+1),evdatp(i,j,1:nPar)
                  	
      				!find the multiplicity
      				d1=wt(i,j)*npost
            			!calculate the integer part of multiplicity
        	    		k=int(d1)
            			!calculate the remaining part of multiplicity
            			d2=d1-dble(k)
            			!increase the multiplicity by one with probability d2
				urv=ranmarns(0)
        	    		if(urv<=d2) k=k+1
				if(outfile) then
      					do i1=1,k
            					write(56,fmt) evdatp(i,j,1:nPar+1)
            				enddo
				endif
			endif
		enddo
	enddo
	
	! global paramater standard deviations
	paramConstr(nPar+1:2*nPar) = sqrt( paramConstr(nPar+1:2*nPar) - paramConstr(1:nPar)**2.0 )
	! global MAP parameters
	paramConstr(nPar*3+1:nPar*4) = posterior(indx,1:nPar)
	
	if(outfile) then
      		close(55)
      		close(56)
      			
		open(unit=57,file=statsFile,form='formatted',status='replace')
      		!stats file
		write(57,'(a,E28.18,a,E28.18)')"Global Evidence:",gzloc,"  +/-",sqrt(ginfoloc/dble(nLpt))
	      		
		!now the separated posterior samples
	      
	      	!generate the point set to be used by the constrained clustering algorithm
	      	nUncon=0
	      	nClst=0
		i=0
		j=1
	      	call genPoints(i,j,nPar,ic_reme)
		
		do i=1,nClst
			temp(i,1)=ic_zloc(clstrdNode(i))
			temp(i,2)=ic_infoloc(clstrdNode(i))
			temp(i,3)=ic_nptloc(clstrdNode(i))
		enddo
		ic_zloc(1:nClst)=temp(1:nClst,1)
		ic_infoloc(1:nClst)=temp(1:nClst,2)
		ic_nptloc(1:nClst)=temp(1:nClst,3)
	      			
		!now arrange the point set, constrained points first, unconstrained later
		pNwt=0.d0
		pwt=0.d0
	      	k=0
	      	do i=1,nClst
			llike(i)=minval(consP(k+1:k+nCon(i),nPar+1))
	      		do j=1,nCon(i)
	      			k=k+1
	          		pts(1:nPar+2,k)=consP(k,1:nPar+2)
				pwt(k,i)=1.d0
			enddo
		enddo
	      	do i=1,nUncon
	      		k=k+1
	      		pts(1:nPar+2,k)=unconsP(i,1:nPar+2)
		enddo
		
		i1=0
	      	do
			if(.false.) then
				i1=i1+1
	      		
				do i=1,k
					pts2(1:ndim,i)=pts(1:ndim,i)
					pts2(ndim+1,i)=pts(nPar+1,i)
				enddo
			
	      			!perform constrained clustering
				i=nCdim
				l1=.true.
	
	      			call GaussMixExpMaxLike(i,nClst,k,nCon(1:nClst),.true.,pts2(1:i,1:k), &
	      			pts(nPar+1:nPar+2,1:k),pwt(1:k,1:nClst),pNwt(1:k,1:nClst),ic_zloc(1:nClst), &
				llike(1:nClst),l1)
			endif
			
			!calculate cluster properties
			do i=1,nClst
				call rGaussProp(k,nCdim,pts(1:nCdim,1:k),pts(nPar+1:nPar+2,1:k), &
				pwt(1:k,i),pNwt(1:k,i),stMu(i,1:nCdim),stSigma(i,1:nCdim),ic_zloc(i))
			enddo
			
			i=sum(nCon(1:nClst))
			j=nCdim
	
			!if(.not.merge(nClst,j,nPar,i,k,nCon(1:nClst),pts(:,1:k),stMu(1:nClst,1:j), &
			!stSigma(1:nClst,1:j),ic_zloc(1:nClst),Ztol,pwt(1:k,1:nClst),pNwt(1:k,1:nClst), &
			!llike(1:nClst))) exit
			exit
		enddo
		
		!open the output file
		if(multimodal .and. outfile) then
			!open(unit=56,file=strictSepFile,form='formatted',status='replace')
	      		open(unit=55,file=sepFile,form='formatted',status='replace')
			write(57,'(a)')
	      		write(57,'(a)')"Local Mode Properties"
	      		write(57,'(a)')"-------------------------------------------"
	      		write(57,'(a)')
			write(57,'(a,i12)')"Total Modes Found:",nClst
		endif
		
		open(unit=58,file=summaryFile,status='unknown')
		call genSepFiles(k,nPar,nClst,Ztol,pts,pNwt(1:k,1:nClst),nCon(1:nClst),ic_zloc(1:nClst), &
		ic_infoloc(1:nClst), ic_nptloc(1:nClst),55,56,57,58,multimodal)
		close(58)
		
	      	if(multimodal) close(55)
	      	close(57)
	endif
	
	!error on global evidence
	logZerr=sqrt(ginfoloc/dble(nLpt))
	
	! call the dumper
	call dumper(nSamples, nLpt, nPar, physLive, posterior, paramConstr, maxLogLike, logZ, logZerr, context)

      	deallocate(branchp,evdatp,wt)
      	deallocate(nbranchp,nPtPerNode)
      	deallocate(pts2,pts,consP,unconsP,pwt,pNwt)
      	deallocate(ncon,nSamp,clstrdNode)
	deallocate(check,ic_zold,temp)
	deallocate(stMu,stSigma,llike)
	deallocate(posterior, physLive, paramConstr)

  end subroutine pos_samp
  
!------------------------------------------------------------------------
      
  recursive subroutine genPoints(br,brNum,nPar,ic_reme)
  
  	implicit none
      
      	integer br !branch to be analyzed
      	integer brNum !branching no. of the branch to be analyzed
	integer nPar !dimensionality
	logical ic_reme(:)
      	!work variables
      	integer i,j,k,i1,i2,node

      	node=int(branchp(br,brNum))
	
      	!find starting node
      	i1=1
      
      	!find out the ending position
      	i2=nPtPerNode(node)
	
      	if(nbranchp(node)==0 .and. .not.ic_reme(node) .and. nPtPerNode(node)>0) then
      		!add the points to the constrained point set if encountered the leaf
		!& calculate the means & sigmas
           	nClst=nClst+1
           	nCon(nClst)=i2-i1+1
            	j=sum(nCon(1:nClst-1))
		stMu(nClst,:)=0.d0
		stSigma(nClst,:)=0.d0
		ClstrdNode(nClst)=node
            	do i=i1,i2
            		j=j+1
            		consP(j,1:nPar+2)=evdatp(node,i,1:nPar+2)
			stMu(nClst,1:nPar)=stMu(nClst,1:nPar)+evdatp(node,i,1:nPar)
			stSigma(nClst,1:nPar)=stSigma(nClst,1:nPar)+evdatp(node,i,1:nPar)**2
            	enddo
		stMu(nClst,1:nPar)=stMu(nClst,1:nPar)/dble(nCon(nClst))
		stSigma(nClst,1:nPar)=stSigma(nClst,1:nPar)/dble(nCon(nClst))
		stSigma(nClst,1:nPar)=sqrt(stSigma(nClst,1:nPar)-stMu(nClst,1:nPar)**2)
	elseif(.not.(nbranchp(node)==0)) then
		if(.not.check(node)) then
            		!add the points to the un-constrained point set
            		do i=i1,i2
            			nUncon=nUncon+1
            			unconsP(nUncon,1:nPar+2)=evdatp(node,i,1:nPar+2)
            		enddo
			check(node)=.true.
		endif
            	!now parse the daughter branches
      		do i=1,nbranchp(node)
            		k=node
            		i1=i
			call genPoints(k,i1,nPar,ic_reme)
		enddo
	endif
      
  end subroutine genPoints
  
!------------------------------------------------------------------------
      
  subroutine genSepFiles(npt,nPar,nCls,Ztol,pt,pwt,nCon,locZ,locInfo,locNpt,funit1,funit2,funit3,funit4,multimodal)
  
  	implicit none
	
	!input variables
	integer npt !total no. of points
	integer nPar !dimensionality
	integer nCls !no. of modes
	double precision pt(nPar+2,npt) !points
	integer nCon(nCls) !no. of constrained points
	double precision pwt(npt,nCls) !ptrobability weights
	double precision locZ(nCls), locInfo(nCls) !local evidence
	integer locNpt(nCls)
	double precision Ztol
	integer funit1 !file having strictly separated samples
	integer funit2 !file having separated samples
	integer funit3 !stats file
	integer funit4 !summary file
	logical multimodal
	!work variables
	integer i,j,k,indx(1),nliveP
	double precision d1,d2,d3,mean(nCls,nPar),sigma(nCls,nPar),maxLike(nPar),MAP(nPar)
	double precision old_slocZ,sinfo,slocZ
	character*30 fmt,stfmt

	nliveP=sum(locNpt(1:nCls))
	write(stfmt,'(a,i4,a)')  '(',nPar*4+2,'E28.18)'

	!calculate the weights including the posterior component
	do i=1,nCls
		k=sum(nCon(1:i-1))
		!first calculate the evidence & info for points strictly lying the cluster
		slocZ=-huge(1.d0)*epsilon(1.d0) !logZero
		sinfo=0.d0
		do j=k+1,k+nCon(i)
			!local evidence
			old_slocZ=slocZ
			slocZ=logSumExp(slocZ,pt(nPar+1,j)+pt(nPar+2,j))
			!local info
			sinfo=exp(pt(nPar+1,j)+pt(nPar+2,j)-slocZ)*pt(nPar+1,j) &
				+exp(old_slocZ-slocZ)*(sinfo+old_slocZ)-slocZ
		enddo
		
		if(locZ(i)<Ztol) cycle
		
		mean(i,:)=0.d0
		sigma(i,:)=0.d0
		nSamp(i)=0
		
		!normalize
		d1=sum(pwt(1:npt,i))
		pwt(1:npt,i)=pwt(1:npt,i)/d1		
		
		!insert two blank line
!		write(funit2,*)
!            	write(funit2,*)
		if(multimodal) then
			write(funit1,*)
            		write(funit1,*)
		endif
		do j=1,npt
			mean(i,1:nPar)=mean(i,1:nPar)+pt(1:nPar,j)*pwt(j,i)
			sigma(i,1:nPar)=sigma(i,1:nPar)+(pt(1:nPar,j)**2)*pwt(j,i)
			
			if(multimodal) then
				!write the strictly separate file
      				write(fmt,'(a,i4,a)')  '(',nPar+2,'E28.18)'
				!strictly separate points
!				if(j>k .and. j<k+nCon(i)+1) then
!					!probability weight
!					swt=exp(pt(nPar+1,j)+pt(nPar+2,j)-slocZ)
!					if(swt>1.d-99) then
!						write(funit2,fmt)swt,-2.d0*pt(nPar+1,j),pt(1:nPar,j)
!					endif
!				endif
			
				!write the separate file
				if(pwt(j,i)>1.d-99) then
					nSamp(i)=nSamp(i)+1
					write(funit1,fmt)pwt(j,i),-2.d0*pt(nPar+1,j),pt(1:nPar,j)
				endif
			endif
		enddo
		sigma(i,1:nPar)=sqrt(sigma(i,1:nPar)-mean(i,1:nPar)**2.)
		
		!stMu(i,:)=mean(i,:)
		!stSigma(i,:)=sigma(i,:)
		
		!find maxLike parameters
		k=sum(nCon(1:i-1))
		indx=maxloc(pt(nPar+1,k+1:k+nCon(i)))
		maxLike(1:nPar)=pt(1:nPar,indx(1)+k)
		d2=pt(nPar+1,indx(1)+k)
		
		!find MAP parameters
		indx=maxloc(pwt(1:npt,i))
		MAP(1:nPar)=pt(1:nPar,indx(1))
		
		!write the stats file
		if(multimodal) then
			write(funit3,*)
			write(funit3,*)
			write(funit3,'(a,i4)')'Mode',i
			d3=(nliveP-locNpt(i))*sinfo/locInfo(i)+locNpt(i)
			write(funit3,'(a,E28.18,a,E28.18)')"Strictly Local Evidence",slocZ," +/-",sqrt(sinfo/locNpt(i))
			write(funit3,'(a,E28.18,a,E28.18)')"Local Evidence",locZ(i)," +/-",sqrt(locInfo(i)/d3)
     		endif
		write(funit3,'(a)')""
		write(funit3,'(a)')"Dim No.       Mean        Sigma"
           	do j=1,nPar
           		!write(funit3,'(i4,2E28.18)')j,stMu(i,j),stSigma(i,j)
           		write(funit3,'(i4,2E28.18)')j,mean(i,j),sigma(i,j)
           	enddo      		
            	write(funit3,'(a)')""
            	write(funit3,'(a)')"Maximum Likelihood Parameters"
            	write(funit3,'(a)')"Dim No.        Parameter"
            	do j=1,nPar
           		write(funit3,'(i4,1E28.18)')j,maxLike(j)
           	enddo
      		write(funit3,'(a)')""
           	write(funit3,'(a)')"MAP Parameters"
           	write(funit3,'(a)')"Dim No.        Parameter"
           	do j=1,nPar
           		write(funit3,'(i4,1E28.18)')j,MAP(j)
           	enddo
		write(funit4,stfmt)mean(i,1:nPar),sigma(i,1:nPar),maxLike(1:nPar),MAP(1:nPar),locZ(i),d2
	enddo
      
  end subroutine genSepFiles
  
!------------------------------------------------------------------------

 logical function merge(n,nCdim,nPar,npt,gnpt,nptx,pt,mean,sigma,locEv,nullEv,wt,nWt,llike)
 
 	implicit none
	
	!input variables
	integer n !no. of clusters
	integer nCdim,nPar !dimensionality
	integer npt !total no. of constrained points
	integer gnpt !total no. of  points
	integer nptx(n) !no. of points in each cluster
	double precision pt(nPar+2,gnpt) !points
	double precision mean(n,nCdim),sigma(n,nCdim)
	double precision locEv(n) !local evidence
	double precision nullEv !null evidence
	double precision wt(gnpt,n)
	double precision nWt(gnpt,n),llike(n)
	
	!work variables
	integer i,j,k,l,m
	double precision tP(nPar+2,npt),tWt(npt,n)
	logical check(n),flag
 	
	
	flag = .false.
	merge=.false.
	check=.false.
	do i=1,n
		if(nptx(i)==0 .or. locEv(i)<nullEv) cycle
		do
			do j=1,n
				if(i==j .or. nptx(j)==0 .or. locEv(j)<nullEv) cycle
				flag=.false.
					
				!merge required?
				l=0
				do k=1,nCdim
					if(abs(mean(i,k)-mean(j,k))<=1.d0*sigma(i,k)) then
						l=l+1
					else
						exit
					endif
				enddo
					
				!yes, then merge the modes
				if(l==nCdim) then
					flag=.true.
					merge=.true.
					
					!re-arrange the constrained points & weights
					wt(:,i)=wt(:,i)+wt(:,j)
					
					!lowest likelihood
					llike(i)=max(llike(i),llike(j))
					
					l=sum(nptx(1:j-1))
					m=sum(nptx(1:i-1))
					tP(:,1:nptx(j))=pt(:,l+1:l+nptx(j))
					tWt(1:nptx(j),:)=wt(l+1:l+nptx(j),:)
					if(i<j) then
						pt(:,m+nptx(i)+nptx(j)+1:l+nptx(j))=pt(:,m+nptx(i)+1:l)
						pt(:,m+nptx(i)+1:m+nptx(i)+nptx(j))=tP(:,1:nptx(j))
						wt(m+nptx(i)+nptx(j)+1:l+nptx(j),:)=wt(m+nptx(i)+1:l,:)
						wt(m+nptx(i)+1:m+nptx(i)+nptx(j),:)=tWt(1:nptx(j),:)
					else
						pt(:,l+1:m-nptx(j))=pt(:,l+nptx(j)+1:m)
						pt(:,m-nptx(j)+1:m)=tP(:,1:nptx(j))
						wt(l+1:m-nptx(j),:)=wt(l+nptx(j)+1:m,:)
						wt(m-nptx(j)+1:m,:)=tWt(1:nptx(j),:)
					endif
					nptx(i)=nptx(i)+nptx(j)
					nptx(j)=0
					
					!recalculate the means & sigmas
					call rGaussProp(gnpt,nCdim,pt(1:nCdim,:),pt(nPar+1:nPar+2,:),wt(:,i), &
					nWt(:,i),mean(i,:),sigma(i,:),locEv(i))
						
					exit
				endif
			enddo
				
			if(.not.flag) exit
		enddo
	enddo
	j=0
	do i=1,n
		if(nptx(i)>0) then
			j=j+1
			nptx(j)=nptx(i)
			locEv(j)=locEv(i)
			mean(j,:)=mean(i,:)
			sigma(j,:)=sigma(i,:)
			wt(:,j)=wt(:,i)
			nWt(:,j)=nWt(:,i)
			llike(j)=llike(i)
		endif
	enddo
	n=j
 
 end function merge
  
!------------------------------------------------------------------------
  
  subroutine rGaussProp(n,d,p,LX,wt,nWt,mean,sigma,Z)
  
	implicit none
      
      	!input variables
      	integer n !no. of points
	integer d !dimensionality
	double precision p(d,n) !points
	double precision LX(2,n) !log-like & log-dX of points
	double precision wt(n)
	double precision nWt(n)
      
      	!output variables
      	double precision mean(d) !mean
      	double precision sigma(d) !standard deviations
      	double precision Z !local evidence
	
      
      	!work variables
      	integer i
  	
      	
	nWt=wt
	!calculate the evidence
!      	Z=-huge(1.d0)*epsilon(1.d0) !logZero
!	do i=1,n
!		if(nWt(i)>0.d0) Z=logSumExp(Z,LX(1,i)+LX(2,i)+log(nWt(i)))
!	enddo
	
	!now calculate the normalized posterior probabilty weights
	do i=1,n
		if(nWt(i)>0.d0) nWt(i)=nWt(i)*exp(LX(1,i)+LX(2,i)-Z)
	enddo
      	
	!mean & sigma
      	mean=0.d0
	sigma=0.d0
      	do i=1,n
            	mean(1:d)=mean(1:d)+p(1:d,i)*nWt(i)
            	sigma(1:d)=sigma(1:d)+(p(1:d,i)**2)*nWt(i)
	enddo
	sigma(1:d)=sqrt(sigma(1:d)-mean(1:d)**2)
	
  end subroutine rGaussProp
      
!----------------------------------------------------------------------


end module posterior
