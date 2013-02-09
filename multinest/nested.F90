! Do nested sampling algorithm to calculate Bayesian evidence
! Jan 2013
! Farhan Feroz

module Nested
  use utils1
  use kmeans_clstr
  use xmeans_clstr
  use posterior
  use priors
  implicit none

#ifdef MPI
  include 'mpif.h'
  integer mpi_status(MPI_STATUS_SIZE), errcode
#endif
  integer my_rank
  integer maxCls,maxeCls
  integer mpi_nthreads !total no. of mpi processors
  integer min_pt,totPar
  integer nCdims !total no. of parameters on which clustering should be done
  integer nlive ! Number of live points
  integer ndims ! Number of dimensions
  integer nsc_def !no. of iterations per every sub-clustering step
  integer updInt !update interval
  double precision Ztol !lowest local evidence for which samples to produce
  double precision tol ! tolerance at end
  double precision ef
  logical multimodal ! multimodal or unimodal sampling
  logical ceff ! constant efficiency?
  integer numlike,globff
  double precision logZero
  integer maxIter
  logical fback,resumeFlag,dlive,genLive,dino
  !output files name
  character(LEN=100)physname,broot,rname,resumename,livename,evname,IS_Files(3)
  !output file units
  integer u_ev,u_resume,u_phys,u_live,u_IS(3)
  double precision gZ,ginfo !total log(evidence) & info
  integer count,sCount
  logical, dimension(:), allocatable :: pWrap
  logical mWrap,aWrap !whether to do wraparound for mode separation
  logical debug, prior_warning, resume, outfile
  !importance sampling
  logical IS

contains
  
  subroutine nestRun(nest_IS,nest_mmodal,nest_ceff,nest_nlive,nest_tol,nest_ef,nest_ndims,nest_totPar,nest_nCdims,maxClst, &
  nest_updInt,nest_Ztol,nest_root,seed,nest_pWrap,nest_fb,nest_resume,nest_outfile,initMPI,nest_logZero,nest_maxIter, &
  loglike,dumper,context)
        
  	implicit none
        
	integer nest_ndims,nest_nlive,nest_updInt,context,seed,i
	integer maxClst,nest_nsc,nest_totPar,nest_nCdims,nest_pWrap(*),nest_maxIter
	logical nest_IS,nest_mmodal,nest_fb,nest_resume,nest_ceff,nest_outfile,initMPI
	character(LEN=100) nest_root
	double precision nest_tol,nest_ef,nest_Ztol,nest_logZero
	
	INTERFACE
    		!the likelihood function
    		subroutine loglike(Cube,n_dim,nPar,lnew,context_pass)
			integer n_dim,nPar,context_pass
			double precision lnew,Cube(nPar)
		end subroutine loglike
    	end INTERFACE
	
	INTERFACE
		!the user dumper function
    		subroutine dumper(nSamples, nlive, nPar, physLive, posterior, paramConstr, maxLogLike, logZ, logZerr, context_pass)
			integer nSamples, nlive, nPar, context_pass
			double precision, pointer :: physLive(:,:), posterior(:,:), paramConstr(:)
			double precision maxLogLike, logZ, logZerr
		end subroutine dumper
	end INTERFACE
	
#ifdef MPI
	if( initMPI ) then
		!MPI initializations
		call MPI_INIT(errcode)
		if (errcode/=MPI_SUCCESS) then
     			write(*,*)'Error starting MPI. Terminating.'
     			call MPI_ABORT(MPI_COMM_WORLD,errcode)
  		end if
	endif
	call MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, errcode)
	call MPI_COMM_SIZE(MPI_COMM_WORLD, mpi_nthreads, errcode)
#else
	mpi_nthreads=1
	my_rank=0
#endif
	IS=nest_IS
	nest_nsc=50
      	nlive=nest_nlive
	Ztol=nest_Ztol
	updInt=nest_updInt
	logZero=nest_logZero
	maxIter=nest_maxIter
	if(maxIter<=0) maxIter=huge(1)
	
	ndims=nest_ndims
      	totPar=nest_totPar
	nCdims=nest_nCdims
	debug=.false.
	prior_warning=.true.
	resume=nest_resume
	outfile=nest_outfile
	if(.not.outfile) resume=.false.
	
      	if(nCdims>ndims) then
		if(my_rank==0) then
			write(*,*)"ERROR: nCdims can not be greater than ndims."
			write(*,*)"Aborting"
		endif
#ifdef MPI
		call MPI_ABORT(MPI_COMM_WORLD,errcode)
#endif
            	stop
	endif
	
	if(my_rank==0) then
		count=0
		sCount=0
      	
      		broot=nest_root
      		rname = trim(broot)
      		fback = nest_fb
		
      		!output file info
      		!setup the output files
      		resumename = trim(rname)//'resume.dat'
      		physname = trim(rname)//'phys_live.points'
      		livename = trim(rname)//'live.points'
      		evname = trim(rname)//'ev.dat'
      		u_ev=55
		u_phys=57
      		u_live=59
		u_resume=61
		
		if( IS ) then
			IS_Files(1) = trim(rname)//'IS.points'
			IS_Files(2) = trim(rname)//'IS.ptprob'
			IS_Files(3) = trim(rname)//'IS.iterinfo'
			u_IS(1) = 62
			u_IS(2) = 63
			u_IS(3) = 64
		endif
	endif
	
	allocate(pWrap(ndims))
	mWrap=.false.
	aWrap=.true.
	do i=1,ndims
		if(nest_pWrap(i)==0) then
			pWrap(i)=.false.
			if(i<=nCdims) aWrap=.false.
		else
			pWrap(i)=.true.
			if(i<=nCdims) mWrap=.true.
		endif
	enddo
      	multimodal=nest_mmodal
	!if( IS ) multimodal = .false.
	ceff=nest_ceff
      	tol=nest_tol
      	ef=nest_ef
	if(ef>=1d6) then
		dino=.false.
		ef=1d0
		min_pt=ndims+1
	elseif(ceff) then
		if(ef>1d0) then
			if(my_rank==0) then
				write(*,*)"ERROR: Can not undersample in constant efficiency mode."
				write(*,*)"Aborting"
			endif
#ifdef MPI
			call MPI_ABORT(MPI_COMM_WORLD,errcode)
#endif
            		stop
		endif
		dino=.true.
		min_pt=2
	else
		dino=.true.
		min_pt=2
	endif
	nsc_def=nest_nsc
	      
      	if(.not.multimodal) then
      		maxCls=1
	else
      		maxCls=maxClst*2
	endif
      
      	maxeCls=nlive/min_pt

      	if(seed<0) then
		!take the seed from system clock
        	call InitRandomNS(mpi_nthreads)
      	else
        	call InitRandomNS(mpi_nthreads,seed)
	endif
	
	if(my_rank==0) then
      		!set the resume flag to true if the resume file exists else to false
      		if(resume) then
			inquire(file=resumename,exist=resumeFlag)
			if(.not.resumeFlag) write(*,*)"MultiNest Warning: no resume file found, starting from scratch"
		else
			resumeFlag=.false.
		endif
      
		write(*,*)"*****************************************************"
		write(*,*)"MultiNest v3.0"
      		write(*,*)"Copyright Farhan Feroz & Mike Hobson"
      		write(*,*)"Release Jan 2013"
		write(*,*)
      		write(*,'(a,i4)')" no. of live points = ",nest_nlive
      		write(*,'(a,i4)')" dimensionality = ",nest_ndims
		if(IS) write(*,'(a)')" using Nested Importance Sampling"
		if(ceff) write(*,'(a)')" running in constant efficiency mode"
      		if(resumeFlag) write(*,'(a)')" resuming from previous job"
      		write(*,*)"*****************************************************"
        
		if (fback) write (*,*) 'Starting MultiNest'
	
	
      		!create the output files
		if(.not.resumeFlag .and. outfile) then
			open(unit=u_ev,file=evname,status='replace')
			close(u_ev)
		endif
	endif
	
	call Nestsample(loglike, dumper, context)
	deallocate(pWrap)
      	call killRandomNS()
#ifdef MPI
	if( initMPI ) call MPI_FINALIZE(errcode)
#endif

  end subroutine nestRun

!----------------------------------------------------------------------

  subroutine Nestsample(loglike, dumper, context)
	
	implicit none
	
	integer context
	double precision, allocatable :: p(:,:), phyP(:,:) !live points
	double precision, allocatable :: l(:) !log-likelihood
	double precision vnow1!current vol
	double precision ltmp(totPar+2)
	character(len=100) fmt
	integer np,i,j,k,ios

	
	INTERFACE
		!the likelihood function
		subroutine loglike(Cube,n_dim,nPar,lnew,context_pass)
			integer n_dim,nPar,context_pass
			double precision lnew,Cube(nPar)
		end subroutine loglike
	end INTERFACE
	
	INTERFACE
		!the user dumper function
    		subroutine dumper(nSamples, nlive, nPar, physLive, posterior, paramConstr, maxLogLike, logZ, logZerr, context_pass)
			integer nSamples, nlive, nPar, context_pass
			double precision, pointer :: physLive(:,:), posterior(:,:), paramConstr(:)
			double precision maxLogLike, logZ, logZerr
		end subroutine dumper
	end INTERFACE
	
	
	allocate( p(ndims,nlive+1), phyP(totPar,nlive+1), l(nlive+1) )

	if(my_rank==0) then
		np=ndims
		globff=0
		numlike=0
		vnow1=1.d0

		write(fmt,'(a,i5,a)')  '(',np+1,'E28.18)'
	
		genLive=.true.
	
		if(resumeflag) then
			!check if the last job was aborted during the live points generation
			open(unit=u_resume,file=resumename,status='old')
			read(u_resume,*)genLive
			
			if( .not.genLive ) then
				read(u_resume,*)i,j,j,j
		
				if( j /= nlive ) then
				  	write(*,*)"ERROR: no. of live points in the resume file is not equal to the the no. passed to nestRun."
					write(*,*)"Aborting"
#ifdef MPI
					call MPI_ABORT(MPI_COMM_WORLD,errcode)
#endif
					stop
				endif
			endif
			close(u_resume)
		
			if( .not.genLive ) then
				j = 0
				open(unit=u_ev,file=evname,status='old') 
				write(fmt,'(a,i2.2,a)')  '(',totPar+2,'E28.18,i3)'
				do
					read(55,*,IOSTAT=ios) ltmp(1:totPar+2),k
				
					!end of file?
					if(ios<0) exit
				
					j = j + 1
				enddo
			
				close(u_ev)
			
				if( j + nlive /= i ) then
					write(*,*)"ERROR: no. of points in ev.dat file is not equal to the no. specified in resume file."
					write(*,*)"Aborting"
#ifdef MPI
					call MPI_ABORT(MPI_COMM_WORLD,errcode)
#endif
					stop
				endif
			endif
		
		endif
	
		gZ=logZero
		ginfo=0.d0
	endif

#ifdef MPI
	call MPI_BARRIER(MPI_COMM_WORLD,errcode)
	call MPI_BCAST(genLive,1,MPI_LOGICAL,0,MPI_COMM_WORLD,errcode)
#endif
	
	if(genLive) then
		if(my_rank==0 .and. fback) write(*,*) 'generating live points'
		
		call gen_initial_live(p,phyP,l,loglike,dumper,context)
	
		if(my_rank==0) then
			globff=nlive
			numlike=nlive
	  		if(fback) write(*,*) 'live points generated, starting sampling'
		endif
	endif

#ifdef MPI
	call MPI_BARRIER(MPI_COMM_WORLD,errcode)
#endif
	
	call clusteredNest(p,phyP,l,loglike,dumper,context)
	
	if(my_rank==0) then
		write(*,*)"ln(ev)=",gZ,"+/-",sqrt(ginfo/dble(nlive))
		write(*,'(a,i12)')' Total Likelihood Evaluations: ', numlike
		write(*,*)"Sampling finished. Exiting MultiNest"
		setBlk=.false.
	endif
	
	deallocate( p, phyP, l )

  end subroutine Nestsample

!----------------------------------------------------------------------

  subroutine gen_initial_live(p,phyP,l,loglike,dumper,context)
    
	implicit none
    
    	integer i,j,iostatus,idum,k,m,nptPerProc,nGen,nstart,nend,context
    	double precision, allocatable :: pnewP(:,:), phyPnewP(:,:), lnewP(:)
    	double precision p(ndims,nlive+1), phyP(totPar,nlive+1), l(nlive+1)
    	integer id
    	character(len=100) fmt,fmt2
#ifdef MPI
	double precision, allocatable ::  tmpl(:), tmpp(:,:), tmpphyP(:,:)
	integer q
#endif
    
    	INTERFACE
    		!the likelihood function
    		subroutine loglike(Cube,n_dim,nPar,lnew,context_pass)
			integer n_dim,nPar,context_pass
			double precision lnew,Cube(nPar)
		end subroutine loglike
    	end INTERFACE
	
	INTERFACE
		!the user dumper function
    		subroutine dumper(nSamples, nlive, nPar, physLive, posterior, paramConstr, maxLogLike, logZ, logZerr, context_pass)
			integer nSamples, nlive, nPar, context_pass
			double precision, pointer :: physLive(:,:), posterior(:,:), paramConstr(:)
			double precision maxLogLike, logZ, logZerr
		end subroutine dumper
	end INTERFACE
	
	allocate( pnewP(ndims,10), phyPnewP(totPar,10), lnewP(10) )
#ifdef MPI
	allocate( tmpl(10), tmpp(ndims,10), tmpphyP(totPar,10) )
#endif

	if(my_rank==0) then
	
		if(outfile) then
    			open(unit=u_resume,file=resumename,form='formatted',status='replace')
    			write(u_resume,'(l2)')genLive
    			close(u_resume)
    			write(fmt,'(a,i5,a)')  '(',ndims+1,'E28.18)'
    			write(fmt2,'(a,i5,a)')  '(',totPar+1,'E28.18,i4)'
		endif

    		id=0
    		i=0
    
    		!resume from previous live points generation?
		if(outfile) then
	    		if(resumeflag) then
	    			!read hypercube-live file
				open(unit=u_live,file=livename,status='old')
				do
	      				i=i+1
					read(u_live,*,IOSTAT=iostatus) p(:,i),l(i)
	            			if(iostatus<0) then
	            				i=i-1
	                  			if(i>nlive) then
							write(*,*)"ERROR: more than ",nlive," points in the live points file."
							write(*,*)"Aborting"
#ifdef MPI
							call MPI_ABORT(MPI_COMM_WORLD,errcode)
#endif
	                        			stop
						endif
	                  			exit
					endif
				enddo
	    			close(u_live)
	
				if(i>0) then
					!read physical live file
					open(unit=u_phys,file=physname,status='old')
					do j=1,i
						read(u_phys,*) phyP(:,j),l(j),idum
					enddo
	    				close(u_phys)
				endif
	    	
	      			open(unit=u_live,file=livename,form='formatted',status='old',position='append')
	    			open(unit=u_phys,file=physname,form='formatted',status='old',position='append')
	    		else
	      			open(unit=u_live,file=livename,form='formatted',status='replace')
	    			open(unit=u_phys,file=physname,form='formatted',status='replace')
				
				if( IS ) then
	      				open(unit=u_IS(2),file=IS_Files(2),form='unformatted',access='sequential',status='replace')
					write(u_IS(2))0,0,0
					close(u_IS(2))
	      				open(unit=u_IS(1),file=IS_Files(1),form='unformatted',access='sequential',status='replace')
					close(u_IS(1))
	      				open(unit=u_IS(3),file=IS_Files(3),form='unformatted',access='sequential',status='replace')
					close(u_IS(3))
				endif
	    		endif
		endif
    
    		j=i
		nend=i
		
		nGen = nlive - j
		nptPerProc = ceiling( dble(nGen) / dble(mpi_nthreads) )
	endif
	
#ifdef MPI
	call MPI_BARRIER(MPI_COMM_WORLD,errcode)
	call MPI_BCAST(nGen,1,MPI_INTEGER,0,MPI_COMM_WORLD,errcode)
	call MPI_BCAST(nptPerProc,1,MPI_INTEGER,0,MPI_COMM_WORLD,errcode)
#endif
      
      	if(nGen==0) then
      		if(my_rank==0) then
            		genLive=.false.
    			resumeFlag=.false.
		endif
		
		deallocate( pnewP, phyPnewP, lnewP )
#ifdef MPI
		deallocate( tmpl, tmpp, tmpphyP )
#endif
		if( outfile ) then
			close(u_live)
			close(u_phys)
		endif
		
            	return
		
	elseif(nGen<0) then
      		if(my_rank==0) then
			write(*,*)"ERROR: live points files have more live points than required."
			write(*,*)"Aborting"
		endif
#ifdef MPI
            	call MPI_ABORT(MPI_COMM_WORLD,errcode)
#endif
            	stop
	endif
	
	k=0
	j=0
	id=my_rank
	do
    		k=k+1
		j=j+1
            	do
			call getrandom(ndims,pnewP(:,j),id)  ! start points
			phyPnewP(1:ndims,j)=pnewP(1:ndims,j)
			lnewP(j)=logZero
			call loglike(phyPnewP(:,j),ndims,totPar,lnewP(j),context)
                  	if(lnewP(j)>logZero) exit
		enddo
		if(k==nptPerProc .or. j==10) then
			if(k==nptPerProc) then
				i=mod(nptPerProc,10)
				if(i==0) i=10
			else
				i=10
				j=0
			endif

#ifdef MPI
			call MPI_BARRIER(MPI_COMM_WORLD,errcode)
#endif
			
			if(id/=0) then
#ifdef MPI
				!send the generated points to the root node
				call MPI_SEND(lnewP(1:i),i,MPI_DOUBLE_PRECISION,0,id,MPI_COMM_WORLD,errcode)
				call MPI_SEND(pnewP(1:ndims,1:i),ndims*i,MPI_DOUBLE_PRECISION,0,id,MPI_COMM_WORLD,errcode)
				call MPI_SEND(phyPnewP(1:totPar,1:i),totPar*i,MPI_DOUBLE_PRECISION,0,id,MPI_COMM_WORLD,errcode)
#endif
			else
				!first write the points generated by the root node
				nstart=nend+1
				p(1:ndims,nstart:nstart+i-1)=pnewP(1:ndims,1:i)
				phyP(1:totPar,nstart:nstart+i-1)=phyPnewP(1:totPar,1:i)
				l(nstart:nstart+i-1)=lnewP(1:i)
				nend=nstart+i-1

#ifdef MPI				
				!receive the points from other nodes
				do m=1,mpi_nthreads-1
					call MPI_RECV(tmpl(1:i),i,MPI_DOUBLE_PRECISION,m,m,MPI_COMM_WORLD,mpi_status,errcode)
					call MPI_RECV(tmpp(1:ndims,1:i),i*ndims,MPI_DOUBLE_PRECISION,m,m,MPI_COMM_WORLD,mpi_status,errcode)
					call MPI_RECV(tmpphyP(1:totPar,1:i),i*totPar,MPI_DOUBLE_PRECISION,m,m,MPI_COMM_WORLD,mpi_status,errcode)
					do q = 1 , i
						if( nend + 1 <= nlive ) then
							l(nend + 1) = tmpl(q)
							p(1 : ndims, nend + 1) = tmpp(1 : ndims, q)
							phyP(1 : totPar, nend + 1) = tmpphyP(1 : totPar, q)
							nend = nend + 1
						endif
					enddo
				enddo
#endif
				
				if( outfile ) then
					!now write this batch to the files
					do m=nstart,nend
						write(u_live,fmt) p(1:ndims,m),l(m)
            					write(u_phys,fmt2) phyP(:,m),l(m),1
					enddo
				endif
			endif
		endif
		if(k==nptPerProc) exit
	enddo
	
	
		
	deallocate( pnewP, phyPnewP, lnewP )
#ifdef MPI
	deallocate( tmpl, tmpp, tmpphyP )
#endif
	
	if( outfile ) then
    		close(u_live)
   		close(u_phys)
	endif
	genLive=.false.
    	resumeFlag=.false.
#ifdef MPI
	call MPI_BARRIER(MPI_COMM_WORLD,errcode)
#endif
    
  end subroutine gen_initial_live

!----------------------------------------------------------------------
  
  subroutine getrandom(n,x,id)
    
    implicit none
    
    integer, intent(in) :: n
    double precision, intent(out) :: x(:)
    integer i,id
    
    ! --- uniform prior ----
    do i = 1,n
       x(i)=ranmarNS(id)
    enddo
    return
  end subroutine getrandom

!----------------------------------------------------------------------
  
  subroutine clusteredNest(p,phyP,l,loglike,dumper,context)
  	
	implicit none
	
	
	!input variables
	
	integer context
	double precision p(ndims,nlive+1) !live points
	double precision phyP(totPar,nlive+1) !physical live points
	double precision l(nlive+1) !log-likelihood
	
	
	!work variables
	
	!misc
	integer i, j, k, m, n, j1, i1, i2, i3, i4, ff, sff, n1, n2, q, nd, nd_i, nd_j, iostatus
	integer num_old
	integer, allocatable :: eswitchff(:), escount(:), dmin(:)
	double precision d1, d2, d3, d4, d5, urv
	double precision h, logX, vprev, vnext, shrink !prior volume
	double precision mar_r !marginal acceptance rate
	double precision gZOld !global evidence & info
	logical eswitch,peswitch,cSwitch !whether to do ellipsoidal sampling or not
	logical remFlag, acpt, flag, flag2
	integer funit1, funit2 !file units
	character(len=100) fName1, fName2 !file names
	character(len=100) fmt,fmt1
	
	!diagnostics for determining when to do eigen analysis
	integer neVol
	parameter(neVol=4)
	double precision, dimension(:), allocatable :: totVol, x1, x2, y1, y2, slope, intcpt, cVolFrac, pVolFrac
	double precision, dimension(:,:,:), allocatable :: eVolFrac
	
	!info for output file update
	double precision, allocatable :: evData(:,:), evDataAll(:), evDataTemp(:)
	
	!isolated cluster info
	integer ic_n !no. of nodes
	integer, allocatable :: ic_sc(:), ic_npt(:)
	logical, allocatable :: ic_done(:)
	integer, dimension(:), allocatable :: ic_fNode, ic_nsc, ic_nBrnch
	double precision, dimension(:,:,:),  allocatable :: ic_brnch, ic_llimits, ic_plimits
	double precision, allocatable :: ic_climits(:,:,:), ic_volFac(:)
	double precision, dimension(:),  allocatable :: ic_Z, ic_Zold, ic_info, ic_vnow, ic_hilike, ic_inc
	double precision, dimension(:,:),  allocatable :: ic_eff
	logical, dimension(:),  allocatable :: ic_reme, ic_rFlag, ic_chk
 	logical modeFound
	
	!means & standard deviations of the live points (for prior edge detection)
	double precision, dimension(:,:), allocatable :: ic_mean, ic_sigma
	double precision lPts(ndims)
	
	!sub-cluster properties
	integer sc_n !no. of sub-clusters
	integer, dimension(:), allocatable :: sc_npt, nptk, nptx ,sc_node, nodek, sck
	double precision, dimension(:,:), allocatable :: meank, sc_eval, evalk
	double precision, dimension(:,:,:), allocatable :: sc_invcov, invcovk, sc_evec, eveck, tMatk
	double precision, dimension(:), allocatable :: kfack, volk, effk
	double precision, allocatable :: sc_mean(:,:), sc_tmat(:,:,:), sc_kfac(:), sc_eff(:), sc_vol(:)
	
	!auxiliary points (to be re-arranged with main points during clustering)
	integer naux !dimensionality of aux points 
	double precision, dimension(:,:), allocatable :: aux, pt
	
	!rejected point info
	double precision lowlike !lowest log-like
	double precision, allocatable :: lowp(:), lowphyP(:) !point with the lowlike
	integer indx(1) !point no. of lowlike
	
	!new point
	double precision lnew
	double precision, allocatable :: pnew(:), phyPnew(:) ! new point
	double precision, dimension(:,:,:), allocatable :: pnewa, phyPnewa
	double precision, dimension(:,:), allocatable :: lnewa
	integer, dimension(:), allocatable :: rIdx
	integer, dimension(:,:), allocatable :: sEll
	logical, dimension(:), allocatable :: remain
	
	!mode separation
	integer nCdim
	
	!importance sampling
	double precision, allocatable :: IS_allpts(:,:), IS_iterinfo(:,:), IS_V(:)
	integer IS_counter(7)
	integer :: IS_nstore = 10000, IS_nMC = 1000
	double precision :: IS_Z
	logical :: IS_CheckAll = .false., IS_betterMC = .true.
	      
	INTERFACE
    		!the likelihood function
    		subroutine loglike(Cube,n_dim,nPar,lnew,context_pass)
			integer n_dim,nPar,context_pass
			double precision lnew,Cube(nPar)
		end subroutine loglike
      	end INTERFACE
	
	INTERFACE
		!the user dumper function
    		subroutine dumper(nSamples, nlive, nPar, physLive, posterior, paramConstr, maxLogLike, logZ, logZerr, context_pass)
			integer nSamples, nlive, nPar, context_pass
			double precision, pointer :: physLive(:,:), posterior(:,:), paramConstr(:)
			double precision maxLogLike, logZ, logZerr
		end subroutine dumper
	end INTERFACE
	
	
	allocate( eswitchff(maxCls), escount(maxCls), dmin(maxCls) )
	allocate( evData(updInt,totPar+3) )
	allocate( ic_sc(maxCls), ic_npt(maxCls) )
	allocate( ic_done(0:maxCls) )
	allocate( ic_climits(maxCls,ndims,2), ic_volFac(maxCls) )
	allocate( sc_mean(maxeCls,ndims), sc_tmat(maxeCls,ndims,ndims), sc_kfac(maxeCls), sc_eff(maxeCls), sc_vol(maxeCls) )
	allocate( lowp(ndims), lowphyP(totPar) )
	allocate( pnew(ndims), phyPnew(totPar) )
	
	
	!initializations
	ic_done=.false.
	ic_npt=nlive
	ic_climits(:,:,2)=1d0
	ic_climits(:,:,1)=0d0
	ic_volFac(:)=1d0
	
	if(my_rank==0) then
		!memory allocation
		if( IS ) then
			allocate(IS_allpts(nlive+IS_nstore,ndims+6), IS_iterinfo(nlive+IS_nstore/10,5), IS_V(maxeCls))
			IS_allpts = 0d0
			IS_iterinfo = 0d0
			
			do i = 1, nlive
				IS_allpts(i,1:ndims) = p(1:ndims,i)	!point
			enddo
			IS_allpts(1:nlive,ndims+1) = l(1:nlive)		!likelihood
			IS_allpts(1:nlive,ndims+2) = dble(nlive)	!p(\Theta) n = \Sum_{i}^{niter} n_{i} E_{i}(\Theta) / V_{i}
			IS_allpts(1:nlive,ndims+3) = 1d0		!check this point for ellipsoid membership in later iterations
			IS_allpts(1:nlive,ndims+4) = 1d0		!which ellipsoid this point lies in
			IS_allpts(1:nlive,ndims+5) = 1d0		!Mahalanobis distance
			IS_allpts(1:nlive,ndims+6) = 1d0		!node
			
			IS_iterinfo(1:nlive,1) = 1d0			!volume
			IS_iterinfo(1:nlive,2) = 1d0			!total no. of points collected at each iteration, excluding points outside prior
			IS_iterinfo(1:nlive,4) = 1d0			!total no. of points collected at each iteration, including points outside prior
			IS_iterinfo(1:nlive,3) = 1d0			!effective no. of points collected at each iteration
			IS_iterinfo(1:nlive,5) = 1d0			!node
			
			IS_counter(1) = nlive				!total no. of points collected so far
			IS_counter(2) = 1				!first point to be checked for ellipsoid membership in later iterations
			IS_counter(3) = nlive+IS_nstore			!total no. of points that can be stored in IS_allpts array
			IS_counter(4) = nlive+IS_nstore/10		!total no. of points that can be stored in IS_iterinfo array
			IS_counter(5) = nlive				!total no. of iterations done so far
			IS_counter(6) = 0				!total no. of IS_iterinfo members written to output file
			IS_counter(7) = 0				!total no. of IS_allpts members written to output file
		endif
		
		allocate(evDataAll(1))
		allocate(sc_npt(maxeCls), nptk(maxeCls), nptx(nlive), meank(maxeCls,ndims), &
		sc_eval(maxeCls,ndims), evalk(maxeCls,ndims), sc_invcov(maxeCls,ndims,ndims), &
		invcovk(maxeCls,ndims,ndims), tMatk(maxeCls,ndims,ndims), &
		sc_evec(maxeCls,ndims,ndims), eveck(maxeCls,ndims,ndims), kfack(maxeCls), &
		effk(maxeCls), volk(maxeCls), &
 		sc_node(maxeCls),nodek(maxeCls),sck(maxCls))
		allocate(pt(ndims,nlive), aux(ndims+totPar+4-nCdims,nlive))
		allocate(ic_fNode(maxCls),ic_nsc(maxCls),ic_nBrnch(maxCls), &
		ic_brnch(maxCls,maxCls,2),ic_reme(maxCls),ic_rFlag(maxCls),ic_z(maxCls),ic_zold(maxCls),ic_info(maxCls), &
		ic_vnow(maxCls),ic_hilike(maxCls),ic_inc(maxCls),ic_chk(maxCls),ic_llimits(maxCls,ndims,2))
		if(multimodal) allocate(ic_plimits(maxCls,nCdims,2))
		if(ceff) then
			allocate(ic_eff(maxCls,4))
			ic_eff(:,1:2)=0d0
			ic_eff(:,3)=1d0
			ic_eff(:,3)=ef
			ic_eff(:,4)=1d0
		endif
		allocate(pnewa(maxCls,mpi_nthreads,ndims),phyPnewa(maxCls,mpi_nthreads,totPar),lnewa(maxCls,mpi_nthreads), &
		sEll(maxCls,mpi_nthreads),remain(maxCls),rIdx(maxCls))
		allocate(totVol(maxCls),eVolFrac(maxCls,neVol*2,neVol),x1(maxCls),x2(maxCls),y1(maxCls), &
		y2(maxCls),slope(maxCls),intcpt(maxCls),cVolFrac(maxCls),pVolFrac(maxCls))
		if( prior_warning ) allocate( ic_mean(maxCls, ndims), ic_sigma(maxCls, ndims) )
	
		!global logZ = log(0)
		gZ=logZero
		gZOld=logZero
		ic_Z=logZero
		ic_zold=logZero
		ginfo=0.d0
		ic_info=0.d0
		
		ic_llimits(:,:,2)=1d0
		ic_llimits(:,:,1)=0d0
		if(multimodal) then
			ic_plimits(:,:,2)=huge(1d0)
			ic_plimits(:,:,1)=-huge(1d0)
		endif
		
		!just one node to sttart with
		ic_n=1
		ic_fNode(1)=0
		ic_sc(1)=0
		ic_nsc=nsc_def
		ic_reme=.false.
		ic_nBrnch=0
	
		!set the prior volume
		ic_vnow=0.d0
      		ic_vnow(1)=1.d0
	
		!no ellipsoidal sampling to start with
		eswitch=.false.
		peswitch=.false.
		escount=0
      	
		!no leftover points
		remain=.false.
		rIdx=1
	
		!no sub-clusters at the beginning
		sc_n=0
		totVol=1.d0
		sc_node=1
	
		!eigen analysis diagnostics
		eVolFrac=0.d0
	
		dmin=ndims+1
		ic_rFlag=.false.
		sff=0
	
		if(resumeFlag) then
			!read the resume file
			funit1=u_resume
	
			open(unit=funit1,file=resumename,status='old')
                	  	
			read(funit1,*)genLive
			read(funit1,*)globff,numlike,ic_n,nlive
			read(funit1,*)gZ,ginfo
			ic_rFlag(1:ic_n)=.true.
			read(funit1,*)eswitch
			peswitch=eswitch
			if(eswitch) totVol=1.d0
			
			!read branching info
	            	do i=1,ic_n
            			read(funit1,*)ic_nBrnch(i)
				if(ic_nBrnch(i)>0) read(funit1,*)ic_brnch(i,1:ic_nBrnch(i),1),ic_brnch(i,1:ic_nBrnch(i),2)
			enddo
				
			!read the node info
			do i=1,ic_n
				read(funit1,*)ic_done(i),ic_reme(i),ic_fNode(i),ic_npt(i)
				read(funit1,*)ic_vnow(i),ic_Z(i),ic_info(i)
				if(ceff) then
					read(funit1,*)ic_eff(i,4)
					ic_eff(i,3)=ic_eff(i,4)
				endif
			enddo
			if(.not.eswitch .or. ic_n==1) ic_npt(1)=nlive
                  	close(funit1)
			eswitchff=0
			
    			!read hypercube-live file
			open(unit=u_live,file=livename,status='old')
			i=0
			do
      				i=i+1
				read(u_live,*,IOSTAT=iostatus) p(1:ndims,i),l(i)
            			if(iostatus<0) then
            				i=i-1
                  			if(i<nlive) then
                  				write(*,*)"ERROR: live points file has less than ",nlive," points."
						write(*,*)"Aborting"
#ifdef MPI
						call MPI_ABORT(MPI_COMM_WORLD,errcode)
#endif
                        			stop
					endif
                  			exit
				endif
				if(i>nlive) then
					write(*,*)"ERROR: live points file has greater than ",nlive," points."
					write(*,*)"Aborting"
#ifdef MPI
					call MPI_ABORT(MPI_COMM_WORLD,errcode)
#endif
                        		stop
				endif
			enddo
			close(u_live)
			
    			!read physical-live file
			open(unit=u_phys,file=physname,status='old')
			i=0
			do
      				i=i+1
				read(u_phys,*,IOSTAT=iostatus) phyP(1:totPar,i),d1,j
            			if(iostatus<0) then
            				i=i-1
                  			if(i<nlive) then
                  				write(*,*)"ERROR: phys live points file has less than ",nlive," points."
						write(*,*)"Aborting"
#ifdef MPI
						call MPI_ABORT(MPI_COMM_WORLD,errcode)
#endif
                        			stop
					endif
                  			exit
				endif
				if(i>nlive) then
					write(*,*)"ERROR: phys live points file has greater than ",nlive," points."
					write(*,*)"Aborting"
#ifdef MPI
					call MPI_ABORT(MPI_COMM_WORLD,errcode)
#endif
                        		stop
				endif
			enddo
			close(u_phys)
			
			!read the IS files
			if( IS ) then
				!read all the points collected so far
				open(unit=u_IS(2),file=IS_Files(2),form='unformatted',access='sequential',status='old')
				read(u_IS(2))j,IS_counter(2),k
				call ExtendArrayIfRequired(IS_counter(1),j-IS_counter(1),IS_counter(3),IS_nstore,ndims+6,IS_allpts)
				IS_counter(1) = j
				IS_counter(7) = j
				call ExtendArrayIfRequired(IS_counter(5),k-IS_counter(5),IS_counter(4),IS_nstore/10,5,IS_iterinfo)
				IS_counter(5) = k
				IS_counter(6) = k
				
				do i = 1, IS_counter(1)
					read(u_IS(2),IOSTAT=iostatus)IS_allpts(i,ndims+2),IS_allpts(i,ndims+3)
					
					!end of file?
					if(iostatus<0) then
						write(*,*)"ERROR: Not enough points in ",IS_Files(2)
						write(*,*)"Aborting"
#ifdef MPI
						call MPI_ABORT(MPI_COMM_WORLD,errcode)
#endif
						stop
					endif
				enddo
				close(u_IS(2))
				
				open(unit=u_IS(1),file=IS_Files(1),form='unformatted',access='sequential',status='old')
				do i = 1, IS_counter(1)
					read(u_IS(1),IOSTAT=iostatus)IS_allpts(i,1:ndims+1),IS_allpts(i,ndims+6)
					
					!end of file?
					if(iostatus<0) then
						write(*,*)"ERROR: Not enough points in ",IS_Files(1)
						write(*,*)"Aborting"
#ifdef MPI
						call MPI_ABORT(MPI_COMM_WORLD,errcode)
#endif
						stop
					endif
				enddo
				close(u_IS(1))
				
				!read the iteration info
				open(unit=u_IS(3),file=IS_Files(3),form='unformatted',access='sequential',status='old')
				do i = 1, IS_counter(5)
					read(u_IS(3),IOSTAT=iostatus)IS_iterinfo(i,1:5)
					
					!end of file?
					if(iostatus<0) then
						write(*,*)"ERROR: Not enough points in ",IS_Files(3)
						write(*,*)"Aborting"
#ifdef MPI
						call MPI_ABORT(MPI_COMM_WORLD,errcode)
#endif
						stop
					endif
				enddo
				close(u_IS(3))
			endif
		endif
		
		!find the highest likelihood for each node
		j=0
		ic_done(0)=.true.
		do i=1,ic_n
			ic_hilike(i)=maxval(l(j+1:j+ic_npt(i)))
			lowlike=minval(l(j+1:j+ic_npt(i)))
			ic_inc(i)=ic_hilike(i)+log(ic_vnow(i))-ic_Z(i)
			if(ic_npt(i)<ndims+1 .or. abs(lowlike-ic_hilike(i))<= 0.0001d0 .or. (ic_inc(i)<log(tol) .and. globff-nlive>50)) then
				ic_done(i)=.true.
			else
				ic_done(i)=.false.
				ic_done(0)=.false.
			endif
			
			!set the limits
			if(.not.ic_done(i)) then
				ic_llimits(i,:,1)=p(:,j+1)
				ic_llimits(i,:,2)=p(:,j+1)
				if(multimodal) then
					ic_plimits(i,1:nCdims,1)=phyP(1:nCdims,j+1)
					ic_plimits(i,1:nCdims,2)=phyP(1:nCdims,j+1)
				endif
				do k=j+2,j+ic_npt(i)
					if(multimodal) then
						call setLimits(multimodal,ndims,nCdims,ic_llimits(i,:,:), &
						ic_plimits(i,:,:),p(:,k),phyP(1:nCdims,k),ic_climits(i,:,:))
					else
						call setLimits(multimodal,ndims,nCdims,ic_llimits(i,:,:), &
						ic_climits(i,:,:),p(:,k),phyP(1:nCdims,k),ic_climits(i,:,:))
					endif
				enddo
			endif
			
			j=j+ic_npt(i)
		enddo
	endif
    	
#ifdef MPI
	call MPI_BARRIER(MPI_COMM_WORLD,errcode)
	call MPI_BCAST(ic_n,1,MPI_INTEGER,0,MPI_COMM_WORLD,errcode)
	call MPI_BCAST(ic_npt(1:ic_n),ic_n,MPI_INTEGER,0,MPI_COMM_WORLD,errcode)
#endif
		
	do ff=1,maxIter

#ifdef MPI
    		call MPI_BARRIER(MPI_COMM_WORLD,errcode)
    		call MPI_BCAST(ic_done(0:ic_n),ic_n+1,MPI_LOGICAL,0,MPI_COMM_WORLD,errcode)
#endif	
		!stopping condition reached
		if(ic_done(0)) then
			if(my_rank==0) then
                        
				if(outfile) then
	                        	!write the resume file
	                        	funit1=u_resume
					fName1=resumename
	                        	write(fmt,'(a,i5,a)')  '(',totPar+1,'E28.18,i4)'
					open(unit=funit1,file=fName1,form='formatted',status='replace')
					write(funit1,'(l2)').false.
					write(funit1,'(4i12)')globff,numlike,ic_n,nlive
					write(funit1,'(2E28.18)')gZ,ginfo
					write(funit1,'(l2)')eswitch
	            			!write branching info
		            		do i=1,ic_n
	            				write(funit1,'(i4)')ic_nBrnch(i)
						if(ic_nBrnch(i)>0) then
							write(fmt,'(a,i5,a)')  '(',2*ic_nBrnch(i),'E28.18)'
							write(funit1,fmt)ic_brnch(i,1:ic_nBrnch(i),1),ic_brnch(i,1:ic_nBrnch(i),2)
						endif
					enddo
					!write the node info
					do i=1,ic_n
						write(funit1,'(2l2,2i6)')ic_done(i),ic_reme(i),ic_fNode(i),ic_npt(i)
						write(funit1,'(3E28.18)')ic_vnow(i),ic_Z(i),ic_info(i)
						if(ceff) write(funit1,'(1E28.18)')ic_eff(i,4)
					enddo
	                  		close(funit1)
				endif
				
				if( IS ) then
					IS_Z = logZero
					do j = 1, IS_counter(1)
						d1 = IS_allpts(j,ndims+1)-log(IS_allpts(j,ndims+2))
						IS_Z = LogSumExp(IS_Z, d1)
					enddo
				endif
				
				!fback
                		if(fback) call gfeedback(gZ,IS,IS_Z,numlike,globff,.false.)
				
				call pos_samp(Ztol,globff,broot,nlive,ndims,nCdims,totPar,multimodal,outfile,gZ,ginfo,ic_n,ic_Z(1:ic_n), &
				ic_info(1:ic_n),ic_reme(1:ic_n),ic_vnow(1:ic_n),ic_npt(1:ic_n),ic_nBrnch(1:ic_n),ic_brnch(1:ic_n,:,1),phyP(:,1:nlive), &
				l(1:nlive),evDataAll,IS,IS_Z,dumper,context)
				
				!if done then add in the contribution to the global evidence from live points
				j=0
				do i1=1,ic_n
					!live point's contribution to the evidence
                			logX=log(ic_vnow(i1)/dble(ic_npt(i1)))
					do i=j+1,j+ic_npt(i1)
						d1=l(i)+logX
						!global evidence
            					gZold=gZ
                  				gZ=LogSumExp(gZ,d1)
						!ginfo=exp(d1-gZ)*l(i)+exp(gZold-gZ)*(ginfo+gZold)-gZ
						ginfo=ginfo*exp(gZold-gZ)+exp(d1-gZ)*l(i)
						!local evidence
            					!gZold=ic_Z(i1)
						ic_Zold(i1)=ic_Z(i1)
                  				ic_Z(i1)=LogSumExp(ic_Z(i1),d1)
						!ic_info(i1)=exp(d1-ic_Z(i1))*l(i)+ &
						!exp(gZold-ic_Z(i1))*(ic_info(i1)+gZold)-ic_Z(i1)
						ic_info(i1)=ic_info(i1)*exp(ic_zold(i1)-ic_Z(i1))+exp(d1-ic_Z(i1))*l(i)
					enddo
					j=j+ic_npt(i1)
                		enddo
				
				ginfo=ginfo-gZ
			
				!memory deallocation
				deallocate(sc_npt, sc_invcov, sc_node, sc_eval, sc_evec)
 				deallocate(nptk, nptx, meank, evalk, invcovk, tMatk, eveck, kfack, effk, volk, &
				nodek ,sck)
				deallocate(pt, aux)
				deallocate(ic_fNode,ic_nsc,ic_nBrnch,ic_brnch,ic_reme,ic_rFlag,ic_Z,ic_zold, &
				ic_info,ic_vnow,ic_hilike,ic_inc,ic_chk,ic_llimits)
				if(multimodal) deallocate(ic_plimits)
				if(ceff) deallocate(ic_eff)
				deallocate(pnewa,phyPnewa,lnewa,sEll,remain,rIdx)
				deallocate(totVol,eVolFrac,x1,x2,y1,y2,slope,intcpt,cVolFrac,pVolFrac)
				if( prior_warning ) deallocate( ic_mean, ic_sigma )
				deallocate( evDataAll )
				if( IS ) deallocate(IS_allpts,IS_iterinfo,IS_V)
			endif
			
			deallocate( eswitchff, escount, dmin, evData, ic_sc, ic_npt, ic_done, ic_climits, &
			ic_volFac, sc_mean, sc_tmat, sc_kfac, sc_eff, sc_vol, lowp, lowphyP, pnew, phyPnew )
			return
		endif
		
		modeFound=.false.
		
		if(my_rank==0) then
			!mode separation
			if(ff/=1 .and. eswitch .and. multimodal .and. mod(ff,15)==0 .and. ic_n<maxCls) then
			
!				if(mod(ff,30)==0) then
					nCdim=nCdims
!				else
!					nCdim=3
!				endif
				
				!aux information to be re-arranged with the live points
				naux=ndims+totPar+1-nCdim
				aux(1,1:nlive)=l(1:nlive)
				aux(2:totPar+1,1:nlive)=phyP(1:totPar,1:nlive)
				aux(totPar+2:naux,1:nlive)=p(nCdim+1:ndims,1:nlive)
				
				!save old no. of modes
				i=ic_n
				nptk(1:i)=ic_npt(1:i)
				
				!decide which nodes to check for mode separation
				ic_chk(1:ic_n)=.not.(ic_done(1:ic_n))
				do j=1,ic_n
					if(ic_chk(j)) then
						if(ic_inc(j)<log(0.5d0)) ic_chk(j)=.false.
					endif
				enddo
				
				modeFound=isolateModes2(nlive,ndims,nCdim,p(1:nCdim,1:nlive),naux,aux(1:naux,:), &
				ic_n,ic_fnode,ic_npt,ic_reme,ic_chk(1:ic_n),ic_vnow(1:ic_n),ic_rFlag,ic_climits(:,1:nCdim,:))

				if(modeFound) then			
					!re-arrange info
					l(1:nlive)=aux(1,1:nlive)
					phyP(1:totPar,1:nlive)=aux(2:totPar+1,1:nlive)
					p(nCdim+1:ndims,1:nlive)=aux(totPar+2:naux,1:nlive)
					
					if(nCdims<ndims .or. .true.) then
						ic_sc(i+1:ic_n)=1
						sc_npt(sc_n+1:sc_n+ic_n-i)=ic_npt(i+1:ic_n)
						sc_n=sc_n+ic_n-i
					endif
					
					!new nodes should inherit nsc from their parents
					do j=i+1,ic_n
						ic_nsc(j)=nsc_def
						ic_nsc(ic_fNode(j))=nsc_def
						
						!branching info
						ic_nBrnch(ic_fNode(j))=ic_nBrnch(ic_fNode(j))+1
						ic_brnch(ic_fNode(j),ic_nBrnch(ic_fNode(j)),1)=dble(j)
						ic_brnch(ic_fNode(j),ic_nBrnch(ic_fNode(j)),2)=dble(ic_npt(j))/dble(nptk(ic_fNode(j)))
					enddo
					
					!assign prior volumes
					nodek(1:i)=0
					!first the new nodes
					q=sum(ic_npt(1:i))
					do j=i+1,ic_n
						if(ic_rFlag(j)) then
							!parent node
							k=ic_fNode(j)
							!divide prior volume
							ic_vnow(j)=ic_vnow(k)*dble(ic_npt(j))/dble(nptk(k))
							!divide the evidences & info
							ic_zold(j)=ic_z(k)
							ic_Z(j)=log(dble(ic_npt(j))/dble(nptk(k)))+ic_Z(k)
							ic_info(j)=ic_info(k)*dble(ic_npt(j))/dble(nptk(k))*exp(ic_zold(j)-ic_z(j))
							nodek(k)=1 !note the parent node
							ic_hilike(j)=maxval(l(q+1:q+ic_npt(j)))
							if(ic_hilike(j)-minval(l(q+1:q+ic_npt(j)))<= 0.0001) ic_done(j)=.true.
							if(ceff) then
								ic_eff(j,1:2)=0d0
								ic_eff(j,3:4)=ic_eff(k,3:4)
							endif
							
							!ellipsoidal decomposition prediction variables
							eVolFrac(j,:,:)=0d0
							eVolFrac(k,:,:)=0d0
							
							!set the limits
							ic_climits(j,:,:)=ic_climits(k,:,:)
							ic_llimits(j,:,:)=ic_llimits(k,:,:)
							ic_plimits(j,:,:)=ic_plimits(k,:,:)
							ic_volFac(j)=ic_volFac(k)
						endif
						q=q+ic_npt(j)
					enddo
					!now the father nodes
					q=0
					m=0
					do j=1,i
						if(nodek(j)==1) then
							ic_vnow(j)=ic_vnow(j)*dble(ic_npt(j))/dble(nptk(j))
							if(ic_npt(j)==0) then
								ic_Z(j)=logZero
								ic_info(j)=0.d0
							else
								ic_zold(j)=ic_z(j)
								ic_Z(j)=log(dble(ic_npt(j))/dble(nptk(j)))+ic_Z(j)
								ic_info(j)=ic_info(j)*dble(ic_npt(j))/dble(nptk(j))*exp(ic_zold(j)-ic_z(j))
							endif
							ic_hilike(j)=maxval(l(q+1:q+ic_npt(j)))
							if(ic_hilike(j)-minval(l(q+1:q+ic_npt(j)))<= 0.0001) ic_done(j)=.true.
							if(ceff) ic_eff(j,1:2)=0d0
							
							sc_npt(m+1)=ic_npt(j)
							sc_npt(m+2:m+ic_sc(j))=0
						endif
						
						!no leftover points
						remain(j)=.false.
						rIdx(j)=1
						
						q=q+ic_npt(j)
						m=m+ic_sc(j)
					enddo
					ic_rFlag(1:ic_n)=.true.
				endif
			endif
			
			num_old=numlike
		
			!ellipsoidal decomposition
			flag2=.false.

			if(peswitch) then
				q=0 !no. of sub-clusters traversed
				m=0 !no. of points traversed
				n=0 !no. of sub-clusters created so far
				do i1=1,ic_n
					if(ic_npt(i1)==0) then
						totvol(i1) = 0d0
						sck(i1)=0
						q=q+ic_sc(i1) !no. of sub-clusters traversed
						cycle
					endif
					
					if(ic_done(i1)) then
						sck(i1)=1
						q=q+ic_sc(i1) !no. of sub-clusters traversed
						m=m+ic_npt(i1)
						nptk(n+1)=ic_npt(i1)
						nodek(n+1)=i1
						n=n+1
						cycle
					endif
					
					if(ic_npt(i1)>=ndims+1 .and. totVol(i1)==0.d0) ic_rFlag(i1)=.true.
					
					if(.not.eswitch) then
						if(mod(ff-1-eswitchff(i1),ic_nsc(i1))==0) then
							flag=.true.
						else
							flag=.false.
						endif
					elseif(ic_npt(i1)<ndims+1 .and. .not.ic_rFlag(i1)) then
						flag=.false.
					elseif(ic_rFlag(i1)) then
						flag=.true.
					else
						flag=.false.
						
						if(ceff) then
							d4=ic_eff(i1,3)
						else
							d4=ef
						endif
						
						!current vol fraction
						cVolFrac(i1)=totVol(i1)*d4/((ic_vnow(i1)*ic_volFac(i1)))
				
						!predicted vol fraction
!						if(dino) then
!							pVolFrac(i1)=max(slope(i1)*dble(globff)+intcpt(i1),1.d0)
!						else
!							pVolFrac(i1)=slope(i1)*dble(globff)+intcpt(i1)
!						endif
						
						pVolFrac(i1)=0d0
						i3=0
						do i2=1,neVol
							if(eVolFrac(i1,i2,1)<=0d0) exit
							pVolFrac(i1)=pVolFrac(i1)+eVolFrac(i1,i2,1)
							i3=i3+1
						enddo
						if(pVolFrac(i1)==0d0) then
							pVolFrac(i1)=1d0
						else
							pVolFrac(i1)=pVolFrac(i1)/dble(i3)
						endif

						if(.not.flag .and. (cVolFrac(i1)>1.1 .or. .not.dino) .and.  &
						(pVolFrac(i1)<cVolFrac(i1) .or. mod(ff-1-eswitchff(i1),ic_nsc(i1))==0)) flag=.true.
				
						if(flag .and. dino) then
							do i=q+1,q+ic_sc(i1)
								if(sc_eff(i)<=1.00001) exit
								if(i==sc_n) flag=.false.
							enddo
						endif
	
						if(.not.flag .and. eVolFrac(i1,2*neVol,1)==0.d0 .and. mod(ff-1-eswitchff(i1),ic_nsc(i1))==0) flag=.true.
						!if(mod(ff-1-eswitchff(i1),20)==0) flag=.true.
					endif
					
					if(flag) then
						!target vol
						if(ceff) then
							d4=ic_eff(i1,3)
						else
							d4=ef
						endif

						d5=1d0
						do i2=1,ndims
							d5=d5/(ic_llimits(i1,i2,2)-ic_llimits(i1,i2,1))
						enddo
						d1=ic_vnow(i1)*d5/d4
				
						!volume threshold
						if(ic_rFlag(i1)) then
							d2=1000.*d1
						else
							d2=totVol(i1)*d5/ic_volFac(i1)
						endif
					
						!volume tolerance
						d3=0.d0
					
						eswitchff(i1)=ff-1
						
						flag=.false.
						
						!aux information to be re-arranged with the live points
						naux=totPar+2
						!rescaling
						do i3=1,ndims
							!rescale back into unit hypercube
							d4=ic_climits(i1,i3,2)-ic_climits(i1,i3,1)
							pt(i3,1:ic_npt(i1))=ic_climits(i1,i3,1)+d4*p(i3,m+1:m+ic_npt(i1))
							
							if(pWrap(i3)) then
								do i4=1,ic_npt(i1)
									call wraparound(pt(i3,i4),pt(i3,i4))
								enddo
							endif
							
							!scale to the new limits
							d4=ic_llimits(i1,i3,2)-ic_llimits(i1,i3,1)
							pt(i3,1:ic_npt(i1))=(pt(i3,1:ic_npt(i1))-ic_llimits(i1,i3,1))/d4
						enddo
						aux(1,1:ic_npt(i1))=l(m+1:m+ic_npt(i1))
						aux(2:totPar+1,1:ic_npt(i1))=phyP(1:totPar,m+1:m+ic_npt(i1))
						lowlike=minval(l(m+1:m+ic_npt(i1)))
						aux(naux,1:ic_npt(i1))=(l(m+1:m+ic_npt(i1))-lowlike)/(ic_hilike(i1)-lowlike)
						
						!max no. of sub-clusters allowed
						n2=min(ceiling(dble(ic_npt(i1))/dble(min_pt)),maxeCls-n)
						
						if(ic_npt(i1)<nlive/1.5) then
							cSwitch=.false.
						else
							cSwitch=.false.
						endif

						if(dino) then
							d5=d1
							do i3=1,50
								count=count+1
							
								!make dinosaur (sub-clustering)
								if(Dinosaur(pt(:,1:ic_npt(i1)),ic_npt(i1),ndims,sck(i1),nptk(n+1:n+n2), &
								naux,aux(1:naux,1:ic_npt(i1)),min_pt,n2,meank(n+1:n+n2,:),invcovk(n+1:n+n2,:,:), &
								tMatk(n+1:n+n2,:,:),evalk(n+1:n+n2,:),eveck(n+1:n+n2,:,:),kfack(n+1:n+n2), &
								effk(n+1:n+n2),volk(n+1:n+n2),d1,d2,neVol,eVolFrac(i1,:,:),globff,d3,.false., &
								ic_rFlag(i1),cSwitch,nCdims)) then
									if(eswitch) then
										ic_nsc(i1)=max(1,ic_nsc(i1)-10)
									else
										eswitch=.true.
										ic_sc(1)=sck(i1)
									endif
								
									scount=scount+1
									totVol(i1)=sum(volk(n+1:n+sck(i1)))
								
									flag=.true.
									flag2=.true.
								
									!no leftover points
									remain(i1)=.false.
									rIdx(i1)=1
									
									!set the limits
									ic_climits(i1,:,:)=ic_llimits(i1,:,:)
									
									ic_volFac(i1)=1d0
									do i2=1,ndims
										ic_volFac(i1)=ic_volFac(i1)/(ic_climits(i1,i2,2)-ic_climits(i1,i2,1))
									enddo
									
									!eVolFrac(i1,1,1)=eVolFrac(i1,1,1)*d1/d5
									
									if(ceff) then
										!eVolFrac(i1,1,1)=eVolFrac(i1,1,1)/ic_eff(i1,3)
										ic_eff(i1,4)=ic_vnow(i1)*ic_volFac(i1)/totVol(i1)
									else
										!eVolFrac(i1,1,1)=eVolFrac(i1,1,1)/ef
									endif
									
									exit
								elseif(.not.ic_rFlag(i1)) then
									if(eswitch) ic_nsc(i1)=ic_nsc(i1)+10
									exit
								endif
								
								eVolFrac(i1,1,1)=eVolFrac(i1,1,1)*d1/d5
							
								if(mod(i3,5)==0) then
									d2=d2*2.d0
									d1=d1*2.d0
								endif
							enddo
						endif
						
						if(.not.Flag .and. ic_rFlag(i1)) then
							sck(i1)=1
							nptk(n+1)=ic_npt(i1)
							kfack(n+1)=0.d0
							volk(n+1)=0.d0
							effk(n+1)=1.d0
							totVol(i1)=0.d0
							flag=.true.
							ic_done(i1)=.true.
							ic_done(0)=.true.
							do i3=1,ic_n
								if(.not.ic_done(i3)) then
									ic_done(0)=.false.
									exit
								endif
							enddo
						else
							!predict the current vol fraction
							x1(i1)=sum(eVolFrac(i1,neVol+1:2*neVol,2))/neVol
							y1(i1)=sum(eVolFrac(i1,neVol+1:2*neVol,1))/neVol
							x2(i1)=sum(eVolFrac(i1,1:neVol,2))/neVol
							y2(i1)=sum(eVolFrac(i1,1:neVol,1))/neVol
							slope(i1)=(y2(i1)-y1(i1))/(x2(i1)-x1(i1))
							intcpt(i1)=(x2(i1)*y1(i1)-x1(i1)*y2(i1))/(x2(i1)-x1(i1))
						endif
					endif
					
					if(ic_rFlag(i1) .and. .not.eswitch) exit
						
					if(flag) then
						!aux information to be re-arranged with the live points
						p(:,m+1:m+ic_npt(i1))=pt(:,1:ic_npt(i1))
						l(m+1:m+ic_npt(i1))=aux(1,1:ic_npt(i1))
						phyP(1:totPar,m+1:m+ic_npt(i1))=aux(2:totPar+1,1:ic_npt(i1))
					else
						sck(i1)=ic_sc(i1)
						meank(n+1:n+sck(i1),:)=sc_mean(q+1:q+sck(i1),:)
						invcovk(n+1:n+sck(i1),:,:)=sc_invcov(q+1:q+sck(i1),:,:)
						tMatk(n+1:n+sck(i1),:,:)=sc_tMat(q+1:q+sck(i1),:,:)
						evalk(n+1:n+sck(i1),:)=sc_eval(q+1:q+sck(i1),:)
						eveck(n+1:n+sck(i1),:,:)=sc_evec(q+1:q+sck(i1),:,:)
						kfack(n+1:n+sck(i1))=sc_kfac(q+1:q+sck(i1))
						effk(n+1:n+sck(i1))=sc_eff(q+1:q+sck(i1))
						volk(n+1:n+sck(i1))=sc_vol(q+1:q+sck(i1))
						nptk(n+1:n+sck(i1))=sc_npt(q+1:q+sck(i1))
						
						slope(i1)=slope(i1)*1.01
					endif
					nodek(n+1:n+sck(i1))=i1
					
					q=q+ic_sc(i1) !no. of sub-clusters traversed
					m=m+ic_npt(i1) !no. of points traversed
					n=n+sck(i1) !no. of sub-clusters created so far
				enddo
				
				if(flag2) then
					!update mode info
					ic_sc(1:ic_n)=sck(1:ic_n)
					!update sub-cluster info
					sc_n=sum(ic_sc(1:ic_n))
					sc_mean(1:sc_n,:)=meank(1:sc_n,:)
					sc_invcov(1:sc_n,:,:)=invcovk(1:sc_n,:,:)
					sc_tMat(1:sc_n,:,:)=tMatk(1:sc_n,:,:)
					sc_eval(1:sc_n,:)=evalk(1:sc_n,:)
					sc_evec(1:sc_n,:,:)=eveck(1:sc_n,:,:)
					sc_kfac(1:sc_n)=kfack(1:sc_n)
					sc_eff(1:sc_n)=effk(1:sc_n)
					sc_vol(1:sc_n)=volk(1:sc_n)
					sc_npt(1:sc_n)=nptk(1:sc_n)
					sc_node(1:sc_n)=nodek(1:sc_n)
					
					!find the index of the low-like points because of re-arrangement
					indx=minloc(l(1:nlive))
				endif
			endif
			
			ic_rFlag(1:ic_n)=.false.
		endif
		
		
            	if(my_rank==0 .and. eswitch .and. sc_n==0) eswitch=.false.
#ifdef MPI
		call MPI_BARRIER(MPI_COMM_WORLD,errcode)
		call MPI_BCAST(eswitch,1,MPI_LOGICAL,0,MPI_COMM_WORLD,errcode)
		call MPI_BCAST(flag2,1,MPI_LOGICAL,0,MPI_COMM_WORLD,errcode)
		call MPI_BCAST(modeFound,1,MPI_LOGICAL,0,MPI_COMM_WORLD,errcode)
		
		if(modeFound) then
			call MPI_BCAST(ic_n,1,MPI_INTEGER,0,MPI_COMM_WORLD,errcode)
			call MPI_BCAST(ic_npt(1:ic_n),ic_n,MPI_INTEGER,0,MPI_COMM_WORLD,errcode)
    			call MPI_BCAST(ic_done(0:ic_n),ic_n+1,MPI_LOGICAL,0,MPI_COMM_WORLD,errcode)
		endif
		
		if(flag2 .and. eswitch) then
			call MPI_BCAST(sc_n,1,MPI_INTEGER,0,MPI_COMM_WORLD,errcode)
			call MPI_BCAST(ic_sc(1:ic_n),ic_n,MPI_INTEGER,0,MPI_COMM_WORLD,errcode)
			call MPI_BCAST(sc_mean(1:sc_n,1:ndims),sc_n*ndims,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,errcode)
			call MPI_BCAST(sc_tmat(1:sc_n,1:ndims,1:ndims),sc_n*ndims*ndims,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,errcode)
			call MPI_BCAST(sc_kfac(1:sc_n),sc_n,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,errcode)
			call MPI_BCAST(sc_eff(1:sc_n),sc_n,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,errcode)
			call MPI_BCAST(sc_vol(1:sc_n),sc_n,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,errcode)
			call MPI_BCAST(ic_climits(1:ic_n,1:ndims,1:2),ic_n*ndims*2,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,errcode)
			call MPI_BCAST(ic_volFac(1:ic_n),ic_n,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,errcode)
		endif
#endif

		if( IS .and. flag2 .and. eswitch .and. my_rank == 0 ) then
			!estimate the total ellipsoidal volume taking the overlap into account
			
			nd_i = 0 !no. of ellipsoids traversed
			nd_j = 0 !no. of points traversed
			do i = 1, ic_n
				if( ic_sc(i) <= 1 .or. ic_done(i) ) then
					j1 = 1
					m = 1
				else
					m = 0
					j1 = ic_npt(i)
					if( IS_betterMC .and. IS_nMC > nlive ) j1 = int(dble(IS_nMC) * dble(ic_npt(i)) / dble(nlive))
					do j = nd_j+1, nd_j+j1
						if( j > nd_j+ic_npt(i) ) then
							!first pick an ellipsoid according to the vol
							do
								!pick a sub-cluster according to the vol
								call selectEll(ic_sc(i), sc_vol(nd_i+1:nd_i+ic_sc(i)), k)
								k = k + nd_i
								if( sc_kfac(k) == 0d0 .or. sc_vol(k) == 0d0 ) cycle
								exit
							enddo
							
							d1=sc_kfac(k)*sc_eff(k)
							call genPtInEll(ndims, sc_mean(k,:), d1, sc_tMat(k,:,:), my_rank, pnew(1:ndims))
						else
							pnew(1:ndims) = p(1:ndims,j)
						endif
						
						i2 = 0
						do k = nd_i+1, nd_i+ic_sc(i)
		            				if( ptIn1Ell(ndims, pnew(1:ndims), sc_mean(k,:), sc_invcov(k,:,:), sc_kfac(k)*sc_eff(k)) ) i2 = i2+1
		            			enddo
						
						if( j > nd_j+ic_npt(i) ) then
							if( i2 > 1 ) then
								urv = ranmarNS(0)
								if( urv > 1d0/dble(i2) ) then
									i2 = 0
									j1 = j1-1
								endif
							endif
						endif
						m = m+i2
					enddo
				endif
				
				IS_V(i) = ( totVol(i) / ic_volFac(i) ) * dble(j1) / dble(m)
				
				nd_j = nd_j+ic_npt(i)
				nd_i = nd_i+ic_sc(i)
			enddo
			
			
			lowlike = minval(l(1:nlive))
			flag = .false.	!found the very first point which is inside the ellipsoidal bound?
			j1 = IS_counter(2)
			if( IS_CheckAll ) j1 = 1
			do j = j1, IS_counter(1)
				d1 = huge(1d0)
				if( IS_CheckAll .or. IS_allpts(j,ndims+3) == 1d0 ) then
					IS_allpts(j,ndims+3) = 0d0
					
					nd_i = 0 !no. of ellipsoids traversed
					do i = 1, ic_n
						if( ic_done(i) .or. totvol(i) == 0d0 .or. ic_npt(i) == 0 .or. ( multimodal .and. .not.isAncestor(int(IS_allpts(j,ndims+6)), i, ic_fnode(1:i)) ) ) then
							nd_i = nd_i+ic_sc(i)
							cycle
						endif
						
						!apply current limits to this point
						call ApplyLimits(0, ic_climits(i,:,:), IS_allpts(j,1:ndims), pt(1,1:ndims))
						
						do k = nd_i+1, nd_i+ic_sc(i)
							!check if this point is inside the ellipsoid
							call ScaleFactor(1, ndims , pt(1,1:ndims), sc_mean(k,:), sc_invcov(k,:,:), d2)
							if( d2  < sc_kfac(k)*sc_eff(k) .and. d2/(sc_kfac(k)*sc_eff(k)) < d1 ) then
								d1 = d2/(sc_kfac(k)*sc_eff(k))
									
								IS_allpts(j,ndims+3) = 1d0
								IS_allpts(j,ndims+4) = dble(k)	!ellipsoid ID
								IS_allpts(j,ndims+5) = d2	!Mahalanobis distance of this point
									
								if( .not.flag ) then
									flag = .true.
									IS_counter(2) = j
								endif
							endif
						enddo
						
						nd_i = nd_i+ic_sc(i)
					enddo
				endif
			enddo
			
			!not a single point inside the current ellipsoidal decomposition
			if( .not.flag ) IS_counter(2) = IS_counter(1)+1
		endif

		nd_i=0 !no. of ellipsoids traversed
		nd_j=0 !no. of points traversed
		do nd=1,ic_n	
			if(ic_done(nd)) then
				nd_i=nd_i+ic_sc(nd)
				nd_j=nd_j+ic_npt(nd)
				cycle
			endif
			
			if(ic_npt(nd)<ndims+1 .or. (ic_sc(nd)==1 .and. sc_vol(nd_i+1)==0.d0)) then
				ic_done(nd)=.true.
				ic_done(0)=.true.
				do i3=1,ic_n
					if(.not.ic_done(i3)) then
						ic_done(0)=.false.
						exit
					endif
				enddo
				nd_i=nd_i+ic_sc(nd)
				nd_j=nd_j+ic_npt(nd)
				cycle
			endif
			
			if(.not.ic_done(0)) then
				!adjust the prior volumes & inc
				shrink=exp(-1.d0/dble(ic_npt(nd)))
				if(my_rank==0) then
	      				vprev=ic_vnow(nd)/shrink
					vnext=ic_vnow(nd)*shrink
					h=(vprev-vnext)/2.d0
		            		ic_inc(nd)=ic_hilike(nd)+log(ic_vnow(nd))-ic_Z(nd)
				
					!find lowlike
					indx=minloc(l(nd_j+1:nd_j+ic_npt(nd)))
					indx(1)=indx(1)+nd_j
		      			lowlike=l(indx(1))
					lowp(:)=p(:,indx(1))
					lowphyP(:)=phyP(:,indx(1))
					
					!set the limits
					do i3=1,ndims
						d4=ic_climits(nd,i3,2)-ic_climits(nd,i3,1)
						d1=ic_climits(nd,i3,1)+d4*lowp(i3)
						
						if( abs( d1 - ic_llimits(nd,i3,1) ) < d4 * 1d-5 ) then
							pt(1,1:ic_npt(nd))=ic_climits(nd,i3,1)+d4*p(i3,nd_j+1:nd_j+ic_npt(nd))
							ic_llimits(nd,i3,1)=minval(pt(1,1:ic_npt(nd)),MASK=pt(1,1:ic_npt(nd))>d1)
						elseif( abs( d1 - ic_llimits(nd,i3,2) ) < d4 * 1d-5 ) then
							pt(1,1:ic_npt(nd))=ic_climits(nd,i3,1)+d4*p(i3,nd_j+1:nd_j+ic_npt(nd))
							ic_llimits(nd,i3,2)=maxval(pt(1,1:ic_npt(nd)),MASK=pt(1,1:ic_npt(nd))<d1)
						endif
					enddo
					if(multimodal) then
						do i3=1,nCdims
							d4=ic_plimits(nd,i3,2)-ic_plimits(nd,i3,1)
							if( abs( lowPhyP(i3) - ic_plimits(nd,i3,1) ) < d4 * 1d-5 ) then
								ic_plimits(nd,i3,1)=minval(phyP(i3,nd_j+1:nd_j+ic_npt(nd)), &
								MASK=phyP(i3,nd_j+1:nd_j+ic_npt(nd))>lowPhyP(i3))
							elseif( abs( lowPhyP(i3) - ic_plimits(nd,i3,2) ) < d4 * 1d-5 ) then
								ic_plimits(nd,i3,2)=maxval(phyP(i3,nd_j+1:nd_j+ic_npt(nd)), &
								MASK=phyP(i3,nd_j+1:nd_j+ic_npt(nd))<lowPhyP(i3))
							endif
						enddo
					endif
				endif
			
				!now find a new point inside the hard constraint
				if(.not.eswitch) then
	            			acpt=.false.
					do
						if(my_rank==0) remFlag=remain(nd)
#ifdef MPI
						call MPI_BARRIER(MPI_COMM_WORLD,errcode)
						call MPI_BCAST(remFlag,1,MPI_LOGICAL,0,MPI_COMM_WORLD,errcode)
#endif
	                  			if(.not.remFlag) then
	            					!generate mpi_nthreads potential points
							call samp(pnew,phyPnew,lnew,sc_mean(1,:),d1,sc_tMat(1,:,:),ic_climits(nd,:,:),loglike,eswitch,lowlike,n,context)
						
							if(my_rank==0) then
								lnewa(nd,1)=lnew
							
								if(lnew>logZero) then
									pnewa(nd,1,:)=pnew(:)
									phyPnewa(nd,1,:)=phyPnew(:)
								endif
							endif
#ifdef MPI
							call MPI_BARRIER(MPI_COMM_WORLD,errcode)
							!now send the points to the root node
							if(my_rank/=0) then
								call MPI_SEND(lnew,1,MPI_DOUBLE_PRECISION,0,my_rank,MPI_COMM_WORLD,errcode)
								if(lnew>logZero) then
									call MPI_SEND(pnew,ndims,MPI_DOUBLE_PRECISION,0,my_rank,MPI_COMM_WORLD,errcode)
									call MPI_SEND(phyPnew,totPar,MPI_DOUBLE_PRECISION,0,my_rank,MPI_COMM_WORLD,errcode)
								endif
							else
								do i2=1,mpi_nthreads-1
									call MPI_RECV(lnewa(nd,i2+1),1,MPI_DOUBLE_PRECISION,i2,i2,MPI_COMM_WORLD,mpi_status,errcode)
									if(lnewa(nd,i2+1)>logZero) then
										call MPI_RECV(pnewa(nd,i2+1,1:ndims),ndims,MPI_DOUBLE_PRECISION,i2,i2,MPI_COMM_WORLD,mpi_status,errcode)
										call MPI_RECV(phyPnewa(nd,i2+1,1:totPar),totPar,MPI_DOUBLE_PRECISION,i2,i2,MPI_COMM_WORLD,mpi_status,errcode)
									endif
								enddo
							endif
							call MPI_BARRIER(MPI_COMM_WORLD,errcode)
#endif
							
							if( IS .and. my_rank == 0 ) then
								call ExtendArrayIfRequired(IS_counter(5),globff+1-IS_counter(5),IS_counter(4),IS_nstore/10,5,IS_iterinfo)
								call ExtendArrayIfRequired(IS_counter(1),mpi_nthreads,IS_counter(3),IS_nstore,ndims+6,IS_allpts)
								
								IS_counter(5) = globff+1
								
								i3 = IS_counter(1)	!no. of points in IS_allpts array
								do i2 = 1, mpi_nthreads
									if( lnewa(nd,i2)>logZero ) then
										IS_counter(1) = IS_counter(1)+1					!total no. of points collected so far
										IS_allpts(IS_counter(1),1:ndims) = pnewa(nd,i2,1:ndims)		!point
										IS_allpts(IS_counter(1),ndims+1) = lnewa(nd,i2)			!likelihood
									else
										if( IS ) IS_iterinfo(globff+1,4) = IS_iterinfo(globff+1,4) + 1
									endif
								enddo
								IS_allpts(i3+1:IS_counter(1),ndims+3) = 1d0					!check this point for ellipsoid membership in later iterations
								IS_allpts(i3+1:IS_counter(1),ndims+6) = dble(nd)				!node
								
								IS_iterinfo(globff+1,1) = 1d0							!current volume, V_{i}
								IS_iterinfo(globff+1,2:4) = IS_iterinfo(globff+1,2:4) + IS_counter(1)-i3	!no. of points collected at current iteration, n_{i}
								IS_iterinfo(globff+1,5) = dble(nd)						!current node
							endif
						endif
					
						if(my_rank==0) then
							!now check if any of them is a good one
	                        			do i=rIdx(nd),mpi_nthreads
	                        				numlike=numlike+1
	                        				if(lnewa(nd,i)>lowlike) then
									pnew(:)=pnewa(nd,i,:)
	                                    				phyPnew(:)=phyPnewa(nd,i,:)
	                                    				lnew=lnewa(nd,i)
	                                    				acpt=.true.
	                                    				!leftover points?
	                                    				if(i==mpi_nthreads) then
	                                    					remain(nd)=.false.
	                                          				rIdx(nd)=1
									else
	                                    					remain(nd)=.true.
	                                    					rIdx(nd)=i+1
									endif
	                                    				exit
								endif
	                              				if(i==mpi_nthreads) then
	                                    				remain(nd)=.false.
									rIdx(nd)=1
								endif
							enddo
						endif
					
#ifdef MPI
						call MPI_BARRIER(MPI_COMM_WORLD,errcode)
						call MPI_BCAST(acpt,1,MPI_LOGICAL,0,MPI_COMM_WORLD,errcode)
#endif
	                        	
						if(acpt) exit
					enddo
				
	                        	if(my_rank==0) then      	
	                  			p(:,indx(1))=pnew(:)
	                  			phyP(:,indx(1))=phyPnew(:)
	                  			l(indx(1))=lnew
						if(lnew>ic_hilike(nd)) ic_hilike(nd)=lnew
						
						!set the limits
						if(multimodal) then
							call setLimits(multimodal,ndims,nCdims,ic_llimits(nd,:,:), &
							ic_plimits(nd,:,:),pnew,phyPnew(1:nCdims),ic_climits(nd,:,:))
						else
							call setLimits(multimodal,ndims,nCdims,ic_llimits(nd,:,:), &
							ic_climits(nd,:,:),pnew,phyPnew(1:nCdims),ic_climits(nd,:,:))
						endif
						
						
						!calculate p(\Theta) n = \Sum_{i}^{niter} n_{i} E_{i}(\Theta) / V_{i}
						if( IS ) then
							n1 = 0; n2 = 0
							if( IS_counter(5) == globff+1 ) then
								n1 = int(IS_iterinfo(globff+1,2))	!n_{i}
								n2 = int(IS_iterinfo(globff+1,4))	!n_{i}
							endif
							
							if( n1 > 0 .and. n2 > 0 ) then
								!(p(\Theta) n1) for points collected in current iteration from current iteration
								IS_allpts(IS_counter(1)-n1+1:IS_counter(1),ndims+2) = dble(n2)
							
								!(p(\Theta) n1) for points collected in previous iterations from current iteration
								IS_allpts(1:IS_counter(1)-n1,ndims+2) = IS_allpts(1:IS_counter(1)-n1,ndims+2) + dble(n2)
							
								!(p(\Theta) n1) for points collected in current iteration from previous iterations
								do i = 1, globff
									IS_allpts(IS_counter(1)-n1+1:IS_counter(1),ndims+2) = IS_allpts(IS_counter(1)-n1+1:IS_counter(1),ndims+2) + IS_iterinfo(i,4)
								enddo
							endif
						endif
					endif
	            		else
					num_old=numlike
					
					acpt=.false.
					do
						if(my_rank==0) remFlag=remain(nd)
#ifdef MPI
						call MPI_BARRIER(MPI_COMM_WORLD,errcode)
						call MPI_BCAST(remFlag,1,MPI_LOGICAL,0,MPI_COMM_WORLD,errcode)
#endif
						if(.not.remFlag) then
							!first pick an ellipsoid according to the vol
							do
								!pick a sub-cluster according to the vol
								call selectEll(ic_sc(nd),sc_vol(nd_i+1:nd_i+ic_sc(nd)),i)
								i=i+nd_i
								if(sc_kfac(i)==0.d0 .or. sc_vol(i)==0.d0) then
									cycle
								else
									exit
								endif
							enddo
						
							!generate mpi_nthreads potential points
							d1=sc_kfac(i)*sc_eff(i)
							call samp(pnew,phyPnew,lnew,sc_mean(i,:),d1,sc_tMat(i,:,:),ic_climits(nd,:,:),loglike,eswitch,lowlike,n,context)
							
							if(my_rank==0) then
								if( IS ) then
									call ExtendArrayIfRequired(IS_counter(5),globff+1-IS_counter(5),IS_counter(4),IS_nstore/10,5,IS_iterinfo)
									call ExtendArrayIfRequired(IS_counter(1),mpi_nthreads,IS_counter(3),IS_nstore,ndims+6,IS_allpts)
									IS_counter(5) = globff+1
								endif
								
								lnewa(nd,1)=lnew
						
								if(lnew>logZero) then
									sEll(nd,1)=i-nd_i
									pnewa(nd,1,:)=pnew(:)
									phyPnewa(nd,1,:)=phyPnew
								endif
								
								if( IS ) IS_iterinfo(globff+1,4) = IS_iterinfo(globff+1,4) + n
							endif
#ifdef MPI
							call MPI_BARRIER(MPI_COMM_WORLD,errcode)
							!now send the points to the root node
							if(my_rank/=0) then
								if( IS ) call MPI_SEND(n,1,MPI_INTEGER,0,my_rank,MPI_COMM_WORLD,errcode)
								call MPI_SEND(lnew,1,MPI_DOUBLE_PRECISION,0,my_rank,MPI_COMM_WORLD,errcode)
								if(lnew>logZero) then
									call MPI_SEND(i-nd_i,1,MPI_INTEGER,0,my_rank,MPI_COMM_WORLD,errcode)
									call MPI_SEND(pnew,ndims,MPI_DOUBLE_PRECISION,0,my_rank,MPI_COMM_WORLD,errcode)
									call MPI_SEND(phyPnew,totPar,MPI_DOUBLE_PRECISION,0,my_rank,MPI_COMM_WORLD,errcode)
								endif
							else
								do i2=1,mpi_nthreads-1
									if( IS ) then
										call MPI_RECV(n,1,MPI_INTEGER,i2,i2,MPI_COMM_WORLD,mpi_status,errcode)
										IS_iterinfo(globff+1,4) = IS_iterinfo(globff+1,4) + n
									endif
									call MPI_RECV(lnewa(nd,i2+1),1,MPI_DOUBLE_PRECISION,i2,i2,MPI_COMM_WORLD,mpi_status,errcode)
									if(lnewa(nd,i2+1)>logZero) then
										call MPI_RECV(sEll(nd,i2+1),1,MPI_INTEGER,i2,i2,MPI_COMM_WORLD,mpi_status,errcode)
										call MPI_RECV(pnewa(nd,i2+1,1:ndims),ndims,MPI_DOUBLE_PRECISION,i2,i2,MPI_COMM_WORLD,mpi_status,errcode)
										call MPI_RECV(phyPnewa(nd,i2+1,1:totPar),totPar,MPI_DOUBLE_PRECISION,i2,i2,MPI_COMM_WORLD,mpi_status,errcode)
									endif
								enddo
							endif
							call MPI_BARRIER(MPI_COMM_WORLD,errcode)
#endif
							if( IS .and. my_rank == 0 ) then
								i3 = IS_counter(1)
								do i2=1, mpi_nthreads
									if( lnewa(nd,i2) > logZero ) then
										IS_counter(1) = IS_counter(1)+1					!total no. of points collected so far
										
										!reverse the limits on this point such that it is in unit hypercube & store it
										call ApplyLimits(1, ic_climits(nd,:,:), pnewa(nd,i2,1:ndims), IS_allpts(IS_counter(1),1:ndims))
										IS_allpts(IS_counter(1),ndims+1) = lnewa(nd,i2)			!likelihood
										
										!check how many ellipsoids this point lies in
										j = 1
	            								do m = nd_i+1, nd_i+ic_sc(nd)
	            									if( m == sEll(nd,i2)+nd_i .or. sc_npt(m) == 0 ) then
												IS_allpts(IS_counter(1),ndims+4) = m
												call ScaleFactor(1, ndims , pnewa(nd,i2,:), sc_mean(m,:), sc_invcov(m,:,:), IS_allpts(IS_counter(1),ndims+5))
												cycle
											endif
											if( ptIn1Ell(ndims,pnewa(nd,i2,:),sc_mean(m,:),sc_invcov(m,:,:),sc_kfac(m)*sc_eff(m))) j=j+1
	            								enddo
										IS_allpts(IS_counter(1),ndims+2) = dble(j) / ( totvol(nd) / ic_volFac(nd) )
										IS_iterinfo(globff+1,3) = IS_iterinfo(globff+1,3) + 1d0/dble(j)	!effective no. of points collected at each iteration
									endif
								enddo
								
								IS_iterinfo(globff+1,1) = IS_V(nd)						!volume, V_{i}
								IS_iterinfo(globff+1,2) = IS_iterinfo(globff+1,2) + IS_counter(1)-i3		!total no. of points collected at each iteration
								IS_iterinfo(globff+1,5) = dble(nd)						!node
								IS_allpts(i3+1:IS_counter(1),ndims+3) = 1d0					!check this point for ellipsoid membership in later iterations					!node
								IS_allpts(i3+1:IS_counter(1),ndims+6) = dble(nd)				!node
							endif
						endif
					
						if(my_rank==0) then
							!check if any of them is inside the hard edge
	                        			do j1=rIdx(nd),mpi_nthreads
								numlike=numlike+1
								if(my_rank==0 .and. ceff) then
									if(.not.remain(nd)) ic_eff(nd,2)=ic_eff(nd,2)+1d0
								endif
								
								if(lnewa(nd,j1)>lowlike) then
	                              					if(sc_npt(sEll(nd,j1)+nd_i)>0) then
	            								!find the no. of ellipsoids, j, the new points lies in
										j=1
	                                   					acpt=.true.
	            								do m=nd_i+1,nd_i+ic_sc(nd)
	            									if(m==sEll(nd,j1)+nd_i .or. sc_npt(m)==0) cycle
											if(ptIn1Ell(ndims,pnewa(nd,j1,:),sc_mean(m,:), &
											sc_invcov(m,:,:),sc_kfac(m)*sc_eff(m))) j=j+1
	            								enddo
									endif
	                                    			
									if(acpt) then
										!accept it with probability 1/j						
										if(j>1) then
	            									urv=ranmarNS(0)
											if(urv<=(1.d0/dble(j))) then
	           	 									acpt=.true.
											else
												acpt=.false.
											endif
										else
											acpt=.true.
										endif
										
										!point accepted
										if(acpt) then
											if(my_rank==0 .and. ceff) then
												if(.not.remain(nd)) ic_eff(nd,1)=ic_eff(nd,1)+1d0
											endif
											!leftover points?
	                                    						if(j1==mpi_nthreads) then
	                                    							remain(nd)=.false.
	                                          						rIdx(nd)=1
											else
	                                    							remain(nd)=.true.
	                                    							rIdx(nd)=j1+1
											endif
	                                    						exit
										endif
									endif
								endif
	                              				
								if(j1==mpi_nthreads) then
	                                    				remain(nd)=.false.
									rIdx(nd)=1
								endif
							enddo
						endif
					
#ifdef MPI
						call MPI_BARRIER(MPI_COMM_WORLD,errcode)
						call MPI_BCAST(acpt,1,MPI_LOGICAL,0,MPI_COMM_WORLD,errcode)
#endif
					
	                        		if(acpt) then
							if(my_rank==0) then
								pnew=pnewa(nd,j1,:)
								phyPnew=phyPnewa(nd,j1,:)
								lnew=lnewa(nd,j1)
								i=sEll(nd,j1)+nd_i
							
								if(lnew>ic_hilike(nd)) ic_hilike(nd)=lnew
						
								!set the limits
								if(multimodal) then
									call setLimits(multimodal,ndims,nCdims,ic_llimits(nd,:,:), &
									ic_plimits(nd,:,:),pnew,phyPnew(1:nCdims),ic_climits(nd,:,:))
								else
									call setLimits(multimodal,ndims,nCdims,ic_llimits(nd,:,:), &
									ic_climits(nd,:,:),pnew,phyPnew(1:nCdims),ic_climits(nd,:,:))
								endif
	                              			endif
							exit
	            				endif
	      				enddo

				
					if(my_rank==0) then
						
						!calculate p(\Theta) n = \Sum_{i}^{niter} n_{i} E_{i}(\Theta) / V_{i}
						if( IS ) then
							n1 = 0; n2 = 0
							if( IS_counter(5) == globff+1 ) then
								n1 = int(IS_iterinfo(globff+1,2))	!n_{coll, i}
								n2 = int(IS_iterinfo(globff+1,4))	!n_{i}
							endif
							
							if( n1 > 0 .and. n2 > 0 ) then
								
								!(p(\Theta) n) for points collected in current iteration from current iteration
								IS_allpts(IS_counter(1)-n1+1:IS_counter(1),ndims+2) = IS_allpts(IS_counter(1)-n1+1:IS_counter(1),ndims+2) * dble(n2)
								
								!(p(\Theta) n) for points collected in previous iterations from current iteration
								i1 = 0
								do j = IS_counter(2), IS_counter(1)-n1
									if( IS_allpts(j,ndims+3) == 0d0 ) cycle
									
									if( multimodal .and. .not.isAncestor(int(IS_allpts(j,ndims+6)), nd, ic_fnode(1:nd)) ) then
										IS_allpts(j,ndims+3) = 0d0
										cycle
									endif
								
									!check if this point lies in current ellipsoidal decomposition
									m = int(IS_allpts(j,ndims+4))
									if( m > nd_i .and. m <= nd_i+ic_sc(nd) .and. IS_allpts(j,ndims+5) <= sc_kfac(m)*sc_eff(m) ) then
										IS_allpts(j,ndims+2) = IS_allpts(j,ndims+2) + dble(n2) / IS_iterinfo(globff+1,1)
										if( i1 == 0 ) IS_counter(2) = j
										i1 = j
									else
										IS_allpts(j,ndims+3) = 0d0
									endif
								enddo
								if( i1 == 0 ) IS_counter(2) = IS_counter(1)-n1+1
								
								!(p(\Theta) n) for points collected in current iteration from previous iterations
								do j = 1, globff
									if( multimodal .and. .not.isAncestor(int(IS_iterinfo(j,5)), nd, ic_fnode(1:nd)) ) cycle
									
									if( int(IS_iterinfo(j,4)) > 0 ) IS_allpts(IS_counter(1)-n1+1:IS_counter(1),ndims+2) = IS_allpts(IS_counter(1)-n1+1:IS_counter(1),ndims+2) + IS_iterinfo(j,4) / IS_iterinfo(j,1)
								enddo
							endif
						endif
				
						if(ceff) then
							if(ic_eff(nd,1)>0d0 .and. mod(int(ic_eff(nd,1)),10)==0) then
								d1=ic_eff(nd,1)/ic_eff(nd,2)
								ic_eff(nd,1)=0d0
								ic_eff(nd,2)=0d0
								
								if(1.2d0*d1<ef) then
									!ic_eff(nd,3)=ic_eff(nd,4)*(1d0+0.2d0*sqrt(1000d0/dble(ic_npt(nd))))
									ic_eff(nd,3)=max(ic_eff(nd,3),ic_eff(nd,4)*(1d0+0.2d0*sqrt(1000d0/dble(ic_npt(nd)))))
								elseif(d1>1.2d0*ef) then
									ic_eff(nd,3)=max(ef,ic_eff(nd,4)/(1d0+0.2d0*sqrt(1000d0/dble(ic_npt(nd)))))
									!ic_eff(nd,3)=ic_eff(nd,4)/(1d0+0.2d0*sqrt(1000d0/dble(ic_npt(nd))))
									ic_eff(nd,4)=ic_eff(nd,3)
								endif
							endif
						endif
						
						!find the sub-cluster in which the rejected point lies
						n1=nd_j	
						do q=nd_i+1,nd_i+ic_sc(nd)
							if(sc_npt(q)==0) cycle
							n1=n1+sc_npt(q)
							if(indx(1)<=n1) then
								n1=n1-sc_npt(q)
	            						exit
	            					endif
						enddo
	
						i1=nlive
		
						!remove the rejected point & its likelihood & update the no. of points in the cluster
						if(indx(1)<nlive) then
							p(:,indx(1):nd_j+ic_npt(nd)-1)=p(:,indx(1)+1:nd_j+ic_npt(nd))
							phyP(:,indx(1):nd_j+ic_npt(nd)-1)=phyP(:,indx(1)+1:nd_j+ic_npt(nd))
							l(indx(1):nd_j+ic_npt(nd)-1)=l(indx(1)+1:nd_j+ic_npt(nd))
						endif
						sc_npt(q)=sc_npt(q)-1
						
						if(sc_npt(q)==0) then
							sc_vol(q)=0.d0
							sc_kfac(q)=0.d0
							sc_eff(q)=1.d0
							if(ic_sc(nd)==1) ic_inc(nd)=log(tol)
						elseif(sc_vol(q)>0.d0 .and. sc_kfac(i)>0.d0 .and. (i/=q .or. (i==q .and. sc_npt(q)>0))) then
				
							!min vol this ellipsoid should occupy
							if(dino) then
								if(ceff) then
									d4=ic_eff(nd,3)
								else
									d4=ef
								endif
								
								d1=(ic_vnow(nd)*ic_volFac(nd)*shrink/d4)*(sc_npt(q))/dble(ic_npt(nd))
							else
								d1=tiny(1.d0)
							endif
	
							!now evolve the ellipsoid with the rejected point
							call evolveEll(0,sc_npt(q),ndims,lowp,p(:,n1+1:n1+sc_npt(q)),sc_mean(q,:), &
							sc_eval(q,:),sc_invcov(q,:,:),sc_kfac(q),sc_eff(q),sc_vol(q),d1)
						endif
	                  	
						n1=sum(sc_npt(1:i-1))
						
						if(i/=q .or. (i==q .and. sc_npt(q)>0)) then
							!min vol this ellipsoid should occupy
							if(dino) then
								if(ceff) then
									d4=ic_eff(nd,3)
								else
									d4=ef
								endif
								
								d1=(ic_vnow(nd)*ic_volFac(nd)*shrink/d4)*(dble(sc_npt(i))+1.d0)/dble(ic_npt(nd))
							else
								d1=tiny(1.d0)
							endif
	
							!now evolve the ellipsoid with the inserted point
							call evolveEll(1,sc_npt(i),ndims,pnew,p(:,n1+1:n1+sc_npt(i)),sc_mean(i,:), &
							sc_eval(i,:),sc_invcov(i,:,:),sc_kfac(i),sc_eff(i),sc_vol(i),d1)
						endif
					endif
				
#ifdef MPI
					call MPI_BARRIER(MPI_COMM_WORLD,errcode)
					call MPI_BCAST(q,1,MPI_INTEGER,0,MPI_COMM_WORLD,errcode)
					call MPI_BCAST(i,1,MPI_INTEGER,0,MPI_COMM_WORLD,errcode)
					call MPI_BCAST(sc_kfac(i),1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,errcode)
					call MPI_BCAST(sc_eff(i),1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,errcode)
					call MPI_BCAST(sc_vol(i),1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,errcode)
					call MPI_BCAST(sc_kfac(q),1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,errcode)
					call MPI_BCAST(sc_eff(q),1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,errcode)
					call MPI_BCAST(sc_vol(q),1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,errcode)
#endif
				
					!now evolve the rest of the sub-clusters
					!only if sub-clustering wasn't done in the current iteration
					do j=nd_i+1,nd_i+ic_sc(nd)
						!sub-clusters with rejected/inserted point already evolved
						if(j==i .or. j==q) cycle
						if(sc_eff(j)>1.d0 .and. sc_vol(j)>0.d0 .and. sc_kfac(i)>0.d0) then
							d1=sc_eff(j)
							sc_eff(j)=max(1.d0,sc_eff(j)*(shrink**(2.d0/dble(ndims))))
							d1=sc_eff(j)/d1
							sc_vol(j)=sc_vol(j)*(d1**(dble(ndims)/2.d0))
						endif						
					enddo
				
					if(my_rank==0) then
						!add the new point & its likelihood & update the no. of points 
						!in the cluster new point's index
						n1=n1+sc_npt(i)+1
						!make room
						p(:,n1+1:nd_j+ic_npt(nd))=p(:,n1:nd_j+ic_npt(nd)-1)
						phyP(:,n1+1:nd_j+ic_npt(nd))=phyP(:,n1:nd_j+ic_npt(nd)-1)
						l(n1+1:nd_j+ic_npt(nd))=l(n1:nd_j+ic_npt(nd)-1)
						!insert the new point
						p(:,n1)=pnew(:)
						phyP(:,n1)=phyPnew(:)
						l(n1)=lnew
						!increment no. of points in sub-cluster i
						sc_npt(i)=sc_npt(i)+1
					endif
				endif
			
				if(my_rank==0) then
					!update evidence, info, prior vol, sampling statsitics
					globff=globff+1
					sff=sff+1
					gZold=gZ
					ic_zold(nd)=ic_z(nd)
					d1=lowlike+log(h)
	      				gZ=LogSumExp(gZ,d1)
	      				ic_Z(nd)=LogSumExp(ic_Z(nd),d1)
!					ic_info(nd)=exp(d1-ic_Z(nd))*lowlike+exp(gZold-ic_Z(nd))* &
!					(ic_info(nd)+gZold)-ic_Z(nd)
					ginfo=ginfo*exp(gzold-gz)+exp(d1-gz)*lowlike
					ic_info(nd)=ic_info(nd)*exp(ic_zold(nd)-ic_z(nd))+exp(d1-ic_z(nd))*lowlike
		
					!data for ev.dat file
					j1=mod(sff-1,updInt)+1
					evData(j1,1:totPar)=lowphyP(1:totPar)
					evData(j1,totPar+1)=lowlike
					evData(j1,totPar+2)=log(h)
					evData(j1,totPar+3)=dble(nd)
					
					lowlike=minval(l(nd_j+1:nd_j+ic_npt(nd)))
				
					if(abs(lowlike-ic_hilike(nd))<= 0.0001 .or. (ic_inc(nd)<log(tol) .and. &
					globff-nlive>50) .or. ff==maxIter) then
						ic_done(nd)=.true.
							
						!check if all done
						ic_done(0)=.true.
						do i=1,ic_n
							if(.not.ic_done(i)) then
								ic_done(0)=.false.
								exit
							endif
						enddo
					endif
				endif
			endif
		
			if(my_rank==0) then
				if(sff>0 .and. (mod(sff,updInt)==0 .or. ic_done(0))) then
				
					if( .not.outfile ) then
						k=0
						if( sff > updInt ) then
							k=size(evDataAll)
							allocate( evDataTemp(k) )
							evDataTemp=evDataAll
							deallocate( evDataAll )
							allocate( evDataAll(k+j1*(totPar+3)) )
							evDataAll(1:k)=evDataTemp(1:k)
							deallocate( evDataTemp )
						else
							deallocate( evDataAll )
							allocate( evDataAll(j1*(totPar+3)) )
						endif
						do i=1,j1
    							evDataAll(k+1:k+totPar+3) = evData(i,1:totPar+3)
							k=k+totPar+3
						enddo
					else
						!write the evidence file
						funit1=u_ev
						fName1=evname
						open(unit=funit1,file=fName1,form='formatted',status='old', position='append')
	    					write(fmt,'(a,i5,a)')  '(',totPar+2,'E28.18,i5)'	
						do i=1,j1
	    						write(funit1,fmt) evData(i,1:totPar+2),int(evData(i,totPar+3))
						enddo
						!close the files
						close(funit1)
					
		                		!write the live file
						funit1=u_phys
						fName1=physname
						funit2=u_live
						fName2=livename
						open(unit=funit1,file=fName1,form='formatted',status='replace')
						open(unit=funit2,file=fName2,form='formatted',status='replace')
	                		
						write(fmt,'(a,i5,a)')  '(',totPar+1,'E28.18,i4)'
	                			write(fmt1,'(a,i5,a)')  '(',ndims+1,'E28.18)'
						k=0
						do i=1,ic_n
							do j=1,ic_npt(i)
								k=k+1
								write(funit1,fmt) phyP(1:totPar,k),l(k),i
								lPts(1:ndims) = ic_climits(i,1:ndims,1)+(ic_climits(i,1:ndims,2)-ic_climits(i,1:ndims,1))*p(1:ndims,k)
								write(funit2,fmt1) lPts(1:ndims),l(k)
							enddo
	                			enddo
						!close the files
						close(funit1)
						close(funit2)
						
                  	
	                  			!write the resume file
						funit1=u_resume
						fName1=resumename
						open(unit=funit1,file=fName1,form='formatted',status='replace')
	                	  	
						write(funit1,'(l2)')genLive
						write(funit1,'(4i12)')globff,numlike,ic_n,nlive
						write(funit1,'(2E28.18)')gZ,ginfo
						write(funit1,'(l2)')eswitch
					
	            				!write branching info
		            			do i=1,ic_n
	            					write(funit1,'(i4)')ic_nBrnch(i)
							if(ic_nBrnch(i)>0) then
								write(fmt,'(a,i5,a)')  '(',2*ic_nBrnch(i),'E28.18)'
								write(funit1,fmt)ic_brnch(i,1:ic_nBrnch(i),1), ic_brnch(i,1:ic_nBrnch(i),2)
							endif
						enddo
					
						!write the node info
						do i=1,ic_n
							write(funit1,'(2l2,2i6)')ic_done(i),ic_reme(i),ic_fNode(i), ic_npt(i)
							write(funit1,'(3E28.18)')ic_vnow(i),ic_Z(i),ic_info(i)
							if(ceff) write(funit1,'(1E28.18)')ic_eff(i,4)
						enddo
	                  			close(funit1)
						
						
						!write the IS files
						if( IS ) then
							!write all the points collected so far
							open(unit=u_IS(1), file=IS_Files(1), form='unformatted', access='sequential', status='old', position='append')
							do i = IS_counter(7)+1, IS_counter(1)
								write(u_IS(1))IS_allpts(i,1:ndims+1),int(IS_allpts(i,ndims+6))
							enddo
							IS_counter(7) = IS_counter(1)
							close(u_IS(1))
							
							!write all the points collected so far
							open(unit=u_IS(2), file=IS_Files(2), form='unformatted', access='sequential', status='replace')
							write(u_IS(2))IS_counter(1:2),IS_counter(5)
							do i = 1, IS_counter(1)
								write(u_IS(2))IS_allpts(i,ndims+2),int(IS_allpts(i,ndims+3))
							enddo
							close(u_IS(2))
							
							!write the iteration info
							open(unit=u_IS(3), file=IS_Files(3), form='unformatted', access='sequential',status='old', position='append')
							do i = IS_counter(6)+1, IS_counter(5)
								write(u_IS(3))IS_iterinfo(i,1),int(IS_iterinfo(i,2:5))
							enddo
							IS_counter(6) = IS_counter(5)
							close(u_IS(3))
						endif
					endif
					
					!check if the parameters are close the prior edges
					if( prior_warning .and. mod(ff,50)== 0 ) then
						flag = .false.
						k=0
						do i=1,ic_n
							if( ic_npt(i) == 0 .or. ic_done(i) ) cycle
								
							ic_mean(i,1:ndims) = 0d0
							ic_sigma(i,1:ndims) = 0d0
							
							do j=1,ic_npt(i)
								k=k+1
								lPts(1:ndims) = ic_climits(i,1:ndims,1)+(ic_climits(i,1:ndims,2)-ic_climits(i,1:ndims,1))*p(1:ndims,k)
								if( prior_warning .and. mod(ff,50)== 0 ) then
									ic_mean(i,1:ndims) = ic_mean(i,1:ndims) + lPts(1:ndims)
									ic_sigma(i,1:ndims) = ic_sigma(i,1:ndims) + lPts(1:ndims) * lPts(1:ndims)
								endif
							enddo
							
							ic_mean(i,1:ndims) = ic_mean(i,1:ndims) / dble(ic_npt(i))
							ic_sigma(i,1:ndims) = sqrt(max(0d0, ic_sigma(i,1:ndims) / dble(ic_npt(i)) + ic_mean(i,1:ndims) * ic_mean(i,1:ndims)))
							
							do j = 1, ndims
								if( ic_sigma(i,j) <= 0.05 .and. ( ic_sigma(i,j) <= 0.05 .or. ic_sigma(i,j) >= 0.95 ) ) then
									if( .not. flag ) then
										write(*,*)
										write(*,*)"MultiNest Warning!"
										flag = .true.
									endif
									write(*,*)"Parameter ", j, " of mode ", i, " is converging towards the edge of the prior."
								endif
							enddo
						enddo
					endif
					
					if(mod(sff,updInt*10)==0 .or. ic_done(0)) call pos_samp(Ztol,globff,broot,nlive,ndims,nCdims,totPar, &
					multimodal,outfile,gZ,ginfo,ic_n,ic_Z(1:ic_n),ic_info(1:ic_n),ic_reme(1:ic_n),ic_vnow(1:ic_n),ic_npt(1:ic_n), &
					ic_nBrnch(1:ic_n),ic_brnch(1:ic_n,:,1),phyP(:,1:nlive),l(1:nlive),evDataAll,IS,IS_Z,dumper,context)
				endif
			endif
			
			nd_i=nd_i+ic_sc(nd)
			
			!prior volume shrinks
			if(my_rank==0) then
				if(.not.ic_done(nd)) ic_vnow(nd)=ic_vnow(nd)*shrink
				nd_j=nd_j+ic_npt(nd)
			endif
		enddo
		
		if(my_rank==0) then
			!update the total volume
			if(eswitch) then
				k = 0
				do i = 1, ic_n
					if( ic_done(i) ) then
						totVol(i) = 0d0
					else
						d1 = totVol(i)
						totVol(i) = sum(sc_vol(k+1:k+ic_sc(i)))
						if( ceff ) ic_eff(i,4) = ic_vnow(i) * ic_volFac(i) / totVol(i)
						
						if( IS ) IS_V(i) = IS_V(i) * totVol(i) / d1
					endif
					k = k + ic_sc(i)
				enddo
			endif
			
			!calculate the global evidence & info
			if(mod(ff,50)==0) then
				
				if(fback) then
					if( IS ) then
						IS_Z = logZero
						do j = 1, IS_counter(1)
							d1 = IS_allpts(j,ndims+1)-log(IS_allpts(j,ndims+2))
							IS_Z = LogSumExp(IS_Z, d1)
						enddo
					endif
					
					call gfeedback(gZ,IS,IS_Z,numlike,globff,.false.)
				
					if(debug) then
						d1=0.d0
						d2=0.d0
						j=0
						do i=1,ic_n
							if(ic_done(i)) then
								j=j+ic_sc(i)
								cycle
							endif
							d1=d1+ic_vnow(i)*ic_volFac(i)
							d2=d2+sum(sc_vol(j+1:j+ic_sc(i)))
							j=j+ic_sc(i)
						enddo
						write(*,*)ic_n,sc_n,d2/d1,count,scount
						if(ceff) write(*,*)ic_eff(1:ic_n,4)
					endif
				endif
			endif
			
			!switch ellipsoidal sampling on if the sum of the volumes of all the
		      	!clusters/sub-cluster is less than 1
			if(.not.peswitch) then
            			!marginal acceptance rate
				if(.not.peswitch) mar_r=1.d0/dble(numlike-num_old)
            			
				if(mar_r<ef) then
					escount(1)=escount(1)+1
					if(escount(1)==5) then
						peswitch=.true.
		                               	eswitchff(1)=ff
						totVol=1.d0
						escount(1)=0
					endif
				else
					escount(1)=0
				endif
			endif
		endif
	enddo
	
	deallocate( eswitchff, escount, dmin, evData, ic_sc, ic_npt, ic_done, &
	ic_climits, ic_volFac, sc_mean, sc_tmat, sc_kfac, sc_eff, sc_vol, lowp, lowphyP, pnew, phyPnew )
	
  end subroutine clusteredNest
  
!----------------------------------------------------------------------
 
 logical function ellIntersect(ndim,mean1,mean2,eval1,evec1,ef1,ef2,inv_cov1,inv_cov2)
 	
      implicit none
      
      integer ndim !dimensionality
      double precision eval1(:), evec1(:,:), mean1(:), mean2(:), inv_cov1(:,:), inv_cov2(:,:)
      ! matrices for ellipsoid interaction detection
      double precision, allocatable :: matA(:,:), matB(:,:), matR(:,:), matAinvB(:,:)
      ! variables for calculating eigenvalues of N*N real non-sym matrix
      double precision, allocatable :: evalR(:), evalI(:), VL(:,:), VR(:,:), WORK(:), delMean(:)
      ! effective enlargement factor
      double precision ef1,ef2
      integer inf,k,i
      integer i1,i2,i3,i4,i5
      
      
      
	if( ef1 == 0.d0 .and. ef2 == 0.d0 ) then
      		ellIntersect=.false.
            	return
      	else if( ef1 == 0.d0 ) then
      		ellIntersect = ptIn1Ell(ndim, mean1, mean2, inv_cov2, ef2)
            	return
	else if( ef2==0.d0 ) then
      		ellIntersect = ptIn1Ell(ndim, mean2, mean1, inv_cov1, ef1)
            	return
	endif
	
	allocate( matA(ndim+1,ndim+1), matB(ndim+1,ndim+1), matR(ndim+1,ndim+1), matAinvB(ndim+1,ndim+1) )
	allocate( evalR(ndim+1),evalI(ndim+1),VL(ndim+1,ndim+1),VR(ndim+1,ndim+1), WORK(4*ndim+4), delMean(ndim) )
      
      	delMean(1:ndim) = mean1(1:ndim) - mean2(1:ndim)
      	matA = 0.d0

	do i = 1, ndim
		matA(i,i)=eval1(i)
		matB(i,1:ndim)=inv_cov2(i,1:ndim)
		matB(ndim+1,i)=sum(-delMean(:)*inv_cov2(:,i))
		matB(i,ndim+1)=matB(ndim+1,i)
      	enddo
      	
	matA(ndim+1,ndim+1)=-1.d0/(ef1)
      	matB(ndim+1,ndim+1)=sum(matB(ndim+1,1:ndim)*(-delMean(:)))-ef2
      	matR=0.d0
      	
	do i=1,ndim
		matR(i,1:ndim)=evec1(:,i)
	enddo
      
      	matR(ndim+1,ndim+1)=1.d0
      	matB=MATMUL(MATMUL(matR,matB),TRANSPOSE(matR))
      	matAinvB=MATMUL(matA,matB)
      	i1=ndim+1
      	i2=ndim+1
      	i3=ndim+1
      	i4=ndim+1
      	i5=4*ndim+4
      	call DGEEV('N','N',i1,matAinvB,i2,evalR,evalI,VL,i3,VR,i4,WORK,i5,INF)
      	k=0
      
      	do i=1,ndim+1
		if(evalI(i)/=0.d0) then
			ellIntersect=.true.
			deallocate( matA, matB, matR, matAinvB, evalR, evalI, VL,VR, WORK, delMean )
			return
		else if(evalR(i)<0.d0) then
            		k=k+1
		endif
      	enddo
      
      	if(k<2) then
		ellIntersect=.true.
		deallocate( matA, matB, matR, matAinvB, evalR, evalI, VL,VR, WORK, delMean )
		return
      	endif
            
      	ellIntersect=.false.
	deallocate( matA, matB, matR, matAinvB, evalR, evalI, VL,VR, WORK, delMean )
 end function ellIntersect

!---------------------------------------------------------------------- 
  !sample a point inside the given ellipsoid with log-likelihood>lboundary
  subroutine samp(pnew,phyPnew,lnew,mean,ekfac,TMat,limits,loglike,eswitch,lowlike,n,context)
	
	implicit none
	double precision lnew
    	double precision pnew(ndims),spnew(ndims),phyPnew(totPar),ekfac
    	double precision mean(ndims),TMat(ndims,ndims)
	double precision limits(ndims,2)
	double precision lowlike	!likelihood threshold
	integer n			!no. of points drawn
    	logical eswitch
    	integer id,i,context
    	logical mark_bad_ellipses_SMF !SMF
    	integer i_cyc_SMF, i_SMF !SMF
    
    	INTERFACE
    		!the likelihood function
    		subroutine loglike(Cube,n_dim,nPar,lnew,context_pass)
			integer n_dim,nPar,context_pass
			double precision lnew,Cube(nPar)
		end subroutine loglike
    	end INTERFACE
    	
	
	n = 0
    	id = my_rank
        i_cyc_SMF = 0 !SMF
        mark_bad_ellipses_SMF = .true. !SMF
    	
	do
		n = n+1
		
	    	if(.not.eswitch) then
			!generate a random point inside unit hypercube
			call getrandom(ndims,pnew(1:ndims),id)
			phyPnew(1:ndims)=pnew(1:ndims)
			lnew=lowlike
	    		call loglike(phyPnew,ndims,totPar,lnew,context)
	    	else
			!generate a point uniformly inside the given ellipsoid
			call genPtInEll(ndims,mean,ekfac,TMat,id,pnew(1:ndims))
			spnew(:)=limits(:,1)+(limits(:,2)-limits(:,1))*pnew(:)
			do i=1,ndims
				if(pWrap(i)) then
					call wraparound(spnew(i),spnew(i))
				endif
			enddo
!SMF			if(.not.inprior(ndims,spnew(1:ndims))) cycle
!begin SMF bugfix
			if(.not.inprior(ndims,spnew(1:ndims))) then
              
              ! count the number of times the code has cycled. we do 
              ! not want to let the code loop infinitely, which can
              ! happen if an ellipse overlaps minimally with the prior
              i_cyc_SMF = i_cyc_SMF + 1
              !print "(I2,1X,I4,5(1X,F9.6))", id, i_cyc_SMF, spnew(1), &
              !                               spnew(2), spnew(3), &
              !                               spnew(4), spnew(5)
              if (i_cyc_SMF .lt. 1000) then
                 
                 cycle
                 
              else
                 
                 if (mark_bad_ellipses_SMF) then
                    
                    ! the code is stuck in an infinite loop: let's try
                    ! to get it out. there is an ellipse that has been 
                    ! drawn which extends well beyond the prior edge.
                    ! DO NOT call the likelihood function: return a 
                    ! very poor likelihood!
                    print*, '1000 CYCLES REACHED: BAD ELLIPSE!'
                    print*, 'auto-returning loglike = -1.d90 for final point'
                    print "(I2,1X,I4,5(1X,F9.6))", id, i_cyc_SMF, spnew(1), &
                                                   spnew(2), spnew(3), &
                                                   spnew(4), spnew(5)
                    lnew = -1.d90
                    exit
                    
                 else
                    
                    ! abort: the code is stuck in an infinite loop.
                    ! dump out the info MN is using, including the 
                    ! ellipse parameters.
                    print*, '1000 CYCLES REACHED: ABORTING!'
                    print*, 'mean:'
                    print*, mean
                    print*, 'ekfac:'
                    print*, ekfac
                    print*, 'TMat:'
                    print*, TMat
                    print*, 'limits:'
                    print*, limits
                    call MPI_ABORT(MPI_COMM_WORLD, MPI_ERR_OTHER)
                    call MPI_FINALIZE(MPI_SUCCESS)
                    stop
                    
                 end if
                 
              end if
              
           end if
!end SMF bug fix
			phyPnew(1:ndims)=spnew(1:ndims)
			lnew=lowlike
	    		call loglike(phyPnew,ndims,totPar,lnew,context)
	      	endif
		if(lnew>logZero) exit
	enddo
        
  end subroutine samp
  
  !----------------------------------------------------------------------
  
  subroutine ApplyLimits(flag,limits,pt,transpt)
  	
	implicit none
	
	!input parameters
	integer flag				!0 => given point in unit hypercube, apply the limits
						!1 => point to be transformed to unit hypercube, reverse the limits
	double precision limits(ndims,2)	!current limits
	double precision pt(ndims)		!point
	
	!output variables
	double precision transpt(ndims)		!final result
	
	
	if( flag == 0 ) then
		!apply the limits to point in unit hypercube
		transpt(:)=(pt(:)-limits(:,1))/(limits(:,2)-limits(:,1))
	elseif( flag == 1 ) then
		!reverse the limits such that the point is in unit hypercube
		transpt(:)=limits(:,1)+(limits(:,2)-limits(:,1))*pt(:)
	else
		write(*,*)'Incorrect value of flag passed to ApplyLimits'
		stop
	endif
	
  end subroutine ApplyLimits
  
  !----------------------------------------------------------------------
 
  subroutine selectEll(n,volx,sEll)
 	
      implicit none
      
      integer n!total no. of ellipsoids
      double precision volx(n)!no. points in each ellipsoid
      integer sEll!answer, the ellipsoid to sample from
      double precision volTree(n),totvol,urv
      integer i
 
 	totvol=0.d0
	volTree=0.d0
	do i=1,n
      		volTree(i)=totvol+volx(i)
            	totvol=volTree(i)
	enddo
	volTree=volTree/totvol
	urv=ranmarNS(0)
	do i=1,n
		if(urv<=volTree(i)) exit
	enddo
	sEll=i
  end subroutine selectEll
  
!----------------------------------------------------------------------
   
   !provide fback to the user
  subroutine gfeedback(logZ,IS,IS_logZ,nlike,nacc,dswitch)
    
	implicit none
    	!input variables
    	double precision logZ !log-evidence
	logical IS !importance sampling being done?
	double precision IS_logZ !importance sampling log-evidence
    	integer nlike !no. of likelihood evaluations
    	integer nacc !no. of accepted samples
	logical dswitch !dynamic live points
    
    	write(*,'(a,F14.6)')	     'Acceptance Rate:                  ',dble(nacc)/dble(nlike)
	write(*,'(a,i14)')   	     'Replacements:                     ',nacc
	write(*,'(a,i14)')   	     'Total Samples:                    ',nlike
	write(*,'(a,F14.6)')	     'Nested Sampling ln(Z):            ',logZ
	if( IS ) write(*,'(a,F14.6)')'Nested Importance Sampling ln(Z): ',IS_logZ
	if(dswitch) write(*,'(a,i5)')'Total No. of Live Points:         ',nlive
    
  end subroutine gfeedback
  
!----------------------------------------------------------------------
 
 !isolate the modes
 logical function isolateModes2(npt,ndim,nCdim,pt,naux,aux,ic_n,ic_fnode,ic_npt,ic_reme,ic_chk,ic_vnow,reCluster,limits)
 	implicit none
	
	!input parameters
	integer npt !total no. of points
	integer ndim !dimensionality
	integer nCdim !clustering dimensionality
	double precision pt(nCdim,npt)
	integer naux !no. of dimensions of the aux array
	double precision aux(naux,npt)
	logical ic_chk(ic_n) !whether to check a node for mode separation
	double precision ic_vnow(ic_n) !prior volumes
	double precision limits(maxCls,nCdim,2) !physical parameter ranges
	
	!input/output parameters
	!mode info
	integer ic_n !input: no. of existing modes. output: updated no. of existing modes
	integer ic_fnode(maxCls),ic_npt(maxCls)
	logical ic_reme(maxCls)
	
	!output variables
	logical reCluster(maxCls)
	
	!work variables
	integer i,j,k,i1,j2,j3,j4,i2,i3,i4,i5,n,n1,n2,n_mode,npt_mode,m,nLost
	integer, allocatable :: order(:), nptx(:), nodex(:)
	double precision d1,d2,d4,ef0, ef1, ef2
	integer nN,sc_n
	logical, allocatable :: gList(:), lList(:), toBeChkd(:), overlapk(:,:)
	logical flag,intFlag
	double precision, allocatable :: ptk(:,:), auxk(:,:), ptx(:,:), auxx(:,:)
	double precision, allocatable :: mMean(:), lMean(:), mStdErr(:), lStdErr(:)
	double precision, allocatable :: mean1(:), mean2(:), mean1w(:), mean2w(:)
	double precision, allocatable :: eval1(:), evec1(:,:), invcov1(:,:), invcov2(:,:)
	logical, allocatable :: wrapEll(:), wrapDim(:,:,:)
	integer, allocatable :: wrapN(:)
	double precision, allocatable :: meanw(:,:),meank(:,:),evalk(:,:),eveck(:,:,:),invcovk(:,:,:),tmatk(:,:,:),kfack(:)
	
	
	allocate( order(nCdim), nptx(npt/(ndim+1)+1), nodex(npt/(ndim+1)+1) )
	allocate( gList(npt/(ndim+1)+1), lList(npt/(ndim+1)+1), toBeChkd(npt/(ndim+1)+1), overlapk(npt/(ndim+1)+1,npt/(ndim+1)+1) )
	allocate( ptk(nCdim,npt), auxk(naux,npt), ptx(nCdim,npt), auxx(naux,npt), mMean(nCdim), lMean(nCdim), mStdErr(nCdim), &
	lStdErr(nCdim), mean1(nCdim), mean2(nCdim), mean1w(nCdim), mean2w(nCdim), eval1(nCdim), evec1(nCdim,nCdim), &
	invcov1(nCdim,nCdim), invcov2(nCdim,nCdim) )
	
	nN=ic_n
	isolateModes2=.false.
	reCluster=.false.
	
	i1=0 !no. of points traversed
	sc_n=0 !no. of clusters created so far
	
	do i=1,nN
		!enlargement factor
		ef0 = max( ( ic_vnow(i) * 100d0 / 111d0 ) + ( 11d0 / 111d0 ) , 0.8d0 )
	
		i3=ic_n
		
		if(ic_npt(i)==0) cycle
		
		if(.not.ic_chk(i)) then
			do j=1,nCdim
				ptk(j,i1+1:i1+ic_npt(i))=pt(j,i1+1:i1+ic_npt(i))
			enddo
			auxk(1:naux,i1+1:i1+ic_npt(i))=aux(1:naux,i1+1:i1+ic_npt(i))
			nptx(sc_n+1)=ic_npt(i)
			nodex(sc_n+1)=i
			i1=i1+ic_npt(i)
			sc_n=sc_n+1
			cycle
		endif
		
		nLost=0
		overlapk=.true.
		gList=.false.
		
		!cluster using G-means
		do j=1,nCdim
			d4=limits(i,j,2)-limits(i,j,1)
			ptk(j,i1+1:i1+ic_npt(i))=pt(j,i1+1:i1+ic_npt(i))
			
			if(pWrap(j)) then
				do i2=i1+1,i1+ic_npt(i)
					!scale to unit hypercube
					d2=limits(i,j,1)+d4*ptk(j,i2)
					
					call wraparound(d2,d1)
					
					!scale back to the limits
					if(d1 /= d2) ptk(j,i2)=(d1-limits(i,j,1))/d4
				enddo
			endif
		enddo
		
		auxk(1:naux,i1+1:i1+ic_npt(i))=aux(1:naux,i1+1:i1+ic_npt(i))
		n1=max(nCdim+1,3) !min no. of points allowed in a cluster
		n2=ic_npt(i)/n1+1 !max no. of clusters possible
		
		call doGmeans(ptk(1:nCdim,i1+1:i1+ic_npt(i)),ic_npt(i),nCdim,k,nptx(sc_n+1:sc_n+n2), &
		naux,auxk(:,i1+1:i1+ic_npt(i)),n1,n2)

		allocate(wrapN(k),wrapEll(k),wrapDim(k,nCdim,2),meank(k,nCdim),evalk(k,nCdim), &
		eveck(k,nCdim,nCdim),invcovk(k,nCdim,nCdim),tmatk(k,nCdim,nCdim),kfack(k),meanw(k,nCdim))
		wrapEll=.false.
		wrapDim=.false.
		
		nodex(sc_n+1:sc_n+k)=i

		!enclose sub-clusters in bounding ellipsoids
		n=i1
		wrapN=1
		do j=1,k
			call CalcBEllInfo(nptx(sc_n+j),nCdim,ptk(1:nCdim,n+1:n+nptx(sc_n+j)),meank(j,:), &
			evalk(j,:),eveck(j,:,:),invcovk(j,:,:),tmatk(j,:,:),kfack(j),n1)
			
			if(mWrap) then
				ef1=kfack(j)*max(2d0,((1.d0+ef0*sqrt(50.d0/dble(nptx(sc_n+j))))**(1d0/nCdim)))
				call wrapEllCheck(nCdim,meank(j,:),tmatk(j,:,:),ef1,limits(i,:,:), &
				wrapEll(j),wrapDim(j,:,:))
				
				!transform the mean to hypercube
				call Scaled2Cube(nCdim,limits(j,:,:),meank(j,:),meanw(j,:))
				
				if(wrapEll(j)) then
					do i2=1,nCdim
						if(wrapDim(j,i2,1)) then
							wrapN(j)=wrapN(j)+1
							meanw(j,i2)=1d0+meanw(j,i2)
						elseif(wrapDim(j,i2,2)) then
							wrapN(j)=wrapN(j)+1
							meanw(j,i2)=meanw(j,i2)-1d0
						endif
					enddo
				endif
				
				!scale the mean
				call Cube2Scaled(nCdim,limits(j,:,:),meanw(j,:),meanw(j,:))
			endif

			n=n+nptx(sc_n+j)
		enddo
		
		!calculate the standard error, required for localization
		if(.not.aWrap) then
			mMean=0.d0
			mStdErr=0.d0
			do j=i1+1,i1+ic_npt(i)
				mMean(:)=mMean(:)+ptk(:,j)
				mStdErr(:)=mStdErr(:)+ptk(:,j)**2
			enddo
			mMean=mMean/dble(ic_npt(i))
			mStdErr=mStdErr/dble(ic_npt(i))
			mStdErr=sqrt(mStdErr-mMean**2)
		endif
		
		do
			npt_mode=0
			lList=.false.
			toBeChkd=.false.
			n_mode=0 !no. of clusters in the mode
		
			!find a starting sub-cluster
			do j=1,k
				if(gList(j)) then
					cycle
				else
					n_mode=1
					npt_mode=nptx(sc_n+j)
					lList(j)=.true.
					gList(j)=.true.
					toBeChkd(j)=.true.
					m=j
					exit
				endif
			enddo
		
			!didn't find a starting position?
			if(n_mode==0) exit
			
			do
				flag=.false.
				do j=1,k
					if(.not.toBeChkd(j)) cycle
					flag=.true.
					toBeChkd(j)=.false.
					mean1(:)=meank(j,:)
					eval1(:)=evalk(j,:)
					ef1=kfack(j)*((1.d0+ef0*sqrt(60.d0/dble(nptx(sc_n+j))))**(1d0/nCdim))
					ef1=kfack(j)*1.5d0
					invcov1(:,:)=invcovk(j,:,:)
					evec1(:,:)=eveck(j,:,:)
					exit
				enddo
				
				if(.not.flag) exit
						
				do n=1,k
					if(lList(n) .or. n==j .or. .not.overlapk(n,j)) cycle
					mean2(:)=meank(n,:)
					ef2=kfack(n)*((1.d0+ef0*sqrt(60.d0/dble(nptx(sc_n+n))))**(1d0/nCdim))
					ef2=kfack(n)*1.5d0
					invcov2(:,:)=invcovk(n,:,:)
					
					intFlag=.false.
					
					if(mWrap .and. (wrapEll(j) .or. wrapEll(n))) then
						do i2=0,2**(wrapN(j)-1)-1
							call returnOrder(wrapN(j)-1,i2,order(1:wrapN(j)-1))
							
							i4=0
							do i5=1,nCdim
								if(wrapDim(j,i5,1) .or. wrapDim(j,i5,2)) then
									i4=i4+1
									if(order(i4)==0) then
										mean1w(i5)=meanw(j,i5)
									else
										mean1w(i5)=mean1(i5)
									endif
								else
									mean1w(i5)=mean1(i5)
								endif
							enddo
							
							
							do j2=0,2**(wrapN(n)-1)-1
								call returnOrder(wrapN(n)-1,j2,order(1:wrapN(n)-1))
								
								j4=0
								do j3=1,nCdim
									if(wrapDim(n,j3,1) .or. wrapDim(n,j3,2)) then
										j4=j4+1
										if(order(j4)==0) then
											mean2w(j3)=meanw(n,j3)
										else
											mean2w(j3)=mean2(j3)
										endif
									else
										mean2w(j3)=mean2(j3)
									endif
								enddo
							
								if(ellIntersect(nCdim,mean1w,mean2w,eval1,evec1,ef1,ef2,invcov1,invcov2)) then
									intFlag=.true.
									exit
								endif
							enddo
							
							if(intFlag) exit
						enddo
					elseif(ellIntersect(nCdim,mean1,mean2,eval1,evec1,ef1,ef2,invcov1,invcov2)) then
						intFlag=.true.
					endif
					
					if(intFlag) then
						lList(n)=.true.
						gList(n)=.true.
						toBeChkd(n)=.true.
						n_mode=n_mode+1
						npt_mode=npt_mode+nptx(sc_n+n)
					else
						overlapk(n,j)=.false.
						overlapk(j,n)=.false.
					endif
				enddo
			enddo
			
			!found a candidate?
			if((n_mode<k .or. ic_reme(i)) .and. npt_mode>=2*(ndim+1) .and. ((ic_npt(i)-npt_mode)>=2*(ndim+1) &
			.or. (ic_reme(i) .and. ic_npt(i)-npt_mode==0))) then
				flag=.true.
			else
				flag=.false.
			endif

			!now check for localization in ndim parameters
			!calculate the standard error
			if(flag .and. n_mode<k .and. npt_mode<ic_npt(i) .and. .not.aWrap) then
				n=i1
				lMean=0.d0
				lStdErr=0.d0
				do j=sc_n+1,sc_n+k
					if(lList(j)) then
						do i2=n+1,n+nptx(j)
							lMean(:)=lMean(:)+ptk(:,i2)
							lStdErr(:)=lStdErr(:)+ptk(:,i2)**2
						enddo
					endif
					n=n+nptx(j)
				enddo
				lMean=lMean/dble(npt_mode)
				lStdErr=lStdErr/dble(npt_mode)
				lStdErr=sqrt(lStdErr-lMean**2)
				
				do j=1,nCdim
					if(lStdErr(j)<mStdErr(j)) exit
					if(j==nCdim) flag=.false.
				enddo
			endif
			
			if(flag) then
				nLost=nLost+npt_mode
				ic_n=ic_n+1
				if(ic_n>maxCls) then
					write(*,*)"ERROR: More modes found than allowed memory."
					write(*,*)"Increase maxmodes in the call to nestrun and run MultiNest again."
					write(*,*)"Aborting"
#ifdef MPI
					call MPI_ABORT(MPI_COMM_WORLD,errcode)
#endif
                        		stop
				endif
				ic_fnode(ic_n)=i
				ic_reme(ic_n)=.false.
				ic_reme(i)=.true.
				isolateModes2=.true.
				reCluster(i)=.true.
				reCluster(ic_n)=.true.
				do j=1,k
					if(lList(j)) nodex(sc_n+j)=ic_n
				enddo
			endif
		enddo
		
		!make the leftover points into a new node
		if(ic_reme(i) .and. ic_npt(i)-nLost>=ndim+1 .and. .false.) then
			ic_n=ic_n+1
			ic_fnode(ic_n)=i
			ic_reme(ic_n)=.false.
			reCluster(ic_n)=.true.
			do j=1,k
				if(nodex(sc_n+j)==i) nodex(sc_n+j)=ic_n
			enddo
		endif
		
		if(ic_n==i3) then
			do j=1,nCdim
				ptk(j,i1+1:i1+ic_npt(i))=pt(j,i1+1:i1+ic_npt(i))
			enddo
			auxk(1:naux,i1+1:i1+ic_npt(i))=aux(1:naux,i1+1:i1+ic_npt(i))
		endif
		
		i1=i1+ic_npt(i)
		sc_n=sc_n+k
		deallocate(wrapN,wrapEll,wrapDim,meank,evalk,eveck,invcovk,tmatk,kfack,meanw)
	enddo
	
	!if modes found then re-arrange everything
	if(isolateModes2) then
		i1=0 !no. of points re-arranged
		do i=1,ic_n
			ic_npt(i)=0
			k=0 !no. of points traversed
			do j=1,sc_n
				if(nodex(j)==i) then
					!arrange the points
					do i2=1,nCdim
						ptx(i2,i1+1:i1+nptx(j))=ptk(i2,k+1:k+nptx(j))
					enddo
					auxx(1:naux,i1+1:i1+nptx(j))=auxk(1:naux,k+1:k+nptx(j))
					i1=i1+nptx(j)
					ic_npt(i)=ic_npt(i)+nptx(j)
				endif
				k=k+nptx(j)
			enddo
		enddo
		pt=ptx
		aux=auxx
	endif
	
	deallocate( order, nptx, nodex )
	deallocate( gList, lList, toBeChkd, overlapk )
	deallocate( ptk, auxk, ptx, auxx, mMean, lMean, mStdErr, lStdErr, mean1, mean2, mean1w, &
	mean2w, eval1, evec1, invcov1, invcov2 )
 
 end function isolateModes2
  
!----------------------------------------------------------------------

 subroutine setLimits(multimodal,ndim,nCdim,llimits,plimits,pnew,phyPnew,climits)
 
 	implicit none
	
	!input variables
	logical multimodal !set clustering limits?
	integer ndim !dimensionality
	integer nCdim !clustering dimension
	double precision pnew(ndim) !new point
	double precision phyPnew(nCdim) !new physical point
	double precision climits(ndim,2) !current scaling limits
	
	!input/output variables
	double precision llimits(ndim,2) !current limits
	double precision plimits(nCdim,2) !current clustering limits
	
	!work variables
	integer i
	double precision pt(ndim)
	
	
	!first scale the point
	do i=1,ndim
		pt(i)=climits(i,1)+(climits(i,2)-climits(i,1))*pnew(i)
		
		!wraparound
		if(pWrap(i)) call wraparound(pt(i),pt(i))
	enddo
					
	do i=1,ndim
		if(pt(i)<llimits(i,1)) then
			llimits(i,1)=pt(i)
		elseif(pt(i)>llimits(i,2)) then
			llimits(i,2)=pt(i)
		endif
	enddo
	
	if(multimodal) then
		do i=1,nCdim
			if(phyPnew(i)<plimits(i,1)) then
				plimits(i,1)=phyPnew(i)
			elseif(phyPnew(i)>plimits(i,2)) then
				plimits(i,2)=phyPnew(i)
			endif
		enddo
	endif
 
 end subroutine setLimits
  
!----------------------------------------------------------------------

 subroutine wraparound(oPt,wPt)
 	
	implicit none
	double precision oPt !actual point
	double precision wPt !wrapped-around point
	
	wPt=oPt
	do
		if(wPt<0.d0 .or. wPt>1.d0) then
			wPt=wPt-floor(wPt)
		else
			exit
		endif
	enddo
	
 end subroutine wraparound
  
!----------------------------------------------------------------------

 subroutine wrapEllCheck(ndim,mean,TMat,ef,limits,wrapEll,wrapDim)
 	
	implicit none
	
	!input variables
	integer ndim !dimensionality
	double precision mean(ndim)
	double precision TMat (ndim,ndim) !transformation matrix
	double precision ef !enlargement factor
	double precision limits(ndim,2) !current scale limits
	
	!output variable
	logical wrapEll
	logical wrapDim(ndim,2) !which dimensions in which directions to wrap
	
	!work variable
	integer i,j,k
	double precision cubeEdge(ndim),sCubeEdge
	double precision pnewM(1,ndim),u(1,ndim)
	
	
	wrapEll=.false.
	wrapDim=.false.
	
	!loop over the principle directions
	do i=1,ndim
		!loop over the two edges
		do k=1,2
			u(1,:)=0d0
			if(k==1) then
				u(1,i)=1d0
			else
				u(1,i)=-1d0
			endif
				
			pnewM=MatMul(u,TMat)
			cubeEdge(:)=sqrt(ef)*pnewM(1,:)+mean(:)
		
			!loop over dimensions
			do j=1,ndim
				if(pWrap(j) .and. (.not.wrapDim(j,1) .or. .not.wrapDim(j,2))) then
					!transform back to unit hypercube
					sCubeEdge=limits(j,1)+(limits(j,2)-limits(j,1))*cubeEdge(j)
				
					if(sCubeEdge<0d0) then
						wrapDim(j,1)=.true.
						wrapEll=.true.
					endif
					if(sCubeEdge>1d0) then
						wrapDim(j,2)=.true.
						wrapEll=.true.
					endif
				endif
			enddo
		enddo
	enddo
	
	!sanity check
	!no wraparound if both edges are outside the hyper cube boundary
	if(wrapEll) then
		wrapEll=.false.
		do i=1,ndim
			if(wrapDim(i,1) .and. wrapDim(i,2)) then
				wrapDim(i,1)=.false.
				wrapDim(i,2)=.false.
			else
				wrapEll=.true.
			endif
		enddo
	endif
	
	
 end subroutine wrapEllCheck
  
!----------------------------------------------------------------------

 subroutine scaled2Cube(ndim,limits,sP,cP)
 
 	implicit none
	
	!input variables
	integer ndim
	double precision limits(ndim,2) !current limits
	double precision sP(ndim) !scaled point
	
	!output variable
	double precision cP(ndim) !point in unit hypercube
	
	!work variables
	integer i
	
	
	do i=1,ndim
		cP(i)=limits(i,1)+(limits(i,2)-limits(i,1))*sP(i)
	enddo
 
 end subroutine scaled2Cube
  
!----------------------------------------------------------------------

 subroutine Cube2Scaled(ndim,limits,sP,cP)
 
 	implicit none
	
	!input variables
	integer ndim
	double precision limits(ndim,2) !current limits
	double precision sP(ndim) !scaled point
	
	!output variable
	double precision cP(ndim) !point in unit hypercube
	
	!work variables
	integer i
	
	
	do i=1,ndim
		sP(i)=(cP(i)-limits(i,1))/(limits(i,2)-limits(i,1))
	enddo
 
 end subroutine Cube2Scaled
  
!----------------------------------------------------------------------

 subroutine returnOrder(nBits,n,order)
 
 	implicit none
	
	!input variables
	integer nBits !no. of bits
	integer n !number under consideration
	
	!output variable
	integer order(nBits) !list with the order
	
	!work variables
	integer i,j
	
	
	order=0
	
	do i=1,nBits
		j=2**(i-1)
		if(mod(n/j,2)==0) then
			order(i)=0
		else
			order(i)=1
		endif
	enddo
 
 end subroutine returnOrder
 
!----------------------------------------------------------------------
 
 subroutine ExtendArrayIfRequired(n1,n2,n3,n4,n5,array)
 
 	implicit none
	
	!input variables
	integer n1					!number of rows stored
	integer n2					!additional rows to be stored
	integer n3					!max rows that can be stored
	integer n4					!no. of rows to be added if required
	integer n5					!no. of columns
	
	!output variables
	double precision, allocatable :: array(:,:)	!array
	
	!work variables
	integer i
	double precision, allocatable :: temp(:,:)
	
	
	if( n2 == 0 ) return
	
	if( n1+n2 > n3 ) then
		allocate(temp(n1,n5))
		temp = array
		deallocate(array)
		allocate(array(n1+max(n2,n4),n5))
		array = 0d0
		do i = 1, n1
			array(i,:) = temp(i,:)
		enddo
		deallocate(temp)
		
		n3 = n1+max(n2,n4)
	endif
 
 end subroutine ExtendArrayIfRequired
  
!----------------------------------------------------------------------
 
 !check if node1 is ancestor of node2
 logical function isAncestor(node1, node2, fnode)
 
 	implicit none
	
	!input variables
	integer node1			!ancestor node to be checked
	integer node2			!child node to be checked
	integer fnode(node2)		!array with parent nodes
	
	!work variables
	integer i, n2
	
	isAncestor = .false.
	if( node1 > node2 ) return
	
	if( node1 == 1 ) then
		isAncestor = .true.
		return
	endif
	
	n2 = node2
	do
		if( n2 == node1 ) then
			isAncestor = .true.
			return
		else
			n2 = fnode(n2)
			if( node1 > n2 ) return
		endif
	enddo
 
 end function isAncestor
  
!----------------------------------------------------------------------

end module Nested
