! The wrapper for MultiNest
! Sep 2010
! Author: Farhan Feroz

module nestwrap
	
	use ParamDef
	use CalcLike
	use propose
	use nested
	
	implicit none
	
	!nested sampling parameters
	logical nest_mmodal !multiple modes expected?
	logical nest_ceff !run in constant efficiency mode?
	integer nest_nlive !no. of live points
	integer nest_nPar !tot no. of parameters to be saved along with the sampling parameters
	integer nest_seed !seed for nested sampler, -ve means take it from sys clock
	double precision nest_tol !evidence tolerance factor
	double precision nest_efr !sampling efficiency
	character*100 nest_root !root for saving posterior files
	integer nest_maxModes !max modes expected, for memory allocation
	logical nest_fb !feedback on the sampling progress?
	logical nest_resume !resume from a previous run?
	integer sdim !dimensionality
	integer, dimension(:), allocatable :: nest_pWrap
	integer nest_updInt !no. of iterations after which to update on the progress
	
	!retain Curparams throughout, for the code to be run in the only fast params mode
	Type(ParamSet), SAVE :: Curparams

contains

!-------------------------------------------------------------------------
	
	subroutine setup
		
		implicit none
 
		double precision cosparams(num_params),lnew
		Type(ParamSet) Trial
		integer i
		
		cosparams = Scales%center
		Trial%P = cosparams
		lnew = -GetLogLike(Trial)
		Curparams=Trial !retain Curparams throughout, for the code to be run in the only fast params mode
		write(*,*) 'LIKELIHOOD TEST', lnew, cosparams
		
		!set up the prior ranges & the dimensionality
		sdim=0
		do i=1,num_params_used
			if(Scales%Pmin(params_used(i))/=Scales%Pmax(params_used(i))) sdim=sdim+1
		enddo
		
		!set dimensionaity
		nest_nPar=num_params
		
		if (Feedback > 1) then
			write(*,*) 'Dimensionality of param space = ', sdim
			write(*,*) 'Number of derived parameters  = ', nest_nPar-sdim
			write(*,*) 'Location of ln(like) in live points file, column = ', nest_nPar+1
		endif
		
	end subroutine setup
	
!-------------------------------------------------------------------------

	! Wrapper around Likelihood Function
	
	subroutine getLogLikeNS(Cube,n_dim,nPar,lnew,context)
	
		implicit none
		
		integer i,j,context,n_dim,nPar
		double precision cosparams(num_params),Cube(nPar),lnew
		logical accept
		Type(ParamSet) Trial
		double precision logZero
		parameter(logZero=-huge(1.d0)*epsilon(1.d0))
		character(LEN=60) fmt
		accept=.false.
		
		cosparams = Scales%center
		j=0
		do i=1,num_params_used
			if(Scales%Pmin(params_used(i))/=Scales%Pmax(params_used(i))) then
		   		j=j+1
			  	Cube(j)=Scales%Pmin(params_used(i))+(Scales%Pmax(params_used(i))-Scales%Pmin(params_used(i)))*Cube(j)
		   		cosparams(params_used(i))=Cube(j)
	  		else
		   		cosparams(params_used(i))=Scales%Pmin(params_used(i))
	  		endif
		enddo
		
		Trial=Curparams
		Trial%P=cosparams
		Cube=Trial%P
		lnew=-GetLogLike(Trial)
		if(lnew<-1.d20) lnew=logZero
		call AcceptReject(accept,CurParams%Info,Trial%Info)


		call WriteParams_nest(Trial%P, nPar, Cube)
		!write(*,*) ""
		!write(*,*) "num_params_used", num_params, num_params_used, nPar
		!write(*,*) "CubeII:", Cube

		fmt = trim(numcat('(2E16.7,',num_params))//'E16.7)'

!		write(109, fmt) lnew, lnew, Trial%P	
	
	end subroutine getLogLikeNS
		  	
!-------------------------------------------------------------------------

	! dumper, called after every updInt*10 iterations
	
	subroutine dumper(nSamples, nlive, nPar, physLive, posterior, paramConstr, maxLogLike, logZ, logZerr, context)
	
		implicit none
	
		integer nSamples				! number of samples in posterior array
		integer nlive					! number of live points
		integer nPar					! number of parameters saved (physical plus derived)
                integer context                                 ! unused
		double precision, pointer :: physLive(:,:)	! array containing the last set of live points
		double precision, pointer :: posterior(:,:)	! array with the posterior distribution
		double precision, pointer :: paramConstr(:)	! array with mean, sigmas, maxlike & MAP parameters
		double precision maxLogLike			! max loglikelihood value
		double precision logZ				! log evidence
                double precision logZerr			! error log evidence
		
	end subroutine dumper

!-----*-----------------------------------------------------------------

subroutine nest_Sample
		
		implicit none
		
	  	integer nclusters,context !total number of clusters found
	  	integer maxNode !variables used by the posterior routine
	
	  	call setup

	  	allocate(nest_pWrap(sdim))
	  	nest_pWrap=0
		
		
		!Calling MultiNest
	  	call nestRun(nest_mmodal,nest_ceff,nest_nlive,nest_tol,nest_efr,sdim,nest_nPar, &
	  	2,nest_maxModes,nest_updInt,-1.d90,nest_root,nest_seed,nest_pWrap, &
	  	nest_fb,nest_resume,.true.,.true.,-1.d90,-1,getLogLikeNS,dumper,context)
	  	
	  	deallocate(nest_pWrap)

end subroutine nest_Sample
		  	
!-------------------------------------------------------------------------

end module nestwrap
