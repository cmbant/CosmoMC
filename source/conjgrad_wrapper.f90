! slb 10aug04

module ConjGradModule
 use ParamDef
 use CalcLike
 use Random
 implicit none
 private 
 real ftol, Bestfit_loglike
 public conjgrad_wrapper, Bestfit_loglike

contains

 subroutine conjgrad_wrapper(Params,ftol_in,exitstatus)
  Type(ParamSet) Params,CurParams
  integer  niter
  real, intent(in) :: ftol_in
  real     vect(num_params_used),scale,tol,frac,fun
  real     change(num_params_used),gar(num_params_used),h(num_params_used)
  real     dfun(num_params_used),dfuntemp(num_params_used)
  real     MaxLike,CurLike,delta(num_params_used)
  integer  i,its,exitstatus,v
  real, parameter :: eps=1.E-10,tiny=1.E-5
  integer ii

  MaxLike = LogZero
  CurLike = StartLike
  CurParams=Params

  ftol = ftol_in

  scale = 0.2 ! change parameters by about 20 percent of their estimated widths. A guess
  tol= 0.1   ! get step length good to 10 per cent accuracy
  frac = 1e-5  ! don't bother tuning distance in param space if it only changes vect by frac
  ! make frac really small, as would like ftol to determine when it stops iterating
  niter = 200 ! maximum no of iterations (but hopefully will converge before this)
  v = Feedback

 ! scale the params so they are all roughly the same order of magnitude
  do i=1, num_params_used
        ii=params_used(i)
        vect(i)=CurParams%P(ii)/Scales%PWidth(ii)
        delta(i)=0.2 ! Step size in scaled parameter to use for numerical diffn. Update this in dffn when needed.
  end do

  call conjgrad(vect,CurParams,delta,scale,tol,frac,niter,ftol,fun,change,gar,h,dfun,dfuntemp,exitstatus,its,v)

! scale the params so they are all roughly the same order of magnitude
  do i=1, num_params_used
        ii=params_used(i)
        CurParams%P(ii)=vect(i)*Scales%PWidth(ii)
  end do
  
  Bestfit_loglike= fun
  write(*,*) 'Found minimum value of -Log(Like) = ',fun
  call AcceptReject(exitstatus < 3,Params%Info,CurParams%Info)
  if (exitstatus.eq.0) then
    write (*,*) 'Incredible! Loglike changed by less than ',eps,' so must have converged'
        Params=CurParams
  elseif (exitstatus.eq.1) then
    write(*,*) 'Converged sucessfully (loglike changes by less than',ftol,')'
        write(*,*) ' after ',its,' iterations'
        Params=CurParams
  elseif (exitstatus.eq.2) then
        write(*,*) 'Gave up trying to improve minimum after ',niter,' iterations'
        Params=CurParams ! presumably is an improvement
!  elseif(exitstatus.eq.10) then
!       write(*,*) 'Exiting because we seem to be going up hill now'
  elseif(exitstatus.eq.11) then
        write(*,*) 'Error: downhill gradient seems to point uphill!'
  else
        write(*,*) 'ERROR! Unknown exitstatus from conjgrad_cosmomc'
  end if


 end subroutine conjgrad_wrapper


 subroutine OutOfBounds
    write(*,*) 'ERROR: Trial parameters excluded by prior or error in likelihood'
    write(*,*) 'Try starting further away from problem regions?'
    call DoStop
 end  subroutine OutOfBounds


 function ffn(vect,OtherParams, keepNew)
  use ParamDef
  use CalcLike
  implicit none
  logical, intent(in) :: KeepNew

  Type(ParamSet) P, OtherParams
  integer i,ii
  real ffn
  real vect(num_params_used)


! OtherParams contains all the other info about the cosmology
! The subroutines only know about vect, and they make changes to this.
! But they carry around OtherParams (never changing it) so it can be used here.

! scale the params so they are all roughly the same order of magnitude
  P = OtherParams

  do i=1, num_params_used
        ii=params_used(i)
        P%P(ii)=vect(i)*Scales%PWidth(ii)
  end do
  ffn = GetLogLike(P)
  call AcceptReject(keepNew,OtherParams%Info,P%Info)
  if (keepNew) OtherParams = P

  if (ffn.eq.LogZero) call OutOfBounds
 
 end function ffn


 subroutine dffn(vect,df,OtherParams,delta)
  use ParamDef
  implicit none

  Type(ParamSet) OtherParams
  real vect(num_params_used)
  real df(num_params_used)
  integer  i
  real     delta(num_params_used),funcvect
  real  dfplus,dfminus
  real b, cc, bb, xhat
!  external ffn

! find the gradient numerically
      if (Feedback>1) write(*,*) 'Finding the gradient numerically...'
      funcvect=ffn(vect,OtherParams,.true.)
      if(Feedback>2) write(*,*) 'function value at initial position is:',funcvect

      do i=1,num_params_used

            call step1d(i,vect,delta,OtherParams,funcvect,dfplus,dfminus) ! step either side
                do while ( (abs(dfplus)<(0.5*ftol)) .and. (abs(dfminus)<(0.5*ftol)) )
                  ! while both are too small then increase step size to avoid numerical problems
                  delta(i)=delta(i)*2
                  call step1d(i,vect,delta,OtherParams,funcvect,dfplus,dfminus)
                end do

        ! Assume is a quadratic f(x) = a + bx + cx^2
        ! and use this to find gradient, minimum and better step
                cc=1/ (2*delta(i)*delta(i)) * ( dfplus+dfminus )
                bb=(dfplus-dfminus)/(2*delta(i)) ! the local gradient
                b= bb - 2*vect(i)*cc 
                df(i)=bb       ! estimate of the local gradient
                xhat=-b/(2*cc)  ! estimate of the position of the minimum
 
                ! refine the step size further if we bracket the minimum
                if ( (xhat<(vect(i)+abs(delta(i)))) .and. (xhat>(vect(i)-abs(delta(i)))) ) then
                  ! probably want either dfplus or dfminus near ftol to be happy
                  do while (.not.( ( (abs(dfplus) < (2*ftol) ) .and. (abs(dfplus)>(0.5*ftol)) ) .or. &
                              ( (abs(dfminus) < (2*ftol) ) .and. (abs(dfminus)>(0.5*ftol)) ) ))
                    delta(i)=abs(  1/(2*cc) *( abs(bb) - sqrt( bb*bb +4*cc*ftol) )  )
                        call step1d(i,vect,delta,OtherParams,funcvect,dfplus,dfminus)
                        cc=1/ (2*delta(i)*delta(i)) * ( dfplus+dfminus )
                    bb=(dfplus-dfminus)/(2*delta(i)) 
                    df(i)=bb
                  end do
                  
        end if

        if (Feedback>1) write(*,*) 'gradient wrt ',i,'th element is: ',df(i)
      end do

 end subroutine dffn



 subroutine step1d(i,vect,delta,OtherParams,funcvect,dfplus,dfminus)
  use ParamDef
  implicit none

  Type(ParamSet) OtherParams
  real vect(num_params_used),vectplus(num_params_used),vectminus(num_params_used)
  integer  i
  real     funcvect,delta(num_params_used),funcvectplus,funcvectminus
  real  dfplus,dfminus
 

    vectplus=vect
    vectplus(i)=vectplus(i) + delta(i)
    funcvectplus=ffn(vectplus,OtherParams,.false.)

    if (funcvectplus.eq.LogZero) call OutOfBounds

    dfplus=funcvectplus-funcvect
                    
    vectminus=vect
    vectminus(i)=vectminus(i) - delta(i)
    funcvectminus=ffn(vectminus,OtherParams,.false.)            

    if (funcvectminus.eq.LogZero) call OutOfBounds

    dfminus=funcvectminus-funcvect

 end subroutine step1d



! -----*-----------------------------------------------------------------
! -----*-----------------------------------------------------------------
! -----*-----------------------------------------------------------------


! Modified version of conjgrad from lensent. Only a few small changes:
!  - minimal changes to convert to f90 syntax
!  - allowed OtherParams to be carried around everywhere so the other info can be stored
!  - also carry around delta, the numerical diff step size.
!  - don't bother carrying len around, instead use Paramdef in the declarations
!  - meaning of ftol is now whether fn values are within a difference ftol (instead of factor)
!  - added the various reqd subroutines from lensent: rmsvect,dotvects,copyvect
!  - removed calls to display
!  - changed references to `mass'
!  - switched to Polak-Ribiere constrained by Fletcher-Reeves (was Polak-Ribiere before)
!    (looks slightly better?) 

! -----*-----------------------------------------------------------------

      subroutine conjgrad(vect,OtherParams,delta,scale,tol,frac,niter,conjgrad_ftol,fun,  &
                           change,gar,h,dfun,dfuntemp,exitstatus,its,v)

! -----*-----------------------------------------------------------------

! conjgrad minimises fun with respect to vect, of length len
! Adapted from Numerical Recipes implementation of the Polak-Ribiere
! variant of the Fletcher Reeves version of the conjugate gradient
! algorithm (see Numerical Recipies in Fortran 2nd Ed. p413 onwards)
! as suggested by David Mackay in
! http://wol.ra.phy.cam.ac.uk/mackay/cc/macopt.html
! slb dec 1997

! The arguments (NB. slightly different from the NR arguments):
! vect    On input: the place to start the minimisation from
!         On output: the lowest place found so far
! len     length of vect
! scale   typical (rms) amount to consider changing vect by (eg. 0.2)
! tol     "keep trying to decide how long a step to take 
!         until the two guesses at a good step length are within a 
!         factor of tol of each other"
! frac    "don't bother quite how small the step length 
!         should be if it is only changing vect by an rms of frac"
! niter   the maximum number of iterations that will be performed
!         unless convergence is reached earlier
! ftol    deemed to have converged sucessfully if fun changes by 
!         less than a fraction ftol in one iteration
! fun     On input: anything/unset
!         On output: the value of ffn at vect
! change,gar,h,dfun,temp are workspaces of length len as far
!         as the calling program is concerned
! exitstatus is an integer, unset on input
!         On output: 0 if exiting because mod(df).lt.eps
!                    1 if exiting because f has changed by less than ftol
!                    2 if exiting because number of iterations.gt.niter
!                    10 if exiting because f has increased on an iteration
!                    11 if exiting because df points uphill
! its is an integer. On output = the number of minimisation iterations.
! v is an integer. If v.ne.0 then stout is verbose, else not

! Assumed functions and subroutines:
! 1. a function called ffn(vect) which takes the argument vect and 
!    returns the value of the function for that vect
! 2. a subroutine called dffn(vect,dfun), where dfun
!    is an output argument containing the gradient of the function wrt
!    each element of vect
!    (see also chkdffn(vect,dfun) which checks your gradient
!    function against the numerically calculated gradient)
! 3. a subroutine called display which has no effect on its arguments
!    but plots and reports on the progress of the minimisation.
!    NB, it is only ever called just after calling ffn.
!    subroutine display(vect,h,change,len,its,dfun,fun,exitstatus)
!      vect    the current position
!      h       the direction just tried
!      change  the actual change that was just made to the postion
!      len     as before
!      its     the current iteration number
!      dfun    the current function gradient
!      fun     the current function value
!      exitstatus exitstatus from display, if .gt.20 then exit from
!              conjgrad with current position and this exitstatus. 
!              (ignored if exitstatus.lt.20)
! 
! NB. if exitstatus.ne.0 then the last time ffn was called it was at
! the output position, and any subsequent calls to dffn or display
! were also at that position.
! if exitstatus.eq.0 then the last time ffn was called it was at the 
! output position, but dffn may have been called at other positions
! between when ffn was called and exiting.

! -----*-----------------------------------------------------------------

      use Paramdef
      implicit none
      Type(ParamSet) OtherParams
      integer  niter
      real     vect(num_params_used),scale,tol,frac,conjgrad_ftol,fun,delta(num_params_used)
      real     change(num_params_used),gar(num_params_used),h(num_params_used),dfun(num_params_used),dfuntemp(num_params_used)

      integer  i,its,exitstatus,v
      real     eps,tiny
      parameter (eps=1.E-10,tiny=1.E-5)
      real gg, dgg_pr,dgg_fr,gam,gam_pr,gam_fr,funprev,df1,moddfun
!      real ffn
!      external ffn

! initialise
      ftol=conjgrad_ftol ! make this value global to the module
      exitstatus=0.

! find initial value of fun, used later for checking for convergence
      funprev=ffn(vect,OtherParams,.true.)
      call dffn(vect,dfun,OtherParams,delta)

! check that dfun.dfun .ne. 0, if so, vect must be at the minimum already
      call rmsvect(dfun,num_params_used,moddfun)
!      write(*,*)'moddfun is: ',moddfun
      if (moddfun.lt.eps) then
        if (Feedback > 1) write(*,*) 'The modulus of the gradient is zero for the vector', &
           ' input to conjgrad. rms value of the gradient is: ',moddfun
        fun=funprev
        return
      end if
      df1=-moddfun**2   ! so it is the dot product of the gradient with the dirn

! initialise vectors
      do i=1,num_params_used
        gar(i) = -dfun(i)
        h(i) = gar(i)
      end do


! main loop
      do its=1,niter

!        call plotmap(h,16,1.0,'x','y','dirn')
! decide how much to change vect in the direction h, and do it
        call slinmin(vect,h,OtherParams,delta,scale,tol,frac,change,dfun,df1,  &
                       dfuntemp,v)
! on return dfun is the gradient at the new position
        fun=ffn(vect,OtherParams, .true.)
        if(v>1)write(*,*) ' for which F =',fun
!        call display(vect,h,change,num_params_used,its,dfun,fun,exitstatus)
! ... check exitstatus from display
!        if (exitstatus.gt.20) then
!          if(v>0)write(*,*) 'Exiting conjgrad because exitstatus from display .gt.20'
!          return
!        else
!          exitstatus=0. ! reset in case changed
!        end if

! test for convergence:
        call rmsvect(dfun,num_params_used,moddfun)
        if (moddfun.lt.eps) then
          write(*,*) 'modulus of the gradient is less than ',eps
          write(*,*) 'mod(gradient) = ',moddfun
          write(*,*) 'Sucessfully converged after ',its,' iterations'
          return
        end if

!        if(2.*abs(fun-funprev).le.ftol*(abs(fun)+abs(funprev)+eps)) then
        if(abs(fun-funprev).le.ftol) then
          write(*,'(a,f12.8,a,i6,a)')   &
       ' Change in F is less than ',ftol,' after ',  &
        its,' iterations'
          if(v>1) write(*,*) 'The modulus of the gradient is ',moddfun
          if(v>1) write(*,*) 'Returning to calling program.'
          exitstatus=1
          return
        end if

! test for sanity
        if (fun.gt.funprev*(1+tiny)) then
          write(*,*)'The new function value is larger than the previous'
          write(*,*)'one. Continuing anyway.'
          write(*,*)'You are moving to a new local minimum.'
          write(*,*)'Reducing the scale suggestion might prevent this'
        end if          

        funprev=fun ! record previous value of fun

! evaluate gam using eqn 10.6.7, Numerical Recipies in Fortran 2nd Ed.
        gg=0.
        dgg_pr=0.
        dgg_fr=0.
        do i=1,num_params_used
          gg=gg+gar(i)**2
          ! g_{i+1} is the previous gradient of the function, stored in dfun
          dgg_pr=dgg_pr+(dfun(i)+gar(i))*dfun(i)
          dgg_fr=dgg_fr+(dfun(i))*dfun(i)
        end do
        gam_pr=dgg_pr/gg
        gam_fr=dgg_fr/gg

! use condition from http://www.cs.nyu.edu/faculty/overton/software/cgqn/cgprfr.m                              
     if (gam_pr < -gam_fr) then     ! ensures beta <= |beta_fr|, allowing proof of
         gam = -gam_fr        ! global convergence, but avoids inefficiency
      elseif (gam_pr > gam_fr) then ! of FR which happens when beta_fr gets stuck near 1
         gam = gam_fr
      else
         gam = gam_pr
      end if

! work out the new vectors
        do i=1,num_params_used
          gar(i) = -dfun(i)
          h(i) = gar(i) + gam*h(i)
        end do

! just check that the suggested direction points downhill
        call dotvects(h,dfun,num_params_used,df1)
        if (df1.gt.0) then
          write(*,*) 'ERROR: suggested direction points uphill ',  &
                      'after ',its,' iterations'
          exitstatus=11
          return
        end if

      end do
      
      write(*,*) 'have done all ',niter,' iterations without converging'
      write(*,*) 'The modulus of the gradient is ',moddfun
      exitstatus=2

      end subroutine conjgrad

! -----*-----------------------------------------------------------------
! note slinmin updates the vect to the new one 
! and doesn't mess up h as the one in NR does (for no reason?)
! and returns the new value of the function, fun
! and the gradient at that point, dfun

      subroutine slinmin(vect,dirn,OtherParams,delta,scale,tol,frac,change,  &
                              dfun,df1,dfuntemp,v)
      use ParamDef
      implicit none
      Type(ParamSet) OtherParams
      integer  i,v
      real     vect(num_params_used),dirn(num_params_used),scale,tol,frac,change(num_params_used)
      real     dfun(num_params_used),dfuntemp(num_params_used)
      real     lambda1,lambda2,df1,df2,lambda,delta(num_params_used)

! dfun is the gradient at the current (soon to be old) position

! bracket the minimum ie output lambda1 and lambda2
      if (v>1) write(*,*)
      if (v>1) write(*,*)' Bracketing the minimum in the new direction...'
      call sbrac(vect,dirn,OtherParams,delta,scale,lambda1,lambda2,df1,df2,  &
                     change,dfun,dfuntemp,v)
! use change as a temporary array in the above subroutine
! dfuntemp is now the gradient at lambda2
! dfun is now the gradient at lambda1

! find the position of the minimum, ie. lambda
! dfun is the gradient at lambda1 since lambda1=0.0 at this point
      if (v>1)write(*,*) '   Locating the minimum in this direction...'
      call sloc(lambda1,lambda2,df1,df2,vect,dirn,OtherParams,delta,tol,frac,  &
                   lambda,change,dfun,dfuntemp,v)
! use change as a temporary array in the above subroutine
! dfun is updated with the gradient at lambda

! update the vect with vect for the miminum
      do i=1,num_params_used
        change(i) = lambda*dirn(i)
        vect(i) = vect(i) + change(i)
      end do

      end subroutine slinmin

! -----*-----------------------------------------------------------------

      subroutine sbrac(vect,dirn,OtherParams,delta,scale,lambda1,lambda2,  &
                          df1,df2,trialvect,dfun1,trialdfun2,v)
      use ParamDef
      implicit none
      Type(ParamSet) OtherParams
      integer  v
      real     scale,lambda1,lambda2,df1,df2,delta(num_params_used)
      real     vect(num_params_used),dirn(num_params_used),trialvect(num_params_used)
      real     dfun1(num_params_used),trialdfun2(num_params_used)
      real     eps
      parameter(eps=1E-6)

      real     moddirn,lambdaadd
      integer  its,itsmax,i
      parameter(itsmax=100)

! find the magnitude of dirn
      call rmsvect(dirn,num_params_used,moddirn)
!      write(*,*) 'modulus of direction to try is:',moddirn

! give the starting values for lambda:
!   dirn points downhill by construction (and see check), 
!   so lambda must be positive
      lambda1=0.0
!   we don't know about lambda2, but scale/mod(dirn) would make
!   the rms of the changes in each vect entry = scale
!   NB we wouldn't get this far if moddirn=0
      lambda2=scale/moddirn
      lambdaadd=lambda2

      call trialdfscalar(vect,OtherParams,delta,lambda2,dirn,df2,trialvect,trialdfun2)
      if (v>1) write(*,*)   &
        'For initial lower limit on step length (lambda1),',lambda1,  &
        ', gradient.dirn is ',df1
      if (v>1) write(*,*) 'For initial upper estimate (lambda2),',  &
             lambda2,', gradient.dirn is ',df2
      if (abs(df2).lt.eps) return

      its=1
      do while (df2.le.0.0)
        lambda1=lambda2
        do i=1,num_params_used
          dfun1(i)=trialdfun2(i)
        end do
        df1=df2
        lambda2=lambda2+lambdaadd
        call trialdfscalar(vect,OtherParams,delta,lambda2,dirn,df2,  &
                                             trialvect,trialdfun2)

!         write (13, '(3I6, 3F12.5, F22.10)') i, j, k,
!      write (*,'(A12,F12.5)') '    gam:',gam           
!        write(*,'(a,i2,a,f,a,f)') 'After ',its,' iterations, lambda2 is , '
        if (v>1) write(*,*) 'After ',its,' iterations, lambda2 is'  &
                  ,lambda2,', where gradient is ',df2
        if (v>1) write(*,*) '              and now lambda1 is'  &
                  ,lambda1,', where gradient is ',df1
        its=its+1
        if (its.ge.itsmax) then
          if (v>1)write(*,*)   &
           'Warning: have not bracketed the minimum after',  &
           its,' increments. ',  &
           'Increase scale?'
        end if
      end do

      end subroutine sbrac

! -----*-----------------------------------------------------------------

      subroutine trialdfscalar(vect,OtherParams,delta,lambda,dirn,dfscalar,   &
                                trialvect,trialdf)
      use ParamDef
      implicit none
      Type(ParamSet) OtherParams
      real     vect(num_params_used),lambda,dirn(num_params_used),dfscalar
!      real     sumdirn
      real delta(num_params_used)
      real     trialvect(num_params_used),trialdf(num_params_used)
      integer i

! make up the trial position in vector space
      do i=1,num_params_used
        trialvect(i)=vect(i)+lambda*dirn(i)
      end do

! find the direction of the gradient there
      call dffn(trialvect,trialdf,OtherParams,delta)

! find component of the gradient in the direction we are investigating
      call dotvects(dirn,trialdf,num_params_used,dfscalar)

      end subroutine trialdfscalar

! -----*-----------------------------------------------------------------

      subroutine sloc(lambda1,lambda2,df1,df2,vect,dirn,OtherParams,delta,tol,frac,  &
                            lambda,trialvect,trialdfun1,trialdfun2,v)
! locate a good value for the step length, lambda, between the limits
! already found, lambda1 and lambda2
      use ParamDef
      implicit none
      Type(ParamSet) OtherParams
      integer  v
      real     lambda1,lambda2,lambda,frac,tol,delta(num_params_used)
      real     vect(num_params_used),dirn(num_params_used),trialvect(num_params_used),trialdfun1(num_params_used)
      real     trialdfun2(num_params_used)

      integer  nbis,nbismax
      parameter(nbismax=30)
      real     lambdabis,df1,df2,dfbis
      real     rmsdirn,rmsv
      real     eps
      parameter(eps=1E-6)   ! machine precision     

! trialdfun1 is ALWAYS the gradient at lambda1
! on entry (ie. now) trialdfun2 is the gradient at lambda2

! this is left as a possibility by sloc:
      if (abs(df2).lt.eps) then
        lambda=lambda2
        call copyvect(trialdfun2,num_params_used,trialdfun1)
        write(*,*) 'larger bracket hits the minimum in this direction'
        return
      end if

! report status:
      call rmsvect(dirn,num_params_used,rmsdirn)
      call rmsvect(vect,num_params_used,rmsv)
!      write(*,*) 'initially: lambda1= ',lambda1
!      write(*,*) 'initially: lambda2, df2= ',lambda2,df2


! do bisections to home in on a good step length lambda:
      nbis=1
      do while ((abs(lambda1-lambda2).gt.(tol*(lambda1+lambda2))).and.  &
                ((lambda1+lambda2)*rmsdirn.gt.frac*rmsv))
! ... while lambda1 and lambda2 differ by a fraction tol
! and lambda1 and lambda2 are not both too very small
        lambdabis=(lambda1+lambda2)/2  ! bisect
        call trialdfscalar(vect,OtherParams,delta,lambdabis,dirn,dfbis,  &
                                                trialvect,trialdfun2)
! trialdfun1 is the gradient at lambda1
! trialdfun2 is the gradient at lambdabis

! replace one bound with the new one such that the minimum is still bracketed
        if (dfbis.lt.0) then
          lambda1=lambdabis
          df1=dfbis
          if(v>1)write(*,*) 'after ',nbis,' bisections, lambda1,df1:',  &
                   lambda1,df1
          call copyvect(trialdfun2,num_params_used,trialdfun1)
! trialdfun1 is the gradient at lambda1 
! trialdfun2 is the gradient at lambda1
        else 
          lambda2=lambdabis
          df2=dfbis
          if(v>1)write(*,*) 'after ',nbis,' bisections, lambda2,df2:',  &
                   lambda2,df2
! trialdfun1 is the gradient at lambda1
! trialdfun2 is the gradient at lambda2
        end if

! increment counter and warn if too many its:
        nbis=nbis+1
        if (nbis.ge.nbismax) then
          write(*,*) 'Warning: no. of bisections to find fmin is'  &
                     ,nbis
        end if

      end do


! report status
      if ((abs(lambda1-lambda2).le.(tol*(lambda1+lambda2)))) then
        if(v>1)write(*,*)   &
           'finishing because lambda1 and lambda2 are the same ',  &
              'to within a fraction tol=',tol
!        write(*,*) 'abs(lambda1-lambda2)  =',abs(lambda1-lambda2)
!        write(*,*) 'tol*(lambda1+lambda2))=',tol*(lambda1+lambda2)
      end if

      if ((lambda1+lambda2)*rmsdirn.le.frac*rmsv) then
        if(v>1)write(*,*) 'not bothering to refine lambda any further because '
        if(v>1)write(*,*) 'it will only change the normalised parameters by a fraction ',frac
        if(v>1)write(*,*)'(lambda1+lambda2)*rmsdirn=',(lambda1+lambda2)*rmsdirn
        if(v>1)write(*,*)'frac*rmsv             =',frac*rmsv
      end if

      lambda=lambda1   ! err on the safe side
      if(v>1)write(*,*) 'using lambda = ',lambda

      end subroutine sloc

! -----*-----------------------------------------------------------------

! -----*-----------------------------------------------------------------
      subroutine rmsvect(vect,len,rms)
      implicit none
      integer len,i
      real    vect(len),rms,sumsq
      sumsq=0.
      do i=1,len
          sumsq=sumsq + vect(i)**2
      end do
      rms=sqrt(sumsq)
      end subroutine rmsvect
! -----*-----------------------------------------------------------------
! find the dot product of two vectors
      subroutine dotvects(vect1,vect2,len,dotprod)
      implicit none
      integer len,i
      real    vect1(len),vect2(len),dotprod
      dotprod=0.0
      do i=1,len
          dotprod=dotprod + vect1(i)*vect2(i)
      end do
      end subroutine dotvects
! -----*-----------------------------------------------------------------
! in f90:  vectout=vectin
      subroutine copyvect(vectin,n,vectout)
      implicit none
      integer  n,i
      real     vectin(n),vectout(n)
      do i=1,n
        vectout(i)=vectin(i)
      end do
      end subroutine copyvect
! -----*-----------------------------------------------------------------



 
end module ConjGradModule
