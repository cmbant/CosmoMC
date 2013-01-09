!Do Metropolis-Hastings and slice sampling algorithms
!Also directional gridding method for fast-slow parameters
!Provisional implementation of Wang-Landau and Multicanonical methods
!August 2006

module MonteCarlo
 use ParamDef
 use CalcLike
 use Random
 use propose
 use IO
 implicit none

 integer :: indep_sample = 0
  !number of iterations between dumping full model data. If zero then skip.


 integer :: directional_grid_steps = 20 
  !for sampling_method = sampling_slowgrid, number of steps per grid

 real(mcp), parameter    :: eps_1 = 1.00001
 real(mcp) :: MaxLike

 logical :: fast_slicing = .false. !slice fast parameters when using Metropolis 
 logical :: slice_stepout = .true.
 logical :: slice_randomsize = .false.

!Multicanonical/W-L sampling parameters
  real(mcp), parameter :: mc_logspace =  6. !1/ steps to use in log-likelihood * number of parameter
  integer :: mc_mini = 0, mc_maxi = 0
  integer :: mc_update_steps = 7000, mc_update_burn = 50
  integer :: mc_steps_inc = 1000
  real(mcp), dimension(:,:), allocatable :: mc_like_counts
  real(mcp), dimension(:), allocatable :: mc_lnweights
!Wang-Landau parameter
  real(mcp) :: WL_f = 1, WL_min = 0.1, WL_minW
  real(mcp) :: WL_maxL=8. !WL_maxL is log like from best to use, in units of number of parameters
  integer :: WL_update_steps = 4000
  real(mcp) :: WL_flat_tol = 1.4 !factor by which smallest histogram can be smaller than mean
contains

 !function UpdateParamsLike(Params, fast, dist, i, freeNew, Lik) result(NewLik) 
 ! !For slice sampling movement and likelihood
 ! Type(ParamSet) tmp, Params
 ! integer, intent(in) :: i
 ! real(mcp), intent(in), optional :: Lik
 ! real(mcp) NewLik
 ! logical, intent(in) :: fast, freeNew
 ! real(mcp), intent(in) :: dist
 !
 ! tmp = Params
 !
 ! call UpdateParamsDirection(tmp, fast, dist, i)
 !
 ! NewLik = GetLogLike(tmp)
 ! if (freeNew) then
 !   call AcceptReject(.true., tmp%Info, Params%Info)
 ! else
 !    if (NewLik <= Lik*eps_1) then 
 !      !accept move
 !      call AcceptReject(.false., tmp%Info, Params%Info)
 !      Params = tmp
 !    else
 !      call AcceptReject(.true., tmp%Info, Params%Info)
 !    endif
 ! endif
 !
 !end function UpdateParamsLike
 !
 !
 !subroutine SliceUpdate(Params, fast, i, Like)
 !!Do slice sampling, see http://www.cs.toronto.edu/~radford/slice-aos.abstract.html
 !!update parameter i, input Like updated to new value.
 !!Use linear stepping for the moment
 ! Type(ParamSet) Params
 ! integer, intent(in) :: i
 ! logical, intent(in) :: fast
 ! real(mcp), intent(inout) :: Like
 ! real(mcp) offset,  P
 ! real(mcp) L, R 
 ! real(mcp) LikL, LikR, Range, LikT, w
 ! integer fevals
 !
 !
 !  w = propose_scale
 !  if (slice_randomsize) w =w * Gaussian1()
 !  Like = Like + randexp1()  !New vertical position (likelihood)  
 !  fevals = 0
 !
 !  offset = ranmar() 
 !
 !  L = -offset
 !  R = 1- offset
 !
 !  !step out
 !  if (slice_stepout) then
 !      do
 !       LikL = UpdateParamsLike(Params, fast, w*L, i, .true.)
 !       fevals = fevals + 1
 !       if (LikL*eps_1 > Like) exit
 !       L = L  - 1
 !      end do
 !
 !      do
 !       LikR = UpdateParamsLike(Params, fast, w*R, i,.true.)
 !       fevals = fevals + 1
 !       if (LikR*eps_1 > Like) exit
 !       R = R + 1
 !      end do
 !  end if
 !    
 !  !stepping in
 !  do
 !   Range =  R - L
 !   P  =  ranmar() * Range + L
 !   LikT = UpdateParamsLike(Params, fast, w*P, i, .false., Like)
 !   fevals = fevals + 1
 !   if (LikT < Like*eps_1) exit
 !   if (P > 0) then
 !       R = P
 !   else
 !       L = P
 !   end if
 !  end do
 !  Like = LikT
 !  if (.not. fast) slow_proposals = slow_proposals + fevals
 !
 !end subroutine SliceUpdate
 !
 !
 !subroutine SliceSampleSlowParam(CurParams, CurLike) 
 ! Type(ParamSet) CurParams
 ! real(mcp) CurLike
 ! integer, save:: loopix = 0
 !
 !  if (mod(loopix,num_slow)==0) then
 !       if (.not. allocated(Rot_slow)) allocate(Rot_slow(num_slow,num_slow))
 !       call RotMatrix(Rot_slow,num_slow)
 !       loopix = 0
 !  end if
 !  loopix = loopix + 1
 !  call SliceUpdate(CurParams, .false., loopix,CurLike) 
 !
 !end subroutine SliceSampleSlowParam
 !
 !subroutine SliceSampleFastParams(CurParams, CurLike)
 ! Type(ParamSet) CurParams
 ! real(mcp) CurLike
 ! integer j 
 !
 !  if (.not. allocated(Rot_fast)) allocate(Rot_fast(num_fast,num_fast))
 !  call RotMatrix(Rot_fast,num_fast)
 !  
 !  do j = 1, num_fast * max(1,oversample_fast)
 !     call SliceUpdate(CurParams, .true., mod(j-1,num_fast)+1,CurLike) 
 !  end do
 !  
 !end subroutine SliceSampleFastParams
 !
 function MetropolisAccept(Like, CurLike)
    real(mcp) Like, CurLike
    logical MetropolisAccept
    
    if (Like /=LogZero) then
      MetropolisAccept = CurLike > Like
      if (.not. MetropolisAccept) MetropolisAccept = randexp1() > Like - CurLike
    else
      MetropolisAccept = .false.
    end if
    
 end function MetropolisAccept
 
  subroutine FastDragging(CurParams, CurLike, mult, MaxLike) 
  !Make proposals in fast-marginalized slow parameters 
  !'drag' fast parameters using method of Neal
  Type(ParamSet) TrialEnd, TrialStart, CurParams, CurEndParams, CurStartParams
  real(mcp) CurLike, MaxLike
  integer mult
  integer numaccpt
  real(mcp) CurIntLike, IntLike, CurEndLike, CurStartLike, EndLike, StartLike
  real(mcp) CurDragLike, DragLike
  logical :: accpt
  real(mcp)  :: likes_start_sum, likes_end_sum
  integer interp_step, interp_steps
  real(mcp) frac, delta(num_params)
  integer, save :: num_fast_calls = 0, num_slow_calls = 0

   call Timer()
   TrialEnd = CurParams
   call Proposer%GetProposalSlow(TrialEnd%P)
   
   CurEndLike = GetLogLike(TrialEnd)
   if (CurEndLike==logZero) then
       call AcceptReject(.false., CurParams%Info, TrialEnd%Info)
       mult = mult + 1
       return
   end if
   CurStartLike = CurLike
   if (Feedback > 1) call Timer('Dragging Slow time')

   num_slow_calls = num_slow_calls + 1

   likes_end_sum = CurEndLike
   likes_start_sum = CurStartLike

   CurStartParams = CurParams
   CurEndParams = TrialEnd

   interp_steps = max(2,nint(dragging_steps * num_fast * Proposer%oversample_fast) + 1)

   numaccpt = 0
   do interp_step = 1, interp_steps-1 
    call Proposer%GetProposalFastDelta(delta)
    TrialEnd = CurEndParams
    TrialEnd%P(1:num_params) = TrialEnd%P(1:num_params) + delta 
    EndLike = GetLogLike(TrialEnd)
    accpt = EndLike /= logZero
    if (accpt) then 
     num_fast_calls = num_fast_calls + 1
     TrialStart = CurStartParams
     TrialStart%P(1:num_params) = TrialStart%P(1:num_params)  + delta
     StartLike = GetLogLike(TrialStart)
     accpt = StartLike/=logZero
     if (accpt) then

      num_fast_calls = num_fast_calls + 1
      if (Feedback > 2) print *,'End,start drag: ', interp_step, EndLike, StartLike

      frac = real(interp_step, mcp)/interp_steps
      CurIntLike = CurStartLike*(1-frac) + frac*CurEndLike
      IntLike = StartLike*(1-frac) + frac*EndLike
      accpt = MetropolisAccept(IntLike, CurIntLike)

     end if
    end if
    
    call AcceptReject(accpt, CurEndParams%Info, TrialEnd%Info)
    call AcceptReject(accpt, CurStartParams%Info, TrialStart%Info)

    if (accpt) then
       CurEndParams = TrialEnd
       CurStartParams = TrialStart
       CurEndLike = EndLike
       CurStartLike = StartLike
       CurIntLike = IntLike
       numaccpt = numaccpt + 1
    end if

   likes_start_sum = likes_start_sum + CurStartLike
   likes_end_sum = likes_end_sum + CurEndLike

   end do

   if (Feedback > 1) call Timer('Dragging time')
   
   if (Feedback > 1) print *,'drag steps accept ratio:', &
        real(numaccpt)/(interp_steps)

   CurDragLike = likes_start_sum/interp_steps !old slow
   DragLike = likes_end_sum/interp_steps !proposed new slow
   if (Feedback > 2) print *,'CurDragLike, DragLike: ', CurDragLike, DragLike

   accpt = MetropolisAccept(DragLike, CurDragLike)

   call AcceptReject(accpt, CurParams%Info, CurEndParams%Info)

   if (accpt) then
       if (num_accept> burn_in) then
          output_lines = output_lines +1
          call WriteParams(CurParams,real(mult,mcp),CurLike)
       end if
       num_accept = num_accept + 1
       CurParams = CurEndParams
       CurLike = CurEndLike
       if (CurLike < MaxLike) MaxLike = CurLike
       mult=1
       if (Feedback > 0) write (*,*) mult, ' accepting drag, accpt:', real(num_accept)/num, &
                         'fast/slow',real(num_fast_calls)/num_slow_calls
   else
       mult = mult + 1
   end if

  end subroutine FastDragging


 !subroutine SampleSlowGrid(CurParams, CurLike, mult, MaxLike, num, num_accept) 
 ! !sampling_method = method_slowgrid
 ! !Make proposals on a grid in each slow direction simultaneously with changes to
 ! !fast parameters. Allows fast paraters to 'adapt' to slow parameter values, thereby
 ! !increasing slow acceptance rate
 ! Type(ParamSet) Trial, CurParams
 ! real(mcp) CurLike, MaxLike
 ! integer num, num_accept, mult
 ! real(mcp) w
 ! integer r, r_min, r_max, last_r
 ! Type(ParamSet), dimension(:), allocatable :: grid
 ! integer i, accpt
 ! real(mcp) Like
 ! integer, save:: loopix = 0
 !
 !   TrialEnd = CurParams
 !   call Proposer%GetProposalSlow(TrialEnd%P)
 !
 !  w = max(1e-5,propose_r(num_slow)) * propose_scale / 2
 !
 !  r = 0
 !  r_min = 0
 !  r_max = 0
 !  accpt = 0
 !
 !  allocate(grid(-directional_grid_steps:directional_grid_steps))
 !  grid(0) = CurParams
 !  
 !  do i = 1, directional_grid_steps*oversample_fast
 !
 !   last_r = r
 !
 !   if (ranmar() < 0.5) then
 !    r =r  + 1
 !   else
 !    r =r  - 1
 !   end if
 ! 
 !   Trial = Grid(last_r)
 !   call UpdateParamsDirection(Trial, .false., (r-last_r)*w, loopix)
 !
 !   if (r > r_max  .or. r < r_min) then
 !    grid(r) = Trial
 !    slow_proposals = slow_proposals + 1
 !    r_min = min(r,r_min)
 !    r_max = max(r,r_max)
 !   else
 !    grid(r)%P(fast_params_used) = Trial%P(fast_params_used)
 !   end if
 !
 !   call GetProposalProjFast(grid(r), grid(r), propose_scale)
 ! 
 !   Like = GetLogLike(grid(r))
 !  
 !   if (Feedback > 1) write (*,*) r, 'Likelihood: ', Like, 'Current Like:', CurLike
 !
 !   if (MetropolisAccept(Like, CurLike)) then
 !   !Accept
 !
 !     if (num_accept> burn_in) then
 !         output_lines = output_lines +1
 !        call WriteParams(grid(last_r),real(mult,mcp),CurLike)
 !     end if
 !  
 !     CurLike = Like
 !     accpt = accpt + 1
 !     if (Like < MaxLike) MaxLike = Like
 !     mult=1
 !
 !   else
 !   !Reject
 !     r = last_r
 !     mult=mult+1
 !   end if  
 !
 !   call AddMPIParams(grid(r)%P,CurLike) 
 !
 !  num = num + 1
 !
 !  end do
 !  
 !  if (Feedback > 1) write(*,*) 'Done slow grid, direction', loopix
 !  if (Feedback > 1) write(*,*) 'grid steps moved:',abs(r), ' acc rate =',&
 !           real(accpt)/(directional_grid_steps*oversample_fast)
 !  
 !  !Have to be careful freeing up here.. (logZero results don't allocate new)
 !  do last_r = r_max, r+1, -1
 !     call AcceptReject(.true., grid(last_r)%Info,grid(last_r-1)%Info)  
 !  end do
 !  do last_r = r_min, r-1
 !     call AcceptReject(.true., grid(last_r)%Info,grid(last_r+1)%Info)  
 !  end do
 !
 !  CurParams = grid(r)
 !  num_accept = num_accept + accpt
 !
 !  deallocate(grid)
 !
 !end subroutine SampleSlowGrid

!Multicanonical/Wang-Landau

 subroutine MC_AddLike(L)
  real(mcp), intent(in) :: L
  integer like_ix
  integer, dimension(:,:), allocatable :: tmp
  integer sc

   sc =0

   like_ix = int(L*mc_logspace/num_params_used)
   if (.not. allocated(mc_like_counts)) then
     mc_mini=like_ix-100
     mc_maxi=like_ix+100
     allocate(mc_like_counts(mc_mini:mc_maxi,2))
     mc_like_counts=sc
   end if
   if (like_ix > mc_maxi) then
     allocate(tmp(mc_mini:like_ix+100,2))
     tmp=sc
     tmp(mc_mini:mc_maxi,:) = mc_like_counts
     deallocate(mc_like_counts)
     mc_maxi = like_ix + 100
     allocate(mc_like_counts(mc_mini:mc_maxi,2))
     mc_like_counts = tmp
     deallocate(tmp)
   end if
   if (like_ix < mc_mini) then
     allocate(tmp(mc_mini-100:mc_maxi,2))
     tmp=sc
     tmp(mc_mini:mc_maxi,:) = mc_like_counts
     deallocate(mc_like_counts)
     mc_mini = like_ix - 100
     allocate(mc_like_counts(mc_mini:mc_maxi,2))
     mc_like_counts = tmp
     deallocate(tmp)
   end if

   mc_like_counts(like_ix,1)=mc_like_counts(like_ix,1)+ WL_f
   mc_like_counts(like_ix,2)=mc_like_counts(like_ix,2)+ 1


 end subroutine MC_AddLike

 subroutine MC_UpdateWeights
    real(mcp), dimension(:), allocatable :: tmp
    real(mcp) win(-2:2)
    integer i
         
    allocate(tmp(mc_mini:mc_maxi))
    tmp=0
    tmp(LBOUND(mc_lnweights,DIM=1):UBOUND(mc_lnweights,DIM=1)) = mc_lnweights   
    if (allocated(mc_lnweights)) deallocate(mc_lnweights)
    allocate(mc_lnweights(mc_mini:mc_maxi))
    mc_lnweights = tmp

    tmp = 0
    do i=-2, 2
     win(i) = exp(-real(i**2)*2)
    end do
    do i = mc_mini+2 , mc_maxi-2
       tmp(i) = sum(mc_like_counts(i-2:i+2,1)*win)
    end do 
    !tmp=mc_like_counts
    tmp = max(tmp,1._mcp)
    mc_lnweights = mc_lnweights + log(tmp/real(maxval(tmp))) 
    deallocate(tmp)
    mc_like_counts = 0

 end  subroutine MC_UpdateWeights
 
 function MC_LnWeight(L)
  real(mcp), intent(in) :: L
  real(mcp) MC_LnWeight
  integer like_ix

  if (.not. allocated(mc_lnweights)) then
       MC_LnWeight=0
  else
      like_ix = int(L*mc_logspace/num_params_used)
      if (like_ix > UBound(mc_lnweights,dim=1)) then
       MC_LnWeight = mc_lnweights(UBound(mc_lnweights,dim=1))    
      elseif (like_ix < LBound(mc_lnweights,dim=1)) then
       MC_LnWeight= mc_lnweights(LBound(mc_lnweights,dim=1))    
      else
       MC_LnWeight = mc_lnweights(like_ix)    
      end if
  endif

 end  function MC_LnWeight


 function MC_WeightLike(L) result (Like)
  real(mcp), intent(in) :: L
  real(mcp) Like

   Like = L + MC_LnWeight(L)

 end function MC_WeightLike


 function MC_Weight(L) result(W)
  real(mcp), intent(in) :: L
  real(mcp) mx, W
  integer like_ix

  if (.not. allocated(mc_lnweights)) then
   W=1
  else
   like_ix =int(L*mc_logspace/num_params_used)

   mx = maxval(mc_lnweights)
   W =  exp((MC_LnWeight(L) - mx))
  end if

 end function MC_Weight


 subroutine WL_UpdateWeights
  integer i,n
  real(mcp) amin,asum,L

   if (WL_F > WL_min) then

     amin = LogZero
     asum = 0
     n=0
     do i=mc_mini,mc_maxi
      L = i*num_params_used/mc_logspace 
      if (mc_like_counts(i,1) /= 0 .and. L <= MaxLike + WL_maxL*num_params_used) then
           n=n+1
           amin = min(amin,mc_like_counts(i,2))
           asum = asum + mc_like_counts(i,2)
      end if
     end do
    WL_update_steps = WL_update_steps + 30*num_params_used 

    if (amin > asum / (WL_flat_tol*n)) then
          WL_f = WL_f*0.5
          mc_like_counts(:,2)=0
          WL_update_steps = WL_update_steps + 300*num_params_used  

          if (WL_F < WL_min) then
             if (Feedback > 0) write(*,*) 'Outputting Markovian samples'
             WL_f = 0
             WL_minW=LogZero
!            call CreateTxtFile('density.txt',1)
             do i=mc_mini,mc_maxi
              L = i*num_params_used/mc_logspace 
              if (mc_like_counts(i,1) /= 0 .and. L <= MaxLike + WL_maxL*num_params_used) then
                     WL_minW = min(WL_minW,L - mc_like_counts(i,1))
  !                   write (1,*) L,mc_like_counts(i,1)
              end if
             end do
 !            close(1)
          end if
    end if
   end if
 end  subroutine WL_UpdateWeights

 function WL_WeightLike(L) result (Like)
  real(mcp), intent(in) :: L
  real(mcp) Like
  integer like_ix

  if (L > MaxLike + WL_maxL*num_params_used) then
    Like = LogZero  !!!!
    return
  end if

  like_ix = int(L*mc_logspace/num_params_used)
  if (like_ix >= mc_maxi .or. like_ix <= mc_mini) then
   Like = 0     
  else
   Like = mc_like_counts(like_ix,1)    
  end if

 end function WL_WeightLike


function WL_Weight(L) result (W)
  real(mcp), intent(in) :: L
  real(mcp) W
  integer like_ix

  if (L > MaxLike + WL_maxL*num_params_used) then
   W=0
   return
  end if

  like_ix = int(L*mc_logspace/num_params_used)
  if (like_ix >= mc_maxi .or. like_ix <= mc_mini) then
   W = 0     
  else
   W  = exp( -L + mc_like_counts(like_ix,1) + WL_minW)    
  end if

 end function WL_Weight


 subroutine MCMCsample(Params, samples_to_get)
   integer samples_to_get
   Type(ParamSet) Trial, Params, CurParams
   real(mcp) Like, CurLike
   logical accpt
   real(mcp) rmult
   integer mult, last_num
   integer num_metropolis,  num_metropolis_accept
   real(mcp) testlike, testCurLike
   integer mc_update_burn_num, mc_updates
   character(LEN=128) logLine
   
   MaxLike = LogZero
   CurLike = StartLike
   testCurLike = StartLike
   CurParams = Params
   last_num = 0
   mult= 1
   num_metropolis = 0
   num_metropolis_accept = 0
   mc_update_burn_num = mc_update_burn
   mc_updates = 0

   do while (num <= samples_to_get)

   num = num + 1
 
   if (indep_sample /= 0) then
      if (num_accept>= burn_in .and. mod(num, indep_sample)==0) then
       if (CurLike /= LogZero) call WriteIndepSample(CurParams, CurLike)
      end if
   end if

   if (CurLike /= LogZero .and. sampling_method == sampling_slice) then 
    !Slice sampling
          !if (num_accept> burn_in) then
          ! output_lines = output_lines +1
          ! call WriteParams(CurParams, real(mult,mcp), CurLike)
          !end if
          !if (Feedback > 1) write (*,*) instance, 'Slicing, Current Like:', CurLike
          !mult = 1
          !if (num_slow /=0) call SliceSampleSlowParam(CurParams, CurLike) 
          !if (num_fast /=0) call SliceSampleFastParams(CurParams, CurLike) 
          !num_accept = num_accept + 1
          !if (CurLike < MaxLike) MaxLike = CurLike

  elseif (sampling_method == sampling_slowgrid .and. CurLike /= LogZero .and. num_fast /= 0 .and. num_slow /=0) then
          if (Feedback > 1) write (*,*) instance, 'Directional gridding, Like: ', CurLike
!          call SampleSlowGrid(CurParams, CurLike,mult, MaxLike, num, num_accept) 
  elseif (sampling_method == sampling_fast_dragging .and. CurLike /= LogZero .and. num_fast/=0 .and. num_slow /=0) then
          if (Feedback > 1) write (*,*) instance, 'Fast dragging, Like: ', CurLike
          call FastDragging(CurParams, CurLike,mult, MaxLike) 
  else

   !Do metropolis, except for (optional) slice sampling on fast parameters
   !(slicing dynamically adjusts for the marginal distribution width at the
   ! expense of additional function evaluations)
    if (sampling_method == sampling_fastslice) then

        !if (CurLike /=LogZero .and. num_fast /=0) then
        !    if (num_accept> burn_in) then
        !      output_lines = output_lines +1
        !      call WriteParams(CurParams,real(mult,mcp), CurLike)
        !    end if
        !    call SliceSampleFastParams(CurParams, CurLike)
        !   if (CurLike < MaxLike) MaxLike = CurLike
        !    mult = 1 
        !end if 
        !call GetProposalProjSlow(CurParams, Trial)
    else
        Trial = CurParams
        call Proposer%GetProposal(Trial%P)
    end if

    Like = GetLogLike(Trial) 

    if (Feedback > 1) write (*,*) instance, 'Likelihood: ', Like, 'Current Like:', CurLike

    if (Like /= logZero) then
       if (sampling_method == sampling_multicanonical) then
        testLike = MC_WeightLike(Like)   
       elseif (sampling_method == sampling_wang_landau) then
        testLike = WL_WeightLike(Like)   
       else
        testLike = Like
       end if
       accpt = MetropolisAccept(testLike, testCurLike)
       if (.not. accpt .and. sampling_method == sampling_wang_landau) testCurLike = testCurLike + WL_f
    else
      accpt = .false.
    end if

   call AcceptReject(accpt, CurParams%Info, Trial%Info)
   num_metropolis = num_metropolis +1

   if (accpt) then


      if (num_accept> burn_in .and. CurLike/=LogZero ) then
         if (sampling_method == sampling_multicanonical) then
          if (num>mc_update_burn_num) then
           rmult = mult*MC_Weight(CurLike)
           else
           rmult = 0.
           end if
         else
          rmult = mult
         end if
         if (sampling_method == sampling_wang_landau) then
           if (WL_f >= WL_min) then
             rmult = 0.
           else
             rmult = mult*WL_Weight(CurLike) 
           end if
         end if

         if (rmult /= 0.) then
          output_lines = output_lines +1
          call WriteParams(CurParams,rmult,CurLike)
         end if
      end if
   
      CurLike = Like
      testCurLike = testLike
      CurParams = Trial
      num_accept = num_accept + 1
      num_metropolis_accept = num_metropolis_accept + 1
      if (Like < MaxLike) MaxLike = Like
      mult=1

      if (Feedback > 1) write (*,*) num_metropolis, ' accepting. ratio:', real(num_metropolis_accept)/num_metropolis

      if (mod(num_metropolis_accept,50) ==0) then
          write (*,*) MPIRank, 'rat:',real(num_metropolis_accept)/num_metropolis, ' in ',num_metropolis, &
               ' (M) best: ',real(MaxLike)
      end  if
      if (logfile_unit /=0 .and. mod(num_metropolis_accept,50) ==0) then
           write (logLine,*) 'rat:',real(num_metropolis_accept)/num_metropolis, ' in ',num_metropolis, &
               ', (M) best: ',real(MaxLike)
           call IO_WriteLog(logfile_unit,logLine)
           write (logLine,*) 'local acceptance ratio:', 50./(num_metropolis - last_num) 
           call IO_WriteLog(logfile_unit,logLine)
           last_num = num_metropolis
      end if
    else
       mult=mult+1
       if (Feedback > 1) write (*,*) num_metropolis,' rejecting. ratio:', real(num_metropolis_accept)/num_metropolis
    end if

    end if !not slicing 

    if (CurLike /= logZero) then
      if (sampling_method /= sampling_slowgrid) call AddMPIParams(CurParams%P,CurLike) 
      if (mod(num,100)==0) call CheckParamChange
      
      if (sampling_method == sampling_multicanonical .and. &
              mod(num,mc_update_steps + mc_updates*mc_steps_inc)==0) then
        call MC_UpdateWeights
        Like=LogZero         
        testCurLike = StartLike
        mc_update_burn_num = num + mc_update_burn
        mc_updates = mc_updates + 1
      end if
     
      if (sampling_method == sampling_wang_landau .and. &
              num > WL_update_steps) then
           call WL_UpdateWeights
      end if


      if (sampling_method == sampling_multicanonical .and. num>mc_update_burn_num .or. &
           sampling_method == sampling_wang_landau) then
         call MC_AddLike(CurLike)
      end if

    else
      if (num > 1000) then
        call DoAbort('MCMC.f90: Couldn''t start after 1000 tries - check starting ranges')
      end if    
    end if

   end do

   if (Feedback > 0) then
     write(*,*) MPIRank, 'Stopping as have ',samples_to_get ,' samples. '
     if (use_fast_slow) write(*,*) 'slow changes', slow_changes, 'power changes', power_changes
   end if

 end subroutine MCMCsample

 
end module MonteCarlo
