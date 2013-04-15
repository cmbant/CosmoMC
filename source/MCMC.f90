    !Do Metropolis-Hastings and slice sampling algorithms
    !Also directional gridding method for fast-slow parameters
    !Provisional implementation of Wang-Landau and Multicanonical methods

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

    subroutine FastDragging(CurParams, CurLike, mult)
    !Make proposals in fast-marginalized slow parameters
    !'drag' fast parameters using method of Neal
    Type(ParamSet) TrialEnd, TrialStart, CurParams, CurEndParams, CurStartParams
    real(mcp) CurLike
    integer mult
    integer numaccpt
    real(mcp) CurIntLike, IntLike, CurEndLike, CurStartLike, EndLike, StartLike
    real(mcp) CurDragLike, DragLike
    logical :: accpt
    real(mcp)  :: likes_start_sum, likes_end_sum
    integer interp_step, interp_steps
    real(mcp) frac, delta(num_params)
    integer, save :: num_fast_calls = 0, num_slow_calls = 0, drag_accpt=0
    integer, save :: counter=0

    counter=counter+1
    if (mod(counter, Proposer%oversample_fast)/=0) then 
        call FastMetropolisHastings(CurParams, CurLike, mult)
        return
    end if

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

    interp_steps = max(2,nint(dragging_steps * num_fast) + 1)

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

    call MoveDone(accpt, CurParams, CurLike, mult)
    if (accpt) then
        CurParams = CurEndParams
        CurLike = CurEndLike
        drag_accpt=drag_accpt+1
        if (Feedback > 0 .and. mod(drag_accpt,30)==0 .or. Feedback>1) &
        write (*,*) trim(concat('Chain:',MpiRank,' drag accpt:')), real(drag_accpt)/(counter/Proposer%oversample_fast), &
        'fast/slow',real(num_fast_calls)/num_slow_calls, 'slow:', num_slow_calls
    end if

    end subroutine FastDragging

    subroutine MetropolisHastings(CurParams, CurLike, mult)
    !Standard metropolis hastings
    Type(ParamSet) CurParams, Trial
    real(mcp) CurLike, Like
    integer mult
    logical :: accpt
    character(LEN=128) logLine
    integer, save :: num_metropolis_accept=0, num_metropolis=0, last_num=0

    Trial = CurParams
    call Proposer%GetProposal(Trial%P)

    Like = GetLogLike(Trial)
    if (Feedback > 1) write (*,*) instance, 'Likelihood: ', Like, 'Current Like:', CurLike

    if (Like /= logZero) then
        accpt = MetropolisAccept(Like, CurLike)
    else
        accpt = .false.
    end if

    call AcceptReject(accpt, CurParams%Info, Trial%Info)
    num_metropolis = num_metropolis + 1

    call MoveDone(accpt, CurParams, CurLike, mult, Proposer%Oversample_fast)

    if (accpt) then
        CurParams = Trial
        CurLike = Like
        num_metropolis_accept = num_metropolis_accept + 1
        if (Feedback > 1) write (*,*) num_metropolis, ' metropolis accept. ratio:', real(num_metropolis_accept)/num_metropolis
        if (logfile_unit /=0 .and. mod(num_metropolis_accept,50*Proposer%Oversample_fast) ==0) then
            write (logLine,*) 'metropolis rat:',real(num_metropolis_accept)/num_metropolis, ' in ',num_metropolis, &
            ', best: ',real(MaxLike)
            call IO_WriteLog(logfile_unit,logLine)
            write (logLine,*) 'local acceptance ratio:', 50./(num_metropolis - last_num)
            call IO_WriteLog(logfile_unit,logLine)
            last_num = num_metropolis
        end if
    end if

    end subroutine MetropolisHastings

    subroutine FastMetropolisHastings(CurParams, CurLike, mult)
    !Metropolis hastings on fast parameters
    Type(ParamSet) CurParams, Trial
    real(mcp) CurLike, Like
    integer mult
    logical :: accpt

    Trial = CurParams
    call Proposer%GetProposalFast(Trial%P)

    Like = GetLogLike(Trial)
    if (Like /= logZero) then
        accpt = MetropolisAccept(Like, CurLike)
    else
        accpt = .false.
    end if

    call AcceptReject(accpt, CurParams%Info, Trial%Info)
    call MoveDone(accpt, CurParams, CurLike, mult)

    if (accpt) then
        CurParams = Trial
        CurLike = Like
    end if

    end subroutine FastMetropolisHastings

    subroutine MoveDone(accpt, CurParams, CurLike, mult, thin_fac)
    logical, intent(in) :: accpt
    Type(ParamSet), intent(in) :: CurParams
    real(mcp) CurLike
    integer, intent(inout):: mult
    integer, intent(in), optional :: thin_fac
    integer, save :: acc = 0, indep_acc=0
    integer thin
    logical want 

    if (accpt) then
        num_accept = num_accept + 1
        thin=1
        if (present(thin_fac)) thin=thin_fac
        want= num_accept> burn_in .and. checkpoint_burn ==0  .and. CurLike /= LogZero
        if (want) then
            acc = acc + mult
            indep_acc= indep_acc + mult
        end if
        if (acc >= thin) then
            if (want) then
                output_lines = output_lines +1
                call WriteParams(CurParams,real(acc/thin,mcp),CurLike)
            end if
            acc = mod(acc, thin)
            if (checkpoint_burn/=0) checkpoint_burn = checkpoint_burn-1
        end if
        if (indep_sample /= 0 .and. indep_acc >= indep_sample*thin .and. want) then
            call WriteIndepSample(CurParams, CurLike,real(indep_acc/(indep_sample*thin),mcp))
            indep_acc = mod(indep_acc, indep_sample*thin)
        end if
        if (CurLike < MaxLike) MaxLike = CurLike
        mult=1
    else
        mult = mult + 1
    end if

    end subroutine MoveDone

    subroutine MCMCsample(Params, samples_to_get)
    integer samples_to_get
    Type(ParamSet) Params, CurParams
    real(mcp)  CurLike
    integer mult

    MaxLike = LogZero
    CurLike = StartLike
    CurParams = Params
    mult= 1

    do while (num <= samples_to_get*Proposer%Oversample_fast)
        num = num + 1

        if (sampling_method == sampling_fast_dragging .and. CurLike /= LogZero .and. num_fast/=0 .and. num_slow /=0) then
            if (Feedback > 1) write (*,*) instance, 'Fast dragging, Like: ', CurLike
            call FastDragging(CurParams, CurLike, mult)
        else
            if(sampling_method/=sampling_metropolis .and. sampling_method /=sampling_fast_dragging) &
            call MpiStop('Sampling method not currently updated to this CosmoMC version')
            call MetropolisHastings(CurParams, CurLike, mult)
        end if

        if (CurLike /= logZero) then
            call AddMPIParams(CurParams%P,CurLike)
            if (mod(num*Proposer%Oversample_fast,100)==0) call CheckParamChange
        else
            if (num > 1000) then
                call DoAbort('MCMC.f90: Couldn''t start after 1000 tries - check starting ranges')
            end if
        end if
    end do

    if (Feedback > 0) then
        write(*,*) MPIRank, 'Stopping as have ',samples_to_get ,' samples. '
        if (use_fast_slow) write(*,*) 'slow changes', slow_changes, 'semi-slow changes', semislow_changes
    end if

    end subroutine MCMCsample


    end module MonteCarlo
