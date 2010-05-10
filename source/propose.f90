!Proposal density 
!Using the covariance matrix to give the proposal directions typically
!significantly increases the acceptance rate and gives faster movement 
!around parameter space.
!We generate a random basis in the eigenvectors, then cycle through them
!proposing changes to each, then generate a new random basis
!The distance proposal in the random direction is given by a two-D Gaussian 
!radial function mixed with an exponential, which is quite robust to wrong width estimates
!See http://cosmologist.info/notes/cosmomc.ps.gz

module propose
  use settings
  use Random
  use ParamDef
  implicit none

 logical :: propose_rand_directions = .true.
 logical :: fast_slicing = .false. !slice fast parameters when using Metropolis 
 logical :: slice_stepout = .true.
 logical :: slice_randomsize = .false.

 real, allocatable, dimension(:,:), save :: Rot_slow, Rot_fast

contains

 subroutine RotMatrix(M,n)
  integer, intent(in) :: n
  real M(n,n)
  integer i
  
   if (propose_rand_directions .and. n > 1) then
      call RandRotation(M, n)      
   else
      M = 0
      do i = 1, n
           M(i,i) = sign(1., real(ranmar())-0.5)
      end do
   end if
 
 end  subroutine RotMatrix


 function Propose_r(in_n) result(r_fac)
  !distance proposal function (scale is 1)
   integer, intent(in) :: in_n
   integer i,n
   real r_fac

      if (ranmar() < 0.33d0) then
       r_fac = randexp1()
      else
       n = min(in_n,2)
       r_fac = 0
       do i = 1, n
        r_fac = r_fac + Gaussian1()**2
       end do
       r_fac = sqrt( r_fac / n )
      end if

 end function Propose_r


 subroutine UpdateParamsDirection(tmp, fast, dist,i)
  !Change parameters in tmp by dist in the direction of the ith random e-vector
   Type(ParamSet) tmp
   real vec(num_params_used)
   integer, intent(in) :: i
   logical, intent(in) :: fast
   real, intent(in) :: dist

  if (has_propose_matrix) then
    vec = 0
    if (fast) then
      vec(1:num_fast) =  Rot_fast(:,i) * dist * propose_diag_fast  
      tmp%P(fast_params_used) =  tmp%P(fast_params_used) + & 
        sigmas(num_slow+1:num_slow+ num_fast) * matmul (propose_matrix_fast, vec(1:num_fast))
    else
      vec(1:num_slow) =  Rot_slow(:,i) * dist * propose_diag(slow_evecs) 
      tmp%P(params_used) =  tmp%P(params_used) + &
            sigmas * matmul (propose_matrix(:,slow_evecs), vec(1:num_slow))
    end if
  else
    if (fast) then
     tmp%P(fast_params_used) = tmp%P(fast_params_used) + &
              Rot_fast(:,i) * dist *  Scales%PWidth(fast_params_used)
    else
     tmp%P(params_used(1:num_slow)) = tmp%P(params_used(1:num_slow)) +  &
        Rot_slow(:,i) * dist *  Scales%PWidth(params_used(1:num_slow))
    end if
  end if
 
 end subroutine UpdateParamsDirection


subroutine GetProposal(In, Out)

  type (ParamSet) In, Out
  integer, save :: fast_ix = 0

  fast_ix = fast_ix + 1
  if (num_slow ==0 .or. num_fast /= 0 .and. &
       mod(fast_ix, 2*(1 + (num_fast*oversample_fast)/num_slow)) /= 0) then
   call GetProposalProjFast(In, Out)
  else
   call GetProposalProjSlow(In, Out)
  end if

end subroutine GetProposal

subroutine GetProposalProjSlow(In, Out)
  use settings
  use Random
  use ParamDef
  implicit none
  type (ParamSet) In, Out
  real  wid
  integer, save :: loopix = 0
 
  slow_proposals = slow_proposals + 1
!Change only slow params, or eigenvectors of covmat which most change the slow params

   Out= In
   wid = propose_scale
   if (mod(loopix,num_slow)==0) then
        if (.not. allocated(Rot_slow)) allocate(Rot_slow(num_slow,num_slow))
        call RotMatrix(Rot_slow, num_slow)      
        loopix = 0
   end if
   loopix = loopix + 1

   call UpdateParamsDirection(Out,.false.,Propose_r(num_slow) * wid, loopix)

end subroutine GetProposalProjSlow


subroutine GetProposalProjFast(In, Out, ascale)
  type (ParamSet) In, Out
  real wid
  real, intent(in), optional :: ascale
  integer, save :: loopix = 0

   Out= In

   wid = propose_scale
   if (present(ascale)) wid = ascale

   if (mod(loopix,num_fast)==0) then
        if (.not. allocated(Rot_fast)) allocate(Rot_fast(num_fast,num_fast))
        call RotMatrix(Rot_fast, num_fast)      
        loopix = 0
   end if
   loopix = loopix + 1

   call UpdateParamsDirection(Out,.true.,Propose_r(num_fast) * wid, loopix)
     
end subroutine GetProposalProjFast

end module propose
