!Proposal density 
!Using the covariance matrix to give the proposal directions typically
!significantly increases the acceptance rate and gives faster movement 
!around parameter space.
!We generate a random basis in the eigenvectors, then cycle through them
!proposing changes to each, then generate a new random basis
!The distance proposal in the random direction is given by a two-D Gaussian 
!radial function mixed with an exponential, which is quite robust to wrong width estimates
!See http://cosmologist.info/notes/cosmomc.pdf

module propose
  use Random
  use settings
  implicit none
 
  Type IndexCycler
      integer n
      integer :: loopix=0
  end Type IndexCycler
  
  Type, extends(IndexCycler) :: CyclicIndexRandomizer
      integer, allocatable :: indices(:)
  contains
     procedure :: Next
  end Type CyclicIndexRandomizer
  
  Type, extends(IndexCycler) :: RandDirectionProposer
    integer :: propose_count = 0
    real(mcp), allocatable, dimension(:,:) :: R
  contains
   procedure :: ProposeVec
   procedure :: Propose_r
  end type RandDirectionProposer
  
  Type, extends(RandDirectionProposer) :: BlockProposer
      integer :: block_start
      integer, allocatable :: used_param_indices(:)
      integer, allocatable :: params_changed(:), used_params_changed(:)
      real(mcp), allocatable :: mapping_matrix(:,:)
  contains
    procedure :: UpdateParams
  end Type BlockProposer
  
  Type BlockedProposer
    integer :: nBlocks
    integer, allocatable :: indices(:), ProposerForIndex(:)
    Type(CyclicIndexRandomizer) :: Slow, Fast, All
    Type(BlockProposer), allocatable :: Proposer(:)
  contains
   procedure :: Init
   procedure :: SetCovariance
   procedure :: GetBlockProposal
   procedure :: GetProposal
   procedure :: GetProposalSlow
   procedure :: GetProposalFastDelta
  end Type BlockedProposer

  type int_arr_pointer
    integer, dimension(:), pointer :: p 
  end type int_arr_pointer

 logical :: propose_rand_directions = .true.

contains

 function Next(this)
  Class(CyclicIndexRandomizer) :: this
  integer Next
   this%loopix = mod(this%loopix, this%n) + 1
   if (this%loopix == 1) then
       if (.not. allocated(this%indices)) allocate(this%indices(this%n))
       call randIndices(this%indices, this%n, this%n)
   end if
   Next = this%indices(this%loopix)

 end function Next
  
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

    
  function ProposeVec(this, ascale) result(vec)
   Class(RandDirectionProposer) :: this
   real(mcp) :: vec(this%n)
   real(mcp), intent(in), optional :: ascale
   real(mcp) ::  wid 
   
   if (present(ascale)) then
    wid = ascale
   else
    wid = propose_scale
   end if 
   
   if (mod(this%loopix,this%n)==0) then
       !Get a new random rotation
        if (.not. allocated(this%R)) allocate(this%R(this%n,this%n))
        call RotMatrix(this%R, this%n)
        this%loopix = 0
   end if
   this%loopix = this%loopix + 1
   this%propose_count= this%propose_count + 1
   vec = this%R(:,this%loopix) * ( this%Propose_r() * wid )
   
  end function ProposeVec

 function Propose_r(this) result(r_fac)
  !distance proposal function (scale is 1)
   Class(RandDirectionProposer), intent(in) :: this
   integer i,n
   real(mcp) r_fac

      if (ranmar() < 0.33d0) then
       r_fac = randexp1()
      else
       n = min(this%n, 2)
       r_fac = 0
       do i = 1, n
        r_fac = r_fac + Gaussian1()**2
       end do
       r_fac = sqrt( r_fac / n )
      end if

 end function Propose_r

 
 subroutine UpdateParams(this, P, vec)
  Class(BlockProposer) :: this
  real(mcp) :: P(:)
  real(mcp):: vec(this%n)

   P(this%params_changed) =  P(this%params_changed) + matmul(this%mapping_matrix, vec)
  
 end subroutine UpdateParams

  subroutine Init(this, parameter_blocks, slow_block_max)
  !slow_block_max determines which parameter blocks are grouped together as being 'slow'
   Class(BlockedProposer), target :: this
   type(int_arr_pointer) :: parameter_blocks(:)
   integer, intent(in) :: slow_block_max
   integer used_blocks(size(parameter_blocks))
   integer i, ix, n, np
   Type(BlockProposer), pointer :: BP
   
   this%nBlocks = size(parameter_blocks) 
   n=0
   this%All%n=0
   this%Slow%n=0
   do i=1, this%nBlocks
       np = size(parameter_blocks(i)%P)
       if (np >0) then
        this%All%n = this%All%n+ np
        if (i <= slow_block_max) this%Slow%n= this%Slow%n + np
        n=n+1
        used_blocks(n)=i
       end if
   end do
   this%Fast%n = this%All%n - this%Slow%n
   this%nBlocks = n
   allocate(this%Proposer(this%nBlocks))
   allocate(this%indices(this%All%n))
   allocate(this%ProposerForIndex(this%All%n))
   this%indices=0
   ix = 1
   do i=1, this%nBlocks
       !loop from slow to fast
       BP => this%Proposer(i)
       BP%block_start = ix
       BP%n = size(parameter_blocks(used_blocks(i))%P)
       allocate(BP%used_param_indices(BP%n))
       BP%used_param_indices = parameter_blocks(used_blocks(i))%P
       this%indices(ix:ix+BP%n-1) = BP%used_param_indices
       this%ProposerForIndex(ix:ix+BP%n-1) = i
       ix = ix + BP%n
   end do
   if (any(this%indices==0)) stop 'DecomposeCovariance: not all used parameters blocked'
   do i=1, this%nBlocks
       BP => this%Proposer(i)
       allocate(BP%used_params_changed(this%All%n-BP%block_start+1))
       allocate(BP%params_changed(this%All%n-BP%block_start+1))
       BP%used_params_changed = this%indices(BP%block_start:this%All%n)
       BP%params_changed = params_used(BP%used_params_changed)
   end do

  end subroutine Init

    subroutine SetCovariance(this, propose_matrix)
  !take covariance of used parameters (propose_matrix), and construct orthonormal parmeters 
  !where orthonormal parameters are grouped in blocks by speed, so changes in slowest block 
  !changes slow and fast parameters, but changes in the fastest block only changes fast parameters
   use MatrixUtils
   Class(BlockedProposer), target :: this
   real(mcp), intent(in) :: propose_matrix(num_params_used,num_params_used)
   real(mcp) corr(num_params_used,num_params_used)
   real(mcp) :: sigmas(num_params_used), L(num_params_used,num_params_used)
   integer i, j
   Type(BlockProposer), pointer :: BP

   do i = 1, num_params_used
     sigmas(i) = sqrt(propose_matrix(i,i))
     corr(i,:) = propose_matrix(i,:) / sigmas(i)
   end do
   do i = 1, num_params_used
     corr(:,i) = corr(:,i) / sigmas(i)
   end do
   L = corr(this%indices,this%indices)
   call Matrix_Cholesky(L,zeroed=.true.)
   do i=1, this%nBlocks
       BP => this%Proposer(i)
       if (.not. allocated(BP%mapping_matrix)) allocate(BP%mapping_matrix(size(BP%used_params_changed), BP%n))
       do j = 1, size(BP%used_params_changed)
           BP%mapping_matrix(j,:)  =  sigmas(BP%used_params_changed(j)) * &
                        L(BP%block_start+j-1,BP%block_start:BP%block_start+BP%n-1) 
       end do
   end do
   !For two blocks, fast and slow, the effect is like this
   !param_transform_slow = L(1:num_params_used,1:num_slow)
   !do i=1, num_slow
   ! param_transform_slow(slow_in_used(i),:) =  sigmas(slow_in_used(i)) * L(i,1:num_slow) 
   !end do
   !do i=1, num_fast
   ! param_transform_slow(fast_in_used(i),:) =  sigmas(fast_in_used(i)) * L(num_slow+i,1:num_slow) 
   !end do 
    end subroutine SetCovariance

    
  subroutine GetBlockProposal(this, P, i)
  class(BlockedProposer) :: this
  real(mcp) :: P(:)
  integer, intent(in) :: i

   call this%Proposer(i)%UpdateParams(P, this%Proposer(i)%Proposevec()) 

  end subroutine GetBlockProposal


  subroutine GetProposal(this, P)
   class(BlockedProposer) :: this
   real(mcp) :: P(:) 

   call this%GetBlockProposal(P, this%ProposerForIndex(this%All%Next()))

  end subroutine GetProposal

  subroutine GetProposalSlow(this, P)
   class(BlockedProposer) :: this
   real(mcp) :: P(:) 

   call this%GetBlockProposal(P, this%ProposerForIndex(this%Slow%Next()))

  end subroutine GetProposalSlow
  
  subroutine GetProposalFastDelta(this, P)
   class(BlockedProposer) :: this
   real(mcp) :: P(:) 

   P=0
   call this%GetBlockProposal(P, this%ProposerForIndex(this%Slow%n + this%Fast%Next()))

  end subroutine GetProposalFastDelta

end module propose
