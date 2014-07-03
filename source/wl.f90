! Module for galaxy weak lensing

module wl
  use settings
  use CosmologyTypes
  use CosmoTheory
  use Calculator_Cosmology
  use Likelihood_Cosmology
  implicit none
  private

  type, extends(TCosmoCalcLikelihood) :: WLLikelihood
     real(mcp), allocatable, dimension(:,:) :: wl_invcov,wl_cov
     integer :: num_z_bins
     integer :: num_theta_bins
     real(mcp), allocatable, dimension(:) :: theta_bins
     real(mcp), allocatable, dimension(:) :: z_bins
     real(mcp), allocatable, dimension(:) :: xi
     real(mcp), allocatable, dimension(:,:) :: p
     integer :: num_z_p
     real(mcp), allocatable, dimension(:) :: z_p
     real(mcp) :: ah_factor
   contains
     procedure :: LogLike => WL_LnLike
     procedure :: ReadIni => WL_ReadIni
     procedure, private :: WL_CFHTLENS_loglike
     procedure, private :: get_convergence
  end type WLLikelihood

  logical :: use_wl_lss  = .false.

  public WLLikelihood, WLLikelihood_Add, use_wl_lss
  contains

  subroutine WLLikelihood_Add(LikeList, Ini)
    class(TLikelihoodList) :: LikeList
    class(TSettingIni) :: ini
    Type(WLLikelihood), pointer :: this
    integer numwlsets, i

    if (Ini%Read_Logical('use_WL',.false.)) then
       use_wl_lss = .true.
       numwlsets = Ini%Read_Int('wl_numdatasets',0)
       do i= 1, numwlsets
          allocate(this)
          this%needs_nonlinear_pk = .true.
          this%kmax=100.0
          call this%ReadDatasetFile(Ini%ReadFileName(numcat('wl_dataset',i)))
          this%LikelihoodType = 'WL'
          this%needs_powerspectra = .true.
          this%num_z = Ini%Read_Int('nz_wl',100)
          this%max_z = Ini%Read_Double('max_z',10.0d0)
       end do
       call LikeList%Add(this)
    end if

  end subroutine WLLikelihood_Add

  subroutine WL_ReadIni(this, Ini)
    class(WLLikelihood) this
    class(TSettingIni) :: Ini
    character(LEN=:), allocatable :: measurements_file, cov_file, window_file
    Type(TTextFile) :: F
    real(mcp) :: dummy1,dummy2
    integer i,iopb, nt

    if (Feedback > 0) write (*,*) 'reading WL data set: '//trim(this%name)

    this%num_z_bins = Ini%Read_Int('num_z_bins')
    this%num_theta_bins = Ini%Read_Int('num_theta_bins')
    allocate(this%theta_bins(this%num_theta_bins))
    allocate(this%z_bins(this%num_z_bins))
    nt = this%num_z_bins*(1+this%num_z_bins)/2
    allocate(this%xi(this%num_theta_bins*nt*2))
    this%num_z_p = Ini%Read_Int('num_z_p')
    allocate(this%z_p(this%num_z_p))
    this%ah_factor = Ini%Read_Double('ah_factor',1.0d0)

    measurements_file  = Ini%ReadFileName('measurements_file')
    window_file  = Ini%ReadFileName('window_file')
    cov_file  = Ini%ReadFileName('cov_file')
   
    !------------------------------------------------------
    ! TODO - ignore theta range if required
    !------------------------------------------------------

    if (this%name == 'CFHTLENS_1bin') then

       allocate(this%p(this%num_z_p,1))
       call F%Open(window_file)
       do i = 1,this%num_z_p
          read (F%unit,*,iostat=iopb) this%z_p(i),this%p(i,1)
       end do
          
       call F%Open(measurements_file)
       do i = 1,this%num_theta_bins
          read (F%unit,*,iostat=iopb) this%theta_bins(i),this%xi(i),dummy1,this%xi(i+this%num_theta_bins),dummy2
       end do

       allocate(this%wl_cov(2*this%num_theta_bins*nt,2*this%num_theta_bins*nt))
       this%wl_cov = 0
       call File%ReadTextMatrix(cov_file,this%wl_cov)

    else  
       write(*,*)'ERROR: Not yet implemented WL dataset: '//trim(this%name)
       call MPIStop()
    end if

  end subroutine WL_ReadIni

  function WL_LnLike(this, CMB, Theory, DataParams)
    Class(WLLikelihood) :: this
    Class(CMBParams) CMB
    Class(TCosmoTheoryPredictions), target :: Theory
    real(mcp) :: DataParams(:)
    real(mcp) WL_LnLike

    WL_LnLike=0

    if (this%name=='CFHTLENS_1bin' .or. this%name=='CFHTLENS_6bin') then
       WL_LnLike = this%WL_CFHTLENS_loglike(CMB,Theory)
    end if
    
  end function WL_LnLike

  ! CFHTLENS likelihood
  function WL_CFHTLENS_loglike(this,CMB,Theory)
    Class(WLLikelihood) :: this
    Class(CMBParams) CMB
    Class(TCosmoTheoryPredictions), target :: Theory
    real(mcp) WL_CFHTLENS_loglike
    WL_CFHTLENS_loglike = 0

    call this%get_convergence(CMB,Theory)

    write(*,*) WL_CFHTLENS_loglike

  end function WL_CFHTLENS_loglike

  subroutine get_convergence(this,CMB,Theory)
    Class(WLLikelihood) :: this
    Class(CMBParams) CMB
    Class(TCosmoTheoryPredictions), target :: Theory
    integer :: num_z,i
    real(mcp) :: z
    real(mcp), allocatable :: r(:),dzodr(:)

    ! Compute comoving distance r and dz/dr for redshifts
    num_z = Theory%NL_MPK%ny
    allocate(r(num_z),dzodr(num_z))
    do i=1,num_z
       z = Theory%NL_MPK%y(i)
       r(i) = this%Calculator%ComovingRadialDistance(z)
       dzodr(i) = this%Calculator%Hofz(z)
       write(*,*) r(i),dzodr(i)
    end do
    
    deallocate(r,dzodr)

  end subroutine get_convergence


end module wl
  
  
