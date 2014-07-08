! Module for galaxy weak lensing
! AJM, some code from Julien Lesgourgues 

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
     real(mcp), allocatable, dimension(:) :: xi_obs ! Observed correlation functions
     real(mcp), allocatable, dimension(:) :: xi ! Theoretical correlation functions
     integer :: num_z_p ! Source galaxy distribution p(z,bin)
     real(mcp), allocatable, dimension(:,:) :: p
     real(mcp), allocatable, dimension(:) :: z_p
     real(mcp) :: ah_factor ! Anderson-Hartlap factor 
   contains
     procedure :: LogLike => WL_LnLike
     procedure :: ReadIni => WL_ReadIni
     procedure, private :: WL_CFHTLENS_loglike
     procedure, private :: get_convergence
  end type WLLikelihood

  integer, parameter :: nlmax=65
  integer, parameter :: nthetamax = 70
  real(mcp), parameter :: dlnl=0.2d0 
  real(mcp), parameter :: dlntheta = 0.1d0
  real(mcp), parameter :: thetamin = 0.5d0
  real(mcp), parameter :: dx = 0.01d0 
  real(mcp), parameter :: xstop=100.0d0

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
          this%kmax=200.0
          call this%ReadDatasetFile(Ini%ReadFileName(numcat('wl_dataset',i)))
          this%LikelihoodType = 'WL'
          this%needs_powerspectra = .true.
          this%num_z = Ini%Read_Int('nz_wl',100)
          this%max_z = Ini%Read_Double('max_z',10.0d0)
          call LikeList%Add(this)
       end do
       if (Feedback>1) write(*,*) 'read WL data sets'
    end if

  end subroutine WLLikelihood_Add

  subroutine WL_ReadIni(this, Ini)
    class(WLLikelihood) this
    class(TSettingIni) :: Ini
    character(LEN=:), allocatable :: measurements_file, cov_file, window_file
    Type(TTextFile) :: F
    real(mcp) :: dummy1,dummy2,pnorm
    real(mcp), allocatable, dimension(:,:) :: temp
    integer i,iz,it,ib,iopb,j,k,nt

    if (Feedback > 0) write (*,*) 'reading WL data set: '//trim(this%name)

    this%num_z_bins = Ini%Read_Int('num_z_bins')
    this%num_theta_bins = Ini%Read_Int('num_theta_bins')
    this%num_z_p = Ini%Read_Int('num_z_p')
    nt = this%num_z_bins*(1+this%num_z_bins)/2

    allocate(this%theta_bins(this%num_theta_bins))
    allocate(this%z_bins(this%num_z_bins))
    allocate(this%xi_obs(this%num_theta_bins*nt*2))
    allocate(this%xi(this%num_theta_bins*nt*2))
    allocate(this%z_p(this%num_z_p))
    allocate(this%p(this%num_z_p,this%num_z_bins))
    this%ah_factor = Ini%Read_Double('ah_factor',1.0d0)
    
    measurements_file  = Ini%ReadFileName('measurements_file')
    window_file  = Ini%ReadFileName('window_file')
    cov_file  = Ini%ReadFileName('cov_file')
   
    !------------------------------------------------------
    ! TODO - ignore theta range if required - scale diagonal of covariance matrix? 
    !------------------------------------------------------

    allocate(this%wl_cov(2*this%num_theta_bins*nt,2*this%num_theta_bins*nt))
    this%wl_cov = 0
    call File%ReadTextMatrix(cov_file,this%wl_cov)

    if (this%name == 'CFHTLENS_1bin') then

       call F%Open(window_file)
       do iz = 1,this%num_z_p
          read (F%unit,*,iostat=iopb) this%z_p(iz),this%p(iz,1)
       end do
       call F%Close()
          
       call F%Open(measurements_file)
       do it = 1,this%num_theta_bins
          read (F%unit,*,iostat=iopb) this%theta_bins(it),this%xi_obs(it),dummy1,this%xi_obs(it+this%num_theta_bins),dummy2
       end do
       call F%Close()
       
    elseif (this%name == 'CFHTLENS_6bin') then
       
       do ib=1,this%num_z_bins
          call F%Open(window_file(1:index(window_file,'BIN_NUMBER')-1)//IntToStr(ib)//window_file(index(window_file,'BIN_NUMBER')+len('BIN_NUMBER'):len(window_file)))
          do iz=1,this%num_z_p
             read (F%unit,*,iostat=iopb) this%z_p(iz),this%p(iz,ib)
          end do
          call F%Close()
       end do

       call F%Open(measurements_file)
       k = 1
       allocate(temp(2*this%num_theta_bins,nt))
       do i=1,2*this%num_theta_bins
          read (F%unit,*, iostat=iopb) dummy1,temp(i,:)
          if (i.le.this%num_theta_bins) this%theta_bins(i)=dummy1
       end do
       do j=1,nt
          do i=1,2*this%num_theta_bins
             this%xi_obs(k) = temp(i,j)
             k = k + 1
          end do
       end do
       deallocate(temp)
       call F%Close()
       
    else  

       write(*,*)'ERROR: Not yet implemented WL dataset: '//trim(this%name)
       call MPIStop()

    end if

    !Normalize window functions p so \int p(z) dz = 1
    do ib=1,this%num_z_bins
       pnorm = 0
       do iz=2,this%num_z_p
          pnorm = pnorm + 0.5d0*(this%p(iz-1,ib)+this%p(iz,ib))*(this%z_p(iz)-this%z_p(iz-1))
       end do
       this%p(:,ib) = this%p(:,ib)/pnorm
    end do

    ! Apply Anderson-Hartap correction 
    this%wl_cov = this%wl_cov/this%ah_factor

  end subroutine WL_ReadIni

  function WL_LnLike(this, CMB, Theory, DataParams)
    use MatrixUtils
    Class(WLLikelihood) :: this
    Class(CMBParams) CMB
    Class(TCosmoTheoryPredictions), target :: Theory
    real(mcp) :: DataParams(:)
    real(mcp) WL_LnLike

    WL_LnLike=0

    call this%get_convergence(CMB,Theory)

    if (this%name=='CFHTLENS_1bin' .or. this%name=='CFHTLENS_6bin') then
       WL_LnLike = this%WL_CFHTLENS_loglike(CMB,Theory)
    end if
    
  end function WL_LnLike

  function WL_CFHTLENS_loglike(this,CMB,Theory)
    use MatrixUtils
    Class(WLLikelihood) :: this
    Class(CMBParams) CMB
    Class(TCosmoTheoryPredictions), target :: Theory
    real(mcp) WL_CFHTLENS_loglike
    real(mcp), allocatable :: vec(:),cov(:,:)
    integer nt
   
    nt = this%num_z_bins*(1+this%num_z_bins)/2
    allocate(vec(this%num_theta_bins*nt*2))
    allocate(cov(this%num_theta_bins*nt*2,this%num_theta_bins*nt*2))
    vec(:) = this%xi(:)-this%xi_obs(:)
    cov(:,:) = this%wl_cov(:,:)
    WL_CFHTLENS_loglike = Matrix_GaussianLogLike(cov,vec) 
    deallocate(cov,vec)
  
  end function WL_CFHTLENS_loglike

  subroutine get_convergence(this,CMB,Theory)
    use Interpolation
    Class(WLLikelihood) :: this
    Class(CMBParams) CMB
    Class(TCosmoTheoryPredictions), target :: Theory
    type(TCubicSpline),  allocatable :: r_z, dzodr_z, P_z, C_l(:,:)
    type(TCubicSpline),  allocatable :: xi1_theta(:,:),xi2_theta(:,:)
    real(mcp) :: h,z,kh
    real(mcp), allocatable :: r(:),dzodr(:)
    real(mcp), allocatable :: rbin(:),gbin(:,:)
    real(mcp), allocatable :: ll(:),PP(:)
    real(mcp), allocatable :: integrand(:)
    real(mcp), allocatable :: Cl(:,:,:)
    real(mcp), allocatable :: theta(:)
    real(mcp), allocatable :: xi1(:,:,:),xi2(:,:,:) 
    real(mcp) :: khmin, khmax, lmin, lmax, xmin, xmax, x, lll
    real(mcp) :: Bessel0, Bessel4, Cval
    real(mcp) :: i1, i2, i1p, i2p, lp
    integer :: i,ib,jb,il,it,iz,nr,nrs,izl,izh,j
    integer :: num_z

    h = CMB%H0/100 
    num_z = Theory%NL_MPK%ny
    khmin = exp(Theory%NL_MPK%x(1))
    khmax = exp(Theory%NL_MPK%x(Theory%NL_MPK%nx))

    !-----------------------------------------------------------------------
    ! Compute comoving distance r and dz/dr
    !-----------------------------------------------------------------------
    
    allocate(r(num_z),dzodr(num_z))
    do iz=1,num_z
       z = Theory%NL_MPK%y(iz)
       r(iz) = this%Calculator%ComovingRadialDistance(z)
       dzodr(iz) = this%Calculator%Hofz(z)
!!!    DEBUG
!!!       write(*,*) z,r(iz),dzodr(iz)
    end do
    allocate(r_z, dzodr_z)
    call r_z%Init(Theory%NL_MPK%y,r,n=num_z)
    call dzodr_z%Init(Theory%NL_MPK%y,dzodr,n=num_z)
    
    !-----------------------------------------------------------------------
    ! Compute lensing efficiency
    !-----------------------------------------------------------------------

    allocate(rbin(this%num_z_p),gbin(this%num_z_p,this%num_z_bins))
    rbin=0
    gbin=0
    do iz=1,this%num_z_p
       rbin(iz) = r_z%Value(this%z_p(iz))
    end do
    do ib=1,this%num_z_bins
       do nr=2,this%num_z_p-1
          do nrs=nr+1,this%num_z_p
             gbin(nr,ib)=gbin(nr,ib)+0.5*(dzodr_z%Value(this%z_p(nrs))*this%p(nrs,ib)*(rbin(nrs)-rbin(nr))/rbin(nrs) &
                  + dzodr_z%Value(this%z_p(nrs-1))*this%p(nrs-1,ib)*(rbin(nrs-1)-rbin(nr))/rbin(nrs-1))*(rbin(nrs)-rbin(nrs-1))
          end do
       end do
    end do
   
!!!    DEBUG 
!!!    write(*,*) gbin
    
    !-----------------------------------------------------------------------
    ! Find convergence power spectrum using Limber approximation
    !-----------------------------------------------------------------------

    allocate(ll(nlmax),PP(num_z))
    allocate(integrand(this%num_z_p))
    allocate(Cl(nlmax,this%num_z_bins,this%num_z_bins))
    Cl = 0
    do il=1,nlmax

       ll(il)=1.*exp(dlnl*(il-1._mcp))
       PP=0
       do iz=1,num_z
          kh = ll(il)/r(iz)/h ! CAMB wants k/h values 
          z = Theory%NL_MPK%y(iz)
          if ((kh .le. khmin) .or. (kh .ge. khmax)) then
             PP(iz)=0.0d0
          else   
             PP(iz)=Theory%NL_MPK%PowerAt(kh,z)
             !-----------------------------------------------------------------------
             ! TODO - need to apply a correction for theories with anisotropic stress
             !-----------------------------------------------------------------------
          end if
       end do
      
       ! Compute integrand over comoving distance 
       allocate(P_z)
       call P_z%Init(r,PP,n=num_z)
       do ib=1,this%num_z_bins
          do jb=1,this%num_z_bins
             integrand = 0
             do nr=1,this%num_z_p
                integrand(nr) = gbin(nr,ib)*gbin(nr,jb)*(1.0+this%z_p(nr))**2.0*P_z%Value(rbin(nr))
             end do             
             do nr=2,this%num_z_p
                Cl(il,ib,jb)=Cl(il,ib,jb)+0.5d0*(integrand(nr)+integrand(nr-1))*(rbin(nr)-rbin(nr-1)) 
             end do
          end do
       end do
       Cl(il,:,:) = Cl(il,:,:)/h**3.0*9._mcp/4._mcp*(h*1e5_mcp/const_c)**4.0*(CMB%omdm+CMB%omb)**2
       deallocate(P_z)

!!!    DEBUG
!!!       write(*,'(10E15.5)') ll(il),Cl(il,1,1),ll(il)**2*Cl(il,1,1)

    end do

    !-----------------------------------------------------------------------
    ! Convert C_l to xi's
    !-----------------------------------------------------------------------

    !----------------------------------------------------------------------
    ! TODO - option for other observable? 
    !-----------------------------------------------------------------------

    allocate(C_l(this%num_z_bins,this%num_z_bins))
    do ib=1,this%num_z_bins
       do jb=1,this%num_z_bins
          call C_l(ib,jb)%Init(ll,Cl(:,ib,jb),n=nlmax)
       end do
    end do
    lmin=ll(1)
    lmax=ll(nlmax)
    allocate(theta(nthetamax))
    allocate(xi1(nthetamax,this%num_z_bins,this%num_z_bins),xi2(nthetamax,this%num_z_bins,this%num_z_bins))
    xi1 = 0
    xi2 = 0
    do it=1,nthetamax
       theta(it) = thetamin*exp(dlntheta*(it-1._mcp))
       xmin=lmin*theta(it)*pi/(180._mcp*60._mcp) ! Convert from arcmin to radians
       xmax=lmax*theta(it)*pi/(180._mcp*60._mcp) 
       x = xmin
       lp = 0
       i1p = 0
       i2p = 0
       do while(x<xmax .and. x<xstop)
          lll=x/(theta(it)*pi/(180._mcp*60._mcp)) 
          if(lll>lmax) then
             write(*,*)'ERROR: l>lmax: '//trim(this%name)
             call MPIStop()
          end if
          Bessel0 = Bessel_J0(x)
          Bessel4 = Bessel_JN(4,x)
          do ib=1,this%num_z_bins
             do jb=1,this%num_z_bins
                Cval = C_l(ib,jb)%Value(lll)*lll/pi/2._mcp
                i1 = Cval*Bessel0
                i2 = Cval*Bessel4
                xi1(it,ib,jb) = xi1(it,ib,jb)+0.5*(i1p+i1)*(lll-lp)
                xi2(it,ib,jb) = xi2(it,ib,jb)+0.5*(i2p+i2)*(lll-lp)
                i1p = i1
                i2p = i2
             end do
          end do
          x = x+dx
          lp = lll
       end do

!!!    DEBUG
!!!    write(*,*) theta(it),xi1(it,1,1),xi2(it,1,1) 

    end do

    !-----------------------------------------------------------------------
    ! Get xi's in column vector format 
    !-----------------------------------------------------------------------
    
    allocate(xi1_theta(this%num_z_bins,this%num_z_bins),xi2_theta(this%num_z_bins,this%num_z_bins))
    do ib=1,this%num_z_bins
       do jb=1,this%num_z_bins
          call xi1_theta(ib,jb)%Init(theta,xi1(:,ib,jb),n=nthetamax)
          call xi2_theta(ib,jb)%Init(theta,xi2(:,ib,jb),n=nthetamax)
       end do
    end do

    iz = 0
    do izl = 1,this%num_z_bins
       do izh = izl,this%num_z_bins
          iz = iz + 1 ! this counts the bin combinations iz=1 =>(1,1), iz=1 =>(1,2) etc
          do i = 1,this%num_theta_bins
             j = (iz-1)*2*this%num_theta_bins
             this%xi(j+i) = xi1_theta(izl,izh)%Value(this%theta_bins(i))      
             this%xi(this%num_theta_bins + j+i) = xi2_theta(izl,izh)%Value(this%theta_bins(i))    
          end do
       end do
    end do

!!!    DEBUG
!!!    write(*,*) this%theta_bins
!!!    write(*,*) this%xi

    deallocate(r,dzodr)
    deallocate(gbin,rbin)
    deallocate(ll,PP,integrand,Cl)
    deallocate(C_l,theta,xi1,xi2)
    deallocate(xi1_theta,xi2_theta)

  end subroutine get_convergence

end module wl
  
  
