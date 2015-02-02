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
        integer :: num_mask
        real(mcp), allocatable, dimension(:) :: theta_bins
        real(mcp), allocatable, dimension(:) :: z_bins
        real(mcp), allocatable, dimension(:) :: xi_obs ! Observed correlation functions
        integer :: num_z_p ! Source galaxy distribution p(z,bin)
        real(mcp), allocatable, dimension(:,:) :: p
        real(mcp), allocatable, dimension(:) :: z_p
        integer, allocatable, dimension(:) :: mask_indices
        real(mcp) :: ah_factor ! Anderson-Hartlap factor
        logical :: cut_theta  ! Cut non-linear scales
        logical :: use_non_linear ! Whether to use non-linear corrections
        logical :: use_weyl ! Whether to use Weyl potential or matter P(k)
    contains
    procedure :: LogLike => WL_LnLike
    procedure :: ReadIni => WL_ReadIni
    procedure, private :: get_convergence
    end type WLLikelihood

    ! Integration accuracy parameters
    integer, parameter :: nlmax = 65
    real(mcp), parameter :: dlnl = 0.2d0
    real(mcp) :: dlntheta = 0.25d0
    real(mcp), parameter :: dx = 0.02d0
    real(mcp), parameter :: xstop = 200.0d0

    public WLLikelihood, WLLikelihood_Add
    contains

    subroutine WLLikelihood_Add(LikeList, Ini)
    class(TLikelihoodList) :: LikeList
    class(TSettingIni) :: ini
    Type(WLLikelihood), pointer :: this
    Type(TSettingIni) :: DataSets
    integer i
    logical :: nonlinear, useweyl

    if (Ini%Read_Logical('use_WL',.false.)) then
        nonlinear = Ini%Read_Logical('wl_use_non_linear',.true.)
        useweyl = Ini%Read_Logical('wl_use_weyl',.true.)
        call Ini%TagValuesForName('wl_dataset', DataSets)
        do i= 1, DataSets%Count
            allocate(this)
            this%needs_nonlinear_pk = nonlinear
            this%use_non_linear = nonlinear
            this%use_weyl = useweyl
            call this%ReadDatasetFile(DataSets%Value(i))
            this%LikelihoodType = 'WL'
            this%tag = DataSets%Name(i)
            this%needs_powerspectra = .true.
            this%needs_weylpower = useweyl
            call LikeList%Add(this)
        end do
        if (Feedback>1) write(*,*) 'read WL data sets'
    end if

    end subroutine WLLikelihood_Add

    subroutine WL_ReadIni(this, Ini)
    use MatrixUtils
    class(WLLikelihood) this
    class(TSettingIni) :: Ini
    character(LEN=:), allocatable :: measurements_file, window_file, measurements_format
    Type(TTextFile) :: F
    real(mcp) :: dummy1,dummy2,pnorm
    real(mcp), allocatable, dimension(:,:) :: temp
    real(mcp), allocatable, dimension(:,:) :: cut_values
    integer, allocatable, dimension(:) :: mask
    integer i,iz,it,ib,j,k,nt,izl,izh
    real(mcp) :: xi_plus_cut, xi_minus_cut

    if (Feedback > 0) write (*,*) 'reading WL data set: '//trim(this%name)

    this%num_z = Ini%Read_Int('nz_wl',100)
    this%max_z = Ini%Read_Double('max_z')
    this%cut_theta = Ini%Read_Logical('cut_theta',.false.)

    this%num_z_bins = Ini%Read_Int('num_z_bins')
    this%num_theta_bins = Ini%Read_Int('num_theta_bins')
    this%num_z_p = Ini%Read_Int('num_z_p')
    nt = this%num_z_bins*(1+this%num_z_bins)/2

    allocate(this%theta_bins(this%num_theta_bins))
    allocate(this%z_bins(this%num_z_bins))
    allocate(this%xi_obs(this%num_theta_bins*nt*2))
    allocate(this%z_p(this%num_z_p))
    allocate(this%p(this%num_z_p,this%num_z_bins))

    this%kmax = Ini%Read_Double('kmax')

    this%ah_factor = Ini%Read_Double('ah_factor',1.0d0)

    measurements_file  = Ini%ReadFileName('measurements_file')
    measurements_format = Ini%Read_String('measurements_format',NotFoundFail=.true.)

    window_file  = Ini%ReadFileName('window_file')

    call File%LoadTxt(Ini%ReadFileName('cov_file'),this%wl_cov)
    if (size(this%wl_cov,1)/= 2*this%num_theta_bins*nt) call MpiStop('WL_ReadIni: wrong cov size')
    
    if (this%cut_theta) then
        call File%LoadTxt(Ini%ReadFileName('cut_file'), cut_values)
        if (size(cut_values,1)/=this%num_z_bins) call MpiStop('WL_ReadIni: bad cut values')
    end if

    do ib=1,this%num_z_bins
        call F%Open(FormatString(window_file,ib,allow_unused=.true.))
        do iz=1,this%num_z_p
            read (F%unit,*) this%z_p(iz),this%p(iz,ib)
        end do
        call F%Close()
    end do

    call F%Open(measurements_file)
    if (measurements_format == '1bin') then
        do it = 1,this%num_theta_bins
            read (F%unit,*) this%theta_bins(it),this%xi_obs(it),dummy1,this%xi_obs(it+this%num_theta_bins),dummy2
        end do
    elseif (measurements_format == 'nbin') then
        k = 1
        allocate(temp(2*this%num_theta_bins,nt))
        do i=1,2*this%num_theta_bins
            read (F%unit,*) dummy1,temp(i,:)
            if (i.le.this%num_theta_bins) this%theta_bins(i)=dummy1
        end do
        do j=1,nt
            do i=1,2*this%num_theta_bins
                this%xi_obs(k) = temp(i,j)
                k = k + 1
            end do
        end do
        deallocate(temp)
    else
        call MPIStop('ERROR: Not yet implemented WL measurements format: '//measurements_format)
    end if
    call F%Close()

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

    ! Compute theta mask
    allocate(mask(this%num_theta_bins*nt*2))
    if (allocated(cut_values)) then
        mask = 0
        iz = 0
        do izl = 1,this%num_z_bins
            do izh = izl,this%num_z_bins
                iz = iz + 1 ! this counts the bin combinations iz=1 =>(1,1), iz=1 =>(1,2) etc
                do i = 1,this%num_theta_bins
                    j = (iz-1)*2*this%num_theta_bins
                    xi_plus_cut = max(cut_values(izl,1),cut_values(izh,1))
                    xi_minus_cut = max(cut_values(izl,2),cut_values(izh,2))
                    if (this%theta_bins(i)>xi_plus_cut) mask(j+i) = 1
                    if (this%theta_bins(i)>xi_minus_cut) mask(this%num_theta_bins + j+i) = 1
                    ! Testing
                    !write(*,'(5i4,3E15.3,2i4)') izl,izh,i,i+j,this%num_theta_bins + j+i,xi_plus_cut,&
                    !     xi_minus_cut,this%theta_bins(i),mask(j+i),mask(this%num_theta_bins + j+i)
                end do
            end do
        end do
    else
        mask = 1
    end if
    this%num_mask = sum(mask)
    allocate(this%mask_indices(this%num_mask))
    j = 1
    do i=1,this%num_theta_bins*nt*2
        if (mask(i) == 1) then
            this%mask_indices(j) = i
            j = j+1
        end if
    end do

    ! Precompute masked inverse
    allocate(this%wl_invcov(this%num_mask,this%num_mask))
    this%wl_invcov = this%wl_cov(this%mask_indices,this%mask_indices)
    call Matrix_Inverse(this%wl_invcov)

    end subroutine WL_ReadIni

    function WL_LnLike(this, CMB, Theory, DataParams)
    use MatrixUtils
    Class(WLLikelihood) :: this
    Class(CMBParams) CMB
    Class(TCosmoTheoryPredictions), target :: Theory
    real(mcp) :: DataParams(:)
    real(mcp) WL_LnLike
    real(mcp) vec(this%num_mask)
    real(mcp), allocatable :: xi(:)

    allocate(xi(size(this%xi_obs)))
    call this%get_convergence(CMB,Theory, xi)
    vec = xi(this%mask_indices)-this%xi_obs(this%mask_indices)
    WL_LnLike = Matrix_QuadForm(this%wl_invcov,vec) / 2

    end function WL_LnLike


    subroutine get_convergence(this,CMB,Theory,xi)
    use Interpolation
    Class(WLLikelihood) :: this
    Class(CMBParams) CMB
    Class(TCosmoTheoryPredictions), target :: Theory
    real(mcp), intent(out) :: xi(:)
    type(TCosmoTheoryPK), pointer :: PK
    type(TCubicSpline) :: r_z, dzodr_z, P_z
    type(TCubicSpline),  allocatable :: C_l(:,:), xi1_theta(:,:),xi2_theta(:,:)
    real(mcp) :: h,z,kh,k
    real(mcp), allocatable :: r(:),dzodr(:)
    real(mcp), allocatable :: rbin(:),gbin(:,:)
    real(mcp), allocatable :: ll(:),PP(:)
    real(mcp), allocatable :: integrand(:)
    real(mcp), allocatable :: Cl(:,:,:)
    real(mcp), allocatable :: theta(:)
    real(mcp), allocatable :: xi1(:,:,:),xi2(:,:,:)
    real(mcp), allocatable :: i1p(:,:),i2p(:,:)
    real(mcp) :: khmin, khmax, lmin, lmax, xmin, xmax, x, lll
    real(mcp) :: Bessel0, Bessel4, Cval
    real(mcp) :: i1, i2, lp
    real(mcp) :: a2r
    real(mcp) :: thetamin, thetamax
    integer :: i,ib,jb,il,it,iz,nr,nrs,izl,izh,j
    integer :: num_z, ntheta, min_iz

    if (this%use_non_linear) then
        if (this%use_weyl) then
            PK => Theory%NL_MPK_WEYL
        else
            PK => Theory%NL_MPK
        end if
    else
        if (this%use_weyl) then
            PK => Theory%MPK_WEYL
        else
            PK => Theory%MPK
        end if
    end if

    h = CMB%H0/100
    num_z = PK%ny
    khmin = exp(PK%x(1))
    khmax = exp(PK%x(PK%nx))
    a2r = pi/(180._mcp*60._mcp)

    !-----------------------------------------------------------------------
    ! Compute comoving distance r and dz/dr
    !-----------------------------------------------------------------------

    allocate(r(num_z),dzodr(num_z))
    min_iz=1
    do iz=1,num_z
        z = PK%y(iz)
        if (z==0) min_iz=iz+1
        r(iz) = this%Calculator%ComovingRadialDistance(z)
        dzodr(iz) = this%Calculator%Hofz(z)
    end do
    call r_z%Init(PK%y,r,n=num_z)
    call dzodr_z%Init(PK%y,dzodr,n=num_z)

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
        do iz=min_iz,num_z
            k = ll(il)/r(iz)
            kh = k/h ! CAMB wants k/h values
            z = PK%y(iz)
            if ((kh .le. khmin) .or. (kh .ge. khmax)) then
                PP(iz)=0.0d0
            else
                PP(iz)= PK%PowerAt(kh,z)
                ! Testing
                !write(*,'(10E15.5)') k,z,PK_WEYL%PowerAt(kh,z)*k,9.0/(8.0*pi**2.0)&
                !*PK%PowerAt(kh,z)/(h**3.0)*(h*1e5_mcp/const_c)**4.0*(CMB%omdm+CMB%omb)**2*(1+z)**2.0
            end if
        end do

        ! Compute integrand over comoving distance
        call P_z%Init(r,PP,n=num_z)
        do ib=1,this%num_z_bins
            do jb=1,this%num_z_bins
                integrand = 0
                do nr=1,this%num_z_p
                    integrand(nr) = gbin(nr,ib)*gbin(nr,jb)*P_z%Value(rbin(nr))
                    if (.not. this%use_weyl) integrand(nr) = integrand(nr) *(1.0+this%z_p(nr))**2.0
                end do
                do nr=2,this%num_z_p
                    Cl(il,ib,jb)=Cl(il,ib,jb)+0.5d0*(integrand(nr)+integrand(nr-1))*(rbin(nr)-rbin(nr-1))
                end do
            end do
        end do
        if (.not. this%use_weyl) then
            Cl(il,:,:) = Cl(il,:,:)/h**3.0*9._mcp/4._mcp*(h*1e5_mcp/const_c)**4.0*(CMB%omdm+CMB%omb)**2
        end if
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

    thetamin = minval(this%theta_bins)*0.8
    thetamax = maxval(this%theta_bins)*1.2
    ntheta = ceiling(log(thetamax/thetamin)/dlntheta) + 1

    lmin=ll(1)
    lmax=ll(nlmax)
    allocate(theta(ntheta))
    allocate(xi1(ntheta,this%num_z_bins,this%num_z_bins),xi2(ntheta,this%num_z_bins,this%num_z_bins))
    allocate(i1p(this%num_z_bins,this%num_z_bins),i2p(this%num_z_bins,this%num_z_bins))
    xi1 = 0
    xi2 = 0

    do it=1,ntheta
        theta(it) = thetamin*exp(dlntheta*(it-1._mcp))
        xmin=lmin*theta(it)*a2r! Convert from arcmin to radians
        xmax=lmax*theta(it)*a2r
        x = xmin
        lp = 0
        i1p = 0
        i2p = 0
        do while(x<xstop .and. x<xmax)
            lll=x/(theta(it)*a2r)
            if(lll>lmax) then
                write(*,*)'ERROR: l>lmax: '//trim(this%name)
                call MPIStop()
            end if
            Bessel0 = Bessel_J0(x)
            Bessel4 = Bessel_JN(4,x)
            do ib=1,this%num_z_bins
                do jb=ib,this%num_z_bins
                    Cval = C_l(ib,jb)%Value(lll)*lll
                    i1 = Cval*Bessel0
                    i2 = Cval*Bessel4
                    xi1(it,ib,jb) = xi1(it,ib,jb)+0.5*(i1p(ib,jb)+i1)*(lll-lp)
                    xi2(it,ib,jb) = xi2(it,ib,jb)+0.5*(i2p(ib,jb)+i2)*(lll-lp)
                    i1p(ib,jb) = i1
                    i2p(ib,jb) = i2
                end do
            end do
            x = x+dx
            lp = lll
        end do

        do ib=1,this%num_z_bins
            do jb=ib,this%num_z_bins
                xi1(it,jb,ib) = xi1(it,ib,jb)
                xi2(it,jb,ib) = xi2(it,ib,jb)
            end do
        end do

    end do

    xi1=xi1/pi/2._mcp
    xi2=xi2/pi/2._mcp

    !-----------------------------------------------------------------------
    ! Get xi's in column vector format
    !-----------------------------------------------------------------------

    allocate(xi1_theta(this%num_z_bins,this%num_z_bins),xi2_theta(this%num_z_bins,this%num_z_bins))
    do ib=1,this%num_z_bins
        do jb=1,this%num_z_bins
            call xi1_theta(ib,jb)%Init(theta,xi1(:,ib,jb),n=ntheta)
            call xi2_theta(ib,jb)%Init(theta,xi2(:,ib,jb),n=ntheta)
        end do
    end do

    iz = 0
    do izl = 1,this%num_z_bins
        do izh = izl,this%num_z_bins
            iz = iz + 1 ! this counts the bin combinations iz=1 =>(1,1), iz=1 =>(1,2) etc
            do i = 1,this%num_theta_bins
                j = (iz-1)*2*this%num_theta_bins
                xi(j+i) = xi1_theta(izl,izh)%Value(this%theta_bins(i))
                xi(this%num_theta_bins + j+i) = xi2_theta(izl,izh)%Value(this%theta_bins(i))
            end do
        end do
    end do

    end subroutine get_convergence

    end module wl


