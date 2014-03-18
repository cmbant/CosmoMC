    !Pseudo-Cl (or other C_l esimator) based likelihood approximation for cut sky with polarization
    !AL Mar 2010 - fixed bug using on E, added support for multiple input simulated Chat for bias testing
    !Apr 2011, added fullsky_exact_fksy
    !Mar 2014  renaned from Planck_like, removed low-L likelihood to separate file
    !Updated file structure, allows for binned HL likelihoods. Everything now in L(L+1)CL/2pi

    module CMBLikes
    use settings
    use CosmologyTypes
    use MatrixUtils
    use CosmoTheory
    use Likelihood_Cosmology
    implicit none
    private

    logical, parameter :: bin_test = .false.

    Type TSqMatrix
        real(mcp), dimension(:,:), allocatable :: M
    end Type TSqMatrix

    Type, extends(TCMBLikelihood) :: TCMBLikes
        logical :: highl_cl = .true.
        integer nfields !number of fields
        logical use_field(3)
        integer field_index(3) !mapping from 1,2,3 to index actually used
        integer fields(3) !indices (1=T,2=E,3=B) of fields to use
        character field_order(3)
        integer ncl !calculated from above = nfields*(nfields+1)/2
        integer ncl_used !Number of C_l actually used in covariance matrix (others assumed zero)
        integer, allocatable :: cl_use_index(:)
        integer pcl_lmin, pcl_lmax !The range of l to use psuedo-cl-based likelihood
        integer like_approx
        real(mcp) :: fullsky_exact_fksy = 1 ! only used for testing with exactly fullsky
        integer ncl_hat !1, or more if using multiple simulations
        real(mcp), dimension(:,:), allocatable :: ClFiducial, ClNoise
        real(mcp), dimension(:,:,:), allocatable :: ClHat
        !These are all L(L+1)C_L/2pi

        integer cl_phi_lmin, cl_phi_lmax !lmax for the lensing reconstruction
        integer lensing_recon_ncl !0 for no lensing recon, 1 for phi-phi spectru, 2 phi-phi and phi-T
        integer phi_like_approx
        real(mcp), dimension(:,:), allocatable :: ClPhiHat, ClPhiNoise, phi_inv_covariance
        !for lensing reconstruction
        !note these are [l(l+1)]^4C_l/2pi

        integer bin_width
        integer vecsize
        logical binned
        integer nbins
        integer, allocatable :: bin_cols_in(:,:), bin_cols_out(:)
        integer bin_min, bin_max
        integer, allocatable :: cov_cl_used(:)
        real(mcp), dimension(:,:), allocatable :: inv_covariance
        real(mcp), dimension(:,:,:), allocatable :: binWindows

        Type(TSqMatrix), dimension(:), allocatable :: sqrt_fiducial, NoiseM
        Type(TSqMatrix), dimension(:,:), allocatable :: ChatM
    contains
    procedure :: LogLike => CMBLikes_LogLike
    procedure :: ReadIni => CMBLikes_ReadIni
    procedure, private :: ReadClArr => CMBLikes_ReadClArr
    procedure, private :: ReadLensingReconData=> CMBLikes_ReadLensingReconData
    procedure, private :: ReadClPhiArr => CMBLikes_ReadClPhiArr
    procedure, private :: UseString_to_cols
    procedure, private :: UseString_to_TEB
    procedure, private :: Transform => CMBLikes_Transform
    procedure, private :: LensRecon_Like => CMBLikes_LensRecon_Like
    procedure, private :: ExactChisq
    procedure, private :: MatrixToElements
    procedure, private :: MatrixToElementsInt
    procedure, private :: ElementsToMatrix
    procedure, private :: SetTopHatWindows
    procedure, private :: GetColsFromOrder
    procedure :: ReadBinWindows
    procedure :: ReadCovmat
    procedure :: GetBinnedTheory
    end Type TCMBLikes

    character(LEN=3), parameter :: field_names = 'TEB'
    integer, parameter :: like_approx_HL=1   !new approximation from Hammimeche & Lewis arXiv: 0801.0554
    integer, parameter :: like_approx_fid_gaussian=2 !fiducial fixed covariance matrix, (X-Xhat)^TC^{-1}(X-Xhat)
    integer, parameter :: like_approx_fullsky_exact=3 !ignore all correlations, use exact full sky likelihood function

    character(LEN=Ini_Enumeration_Len), parameter :: &
    & like_Names(3) = [character(Ini_Enumeration_Len)::'HL','gaussian','exact']

    public TCMBLikes
    contains

    function TypeIndex(C)
    character, intent(in) :: C
    integer TypeIndex
    !Get order T, E, B -> 1,2,3
    if (C=='T') then
        TypeIndex=1
    elseif (C=='E') then
        TypeIndex=2
    elseif (C=='B') then
        TypeIndex=3
    else
        call mpiStop('Invalid C_l part, must be T E or B')
    end if

    end function TypeIndex

    subroutine CMBLikes_ReadClArr(this, filename, order, Cl, lmin)
    class(TCMBLikes) :: this
    character(LEN=*), intent(in) :: filename, order
    integer, intent(in) :: lmin
    real(mcp) :: Cl(:,lmin:)
    character(LEN=:), allocatable :: tmp
    integer ix, i, j,i1,l,ll
    integer, allocatable :: cols(:)
    real(mcp), allocatable :: tmp_ar(:)
    integer status, norder
    Type(TTextFile) :: F

    norder = this%GetColsFromOrder(Order, cols)
    allocate(tmp_ar(norder))

    Cl=0
    do while (F%ReadNextContentLine(filename,tmp))
        read(tmp,*, iostat=status) l, tmp_ar
        if (status/=0) call MpiStop('CMBLikes_ReadClArr: error reading line '//trim(filename))
        ll=l
        if (l >=this%bin_min .and. l<=this%bin_max) then
            do ix=1,this%ncl
                if (cols(ix)/=0) Cl(ix,l) = tmp_ar(cols(ix))
            end do
        end if
    end do
    if (ll < this%bin_max) then
        write(*,*) 'CMBLikes_ReadClArr: C_l file does not go up to lmax:', this%pcl_lmax
        write (*,*) trim(filename)
        call MpiStop()
    end if

    end subroutine CMBLikes_ReadClArr

    subroutine UseString_to_cols(this, S, cols)
    class(TCMBLikes) :: this
    character(LEN=*), intent(in) :: S
    integer, allocatable, intent(out) :: cols(:)
    Type(TStringList) :: L
    integer tmp,i,i1,i2,ii,jj,ix
    character(LEN=:), pointer :: P

    call L%SetFromString(S,field_names)
    allocate(cols(L%Count), source=0)

    do i=1, L%Count
        P=> L%Item(i)
        if (len(P)/=2) call mpiStop('Invalid C_l order')
        i1= this%field_index(TypeIndex(L%CharAt(i,1)))
        i2= this%field_index(TypeIndex(L%CharAt(i,2)))
        if (i1==0 .or. i2==0) cycle
        if (i2>i1) then
            tmp=i1
            i1=i2
            i2=tmp
        end if

        ix=0
        do ii=1, this%nfields
            do jj=1,ii
                ix = ix+1
                if (ii==i1 .and. jj==i2) then
                    cols(i) = ix
                end if
            end do
        end do
    end do

    call L%Clear()

    end subroutine UseString_to_cols

    subroutine UseString_to_TEB(this, S, TEB)
    class(TCMBLikes) :: this
    character(LEN=*), intent(in) :: S
    integer, allocatable, intent(out) :: TEB(:,:)
    Type(TStringList) :: L
    integer tmp,i,i1,i2,ii,jj,ix
    character(LEN=:), pointer :: P

    call L%SetFromString(S,field_names)
    allocate(TEB(2,L%Count), source=0)

    do i=1, L%Count
        P=> L%Item(i)
        if (len(P)/=2) call mpiStop('Invalid C_l order')
        i1= TypeIndex(L%CharAt(i,1))
        i2= TypeIndex(L%CharAt(i,2))
        if (i2>i1) then
            tmp=i1
            i1=i2
            i2=tmp
        end if
        TEB(:,i) = [i1,i2]
    end do
    call L%Clear()

    end subroutine UseString_to_TEB

    subroutine SetTopHatWindows(this)
    class(TCMBLikes) :: this
    integer i

    allocate(this%binWindows(this%pcl_lmin:this%pcl_lmax,1,this%nbins), source=0._mcp)
    !Internal CL are now L(L+1)C_L
    do i=this%pcl_lmin,this%pcl_lmin+this%nbins*this%bin_width -1
        this%binWindows(i,1,(i-this%pcl_lmin)/this%bin_width+1)=real(2*i+1,mcp)/(i*(i+1))
    end do
    do i=this%pcl_lmin+this%nbins*this%bin_width,this%pcl_lmax
        this%binWindows(i,1,this%nbins)=real(2*i+1,mcp)/(i*(i+1))
    end do
    do i=1, this%nbins
        this%binWindows(:,1,i) = this%binWindows(:,1,i)/sum(this%binWindows(:,1,i))
    end do

    end subroutine SetTopHatWindows

    function GetColsFromOrder(this, Order, cols) result(num)
    !Converts string Order = TT TE EE XY... into indices into array of power spectra (and zero if not present)
    class(TCMBLikes) :: this
    character(LEN=*), intent(in) :: Order
    integer, allocatable :: cols(:)
    Type(TStringList) :: Li
    integer i1,i,j,ix,num

    call Li%SetFromString(Order,field_names)
    allocate(cols(this%ncl), source=0)
    ix=0
    do i=1,this%nfields
        do j=1,i
            ix = ix +1
            i1 =Li%IndexOf(this%field_order(i)//this%field_order(j))
            if (i1==-1 .and. i/=j) i1 = Li%IndexOf(this%field_order(j)//this%field_order(i))
            if (i1/=-1) then
                if (cols(ix)>0) call MpiStop('GetColsFromOrder: duplicate CL type')
                cols(ix) = i1
            end if
        end do
    end do
    num = Li%Count

    end function GetColsFromOrder

    subroutine ReadBinWindows(this,Ini)
    !For example:  bin_window_in_order = TT TE TE BB, bin_window_out_order = TT EE TE BB
    !has window file with five columns giving L, TT->TT, TE->EE, TE->TE, BB->BB
    class(TCMBLikes) :: this
    class(TSettingIni) :: Ini
    integer i
    character(LEN=:), allocatable :: filename,  S, Order1, Order2, InLine
    real(mcp), allocatable :: tmp_ar(:)
    integer L, status, norder
    Type(TTextFile) :: F

    filename = Ini%ReadFileName('bin_window_files', NotFoundFail=.true.,relative=.true.)
    Order1 = Ini%Read_String('bin_window_in_order', .true.)
    Order2 = Ini%Read_String_Default('bin_window_out_order', Order1)

    call this%UseString_to_TEB(Order1, this%bin_cols_in)
    call this%UseString_to_cols(Order2, this%bin_cols_out)
    norder = size(this%bin_cols_in,2)
    if (norder/=size(this%bin_cols_out)) &
    & call MpiStop('bin_window_in_order and bin_window_out_order must have same numebr of CL')

    allocate(tmp_ar(norder))
    allocate(this%binWindows(this%pcl_lmin:this%pcl_lmax,norder,this%nbins), source=0._mcp)
    do i=1, this%nbins
        S = FormatString(filename, i)
        do while (F%ReadNextContentLine(S,InLine))
            read(InLine,*, iostat=status) l, tmp_ar
            if (status/=0) call MpiStop('ReadBinWindows: error reading line '//trim(filename))
            if (l>=this%pcl_lmin .and. l <=this%pcl_lmax) then
                this%binWindows(l,:,i) = tmp_ar
            else
                write(*,*) 'WARNING: bin window outside cl_lmin--cl_max range'
            end if
        end do
    end do

    end subroutine ReadBinWindows

    subroutine CMBLikes_ReadIni(this, Ini)
    class(TCMBLikes) :: this
    class(TSettingIni) :: Ini
    real(mcp), dimension(:,:), allocatable, target :: Cov, fullcov
    integer ix, i
    character(LEN=:), allocatable :: S, S_order, covmat_cl
    integer l,j, x,y, clix
    real(mcp) :: asum
    real(mcp), allocatable :: avec(:)
    logical :: cl_fiducial_includes_noise

    if (Ini%TestEqual('dataset_format','CMBLike')) &
    & call MpiStop('CMBLikes dataset_format now CMBLike2 (e.g. covmat are for L(L+1)CL/2pi)')
    if (.not. Ini%TestEqual('dataset_format','CMBLike2',EmptyOK=.true.)) call MpiStop('CMBLikes wrong dataset_format')

    S = Ini%Read_String('fields_use', .true.)
    this%use_field = .false.
    do i=1, len_trim(S)
        if (trim(S(i:i))/='') this%use_field(TypeIndex(S(i:i))) = .true.
    end do

    this%nfields=0
    this%vecsize =0
    this%ncl=0
    this%ncl_used=0

    this%like_approx = Ini%Read_Enumeration('like_approx',Like_Names)
    this%nfields = count(this%use_field)
    ix=0
    this%field_index=0
    do i=1,3
        if (this%use_field(i)) then
            ix=ix+1
            this%field_index(i)=ix
            this%fields(ix) =i
            this%field_order(ix) = field_names(i:i)
        end if
    end do
    this%ncl = (this%nfields*(this%nfields+1))/2

    this%pcl_lmin = Ini%Read_Int('cl_lmin')
    this%pcl_lmax = Ini%Read_Int('cl_lmax')
    this%binned = Ini%Read_Logical('binned')

    if (this%binned) then
        this%nbins = Ini%Read_Int('nbins',0)
    else
        this%nbins = this%pcl_lmax - this%pcl_lmin + 1
        Write(*,*) 'WARNING: Unbinned likelihoods untested in this version'
    end if

    if (bin_test) then
        call MpiStop('bin_test not updated/tested yet')
        if (this%binned) call MpiStop('nbins/=0 with bin_test')
        this%bin_width = Ini%Read_Int('bin_width')
        this%nbins = (this%pcl_lmax-this%pcl_lmin+1)/this%bin_width !Make last bin bigger if not exact multiple
    end if

    if (this%binned) then
        if (bin_test) then
            call this%SetTopHatWindows()
        else
            call this%ReadBinWindows(Ini)
        end if
        this%bin_min=1
        this%bin_max=this%nbins
    else
        this%bin_min=this%pcl_lmin
        this%bin_max=this%pcl_lmax
    end if

    this%vecsize = (this%bin_max-this%bin_min+1)

    this%ncl_hat = Ini%Read_Int('ncl_hat', 1) !only >1 if multiple sims for testing

    allocate(this%ClHat(this%ncl,this%bin_min:this%bin_max, this%ncl_hat))
    S = Ini%ReadFileName('cl_hat_file',NotFoundFail=.true., relative=.true.)
    S_order = Ini%read_String('cl_hat_order', .true.)
    call this%ReadClArr(S, S_order,this%ClHat(:,:,1),this%bin_min)
    do j=2, this%ncl_hat
        !for simulated with multiple realizations with same covariance and noise
        call this%ReadClArr(Ini%ReadFileName(numcat('cl_hat_file',j), relative=.true.), S_order,this%ClHat(:,:,j),this%bin_min)
    end do

    if (this%like_approx /= like_approx_fullsky_exact) then
        allocate(this%ClFiducial(this%ncl,this%bin_min:this%bin_max))
        S =Ini%ReadFileName('cl_fiducial_file',NotFoundFail=.true.,relative=.true.)
        S_order = Ini%Read_String('cl_fiducial_order', .true.)
        call this%ReadClArr(S, S_order,this%ClFiducial,this%bin_min)
    else
        !Exact like
        call Ini%Read('fullsky_exact_fksy', this%fullsky_exact_fksy)
    end if

    allocate(this%ClNoise(this%ncl,this%bin_min:this%bin_max))
    S = Ini%ReadFileName('cl_noise_file',relative=.true.)
    S_order = Ini%read_String('cl_noise_order', .true.)
    call this%ReadClArr(S, S_order,this%ClNoise,this%bin_min)

    this%lensing_recon_ncl = Ini%Read_Int('lensing_recon_ncl', 0)
    if (this%lensing_recon_ncl > 0) call this%ReadLensingReconData(Ini)

    if (this%lensing_recon_ncl >0) then
        allocate(this%cl_lmax(CL_Phi,CL_Phi), source=0)
        this%cl_lmax(CL_Phi,CL_Phi) =  this%pcl_lmax
        this%cl_lmax(CL_Phi,1:this%lensing_recon_ncl) =  this%pcl_lmax
    else
        allocate(this%cl_lmax(CL_B,CL_B), source=0)
    end if

    do i=1,CL_E
        if (this%use_field(i)) then
            do j=1,i
                if (this%use_field(j)) this%cl_lmax(i,j) = this%pcl_lmax
            end do
        end if
    end do
    if (this%use_field(CL_B)) this%cl_lmax(CL_B,CL_B) = this%pcl_lmax

    if (.not. Ini%Read_Logical('cl_hat_includes_noise',.false.)) then
        do j=1,this%ncl_hat
            this%ClHat(:,:,j) =  this%ClHat(:,:,j) + this%ClNoise
        end do
    end if

    cl_fiducial_includes_noise = Ini%Read_Logical('cl_fiducial_includes_noise',.false.)

    if (Ini%HasKey('point_source_cl') .or. Ini%HasKey('beam_modes_file')) &
    & call MpiStop('dataset uses keywords no longer supported')

    allocate(this%ChatM(this%bin_min:this%bin_max, this%ncl_hat))
    allocate(this%NoiseM(this%bin_min:this%bin_max))
    allocate(this%sqrt_fiducial(this%bin_min:this%bin_max))

    !if (this%nbins/=0) then
    !allocate(avec(this%ncl))
    !do i=1, this%nbins
    !    allocate(this%NoiseM(i)%M(this%nfields,this%nfields))
    !    do j=1,this%ncl
    !        avec(j) = sum(this%binWindows(:,i)*(this%ClNoise(j,:)))
    !    end do
    !    call this%ElementsToMatrix(avec, this%NoiseM(i)%M)
    !
    !    if (allocated(this%ClFiducial)) then
    !        allocate(this%sqrt_fiducial(i)%M(this%nfields,this%nfields))
    !        do j=1,this%ncl
    !            avec(j) = sum(this%binWindows(:,i)*this%ClFiducial(j,:))
    !        end do
    !        call this%ElementsToMatrix(avec, this%sqrt_fiducial(i)%M)
    !        if (.not. cl_fiducial_includes_noise) &
    !        & this%sqrt_fiducial(i)%M= this%sqrt_fiducial(i)%M + this%NoiseM(i)%M
    !        call Matrix_Root(this%sqrt_fiducial(i)%M, this%nfields, 0.5_mcp)
    !    end if
    !    do clix =1, this%ncl_hat
    !        allocate(this%ChatM(i,clix)%M(this%nfields,this%nfields))
    !        do j=1,this%ncl
    !            avec(j) = sum(this%binWindows(:,i)*this%ClHat(j,:,clix))
    !        end do
    !        call this%ElementsToMatrix(avec, this%ChatM(i,clix)%M)
    !    end do
    !end do
    !deallocate(avec)
    do l=this%bin_min,this%bin_max
        do clix = 1, this%ncl_hat
            allocate(this%ChatM(l,clix)%M(this%nfields,this%nfields))
            call this%ElementsToMatrix(this%ClHat(:,l,clix), this%ChatM(l,clix)%M)
        end do
        allocate(this%NoiseM(l)%M(this%nfields,this%nfields))
        call this%ElementsToMatrix(this%ClNoise(:,l), this%NoiseM(l)%M)
        if (allocated(this%ClFiducial)) then
            allocate(this%sqrt_fiducial(l)%M(this%nfields,this%nfields))
            if (.not. cl_fiducial_includes_noise) this%ClFiducial(:,l)=this%ClFiducial(:,l)+this%ClNoise(:,l)
            call this%ElementsToMatrix(this%ClFiducial(:,l), this%sqrt_fiducial(l)%M)
            call Matrix_Root(this%sqrt_fiducial(l)%M, this%nfields, 0.5_mcp)
        end if
    end do

    if (this%like_approx /= like_approx_fullsky_exact) then
        call this%ReadCovmat(Ini)
    end if
    if (Ini%HasKey('lowl_exact')) call MpiStop('lowl_exact has been separated out as not currently used')

    call this%TCMBLikelihood%ReadIni(Ini)

    end subroutine CMBLikes_ReadIni

    subroutine ReadCovmat(this, Ini)
    class(TCMBLikes) :: this
    class(TSettingIni) :: Ini
    character(LEN=:), allocatable :: filename, covmat_cl
    integer num_in, ix, binx, biny
    real(mcp), allocatable :: Cov(:,:)
    integer, allocatable :: cl_in_index(:)
    integer lmin_covmat,lmax_covmat, vecsize_in
    integer i, j
    real(mcp) :: covmat_scale = 1

    covmat_cl = Ini%Read_String('covmat_cl', .true.)
    filename = Ini%ReadFileName('covmat_fiducial', NotFoundFail=.true.,relative=.true.)
    call Ini%Read('covmat_scale',covmat_scale)

    call this%UseString_to_cols(covmat_cl, cl_in_index)
    num_in = size(cl_in_index)
    this%ncl_used = count(cl_in_index /=0)
    allocate(this%cl_use_index(this%ncl_used))
    allocate(this%cov_cl_used(this%ncl_used))
    ix = 0
    do i=1, num_in
        if (cl_in_index(i)/=0) then
            ix = ix+1
            this%cl_use_index(ix) = cl_in_index(i)
            this%cov_cl_used(ix) = i
        end if
    end do


    if (this%binned) then
        allocate(Cov(num_in*this%nbins,num_in*this%nbins))
        call File%ReadTextMatrix(filename, Cov)
        allocate(this%inv_covariance(this%nbins*this%ncl_used,this%nbins*this%ncl_used))
        do binx=1, this%nbins
            do biny=1, this%nbins
                this%inv_covariance( (binx-1)*this%ncl_used+1:binx*this%ncl_used, (biny-1)*this%ncl_used+1:biny*this%ncl_used) = &
                & covmat_scale*Cov( (binx-1)*num_in+this%cov_cl_used, (biny-1)*num_in+this%cov_cl_used)
            end do
        end do
        call Matrix_Inverse(this%inv_covariance)
    else
        lmax_covmat = Ini%Read_Int('covmat_lmax')
        lmin_covmat = Ini%Read_Int('covmat_lmin')
        if (lmin_covmat > this%pcl_lmin) call MpiStop('lmin_covmat must be  <= cl_lmin')
        if (lmax_covmat < this%pcl_lmax) call MpiStop('lmax_covmat must be  >= cl_lmax')

        vecsize_in =  (lmax_covmat-lmin_covmat+1)
        if (IsMainMPI()) then
            allocate(Cov(vecsize_in*num_in,vecsize_in*num_in))
            call MatrixSym_Read_Binary(filename, Cov)
            allocate(this%inv_covariance(this%vecsize*this%ncl_used, this%vecsize*this%ncl_used))
            do i=1, this%ncl_used
                do j=1,this%ncl_used
                    this%inv_covariance((i-1)*this%vecsize+1:i*this%vecsize,(j-1)*this%vecsize+1:j*this%vecsize)= covmat_scale &
                    *Cov((this%cov_cl_used(i)-1)*vecsize_in+(this%pcl_lmin-lmin_covmat+1):(this%cov_cl_used(i)-1)*vecsize_in + &
                    (this%pcl_lmax-lmin_covmat+1),  (this%cov_cl_used(j)-1)*vecsize_in+(this%pcl_lmin-lmin_covmat+1): &
                    (this%cov_cl_used(j)-1)*vecsize_in +(this%pcl_lmax-lmin_covmat+1))
                end do
            end do
            deallocate(Cov)
            !if (allocated(this%ClPointsources) .and. this%pointsource_MCMC_modes==0) then
            !    do i=1, this%ncl_used
            !        do j=1,this%ncl_used
            !            do l1= this%pcl_lmin, this%pcl_lmax
            !                do l2= this%pcl_lmin, this%pcl_lmax
            !                    fullcov((i-1)*this%vecsize+ l1 -this%pcl_lmin+1 ,(j-1)*this%vecsize+ l2 -this%pcl_lmin+1 ) = &
            !                    fullcov((i-1)*this%vecsize+ l1 -this%pcl_lmin+1,(j-1)*this%vecsize+ l2 -this%pcl_lmin+1 ) +  &
            !                    this%point_source_error**2*this%ClPointsources(i,l1)*this%ClPointsources(j,l2)
            !                end do
            !            end do
            !        end do
            !    end do
            !end if
            !if (allocated(this%beammodes)) then
            !    do i=this%beam_MCMC_modes+1, nmodes
            !        this%BeamModes(this%pcl_lmin:this%pcl_lmax,i)=this%BeamModes(this%pcl_lmin:this%pcl_lmax,i)*&
            !        & this%ClFiducial(1,this%pcl_lmin:this%pcl_lmax)
            !    end do
            !    do i=this%beam_MCMC_modes+1, nmodes
            !        do l1= this%pcl_lmin, this%pcl_lmax
            !            do l2= this%pcl_lmin, this%pcl_lmax
            !                fullcov(l1 -this%pcl_lmin+1, l2 -this%pcl_lmin+1 ) = &
            !                fullcov(l1 -this%pcl_lmin+1, l2 -this%pcl_lmin+1 ) +  this%BeamModes(l1,i)*this%BeamModes(l2,i)
            !            end do
            !        end do
            !    end do
            !end if

            call Matrix_inverse(this%inv_covariance)
        else !Not mainMPI
            allocate(this%inv_covariance(this%vecsize*this%ncl_used, this%vecsize*this%ncl_used))
        end if !MainMPI
#ifdef MPI
        call MPI_BCAST(this%inv_covariance,Size(this%inv_covariance),MPI_real, 0, MPI_COMM_WORLD, i)
#endif
    end if
    end subroutine ReadCovmat


    subroutine CMBLikes_ReadLensingReconData(this, Ini)
    class(TCMBLikes) :: this
    class(TSettingIni) :: Ini
    character(LEN=:), allocatable :: fname

    if (Feedback > 1) print *,'CMBLikes_ReadLensingReconData'

    this%cl_phi_lmin = Ini%Read_Int('cl_phi_lmin', this%pcl_lmin)
    this%cl_phi_lmax = Ini%Read_Int('cl_phi_lmax', this%pcl_lmax)
    this%phi_like_approx = Ini%Read_Int('phi_like_approx',this%like_approx)

    allocate(this%ClPhiHat(this%lensing_recon_ncl,this%cl_phi_lmin:this%cl_phi_lmax))
    allocate(this%ClPhiNoise(this%lensing_recon_ncl,this%cl_phi_lmin:this%cl_phi_lmax))
    call this%ReadClPhiArr(Ini%ReadFileName('cl_hat_phi_file',relative=.true.),this%ClPhiHat)
    call this%ReadClPhiArr(Ini%ReadFileName('cl_noise_phi_file',relative=.true.),this%ClPhiNoise)

    if (.not. Ini%Read_Logical('cl_hat_includes_noise')) then
        this%ClPhiHat = this%ClPhiHat + this%ClPhiNoise
    end if

    if (this%phi_like_approx /= like_approx_fullsky_exact) then
        fname = Ini%ReadFileName('covmat_phi_fiducial',relative=.true.)
        if (fname /='') then
            allocate(this%phi_inv_covariance(this%cl_phi_lmin:this%cl_phi_lmax,this%cl_phi_lmin:this%cl_phi_lmax))
            call MatrixSym_Read_Binary(fname, this%phi_inv_covariance)
            call Matrix_Inverse(this%phi_inv_covariance)
        end if
    end if

    if (Feedback > 1) print *, 'CMBLikes_ReadLensingReconData done'

    end subroutine CMBLikes_ReadLensingReconData

    subroutine CMBLikes_ReadClPhiArr(this, aname, Cl)
    class(TCMBLikes) :: this
    character(LEN=*), intent(in) :: aname
    real(mcp) :: Cl(:,this%cl_phi_lmin:),tmp_arr(this%lensing_recon_ncl)
    character(LEN=:), allocatable :: InLine
    integer l, ll, status
    Type(TTextFile) :: F

    call F%Open(aname)
    Cl=0
    ll=0
    do while(F%ReadLine(InLine))
        read(InLine,*, iostat=status) l, tmp_arr
        if (status/=0) exit
        if (l>=this%cl_phi_lmin .and. l <=this%cl_phi_lmax) then
            ll=l
            Cl(1,l) = tmp_arr(1)
            if (this%lensing_recon_ncl>1) call MpiStop('CMBLikes_ReadClPhiArr: change for n>1')
        end if
    end do
    if (ll<this%cl_phi_lmax) then
        write(*,*) 'CMBLikes_ReadClPhiArr: C_l file does not go up to phi lmax:', this%cl_phi_lmax
        write (*,*) trim(aname)
        call MpiStop()
    end if
    call F%Close()

    end subroutine CMBLikes_ReadClPhiArr

    subroutine CMBLikes_Transform(this, C, Chat, CfHalf, COffset)
    !Get  C = C_s^{1/2}  U f(this) U^T C_s^{1/2} where C^{-1/2} CHat C^{-1/2} = U this U^T

    !Get  C = C_f^{1/2} C^{-1/2} C^{+1/2} U f(this) U^T C^{+1/2} C^{-1/2} C_f^{1/2} where C^{-1/2} CHat C^{-1/2} = U this U^T

    class(TCMBLikes) :: this
    real(mcp) C(this%nfields,this%nfields)
    real(mcp), intent(in), optional :: COffset(this%nfields,this%nfields)
    real(mcp), intent(in) :: CHat(this%nfields,this%nfields), CfHalf(this%nfields,this%nfields)
    real(mcp) :: U(this%nfields,this%nfields), Rot(this%nfields,this%nfields)
    real(mcp) :: roots(this%nfields)
    real(mcp) :: diag(this%nfields)
    integer i

    if (present(COffset)) then
        U = C + Coffset*C
        call Matrix_Diagonalize(U,Diag,this%nfields)
        Rot= matmul(matmul(transpose(U),CHat+ COffset*C),U)
    else
        U = C
        call Matrix_Diagonalize(U,Diag,this%nfields)
        Rot= matmul(matmul(transpose(U),CHat),U)
    end if
    roots = sqrt(Diag)

    do i=1, this%nfields
        Rot(i,:)=Rot(i,:)/roots(i)
        Rot(:,i)=Rot(:,i)/roots(i)
    end do

    Rot = matmul(U,matmul(Rot,transpose(U)))
    call Matrix_Diagonalize(Rot,Diag,this%nfields)

    Diag = sign(sqrt(2*max(0._mcp,Diag-log(Diag)-1)),Diag-1)
    !want f(this)-1 to save calculating X-X_s

    if (present(COffset)) then
        Rot = MatMul(transpose(U),Rot)
        do i=1, this%nfields
            Rot(i,:)=Rot(i,:)*roots(i)
        end do
        Rot = MatMul(U,Rot)
        call Matrix_Root(C,this%nfields, -0.5_mcp)
        Rot = MatMul(C,Rot)
    end if

    U = matmul(CfHalf,Rot)
    C = U
    do i=1, this%nfields
        C(:,i) = C(:,i)*Diag(i)
    end do
    C = MatMul(C,transpose(U))

    end subroutine CMBLikes_Transform


    subroutine MatrixToElements(this, M, X)
    class(TCMBLikes) :: this
    real(mcp) :: M(this%nfields,this%nfields)
    real(mcp) :: X(this%ncl)
    integer ix,i,j

    ix=0
    do i=1, this%nfields
        do j=1,i
            ix = ix+1
            X(ix) = M(i,j)
        end do
    end do

    end subroutine MatrixToElements

    subroutine MatrixToElementsInt(this, M, X)
    class(TCMBLikes) :: this
    integer, intent(in) :: M(this%nfields,this%nfields)
    integer,intent(out) :: X(this%ncl)
    integer ix,i,j

    ix=0
    do i=1, this%nfields
        do j=1,i
            ix = ix+1
            X(ix) = M(i,j)
        end do
    end do

    end subroutine MatrixToElementsInt


    subroutine ElementsToMatrix(this, X, M)
    class(TCMBLikes) :: this
    real(mcp), intent(out) :: M(this%nfields,this%nfields)
    real(mcp), intent(in) :: X(this%ncl)
    integer ix,i,j

    ix=0
    do i=1, this%nfields
        do j=1,i
            ix = ix+1
            M(i,j) = X(ix)
            M(j,i) = M(i,j)
        end do
    end do

    end subroutine ElementsToMatrix

    function ExactChiSq(this, C,Chat,l)
    class(TCMBLikes) :: this
    real(mcp), intent(in) :: C(this%nfields,this%nfields), Chat(this%nfields,this%nfields)
    integer, intent(in) :: l
    real(mcp) ExactChiSq
    real(mcp) M(this%nfields,this%nfields)

    M = C
    call Matrix_root(M,this%nfields,-0.5_mcp)
    M = matmul(M,matmul(Chat,M))
    ExactChiSq = (2*l+1)*this%fullsky_exact_fksy*(Matrix_Trace(M) - this%nfields - MatrixSym_LogDet(M) )

    end function ExactChiSq

    function CMBLikes_LensRecon_Like(this, Theory) result (chisq)
    class(TCMBLikes) :: this
    Class(TCosmoTheoryPredictions) :: Theory
    real(mcp) vec(this%cl_phi_lmin:this%cl_phi_lmax)
    integer l
    real(mcp) chisq, Cphi, CPhiHat

    if (Feedback > 1) print *,'CMBLikes_LensRecon_CMBLike'

    if (this%bin_width/=1) call MpiStop('CMBLikes_LensRecon_CMBLike: unsupported option')
    if (this%lensing_recon_ncl /=1) call MpiStop('CMBLikes_LensRecon_CMBLike: unsupported ncl')
    !not implemented cross-correlation
    chisq = 0
    do l=this%cl_phi_lmin, this%cl_phi_lmax
        Cphi = Theory%Cls(CL_Phi,CL_phi)%CL(L) + this%ClPhiNoise(1,l)
        CPhihat = this%ClPhihat(1,l)
        if (this%phi_like_approx == like_approx_fullsky_exact) then
            chisq = chisq + (2*l+1)*this%fullsky_exact_fksy*( CPhiHat/CPhi + log(CPhi/CPhiHat) -1)
        else if (this%phi_like_approx == like_approx_fid_gaussian) then
            vec(l) = CPhiHat - CPhi
        else
            call MpiStop('only implemented lensing recon exact like')
        end if
    end do

    if (this%phi_like_approx /= like_approx_fullsky_exact) then
        chisq = chisq + Matrix_QuadForm(this%phi_inv_covariance,vec)
    end if

    if (Feedback > 1) print *,'CMBLikes_LensRecon_CMBLike done'

    end function CMBLikes_LensRecon_Like

    subroutine GetBinnedTheory(this, Theory, C, bin)
    class(TCMBLikes) :: this
    Class(TCosmoTheoryPredictions) :: Theory
    real(mcp) Cls(this%ncl), C(this%nfields,this%nfields)
    integer bin
    integer win_ix,ix_in(2), ix_out

    cls=0
    do win_ix = 1, size(this%bin_cols_in,2)
        ix_in = this%bin_cols_in(:,win_ix)
        if (this%cl_lmax(ix_in(1),ix_in(2))>0) then
            ix_out = this%bin_cols_out(win_ix)
            if (ix_out>0) Cls(ix_out) = Cls(ix_out) + &
            & dot_product(this%BinWindows(:,win_ix,bin),Theory%CLs(ix_in(1),ix_in(2))%CL(this%pcl_lmin:this%pcl_lmax))
        end if
    end do
    call this%ElementsToMatrix(Cls, C)
    end subroutine GetBinnedTheory

    function CMBLikes_LogLike(this, CMB, Theory, DataParams)  result (LogLike)
    real(mcp) logLike
    class(TCMBLikes) :: this
    Class(CMBParams) CMB
    Class(TCosmoTheoryPredictions), target :: Theory
    real(mcp) DataParams(:)
    real(mcp) chisq
    real(mcp) C(this%nfields,this%nfields)
    real(mcp) vecp(this%ncl)
    real(mcp) bigX((this%bin_max-this%bin_min+1)*this%ncl_used)
    integer  i, Ti,Ei,Bi, bin,clix
    logical :: quadratic

    chisq =0

    Ti = this%field_index(1)
    Ei = this%field_index(2)
    Bi = this%field_index(3)

    do clix = 1, this%ncl_hat ! 1 or sum over chi-squareds of simulations
        !removed nuisance parameters for the moment
        !if (this%num_nuisance_parameters/=0 .and. Ti /= 0) then
        !  mode=0
        !  if (this%pointsource_MCMC_modes>0) then
        !    mode=mode+1
        !    cl(this%pcl_lmin:this%pcl_lmax,1) = cl(this%pcl_lmin:this%pcl_lmax,1) + &
        !      this%point_source_error*nuisance_params(mode)*this%ClPointsources(1,this%pcl_lmin:this%pcl_lmax)
        !  end if
        !  beamC=0
        !  do i=1, this%beam_MCMC_modes
        !     mode = mode + 1
        !     BeamC = BeamC + nuisance_params(mode) *cl(this%pcl_lmin:this%pcl_lmax,1)*this%beammodes( this%pcl_lmin:this%pcl_lmax,i )
        !  end do
        !   cl(this%pcl_lmin:this%pcl_lmax,1) = cl(this%pcl_lmin:this%pcl_lmax,1) + beamC
        !
        ! end if

        do bin = this%bin_min, this%bin_max
            if (this%binned .or. bin_test) then
                if (this%like_approx == like_approx_fullsky_exact) call mpiStop('CMBLikes: exact like cannot be binned!')

                call this%GetBinnedTheory(Theory, C, bin)

            else
                C=0
                if (Ti/=0) C(Ti,Ti) = Theory%CLs(CL_T,CL_T)%CL(bin)
                if (Ei/=0) then
                    C(Ei,Ei) = Theory%CLs(CL_E,CL_E)%CL(bin)
                    if (Ti/=0) then
                        C(Ei,Ti) = Theory%CLs(CL_E,CL_T)%CL(bin)
                        C(Ti,Ei) = C(Ei,Ti)
                    end if
                end if
                if (Bi/=0) C(Bi,Bi) =  Theory%CLs(CL_B,CL_B)%CL(bin)
            end if
            C = C + this%NoiseM(bin)%M

            if (this%like_approx == like_approx_HL) then
                call this%Transform(C, this%ChatM(bin,clix)%M, this%sqrt_fiducial(bin)%M)
                call this%MatrixToElements(C, vecp)
                quadratic = .true.
            else if (this%like_approx == like_approx_fid_gaussian) then
                call this%MatrixToElements(C- this%ChatM(bin,clix)%M, vecp)
                quadratic = .true.
            else if (this%like_approx == like_approx_fullsky_exact) then
                quadratic = .false.
                chisq = chisq  + this%ExactChisq(C,this%ChatM(bin,clix)%M,bin)
            else
                call MpiStop('Unknown like_approx')
            end if

            if (quadratic) then
                if (this%binned) then
                    bigX( (bin-1)*this%ncl_used + 1:bin*this%ncl_used) = vecp(this%cl_use_index)
                else
                    bigX( (this%cl_use_index-1)*this%vecsize + bin-this%pcl_lmin+1) = vecp(this%cl_use_index)
                end if
            end if
        end do

        if (quadratic) chisq = chisq + Matrix_QuadForm(this%inv_covariance,BigX)

        !        if (this%lowl_exact) then
        !            chisq = chisq + CMBLikes_lowl_CMBLike(this, cl)
        !        end if

        if (this%lensing_recon_ncl>0) then
            chisq = chisq + this%LensRecon_Like(Theory)
        end if
    end do

    LogLike = chisq/2

    end function CMBLikes_LogLike


    end module CMBLikes
