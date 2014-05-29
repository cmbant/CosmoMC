    !Pseudo-Cl (or other C_l esimator) based likelihood approximation for cut sky with polarization
    !AL Mar 2010 - fixed bug using on E, added support for multiple input simulated Chat for bias testing
    !Apr 2011, added fullsky_exact_fksy
    !Mar 2014  renamed from Planck_like, removed low-L likelihood to separate file
    !Updated file structure, allows for binned HL likelihoods. Everything now in L(L+1)CL/2pi

    module CMBLikes
    use FileUtils
    use settings
    use CosmologyTypes
    use MatrixUtils
    use CosmoTheory
    use Likelihood_Cosmology
    use Interpolation
    implicit none
    private

    logical, parameter :: bin_test = .false.

    Type TSqMatrix
        real(mcp), allocatable :: M(:,:)
    end Type TSqMatrix

    Type, extends(TSqMatrix) :: TIndexedMatrix
        integer, allocatable :: Left_indices(:), Right_indices(:)
    end type TIndexedMatrix

    Type, extends(TIndexedMatrix) :: TCorrectionMatrix
        real(mcp), allocatable :: Fiducial(:)
    contains
    procedure :: InterpProduct
    procedure :: LoadFromFile
    end Type TCorrectionMatrix

    Type LensReconData
        Type(TCorrectionMatrix) N1_matrix_dphi
        Type(TCorrectionMatrix) renorm_matrix
        Type(TCorrectionMatrix) renorm_N1_matrix
        real(mcp), allocatable :: FiducialCl(:,:)
        real(mcp), allocatable :: FiducialPhi(:)
    end Type LensReconData

    integer, parameter :: tot_fields =4

    Type, extends(TCMBLikelihood) :: TCMBLikes
        integer nfields !number of fields
        logical use_field(tot_fields) !Actually directly used for likelihood
        logical required_field(tot_fields) !Used for likelihood, or used for something else like renormalization/derived parameters
        integer field_index(tot_fields) !mapping from 1,2,3 to index actually used
        integer fields(tot_fields) !indices (1=T,2=E,3=B,4=B) of fields to use
        character field_order(tot_fields)
        integer ncl !calculated from above = nfields*(nfields+1)/2
        integer ncl_used !Number of C_l actually used in covariance matrix (others assumed zero)
        integer, allocatable :: cl_use_index(:)
        integer pcl_lmin, pcl_lmax !The range of l in data files (esp. covmat and/or window files)
        integer like_approx
        real(mcp) :: fullsky_exact_fksy = 1 ! only used for testing with exactly fullsky
        integer ncl_hat !1, or more if using multiple simulations
        real(mcp), dimension(:,:), allocatable :: ClFiducial, ClNoise
        real(mcp), dimension(:,:,:), allocatable :: ClHat
        !These are all L(L+1)C_L/2pi

        logical has_lensing
        integer bin_width
        logical binned
        integer nbins, nbins_used
        integer, allocatable :: bin_cols_in(:,:), bin_cols_out(:)
        integer bin_min, bin_max
        integer, allocatable :: cov_cl_used(:)
        real(mcp), dimension(:,:), allocatable :: inv_covariance
        real(mcp), dimension(:,:,:), allocatable :: binWindows

        Type(TSqMatrix), dimension(:), allocatable :: sqrt_fiducial, NoiseM
        Type(TSqMatrix), dimension(:,:), allocatable :: ChatM

        Type(LensReconData) :: Lensing
    contains
    procedure :: LogLike => CMBLikes_LogLike
    procedure :: ReadIni => CMBLikes_ReadIni
    procedure, private :: ReadClArr => CMBLikes_ReadClArr
    procedure, private :: UseString_to_cols
    procedure, private :: UseString_to_TEB
    procedure, private :: Transform => CMBLikes_Transform
    procedure, private :: ExactChisq
    procedure, private :: MatrixToElements
    procedure, private :: MatrixToElementsInt
    procedure, private :: ElementsToMatrix
    procedure, private :: SetTopHatWindows
    procedure, private :: GetColsFromOrder
    procedure, private :: ReadLensing => CMBLikes_ReadLensing
    procedure :: ReadBinWindows
    procedure :: ReadCovmat
    procedure :: GetBinnedTheory
    procedure :: GetObservedTheoryCls !Add foregrounds, or renorm lensing
    end Type TCMBLikes

    character(LEN=tot_fields), parameter :: field_names = 'TEBP' !P is CMB lensing
    integer, parameter :: like_approx_HL=1   !approximation from Hammimeche & Lewis arXiv: 0801.0554
    integer, parameter :: like_approx_fid_gaussian=2 !fiducial fixed covariance matrix, (X-Xhat)^TC^{-1}(X-Xhat)
    integer, parameter :: like_approx_fullsky_exact=3 !ignore all correlations, use exact full sky likelihood function

    character(LEN=Ini_Enumeration_Len), parameter :: &
    & like_Names(3) = [character(Ini_Enumeration_Len)::'HL','gaussian','exact']

    public TCMBLikes
    contains

    function TypeIndex(C)
    character, intent(in) :: C
    integer TypeIndex
    !Get order T, E, B, P -> 1,2,3,4

    TypeIndex = index(field_names,C)
    if (TypeIndex==0) then
        call mpiStop('Invalid C_l part, must be one of: '//field_names)
    end if

    end function TypeIndex

    subroutine CMBLikes_ReadClArr(this, filename, order, Cl, lmin)
    class(TCMBLikes) :: this
    character(LEN=*), intent(in) :: filename, order
    integer, intent(in) :: lmin
    real(mcp) :: Cl(:,lmin:)
    character(LEN=:), allocatable :: tmp, incols
    integer ix, l,ll
    integer, allocatable :: cols(:)
    real(mcp), allocatable :: tmp_ar(:)
    integer status, norder
    Type(TTextFile) :: F

    if (Order=='') then
        incols = File%TopCommentLine(filename)
        if (incols=='') then
            call MpiStop('No column order given for '//filename)
        else
            incols = trim(incols(2:))
        end if
    else
        incols = 'L '//order
    end if
    norder = this%GetColsFromOrder(incols, cols)-1
    allocate(tmp_ar(norder))

    Cl=0
    do while (F%ReadNextContentLine(filename,tmp))
        read(tmp,*, iostat=status) l, tmp_ar
        if (status/=0) call MpiStop('CMBLikes_ReadClArr: error reading line '//trim(filename))
        ll=l
        if (l >=this%bin_min .and. l<=this%bin_max) then
            do ix=1,this%ncl
                if (cols(ix)/=0) Cl(ix,l) = tmp_ar(cols(ix)-1)
            end do
        end if
    end do
    if (ll < this%bin_max) then
        write(*,*) 'CMBLikes_ReadClArr: C_l file does not go up to maximum used:', this%bin_max
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

    subroutine UseString_to_TEB(this, S, TEBP)
    class(TCMBLikes) :: this
    character(LEN=*), intent(in) :: S
    integer, allocatable, intent(out) :: TEBP(:,:)
    Type(TStringList) :: L
    integer tmp,i,i1,i2
    character(LEN=:), pointer :: P

    call L%SetFromString(S,field_names)
    allocate(TEBP(2,L%Count), source=0)

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
        TEBP(:,i) = [i1,i2]
    end do
    call L%Clear()

    end subroutine UseString_to_TEB

    subroutine SetTopHatWindows(this)
    class(TCMBLikes) :: this
    integer i
    !Untested
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

    call Li%SetFromString(Order)
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
    allocate(this%binWindows(this%pcl_lmin:this%pcl_lmax,norder,this%bin_min:this%bin_max), source=0._mcp)
    do i=this%bin_min, this%bin_max
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
    integer ix, i
    character(LEN=:), allocatable :: S, S_order
    integer l,j, clix
    logical :: cl_fiducial_includes_noise, includes_noise

    if (Ini%TestEqual('dataset_format','CMBLike')) &
    & call MpiStop('CMBLikes dataset_format now CMBLike2 (e.g. covmat are for L(L+1)CL/2pi)')
    if (.not. Ini%TestEqual('dataset_format','CMBLike2',EmptyOK=.true.)) call MpiStop('CMBLikes wrong dataset_format')

    S = Ini%Read_String('fields_use', .true.)
    this%use_field = .false.
    do i=1, len_trim(S)
        if (trim(S(i:i))/='') this%use_field(TypeIndex(S(i:i))) = .true.
    end do

    S = Ini%Read_String('fields_required')
    this%required_field = this%use_field
    if (S /= '') then
        do i=1, len_trim(S)
            if (trim(S(i:i))/='') this%required_field(TypeIndex(S(i:i))) = .true.
        end do
    end if

    this%nfields=0
    this%ncl=0
    this%ncl_used=0

    this%like_approx = Ini%Read_Enumeration('like_approx',Like_Names)
    this%nfields = count(this%use_field)
    ix=0
    this%field_index=0
    do i=1,tot_fields
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
        this%bin_min=Ini%Read_Int('use_min',1,min=1,max=this%nbins)
        this%bin_max=Ini%Read_Int('use_max',this%nbins,min=this%bin_min,max=this%nbins)
        if (bin_test) then
            call this%SetTopHatWindows()
        else
            call this%ReadBinWindows(Ini)
        end if
    else
        this%bin_min=Ini%Read_Int('use_min',this%pcl_lmin,min=this%pcl_lmin,max=this%pcl_lmax)
        this%bin_max=Ini%Read_Int('use_max',this%pcl_lmax,min=this%bin_min,max=this%pcl_lmax)
    end if
    this%nbins_used = this%bin_max - this%bin_min + 1

    this%ncl_hat = Ini%Read_Int('ncl_hat', 1) !only >1 if multiple sims for testing

    allocate(this%ClHat(this%ncl,this%bin_min:this%bin_max, this%ncl_hat))
    S = Ini%ReadFileName('cl_hat_file',NotFoundFail=.true., relative=.true.)
    S_order = Ini%read_String('cl_hat_order')
    call this%ReadClArr(S, S_order,this%ClHat(:,:,1),this%bin_min)
    do j=2, this%ncl_hat
        !for simulated with multiple realizations with same covariance and noise
        call this%ReadClArr(Ini%ReadFileName(numcat('cl_hat_file',j), relative=.true.), S_order,this%ClHat(:,:,j),this%bin_min)
    end do

    if (this%like_approx == like_approx_HL) then
        S =Ini%ReadFileName('cl_fiducial_file', NotFoundFail=.true. ,relative=.true.)
        allocate(this%ClFiducial(this%ncl,this%bin_min:this%bin_max))
        S_order = Ini%Read_String('cl_fiducial_order')
        call this%ReadClArr(S, S_order,this%ClFiducial,this%bin_min)
    else if (this%like_approx == like_approx_fullsky_exact) then
        !Exact like
        call Ini%Read('fullsky_exact_fksy', this%fullsky_exact_fksy)
    end if

    includes_noise = Ini%Read_Logical('cl_hat_includes_noise',.false.)
    if (this%like_approx/=like_approx_fid_gaussian .or. includes_noise) then
        S = Ini%ReadFileName('cl_noise_file',relative=.true.,NotFoundFail=.true.)
        allocate(this%ClNoise(this%ncl,this%bin_min:this%bin_max))
        S_order = Ini%Read_String('cl_noise_order')
        call this%ReadClArr(S, S_order,this%ClNoise,this%bin_min)
        if (.not. includes_noise) then
            do j=1,this%ncl_hat
                this%ClHat(:,:,j) =  this%ClHat(:,:,j) + this%ClNoise
            end do
        else if (this%like_approx==like_approx_fid_gaussian) then
            do j=1,this%ncl_hat
                this%ClHat(:,:,j) =  this%ClHat(:,:,j) - this%ClNoise
            end do
            deallocate(this%ClNoise)
        end if
    end if

    allocate(this%cl_lmax(CL_Phi,CL_Phi), source=0)

    do i=1, tot_fields
        if (this%required_field(i)) this%cl_lmax(i,i) = this%pcl_lmax
    end do

    if (this%required_field(CL_T) .and. this%required_field(CL_E)) &
    & this%cl_lmax(CL_E,CL_T) = this%pcl_lmax


    if (Ini%HasKey('point_source_cl') .or. Ini%HasKey('beam_modes_file')) &
    & call MpiStop('dataset uses keywords no longer supported')

    allocate(this%ChatM(this%bin_min:this%bin_max, this%ncl_hat))

    if (this%like_approx /= like_approx_fid_gaussian) then
        cl_fiducial_includes_noise = Ini%Read_Logical('cl_fiducial_includes_noise',.false.)
        allocate(this%NoiseM(this%bin_min:this%bin_max))
        allocate(this%sqrt_fiducial(this%bin_min:this%bin_max))
    end if

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
        if (allocated(this%ClNoise)) then
            allocate(this%NoiseM(l)%M(this%nfields,this%nfields))
            call this%ElementsToMatrix(this%ClNoise(:,l), this%NoiseM(l)%M)
        end if
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

    this%has_lensing  = this%cl_lmax(CL_Phi,CL_Phi) > 0
    if (this%has_lensing) then
        call this%ReadLensing(Ini)
    end if

    call this%TCMBLikelihood%ReadIni(Ini)

    end subroutine CMBLikes_ReadIni

    subroutine CMBLikes_ReadLensing(this, Ini)
    class(TCMBLikes) :: this
    class(TSettingIni) :: Ini
    integer i
    character(LEN=:), allocatable :: S, Comment
    Type(TStringList) :: L, RenormCl
    real(mcp), allocatable :: Fid(:,:)
    integer, allocatable :: ls(:)
    !CMB lensing likelihood

    S = Ini%ReadFileName('lensing_fiducial_cl',relative=.true.)
    if (S/='') then
        call File%LoadTxt(S, Fid, comment = Comment)
        call L%SetFromString(Comment)
        allocate(ls(size(Fid,1)),source=int(Fid(:,1)))
        allocate(this%Lensing%FiducialPhi(0:ls(size(ls))), source=0._mcp)
        this%Lensing%FiducialPhi(ls) = Fid(:,L%IndexOf('PP'))
    end if

    S = Ini%ReadFileName('lensing_N1_matrix_fiducial_dphi',relative=.true.)
    if (S/='') then
        call this%Lensing%N1_matrix_dphi%LoadFromFile(S, this%Lensing%FiducialPhi(1:))
    end if
    S = Ini%Read_String('lensing_renorm_cl')
    if (S/='') then
        if (.not. allocated(Fid)) call MpiStop('lensing_renorm_fields but no lensing_fiducial_cl')
        call RenormCl%SetFromString(S)
        allocate(this%Lensing%FiducialCl(0:ls(size(ls)),RenormCl%Count), source=0._mcp)
        do i=1, RenormCl%Count
            this%Lensing%FiducialCl(ls,i) = Fid(:,L%IndexOf(RenormCl%Item(i)))
        end do

        S = Ini%ReadFileName('lensing_renorm_matrix ',relative=.true.)
        if (S/='') then
            do i=1, RenormCl%Count
                if (RenormCl%Item(i)/='TT') call Mpistop('Pol renormalization not done yet')
                call this%Lensing%renorm_matrix%LoadFromFile(FormatString(S,RenormCl%Item(i)), &
                & this%Lensing%FiducialCL(1:,i))
            end do
        end if
        S = Ini%ReadFileName('lensing_renorm_N1_matrix ',relative=.true.)
        if (S/='') then
            do i=1, RenormCl%Count
                if (RenormCl%Item(i)/='TT') call Mpistop('Pol renormalization not done yet')
                call this%Lensing%renorm_N1_matrix%LoadFromFile(FormatString(S,RenormCl%Item(i)), &
                &this%Lensing%FiducialCL(1:,i))
            end do
        end if
    end if

    end subroutine CMBLikes_ReadLensing

    subroutine ReadCovmat(this, Ini)
    class(TCMBLikes) :: this
    class(TSettingIni) :: Ini
    character(LEN=:), allocatable :: filename, covmat_cl
    integer num_in, ix, binx, biny
    real(mcp), allocatable :: Cov(:,:)
    integer, allocatable :: cl_in_index(:)
    integer vecsize_in
    integer i, j, L1, L2
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
        allocate(this%inv_covariance(this%nbins_used*this%ncl_used,this%nbins_used*this%ncl_used))
        do binx=this%bin_min, this%bin_max
            do biny=this%bin_min, this%bin_max
                this%inv_covariance( (binx-this%bin_min)*this%ncl_used+1:(binx-this%bin_min+1)*this%ncl_used, &
                & (biny-this%bin_min)*this%ncl_used+1:(biny-this%bin_min+1)*this%ncl_used) = &
                & covmat_scale*Cov( (binx-1)*num_in+this%cov_cl_used, (biny-1)*num_in+this%cov_cl_used)
            end do
        end do
        call Matrix_Inverse(this%inv_covariance)
    else
        vecsize_in =  (this%pcl_lmax-this%pcl_lmin+1)
        if (IsMainMPI()) then
            allocate(Cov(vecsize_in*num_in,vecsize_in*num_in))
            call MatrixSym_Read_Binary(filename, Cov)
            allocate(this%inv_covariance(this%nbins_used*this%ncl_used, this%nbins_used*this%ncl_used))
            do i=1, this%ncl_used
                do j=1,this%ncl_used
                    do L1=this%bin_min, this%bin_max
                        do L2 = this%bin_min, this%bin_max
                            !Assume L matrices are ordered the other way in blocks of different L for each CL type
                            this%inv_covariance( (L1-this%bin_min)*this%ncl_used+i, (L2-this%bin_min)*this%ncl_used+j) = &
                            & covmat_scale*Cov((this%cov_cl_used(i)-1)*vecsize_in + L1 - this%pcl_lmin+1, &
                            & (this%cov_cl_used(j)-1)*vecsize_in + L2 - this%pcl_lmin+1)
                        end do
                    end do
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
            allocate(this%inv_covariance(this%nbins_used*this%ncl_used, this%nbins_used*this%ncl_used))
        end if !MainMPI
#ifdef MPI
        call MPI_BCAST(this%inv_covariance,Size(this%inv_covariance),MPI_real, 0, MPI_COMM_WORLD, i)
#endif
    end if
    end subroutine ReadCovmat

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

    subroutine GetBinnedTheory(this, TheoryCls, C, bin)
    class(TCMBLikes) :: this
    Type(TSkyPowerSpectrum) :: TheoryCls(:,:)
    real(mcp) Cls(this%ncl), C(this%nfields,this%nfields)
    integer bin
    integer win_ix,ix_in(2), ix_out

    cls=0
    do win_ix = 1, size(this%bin_cols_in,2)
        ix_in = this%bin_cols_in(:,win_ix)
        if (this%cl_lmax(ix_in(1),ix_in(2))>0) then
            ix_out = this%bin_cols_out(win_ix)
            if (ix_out>0) then
                Cls(ix_out) = Cls(ix_out) + &
                & dot_product(this%BinWindows(:,win_ix,bin),TheoryCls(ix_in(1),ix_in(2))%CL(this%pcl_lmin:this%pcl_lmax))
            end if
        end if
    end do
    call this%ElementsToMatrix(Cls, C)

    end subroutine GetBinnedTheory


    subroutine GetObservedTheoryCls(this, Theory, Cls)
    class(TCMBLikes) :: this
    class(TCosmoTheoryPredictions) :: Theory
    Type(TSkyPowerSpectrum), allocatable, intent(out) :: Cls(:,:)
    real(mcp), allocatable :: N1_Interp(:), renorm(:)
    integer lmax, lmin, lmaxN1, lminN1

    allocate(Cls, source= Theory%Cls)

    if (this%has_lensing) then
        !Correct for normalization and actual N1
        associate(CL=>Cls(CL_Phi,CL_Phi)%CL)
            if (allocated(this%Lensing%N1_matrix_dphi%M)) then
                call this%Lensing%N1_matrix_dphi%InterpProduct(CL(1:), N1_interp, lminN1, lmaxN1, this%pcl_lmax, sub_fiducial = .true.)
            end if

            if (allocated(this%Lensing%renorm_matrix%M)) then
                call this%Lensing%renorm_matrix%InterpProduct(Theory%Cls(CL_T,CL_T)%Cl(1:), renorm, lmin, lmax, &
                & this%pcl_lmax, sub_fiducial = .true.)
                CL(lmin:lmax) = CL(lmin:lmax)*(1 + 2*renorm)
            end if

            if (allocated(N1_interp)) CL(lminN1:lmaxN1) =  CL(lminN1:lmaxN1) + N1_interp

            if (allocated(this%Lensing%renorm_N1_matrix%M)) then
                call this%Lensing%renorm_N1_matrix%InterpProduct(Theory%Cls(CL_T,CL_T)%Cl(1:), &
                & N1_interp, lminN1, lmaxN1, this%pcl_lmax, sub_fiducial = .true.)
                CL(lminN1:lmaxN1) =  CL(lminN1:lmaxN1) + 2*N1_interp
            end if
            end associate
    end if

    end subroutine GetObservedTheoryCls

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
    integer  Ti,Ei,Bi,Phii, bin,clix
    logical :: quadratic
    Type(TSkyPowerSpectrum), allocatable :: TheoryCls(:,:)

    chisq =0

    Ti = this%field_index(1)
    Ei = this%field_index(2)
    Bi = this%field_index(3)
    Phii = this%field_index(4)
    call this%GetObservedTheoryCls(Theory, TheoryCls)

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
                call this%GetBinnedTheory(TheoryCls, C, bin)
            else
                C=0
                if (Ti/=0) C(Ti,Ti) = TheoryCLs(CL_T,CL_T)%CL(bin)
                if (Ei/=0) then
                    C(Ei,Ei) = TheoryCLs(CL_E,CL_E)%CL(bin)
                    if (Ti/=0) then
                        C(Ei,Ti) = TheoryCLs(CL_E,CL_T)%CL(bin)
                        C(Ti,Ei) = C(Ei,Ti)
                    end if
                end if
                if (Bi/=0) C(Bi,Bi) =  TheoryCLs(CL_B,CL_B)%CL(bin)
                if (Phii/=0) then
                    C(Phii,Phii) = TheoryCLs(CL_Phi,CL_Phi)%CL(bin)
                    !No cross for the moment
                    !if (Ti/=0) C(Phii,Ti) = TheoryCLs(CL_Phi,CL_T)%CL(bin)
                end if
            end if

            if (allocated(this%NoiseM)) then
                C = C + this%NoiseM(bin)%M
            end if

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
                bigX( (bin-this%bin_min)*this%ncl_used + 1:(bin-this%bin_min+1)*this%ncl_used) = vecp(this%cl_use_index)
            end if
        end do

        if (quadratic) chisq = chisq + Matrix_QuadForm(this%inv_covariance,BigX)

    end do

    LogLike = chisq/2

    end function CMBLikes_LogLike


    subroutine LoadFromFile(this, fname, fidCL)
    class(TCorrectionMatrix) :: this
    character(LEN=*), intent(in) :: fname
    real(mcp), allocatable :: Mat(:,:)
    real(mcp), intent(in) :: fidCL(:)
    integer m,n

    call File%LoadTxt(fname,Mat, m,n)
    allocate(this%M(m-1, n-1), source = Mat(2:m,2:n))
    allocate(this%left_indices(m-1), this%right_indices(n-1))
    this%left_indices = int(Mat(2:m,1))
    this%right_indices = int(Mat(1,2:n))

    allocate(this%Fiducial(m-1))
    this%Fiducial = matmul(this%M, fidCL(this%Right_indices))

    end subroutine LoadFromFile

    subroutine InterpProduct(this, RightInput, InterpOutput, lmin, lmax, lmaxmax, sub_fiducial)
    class (TCorrectionMatrix) :: this
    real(mcp), intent(in) :: RightInput(:)
    real(mcp), allocatable, intent(out) :: InterpOutput(:)
    integer, intent(out) :: lmin, lmax
    integer, intent(in) :: lmaxmax
    logical sub_fiducial
    Type(TCubicSpline) :: Interp
    real(mcp), allocatable :: LeftOutput(:)

    allocate(LeftOutput(size(this%Left_Indices)))
    Leftoutput = matmul(this%M, RightInput(this%Right_indices))
    if (sub_fiducial) LeftOutput = LeftOutput - this%Fiducial
    call Interp%Init(this%Left_indices, LeftOutput)
    lmin = this%Left_indices(1)
    lmax = min(lmaxmax,this%Left_indices(size(this%Left_Indices)))
    allocate(InterpOutput(lmin:lmax))
    call Interp%Array(lmin, lmax, InterpOutput)

    end subroutine InterpProduct

    end module CMBLikes
