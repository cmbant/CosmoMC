    module CMBlikes
    !Pseudo-Cl (or other C_l esimator) based likelihood approximation for cut sky with polarization
    !AL Mar 2010 - fixed bug using on E, added support for multiple input simulated Chat for bias testing
    !Apr 2011, added fullsky_exact_fksy
    !Mar 2014  renamed from Planck_like, removed low-L likelihood to separate file
    !Updated file structure, allows for binned HL likelihoods. Everything now in L(L+1)CL/2pi
    !2014 linear correction bin windows (e.g. for Planck lensing likelihood)
    !Nov 2014: added generalization for more general map cross spectra
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

    character(LEN=*), parameter :: cross_separators = 'x'

    Type TSqMatrix
        real(mcp), allocatable :: M(:,:)
    end Type TSqMatrix

    Type TBinWindows
        integer :: lmin, lmax
        real(mcp), dimension(:,:,:), allocatable :: W
        integer, allocatable :: bin_cols_in(:,:), bin_cols_out(:)
        Type(TSkyPowerSpectrum), allocatable :: fixCls(:,:) !to force particular spectra to be fixed (if allocated)
    contains
    procedure :: bin => TBinWindows_bin
    end type

    Type, extends(TSkyPowerSpectrum) :: TMapCrossPowerSpectrum
        integer map_i, map_j
        integer theory_i, theory_j
    end type TMapCrossPowerSpectrum

    integer, parameter :: tot_theory_fields = len(CMB_CL_Fields)

    Type, extends(TCMBLikelihood) :: TCMBLikes
        integer nmaps !number of maps used in the likelihood
        integer nmaps_required !number of maps used before binning or corrections
        logical, allocatable :: use_map(:) !which maps are used directly in likelihood (others ignored)
        logical, allocatable :: require_map(:) !which maps are used directly in likelihood (others ignored)
        logical required_theory_field(tot_theory_fields)
        !Used for likelihood, or used for something else like renormalization/derived parameters
        logical :: has_map_names = .false.
        Type(TStringList) :: map_names !strings labelling each distinct map, e.g. T143, B353, .. or just T E etc.
        integer, allocatable :: map_fields(:) !the theory field that each map is measuring

        integer, allocatable :: map_used_index(:) !for each map_name, zero if not used, otherwise the index into the used map list
        integer, allocatable :: map_required_index(:)  !same for required maps
        integer, allocatable :: required_order(:)
        logical :: has_foregrounds = .false.
        Type(TStringList) :: used_map_order !names of the maps actually used
        integer ncl !calculated from above = nmaps*(nmaps+1)/2
        integer ncl_used !Number of C_l actually used in covariance matrix (others assumed zero)
        integer, allocatable :: cl_use_index(:)
        integer pcl_lmin, pcl_lmax !The range of l in data files (esp. covmat and/or window files)
        integer like_approx
        real(mcp) :: fullsky_exact_fksy = 1 ! only used for testing with exactly fullsky
        integer :: calibration_index = 0 !if non-zero, index into nuisance parameter array of global calibration of TEB
        real(mcp), dimension(:,:), allocatable :: ClFiducial, ClNoise, ClHat
        !These are all L(L+1)C_L/2pi

        real(mcp), dimension(:,:), allocatable :: inv_covariance

        logical has_lensing

        !Binning parameters
        integer bin_width
        logical binned
        integer nbins, nbins_used
        integer bin_min, bin_max
        Type(TBinWindows) binWindows

        !For correcting theory bandpowers using linear correction (which is zero at fiducial model):
        !   CL_bin --> CL_bin + dot(binCorrectionWindows, CL) - FiducialCorrection
        !Could be incorporated in bandpowers and bin window, but allows cleaner separation of effect of correction for testing
        Type(TBinWindows) binCorrectionWindows
        real(mcp), allocatable :: FiducialCorrection(:,:)

        Type(TSqMatrix), dimension(:), allocatable :: sqrt_fiducial, NoiseM
        Type(TSqMatrix), dimension(:), allocatable :: ChatM
        class(TMapCrossPowerSpectrum), allocatable :: MapCls(:,:)
    contains
    procedure :: LogLike => CMBLikes_LogLike
    procedure :: ReadIni => CMBLikes_ReadIni
    procedure, private :: ReadClArr => CMBLikes_ReadClArr
    procedure, private :: UseString_to_cols
    procedure, private :: UseString_to_Cl_i_j
    procedure, private :: Transform => CMBLikes_Transform
    procedure, private :: ExactChisq
    procedure, private :: MatrixToElements
    procedure, private :: MatrixToElementsInt
    procedure, private :: ElementsToMatrix
    !    procedure, private :: SetTopHatWindows
    procedure, private :: GetColsFromOrder
    procedure, private :: Cl_used_i_j_name
    procedure, private :: Cl_i_j_name
    procedure, private :: MapPair_to_Theory_i_j
    procedure, private :: PairStringToUsedMapIndices
    procedure, private :: PairStringToMapIndices
    procedure, nopass :: TypeIndex
    procedure :: ReadBinWindows
    procedure :: ReadCovmat
    procedure :: InitMapCls
    procedure :: GetBinnedMapCls
    procedure :: GetTheoryMapCls !take theory calcualtion and add foregrounds etc using sub below
    procedure :: AdaptTheoryForMaps !call AddForegrounds etc for calibrations, beams, etc
    procedure :: AddForegrounds
    procedure :: WriteLikelihoodData => CMBLikes_WriteLikelihoodData
    end Type TCMBLikes

    integer, parameter :: like_approx_HL=1   !approximation from Hammimeche & Lewis arXiv: 0801.0554
    integer, parameter :: like_approx_fid_gaussian=2 !fiducial fixed covariance matrix, (X-Xhat)^TC^{-1}(X-Xhat)
    integer, parameter :: like_approx_fullsky_exact=3 !ignore all correlations, use exact full sky likelihood function

    character(LEN=Ini_Enumeration_Len), parameter :: &
        & CMBLikes_like_Names(3) = [character(Ini_Enumeration_Len)::'HL','gaussian','exact']

    public TCMBLikes, TMapCrossPowerSpectrum, TBinWindows, CMBLikes_like_Names
    contains

    function TypeIndex(C)
    character, intent(in) :: C
    integer TypeIndex
    !Get order T, E, B, P -> 1,2,3,4

    TypeIndex = index(CMB_CL_Fields,C)
    if (TypeIndex==0) then
        call mpiStop('Invalid C_l part, must be one of: '//CMB_CL_Fields)
    end if

    end function TypeIndex

    subroutine CMBLikes_ReadClArr(this, Ini, basename,  Cl, hasKey)
    class(TCMBLikes) :: this
    class(TSettingIni), intent(in) :: Ini
    character(LEN=*), intent(in) :: basename
    logical, intent(inout), optional :: hasKey
    real(mcp), allocatable,intent(out) :: Cl(:,:)
    character(LEN=:), allocatable :: tmp, incols, order, filename
    integer ix, l,ll
    integer, allocatable :: cols(:)
    real(mcp), allocatable :: tmp_ar(:)
    integer status, norder
    Type(TTextFile) :: F

    filename =Ini%ReadRelativeFileName(basename//'_file', NotFoundFail=.not. present(hasKey))
    if (present(hasKey)) then
        hasKey = filename /=''
        if (.not. hasKey) return
    end if
    allocate(Cl(this%ncl,this%bin_min:this%bin_max))
    Order = Ini%Read_String(basename//'_order')

    if (Order=='') then
        incols = File%LastTopComment(filename)
        if (incols=='') call MpiStop('No column order given for '//filename)
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

    subroutine PairStringToMapIndices(this,S, i1,i2)
    class(TCMBLikes) :: this
    character(LEN=*), intent(in) :: S
    integer, intent(out) :: i1,i2
    integer ix

    if (len(S)==2) then
        if (this%has_map_names) call MpiStop('CMBlikes: CL names must use MAP1xMAP2 names')
        i1= this%map_names%indexOf(S(1:1))
        i2= this%map_names%indexOf(S(2:2))
    else
        ix = scan(S, cross_separators)
        if (ix==0) call MpiStop('CMBLikes: invalid spectrum name '//S)
        i1 = this%map_names%indexOf(S(1:ix-1))
        i2 = this%map_names%indexOf(S(ix+1:))
    end if
    if (i1==-1 .or. i2==-1) call MpiStop('CMBLikes: unrecognised map name '//S)

    end subroutine PairStringToMapIndices

    subroutine PairStringToUsedMapIndices(this,used_index, S, i1,i2)
    class(TCMBLikes) :: this
    character(LEN=*), intent(in) :: S
    integer, intent(in) :: used_index(:)
    integer, intent(out) :: i1,i2
    integer tmp

    call this%PairStringToMapIndices(S,i1,i2)
    i1 = used_index(i1)
    i2 = used_index(i2)
    if (i2>i1) then
        tmp=i1
        i1=i2
        i2=tmp
    end if

    end subroutine PairStringToUsedMapIndices


    subroutine UseString_to_cols(this, S, cols)
    class(TCMBLikes) :: this
    character(LEN=*), intent(in) :: S
    integer, allocatable, intent(out) :: cols(:)
    integer i,i1,i2,ii,jj,ix
    integer, allocatable :: cl_i_j(:,:)

    call this%UseString_to_Cl_i_j(S, this%map_used_index, cl_i_j)

    allocate(cols(size(cl_i_j,2)), source=0)

    do i=1, size(cl_i_j,2)
        i1 = cl_i_j(1,i)
        i2 = cl_i_j(2,i)
        if (i1==0 .or. i2==0) cycle
        ix=0
        do ii=1, this%nmaps
            do jj=1,ii
                ix = ix+1
                if (ii==i1 .and. jj==i2) then
                    cols(i) = ix
                end if
            end do
        end do
    end do

    end subroutine UseString_to_cols

    subroutine UseString_to_Cl_i_j(this, S, used_index, cl_i_j)
    class(TCMBLikes) :: this
    character(LEN=*), intent(in) :: S
    integer, allocatable, intent(out) :: cl_i_j(:,:)
    integer, intent(in) :: used_index(:)
    Type(TStringList) :: L
    integer i,i1,i2
    character(LEN=:), pointer :: P

    call L%SetFromString(S)
    allocate(cl_i_j(2,L%Count), source=0)

    do i=1, L%Count
        P=> L%Item(i)
        call this%PairStringToUsedMapIndices(used_index,P,i1,i2)
        cl_i_j(:,i) = [i1,i2]
    end do
    call L%Clear()

    end subroutine UseString_to_Cl_i_j


    subroutine MapPair_to_Theory_i_j(this, order, i1,i2, i, j)
    class(TCMBLikes) :: this
    integer, intent(in) :: order(:)
    integer, intent(in) :: i1,i2
    integer, intent(out) :: i,j
    integer tmp

    i = this%map_fields(order(i1))
    j = this%map_fields(order(i2))
    if (j>i) then
        tmp = i
        i = j
        j=tmp
    end if

    end subroutine MapPair_to_Theory_i_j

    !subroutine SetTopHatWindows(this)
    !class(TCMBLikes) :: this
    !integer i
    !!Untested
    !allocate(this%binWindows(this%pcl_lmin:this%pcl_lmax,1,this%nbins), source=0._mcp)
    !!Internal CL are now L(L+1)C_L
    !do i=this%pcl_lmin,this%pcl_lmin+this%nbins*this%bin_width -1
    !    this%binWindows(i,1,(i-this%pcl_lmin)/this%bin_width+1)=real(2*i+1,mcp)/(i*(i+1))
    !end do
    !do i=this%pcl_lmin+this%nbins*this%bin_width,this%pcl_lmax
    !    this%binWindows(i,1,this%nbins)=real(2*i+1,mcp)/(i*(i+1))
    !end do
    !do i=1, this%nbins
    !    this%binWindows(:,1,i) = this%binWindows(:,1,i)/sum(this%binWindows(:,1,i))
    !end do
    !
    !end subroutine SetTopHatWindows

    function Cl_used_i_j_name(this,i,j) result(ClName)
    class(TCMBLikes) :: this
    integer, intent(in) :: i,j
    character(LEN=:), allocatable :: ClName

    ClName = this%CL_i_j_name(this%used_map_order,i,j)

    end function Cl_used_i_j_name

    function Cl_i_j_name(this,names,i,j) result(ClName)
    class(TCMBLikes) :: this
    class(TStringList), intent(in) :: names
    character(LEN=:), allocatable :: ClName
    integer, intent(in) :: i,j
    character(LEN=:), allocatable :: name1,name2

    name1 = names%Item(i)
    name2 = names%Item(j)
    if (this%has_map_names) then
        ClName = name1//cross_separators(1:1)//name2
    else
        ClName = name1//name2
    end if

    end function Cl_i_j_name

    function GetColsFromOrder(this, Order, cols) result(num)
    !Converts string Order = TT TE EE XY... or AAAxBBB AAAxCCC BBxCC into indices into array of power spectra (and zero if not present)
    class(TCMBLikes) :: this
    character(LEN=*), intent(in) :: Order
    integer, allocatable :: cols(:)
    Type(TStringList) :: Li
    integer i1,i,j,ix,num

    call Li%SetFromString(Order)
    allocate(cols(this%ncl), source=0)
    ix=0
    do i=1,this%nmaps
        do j=1,i
            ix = ix +1
            i1 =Li%IndexOf(this%Cl_used_i_j_name(i,j))
            if (i1==-1 .and. i/=j) i1 = Li%IndexOf(this%Cl_used_i_j_name(j,i))
            if (i1/=-1) then
                if (cols(ix)>0) call MpiStop('GetColsFromOrder: duplicate CL type')
                cols(ix) = i1
            end if
        end do
    end do
    num = Li%Count

    end function GetColsFromOrder

    subroutine ReadBinWindows(this,Ini, bin_type, binWindows)
    !For example:  bin_window_in_order = TT TE TE BB, bin_window_out_order = TT EE TE BB
    !has window file with five columns giving L, TT->TT, TE->EE, TE->TE, BB->BB
    class(TCMBLikes) :: this
    class(TSettingIni) :: Ini
    character(LEN=*), intent(in) :: bin_type
    class(TBinWindows) binWindows
    integer i, j
    character(LEN=:), allocatable :: filename,  S, Order1, Order2, InLine
    real(mcp), allocatable :: tmp_ar(:)
    integer L, status, norder
    Type(TTextFile) :: F
    logical Err
    integer, allocatable :: MapPairsFile(:,:), MapPairsUse(:,:), fixed_index(:)
    integer nfixed

    filename = Ini%ReadRelativeFileName(bin_type // '_files', NotFoundFail=.true.)
    Order1 = Ini%Read_String(bin_type // '_in_order', .true.)
    Order2 = Ini%Read_String_Default(bin_type // '_out_order', Order1)
    call this%UseString_to_CL_i_j(Order1, this%map_required_index, binWindows%bin_cols_in)
    call this%UseString_to_cols(Order2, binWindows%bin_cols_out)
    norder = size(binWindows%bin_cols_in,2)
    if (norder/=size(binWindows%bin_cols_out)) &
        & call MpiStop(bin_type // '_in_order and '//bin_type // '_out_order must have same numebr of CL')

    allocate(tmp_ar(norder))
    binWindows%lmin = this%pcl_lmin
    binWindows%lmax = this%pcl_lmax
    allocate(binWindows%W(this%pcl_lmin:this%pcl_lmax,norder,this%bin_min:this%bin_max), source=0._mcp)
    do i=this%bin_min, this%bin_max
        S = FormatString(filename, i)
        Err = .false.
        do while (F%ReadNextContentLine(S,InLine))
            read(InLine,*, iostat=status) l, tmp_ar
            if (status/=0) call MpiStop('ReadBinWindows: error reading line '//trim(filename))
            if (l>=this%pcl_lmin .and. l <=this%pcl_lmax) then
                binWindows%W(l,:,i) = tmp_ar
            else
                Err = Err .or. any(tmp_ar/=0)
            end if
        end do
        if (Err) then
            write(*,*) FormatString( 'WARNING: %s %u outside cl_lmin-cl_max range: %s', bin_type, i, S)
        end if
    end do

    S = Ini%ReadRelativeFileName(bin_type//'_fix_cl_file')
    if (S/='') then
        !Want to force some of the component CL to a fixed model
        Order1 = Ini%Read_String(bin_type // '_fix_cl_file_order')
        if (order1=='') then
            Order1 = File%LastTopComment(S)
            if (Order1=='') call MpiStop('No column order given for '//S)
            do while(trim(Order1(1:1))=='' .or. Order1(1:1)=='L')
                Order1= Order1(2:)
            end do
        end if
        call this%UseString_to_Cl_i_j(Order1, this%map_required_index, MapPairsFile)
        norder = size(MapPairsFile,2)
        deallocate(tmp_ar)
        allocate(tmp_ar(norder))
        Order2 = Ini%Read_String(bin_type // '_fix_cl', .true.)
        call this%UseString_to_Cl_i_j(Order2,this%map_required_index,  MapPairsUse)
        nfixed = size(MapPairsUse,2)
        allocate(fixed_index(nfixed))
        allocate(binWindows%fixCls(this%nmaps_required,this%nmaps_required))
        do i=1,nfixed
            allocate(binWindows%fixCls(MapPairsUse(1,i),MapPairsUse(2,i))%CL(this%pcl_lmin:this%pcl_lmax), source=0._mcp)
            do j=1, norder+1
                if (j>norder) call MpiStop('ReadBinWindows: fix_cl uses CL not in the fix_cl_file')
                if (all(MapPairsFile(:,j)==MapPairsUse(:,i))) then
                    fixed_index(i) = j
                    exit
                end if
            end do
        end do
        do while (F%ReadNextContentLine(S,inLine))
            read(InLine,*, iostat=status) l, tmp_ar
            if (status/=0) call MpiStop(' error fix_cl_file line '//trim(filename))
            if (l >=this%pcl_lmin .and. l<=this%pcl_lmax) then
                do i=1, nfixed
                    binWindows%fixCls(MapPairsUse(1,i),MapPairsUse(2,i))%CL(l) = tmp_ar(fixed_index(i))
                end do
            end if
        end do
    end if

    end subroutine ReadBinWindows

    subroutine CMBLikes_ReadIni(this, Ini)
    class(TCMBLikes) :: this
    class(TSettingIni) :: Ini
    integer ix, i
    character(LEN=:), allocatable :: S
    integer l,j
    logical :: cl_fiducial_includes_noise, includes_noise
    Type(TStringList) :: map_fields, fields_use, maps_use
    logical use_theory_field(tot_theory_fields)
    character(LEN=:), allocatable :: tmp
    logical :: hasKey

    if (Ini%TestEqual('dataset_format','CMBLike')) &
        & call MpiStop('CMBLikes dataset_format now CMBLike2 (e.g. covmat are for L(L+1)CL/2pi)')
    if (.not. Ini%TestEqual('dataset_format','CMBLike2',EmptyOK=.true.)) call MpiStop('CMBLikes wrong dataset_format')

    S = Ini%Read_String('map_names')
    this%has_map_names = S/=''
    if (this%has_map_names) then
        !E.g. have multiple frequencies for given measurement
        call this%map_names%SetFromString(S)
        S = Ini%Read_String('map_fields',.true.)
        call map_fields%SetFromString(S)
        if (map_fields%Count/=this%map_names%Count) &
            call MpiStop('CMBLikes: number of map_fields does not match map_names')
        allocate(this%map_fields(this%map_names%Count))
        do i=1, map_fields%Count
            this%map_fields(i) = this%TypeIndex(map_fields%CharAt(i,1))
        end do
    else
        !map fields are just directly T E B or P (one each at most)
        allocate(this%map_fields(tot_theory_fields))
        do i=1, tot_theory_fields
            tmp = CMB_CL_Fields(i:i) !just avoid ifort 15.01 bug on adding directly
            call this%map_names%Add(tmp)
            this%map_fields(i) =i
        end do
    end if

    S = Ini%Read_String('fields_use')
    if (S/='') then
        use_theory_field = .false.
        call fields_use%SetFromString(S)
        do i=1, fields_use%count
            use_theory_field(this%TypeIndex(fields_use%CharAt(i,1))) = .true.
        end do
    else
        if (.not. this%has_map_names) call MpiStop('CMBlikes: must have fields_use or map_names')
        use_theory_field = .true.
    end if

    allocate(this%use_map(this%map_names%Count))
    S = Ini%Read_String('maps_use')
    if (S/='') then
        if (any(.not. use_theory_field)) &
            print *, 'CMBlikes WARNING: maps_use overrides fields_use'
        this%use_map = .false.
        call maps_use%SetFromString(S)
        do i=1,maps_use%Count
            j = this%map_names%IndexOf(maps_use%Item(i))
            if (j/=-1) then
                this%use_map(j) = .true.
            else
                call MpiStop('CMBlikes: maps_use item not found - '//maps_use%Item(i))
            end if
        end do
    else
        this%use_map = .true.
        do i=1,this%map_names%count
            this%use_map(i) = use_theory_field(this%map_fields(i))
        end do
    end if

    allocate(this%require_map(this%map_names%Count))
    this%require_map = this%use_map
    !Bandpowers can depend on more fields than are actually used in likelihood
    !e.g. for correcting leakage or other linear corrections
    if (this%has_map_names) then
        S = Ini%Read_String('maps_required')
        if (Ini%HasKey('fields_required')) call MpiStop('CMBLikes: use maps_required not fields_required')
    else
        S = Ini%Read_String('fields_required')
    end if
    if (S/='') then
        call maps_use%SetFromString(S)
        do i=1,maps_use%Count
            j = this%map_names%IndexOf(maps_use%Item(i))
            if (j/=-1) then
                this%require_map(j) = .true.
            else
                call MpiStop('CMBlikes: required item not found - '//maps_use%Item(i))
            end if
        end do
    end if
    this%required_theory_field = .false.
    do i=1, this%map_names%Count
        if (this%require_map(i)) then
            this%required_theory_field(this%map_fields(i)) = .true.
        end if
    end do

    this%ncl_used=0

    this%like_approx = Ini%Read_Enumeration('like_approx',CMBLikes_like_Names)
    this%nmaps = count(this%use_map)
    this%nmaps_required = count(this%require_map)

    allocate(this%required_order(this%nmaps_required))
    allocate(this%map_required_index(this%map_names%Count), source=0)
    ix =0
    do i=1, this%map_names%Count
        if (this%require_map(i)) then
            ix=ix+1
            this%map_required_index(i)=ix
            this%required_order(ix) = i
        end if
    end do

    allocate(this%map_used_index(this%map_names%Count), source=0)
    ix=0
    do i=1, this%map_names%Count
        if (this%use_map(i)) then
            ix=ix+1
            this%map_used_index(i)=ix
            call this%used_map_order%Add(this%map_names%Item(i))
        end if
    end do
    this%ncl = (this%nmaps*(this%nmaps+1))/2

    this%pcl_lmin = Ini%Read_Int('cl_lmin')
    this%pcl_lmax = Ini%Read_Int('cl_lmax')
    this%binned = Ini%Read_Logical('binned')

    if (this%binned) then
        this%nbins = Ini%Read_Int('nbins',0)
    else
        this%nbins = this%pcl_lmax - this%pcl_lmin + 1
        if (this%like_approx/=like_approx_fullsky_exact) &
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
            call MpiStop('bin_test not implemented')
            !           call this%SetTopHatWindows()
        else
            call this%ReadBinWindows(Ini, 'bin_window', this%binWindows)
        end if
    else
        this%bin_min=Ini%Read_Int('use_min',this%pcl_lmin,min=this%pcl_lmin,max=this%pcl_lmax)
        this%bin_max=Ini%Read_Int('use_max',this%pcl_lmax,min=this%bin_min,max=this%pcl_lmax)
    end if
    this%nbins_used = this%bin_max - this%bin_min + 1

    call this%ReadClArr(Ini, 'cl_hat', this%ClHat)

    if (this%like_approx == like_approx_HL) then
        call this%ReadClArr(Ini,'cl_fiducial',this%ClFiducial)
    else if (this%like_approx == like_approx_fullsky_exact) then
        !Exact like
        call Ini%Read('fullsky_exact_fksy', this%fullsky_exact_fksy)
    end if

    includes_noise = Ini%Read_Logical('cl_hat_includes_noise',.false.)
    if (this%like_approx/=like_approx_fid_gaussian .or. includes_noise) then
        call this%ReadClArr(Ini, 'cl_noise',this%ClNoise)
        if (.not. includes_noise) then
            this%ClHat =  this%ClHat + this%ClNoise
        else if (this%like_approx==like_approx_fid_gaussian) then
            this%ClHat =  this%ClHat - this%ClNoise
            deallocate(this%ClNoise)
        end if
    end if

    allocate(this%cl_lmax(tot_theory_fields,tot_theory_fields), source=0)

    do i=1, tot_theory_fields
        if (this%required_theory_field(i)) this%cl_lmax(i,i) = this%pcl_lmax
    end do

    if (this%required_theory_field(CL_T) .and. this%required_theory_field(CL_E)) &
        & this%cl_lmax(CL_E,CL_T) = this%pcl_lmax


    if (Ini%HasKey('point_source_cl') .or. Ini%HasKey('beam_modes_file')) &
        & call MpiStop('dataset uses keywords no longer supported')

    allocate(this%ChatM(this%bin_min:this%bin_max))

    if (this%like_approx /= like_approx_fid_gaussian) then
        cl_fiducial_includes_noise = Ini%Read_Logical('cl_fiducial_includes_noise',.false.)
        allocate(this%NoiseM(this%bin_min:this%bin_max))
        allocate(this%sqrt_fiducial(this%bin_min:this%bin_max))
    end if

    !if (this%nbins/=0) then
    !allocate(avec(this%ncl))
    !do i=1, this%nbins
    !    allocate(this%NoiseM(i)%M(this%nmaps,this%nmaps))
    !    do j=1,this%ncl
    !        avec(j) = sum(this%binWindows(:,i)*(this%ClNoise(j,:)))
    !    end do
    !    call this%ElementsToMatrix(avec, this%NoiseM(i)%M)
    !
    !    if (allocated(this%ClFiducial)) then
    !        allocate(this%sqrt_fiducial(i)%M(this%nmaps,this%nmaps))
    !        do j=1,this%ncl
    !            avec(j) = sum(this%binWindows(:,i)*this%ClFiducial(j,:))
    !        end do
    !        call this%ElementsToMatrix(avec, this%sqrt_fiducial(i)%M)
    !        if (.not. cl_fiducial_includes_noise) &
    !        & this%sqrt_fiducial(i)%M= this%sqrt_fiducial(i)%M + this%NoiseM(i)%M
    !        call Matrix_Root(this%sqrt_fiducial(i)%M, this%nmaps, 0.5_mcp)
    !    end if
    !    do clix =1, this%ncl_hat
    !        allocate(this%ChatM(i,clix)%M(this%nmaps,this%nmaps))
    !        do j=1,this%ncl
    !            avec(j) = sum(this%binWindows(:,i)*this%ClHat(j,:,clix))
    !        end do
    !        call this%ElementsToMatrix(avec, this%ChatM(i,clix)%M)
    !    end do
    !end do
    !deallocate(avec)

    do l=this%bin_min,this%bin_max
        allocate(this%ChatM(l)%M(this%nmaps,this%nmaps))
        call this%ElementsToMatrix(this%ClHat(:,l), this%ChatM(l)%M)
        if (allocated(this%ClNoise)) then
            allocate(this%NoiseM(l)%M(this%nmaps,this%nmaps))
            call this%ElementsToMatrix(this%ClNoise(:,l), this%NoiseM(l)%M)
        end if
        if (allocated(this%ClFiducial)) then
            allocate(this%sqrt_fiducial(l)%M(this%nmaps,this%nmaps))
            if (.not. cl_fiducial_includes_noise) this%ClFiducial(:,l)=this%ClFiducial(:,l)+this%ClNoise(:,l)
            call this%ElementsToMatrix(this%ClFiducial(:,l), this%sqrt_fiducial(l)%M)
            call Matrix_Root(this%sqrt_fiducial(l)%M, this%nmaps, 0.5_mcp)
        end if
    end do

    if (this%like_approx /= like_approx_fullsky_exact) then
        call this%ReadCovmat(Ini)
    end if
    if (Ini%HasKey('lowl_exact')) call MpiStop('lowl_exact has been separated out as not currently used')

    call this%ReadClArr(Ini, 'linear_correction_fiducial',this%FiducialCorrection, hasKey=hasKey)
    if (hasKey) then
        call this%ReadBinWindows(Ini, 'linear_correction_bin_window', this%binCorrectionWindows)
    end if

    this%has_lensing  = this%cl_lmax(CL_Phi,CL_Phi) > 0

    S = Ini%ReadRelativeFileName('calibration_param')
    if (S/='') then
        call this%loadParamNames(S)
        this%calibration_index = this%nuisance_params%nnames
    end if

    call this%TCMBLikelihood%ReadIni(Ini)

    call this%InitMapCls(this%MapCls, this%nmaps_required, this%required_order)

    end subroutine CMBLikes_ReadIni


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
    integer, allocatable :: cov_cl_used(:)
    character(LEN=:), allocatable :: covmat_format


    covmat_cl = Ini%Read_String('covmat_cl', .true.)
    filename = Ini%ReadRelativeFileName('covmat_fiducial', NotFoundFail=.true.)
    call Ini%Read('covmat_scale',covmat_scale)

    call this%UseString_to_cols(covmat_cl, cl_in_index)
    num_in = size(cl_in_index)
    this%ncl_used = count(cl_in_index /=0)
    allocate(this%cl_use_index(this%ncl_used))
    allocate(cov_cl_used(this%ncl_used))
    ix = 0
    do i=1, num_in
        if (cl_in_index(i)/=0) then
            ix = ix+1
            this%cl_use_index(ix) = cl_in_index(i)
            cov_cl_used(ix) = i
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
                    & covmat_scale*Cov( (binx-1)*num_in+cov_cl_used, (biny-1)*num_in+cov_cl_used)
            end do
        end do
        call Matrix_Inverse(this%inv_covariance)
    else
        vecsize_in =  (this%pcl_lmax-this%pcl_lmin+1)
        if (IsMainMPI()) then
            allocate(Cov(vecsize_in*num_in,vecsize_in*num_in))
            call Ini%Read('covmat_scale',covmat_scale)
            covmat_format = Ini%Read_String('covmat_format')
            if (covmat_format == 'symmetric_fortran_binary') then
                call MatrixSym_Read_Binary(filename, Cov)
            else if (covmat_format == 'text') then
                call Matrix_Read(filename, Cov)
            else
                call MpiStop('Unknown covmat_format')
            end if
            allocate(this%inv_covariance(this%nbins_used*this%ncl_used, this%nbins_used*this%ncl_used))
            do i=1, this%ncl_used
                do j=1,this%ncl_used
                    do L1=this%bin_min, this%bin_max
                        do L2 = this%bin_min, this%bin_max
                            !Assume L matrices are ordered the other way in blocks of different L for each CL type
                            this%inv_covariance( (L1-this%bin_min)*this%ncl_used+i, (L2-this%bin_min)*this%ncl_used+j) = &
                                & covmat_scale*Cov((cov_cl_used(i)-1)*vecsize_in + L1 - this%pcl_lmin+1, &
                                & (cov_cl_used(j)-1)*vecsize_in + L2 - this%pcl_lmin+1)
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
        call MPI_BCAST(this%inv_covariance,Size(this%inv_covariance),MPI_real_mcp, 0, MPI_COMM_WORLD, i)
#endif
    end if
    end subroutine ReadCovmat

    subroutine CMBLikes_Transform(this, C, Chat, CfHalf, COffset)
    !Get  C = C_s^{1/2}  U f(this) U^T C_s^{1/2} where C^{-1/2} CHat C^{-1/2} = U this U^T

    !Get  C = C_f^{1/2} C^{-1/2} C^{+1/2} U f(this) U^T C^{+1/2} C^{-1/2} C_f^{1/2} where C^{-1/2} CHat C^{-1/2} = U this U^T

    class(TCMBLikes) :: this
    real(mcp) C(this%nmaps,this%nmaps)
    real(mcp), intent(in), optional :: COffset(this%nmaps,this%nmaps)
    real(mcp), intent(in) :: CHat(this%nmaps,this%nmaps), CfHalf(this%nmaps,this%nmaps)
    real(mcp) :: U(this%nmaps,this%nmaps), Rot(this%nmaps,this%nmaps)
    real(mcp) :: roots(this%nmaps)
    real(mcp) :: diag(this%nmaps)
    integer i

    if (present(COffset)) then
        U = C + Coffset*C
        call Matrix_Diagonalize(U,Diag,this%nmaps)
        Rot= matmul(matmul(transpose(U),CHat+ COffset*C),U)
    else
        U = C
        call Matrix_Diagonalize(U,Diag,this%nmaps)
        Rot= matmul(matmul(transpose(U),CHat),U)
    end if
    roots = sqrt(Diag)

    do i=1, this%nmaps
        Rot(i,:)=Rot(i,:)/roots(i)
        Rot(:,i)=Rot(:,i)/roots(i)
    end do

    Rot = matmul(U,matmul(Rot,transpose(U)))
    call Matrix_Diagonalize(Rot,Diag,this%nmaps)

    Diag = sign(sqrt(2*max(0._mcp,Diag-log(Diag)-1)),Diag-1)
    !want f(this)-1 to save calculating X-X_s

    if (present(COffset)) then
        Rot = MatMul(transpose(U),Rot)
        do i=1, this%nmaps
            Rot(i,:)=Rot(i,:)*roots(i)
        end do
        Rot = MatMul(U,Rot)
        call Matrix_Root(C,this%nmaps, -0.5_mcp)
        Rot = MatMul(C,Rot)
    end if

    U = matmul(CfHalf,Rot)
    C = U
    do i=1, this%nmaps
        C(:,i) = C(:,i)*Diag(i)
    end do
    C = MatMul(C,transpose(U))

    end subroutine CMBLikes_Transform


    subroutine MatrixToElements(this, M, X)
    class(TCMBLikes) :: this
    real(mcp) :: M(this%nmaps,this%nmaps)
    real(mcp) :: X(this%ncl)
    integer ix,i,j

    ix=0
    do i=1, this%nmaps
        do j=1,i
            ix = ix+1
            X(ix) = M(i,j)
        end do
    end do

    end subroutine MatrixToElements

    subroutine MatrixToElementsInt(this, M, X)
    class(TCMBLikes) :: this
    integer, intent(in) :: M(this%nmaps,this%nmaps)
    integer,intent(out) :: X(this%ncl)
    integer ix,i,j

    ix=0
    do i=1, this%nmaps
        do j=1,i
            ix = ix+1
            X(ix) = M(i,j)
        end do
    end do

    end subroutine MatrixToElementsInt


    subroutine ElementsToMatrix(this, X, M)
    class(TCMBLikes) :: this
    real(mcp), intent(out) :: M(this%nmaps,this%nmaps)
    real(mcp), intent(in) :: X(this%ncl)
    integer ix,i,j

    ix=0
    do i=1, this%nmaps
        do j=1,i
            ix = ix+1
            M(i,j) = X(ix)
            M(j,i) = M(i,j)
        end do
    end do

    end subroutine ElementsToMatrix

    function ExactChiSq(this, C,Chat,l)
    class(TCMBLikes) :: this
    real(mcp), intent(in) :: C(this%nmaps,this%nmaps), Chat(this%nmaps,this%nmaps)
    integer, intent(in) :: l
    real(mcp) ExactChiSq
    real(mcp) M(this%nmaps,this%nmaps)

    M = C
    call Matrix_root(M,this%nmaps,-0.5_mcp)
    M = matmul(M,matmul(Chat,M))
    ExactChiSq = (2*l+1)*this%fullsky_exact_fksy*(Matrix_Trace(M) - this%nmaps - MatrixSym_LogDet(M) )

    end function ExactChiSq

    subroutine GetBinnedMapCls(this, MapCls, C, bin)
    class(TCMBLikes) :: this
    class(TMapCrossPowerSpectrum) :: MapCls(:,:)
    real(mcp) Cls(this%ncl), C(this%nmaps,this%nmaps)
    real(mcp) correctionCl(this%ncl)
    integer, intent(in) :: bin

    call this%BinWindows%Bin(MapCls, Cls, bin)
    if (allocated(this%binCorrectionWindows%W)) then
        call this%BinCorrectionWindows%Bin(MapCls, correctionCl, bin)
        Cls = Cls + (correctionCl - this%FiducialCorrection(:,bin))
    end if
    call this%ElementsToMatrix(Cls, C)

    end subroutine GetBinnedMapCls

    subroutine InitMapCls(this, Cls, nmaps, order)
    class(TCMBLikes) :: this
    integer, intent(in) :: nmaps, order(:)
    class(TMapCrossPowerSpectrum), allocatable, target, intent(out) :: Cls(:,:)
    class(TMapCrossPowerSpectrum), pointer :: CL
    integer i,j
    integer f1, f2

    if (nmaps /= size(order)) call MpiStop('CMBLikes InitMapCls: size mismatch')
    allocate(TMapCrossPowerSpectrum::Cls(nmaps, nmaps))
    do i=1, nmaps
        do j=1, i
            CL =>Cls(i,j)
            CL%map_i = order(i)
            CL%map_j = order(j)
            call this%MapPair_to_Theory_i_j(order,i,j,f1,f2)
            CL%theory_i = f1
            CL%theory_j = f2
            allocate(CL%CL(this%pcl_lmin:this%pcl_lmax), source=0._mcp)
        end do
    end do

    end subroutine InitMapCls


    subroutine GetTheoryMapCls(this, Theory, Cls, DataParams)
    class(TCMBLikes), target :: this
    class(TCosmoTheoryPredictions) :: Theory
    class(TMapCrossPowerSpectrum), pointer, intent(out) :: Cls(:,:)
    real(mcp), intent(in) :: DataParams(:)
    integer i,j

    Cls => this%MapCls
    do i=1, this%nmaps_required
        do j=1, i
            !Neat version, buggy in gfortran
            !associate(CL => Cls(i,j))
            !    associate(Th => Theory%Cls(CL%theory_i ,CL%theory_j))
            !        if (allocated(Th%CL)) then
            !            CL%CL(this%pcl_lmin:this%pcl_lmax) = Th%CL(this%pcl_lmin:this%pcl_lmax)
            !        else
            !            CL%CL(this%pcl_lmin:this%pcl_lmax) = 0
            !        end if
            !    end associate
            !end associate
            if (allocated(Theory%Cls(Cls(i,j)%theory_i ,Cls(i,j)%theory_j)%CL)) then
                Cls(i,j)%CL(this%pcl_lmin:this%pcl_lmax) = &
                    Theory%Cls(Cls(i,j)%theory_i ,Cls(i,j)%theory_j)%CL(this%pcl_lmin:this%pcl_lmax)
            else
                Cls(i,j)%CL(this%pcl_lmin:this%pcl_lmax) = 0
            end if
        end do
    end do
    call this%AdaptTheoryForMaps(Cls,DataParams)

    end subroutine GetTheoryMapCls

    subroutine AddForegrounds(this,Cls,DataParams)
    class(TCMBLikes) :: this
    class(TMapCrossPowerSpectrum), target, intent(inout) :: Cls(:,:)
    real(mcp), intent(in) :: DataParams(:)


    end subroutine AddForegrounds

    subroutine AdaptTheoryForMaps(this,Cls,DataParams)
    class(TCMBLikes) :: this
    class(TMapCrossPowerSpectrum), intent(inout) :: Cls(:,:)
    real(mcp), intent(in) :: DataParams(:)
    integer i,j

    call this%AddForegrounds(Cls, DataParams)
    if (this%calibration_index > 0) then
        !Scale T, E, B spectra by the calibration parameter
        do i=1, this%nmaps_required
            do j=1, i
                if (allocated(Cls(i,j)%CL) ) then
                    if (Cls(i,j)%theory_i<=CL_B .and. Cls(i,j)%theory_j<=CL_B) then
                        Cls(i,j)%CL = Cls(i,j)%CL / DataParams(this%calibration_index)**2
                    end if
                end if
            end do
        end do
    end if

    end subroutine AdaptTheoryForMaps


    subroutine CMBLikes_WriteLikelihoodData(this,Theory,DataParams, root)
    Class(TCMBLikes) :: this
    class(TTheoryPredictions) :: Theory
    real(mcp), intent(in) :: DataParams(:)
    class(TMapCrossPowerSpectrum), allocatable :: Cls(:,:)
    character(LEN=*), intent(in) :: root
    integer L, i, j
    Type(TTextFile) F

    if (.not. this%has_foregrounds) return
    call this%InitMapCls(Cls, this%nmaps_required, this%required_order)
    call this%AddForegrounds(Cls, DataParams)

    F%IntegerFormat = '(*(I6))'
    call F%CreateFile(trim(root)//'.'//this%getTag()//'_foregrounds')
    call F%WriteInLine('#    L')
    do i=1, this%nmaps_required
        do j=1,i
            call F%WriteInLine(this%CL_i_j_Name(this%map_names,Cls(i,j)%map_i,Cls(i,j)%map_j), '(a17)')
        end do
    end do
    call F%NewLine()
    do L = this%pcl_lmin, this%pcl_lmax
        call F%WriteInLine(L)
        do i=1, this%nmaps_required
            do j=1,i
                call F%WriteInLine( Cls(i,j)%CL(L) )
            end do
        end do
        call F%Newline()
    end do
    call F%Close()

    end subroutine CMBLikes_WriteLikelihoodData


    function CMBLikes_LogLike(this, CMB, Theory, DataParams)  result (LogLike)
    real(mcp) logLike
    class(TCMBLikes) :: this
    Class(CMBParams) CMB
    Class(TCosmoTheoryPredictions), target :: Theory
    real(mcp) DataParams(:)
    real(mcp) chisq
    real(mcp) C(this%nmaps,this%nmaps)
    real(mcp) vecp(this%ncl)
    real(mcp) bigX((this%bin_max-this%bin_min+1)*this%ncl_used)
    integer  i,j, bin
    logical :: quadratic
    class(TMapCrossPowerSpectrum), pointer :: TheoryCls(:,:)

    chisq =0

    call this%GetTheoryMapCls(Theory, TheoryCls, DataParams)

    do bin = this%bin_min, this%bin_max
        if (this%binned .or. bin_test) then
            if (this%like_approx == like_approx_fullsky_exact) call mpiStop('CMBLikes: exact like cannot be binned!')
            call this%GetBinnedMapCls(TheoryCls, C, bin)
        else
            if (this%nmaps /= this%nmaps_required) call MpiStop('CMBlikes: Unbinned must have required==used')
            do i=1, this%nmaps
                do j=1, i
                    if (allocated(TheoryCls(i,j)%CL)) then
                        C(i,j) =TheoryCLs(i,j)%CL(bin)
                    else
                        C(i,j)=0
                    end if
                    C(j,i) = C(i,j)
                end do
            end do

        end if

        if (allocated(this%NoiseM)) then
            C = C + this%NoiseM(bin)%M
        end if

        if (this%like_approx == like_approx_HL) then
            call this%Transform(C, this%ChatM(bin)%M, this%sqrt_fiducial(bin)%M)
            call this%MatrixToElements(C, vecp)
            quadratic = .true.
        else if (this%like_approx == like_approx_fid_gaussian) then
            call this%MatrixToElements(C- this%ChatM(bin)%M, vecp)
            quadratic = .true.
        else if (this%like_approx == like_approx_fullsky_exact) then
            quadratic = .false.
            chisq = chisq  + this%ExactChisq(C,this%ChatM(bin)%M,bin)
        else
            call MpiStop('Unknown like_approx')
        end if

        if (quadratic) then
            bigX( (bin-this%bin_min)*this%ncl_used + 1:(bin-this%bin_min+1)*this%ncl_used) = vecp(this%cl_use_index)
        end if
    end do

    if (quadratic) chisq = chisq + Matrix_QuadForm(this%inv_covariance,BigX)


    LogLike = chisq/2

    end function CMBLikes_LogLike


    subroutine TBinWindows_bin(this, TheoryCls, Cls, bin)
    class(TBinWindows) :: this
    class(TSkyPowerSpectrum) :: TheoryCls(:,:)
    real(mcp) Cls(:)
    integer bin
    integer win_ix,ix_in(2), ix_out

    cls=0
    do win_ix = 1, size(this%bin_cols_in,2)
        ix_in = this%bin_cols_in(:,win_ix)
        ix_out = this%bin_cols_out(win_ix)
        if (ix_out>0) then
            if (allocated(this%fixCls)) then
                if (allocated(this%fixCls(ix_in(1),ix_in(2))%CL)) then
                    Cls(ix_out) = Cls(ix_out) + &
                        & dot_product(this%W(:,win_ix,bin),this%fixCls(ix_in(1),ix_in(2))%CL(this%lmin:this%lmax))
                    cycle
                end if
            end if
            if (allocated(TheoryCls(ix_in(1),ix_in(2))%CL)) then
                Cls(ix_out) = Cls(ix_out) + &
                    & dot_product(this%W(:,win_ix,bin),TheoryCls(ix_in(1),ix_in(2))%CL(this%lmin:this%lmax))
            end if
        end if
    end do

    end subroutine TBinWindows_bin


    end module CMBLikes
