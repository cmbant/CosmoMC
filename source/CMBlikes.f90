    !Pseudo-Cl (or other C_l esimator) based likelihood approximation for cut sky with polarization
    !Simple harmonic low-l likelihood
    !Obviously this is not a realistic Planck likelihood code
    !AL Mar 2010 - fixed bug using on E, added support for multiple input simulated Chat for bias testing
    !Apr 2011, added fullsky_exact_fksy
    ! - renaned from Planck_like

    module CMBLikes
    use settings
    use CosmologyTypes
    use MatrixUtils
    use Likelihood_Cosmology
    implicit none
    private

    integer :: cl_ix_E = 3, cl_ix_B=4
    logical, parameter :: bin_test = .false.

    Type TSqMatrix
        real(mcp), dimension(:,:), allocatable :: M
    end Type TSqMatrix

    Type TLowlLike
        integer nmodes
        integer tmodes, almmodes
        integer polalmmodes, EBmodes
        integer lexact, lmax
        double precision, dimension(:), pointer :: ModeDataVector
        double precision, dimension(:,:), pointer :: TheoryProj, ReProj, ImProj
        double precision highlScaleT, highlScaleE, highlScaleC, highlScaleB
        double precision, dimension(:,:) , pointer :: NoiseCov, HighlCov
    end Type TLowlLike

    Type, extends(TCMBLikelihood) :: TCMBLikes
        logical highl_cl, lowl_exact !Which bits of likelihood to include

        integer nfields !number of fields
        logical use_field(3)
        integer field_index(3) !mapping from 1,2,3 to index actually used
        integer fields(3) !indices (1=T,2=E,3=B) of fields to use
        character field_order(3)
        integer ncl !calculated from above = nfields*(nfields+1)/2
        integer ncl_used !Number of C_l actually used in covariance matrix (others assumed zero)
        integer cl_use_index(6)
        integer pcl_lmin, pcl_lmax !The range of l to use psuedo-cl-based likelihood
        integer bin_width
        integer vecsize
        integer nbins
        integer like_approx
        real(mcp) fullsky_exact_fksy ! only used for testing with exactly fullsky
        integer ncl_hat !1, or more if using multiple simulations
        real(mcp), dimension(:,:), allocatable :: ClFiducial, ClNoise, ClPointsources, ClOffset
        real(mcp), dimension(:,:,:), allocatable :: ClHat
        !ClOffset is the alpha parameter determining skewness

        integer cl_phi_lmin, cl_phi_lmax !lmax for the lensing reconstruction
        integer lensing_recon_ncl !0 for no lensing recon, 1 for phi-phi spectru, 2 phi-phi and phi-T
        integer phi_like_approx
        real(mcp), dimension(:,:), allocatable :: ClPhiHat, ClPhiNoise, phi_inv_covariance
        !for lensing reconstruction
        !note these are [l(l+1)]^4C_l/2pi

        real(mcp), dimension(:,:), allocatable :: inv_covariance
        real(mcp), dimension(:,:), allocatable :: binWindows, beammodes
        integer beam_MCMC_modes
        real(mcp) point_source_error !fractional error in A
        integer pointsource_MCMC_modes
        integer num_nuisance_parameters
        Type(TSqMatrix) ,dimension(:), allocatable :: sqrt_fiducial, NoiseM, OffsetM
        Type(TSqMatrix) ,dimension(:,:), allocatable :: ChatM
        Type(TLowlLike) :: LowL
    contains
    procedure :: CMBLike => CMBLikes_CMBLike
    procedure :: ReadIni => CMBLikes_ReadIni
    end Type TCMBLikes

    character(LEN=3), parameter :: field_names = 'TEB'
    integer, parameter :: like_approx_diag=1   !new approximation from Hammimeche & Lewis arXiv: 0801.0554
    integer, parameter :: like_approx_fid_gaussian=2 !fixed covariance matrix, (X-Xhat)^TC^{-1}(X-Xhat)
    integer, parameter :: like_approx_fullsky_exact=3 !ignore all correlations, use exact full sky likelihood function
    integer, parameter :: like_approx_gaussian=4 !includes theory dependent determinant

    public TCMBLikes
    contains


    subroutine CMBLikes_ReadLowlFile(this,aname)
    Type(TCMBLikes) :: this
    character(LEN=*), intent(in) :: aname
    integer  filemodes, nmodes, i,j
    double precision, allocatable :: coupling_row(:)
    double precision, dimension(:), allocatable :: TModeData, EModeData, BModeData
    Type(TBinaryFile) :: F

    !Note this doesn't currently support chopping to requested fields - always uses TEB
    !Also uses binary fortran files rather than e.g. FITS or other endian-independent standard

    call F%Open(aname)

    read (F%unit) filemodes, this%LowL%tmodes
    if (filemodes /= (this%Lowl%lmax+1)**2) call MpiStop('lowl likelihood lmax mismatch')
    this%LowL%almmodes = (this%Lowl%lexact+1)**2
    allocate(TModeData(this%LowL%tmodes))
    read(F%unit) TModeData
    allocate(coupling_row(filemodes))
    allocate(this%Lowl%TheoryProj(this%Lowl%almmodes , this%LowL%tmodes))
    do i=1, this%LowL%tmodes
        read(F%unit) coupling_row
        this%Lowl%TheoryProj(1:this%Lowl%almmodes,i) = coupling_row(1:this%Lowl%almmodes)
    end do
    deallocate(coupling_row)

    Read(F%unit) this%LowL%highlScaleT
    nmodes = this%LowL%tmodes

    read(F%unit) this%Lowl%highlScaleE, this%Lowl%highlScaleC, this%Lowl%highlScaleB

    read(F%unit) filemodes, this%LowL%EBmodes
    this%Lowl%polalmmodes  =  (this%Lowl%lexact+1)**2-4

    allocate(EModeData(this%LowL%EBmodes))
    allocate(BModeData(this%LowL%EBmodes))
    read(F%unit) EModeData, BModeData

    allocate(this%Lowl%ReProj(this%Lowl%polalmmodes, this%LowL%EBmodes))
    allocate(this%Lowl%ImProj(this%Lowl%polalmmodes, this%LowL%EBmodes))
    allocate(coupling_row(filemodes))
    do i=1, this%LowL%EBmodes
        read(F%unit) coupling_row
        this%Lowl%ReProj(1:this%Lowl%polalmmodes,i) = coupling_row(1:this%Lowl%polalmmodes)
        read(F%unit) coupling_row
        this%Lowl%ImProj(1:this%Lowl%polalmmodes,i) = coupling_row(1:this%Lowl%polalmmodes)
    end do
    deallocate(coupling_row)
    nmodes= nmodes + this%LowL%EBmodes*2


    allocate(this%LowL%NoiseCov(nmodes,nmodes))
    allocate(this%LowL%HighlCov(nmodes,nmodes))
    do i=1,nmodes
        read(F%unit) this%LowL%NoiseCov(1:i,i)
        read(F%unit) this%LowL%HighlCov(1:i,i)
    end do
    do i=1,nmodes
        do j=i+1, nmodes
            this%LowL%NoiseCov(j,i) =  this%LowL%NoiseCov(i,j)
            this%LowL%HighlCov(j,i) =  this%LowL%HighlCov(i,j)
        end do
    end do

    read(F%unit) i
    if (i/=252353) call MpiStop('Bad low l likelihood data file')

    close(F%unit)

    allocate( this%LowL%ModeDataVector(nmodes))
    this%LowL%ModeDataVector(1:this%LowL%tmodes) = TModeData
    this%LowL%ModeDataVector(this%LowL%tmodes+1:this%LowL%tmodes+this%LowL%EBmodes) = EModeData
    this%LowL%ModeDataVector(this%LowL%tmodes+this%LowL%EBmodes+1:this%LowL%tmodes+2*this%LowL%EBmodes) = BModeData
    deallocate(TModeData,BModeData,EModeData)

    this%LowL%nmodes = nmodes

    end subroutine CMBLikes_ReadLowlFile

    function idx_T(l,m)
    !Conversion from l,m into vectors
    integer, intent(in) :: l,m
    integer idx_T

    idx_T = l*(l+1) + m +1

    end function

    function idx_P(l,m)
    !Conversion from l,m into vectors
    integer, intent(in) :: l,m
    integer idx_P

    idx_P = l*(l+1) + m +1 -4

    end function


    subroutine CMBLikes_lowl_GetFullCovariance(this, Cov, cl, lmin, lsum)
    Type(TCMBLikes) :: this
    double precision Cov(:,:)
    real(mcp) cl(:,:)
    integer, intent (in) :: lmin, lsum
    integer i,j
    double precision :: sum1,sum2,tmp,tmpEE,tmpBB, tmpEB, tmpBE
    integer l, mix, mixT

    Cov=0

    do i=1, this%Lowl%tmodes
        do j=1, i
            tmp = 0
            do l=lmin,lsum
                mix = idx_T(l,-l)
                tmp =tmp + dot_product(this%Lowl%TheoryProj(mix:mix+2*l,i),this%Lowl%TheoryProj(mix:mix+2*l,j))*cl(l,1)
            end do
            Cov(i,j) = tmp
            if (i/=j) Cov(j,i) =  tmp
        end do
    end do

    do i=1, this%Lowl%tmodes
        do j=1, this%Lowl%EBmodes
            !TE
            tmp = 0
            do l=lmin,lsum
                mix = idx_P(l,-l)
                mixT = idx_T(l,-l)
                tmp =tmp + dot_product(this%Lowl%TheoryProj(mixT:mixT+2*l,i),&
                this%Lowl%ReProj(mix:mix+2*l,j))*cl(l,2)
            end do
            Cov(i,this%Lowl%tmodes+j) = tmp
            Cov(this%Lowl%tmodes+j,i) =  tmp

            !TB
            tmp = 0
            do l=lmin,lsum
                mix = idx_P(l,-l)
                mixT = idx_T(l,-l)
                tmp =tmp - dot_product(this%Lowl%TheoryProj(mixT:mixT+2*l,i),&
                this%Lowl%ImProj(mix:mix+2*l,j))*cl(l,2)
            end do
            Cov(i,this%Lowl%tmodes+this%Lowl%EBmodes+j) = tmp
            Cov(this%Lowl%tmodes+this%Lowl%EBmodes+j,i) =  tmp
        end do
    end do


    do i=1, this%Lowl%EBmodes
        do j=1, i
            !EE
            tmpEE = 0
            tmpBB = 0
            do l=lmin,lsum
                mix = idx_P(l,-l)
                sum1=dot_product(this%Lowl%ReProj(mix:mix+2*l,i),this%Lowl%ReProj(mix:mix+2*l,j))
                sum2=dot_product(this%Lowl%ImProj(mix:mix+2*l,i),this%Lowl%ImProj(mix:mix+2*l,j))
                tmpEE =tmpEE + sum1*cl(l,cl_ix_E)  + sum2*cl(l,cl_ix_B)
                tmpBB =tmpBB + sum1*cl(l,cl_ix_B)  + sum2*cl(l,cl_ix_E)
            end do
            Cov(this%Lowl%tmodes+i,this%Lowl%tmodes+j) = tmpEE
            if (i/=j) Cov(this%Lowl%tmodes+j,this%Lowl%tmodes+i) =  tmpEE
            Cov(this%Lowl%tmodes+this%Lowl%EBmodes+i,this%Lowl%tmodes+this%Lowl%EBmodes+j) = tmpBB
            if (i/=j) Cov(this%Lowl%tmodes+this%Lowl%EBmodes+j,this%Lowl%tmodes+this%Lowl%EBmodes+i) =  tmpBB

            !EB/BE
            tmpEB = 0
            tmpBE=0
            do l=lmin,lsum
                mix = idx_P(l,-l)
                sum1=dot_product(this%Lowl%ReProj(mix:mix+2*l,i),this%Lowl%ImProj(mix:mix+2*l,j))
                sum2=dot_product(this%Lowl%ImProj(mix:mix+2*l,i),this%Lowl%ReProj(mix:mix+2*l,j))
                tmpEB =tmpEB - sum1*cl(l,cl_ix_E) + sum2*cl(l,cl_ix_B)
                tmpBE =tmpBE + sum1*cl(l,cl_ix_B) - sum2*cl(l,cl_ix_E)
            end do
            Cov(this%Lowl%tmodes+i,this%Lowl%tmodes+this%Lowl%EBmodes+j) = tmpEB
            Cov(this%Lowl%tmodes+this%Lowl%EBmodes+j,this%Lowl%tmodes+i) =  tmpEB
            Cov(this%Lowl%tmodes+this%Lowl%EBmodes+i,this%Lowl%tmodes+j) = tmpBE
            Cov(this%Lowl%tmodes+j,this%Lowl%tmodes+this%Lowl%EBmodes+i) =  tmpBE
        end do
    end do


    end subroutine CMBLikes_lowl_GetFullCovariance


    function CMBLikes_lowl_CMBLike(this, cl) result (chisq)
    Type(TCMBLikes) :: this
    real(mcp) cl(:,:)
    real(mcp) chisq
    double precision, allocatable :: Cov(:,:)
    integer j


    print *,'getting low l'

    allocate(Cov(this%Lowl%nmodes,this%Lowl%nmodes))

    call CMBLikes_lowl_GetFullCovariance(this, Cov, cl, 2, this%Lowl%lexact)

    do j=1,this%Lowl%nmodes
        Cov(:,j) = Cov(:,j)+this%Lowl%NoiseCov(:,j)
    end do

    !Scale high l
    Cov(1:this%Lowl%tmodes,1:this%Lowl%tmodes) = &
    Cov(1:this%Lowl%tmodes,1:this%Lowl%tmodes)  &
    + cl(this%Lowl%lexact+1,1)/this%Lowl%highlScaleT*this%Lowl%HighlCov(1:this%Lowl%tmodes,1:this%Lowl%tmodes)


    !TE
    Cov(1:this%Lowl%tmodes,this%Lowl%tmodes+1:this%Lowl%tmodes+this%Lowl%EBmodes) = &
    Cov(1:this%Lowl%tmodes,this%Lowl%tmodes+1:this%Lowl%tmodes+this%Lowl%EBmodes) + &
    sqrt(cl(this%Lowl%lexact+1,1)/this%Lowl%highlScaleT*cl(this%Lowl%lexact+1,3)/this%Lowl%highlScaleE) * &
    this%Lowl%HighlCov(1:this%Lowl%tmodes,this%Lowl%tmodes+1:this%Lowl%tmodes+this%Lowl%EBmodes)

    Cov(this%Lowl%tmodes+1:this%Lowl%tmodes+this%Lowl%EBmodes,1:this%Lowl%tmodes) = &
    transpose(Cov(1:this%Lowl%tmodes,this%Lowl%tmodes+1:this%Lowl%tmodes+this%Lowl%EBmodes))


    !EE
    Cov(this%Lowl%tmodes+1:this%Lowl%tmodes+this%Lowl%EBmodes,this%Lowl%tmodes+1:this%Lowl%tmodes+this%Lowl%EBmodes) = &
    Cov(this%Lowl%tmodes+1:this%Lowl%tmodes+this%Lowl%EBmodes,this%Lowl%tmodes+1:this%Lowl%tmodes+this%Lowl%EBmodes) + &
    cl(this%Lowl%lexact+1,3)/this%Lowl%highlScaleE* &
    this%Lowl%HighlCov(this%Lowl%tmodes+1:this%Lowl%tmodes+this%Lowl%EBmodes, &
    this%Lowl%tmodes+1:this%Lowl%tmodes+this%Lowl%EBmodes)


    if (CosmoSettings%num_cls > 3) then
        !BB
        Cov(this%Lowl%tmodes+this%Lowl%EBmodes+1:this%Lowl%tmodes+2*this%Lowl%EBmodes,&
        this%Lowl%tmodes+this%Lowl%EBmodes+1:this%Lowl%tmodes+2*this%Lowl%EBmodes) = &
        Cov(this%Lowl%tmodes+this%Lowl%EBmodes+1:this%Lowl%tmodes+2*this%Lowl%EBmodes, &
        this%Lowl%tmodes+this%Lowl%EBmodes+1:this%Lowl%tmodes+2*this%Lowl%EBmodes) + &
        cl(this%Lowl%lexact+1,cl_ix_B)/this%Lowl%highlScaleB* &
        this%Lowl%HighlCov(this%Lowl%tmodes+this%Lowl%EBmodes+1:this%Lowl%tmodes+2*this%Lowl%EBmodes, &
        this%Lowl%tmodes+this%Lowl%EBmodes+1:this%Lowl%tmodes+2*this%Lowl%EBmodes)


        Cov(this%Lowl%tmodes+this%Lowl%EBmodes+1:this%Lowl%tmodes+2*this%Lowl%EBmodes,&
        this%Lowl%tmodes+1:this%Lowl%tmodes+this%Lowl%EBmodes) = &
        Cov(this%Lowl%tmodes+this%Lowl%EBmodes+1:this%Lowl%tmodes+2*this%Lowl%EBmodes,&
        this%Lowl%tmodes+1:this%Lowl%tmodes+this%Lowl%EBmodes) + &
        cl(this%Lowl%lexact+1,cl_ix_E)/this%Lowl%highlScaleE* &
        this%Lowl%HighlCov(this%Lowl%tmodes+this%Lowl%EBmodes+1:this%Lowl%tmodes+2*this%Lowl%EBmodes,&
        this%Lowl%tmodes+1:this%Lowl%tmodes+this%Lowl%EBmodes)
        !EB
        Cov(this%Lowl%tmodes+1:this%Lowl%tmodes+this%Lowl%EBmodes, &
        this%Lowl%tmodes+this%Lowl%EBmodes+1:this%Lowl%tmodes+2*this%Lowl%EBmodes) = &
        transpose(Cov(this%Lowl%tmodes+this%Lowl%EBmodes+1:this%Lowl%tmodes+2*this%Lowl%EBmodes,&
        this%Lowl%tmodes+1:this%Lowl%tmodes+this%Lowl%EBmodes))
    end if

    chisq = 2*Matrix_GaussianLogLikeDouble(Cov,this%LowL%ModeDataVector)

    if (Feedback > 1) print *,'lowl chisq = ', chisq
    deallocate(Cov)

    end function CMBLikes_lowl_CMBLike

    function TypeIndex(C)
    character(LEN=*), intent(in) :: C
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

    subroutine CMBLikes_ReadClArr(this, aname, order, Cl, lmin, keepnorm)
    Type(TCMBLikes) :: this
    character(LEN=*), intent(in) :: aname, order
    logical, intent(in), optional :: keepnorm
    integer, intent(in) :: lmin
    real(mcp) :: Cl(:,lmin:)
    character(LEN=:), allocatable :: tmp
    integer ix, i, j,i1,l,ll
    integer cols(6)
    logical donorm
    real(mcp) norm,tmp_ar(6)
    Type(TStringList) :: Li
    integer status
    Type(TTextFile) :: F

    if (present(keepnorm)) then
        donorm = .not. keepnorm
    else
        donorm = .true.
    end if
    call Li%SetFromString(order,field_names)
    ix=0
    cols=0
    do i=1,this%nfields
        do j=1,i
            ix = ix +1
            i1 =Li%IndexOf(this%field_order(i)//this%field_order(j))
            if (i1==-1) i1 = Li%IndexOf(this%field_order(j)//this%field_order(i))
            if (i1/=-1) then
                cols(ix) = i1
            end if
        end do
    end do

    call F%Open(aname)
    Cl=0
    do while (F%ReadLine(tmp))
        read(tmp,*, iostat=status) l, tmp_ar(1:Li%Count)
        if (status/=0) call MpiStop('CMBLikes_ReadClArr: error reading line '//trim(aname))
        ll=l
        if (l>=this%pcl_lmin .and. l <=this%pcl_lmax) then
            if (donorm) then
                norm = l*(l+1)/twopi
            else
                norm =1
            end if
            do ix=1,this%ncl
                if (cols(ix)/=0) Cl(ix,l) = tmp_ar(cols(ix))/norm
            end do
        end if
    end do
    if (ll<this%pcl_lmax) then
        write(*,*) 'CMBLikes_ReadClArr: C_l file does not go up to lmax:', this%pcl_lmax
        write (*,*) trim(aname)
        call MpiStop()
    end if
    call F%Close()

    call Li%Clear()

    end subroutine CMBLikes_ReadClArr


    subroutine UseString_to_colIx(this, S, C, totnum)
    Type(TCMBLikes) :: this
    character(LEN=*), intent(in) :: S
    integer, intent(inout), optional :: totnum
    integer :: C(:,:)
    Type(TStringList) :: L
    integer i,i1,i2
    character(LEN=:), pointer :: P

    call L%SetFromString(S, 'TEB')
    C=0

    do i=1, L%Count
        P=> L%Item(i)
        if (len(P)/=2) call mpiStop('Invalid C_l order')
        i1= this%field_index(TypeIndex(L%CharAt(i,1)))
        i2= this%field_index(TypeIndex(L%CharAt(i,2)))
        if (i1/=0 .and. i2/=0) then
            C(i1,i2) = i
            C(i2,i1)= i
        end if
    end do

    if (present(totnum)) totnum = L%Count
    call L%Clear()

    end subroutine UseString_to_colIx


    subroutine CMBLike_ReadModes(this,modes, fname, nmodes)
    Type(TCMBLikes) :: this
    character(LEN=*), intent(in) :: fname
    integer nmodes
    real(mcp) x, modes(this%pcl_lmin:this%pcl_lmax,nmodes)
    integer stat, l, i
    Type(TTextFile) :: F

    call F%Open(fname)
    do
        read(F%unit,*,iostat=stat) i, l, x
        if ( stat /= 0 ) exit
        if ( i <= nmodes .and. l>= this%pcl_lmin .and. l<=this%pcl_lmax ) then
            modes(l,i) = x
        end if
    end do
    call F%Close()

    end subroutine CMBLike_ReadModes

    subroutine CMBLikes_ReadIni(this, Ini)
    class(TCMBLikes) :: this
    class(TSettingIni) :: Ini
    real(mcp), dimension(:,:), allocatable, target :: Cov, fullcov
    integer ix, i
    character(LEN=:), allocatable :: S, S_order
    integer l,j, x,y, clix
    integer l1,l2
    integer, dimension(:,:), allocatable :: indices
    integer lmin_covmat,lmax_covmat, vecsize_in
    integer nmodes,cov_num_cls
    real(mcp) covmat_scale
    double precision :: asum
    real(mcp), allocatable :: avec(:)

    S = Ini%Read_String('fields_use', .true.)
    this%use_field = .false.
    do i=1, len_trim(S)
        if (trim(S(i:i))/='') this%use_field(TypeIndex(S(i:i))) = .true.
    end do

    this%highl_cl = Ini%Read_Logical('highl_cl')

    this%nfields=0
    this%vecsize =0
    this%ncl=0
    this%ncl_used=0
    this%num_nuisance_parameters=0

    if (this%highl_cl) then
        this%like_approx = Ini%read_Int('like_approx')
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

        this%pcl_lmin = Ini%Read_Int('pcl_lmin')
        this%pcl_lmax = Ini%Read_Int('pcl_lmax')
        this%vecsize = (this%pcl_lmax-this%pcl_lmin+1)
        this%ncl_used = 0
        this%bin_width = Ini%Read_Int('bin_width',1)
        this%nbins = this%vecsize/this%bin_width !Make last bin bigger if not exact multiple


        if (this%bin_width/=1 .or. bin_test) then
            allocate(this%binWindows(this%pcl_lmin:this%pcl_lmax,this%nbins))
            this%binWindows=0
            do i=this%pcl_lmin,this%pcl_lmin+this%nbins*this%bin_width -1
                this%binWindows(i,(i-this%pcl_lmin)/this%bin_width+1)=2*i+1
            end do
            do i=this%pcl_lmin+this%nbins*this%bin_width,this%pcl_lmax
                this%binWindows(i,this%nbins)=2*i+1
            end do
            do i=1, this%nbins
                this%binWindows(:,i) = this%binWindows(:,i)/sum(this%binWindows(:,i))
            end do
        end if

        this%ncl_hat = Ini%Read_Int('ncl_hat', 1)

        allocate(this%ClHat(this%ncl,this%pcl_lmin:this%pcl_lmax, this%ncl_hat))
        allocate(this%ClNoise(this%ncl,this%pcl_lmin:this%pcl_lmax))

        S = Ini%ReadFileName('cl_hat_file')
        S_order = Ini%read_String('cl_hat_order', .true.)
        call CMBLikes_ReadClArr(this, S,S_order,this%ClHat(:,:,1),this%pcl_lmin)
        do j=2, this%ncl_hat
            !for simulated with multiple realizations with same covariance and noise
            call CMBLikes_ReadClArr(this, Ini%ReadFileName(numcat('cl_hat_file',j)),&
            S_order,this%ClHat(:,:,j),this%pcl_lmin)
        end do

        if (this%like_approx /= like_approx_fullsky_exact) then
            allocate(this%ClFiducial(this%ncl,this%pcl_lmin:this%pcl_lmax))
            S =Ini% ReadFileName('cl_fiducial_file')
            S_order = Ini%read_String('cl_fiducial_order', .true.)
            call CMBLikes_ReadClArr(this, S,S_order,this%ClFiducial,this%pcl_lmin)
        else
            !Exact like
            this%fullsky_exact_fksy = Ini%Read_Real('fullsky_exact_fksy', 1.)
        end if

        S = Ini%ReadFileName('cl_noise_file')
        S_order = Ini%read_String('cl_noise_order', .true.)
        call CMBLikes_ReadClArr(this, S,S_order,this%ClNoise,this%pcl_lmin)

        this%lensing_recon_ncl = Ini%Read_Int('lensing_recon_ncl', 0)
        if (this%lensing_recon_ncl > 0) call CMBLikes_ReadLensingReconData(this, Ini)

        if (.not. Ini%Read_Logical('cl_hat_includes_noise')) then
            do j=1,this%ncl_hat
                this%ClHat(:,:,j) =  this%ClHat(:,:,j) + this%ClNoise
            end do
        end if

        S = Ini%ReadFileName('cl_offset_file',NotFoundFail=.false.)
        if (S/='') then
            S_order = Ini%read_String('cl_offset_order', .true.)
            allocate(this%ClOffset(this%ncl,this%pcl_lmin:this%pcl_lmax))
            call CMBLikes_ReadClArr(this, S,S_order,this%ClOffset,this%pcl_lmin, .true.)
        end if

        S = Ini%ReadFileName('point_source_cl', NotFoundFail=.false.)
        if (S/='') then
            if (Feedback > 1 .and. IsMainMPI()) print *,'Using point source uncertainty'
            S_order = Ini%read_String('point_source_cl_order',.true.)
            allocate(this%ClPointsources(this%ncl,this%pcl_lmin:this%pcl_lmax))
            call CMBLikes_ReadClArr(this, S,S_order,this%ClPointsources,this%pcl_lmin)
            this%pointsource_MCMC_modes = Ini%Read_Int('pointsource_MCMC_modes');
            this%num_nuisance_parameters = this%num_nuisance_parameters + this%pointsource_MCMC_modes

            this%point_source_error = Ini%Read_Real('point_source_error')
            if (.not. Ini%Read_Logical('cl_noise_includes_pointsources')) then
                this%ClNoise = this%ClNoise + this%ClPointsources
                if (.not. Ini%Read_Logical('cl_hat_includes_pointsources')) then
                    do j=1,this%ncl_hat
                        this%ClHat(:,:,j) =  this%ClHat(:,:,j) + this%ClPointsources
                    end do
                end if
            else if (Ini%Read_Logical('cl_hat_includes_pointsources')) then
                !Already added as part of the noise
                do j=1,this%ncl_hat
                    this%ClHat(:,:,j) =  this%ClHat(:,:,j) - this%ClPointsources
                end do
            end if
        end if

        S = Ini%ReadFileName('beam_modes_file', NotFoundFail=.false.)
        if (S/='') then
            if (Feedback > 1 .and. IsMainMPI()) print *,'Using beam uncertainty modes'
            if (this%ncl/=1) call MpiStop('Planck_like: beam modes currently only for temperature a la WMAP')
            nmodes = Ini%Read_Int('beam_modes_number')
            allocate(this%beammodes(this%pcl_lmin:this%pcl_lmax,nmodes))
            call CMBLike_ReadModes(this, this%beammodes, S, nmodes)
            this%beam_MCMC_modes = Ini%Read_Int('beam_MCMC_modes');
            if (this%beam_MCMC_modes > nmodes) call MpiStop('Planck_like: beam_MCMC_modes > beam_modes number')
            this%num_nuisance_parameters = this%num_nuisance_parameters + this%beam_MCMC_modes
        end if

        if (this%bin_width/=1 .or. bin_test) then
            allocate(this%ChatM(this%nbins, this%ncl_hat))
            allocate(this%NoiseM(this%nbins))
            allocate(this%sqrt_fiducial(this%nbins))
            allocate(avec(this%ncl))
            do i=1, this%nbins
                allocate(this%NoiseM(i)%M(this%nfields,this%nfields))
                do j=1,this%ncl
                    avec(j) = sum(this%binWindows(:,i)*(this%ClNoise(j,:)))
                end do
                call ElementsToMatrix(this, avec, this%NoiseM(i)%M)

                if (allocated(this%ClFiducial)) then
                    allocate(this%sqrt_fiducial(i)%M(this%nfields,this%nfields))
                    do j=1,this%ncl
                        avec(j) = sum(this%binWindows(:,i)*this%ClFiducial(j,:))
                    end do
                    call ElementsToMatrix(this, avec, this%sqrt_fiducial(i)%M)
                    this%sqrt_fiducial(i)%M= this%sqrt_fiducial(i)%M + this%NoiseM(i)%M
                    call Matrix_Root(this%sqrt_fiducial(i)%M, this%nfields, 0.5_mcp)
                end if
                do clix =1, this%ncl_hat
                    allocate(this%ChatM(i,clix)%M(this%nfields,this%nfields))
                    do j=1,this%ncl
                        avec(j) = sum(this%binWindows(:,i)*this%ClHat(j,:,clix))
                    end do
                    call ElementsToMatrix(this, avec, this%ChatM(i,clix)%M)
                end do
            end do
            deallocate(avec)
        else
            allocate(this%sqrt_fiducial(this%pcl_lmin:this%pcl_lmax))
            allocate(this%ChatM(this%pcl_lmin:this%pcl_lmax,this%ncl_hat))
            allocate(this%NoiseM(this%pcl_lmin:this%pcl_lmax))
            if (allocated(this%ClOffset)) allocate(this%OffsetM(this%pcl_lmin:this%pcl_lmax))
            do l=this%pcl_lmin,this%pcl_lmax
                do clix = 1, this%ncl_hat
                    allocate(this%ChatM(l,clix)%M(this%nfields,this%nfields))
                    call ElementsToMatrix(this, this%ClHat(:,l,clix), this%ChatM(l,clix)%M)
                end do
                allocate(this%NoiseM(l)%M(this%nfields,this%nfields))
                call ElementsToMatrix(this, this%ClNoise(:,l), this%NoiseM(l)%M)
                if (allocated(this%ClOffset)) then
                    allocate(this%OffsetM(l)%M(this%nfields,this%nfields))
                    call ElementsToMatrix(this, this%ClOffset(:,l), this%OffsetM(l)%M)
                end if
                if (allocated(this%ClFiducial)) then
                    allocate(this%sqrt_fiducial(l)%M(this%nfields,this%nfields))
                    call ElementsToMatrix(this, this%ClFiducial(:,l)+this%ClNoise(:,l), this%sqrt_fiducial(l)%M)
                    call Matrix_Root(this%sqrt_fiducial(l)%M, this%nfields, 0.5_mcp)
                end if
            end do
        end if

        if (this%like_approx /= like_approx_fullsky_exact) then

        lmax_covmat = Ini%Read_Int('covmat_lmax')
        lmin_covmat = Ini%Read_Int('covmat_lmin')
        if (lmin_covmat > this%pcl_lmin) call MpiStop('lmin_covmat must be  <= pcl_lmin')
        if (lmax_covmat < this%pcl_lmax) call MpiStop('lmax_covmat must be  >= pcl_lmax')
        covmat_scale = Ini%Read_Real('covmat_scale',1.0)
        S = Ini%Read_String('covmat_cl', .true.)

        allocate(indices(this%nfields,this%nfields))
        call UseString_to_colIx(this, S, indices, cov_num_cls)
        call MatrixToElementsInt(this,indices,this%cl_use_index)
        deallocate(indices)

        this%ncl_used = count(this%cl_use_index(1:this%ncl) /=0)

        vecsize_in =  (lmax_covmat-lmin_covmat+1)

        S = Ini%ReadFileName('covmat_fiducial')

        if (IsMainMPI()) then
            allocate(fullcov(this%vecsize*this%ncl_used, this%vecsize*this%ncl_used))
            allocate(Cov(vecsize_in*cov_num_cls,vecsize_in*cov_num_cls))
            call MatrixSym_Read_Binary(S, Cov)
            do i=1, this%ncl_used
                do j=1,this%ncl_used
                    fullcov((i-1)*this%vecsize+1:i*this%vecsize,(j-1)*this%vecsize+1:j*this%vecsize) &
                    = covmat_scale*Cov((i-1)*vecsize_in+(this%pcl_lmin-lmin_covmat+1):(i-1)*vecsize_in + &
                    (this%pcl_lmax-lmin_covmat+1), &
                    (j-1)*vecsize_in+(this%pcl_lmin-lmin_covmat+1):(j-1)*vecsize_in +(this%pcl_lmax-lmin_covmat+1))
                end do
            end do
            deallocate(Cov)
            if (allocated(this%ClPointsources) .and. this%pointsource_MCMC_modes==0) then
                do i=1, this%ncl_used
                    do j=1,this%ncl_used
                        do l1= this%pcl_lmin, this%pcl_lmax
                            do l2= this%pcl_lmin, this%pcl_lmax
                                fullcov((i-1)*this%vecsize+ l1 -this%pcl_lmin+1 ,(j-1)*this%vecsize+ l2 -this%pcl_lmin+1 ) = &
                                fullcov((i-1)*this%vecsize+ l1 -this%pcl_lmin+1,(j-1)*this%vecsize+ l2 -this%pcl_lmin+1 ) +  &
                                this%point_source_error**2*this%ClPointsources(i,l1)*this%ClPointsources(j,l2)
                            end do
                        end do
                    end do
                end do
            end if
            if (allocated(this%beammodes)) then
                do i=this%beam_MCMC_modes+1, nmodes
                    this%BeamModes(this%pcl_lmin:this%pcl_lmax,i)=this%BeamModes(this%pcl_lmin:this%pcl_lmax,i)*&
                    & this%ClFiducial(1,this%pcl_lmin:this%pcl_lmax)
                end do
                do i=this%beam_MCMC_modes+1, nmodes
                    do l1= this%pcl_lmin, this%pcl_lmax
                        do l2= this%pcl_lmin, this%pcl_lmax
                            fullcov(l1 -this%pcl_lmin+1, l2 -this%pcl_lmin+1 ) = &
                            fullcov(l1 -this%pcl_lmin+1, l2 -this%pcl_lmin+1 ) +  this%BeamModes(l1,i)*this%BeamModes(l2,i)
                        end do
                    end do
                end do
            end if
            if (this%bin_width==1 .and. .not. bin_test) then
                allocate(this%inv_covariance(this%vecsize*this%ncl_used, this%vecsize*this%ncl_used))
                this%inv_covariance = fullcov
                deallocate(fullcov)
                !        this%inv_covariance => fullcov
            else
                allocate(this%inv_covariance(this%nbins*this%ncl_used, this%nbins*this%ncl_used))
                do i=1, this%ncl_used
                    do j=1,this%ncl_used
                        do x=1, this%nbins
                            do y=1, this%nbins
                                if (i==j .and. y>x) exit
                                asum=0
                                do ix = 1, this%vecsize
                                    if (this%binWindows(this%pcl_lmin+ix-1,y)/=0) then
                                        asum = asum + sum(this%binWindows(:,x)*fullcov((i-1)*this%vecsize+1: &
                                        & i*this%vecsize,(j-1)*this%vecsize+ix))* this%binWindows(this%pcl_lmin+ix-1,y)
                                    end if
                                end do
                                this%inv_covariance((i-1)*this%nbins+x,(j-1)*this%nbins+y) = asum
                                if (i==j) this%inv_covariance((i-1)*this%nbins+y,(j-1)*this%nbins+x) = asum
                            end do
                        end do
                    end do
                end do
                deallocate(fullcov)
            end if

            call Matrix_inverse(this%inv_covariance)
        else !Not mainMPI
            if (this%bin_width==1 .and. .not. bin_test) then
                allocate(this%inv_covariance(this%vecsize*this%ncl_used, this%vecsize*this%ncl_used))
            else
                allocate(this%inv_covariance(this%nbins*this%ncl_used, this%nbins*this%ncl_used))
            end if
        end if !MainMPI
#ifdef MPI
        call MPI_BCAST(this%inv_covariance,Size(this%inv_covariance),MPI_real, 0, MPI_COMM_WORLD, i)
#endif
        end if

    end if !want high l like

    this%lowl_exact = Ini%Read_Logical('lowl_exact')
    if (this%lowl_exact) then
        if (CosmoSettings%num_cls==3) call MpiStop('CMBLikes current untested for only 3 C_l')
        S = Ini%ReadFileName('lowl_datafile')
        this%Lowl%lexact = Ini%Read_Int('lowl_lexact')
        this%Lowl%lmax = Ini%Read_Int('lowl_lmax')
        call CMBLikes_ReadLowlFile(this,S)
    end if

    end subroutine CMBLikes_ReadIni


    subroutine CMBLikes_ReadLensingReconData(this, Ini)
    Type(TCMBLikes) :: this
    class(TSettingIni) :: Ini
    character(LEN=:), allocatable :: fname

    if (CosmoSettings%num_cls <4) &
    call MpiStop('CMBLikes_ReadLensingReconData: must be run with lensing')

    if (Feedback > 1) print *,'CMBLikes_ReadLensingReconData'

    this%cl_phi_lmin = Ini%Read_Int('cl_phi_lmin', this%pcl_lmin)
    this%cl_phi_lmax = Ini%Read_Int('cl_phi_lmax', this%pcl_lmax)
    this%phi_like_approx = Ini%Read_Int('phi_like_approx',this%like_approx)

    allocate(this%ClPhiHat(this%lensing_recon_ncl,this%cl_phi_lmin:this%cl_phi_lmax))
    allocate(this%ClPhiNoise(this%lensing_recon_ncl,this%cl_phi_lmin:this%cl_phi_lmax))
    call CMBLikes_ReadClPhiArr(this, Ini%ReadFileName('cl_hat_phi_file'),this%ClPhiHat)
    call CMBLikes_ReadClPhiArr(this, Ini%ReadFileName('cl_noise_phi_file'),this%ClPhiNoise)

    if (.not. Ini%Read_Logical('cl_hat_includes_noise')) then
        this%ClPhiHat = this%ClPhiHat + this%ClPhiNoise
    end if

    if (this%phi_like_approx /= like_approx_fullsky_exact) then
        fname = Ini%ReadFileName('covmat_phi_fiducial')
        if (fname /='') then
            allocate(this%phi_inv_covariance(this%cl_phi_lmin:this%cl_phi_lmax,this%cl_phi_lmin:this%cl_phi_lmax))
            call MatrixSym_Read_Binary(fname, this%phi_inv_covariance)
            call Matrix_Inverse(this%phi_inv_covariance)
        end if
    end if

    if (Feedback > 1) print *, 'CMBLikes_ReadLensingReconData done'

    end subroutine CMBLikes_ReadLensingReconData

    subroutine CMBLikes_ReadClPhiArr(this, aname, Cl)
    Type(TCMBLikes) :: this
    character(LEN=*), intent(in) :: aname
    real(mcp) :: Cl(:,this%cl_phi_lmin:), tmp_arr(this%lensing_recon_ncl)
    character(LEN=:), allocatable :: tmp
    integer l, ll, status
    Type(TTextFile) :: F

    call F%Open(aname)
    Cl=0
    ll=0
    do while(F%ReadLine(tmp))
        read(tmp,*, iostat=status) l, tmp_arr
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

    Type(TCMBLikes) :: this
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
    Type(TCMBLikes) :: this
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
    Type(TCMBLikes) :: this
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
    Type(TCMBLikes) :: this
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
    Type(TCMBLikes) :: this
    real(mcp), intent(in) :: C(this%nfields,this%nfields), Chat(this%nfields,this%nfields)
    integer, intent(in) :: l
    real(mcp) ExactChiSq
    real(mcp) M(this%nfields,this%nfields)

    M = C
    call Matrix_root(M,this%nfields,-0.5_mcp)
    M = matmul(M,matmul(Chat,M))
    ExactChiSq = (2*l+1)*this%fullsky_exact_fksy*(Matrix_Trace(M) - this%nfields - MatrixSym_LogDet(M) )

    end function ExactChiSq

    function CMBLikes_LensRecon_Like(this, cl_in) result (chisq)
    Type(TCMBLikes) :: this
    real(mcp), intent(in) :: cl_in(:,:)
    real(mcp) vec(this%cl_phi_lmin:this%cl_phi_lmax)
    integer l, phi_ix
    real(mcp) chisq, Cphi, CPhiHat

    if (Feedback > 1) print *,'CMBLikes_LensRecon_CMBLike'

    if (this%bin_width/=1) call MpiStop('CMBLikes_LensRecon_CMBLike: unsupported option')
    if (this%lensing_recon_ncl /=1) call MpiStop('CMBLikes_LensRecon_CMBLike: unsupported ncl')
    phi_ix = CosmoSettings%num_cls + 1 !put here to avoid compile-time bounds problems
    !not implemented cross-correlation
    chisq = 0
    do l=this%cl_phi_lmin, this%cl_phi_lmax
        Cphi = cl_in(l,phi_ix) + this%ClPhiNoise(1,l)
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

    function CMBLikes_CMBLike(this, cl_in) result (chisq)
    class(TCMBLikes) :: this
    real(mcp), intent(in) :: cl_in(:,:)
    real(mcp)  :: cl(CosmoSettings%lmax,CosmoSettings%num_cls)

    real(mcp) chisq
    real(mcp) C(this%nfields,this%nfields)
    real(mcp) vecp(this%ncl)
    real(mcp) bigX(this%nbins*this%ncl_used)
    integer l,  i, Ti,Ei,Bi, bin,clix
    !  integer mode
    logical :: quadratic

    chisq =0

    do clix = 1, this%ncl_hat ! 1 or sum over chi-squareds of simulations
        cl = cl_in(:,1:CosmoSettings%num_cls) !For the moment have not actually implemented lensing likelihood

        if (this%highl_cl) then
            Ti = this%field_index(1)
            Ei = this%field_index(2)
            Bi = this%field_index(3)

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


            if (Bi/=0 .and. CosmoSettings%num_cls<3) call MpiStop('CMBLikes_CMBLike: Need num_cls =4 to use B modes')

            if (this%bin_width/=1 .or. bin_test) then
                if (this%like_approx == like_approx_fullsky_exact) call mpiStop('Planck_like: exact like cannot be binned!')

                do bin = 1, this%nbins
                    C=0
                    if (Ti/=0) C(Ti,Ti) = sum(this%BinWindows(:,bin)*cl(this%pcl_lmin:this%pcl_lmax,1))
                    if (Ei/=0) then
                        C(Ei,Ei) = sum(this%BinWindows(:,bin)*cl(this%pcl_lmin:this%pcl_lmax,3))
                        if (Ti/=0) then
                            C(Ei,Ti) = sum(this%BinWindows(:,bin)*cl(this%pcl_lmin:this%pcl_lmax,2))
                            C(Ti,Ei) = C(Ei,Ti)
                        end if
                    end if
                    if (Bi/=0) C(Bi,Bi) = sum(this%BinWindows(:,bin)*cl(this%pcl_lmin:this%pcl_lmax,CosmoSettings%num_cls))
                    C = C + this%NoiseM(bin)%M

                    if (this%like_approx == like_approx_diag) then
                        call CMBLikes_Transform(this, C, this%ChatM(bin,clix)%M, this%sqrt_fiducial(bin)%M)
                    else if (this%like_approx == like_approx_fid_gaussian) then
                        C = C - this%ChatM(bin,clix)%M
                    else if (this%like_approx == like_approx_gaussian) then
                        if (this%nfields/=1) call MpiStop('not done Gaussian for polarization')
                        !Det term is generally (n+1)[ log C - log C_f])
                        chisq = chisq + 2*log( C(1,1) / (this%sqrt_fiducial(bin)%M(1,1))**2 )
                        C = (C- this%ChatM(bin,clix)%M) *  (this%sqrt_fiducial(bin)%M)**2 / C
                        quadratic = .true.
                    else
                        call MpiStop('Unknown like_approx')
                    end if

                    quadratic = .true.
                    call MatrixToElements(this, C, vecp)

                    do i=1,this%ncl
                        if (this%cl_use_index(i)/=0) then
                            !            bigX( (this%cl_use_index(i)-1)*this%nbins + bin) = vecp(i)
                            bigX( (i-1)*this%nbins + bin) = vecp(i)
                        end if
                    end do
                end do
            else
                do l = this%pcl_lmin, this%pcl_lmax
                    C=0
                    if (Ti/=0) C(Ti,Ti) = cl(l,1)
                    if (Ei/=0) then
                        C(Ei,Ei) = cl(l,3)
                        if (Ti/=0) then
                            C(Ei,Ti) = cl(l,2)
                            C(Ti,Ei) = C(Ei,Ti)
                        end if
                    end if
                    if (Bi/=0) C(Bi,Bi) =  cl(l,CosmoSettings%num_cls)

                    C =C + this%NoiseM(l)%M

                    if (this%like_approx == like_approx_diag) then
                        if (allocated(this%OffsetM)) then
                            chisq = chisq +2*log( (this%ChatM(l,clix)%M(1,1)+this%OffsetM(l)%M(1,1)*C(1,1)) &
                            & /this%ChatM(l,clix)%M(1,1)/(1+this%OffsetM(l)%M(1,1)))
                            call CMBLikes_Transform(this, C, this%ChatM(l,clix)%M, this%sqrt_fiducial(l)%M, this%OffsetM(l)%M )
                        else
                            call CMBLikes_Transform(this, C, this%ChatM(l,clix)%M, this%sqrt_fiducial(l)%M)
                        end if
                        call MatrixToElements(this, C, vecp)
                        quadratic = .true.
                    else if (this%like_approx == like_approx_fid_gaussian) then
                        call MatrixToElements(this, C- this%ChatM(l,clix)%M, vecp)
                        quadratic = .true.
                    else if (this%like_approx == like_approx_gaussian) then
                        if (this%nfields/=1) call MpiStop('not done Gaussian for polarization')
                        chisq = chisq + 2*log( C(1,1) / (this%sqrt_fiducial(l)%M(1,1))**2 )
                        C = (C- this%ChatM(l,clix)%M) *  (this%sqrt_fiducial(l)%M)**2 / C
                        call MatrixToElements(this, C, vecp)
                        !Det term is generally (n+1)[ log C - log C_f])
                        quadratic = .true.
                    else if (this%like_approx == like_approx_fullsky_exact) then
                        quadratic = .false.
                        chisq = chisq  + ExactChisq(this, C,this%ChatM(l,clix)%M,l)
                    else
                        call MpiStop('Unknown like_approx')
                    end if

                    if (quadratic) then
                        do i=1,this%ncl
                            if (this%cl_use_index(i)/=0) then
                                !      bigX( (this%cl_use_index(i)-1)*this%vecsize + l-this%pcl_lmin+1) = vecp(i)
                                bigX( (i-1)*this%vecsize + l-this%pcl_lmin+1) = vecp(i)
                            end if
                        end do
                    end if
                end do

            end if !binned

            if (quadratic) chisq = chisq + Matrix_QuadForm(this%inv_covariance,BigX)
        end if

        if (this%lowl_exact) then
            chisq = chisq + CMBLikes_lowl_CMBLike(this, cl)
        end if

        if (this%lensing_recon_ncl>0) then
            chisq = chisq + CMBLikes_LensRecon_Like(this, cl_in)
        end if
    end do

    end function CMBLikes_CMBLike


    end module CMBLikes
