    !Pseudo-Cl (or other C_l esimator) based likelihood approximation for cut sky with polarization
    !Simple harmonic low-l likelihood
    !Obviously this is not a realistic Planck likelihood code
    !AL Mar 2010 - fixed bug using on E, added support for multiple input simulated Chat for bias testing
    !Apr 2011, added fullsky_exact_fksy

    module CMBLikes
    use settings
    use cmbtypes
    use IniObjects
    use AMLutils
    use MatrixUtils
    implicit none

    integer :: cl_E = 3, cl_B=4 ! stop compile time errors with num_cls=3
    logical, parameter :: bin_test = .false.

    Type TSqMatrix
        real(mcp), dimension(:,:), pointer :: M
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

    Type TCMBLikes
        logical highl_cl, lowl_exact !Which bits of likelihood to include

        integer nfields !number of fields
        logical use_field(3)
        integer field_index(3) !mapping from 1,2,3 to index actually used
        integer fields(3) !indices (1=T,2=E,3=B) of fields to use
        character field_order(3)
        integer ncl !calculated from above = nfields*(nfields+1)/2
        integer ncl_used !Number of C_l actually used in covariance matrix (others assumed zero)
        integer cl_use_index(6)
        integer cl_lmin, cl_lmax !The range of l to use psuedo-cl-based likelihood
        integer bin_width
        integer vecsize
        integer nbins
        integer like_approx
        real(mcp) fullsky_exact_fksy ! only used for testing with exactly fullsky
        integer ncl_hat !1, or more if using multiple simulations
        real(mcp), dimension(:,:), pointer :: ClFiducial, ClNoise, ClPointsources, ClOffset
        real(mcp), dimension(:,:,:), pointer :: ClHat
        !ClOffset is the alpha parameter determining skewness

        integer cl_phi_lmin, cl_phi_lmax !lmax for the lensing reconstruction
        integer lensing_recon_ncl !0 for no lensing recon, 1 for phi-phi spectru, 2 phi-phi and phi-T
        integer phi_like_approx
        real(mcp), dimension(:,:), pointer :: ClPhiHat, ClPhiNoise, phi_inv_covariance
        !for lensing reconstruction
        !note these are [l(l+1)]^4C_l/2pi

        real(mcp), dimension(:,:), pointer :: inv_covariance
        real(mcp), dimension(:,:), pointer :: binWindows, beammodes
        integer beam_MCMC_modes
        real(mcp) point_source_error !fractional error in A
        integer pointsource_MCMC_modes
        integer num_nuisance_parameters
        Type(TSqMatrix) ,dimension(:), pointer :: sqrt_fiducial, NoiseM, OffsetM
        Type(TSqMatrix) ,dimension(:,:), pointer :: ChatM
        Type(TLowlLike) :: LowL
    end Type TCMBLikes

    character(LEN=3), parameter :: field_names = 'TEB'
    integer, parameter :: like_approx_diag=1   !new approximation from Hammimeche & Lewis arXiv: 0801.0554
    integer, parameter :: like_approx_fid_gaussian=2 !fixed covariance matrix, (X-Xhat)^TC^{-1}(X-Xhat)
    integer, parameter :: like_approx_fullsky_exact=3 !ignore all correlations, use exact full sky likelihood function
    integer, parameter :: like_approx_gaussian=4 !includes theory dependent determinant


    contains


    subroutine CMBLikes_ReadLowlFile(D,aname)
    Type(TCMBLikes) :: D
    character(LEN=*), intent(in) :: aname
    integer  filemodes, nmodes, file_unit,i,j
    double precision, allocatable :: coupling_row(:)
    double precision, dimension(:), allocatable :: TModeData, EModeData, BModeData

    !Note this doesn't currently support chopping to requested fields - always uses TEB
    !Also uses binary fortran files rather than e.g. FITS or other endian-independent standard
    file_unit = new_file_unit()

    call OpenFile(aname,file_unit,'unformatted')

    read (file_unit) filemodes, D%LowL%tmodes
    if (filemodes /= (D%Lowl%lmax+1)**2) call MpiStop('lowl likelihood lmax mismatch')
    D%LowL%almmodes = (D%Lowl%lexact+1)**2
    allocate(TModeData(D%LowL%tmodes))
    read(file_unit) TModeData
    allocate(coupling_row(filemodes))
    allocate(D%Lowl%TheoryProj(D%Lowl%almmodes , D%LowL%tmodes))
    do i=1, D%LowL%tmodes
        read(file_unit) coupling_row
        D%Lowl%TheoryProj(1:D%Lowl%almmodes,i) = coupling_row(1:D%Lowl%almmodes)
    end do
    deallocate(coupling_row)

    Read(file_unit) D%LowL%highlScaleT
    nmodes = D%LowL%tmodes


    read(file_unit) D%Lowl%highlScaleE, D%Lowl%highlScaleC, D%Lowl%highlScaleB

    read(file_unit) filemodes, D%LowL%EBmodes
    D%Lowl%polalmmodes  =  (D%Lowl%lexact+1)**2-4

    allocate(EModeData(D%LowL%EBmodes))
    allocate(BModeData(D%LowL%EBmodes))
    read(file_unit) EModeData, BModeData

    allocate(D%Lowl%ReProj(D%Lowl%polalmmodes, D%LowL%EBmodes))
    allocate(D%Lowl%ImProj(D%Lowl%polalmmodes, D%LowL%EBmodes))
    allocate(coupling_row(filemodes))
    do i=1, D%LowL%EBmodes
        read(file_unit) coupling_row
        D%Lowl%ReProj(1:D%Lowl%polalmmodes,i) = coupling_row(1:D%Lowl%polalmmodes)
        read(file_unit) coupling_row
        D%Lowl%ImProj(1:D%Lowl%polalmmodes,i) = coupling_row(1:D%Lowl%polalmmodes)
    end do
    deallocate(coupling_row)
    nmodes= nmodes + D%LowL%EBmodes*2


    allocate(D%LowL%NoiseCov(nmodes,nmodes))
    allocate(D%LowL%HighlCov(nmodes,nmodes))
    do i=1,nmodes
        read(file_unit) D%LowL%NoiseCov(1:i,i)
        read(file_unit) D%LowL%HighlCov(1:i,i)
    end do
    do i=1,nmodes
        do j=i+1, nmodes
            D%LowL%NoiseCov(j,i) =  D%LowL%NoiseCov(i,j)
            D%LowL%HighlCov(j,i) =  D%LowL%HighlCov(i,j)
        end do
    end do

    read(file_unit) i
    if (i/=252353) call MpiStop('Bad low l likelihood data file')

    call CloseFile(file_unit)

    allocate( D%LowL%ModeDataVector(nmodes))
    D%LowL%ModeDataVector(1:D%LowL%tmodes) = TModeData
    D%LowL%ModeDataVector(D%LowL%tmodes+1:D%LowL%tmodes+D%LowL%EBmodes) = EModeData
    D%LowL%ModeDataVector(D%LowL%tmodes+D%LowL%EBmodes+1:D%LowL%tmodes+2*D%LowL%EBmodes) = BModeData
    deallocate(TModeData,BModeData,EModeData)

    D%LowL%nmodes = nmodes

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


    subroutine CMBLikes_lowl_GetFullCovariance(D, Cov, cl, lmin, lsum)
    Type(TCMBLikes) :: D
    double precision Cov(:,:)
    real(mcp) cl(lmax,num_cls)
    integer, intent (in) :: lmin, lsum
    integer i,j
    double precision :: sum1,sum2,tmp,tmpEE,tmpBB, tmpEB, tmpBE
    integer l, mix, mixT

    Cov=0

    do i=1, D%Lowl%tmodes
        do j=1, i
            tmp = 0
            do l=lmin,lsum
                mix = idx_T(l,-l)
                tmp =tmp + dot_product(D%Lowl%TheoryProj(mix:mix+2*l,i),D%Lowl%TheoryProj(mix:mix+2*l,j))*cl(l,1)
            end do
            Cov(i,j) = tmp
            if (i/=j) Cov(j,i) =  tmp
        end do
    end do

    do i=1, D%Lowl%tmodes
        do j=1, D%Lowl%EBmodes

        !TE
        tmp = 0
        do l=lmin,lsum
            mix = idx_P(l,-l)
            mixT = idx_T(l,-l)
            tmp =tmp + dot_product(D%Lowl%TheoryProj(mixT:mixT+2*l,i),&
            D%Lowl%ReProj(mix:mix+2*l,j))*cl(l,2)
        end do
        Cov(i,D%Lowl%tmodes+j) = tmp
        Cov(D%Lowl%tmodes+j,i) =  tmp

        !TB
        tmp = 0
        do l=lmin,lsum
            mix = idx_P(l,-l)
            mixT = idx_T(l,-l)
            tmp =tmp - dot_product(D%Lowl%TheoryProj(mixT:mixT+2*l,i),&
            D%Lowl%ImProj(mix:mix+2*l,j))*cl(l,2)

        end do
        Cov(i,D%Lowl%tmodes+D%Lowl%EBmodes+j) = tmp
        Cov(D%Lowl%tmodes+D%Lowl%EBmodes+j,i) =  tmp


        end do
    end do


    do i=1, D%Lowl%EBmodes
        do j=1, i

        !EE
        tmpEE = 0
        tmpBB = 0
        do l=lmin,lsum
            mix = idx_P(l,-l)
            sum1=dot_product(D%Lowl%ReProj(mix:mix+2*l,i),D%Lowl%ReProj(mix:mix+2*l,j))
            sum2=dot_product(D%Lowl%ImProj(mix:mix+2*l,i),D%Lowl%ImProj(mix:mix+2*l,j))
            tmpEE =tmpEE + sum1*cl(l,cl_E)  + sum2*cl(l,cl_B)
            tmpBB =tmpBB + sum1*cl(l,cl_B)  + sum2*cl(l,cl_E)
        end do
        Cov(D%Lowl%tmodes+i,D%Lowl%tmodes+j) = tmpEE
        if (i/=j) Cov(D%Lowl%tmodes+j,D%Lowl%tmodes+i) =  tmpEE
        Cov(D%Lowl%tmodes+D%Lowl%EBmodes+i,D%Lowl%tmodes+D%Lowl%EBmodes+j) = tmpBB
        if (i/=j) Cov(D%Lowl%tmodes+D%Lowl%EBmodes+j,D%Lowl%tmodes+D%Lowl%EBmodes+i) =  tmpBB

        !EB/BE
        tmpEB = 0
        tmpBE=0
        do l=lmin,lsum
            mix = idx_P(l,-l)
            sum1=dot_product(D%Lowl%ReProj(mix:mix+2*l,i),D%Lowl%ImProj(mix:mix+2*l,j))
            sum2=dot_product(D%Lowl%ImProj(mix:mix+2*l,i),D%Lowl%ReProj(mix:mix+2*l,j))
            tmpEB =tmpEB - sum1*cl(l,cl_E) + sum2*cl(l,cl_B)
            tmpBE =tmpBE + sum1*cl(l,cl_B) - sum2*cl(l,cl_E)
        end do
        Cov(D%Lowl%tmodes+i,D%Lowl%tmodes+D%Lowl%EBmodes+j) = tmpEB
        Cov(D%Lowl%tmodes+D%Lowl%EBmodes+j,D%Lowl%tmodes+i) =  tmpEB
        Cov(D%Lowl%tmodes+D%Lowl%EBmodes+i,D%Lowl%tmodes+j) = tmpBE
        Cov(D%Lowl%tmodes+j,D%Lowl%tmodes+D%Lowl%EBmodes+i) =  tmpBE

        end do
    end do


    end subroutine CMBLikes_lowl_GetFullCovariance


    function CMBLikes_lowl_CMBLike(D, cl) result (chisq)
    Type(TCMBLikes) :: D
    real(mcp) cl(lmax,num_cls)
    real(mcp) chisq
    double precision, allocatable :: Cov(:,:)
    integer j


    print *,'getting low l'

    allocate(Cov(D%Lowl%nmodes,D%Lowl%nmodes))

    call CMBLikes_lowl_GetFullCovariance(D, Cov, cl, 2, D%Lowl%lexact)

    do j=1,D%Lowl%nmodes
        Cov(:,j) = Cov(:,j)+D%Lowl%NoiseCov(:,j)
    end do

    !Scale high l
    Cov(1:D%Lowl%tmodes,1:D%Lowl%tmodes) = &
    Cov(1:D%Lowl%tmodes,1:D%Lowl%tmodes)  &
    + cl(D%Lowl%lexact+1,1)/D%Lowl%highlScaleT*D%Lowl%HighlCov(1:D%Lowl%tmodes,1:D%Lowl%tmodes)


    !TE
    Cov(1:D%Lowl%tmodes,D%Lowl%tmodes+1:D%Lowl%tmodes+D%Lowl%EBmodes) = &
    Cov(1:D%Lowl%tmodes,D%Lowl%tmodes+1:D%Lowl%tmodes+D%Lowl%EBmodes) + &
    sqrt(cl(D%Lowl%lexact+1,1)/D%Lowl%highlScaleT*cl(D%Lowl%lexact+1,3)/D%Lowl%highlScaleE) * &
    D%Lowl%HighlCov(1:D%Lowl%tmodes,D%Lowl%tmodes+1:D%Lowl%tmodes+D%Lowl%EBmodes)

    Cov(D%Lowl%tmodes+1:D%Lowl%tmodes+D%Lowl%EBmodes,1:D%Lowl%tmodes) = &
    transpose(Cov(1:D%Lowl%tmodes,D%Lowl%tmodes+1:D%Lowl%tmodes+D%Lowl%EBmodes))


    !EE
    Cov(D%Lowl%tmodes+1:D%Lowl%tmodes+D%Lowl%EBmodes,D%Lowl%tmodes+1:D%Lowl%tmodes+D%Lowl%EBmodes) = &
    Cov(D%Lowl%tmodes+1:D%Lowl%tmodes+D%Lowl%EBmodes,D%Lowl%tmodes+1:D%Lowl%tmodes+D%Lowl%EBmodes) + &
    cl(D%Lowl%lexact+1,3)/D%Lowl%highlScaleE* &
    D%Lowl%HighlCov(D%Lowl%tmodes+1:D%Lowl%tmodes+D%Lowl%EBmodes, &
    D%Lowl%tmodes+1:D%Lowl%tmodes+D%Lowl%EBmodes)


    if (num_cls > 3) then
        !BB
        Cov(D%Lowl%tmodes+D%Lowl%EBmodes+1:D%Lowl%tmodes+2*D%Lowl%EBmodes,&
        D%Lowl%tmodes+D%Lowl%EBmodes+1:D%Lowl%tmodes+2*D%Lowl%EBmodes) = &
        Cov(D%Lowl%tmodes+D%Lowl%EBmodes+1:D%Lowl%tmodes+2*D%Lowl%EBmodes, &
        D%Lowl%tmodes+D%Lowl%EBmodes+1:D%Lowl%tmodes+2*D%Lowl%EBmodes) + &
        cl(D%Lowl%lexact+1,cl_B)/D%Lowl%highlScaleB* &
        D%Lowl%HighlCov(D%Lowl%tmodes+D%Lowl%EBmodes+1:D%Lowl%tmodes+2*D%Lowl%EBmodes, &
        D%Lowl%tmodes+D%Lowl%EBmodes+1:D%Lowl%tmodes+2*D%Lowl%EBmodes)


        Cov(D%Lowl%tmodes+D%Lowl%EBmodes+1:D%Lowl%tmodes+2*D%Lowl%EBmodes,&
        D%Lowl%tmodes+1:D%Lowl%tmodes+D%Lowl%EBmodes) = &
        Cov(D%Lowl%tmodes+D%Lowl%EBmodes+1:D%Lowl%tmodes+2*D%Lowl%EBmodes,&
        D%Lowl%tmodes+1:D%Lowl%tmodes+D%Lowl%EBmodes) + &
        cl(D%Lowl%lexact+1,cl_E)/D%Lowl%highlScaleE* &
        D%Lowl%HighlCov(D%Lowl%tmodes+D%Lowl%EBmodes+1:D%Lowl%tmodes+2*D%Lowl%EBmodes,&
        D%Lowl%tmodes+1:D%Lowl%tmodes+D%Lowl%EBmodes)
        !EB
        Cov(D%Lowl%tmodes+1:D%Lowl%tmodes+D%Lowl%EBmodes, &
        D%Lowl%tmodes+D%Lowl%EBmodes+1:D%Lowl%tmodes+2*D%Lowl%EBmodes) = &
        transpose(Cov(D%Lowl%tmodes+D%Lowl%EBmodes+1:D%Lowl%tmodes+2*D%Lowl%EBmodes,&
        D%Lowl%tmodes+1:D%Lowl%tmodes+D%Lowl%EBmodes))

    end if

    chisq = 2*Matrix_GaussianLogLikeDouble(Cov,D%LowL%ModeDataVector)

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

    subroutine CMBLikes_ReadClArr(D, aname, order, Cl, lmin, keepnorm)
    Type(TCMBLikes) :: D
    character(LEN=*), intent(in) :: aname, order
    logical, intent(in), optional :: keepnorm
    integer, intent(in) :: lmin
    real(mcp) :: Cl(:,lmin:)
    character(LEN=1024) :: tmp
    integer ix, i, j,i1,l,ll
    integer cols(6)
    logical donorm
    real(mcp) norm,tmp_ar(6)
    Type (TStringList) :: Li
    integer file_unit

    if (present(keepnorm)) then
        donorm = .not. keepnorm
    else
        donorm = .true.
    end if
    call TStringList_Init(Li)
    call TStringList_SetFromString(Li,order,field_names)
    ix=0
    cols=0
    do i=1,D%nfields
        do j=1,i
            ix = ix +1
            i1 = TStringList_IndexOf(Li, D%field_order(i)//D%field_order(j))
            if (i1==-1) i1 = TStringList_IndexOf(Li, D%field_order(j)//D%field_order(i))
            if (i1/=-1) then
                cols(ix) = i1
            end if
        end do
    end do

    file_unit = new_file_unit()
    call OpenTxtFile(aname, file_unit)
    Cl=0
    do
        read(file_unit,'(a)',end=1) tmp
        read(tmp,*, end=1) l, tmp_ar(1:Li%Count)
        ll=l
        if (l>=D%cl_lmin .and. l <=D%cl_lmax) then
            if (donorm) then
                norm = l*(l+1)/twopi
            else
                norm =1
            end if
            do ix=1,D%ncl
                if (cols(ix)/=0) Cl(ix,l) = tmp_ar(cols(ix))/norm
            end do
        end if
    end do
    if (ll<D%cl_lmax) then
        write(*,*) 'CMBLikes_ReadClArr: C_l file does not go up to lmax:', D%cl_lmax
        write (*,*) trim(aname)
        call MpiStop()
    end if
1   call CloseFile(file_unit)

    call TStringList_Clear(Li)

    end subroutine CMBLikes_ReadClArr


    subroutine UseString_to_colIx(D, S, C, totnum)
    Type(TCMBLikes) :: D
    character(LEN=*), intent(in) :: S
    integer, intent(inout), optional :: totnum
    integer :: C(:,:)
    Type(TStringList) :: L
    integer i,i1,i2

    call TStringList_Init(L)
    call TStringList_SetFromString(L, S, 'TEB')
    C=0

    do i=1, L%Count
        if (size(L%Items(i)%P)/=2) call mpiStop('Invalid C_l order')
        i1= D%field_index(TypeIndex(L%Items(i)%P(1)))
        i2= D%field_index(TypeIndex(L%Items(i)%P(2)))
        if (i1/=0 .and. i2/=0) then
            C(i1,i2) = i
            C(i2,i1)= i
        end if
    end do

    if (present(totnum)) totnum = L%Count
    call TStringList_Clear(L)

    end subroutine UseString_to_colIx


    subroutine CMBLikes_ReadDataFile(D, aname)
    Type(TCMBLikes) :: D
    character(LEN=*) :: aname
    Type(TIniFile) :: Ini
    logical bad

    call Ini%Open(aname, bad, .false.)
    if (bad) then
        call MpiStop('Error opening dataset file '//trim(aname))
    end if

    call CMBLikes_ReadData(D, Ini, ExtractFilePath(aname))

    call Ini%Close()

    end  subroutine CMBLikes_ReadDataFile

    subroutine CMBLike_ReadModes(D,modes, fname, nmodes)
    Type(TCMBLikes) :: D
    character(LEN=*), intent(in) :: fname
    integer nmodes
    real(mcp) x, modes(D%cl_lmin:D%cl_lmax,nmodes)
    integer stat, file_unit, l, i

    file_unit = new_file_unit()
    call OpenTxtFile(fname,file_unit)
    do
        read(file_unit,*,iostat=stat) i, l, x
        if ( stat /= 0 ) exit
        if ( i <= nmodes .and. l>= D%cl_lmin .and. l<=D%cl_lmax ) then
            modes(l,i) = x
        end if
    end do

    call CloseFile(file_unit)

    end subroutine CMBLike_ReadModes

    subroutine CMBLikes_ReadData(D, Ini,dataset_dir)
    Type(TCMBLikes) :: D
    class(TIniFile) :: Ini
    real(mcp), dimension(:,:), allocatable, target :: Cov, fullcov
    character(LEN=*), intent(in) :: dataset_dir
    integer ix, i
    character(LEN=Ini_max_string_len) :: S, S_order
    integer l,j, x,y, clix
    integer l1,l2
    integer, dimension(:,:), allocatable :: indices
    integer lmin_covmat,lmax_covmat, vecsize_in
    integer nmodes,cov_num_cls
    real(mcp) covmat_scale
    double precision :: asum
    real(mcp), allocatable :: avec(:)
    !  character(LEN=Ini_max_string_len) cache_name
    Ini_fail_on_not_found = .true.

    S = Ini%Read_String('fields_use')
    D%use_field = .false.
    do i=1, len_trim(S)
        if (trim(S(i:i))/='') D%use_field(TypeIndex(S(i:i))) = .true.
    end do

    D%highl_cl = Ini%Read_Logical('highl_cl')

    D%nfields=0
    D%vecsize =0
    D%ncl=0
    D%ncl_used=0
    D%num_nuisance_parameters=0

    if (D%highl_cl) then

    D%like_approx = Ini%read_Int('like_approx')
    D%nfields = count(D%use_field)
    ix=0
    D%field_index=0
    do i=1,3
        if (D%use_field(i)) then
            ix=ix+1
            D%field_index(i)=ix
            D%fields(ix) =i
            D%field_order(ix) = field_names(i:i)
        end if
    end do
    D%ncl = (D%nfields*(D%nfields+1))/2

    D%cl_lmin = Ini%Read_Int('cl_lmin')
    D%cl_lmax = Ini%Read_Int('cl_lmax')
    D%vecsize = (D%cl_lmax-D%cl_lmin+1)
    D%ncl_used = 0
    D%bin_width = Ini%Read_Int('bin_width',1)
    D%nbins = D%vecsize/D%bin_width !Make last bin bigger if not exact multiple


    if (D%bin_width/=1 .or. bin_test) then
        allocate(D%binWindows(D%cl_lmin:D%cl_lmax,D%nbins))
        D%binWindows=0
        do i=D%cl_lmin,D%cl_lmin+D%nbins*D%bin_width -1
            D%binWindows(i,(i-D%cl_lmin)/D%bin_width+1)=2*i+1
        end do
        do i=D%cl_lmin+D%nbins*D%bin_width,D%cl_lmax
            D%binWindows(i,D%nbins)=2*i+1
        end do
        do i=1, D%nbins
            D%binWindows(:,i) = D%binWindows(:,i)/sum(D%binWindows(:,i))
        end do
    end if

    D%ncl_hat = Ini%Read_Int('ncl_hat', 1)

    allocate(D%ClHat(D%ncl,D%cl_lmin:D%cl_lmax, D%ncl_hat))
    allocate(D%ClNoise(D%ncl,D%cl_lmin:D%cl_lmax))

    S = ReadIniFileName(Ini,'cl_hat_file',dataset_dir)
    S_order = Ini%read_String('cl_hat_order')
    call CMBLikes_ReadClArr(D, S,S_order,D%ClHat(:,:,1),D%cl_lmin)
    do j=2, D%ncl_hat
        !for simulated with multiple realizations with same covariance and noise
        call CMBLikes_ReadClArr(D, ReadIniFileName(Ini,numcat('cl_hat_file',j),dataset_dir),&
        S_order,D%ClHat(:,:,j),D%cl_lmin)
    end do

    if (D%like_approx /= like_approx_fullsky_exact) then
        allocate(D%ClFiducial(D%ncl,D%cl_lmin:D%cl_lmax))
        S = ReadIniFileName(Ini,'cl_fiducial_file',dataset_dir)
        S_order = Ini%read_String('cl_fiducial_order')
        call CMBLikes_ReadClArr(D, S,S_order,D%ClFiducial,D%cl_lmin)
    else
        !Exact like
        D%fullsky_exact_fksy = Ini%Read_Real('fullsky_exact_fksy', 1.)
        nullify(D%ClFiducial)
    end if

    S = ReadIniFileName(Ini,'cl_noise_file',dataset_dir)
    S_order = Ini%read_String('cl_noise_order')
    call CMBLikes_ReadClArr(D, S,S_order,D%ClNoise,D%cl_lmin)

    D%lensing_recon_ncl = Ini%Read_Int('lensing_recon_ncl', 0)
    if (D%lensing_recon_ncl > 0) call CMBLikes_ReadLensingReconData(D, Ini,dataset_dir)

    if (.not. Ini%Read_Logical('cl_hat_includes_noise')) then
        do j=1,D%ncl_hat
            D%ClHat(:,:,j) =  D%ClHat(:,:,j) + D%ClNoise
        end do
    end if

    S = ReadIniFileName(Ini,'cl_offset_file', dataset_dir,.false.)
    if (S/='') then
        S_order = Ini%read_String('cl_offset_order')
        allocate(D%ClOffset(D%ncl,D%cl_lmin:D%cl_lmax))
        call CMBLikes_ReadClArr(D, S,S_order,D%ClOffset,D%cl_lmin, .true.)
    else
        nullify(D%ClOffset)
    end if

    S = ReadIniFileName(Ini,'point_source_cl', dataset_dir,.false.)
    if (S/='') then
        if (Feedback > 1 .and. IsMainMPI()) print *,'Using point source uncertainty'
        S_order = Ini%read_String('point_source_cl_order')
        allocate(D%ClPointsources(D%ncl,D%cl_lmin:D%cl_lmax))
        call CMBLikes_ReadClArr(D, S,S_order,D%ClPointsources,D%cl_lmin)
        D%pointsource_MCMC_modes = Ini%Read_Int('pointsource_MCMC_modes');
        D%num_nuisance_parameters = D%num_nuisance_parameters + D%pointsource_MCMC_modes

        D%point_source_error = Ini%Read_Real('point_source_error')
        if (.not. Ini%Read_Logical('cl_noise_includes_pointsources')) then
            D%ClNoise = D%ClNoise + D%ClPointsources
            if (.not. Ini%Read_Logical('cl_hat_includes_pointsources')) then
                do j=1,D%ncl_hat
                    D%ClHat(:,:,j) =  D%ClHat(:,:,j) + D%ClPointsources
                end do
            end if
        else if (Ini%Read_Logical('cl_hat_includes_pointsources')) then
            !Already added as part of the noise
            do j=1,D%ncl_hat
                D%ClHat(:,:,j) =  D%ClHat(:,:,j) - D%ClPointsources
            end do
        end if
    else
        nullify(D%ClPointsources)
    end if

    S = ReadIniFileName(Ini,'beam_modes_file', dataset_dir,.false.)
    if (S/='') then
        if (Feedback > 1 .and. IsMainMPI()) print *,'Using beam uncertainty modes'
        if (D%ncl/=1) call MpiStop('Planck_like: beam modes currently only for temperature a la WMAP')
        nmodes = Ini%Read_Int('beam_modes_number')
        allocate(D%beammodes(D%cl_lmin:D%cl_lmax,nmodes))
        call CMBLike_ReadModes(D, D%beammodes, S, nmodes)
        D%beam_MCMC_modes = Ini%Read_Int('beam_MCMC_modes');
        if (D%beam_MCMC_modes > nmodes) call MpiStop('Planck_like: beam_MCMC_modes > beam_modes number')
        D%num_nuisance_parameters = D%num_nuisance_parameters + D%beam_MCMC_modes
    else
        nullify(D%beammodes)
    end if

    nullify(D%OffsetM)
    if (D%bin_width/=1 .or. bin_test) then
        allocate(D%ChatM(D%nbins, D%ncl_hat))
        allocate(D%NoiseM(D%nbins))
        allocate(D%sqrt_fiducial(D%nbins))
        allocate(avec(D%ncl))
        do i=1, D%nbins
            allocate(D%NoiseM(i)%M(D%nfields,D%nfields))
            do j=1,D%ncl
                avec(j) = sum(D%binWindows(:,i)*(D%ClNoise(j,:)))
            end do
            call ElementsToMatrix(D, avec, D%NoiseM(i)%M)

            if (associated(D%ClFiducial)) then
                allocate(D%sqrt_fiducial(i)%M(D%nfields,D%nfields))
                do j=1,D%ncl
                    avec(j) = sum(D%binWindows(:,i)*D%ClFiducial(j,:))
                end do
                call ElementsToMatrix(D, avec, D%sqrt_fiducial(i)%M)
                D%sqrt_fiducial(i)%M= D%sqrt_fiducial(i)%M + D%NoiseM(i)%M
                call Matrix_Root(D%sqrt_fiducial(i)%M, D%nfields, 0.5_mcp)
            end if
            do clix =1, D%ncl_hat
                allocate(D%ChatM(i,clix)%M(D%nfields,D%nfields))
                do j=1,D%ncl
                    avec(j) = sum(D%binWindows(:,i)*D%ClHat(j,:,clix))
                end do
                call ElementsToMatrix(D, avec, D%ChatM(i,clix)%M)
            end do
        end do
        deallocate(avec)
    else
        nullify(D%NoiseM)
        allocate(D%sqrt_fiducial(D%cl_lmin:D%cl_lmax))
        allocate(D%ChatM(D%cl_lmin:D%cl_lmax,D%ncl_hat))
        allocate(D%NoiseM(D%cl_lmin:D%cl_lmax))
        if (associated(D%ClOffset)) allocate(D%OffsetM(D%cl_lmin:D%cl_lmax))
        do l=D%cl_lmin,D%cl_lmax
            do clix = 1, D%ncl_hat
                allocate(D%ChatM(l,clix)%M(D%nfields,D%nfields))
                call ElementsToMatrix(D, D%ClHat(:,l,clix), D%ChatM(l,clix)%M)
            end do
            allocate(D%NoiseM(l)%M(D%nfields,D%nfields))
            call ElementsToMatrix(D, D%ClNoise(:,l), D%NoiseM(l)%M)
            if (associated(D%ClOffset)) then
                allocate(D%OffsetM(l)%M(D%nfields,D%nfields))
                call ElementsToMatrix(D, D%ClOffset(:,l), D%OffsetM(l)%M)
            end if
            if (associated(D%ClFiducial)) then
                allocate(D%sqrt_fiducial(l)%M(D%nfields,D%nfields))
                call ElementsToMatrix(D, D%ClFiducial(:,l)+D%ClNoise(:,l), D%sqrt_fiducial(l)%M)
                call Matrix_Root(D%sqrt_fiducial(l)%M, D%nfields, 0.5_mcp)
            end if
        end do

    end if

    if (D%like_approx /= like_approx_fullsky_exact) then


    lmax_covmat = Ini%Read_Int('covmat_lmax')
    lmin_covmat = Ini%Read_Int('covmat_lmin')
    if (lmin_covmat > D%cl_lmin) call MpiStop('lmin_covmat must be  <= cl_lmin')
    if (lmax_covmat < D%cl_lmax) call MpiStop('lmax_covmat must be  >= cl_lmax')
    covmat_scale = Ini%Read_Real('covmat_scale',1.0)
    S = Ini%Read_String('covmat_cl')

    allocate(indices(D%nfields,D%nfields))
    call UseString_to_colIx(D, S, indices, cov_num_cls)
    call MatrixToElementsInt(D,indices,D%cl_use_index)
    deallocate(indices)

    D%ncl_used = count(D%cl_use_index(1:D%ncl) /=0)

    vecsize_in =  (lmax_covmat-lmin_covmat+1)


    S = ReadIniFileName(Ini,'covmat_fiducial',dataset_dir)

    if (IsMainMPI()) then
        allocate(fullcov(D%vecsize*D%ncl_used, D%vecsize*D%ncl_used))
        allocate(Cov(vecsize_in*cov_num_cls,vecsize_in*cov_num_cls))
        call MatrixSym_Read_Binary(S, Cov)
        do i=1, D%ncl_used
            do j=1,D%ncl_used
                fullcov((i-1)*D%vecsize+1:i*D%vecsize,(j-1)*D%vecsize+1:j*D%vecsize) &
                = covmat_scale*Cov((i-1)*vecsize_in+(D%cl_lmin-lmin_covmat+1):(i-1)*vecsize_in +(D%cl_lmax-lmin_covmat+1), &
                (j-1)*vecsize_in+(D%cl_lmin-lmin_covmat+1):(j-1)*vecsize_in +(D%cl_lmax-lmin_covmat+1))
            end do
        end do
        deallocate(Cov)
        if (associated(D%ClPointsources) .and. D%pointsource_MCMC_modes==0) then
            do i=1, D%ncl_used
                do j=1,D%ncl_used
                    do l1= D%cl_lmin, D%cl_lmax
                        do l2= D%cl_lmin, D%cl_lmax
                            fullcov((i-1)*D%vecsize+ l1 -D%cl_lmin+1 ,(j-1)*D%vecsize+ l2 -D%cl_lmin+1 ) = &
                            fullcov((i-1)*D%vecsize+ l1 -D%cl_lmin+1,(j-1)*D%vecsize+ l2 -D%cl_lmin+1 ) +  &
                            D%point_source_error**2*D%ClPointsources(i,l1)*D%ClPointsources(j,l2)
                        end do
                    end do
                end do
            end do
        end if
        if (associated(D%beammodes)) then
            do i=D%beam_MCMC_modes+1, nmodes
                D%BeamModes(D%cl_lmin:D%cl_lmax,i)=D%BeamModes(D%cl_lmin:D%cl_lmax,i)*D%ClFiducial(1,D%cl_lmin:D%cl_lmax)
            end do
            do i=D%beam_MCMC_modes+1, nmodes
                do l1= D%cl_lmin, D%cl_lmax
                    do l2= D%cl_lmin, D%cl_lmax
                        fullcov(l1 -D%cl_lmin+1, l2 -D%cl_lmin+1 ) = &
                        fullcov(l1 -D%cl_lmin+1, l2 -D%cl_lmin+1 ) +  D%BeamModes(l1,i)*D%BeamModes(l2,i)
                    end do
                end do
            end do
        end if
        if (D%bin_width==1 .and. .not. bin_test) then
            allocate(D%inv_covariance(D%vecsize*D%ncl_used, D%vecsize*D%ncl_used))
            D%inv_covariance = fullcov
            deallocate(fullcov)
            !        D%inv_covariance => fullcov
        else

        allocate(D%inv_covariance(D%nbins*D%ncl_used, D%nbins*D%ncl_used))
        do i=1, D%ncl_used
            do j=1,D%ncl_used
                do x=1, D%nbins
                    do y=1, D%nbins
                        if (i==j .and. y>x) exit
                        asum=0
                        do ix = 1, D%vecsize
                            if (D%binWindows(D%cl_lmin+ix-1,y)/=0) then
                                asum = asum + sum(D%binWindows(:,x)*fullcov((i-1)*D%vecsize+1:i*D%vecsize,(j-1)*D%vecsize+ix)) &
                                * D%binWindows(D%cl_lmin+ix-1,y)
                            end if
                        end do
                        D%inv_covariance((i-1)*D%nbins+x,(j-1)*D%nbins+y) = asum
                        if (i==j) D%inv_covariance((i-1)*D%nbins+y,(j-1)*D%nbins+x) = asum
                    end do
                end do
            end do
        end do
        deallocate(fullcov)
        end if

        call Matrix_inverse(D%inv_covariance)
    else !Not mainMPI
        if (D%bin_width==1 .and. .not. bin_test) then
            allocate(D%inv_covariance(D%vecsize*D%ncl_used, D%vecsize*D%ncl_used))
        else
            allocate(D%inv_covariance(D%nbins*D%ncl_used, D%nbins*D%ncl_used))
        end if
    end if !MainMPI
#ifdef MPI
    call MPI_BCAST(D%inv_covariance,Size(D%inv_covariance),MPI_real, 0, MPI_COMM_WORLD, i)
#endif
    end if

    end if !want high l like

    D%lowl_exact = Ini%Read_Logical('lowl_exact')
    if (D%lowl_exact) then
        if (num_cls==3) call MpiStop('CMBLikes current untested for only 3 C_l')
        S = ReadIniFileName(Ini,'lowl_datafile',dataset_dir)
        D%Lowl%lexact = Ini%Read_Int('lowl_lexact')
        D%Lowl%lmax = Ini%Read_Int('lowl_lmax')
        call CMBLikes_ReadLowlFile(D,S)
    end if

    end subroutine CMBLikes_ReadData


    subroutine CMBLikes_ReadLensingReconData(D, Ini,dataset_dir)
    Type(TCMBLikes) :: D
    class(TIniFile) :: Ini
    character(LEN=*), intent(in) :: dataset_dir
    character(LEN=1024) fname

    if (num_cls_ext ==0) &
    call MpiStop('CMBLikes_ReadLensingReconData: must be compiled with num_cls_ext>0')

    if (Feedback > 1) print *,'CMBLikes_ReadLensingReconData'

    D%cl_phi_lmin = Ini%Read_Int('cl_phi_lmin', D%cl_lmin)
    D%cl_phi_lmax = Ini%Read_Int('cl_phi_lmax', D%cl_lmax)
    D%phi_like_approx = Ini%Read_Int('phi_like_approx',D%like_approx)

    allocate(D%ClPhiHat(D%lensing_recon_ncl,D%cl_phi_lmin:D%cl_phi_lmax))
    allocate(D%ClPhiNoise(D%lensing_recon_ncl,D%cl_phi_lmin:D%cl_phi_lmax))
    call CMBLikes_ReadClPhiArr(D, ReadIniFileName(Ini,'cl_hat_phi_file',dataset_dir),D%ClPhiHat)
    call CMBLikes_ReadClPhiArr(D, ReadIniFileName(Ini,'cl_noise_phi_file',dataset_dir),D%ClPhiNoise)

    if (.not. Ini%Read_Logical('cl_hat_includes_noise')) then
        D%ClPhiHat = D%ClPhiHat + D%ClPhiNoise
    end if

    if (D%phi_like_approx /= like_approx_fullsky_exact) then
        fname = ReadIniFileName(Ini,'covmat_phi_fiducial',dataset_dir)
        if (fname /='') then
            allocate(D%phi_inv_covariance(D%cl_phi_lmin:D%cl_phi_lmax,D%cl_phi_lmin:D%cl_phi_lmax))
            call MatrixSym_Read_Binary(fname, D%phi_inv_covariance)
            call Matrix_Inverse(D%phi_inv_covariance)
        end if
    end if

    if (Feedback > 1) print *, 'CMBLikes_ReadLensingReconData done'

    end subroutine CMBLikes_ReadLensingReconData

    subroutine CMBLikes_ReadClPhiArr(D, aname, Cl)
    Type(TCMBLikes) :: D
    character(LEN=*), intent(in) :: aname
    real(mcp) :: Cl(:,D%cl_phi_lmin:), tmp_arr(D%lensing_recon_ncl)
    integer file_unit
    character(LEN=1024) :: tmp
    integer l, ll

    file_unit = new_file_unit()
    call OpenTxtFile(aname, file_unit)
    Cl=0
    do
        read(file_unit,'(a)',end=1) tmp
        read(tmp,*, end=1) l, tmp_arr
        if (l>=D%cl_phi_lmin .and. l <=D%cl_phi_lmax) then
            ll=l
            Cl(1,l) = tmp_arr(1)
            if (D%lensing_recon_ncl>1) call MpiStop('CMBLikes_ReadClPhiArr: change for n>1')
        end if
    end do
    if (ll<D%cl_phi_lmax) then
        write(*,*) 'CMBLikes_ReadClPhiArr: C_l file does not go up to phi lmax:', D%cl_phi_lmax
        write (*,*) trim(aname)
        call MpiStop()
    end if
1   call CloseFile(file_unit)


    end subroutine CMBLikes_ReadClPhiArr

    subroutine CMBLikes_Transform(D, C, Chat, CfHalf, COffset)
    !Get  C = C_s^{1/2}  U f(D) U^T C_s^{1/2} where C^{-1/2} CHat C^{-1/2} = U D U^T

    !Get  C = C_f^{1/2} C^{-1/2} C^{+1/2} U f(D) U^T C^{+1/2} C^{-1/2} C_f^{1/2} where C^{-1/2} CHat C^{-1/2} = U D U^T

    Type(TCMBLikes) :: D
    real(mcp) C(D%nfields,D%nfields)
    real(mcp), intent(in), optional :: COffset(D%nfields,D%nfields)
    real(mcp), intent(in) :: CHat(D%nfields,D%nfields), CfHalf(D%nfields,D%nfields)
    real(mcp) :: U(D%nfields,D%nfields), Rot(D%nfields,D%nfields)
    real(mcp) :: roots(D%nfields)
    real(mcp) :: diag(D%nfields)
    integer i

    if (present(COffset)) then
        U = C + Coffset*C
        call Matrix_Diagonalize(U,Diag,D%nfields)
        Rot= matmul(matmul(transpose(U),CHat+ COffset*C),U)
    else
        U = C
        call Matrix_Diagonalize(U,Diag,D%nfields)
        Rot= matmul(matmul(transpose(U),CHat),U)
    end if
    roots = sqrt(Diag)

    do i=1, D%nfields
        Rot(i,:)=Rot(i,:)/roots(i)
        Rot(:,i)=Rot(:,i)/roots(i)
    end do

    Rot = matmul(U,matmul(Rot,transpose(U)))
    call Matrix_Diagonalize(Rot,Diag,D%nfields)

    Diag = sign(sqrt(2*max(0._mcp,Diag-log(Diag)-1)),Diag-1)
    !want f(D)-1 to save calculating X-X_s

    if (present(COffset)) then
        Rot = MatMul(transpose(U),Rot)
        do i=1, D%nfields
            Rot(i,:)=Rot(i,:)*roots(i)
        end do
        Rot = MatMul(U,Rot)
        call Matrix_Root(C,D%nfields, -0.5_mcp)
        Rot = MatMul(C,Rot)
    end if

    U = matmul(CfHalf,Rot)
    C = U
    do i=1, D%nfields
        C(:,i) = C(:,i)*Diag(i)
    end do
    C = MatMul(C,transpose(U))

    end subroutine CMBLikes_Transform


    subroutine MatrixToElements(D, M, X)
    Type(TCMBLikes) :: D
    real(mcp) :: M(D%nfields,D%nfields)
    real(mcp) :: X(D%ncl)
    integer ix,i,j

    ix=0
    do i=1, D%nfields
        do j=1,i
            ix = ix+1
            X(ix) = M(i,j)
        end do
    end do

    end subroutine MatrixToElements

    subroutine MatrixToElementsInt(D, M, X)
    Type(TCMBLikes) :: D
    integer, intent(in) :: M(D%nfields,D%nfields)
    integer,intent(out) :: X(D%ncl)
    integer ix,i,j

    ix=0
    do i=1, D%nfields
        do j=1,i
            ix = ix+1
            X(ix) = M(i,j)
        end do
    end do

    end subroutine MatrixToElementsInt


    subroutine ElementsToMatrix(D, X, M)
    Type(TCMBLikes) :: D
    real(mcp), intent(out) :: M(D%nfields,D%nfields)
    real(mcp), intent(in) :: X(D%ncl)
    integer ix,i,j

    ix=0
    do i=1, D%nfields
        do j=1,i
            ix = ix+1
            M(i,j) = X(ix)
            M(j,i) = M(i,j)
        end do
    end do

    end subroutine ElementsToMatrix

    function ExactChiSq(D, C,Chat,l)
    Type(TCMBLikes) :: D
    real(mcp), intent(in) :: C(D%nfields,D%nfields), Chat(D%nfields,D%nfields)
    integer, intent(in) :: l
    real(mcp) ExactChiSq
    real(mcp) M(D%nfields,D%nfields)

    M = C
    call Matrix_root(M,D%nfields,-0.5_mcp)
    M = matmul(M,matmul(Chat,M))
    ExactChiSq = (2*l+1)*D%fullsky_exact_fksy*(Matrix_Trace(M) - D%nfields - MatrixSym_LogDet(M) )

    end function ExactChiSq

    function CMBLikes_LensRecon_Like(D, cl_in) result (chisq)
    Type(TCMBLikes) :: D
    real(mcp), intent(in) :: cl_in(lmax,num_cls_tot)
    real(mcp) vec(D%cl_phi_lmin:D%cl_phi_lmax)
    integer l, phi_ix
    real(mcp) chisq, Cphi, CPhiHat

    if (Feedback > 1) print *,'CMBLikes_LensRecon_CMBLike'

    if (D%bin_width/=1) call MpiStop('CMBLikes_LensRecon_CMBLike: unsupported option')
    if (D%lensing_recon_ncl /=1) call MpiStop('CMBLikes_LensRecon_CMBLike: unsupported ncl')
    phi_ix = num_cls + 1 !put here to avoid compile-time bounds problems
    !not implemented cross-correlation
    chisq = 0
    do l=D%cl_phi_lmin, D%cl_phi_lmax
        Cphi = cl_in(l,phi_ix) + D%ClPhiNoise(1,l)
        CPhihat = D%ClPhihat(1,l)
        if (D%phi_like_approx == like_approx_fullsky_exact) then
            chisq = chisq + (2*l+1)*D%fullsky_exact_fksy*( CPhiHat/CPhi + log(CPhi/CPhiHat) -1)
        else if (D%phi_like_approx == like_approx_fid_gaussian) then
            vec(l) = CPhiHat - CPhi
        else
            call MpiStop('only implemented lensing recon exact like')
        end if
    end do

    if (D%phi_like_approx /= like_approx_fullsky_exact) then
        chisq = chisq + Matrix_QuadForm(D%phi_inv_covariance,vec)
    end if

    if (Feedback > 1) print *,'CMBLikes_LensRecon_CMBLike done'

    end function CMBLikes_LensRecon_Like

    function CMBLikes_CMBLike(D, cl_in) result (chisq)
    Type(TCMBLikes) :: D
    real(mcp), intent(in) :: cl_in(lmax,num_cls_tot)
    real(mcp)  :: cl(lmax,num_cls)

    real(mcp) chisq
    real(mcp) C(D%nfields,D%nfields)
    real(mcp) vecp(D%ncl)
    real(mcp) bigX(D%nbins*D%ncl_used)
    integer l,  i, Ti,Ei,Bi, bin,clix
    !  integer mode
    logical :: quadratic

    chisq =0

    do clix = 1, D%ncl_hat ! 1 or sum over chi-squareds of simulations

    cl = cl_in(:,1:num_cls) !For the moment have not actually implemented lensing likelihood

    if (D%highl_cl) then

    Ti = D%field_index(1)
    Ei = D%field_index(2)
    Bi = D%field_index(3)

    !removed nuisance parameters for the moment
    !if (D%num_nuisance_parameters/=0 .and. Ti /= 0) then
    !  mode=0
    !  if (D%pointsource_MCMC_modes>0) then
    !    mode=mode+1
    !    cl(D%cl_lmin:D%cl_lmax,1) = cl(D%cl_lmin:D%cl_lmax,1) + &
    !      D%point_source_error*nuisance_params(mode)*D%ClPointsources(1,D%cl_lmin:D%cl_lmax)
    !  end if
    !  beamC=0
    !  do i=1, D%beam_MCMC_modes
    !     mode = mode + 1
    !     BeamC = BeamC + nuisance_params(mode) *cl(D%cl_lmin:D%cl_lmax,1)*D%beammodes( D%cl_lmin:D%cl_lmax,i )
    !  end do
    !   cl(D%cl_lmin:D%cl_lmax,1) = cl(D%cl_lmin:D%cl_lmax,1) + beamC
    !
    ! end if


    if (Bi/=0 .and. num_cls<3) call MpiStop('CMBLikes_CMBLike: Need num_cls =4 to use B modes')

    if (D%bin_width/=1 .or. bin_test) then

    if (D%like_approx == like_approx_fullsky_exact) call mpiStop('Planck_like: exact like cannot be binned!')

    do bin = 1, D%nbins
        C=0
        if (Ti/=0) C(Ti,Ti) = sum(D%BinWindows(:,bin)*cl(D%cl_lmin:D%Cl_lmax,1))
        if (Ei/=0) then
            C(Ei,Ei) = sum(D%BinWindows(:,bin)*cl(D%cl_lmin:D%Cl_lmax,3))
            if (Ti/=0) then
                C(Ei,Ti) = sum(D%BinWindows(:,bin)*cl(D%cl_lmin:D%Cl_lmax,2))
                C(Ti,Ei) = C(Ei,Ti)
            end if
        end if
        if (Bi/=0) C(Bi,Bi) = sum(D%BinWindows(:,bin)*cl(D%cl_lmin:D%Cl_lmax,num_cls))
        C = C + D%NoiseM(bin)%M

        if (D%like_approx == like_approx_diag) then
            call CMBLikes_Transform(D, C, D%ChatM(bin,clix)%M, D%sqrt_fiducial(bin)%M)
        else if (D%like_approx == like_approx_fid_gaussian) then
            C = C - D%ChatM(bin,clix)%M
        else if (D%like_approx == like_approx_gaussian) then
            if (D%nfields/=1) call MpiStop('not done Gaussian for polarization')
            !Det term is generally (n+1)[ log C - log C_f])
            chisq = chisq + 2*log( C(1,1) / (D%sqrt_fiducial(bin)%M(1,1))**2 )
            C = (C- D%ChatM(bin,clix)%M) *  (D%sqrt_fiducial(bin)%M)**2 / C
            quadratic = .true.
        else
            call MpiStop('Unknown like_approx')
        end if

        quadratic = .true.
        call MatrixToElements(D, C, vecp)

        do i=1,D%ncl
            if (D%cl_use_index(i)/=0) then
                !            bigX( (D%cl_use_index(i)-1)*D%nbins + bin) = vecp(i)
                bigX( (i-1)*D%nbins + bin) = vecp(i)

            end if
        end do

    end do


    else

    do l = D%cl_lmin, D%cl_lmax
        C=0
        if (Ti/=0) C(Ti,Ti) = cl(l,1)
        if (Ei/=0) then
            C(Ei,Ei) = cl(l,3)
            if (Ti/=0) then
                C(Ei,Ti) = cl(l,2)
                C(Ti,Ei) = C(Ei,Ti)
            end if
        end if
        if (Bi/=0) C(Bi,Bi) =  cl(l,num_cls)

        C =C + D%NoiseM(l)%M

        if (D%like_approx == like_approx_diag) then
            if (Associated(D%OffsetM)) then
                chisq = chisq +2*log( (D%ChatM(l,clix)%M(1,1)+D%OffsetM(l)%M(1,1)*C(1,1))/D%ChatM(l,clix)%M(1,1)/(1+D%OffsetM(l)%M(1,1)))
                call CMBLikes_Transform(D, C, D%ChatM(l,clix)%M, D%sqrt_fiducial(l)%M, D%OffsetM(l)%M )
            else
                call CMBLikes_Transform(D, C, D%ChatM(l,clix)%M, D%sqrt_fiducial(l)%M)
            end if
            call MatrixToElements(D, C, vecp)
            quadratic = .true.
        else if (D%like_approx == like_approx_fid_gaussian) then
            call MatrixToElements(D, C- D%ChatM(l,clix)%M, vecp)
            quadratic = .true.
        else if (D%like_approx == like_approx_gaussian) then
            if (D%nfields/=1) call MpiStop('not done Gaussian for polarization')
            chisq = chisq + 2*log( C(1,1) / (D%sqrt_fiducial(l)%M(1,1))**2 )
            C = (C- D%ChatM(l,clix)%M) *  (D%sqrt_fiducial(l)%M)**2 / C
            call MatrixToElements(D, C, vecp)
            !Det term is generally (n+1)[ log C - log C_f])
            quadratic = .true.

        else if (D%like_approx == like_approx_fullsky_exact) then
            quadratic = .false.
            chisq = chisq  + ExactChisq(D, C,D%ChatM(l,clix)%M,l)
        else
            call MpiStop('Unknown like_approx')
        end if

        if (quadratic) then
            do i=1,D%ncl
                if (D%cl_use_index(i)/=0) then
                    !      bigX( (D%cl_use_index(i)-1)*D%vecsize + l-D%cl_lmin+1) = vecp(i)
                    bigX( (i-1)*D%vecsize + l-D%cl_lmin+1) = vecp(i)

                end if
            end do
        end if

    end do

    end if !binned

    if (quadratic) chisq = chisq + Matrix_QuadForm(D%inv_covariance,BigX)

    end if

    if (D%lowl_exact) then
        chisq = chisq + CMBLikes_lowl_CMBLike(D, cl)
    end if

    if (D%lensing_recon_ncl>0) then
        chisq = chisq + CMBLikes_LensRecon_Like(D, cl_in)
    end if

    end do

    end function CMBLikes_CMBLike


    end module CMBLikes