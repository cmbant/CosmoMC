
    module CMBLikelihoods
    use Likelihood_Cosmology
    use CosmologyTypes
    use CosmoTheory
    use CMBLikes
    implicit none
    private

    type, extends(TCMBLikelihood) :: TCMBSZLikelihood
        ! (inherited) tag set from "cmb_dataset[tag] =" in input file
        real(mcp), pointer, dimension(:) :: sz_template
    contains
    procedure :: ReadSZTemplate
    procedure :: ReadParams => TCMBSZLikelihood_ReadParams
    end type TCMBSZLikelihood

#ifdef WMAP
    type, extends(TCMBSZLikelihood) :: TWMAPLikelihood
    contains
    procedure :: ReadParams => TWMAPLikelihood_ReadParams
    procedure :: LogLike => TWMAPLikelihood_LogLike
    end type TWMAPLikelihood
#endif

    type index_array
        integer, allocatable :: bins(:)
    end type index_array

    type, extends(TCMBLikelihood) :: TPlikLiteLikelihood
        !Plik-lite Planck likelihood (unofficial native cosmomc version)
        !Adapted from code by Erminia Calabrese
        integer :: plmin = 30
        integer :: lmax = 2508
        integer :: nbins = 613
        integer, dimension(3) :: nbincl = (/ 215, 199, 199/)
        integer :: maxbin, nused
        integer, allocatable :: blmin(:), blmax(:), ls(:)
        logical :: used(3)
        character(LEN=2), dimension(3) :: cl_names = ['TT','TE','EE']
        real(mcp), allocatable :: weights(:), invcov(:,:), X_data(:)
        Type(index_array) :: used_bins(3)
        integer pairs(3,2)
    contains
    procedure :: ReadIni => TPlikLiteLikelihood_ReadIni
    procedure :: LogLike => TPlikLiteLikelihood_LogLike
    end type TPlikLiteLikelihood


    public CMBLikelihood_Add
    contains


    subroutine CMBLikelihood_Add(LikeList, Ini)
#ifdef CLIK
    use cliklike
#endif
#ifdef NONCLIK
    use noncliklike
#endif
    use BK_planck
    use smica_planck
    class(TLikelihoodList) :: LikeList
    class(TSettingIni) :: ini
    class(TCMBLikelihood), pointer  :: like
    integer  i
    Type(TSettingIni) :: DataSets, OverrideSettings

    call Ini%TagValuesForName('cmb_dataset', DataSets, filename=.true.)

    do i= 1, DataSets%Count
        call Ini%SettingValuesForTagName('cmb_dataset',DataSets%Name(i),OverrideSettings)
        if (DataSets%Name(i) == 'WMAP') then
#ifdef WMAP
            allocate(TWMAPLikelihood::like)
            like%name = Datasets%Value(i)
            select type(like)
            class is (TWMAPLikelihood)
            end select
#else
            call MpiStop('Set WMAP directory in Makefile to compile with WMAP')
#endif
        else
            if (DataSets%Name(i) == 'BKPLANCK') then
                allocate(TBK_planck::like)
            else if (DataSets%Name(i) == 'SMICA') then
                allocate(TSmica_planck::like)
            else if (DataSets%Name(i) == 'PLIK_LITE') then
                allocate(TPlikLiteLikelihood::like)
            else
                allocate(TCMBLikes::like)
            end if
            call like%ReadDatasetFile(Datasets%Value(i),OverrideSettings)
        end if
        like%Tag = DataSets%Name(i)
        call like%ReadParams(Ini)
        call Ini%Read(Ini%NamedKey('cmb_dataset_speed',DataSets%Name(i)),like%speed)

        call LikeList%Add(like)
    end do
    if (Feedback > 1 .and. DataSets%Count>0 ) write (*,*) 'read CMB data sets'

#ifdef CLIK
    Use_clik = Ini%Read_Logical('use_clik',.false.)
    if (Use_clik) then
        call clik_readParams(LikeList, Ini)
    end if
#else
    if (Ini%Read_Logical('use_clik',.false.)) call MpiStop('compile with CLIK to use clik - see Makefile')
#endif
#ifdef NONCLIK
    call nonclik_readParams(LikeList, Ini)
#endif

    end subroutine CMBLikelihood_Add


    subroutine TCMBSZLikelihood_ReadParams(this, Ini)
    class(TCMBSZLikelihood) :: this
    class(TSettingIni) :: Ini
    character(LEN=:), allocatable :: SZTemplate
    real(mcp) :: SZscale = 1

    SZTemplate = Ini%Read_String(Ini%NamedKey('cmb_dataset_SZ',this%tag))
    if (SZTemplate/='') then
        call Ini%Read(Ini%NamedKey('cmb_dataset_SZ_scale',this%tag),SZScale)
        call this%ReadSZTemplate(SZTemplate,SZScale)
        call this%loadParamNames(trim(DataDir)//'WMAP.paramnames')
    end if

    call this%TCMBLikelihood%ReadParams(Ini)

    end subroutine TCMBSZLikelihood_ReadParams

    subroutine ReadSZTemplate(this, aname, ascale)
    class(TCMBSZLikelihood) :: this
    character(LEN=*), intent(IN) :: aname
    real(mcp), intent(in) :: ascale
    integer l, status
    real(mcp) sz
    Type(TTextFile) :: F

    allocate(this%sz_template(2:this%cl_lmax(CL_T,CL_T)))
    this%sz_template = 0
    call F%Open(aname)
    do
        read(F%unit,*,iostat=status) l, sz
        if (status/=0) exit
        if (l>=2 .and. l<=this%cl_lmax(CL_T,CL_T)) this%sz_template(l) = ascale * sz
    end do
    call F%Close()

    end subroutine ReadSZTemplate

#ifdef WMAP
    subroutine TWMAPLikelihood_ReadParams(this, Ini)
    use WMAP_OPTIONS
    class(TWMAPLikelihood) :: this
    class(TSettingIni) :: ini

    use_TT_beam_ptsrc = Ini%read_Logical('use_WMAP_TT_beam_ptsrc', .true.)
    use_TE = Ini%read_Logical('use_WMAP_TE',.true.)
    use_TT = Ini%read_Logical('use_WMAP_TT',.true.)
    if (MPIRank==0) print *, 'WMAP options (beam TE TT)', use_TT_beam_ptsrc, use_TE, use_TT
    allocate(this%cl_lmax(CL_B,CL_B), source=0)
    this%cl_lmax(CL_T,CL_T) = ttmax
    this%cl_lmax(CL_E,CL_T) = temax
    this%cl_lmax(CL_E,CL_E) = max(gibbs_ell_max,lowl_max)
    this%cl_lmax(CL_B,CL_B) = max(gibbs_ell_max,lowl_max)

    call this%TCMBSZLikelihood%ReadParams(Ini)

    end subroutine TWMAPLikelihood_ReadParams

    function TWMAPLikelihood_LogLike(this, CMB, Theory, DataParams) result(logLike)
    use wmap_likelihood_9yr
    use WMAP_OPTIONS
    use WMAP_UTIL
    Class(TWMAPLikelihood) :: this
    Class (CMBParams) CMB
    Class(TCosmoTheoryPredictions), target :: Theory
    real(mcp) DataParams(:)
    real(mcp) logLike
    real(8) :: likes(num_WMAP),like_tot
    real(mcp) CLTT(2:this%cl_lmax(1,1))

    CLTT = Theory%Cls(1,1)%CL(2:this%cl_lmax(1,1)) + DataParams(1)*this%sz_template
    likes=0
    call wmap_likelihood_compute(CLTT,Theory%Cls(2,1)%CL(2:),Theory%Cls(2,2)%CL(2:),Theory%Cls(3,3)%CL(2:),likes)
    !call wmap_likelihood_error_report

    if (wmap_likelihood_ok) then
        LogLike = sum(likes)
    else
    endif
    end function TWMAPLikelihood_LogLike
#endif


    subroutine TPlikLiteLikelihood_ReadIni(this, Ini)
    use MatrixUtils
    class(TPlikLiteLikelihood) :: this
    class(TSettingIni) :: Ini
    real(mcp), allocatable :: cov(:,:), dat(:,:), ls(:), weights(:)
    integer i,j, lun, maxbin, rangemin, rangemax
    character(LEN=:), allocatable :: bins_for_L_range, use_cl
    real(mcp) bin_centre
    integer, allocatable :: usebins(:), used_indices(:)
    integer nused, offset, offused
    Type(TStringList) :: used_cls

    call this%loadParamNames(Ini%ReadRelativeFileName('calibration_param'))

    use_cl = Ini%Read_String('use_cl')
    call File%LoadTxt(Ini%ReadRelativeFilename('data'), dat)
    call File%LoadTxt(Ini%ReadRelativeFilename('blmin'),this%blmin)
    this%blmin=this%blmin+this%plmin
    call File%LoadTxt(Ini%ReadRelativeFilename('blmax'),this%blmax)
    this%blmax=this%blmax+this%plmin
    call File%LoadTxt(Ini%ReadRelativeFilename('weights'), weights)
    allocate(ls(size(weights)))
    do i=1, size(weights)
        ls(i) = this%plmin + i -1
        weights(i) = weights(i)*twopi/ls(i)/(ls(i)+1)
    end do
    allocate(this%weights(this%plmin:this%plmin+size(weights)-1), source=weights)

    if (Ini%Read_String('cov_file_binary')/= '') then
        allocate(cov(this%nbins,this%nbins))
        open(newunit=lun,file=Ini%ReadRelativeFilename('cov_file_binary'),form='unformatted',status='old')
        read(lun) cov
        close(lun)
        do i=1,this%nbins
            do j=i+1,this%nbins
                cov(j,i)=cov(i,j)
            enddo
        enddo
    else
        call File%LoadTxt(Ini%ReadRelativeFilename('cov_file'), cov)
    end if
    maxbin = maxval(this%nbincl)
    bins_for_L_range = Ini%Read_string('bins_for_L_range')
    if (bins_for_L_range/='') then
        allocate(usebins(maxbin))
        nused=0
        read(bins_for_L_range, *) rangemin, rangemax
        do i =1, maxbin
            bin_centre = (this%blmin(i)+this%blmax(i))/2.0_mcp
            if (rangemin <= bin_centre .and. bin_centre<= rangemax) then
                nused = nused+1
                usebins(nused) = i
            end if
        end do
        print *, 'bins_for_L_range: actual L range used: ', this%blmin(usebins(1)),this%blmax(usebins(nused))
    end if

    call used_cls%SetFromString(use_cl)
    offset = 0
    this%nused =0
    do i=1,3
        this%used(i) = used_cls%IndexOf(this%cl_names(i))/=-1
        if (this%used(i)) then
            if (allocated(usebins)) then
                allocate(this%used_bins(i)%bins, source =  pack(usebins(1:nused), usebins(1:nused) <= this%nbincl(i)))
            else
                allocate(this%used_bins(i)%bins, source = (/ (j,j=1, this%nbincl(i))/))
            end if
            this%nused = this%nused + size(this%used_bins(i)%bins)
        end if
        offset = offset + this%nbincl(i)
    end do

    this%pairs(1,:) = (/ 1,1/)
    this%pairs(2,:) = (/ 2,1/)
    this%pairs(3,:) = (/ 2,2/)

    allocate(used_indices(this%nused))
    offset =0
    offused = 1

    allocate(this%cl_lmax(CL_B,CL_B), source=0)
    do i=1,3
        if (this%used(i)) then
            used_indices(offused:offused + size(this%used_bins(i)%bins)-1)=this%used_bins(i)%bins+offset
            offused=offused + size(this%used_bins(i)%bins)
            this%cl_lmax(this%pairs(i,1),this%pairs(i,2)) = maxval(this%blmax(this%used_bins(i)%bins))
        end if
        offset = offset + this%nbincl(i)
    end do
    allocate(this%X_data, source = dat(used_indices,2))
    allocate(this%invcov, source = cov(used_indices,used_indices))
    call Matrix_Inverse(this%invcov)

    call this%TCMBLikelihood%ReadIni(Ini)
    end subroutine TPlikLiteLikelihood_ReadIni

    function TPlikLiteLikelihood_LogLike(this, CMB, Theory, DataParams) result(logLike)
    use MatrixUtils
    Class(TPlikLiteLikelihood) :: this
    Class(CMBParams) CMB
    Class(TCosmoTheoryPredictions), target :: Theory
    real(mcp) DataParams(:)
    real(mcp) logLike
    real(mcp) cl(this%nused)
    integer ix, i, j, bini

    ix=0
    do i=1,3
        if (this%used(i)) then
            do j=1, size(this%used_bins(i)%bins)
                bini = this%used_bins(i)%bins(j)
                ix=ix+1
                cl(ix) = dot_product(Theory%Cls(this%pairs(i,1),this%pairs(i,2))%CL(this%blmin(bini):this%blmax(bini)), &
                    this%weights(this%blmin(bini):this%blmax(bini)))
            end do
        end if
    end do
    cl = cl/DataParams(1)**2
    logLike = Matrix_quadform(this%invcov,this%X_data - cl)/2

    end function TPlikLiteLikelihood_LogLike


    end module CMBLikelihoods
