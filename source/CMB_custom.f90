    !Things that used to be in CMBdata, for dated likelihoods
    module CMBCustom
    use Likelihood_Cosmology
    use CosmologyTypes
    use CosmoTheory
    implicit none
    private

    type, extends(TCMBLikelihood) :: TCMBDatasetLikelihood
        character(LEN=:), allocatable :: tag ! from "cmb_dataset[tag] =" in input file
        real(mcp), pointer, dimension(:) :: sz_template
    contains
    procedure :: ReadSZTemplate
    procedure :: ReadParams => TCMBDatasetLikelihood_ReadParams
    end type TCMBDatasetLikelihood

#ifdef WMAP
    type, extends(TCMBDatasetLikelihood) :: TWMAPLikelihood
    contains
    procedure :: ReadParams => TWMAPLikelihood_ReadParams
    procedure :: LogLike => TWMAPLikelihood_LogLike
    end type TWMAPLikelihood
#endif

    public CMBCustomLikelihoods_Add
    contains


    subroutine CMBCustomLikelihoods_Add(LikeList, Ini)
#ifdef CLIK
    use cliklike
#endif
#ifdef NONCLIK
    use noncliklike
#endif
    class(TLikelihoodList) :: LikeList
    class(TSettingIni) :: ini

    class(TCMBDatasetLikelihood), pointer  :: like
    integer  i
    character(LEN=:), allocatable :: tag, SZTemplate
    real(mcp) SZScale
    Type(TSettingIni) :: DataSets

    call Ini%TagValuesForName('cmb_dataset', DataSets)

    do i= 1, DataSets%Count
        if (DataSets%Name(i) == 'WMAP') then
#ifdef WMAP
            allocate(TWMAPLikelihood::like)
            like%name = Datasets%Value(i)
#else
            call MpiStop('Set WMAP directory in Makefile to compile with WMAP')
#endif
        else
            call MpiStop('Support for general CMB data sets not added yet')
            allocate(TCMBDatasetLikelihood::like)
            call like%ReadDatasetFile(Datasets%Value(i))
        end if
        like%Tag = DataSets%Name(i)
        call like%ReadParams(Ini)
        call LikeList%Add(like)
    end do
    if (Feedback > 1) write (*,*) 'CMB datasets read'

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


    end subroutine CMBCustomLikelihoods_Add


    subroutine TCMBDatasetLikelihood_ReadParams(this, Ini)
    class(TCMBDatasetLikelihood) :: this
    class(TSettingIni) :: ini
    integer i
    character(LEN=:), allocatable :: SZTemplate
    real(mcp) :: SZscale = 1

    SZTemplate = Ini%Read_String(Ini%NamedKey('cmb_dataset_SZ',this%tag))
    if (SZTemplate/='') then
        call Ini%Read(Ini%NamedKey('cmb_dataset_SZ_scale',this%tag),SZScale)
        call this%ReadSZTemplate(SZTemplate,SZScale)
        call this%loadParamNames(trim(DataDir)//'WMAP.paramnames')
    end if

    call this%TCMBLikelihood%ReadParams(Ini)

    end subroutine TCMBDatasetLikelihood_ReadParams

    subroutine ReadSZTemplate(this, aname, ascale)
    class(TCMBDatasetLikelihood) :: this
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

    call this%TCMBDatasetLikelihood%ReadParams(Ini)

    end subroutine TWMAPLikelihood_ReadParams

    function TWMAPLikelihood_LogLike(like, CMB, Theory, DataParams) result(logLike)
    use wmap_likelihood_9yr
    use WMAP_OPTIONS
    use WMAP_UTIL
    Class(TWMAPLikelihood) :: like
    Class (CMBParams) CMB
    Class(TCosmoTheoryPredictions), target :: Theory
    real(mcp) DataParams(:)
    real(mcp) logLike
    real(8)                     :: likes(num_WMAP),like_tot
    real(mcp) CLTT(2:like%cl_lmax(1,1))

    CLTT = Theory%Cls(1,1)%CL(2:like%cl_lmax(1,1)) + DataParams(1)*like%sz_template
    likes=0
    call wmap_likelihood_compute(CLTT,Theory%Cls(2,1)%CL(2:),Theory%Cls(2,2)%CL(2:),Theory%Cls(3,3)%CL(2:),likes)
    !call wmap_likelihood_error_report

    if (wmap_likelihood_ok) then
        LogLike = sum(likes)
    else
        LogLike = LogZero
    endif
    end function TWMAPLikelihood_LogLike
#endif



    end module CMBCustom