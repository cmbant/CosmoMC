    module noncliklike
    use CosmologyTypes
    use settings
    use temp_like_camspec4
    use Likelihood_Cosmology
    use CosmoTheory
#ifdef highL
    use HIGHELL_OPTIONS
    use highell_likelihood
#endif
    implicit none

    type, extends(TCosmologyLikelihood) :: CamSpeclikelihood
    contains
    procedure :: LogLike => CamSpecLogLike
    procedure :: WriteLikelihoodData => CamSpec_WriteLikelihoodData
    procedure :: derivedParameters => CamSpec_DerivedParameters
    end type CamSpeclikelihood

#ifdef highL
    type, extends(TCosmologyLikelihood) :: highLLikelihood
    contains
    procedure :: LogLike => highLLogLike
    end type highLLikelihood
#endif

    integer, parameter :: dp = kind(1.d0)

    private
    public :: nonclik_readParams

    contains


    subroutine nonclik_readParams(LikeLIst, Ini)
    class(TSettingIni) Ini
    class(TLikelihoodList) :: LikeList
    character (LEN=:), allocatable :: likefilename,sz143filename,&
        beamfilename, kszfilename,tszxcibfilename, tmp, polfilename,&
        data_vector, cib217filename, dust100filename, dust143filename,&
        dust217filename, dust143x217filename
    Class(TCosmologyLikelihood), pointer :: Like
    logical :: use_CAMspec, use_CAMpol
    logical :: use_highL
    logical :: pre_marged
    integer i

    use_CAMspec = Ini%Read_Logical('use_CAMspec',.false.)

    if (use_CAMspec) then
        print *,' using non-clik CamSpec'
        allocate(CamSpeclikelihood::Like)
        call LikeList%Add(Like)
        Like%needs_powerspectra =.true.
        Like%LikelihoodType = 'CMB'
        Like%name='CamSpec'
        Like%speed = 5
        allocate(like%cl_lmax(CL_E,CL_E),source=0)
        Like%cl_lmax(CL_T,CL_T) = 2500
        Like%cl_lmax(CL_E,CL_T) = 2500
        Like%cl_lmax(CL_E,CL_E) = 2500

        apply_tight_sz_prior=Ini%Read_Logical('apply_tight_sz_prior',.false.)
        make_cov_marged = Ini%Read_Logical('make_cov_marged',.false.)
        if (make_cov_marged) marge_file_variant = Ini%Read_string('marge_file_variant')
        pre_marged= .not. make_cov_marged .and. Ini%Read_Logical('pre_marged',.false.)
        if (pre_marged) then
            likefilename=Ini%ReadFileName('margelikefile',NotFoundFail = .true.)
            if (camspec_beam_mcmc_num/=1) call MpiStop('camspec_beam_mcmc_num must be one for precomputed')
        else
            tmp = Ini%Read_String('want_spec')
            if (tmp/='') read(tmp,*) want_spec
            tmp = Ini%Read_String('camspec_lmin')
            if (tmp/='') read(tmp,*) camspec_lmins
            tmp = Ini%Read_String('camspec_lmax')
            if (tmp/='') read(tmp,*) camspec_lmaxs
            camspec_ranges = Ini%Read_String('camspec_ranges')
           
            likefilename=Ini%ReadFileName('likefile',NotFoundFail = .true.)
            camspec_fiducial_foregrounds = Ini%ReadFileName('camspec_fiducial_foregrounds',NotFoundFail = .true.)
            camspec_fiducial_cl = Ini%ReadFileName('camspec_fiducial_cl',NotFoundFail = .true.)
        end if
        call Like%loadParamNames(trim(DataDir)//'camspec_fullbeam.paramnames')
        Like%version = File%ExtractName(likefilename)

        sz143filename=Ini%ReadFileName('sz143file',NotFoundFail = .true.)
        tszxcibfilename=Ini%ReadFileName('tszxcibfile',NotFoundFail = .true.)
        kszfilename=Ini%ReadFileName('kszfile',NotFoundFail = .true.)
        beamfilename=Ini%ReadFileName('beamfile',NotFoundFail = .true.)
        data_vector = Ini%ReadFileName('camspec_data_vector')
        if (data_vector/='') Like%version = File%ExtractName(data_vector)
        cib217filename=Ini%ReadFileName('cib217file')
        dust100filename=Ini%ReadFileName('dust100file')
        dust143filename=Ini%ReadFileName('dust143file')
        dust217filename=Ini%ReadFileName('dust217file')
        dust143x217filename=Ini%ReadFileName('dust143x217file')
        dump_spectra = Ini%Read_Logical('dump_spectra',.false.)

        call Ini%Read('camspec_beam_mcmc_num',camspec_beam_mcmc_num)

        call like_init(pre_marged,likefilename,sz143filename,tszxcibfilename,kszfilename,beamfilename,data_vector,&
            cib217filename,dust100filename,dust143filename,dust217filename,dust143x217filename)

    end if

    use_highL = Ini%Read_Logical('use_highL',.false.)

    if (use_highL) then
#ifdef highL
        print *,' using non-clik highL'
        allocate(highLLikelihood::Like)
        call LikeList%Add(Like)
        Like%LikelihoodType = 'CMB'
        Like%name='highL'
        Like%version = CAMSpec_like_version
        Like%speed = 6
        Like%needs_powerspectra =.true.
        allocate(like%cl_lmax(1,1))
        Like%cl_lmax(CL_T,CL_T) = 6000

        call Like%loadParamNames(trim(DataDir)//'highL.paramnames')

        if (lmax < tt_lmax_mc) call MpiStop('set lmax>=tt_lmax_mc in settings to use highL data')
        data_dir = File%CheckTrailingSlash(Ini%ReadFileName('highL_data_dir'))
        SPT_data_dir = trim(data_dir) // 'data_spt/'
        ACT_data_dir = trim(data_dir) // 'data_act/'
        if (Feedback>0) write(*,*) 'reading High ell data'
        call highell_likelihood_init
#else
        call MpiStop('must compile with -DhighL to use highL')
#endif
    end if

    end subroutine nonclik_readParams

    real(mcp) function CamspecLogLike(this, CMB, Theory, DataParams)
    Class(CamSpeclikelihood) :: this
    Class (CMBParams) CMB
    Class(TCosmoTheoryPredictions), target :: Theory
    real(mcp) DataParams(:)
    real(dp) zlike, cell_cmb(0:this%cl_lmax(CL_T,CL_T)),cell_cmbx(0:this%cl_lmax(CL_E,CL_T)), cell_cmbe(0:this%cl_lmax(CL_E,CL_E))
    integer L

    !Assuming CAMspec nuisance parameters are set as freq_params(2:34), PLik nuisance parameters as
    !freq_params(35:44), ACT/SPT as freq_params(45:65)

    do L=2,this%cl_lmax(CL_T,CL_T)
        cell_cmb(L)=Theory%Cls(CL_T,CL_T)%Cl(L)/(L*(L+1)) !normalization is a georgeism
    enddo
    do L=2,this%cl_lmax(CL_E,CL_T)
        cell_cmbx(L)=Theory%Cls(2,1)%Cl(L)/(L*(L+1)) !this is a georgeism
    enddo
    do L=2,this%cl_lmax(CL_E,CL_E)
        cell_cmbe(L)=Theory%Cls(2,2)%Cl(L)/(L*(L+1)) !this is a georgeism
    enddo

    call calc_like(zlike,  cell_cmb, cell_cmbx, cell_cmbe, DataParams)

    CamspecLogLike = zlike/2

    if (Feedback>2) Print*,'CamSpec lnlike = ',CamspecLogLike

    end function CamspecLogLike


    subroutine CAMSpec_WriteLikelihoodData(this,Theory,DataParams, root)
    Class(CamSpeclikelihood) :: this
    class(TTheoryPredictions) :: Theory
    real(mcp), intent(in) :: DataParams(:)
    character(LEN=*), intent(in) :: root
    real(mcp) , allocatable ::  C_foregrounds(:,:)
    integer L
    Type(TTextFile) F

    if (any(want_spec(1:4))) then
        allocate(C_foregrounds(CAMSpec_lmax_foreground,Nspec))
        C_foregrounds=0
        call compute_fg(C_foregrounds,DataParams, CAMSpec_lmax_foreground)
        call F%CreateFile(trim(root)//'.camspec_foregrounds')
        do l = 2, CAMSpec_lmax_foreground
            write (F%unit,'(1I6,*(E15.5))') l, C_foregrounds(l,:)*l*(l+1)
        end do
        call F%Close()
    end if

    end subroutine CAMSpec_WriteLikelihoodData

    function CamSpec_derivedParameters(this, Theory, DataParams) result(derived)
    class(CamSpeclikelihood) :: this
    class(TTheoryPredictions) :: Theory
    real(mcp) :: derived(this%nuisance_params%num_derived)
    real(mcp) :: DataParams(:)
    real(mcp) , allocatable ::  C_foregrounds(:,:)
    integer l

    l=2000
    allocate(C_foregrounds(l,Nspec))
    C_foregrounds=0
    call compute_fg(C_foregrounds,DataParams, l)
    derived(1) = C_foregrounds(2000,2)*l*(l+1)
    derived(2) = C_foregrounds(2000,3)*l*(l+1)
    derived(3) = C_foregrounds(2000,4)*l*(l+1)

    end function CamSpec_derivedParameters


#ifdef highL
    real(mcp) function highLLogLike(this, CMB, Theory, DataParams)
    Class(highLLikelihood) :: this
    Class (CMBParams) CMB
    Class(TheoryPredictions) Theory
    real(mcp) acl(lmax,num_cls_tot)
    real(mcp) DataParams(:)

    call ClsFromTheoryData(Theory, acl)
    highLLogLike = nonclik_lnlike_highL(acl,DataParams)

    end function highLLogLike


    function nonclik_lnlike_highL(cl,freq_params)
    real(dp) :: nonclik_lnlike_highL
    real(mcp), intent(in) :: cl(lmax,num_cls_tot)
    real(mcp), intent(in)  :: freq_params(:)
    real(dp) like_tot
    integer, parameter :: lmin=2
    real(dp)  cl_tt(2:tt_lmax)
    integer L, offset
    real(dp) A_ps_100, A_ps_143, A_ps_217, A_cib_143, A_cib_217, A_sz_143, r_ps, r_cib, &
        cal0, cal1, cal2, xi, A_ksz, ncib
    real(dp) a_ps_act_148,a_ps_act_217,a_ps_spt_95,a_ps_spt_150,a_ps_spt_220, &
        r_ps_spt_95x150,r_ps_spt_95x220,r_ps_150x220, &
        cal_acts_148,cal_acts_217,cal_acte_148,cal_acte_217,cal_spt_95,cal_spt_150,cal_spt_220
    real(dp) act_dust_s,act_dust_e

    !asz is already removed, start at second feq param
    A_sz_143=freq_params(1)
    A_ksz = freq_params(2)
    xi = freq_params(3)
    a_ps_act_148=freq_params(4)
    a_ps_act_217=freq_params(5)
    a_ps_spt_95=freq_params(6)
    a_ps_spt_150 =freq_params(7)
    a_ps_spt_220 = freq_params(8)
    A_cib_143=freq_params(9)
    A_cib_217=freq_params(10)
    ncib = freq_params(11)
    r_ps_spt_95x150=freq_params(12)
    r_ps_spt_95x220=freq_params(13)
    r_ps_150x220=freq_params(14)
    r_cib=freq_params(15)
    act_dust_s = freq_params(16)
    act_dust_e = freq_params(17)
    cal_acts_148  =freq_params(18)
    cal_acts_217=freq_params(19)
    cal_acte_148=freq_params(20)
    cal_acte_217 =freq_params(21)
    cal_spt_95 =freq_params(22)
    cal_spt_150 =freq_params(23)
    cal_spt_220 =freq_params(24)

    do l =2, tt_lmax
        if (l.le.tt_lmax_mc) then
            cl_tt(l) = cl(l,1)*l*(l+1)/twopi
        else
            cl_tt(l) = 0.0d0
        endif
    enddo

    like_tot = 0.d0
    call highell_likelihood_compute(cl_tt,A_sz_143,A_ksz,xi,a_ps_act_148,a_ps_act_217,a_ps_spt_95,a_ps_spt_150,a_ps_spt_220, &
        A_cib_143,A_cib_217, ncib, r_ps_spt_95x150,r_ps_spt_95x220,r_ps_150x220,r_cib,act_dust_s,act_dust_e, &
        cal_acts_148,cal_acts_217,cal_acte_148,cal_acte_217,cal_spt_95,cal_spt_150,cal_spt_220,like_tot)
    nonclik_lnlike_highL = like_tot

    if (Feedback>2) Print*,'highL lnlike = ',nonclik_lnlike_highL

    end function nonclik_lnlike_highL
#endif

    end module noncliklike
