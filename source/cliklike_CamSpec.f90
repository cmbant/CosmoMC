    module cliklike
    use cmbtypes
    use settings
    use temp_like
    use Likelihood
#ifdef highL
    use HIGHELL_OPTIONS
    use highell_likelihood
#endif
    implicit none

    logical :: use_CAMspec = .true.
    logical :: use_highL = .false.

   type, extends(CosmologyLikelihood) :: CamSpecLikelihood
    contains
    procedure :: LogLike => CamSpecLogLike
    end type CamSpecLikelihood

#ifdef highL
   type, extends(CosmologyLikelihood) :: highLLikelihood
    contains
    procedure :: LogLike => highLLogLike
   end type highLLikelihood
#endif

    integer, parameter :: dp = kind(1.d0)

    private
    public :: clik_readParams

    contains

    
    subroutine clik_readParams(LikeLIst, Ini)
    use settings
    Type(TIniFile) Ini
    class(LikelihoodList) :: LikeList
    character (LEN=Ini_max_string_len) :: likefilename,sz143filename,&
      beamfilename, kszfilename,tszxcibfilename
    Class(CosmologyLikelihood), pointer :: Like

    print *,' using cliklike_CamSpec'
    use_CAMspec = Ini_Read_Logical_File(Ini,'use_CAMspec',.true.)

    if (use_CAMspec) then
        allocate(CamSpecLikelihood::Like)
        call LikeList%Add(Like) 
        Like%needs_powerspectra =.true.
        Like%LikelihoodType = 'CMB'
        Like%name='CamSpec'
        Like%version = CAMSpec_like_version
        Like%speed = 5
        call Like%loadParamNames(trim(DataDir)//'camspec.paramnames')

        likefilename=ReadIniFileName(Ini,'likefile',NotFoundFail = .true.)
        sz143filename=ReadIniFileName(Ini,'sz143file',NotFoundFail = .true.)
        tszxcibfilename=ReadIniFileName(Ini,'tszxcibfile',NotFoundFail = .true.)
        kszfilename=ReadIniFileName(Ini,'kszfile',NotFoundFail = .true.)
        beamfilename=ReadIniFileName(Ini,'beamfile',NotFoundFail = .true.)
        call like_init(likefilename,sz143filename,tszxcibfilename,kszfilename,beamfilename)

    end if 

    use_highL = Ini_Read_Logical_File(Ini,'use_highL',.false.)

    if (use_highL) then
#ifdef highL
        allocate(highLLikelihood::Like)
        call LikeList%Add(Like) 
        Like%LikelihoodType = 'CMB'
        Like%name='highL'
        Like%version = CAMSpec_like_version
        Like%speed = 6
        Like%needs_powerspectra =.true.

        call Like%loadParamNames(trim(DataDir)//'highL.paramnames')
        
      if (lmax < tt_lmax_mc) call MpiStop('set lmax>=tt_lmax_mc in settings to use highL data')
      data_dir = CheckTrailingSlash(ReadIniFileName(Ini,'highL_data_dir'))
      SPT_data_dir = trim(data_dir) // 'data_spt/' 
      ACT_data_dir = trim(data_dir) // 'data_act/'
      if (Feedback>0) write(*,*) 'reading High ell data'
      call highell_likelihood_init
#else
      call MpiStop('must compile with -DhighL to use highL')
#endif
    end if
    
    end subroutine clik_readParams


  real(mcp) function CamspecLogLike(like, CMB, Theory, DataParams) 
     Class(CamspecLikelihood) :: like
     Class (CMBParams) CMB
     Class(TheoryPredictions) Theory
     real(mcp) acl(lmax,num_cls_tot)
     real(mcp) DataParams(:)

     call ClsFromTheoryData(Theory, acl)
!Assuming CAMspec nuisance parameters are set as freq_params(2:34), PLik nuisance parameters as 
!freq_params(35:44), ACT/SPT as freq_params(45:65)
      CamspecLogLike = clik_lnlike_camSpec(acl,DataParams)
 end function CamspecLogLike

 
    function clik_lnlike_camSpec(cl,freq_params)
    real(dp) :: clik_lnlike_camSpec
    real(mcp), intent(in) :: cl(lmax,num_cls_tot)
    real(mcp), intent(in)  :: freq_params(:)
    real(dp) zlike, A_ps_100, A_ps_143, A_ps_217, A_cib_143, A_cib_217, A_sz_143, r_ps, r_cib, &
    cal0, cal1, cal2, xi, A_ksz, ncib

    integer, parameter :: lmin=2
    real(dp) cell_cmb(0:lmax)
    integer, parameter :: nbeams=4, nmodesperbeam=5
    integer, parameter :: nbeammodes = nbeams*nmodesperbeam
    real(mcp) beamcoeffs(nbeammodes)
    real(dp) beam_coeffs(nbeams,nmodesperbeam)
    integer ii,jj,L

    !asz is already removed, start at second feq param
    A_ps_100=freq_params(1)
    A_ps_143 = freq_params(2)
    A_ps_217 = freq_params(3)
    A_cib_143=freq_params(4)
    A_cib_217=freq_params(5) 
    A_sz_143=freq_params(6)  !143
    r_ps = freq_params(7)
    r_cib=freq_params(8)
    ncib = freq_params(9)
    cal0=freq_params(10)
    cal1=freq_params(11) 
    cal2=freq_params(12)
    xi=freq_params(13)
    A_ksz=freq_params(14)
    beamcoeffs = freq_params(15:15+nbeammodes-1)

    do L=lmin,lmax
        cell_cmb(L)=cl(L,1)/twopi !this is a georgeism
    enddo

    do ii=1,nbeams
        do jj=1,nmodesperbeam
            beam_coeffs(ii,jj)=beamcoeffs(jj+nmodesperbeam*(ii-1))
        enddo
    enddo

    call calc_like(zlike,  cell_cmb, A_ps_100,  A_ps_143, A_ps_217, A_cib_143, A_cib_217, A_sz_143,  &
    r_ps, r_cib, xi, A_ksz, ncib, cal0, cal1, cal2, beam_coeffs)

    clik_lnlike_camSpec = zlike/2

    if (Feedback>2) Print*,'CamSpec lnlike = ',clik_lnlike_camSpec

    end function clik_lnlike_camSpec

#ifdef highL
  real(mcp) function highLLogLike(like, CMB, Theory, DataParams) 
     Class(highLLikelihood) :: like
     Class (CMBParams) CMB
     Class(TheoryPredictions) Theory
     real(mcp) acl(lmax,num_cls_tot)
     real(mcp) DataParams(:)

     call ClsFromTheoryData(Theory, acl)
     highLLogLike = clik_lnlike_highL(acl,DataParams)

  end function highLLogLike


  function clik_lnlike_highL(cl,freq_params)
    real(dp) :: clik_lnlike_highL
    real(mcp), intent(in) :: cl(lmax,num_cls_tot)
    real(mcp), intent(in)  :: freq_params(:)
    real(dp) like_tot
    integer, parameter :: lmin=2
    real(dp)  cl_tt(tt_lmax)
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
    clik_lnlike_highL = like_tot

    if (Feedback>2) Print*,'highL lnlike = ',clik_lnlike_highL

    end function clik_lnlike_highL
#endif

   end module cliklike

