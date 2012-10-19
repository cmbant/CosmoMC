    module cliklike
    use cmbtypes
    use settings
    use temp_like
#ifdef highL
    use HIGHELL_OPTIONS
    use highell_likelihood
#endif
    implicit none

    logical :: use_clik = .false.  
    logical :: use_CAMspec = .true.
    logical :: use_highL = .false.

    integer, parameter :: dp = kind(1.d0)

    integer, parameter :: num_camSpec =33
    integer, parameter :: num_plik=10
    integer, parameter :: num_actSpt = 21

    
    private
    public :: clik_readParams,clik_lnlike, use_clik

    contains

    subroutine clik_readParams(Ini)
    use settings
    Type(TIniFile) Ini
    character (LEN=Ini_max_string_len) :: likefilename,sz143filename,&
    beamfilename, kszfilename,tszxcibfilename, highLdir

    print *,' using cliklike_CamSpec'
    use_CAMspec = Ini_Read_Logical_File(Ini,'use_CAMspec',.true.)

    if (use_CAMspec) then
        
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
      if (lmax < tt_lmax_mc) call MpiStop('set lmax>=tt_lmax_mc in settings to use highL data')
      data_dir = AddTrailingSlash(ReadIniFileName(Ini,'highL_data_dir'))
      SPT_data_dir = trim(data_dir) // 'data_act/' 
      ACT_data_dir = trim(data_dir) // 'data_spt/'
      if (Feedback>0) write(*,*) 'reading High ell data'
      call highell_likelihood_init
#else
      call MpiStop('must compile with -DhighL to use highL')
#endif
    end if
    
    end subroutine clik_readParams

    function clik_lnlike(cl,freq_params)
    real(dp) :: clik_lnlike
    real(dp), intent(in) :: cl(lmax,num_cls_tot)
    real(dp), intent(in) :: freq_params(1:num_freq_params-1)   
    
    clik_lnlike = 0
    if (use_CAMspec)  clik_lnlike=clik_lnlike+clik_lnlike_camSpec(cl,freq_params)
#ifdef highL
    if (use_highL)    clik_lnlike=clik_lnlike+clik_lnlike_highL(cl,freq_params)
#endif
    end function clik_lnlike

    function clik_lnlike_camSpec(cl,freq_params)
    real(dp) :: clik_lnlike_camSpec
    real(dp), intent(in) :: cl(lmax,num_cls_tot)
    real(dp), intent(in)  :: freq_params(1:num_freq_params-1)   
    real(dp) zlike, A_ps_100, A_ps_143, A_ps_217, A_cib_143, A_cib_217, A_sz_143, r_ps, r_cib, &
    cal0, cal1, cal2, xi, A_ksz

    integer, parameter :: lmin=2
    real(dp) cell_cmb(0:10000)
    integer, parameter :: nbeams=4, nmodesperbeam=5
    integer, parameter :: nbeammodes = nbeams*nmodesperbeam
    real beamcoeffs(nbeammodes)
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
    cal0=freq_params(9)
    cal1=freq_params(10) 
    cal2=freq_params(11)
    xi=freq_params(12)
    A_ksz=freq_params(13)

    do L=lmin,lmax
        cell_cmb(L)=cl(L,1)/twopi !this is a georgeism
    enddo

    do ii=1,nbeams
        do jj=1,nmodesperbeam
            beam_coeffs(ii,jj)=beamcoeffs(jj+nmodesperbeam*(ii-1))
        enddo
    enddo

    call calc_like(zlike,  cell_cmb, A_ps_100,  A_ps_143, A_ps_217, A_cib_143, A_cib_217, A_sz_143,  &
    r_ps, r_cib, xi, A_ksz, cal0, cal1, cal2, beam_coeffs)

    clik_lnlike_camSpec = zlike/2

    if (Feedback>0) Print*,'CamSpec lnlike = ',clik_lnlike_camSpec

    end function clik_lnlike_camSpec

#ifdef highL
    function clik_lnlike_highL(cl,freq_params)
    real(dp) :: clik_lnlike_highL
    real(dp), intent(in) :: cl(lmax,num_cls_tot)
    real(dp), intent(in)  :: freq_params(1:num_freq_params-1)   
    real(dp) like_tot
    integer, parameter :: lmin=2
    real(dp)  cl_tt(tt_lmax)
    integer L, offset
    real(dp) A_ps_100, A_ps_143, A_ps_217, A_cib_143, A_cib_217, A_sz_143, r_ps, r_cib, &
      cal0, cal1, cal2, xi, A_ksz
    real(dp) a_ps_act_148,a_ps_act_217,a_ps_spt_95,a_ps_spt_150,a_ps_spt_220, &
       r_ps_spt_95x150,r_ps_spt_95x220,r_ps_150x220, &
       cal_acts_148,cal_acts_217,cal_acte_148,cal_acte_217,cal_spt_95,cal_spt_150,cal_spt_220
    
    !asz is already removed, start at second feq param
    A_ps_100=freq_params(1)
    A_ps_143 = freq_params(2)
    A_ps_217 = freq_params(3)
    A_cib_143=freq_params(4)
    A_cib_217=freq_params(5) 
    A_sz_143=freq_params(6)  !143
    r_ps = freq_params(7)
    r_cib=freq_params(8)
    cal0=freq_params(9)
    cal1=freq_params(10) 
    cal2=freq_params(11)
    xi=freq_params(12)
    A_ksz=freq_params(13)

     !skip ACT/SPT parameters 1,2,3,9,10,14 (already included in CAMspec nuisance params)
    offset = num_camSpec + num_plik
    if (.not. use_Camspec) then
     A_sz_143 = freq_params(offset+1)
     A_ksz= freq_params(offset+2)
     xi= freq_params(offset+3)
    end if    
    a_ps_act_148= freq_params(offset+4)
    a_ps_act_217= freq_params(offset+5)
    a_ps_spt_95= freq_params(offset+6)
    a_ps_spt_150= freq_params(offset+7)
    a_ps_spt_220= freq_params(offset+8)
    if (.not. use_Camspec) then
     A_cib_143= freq_params(offset+9)
     A_cib_217= freq_params(offset+10)
    end if
    r_ps_spt_95x150= freq_params(offset+11)
    r_ps_spt_95x220= freq_params(offset+12)
    r_ps_150x220= freq_params(offset+13)
    if (.not. use_Camspec) r_cib = freq_params(offset+14)
    cal_acts_148= freq_params(offset+15)
    cal_acts_217= freq_params(offset+16)
    cal_acte_148= freq_params(offset+17)
    cal_acte_217= freq_params(offset+18)
    cal_spt_95= freq_params(offset+19)
    cal_spt_150= freq_params(offset+20)
    cal_spt_220= freq_params(offset+21)

    do l =2, tt_lmax
     if (l.le.tt_lmax_mc) then
         cl_tt(l) = cl(l,1)*l*(l+1)/twopi
     else
         cl_tt(l) = 0.0d0
     endif
    enddo

    like_tot = 0.d0
    call highell_likelihood_compute(cl_tt,A_sz_143,A_ksz,xi,a_ps_act_148,a_ps_act_217,a_ps_spt_95,a_ps_spt_150,a_ps_spt_220, &
       A_cib_143,A_cib_217,  r_ps_spt_95x150,r_ps_spt_95x220,r_ps_150x220,r_cib, &
       cal_acts_148,cal_acts_217,cal_acte_148,cal_acte_217,cal_spt_95,cal_spt_150,cal_spt_220,like_tot)
    clik_lnlike_highL = like_tot

    if (Feedback>0) Print*,'highL lnlike = ',clik_lnlike_highL

    end function clik_lnlike_highL
#endif

   end module cliklike

