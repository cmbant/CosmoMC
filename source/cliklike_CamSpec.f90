    module cliklike
    use cmbtypes
    use settings
    use temp_like
    implicit none

    logical :: use_clik = .false.  
    logical :: using_CAMspec = .true.
    logical :: using_PLik = .false.
    integer, parameter :: dp = kind(1.d0)


    private
    public :: clik_readParams,clik_lnlike, use_clik, using_CAMspec, using_PLik

    contains

    subroutine clik_readParams(Ini)
    use settings
    Type(TIniFile) Ini
    character (LEN=Ini_max_string_len) :: likefilename,sz143filename,&
    beamfilename, kszfilename,tszxcibfilename

    print *,' using cliklike_CamSpec'
    likefilename=ReadIniFileName(Ini,'likefile',NotFoundFail = .true.)

    sz143filename=ReadIniFileName(Ini,'sz143file',NotFoundFail = .true.)

    tszxcibfilename=ReadIniFileName(Ini,'tszxcibfile',NotFoundFail = .true.)

    kszfilename=ReadIniFileName(Ini,'kszfile',NotFoundFail = .true.)

    beamfilename=ReadIniFileName(Ini,'beamfile',NotFoundFail = .true.)

    call like_init(likefilename,sz143filename,tszxcibfilename,kszfilename,beamfilename)

    end subroutine clik_readParams

    function clik_lnlike(cl,freq_params)
    real(dp) :: clik_lnlike
    real(dp), intent(in) :: cl(lmax,num_cls_tot)
    real(dp), intent(in), optional :: freq_params(1:num_freq_params-1)   
    real(dp) zlike, A_ps_100, A_ps_143, A_ps_217, A_cib_143, A_cib_217, A_sz, r_ps, r_cib, &
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
    A_sz=freq_params(6)  !143
    r_ps = freq_params(7)
    r_cib=freq_params(8)
    cal0=freq_params(9)
    cal1=freq_params(10) 
    cal2=freq_params(11)
    xi=freq_params(12)
    A_ksz=freq_params(13)
    beamcoeffs=freq_params(14:14+nbeammodes-1)

    do L=lmin,lmax
        cell_cmb(L)=cl(L,1)/twopi !this is a georgeism
    enddo

    do ii=1,nbeams
        do jj=1,nmodesperbeam
            beam_coeffs(ii,jj)=beamcoeffs(jj+nmodesperbeam*(ii-1))
        enddo
    enddo

    call calc_like(zlike,  cell_cmb, A_ps_100,  A_ps_143, A_ps_217, A_cib_143, A_cib_217, A_sz,  &
    r_ps, r_cib, xi, A_ksz, cal0, cal1, cal2, beam_coeffs)

    clik_lnlike = zlike/2

    if (Feedback>0) Print*,'CamSpec lnlike = ',clik_lnlike

    end function clik_lnlike


    end module cliklike

