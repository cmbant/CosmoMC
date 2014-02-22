    !!! Generalized BAO module added by J. Dossett
    ! Copied structure from mpk.f90 and Reid BAO code
    !
    ! When using WiggleZ data set cite Blake et al. arXiv:1108.2635

    !for SDSS data set: http://www.sdss3.org/science/boss_publications.php

    ! For rescaling method see Hamann et al, http://arxiv.org/abs/1003.3999

    !AL/JH Oct 2012: encorporate DR9 data into something close to new cosmomc format
    !Dec 2013: merged in DR11 patch (Antonio J. Cuesta, for the BOSS collaboration)

    module bao
    use MatrixUtils
    use settings
    use cmbtypes
    use CosmoTheory
    use Calculator_Cosmology
    use Likelihood_Cosmology
    use IniObjects
    implicit none

    private

    type, extends(TCosmoCalcLikelihood) :: BAOLikelihood
        !    type, extends(TCosmoCalcLikelihood) :: BAOLikelihood
        integer :: num_bao ! total number of points used
        integer :: type_bao
        !what type of bao data is used
        !1: old sdss, no longer
        !2: A(z) =2 ie WiggleZ
        !3 D_V/rs in fitting forumla appox (DR9)
        !4: D_v - 6DF
        real(mcp), allocatable, dimension(:) :: bao_z, bao_obs, bao_err
        real(mcp), allocatable, dimension(:,:) :: bao_invcov
    contains
    procedure :: LogLike => BAO_LnLike
    procedure :: ReadIni => BAO_ReadIni
    procedure, private :: SDSS_dvtors
    procedure, private :: Acoustic
    procedure, private :: BAO_DR7_loglike
    procedure, private :: BAO_DR11_loglike
    end type BAOLikelihood

    integer,parameter :: DR11_alpha_npoints=280
    real(mcp), dimension (DR11_alpha_npoints) :: DR11_alpha_perp_file,DR11_alpha_plel_file
    real(mcp), dimension (DR11_alpha_npoints,DR11_alpha_npoints) ::   DR11_prob_file
    real(mcp) DR11_dalpha_perp, DR11_dalpha_plel
    real(mcp), dimension (10000) :: DR7_alpha_file, DR7_prob_file
    real(mcp) DR7_dalpha
    real rsdrag_theory
    real(mcp) :: BAO_fixed_rs = -1._mcp
    integer DR7_alpha_npoints

    public BAOLikelihood, BAOLikelihood_Add
    contains

    subroutine BAOLikelihood_Add(LikeList, Ini)
    class(TLikelihoodList) :: LikeList
    class(TSettingIni) :: ini
    Type(BAOLikelihood), pointer :: like
    integer numbaosets, i

    if (Ini%Read_Logical('use_BAO',.false.)) then
        numbaosets = Ini%Read_Int('bao_numdatasets',0)
        if (numbaosets<1) call MpiStop('Use_BAO but numbaosets = 0')
        if (Ini%Haskey('BAO_fixed_rs')) then
            BAO_fixed_rs= Ini%Read_Double('BAO_fixed_rs',-1._mcp)
        end if
        do i= 1, numbaosets
            allocate(like)
            call like%ReadDatasetFile(Ini%ReadFileName(numcat('bao_dataset',i)))
            like%LikelihoodType = 'BAO'
            like%needs_background_functions = .true.
            call LikeList%Add(like)
        end do
        if (Feedback>1) write(*,*) 'read bao datasets'
    end if

    end subroutine BAOLikelihood_Add

    subroutine BAO_ReadIni(like, Ini)
    class(BAOLikelihood) like
    class(TSettingIni) :: Ini
    character(LEN=:), allocatable :: bao_measurements_file, bao_invcov_file
    integer i,iopb, unit

    if (Feedback > 0) write (*,*) 'reading BAO data set: '//trim(like%name)
    like%num_bao = Ini%Read_Int('num_bao',0)
    if (like%num_bao.eq.0) write(*,*) ' ERROR: parameter num_bao not set'
    like%type_bao = Ini%Read_Int('type_bao',1)
    if(like%type_bao /= 3 .and. like%type_bao /=2 .and. like%type_bao /=4) then
        write(*,*) like%type_bao
        write(*,*)'ERROR: Invalid bao type specified in BAO dataset: '//trim(like%name)
        call MPIStop()
    end if

    allocate(like%bao_z(like%num_bao))
    allocate(like%bao_obs(like%num_bao))
    allocate(like%bao_err(like%num_bao))

    bao_measurements_file = Ini%ReadFileName('bao_measurements_file')
    unit = OpenNewTxtFile(bao_measurements_file)
    do i=1,like%num_bao
        read (unit,*, iostat=iopb) like%bao_z(i),like%bao_obs(i),like%bao_err(i)
    end do
    close(unit)

    if (like%name == 'DR7') then
        !don't used observed value, probabilty distribution instead
        call BAO_DR7_init(Ini%ReadFileName('prob_dist'))
    elseif (like%name == 'DR11CMASS') then
        !don't used observed value, probabilty distribution instead
        call BAO_DR11_init(Ini%ReadFileName('prob_dist'))
    else
        allocate(like%bao_invcov(like%num_bao,like%num_bao))
        like%bao_invcov=0

        if (Ini%HasKey('bao_invcov_file')) then
            bao_invcov_file  = Ini%ReadFileName('bao_invcov_file')
            unit= OpenNewTxtFile(bao_invcov_file)
            do i=1,like%num_bao
                read (unit,*, iostat=iopb) like%bao_invcov(i,:)
                if (iopb /= 0) then
                    call MpiStop('Error reading bao file '// bao_invcov_file)
                endif
            end do
            close(unit)
        else
            do i=1,like%num_bao
                !diagonal, or actually just 1..
                like%bao_invcov(i,i) = 1/like%bao_err(i)**2
            end do
        end if
    end if

    end subroutine BAO_ReadIni

    function Acoustic(like,CMB,z)
    class(BAOLikelihood) :: like
    class(CMBParams) CMB
    real(mcp) Acoustic
    real(mcp), intent(IN) :: z
    real(mcp) omh2,ckm,omegam,h

    omegam = 1.d0 - CMB%omv - CMB%omk
    h = CMB%h0/100
    ckm = const_c/1e3_mcp !JD c in km/s

    omh2 = omegam*h**2.d0
    Acoustic = 100*like%Calculator%BAO_D_v(z)*sqrt(omh2)/(ckm*z)
    end function Acoustic

    function SDSS_dvtors(like, CMB,z)
    !This uses numerical value of D_v/r_s, but re-scales it to match definition of SDSS
    !paper fitting at the fiducial model. Idea being it is also valid for e.g. varying N_eff
    class(BAOLikelihood) :: like
    class(CMBParams) CMB
    real(mcp) SDSS_dvtors
    real(mcp), intent(IN)::z
    real(mcp) rs
    real(mcp), parameter :: rs_rescale = 153.017d0/148.92 !149.0808

    !    rs = SDSS_CMBToBAOrs(CMB)
    rs = rsdrag_theory*rs_rescale !rescaled to match fitting formula for LCDM
    SDSS_dvtors = like%Calculator%BAO_D_v(z)/rs

    end function SDSS_dvtors


    !===================================================================================

    function BAO_LnLike(like, CMB, Theory, DataParams)
    Class(BAOLikelihood) :: like
    Class(CMBParams) CMB
    Class(TCosmoTheoryPredictions) Theory
    real(mcp) :: DataParams(:)
    integer j,k
    real(mcp) BAO_LnLike
    real(mcp), allocatable :: BAO_theory(:)

    if (BAO_fixed_rs>0) then
        !this is just for use for e.g. BAO 'only' constraints
        rsdrag_theory =  BAO_fixed_rs
    else
        rsdrag_theory =  Theory%derived_parameters( derived_rdrag )
    end if
    BAO_LnLike=0
    if (like%name=='DR7') then
        BAO_LnLike = like%BAO_DR7_loglike(CMB,like%bao_z(1))
    elseif (like%name=='DR11CMASS') then
        BAO_LnLike = like%BAO_DR11_loglike(CMB,like%bao_z(1))
    else
        allocate(BAO_theory(like%num_bao))

        if(like%type_bao ==3)then
            do j=1, like%num_bao
                BAO_theory(j) = like%SDSS_dvtors(CMB,like%bao_z(j))
            end do
        else if(like%type_bao ==2)then
            do j=1, like%num_bao
                BAO_theory(j) = like%Acoustic(CMB,like%bao_z(j))
            end do
        else if(like%type_bao ==4)then
            do j=1, like%num_bao
                BAO_theory(j) = like%Calculator%BAO_D_v(like%bao_z(j))
            end do
        end if

        do j=1, like%num_bao
            do k=1, like%num_bao
                BAO_LnLike = BAO_LnLike +&
                (BAO_theory(j)-like%bao_obs(j))*like%bao_invcov(j,k)*&
                (BAO_theory(k)-like%bao_obs(k))
            end do
        end do
        BAO_LnLike = BAO_LnLike/2.d0

        deallocate(BAO_theory)
    end if

    if(feedback>1) write(*,*) trim(like%name)//' BAO likelihood = ', BAO_LnLike

    end function BAO_LnLike


    subroutine BAO_DR7_init(fname)
    character(LEN=*), intent(in) :: fname
    real(mcp) :: tmp0,tmp1
    real(mcp) :: DR7_alpha =0
    integer ios,ii

    open(unit=7,file=fname,status='old')
    !Read data file
    ios = 0
    ii  = 0
    do while (ios.eq.0)
        read (7,*,iostat=ios) tmp0,tmp1
        if (ios .ne. 0) cycle
        if((ii.gt.1).and.(abs(DR7_dalpha-(tmp0-DR7_alpha)).gt.1e-6)) then
            stop 'binning should be uniform in sdss_baoDR7.txt'
        endif
        ii = ii+1
        DR7_alpha_file(ii) = tmp0
        DR7_prob_file (ii) = tmp1
        DR7_dalpha = tmp0-DR7_alpha
        DR7_alpha  = tmp0
    enddo
    DR7_alpha_npoints = ii
    if (ii.eq.0) call MpiStop('ERROR : reading file')
    close(7)
    !Normalize distribution (so that the peak value is 1.0)
    tmp0=0.0
    do ii=1,DR7_alpha_npoints
        if(DR7_prob_file(ii).gt.tmp0) then
            tmp0=DR7_prob_file(ii)
        endif
    enddo
    DR7_prob_file=DR7_prob_file/tmp0

    end subroutine BAO_DR7_init

    function BAO_DR7_loglike(like,CMB,z)
    Class(BAOLikelihood) :: like
    Class(CMBParams) CMB
    real (mcp) z, BAO_DR7_loglike, alpha_chain, prob
    real,parameter :: rs_wmap7=152.7934d0,dv1_wmap7=1340.177  !r_s and D_V computed for wmap7 cosmology
    integer ii
    alpha_chain = (like%SDSS_dvtors(CMB,z))/(dv1_wmap7/rs_wmap7)
    if ((alpha_chain.gt.DR7_alpha_file(DR7_alpha_npoints-1)).or.(alpha_chain.lt.DR7_alpha_file(1))) then
        BAO_DR7_loglike = logZero
    else
        ii=1+floor((alpha_chain-DR7_alpha_file(1))/DR7_dalpha)
        prob=DR7_prob_file(ii)+(DR7_prob_file(ii+1)-DR7_prob_file(ii))/ &
        (DR7_alpha_file(ii+1)-DR7_alpha_file(ii))*(alpha_chain-DR7_alpha_file(ii))
        BAO_DR7_loglike = -log( prob )
    endif

    end function BAO_DR7_loglike

    subroutine BAO_DR11_init(fname)
    character(LEN=*), intent(in) :: fname
    real(mcp) :: tmp0,tmp1,tmp2
    integer ios,ii,jj,nn

    open(unit=7,file=fname,status='old')
    ios = 0
    nn=0
    do while (ios.eq.0)
        read (7,*,iostat=ios) tmp0,tmp1,tmp2
        if (ios .ne. 0) cycle
        nn = nn + 1
        ii = 1 +     (nn-1)/DR11_alpha_npoints
        jj = 1 + mod((nn-1),DR11_alpha_npoints)
        DR11_alpha_perp_file(ii)   = tmp0
        DR11_alpha_plel_file(jj)   = tmp1
        DR11_prob_file(ii,jj)      = tmp2
    enddo
    close(7)
    DR11_dalpha_perp=DR11_alpha_perp_file(2)-DR11_alpha_perp_file(1)
    DR11_dalpha_plel=DR11_alpha_plel_file(2)-DR11_alpha_plel_file(1)
    !Normalize distribution (so that the peak value is 1.0)
    tmp0=0.0
    do ii=1,DR11_alpha_npoints
        do jj=1,DR11_alpha_npoints
            if(DR11_prob_file(ii,jj).gt.tmp0) then
                tmp0=DR11_prob_file(ii,jj)
            endif
        enddo
    enddo
    DR11_prob_file=DR11_prob_file/tmp0

    end subroutine BAO_DR11_init

    function BAO_DR11_loglike(like,CMB,z)
    Class(BAOLikelihood) :: like
    Class(CMBParams) CMB
    real (mcp) z, BAO_DR11_loglike, alpha_perp, alpha_plel, prob
    real,parameter :: rd_fid=149.28,H_fid=93.558,DA_fid=1359.72 !fiducial parameters
    integer ii,jj
    alpha_perp=(like%Calculator%AngularDiameterDistance(z)/rsdrag_theory)/(DA_fid/rd_fid)
    alpha_plel=(H_fid*rd_fid)/((const_c*like%Calculator%Hofz(z)/1.d3)*rsdrag_theory)
    if ((alpha_perp.lt.DR11_alpha_perp_file(1)).or.(alpha_perp.gt.DR11_alpha_perp_file(DR11_alpha_npoints-1)).or. &
    &   (alpha_plel.lt.DR11_alpha_plel_file(1)).or.(alpha_plel.gt.DR11_alpha_plel_file(DR11_alpha_npoints-1))) then
        BAO_DR11_loglike = logZero
    else
        ii=1+floor((alpha_perp-DR11_alpha_perp_file(1))/DR11_dalpha_perp)
        jj=1+floor((alpha_plel-DR11_alpha_plel_file(1))/DR11_dalpha_plel)
        prob=(1./((DR11_alpha_perp_file(ii+1)-DR11_alpha_perp_file(ii))*(DR11_alpha_plel_file(jj+1)-DR11_alpha_plel_file(jj))))*  &
        &       (DR11_prob_file(ii,jj)*(DR11_alpha_perp_file(ii+1)-alpha_perp)*(DR11_alpha_plel_file(jj+1)-alpha_plel) &
        &       -DR11_prob_file(ii+1,jj)*(DR11_alpha_perp_file(ii)-alpha_perp)*(DR11_alpha_plel_file(jj+1)-alpha_plel) &
        &       -DR11_prob_file(ii,jj+1)*(DR11_alpha_perp_file(ii+1)-alpha_perp)*(DR11_alpha_plel_file(jj)-alpha_plel) &
        &       +DR11_prob_file(ii+1,jj+1)*(DR11_alpha_perp_file(ii)-alpha_perp)*(DR11_alpha_plel_file(jj)-alpha_plel))
        if  (prob.gt.0.) then
            BAO_DR11_loglike = -log( prob )
        else
            BAO_DR11_loglike = logZero
        endif
    endif

    end function BAO_DR11_loglike

    function SDSS_CMBToBAOrs(CMB)
    Type(CMBParams) CMB
    real(mcp) ::  rsdrag
    real(mcp) :: SDSS_CMBToBAOrs
    real(mcp) :: zeq,zdrag,omh2,obh2,b1,b2
    real(mcp) :: rd,req,wkeq

    obh2=CMB%ombh2
    omh2=CMB%ombh2+CMB%omdmh2

    b1     = 0.313*omh2**(-0.419)*(1+0.607*omh2**0.674)
    b2     = 0.238*omh2**0.223
    zdrag  = 1291.*omh2**0.251*(1.+b1*obh2**b2)/(1.+0.659*omh2**0.828)
    zeq    = 2.50e4*omh2*(2.726/2.7)**(-4.)
    wkeq   = 7.46e-2*omh2*(2.726/2.7)**(-2)
    req    = 31.5*obh2*(2.726/2.7)**(-4)*(1e3/zeq)
    rd     = 31.5*obh2*(2.726/2.7)**(-4)*(1e3/zdrag)
    rsdrag = 2./(3.*wkeq)*sqrt(6./req)*log((sqrt(1.+rd)+sqrt(rd+req))/(1.+sqrt(req)))

    SDSS_CMBToBAOrs = rsdrag

    end function SDSS_CMBToBAOrs


    end module bao