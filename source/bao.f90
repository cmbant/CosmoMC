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
    use CosmologyTypes
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
        !5: D_A/rs (DR8)
        real(mcp), allocatable, dimension(:) :: bao_z, bao_obs, bao_err
        real(mcp), allocatable, dimension(:,:) :: bao_invcov
        ! 5: (D_V/rs, F_AP, f sigma_8) 
        real(mcp), allocatable, dimension(:) :: bao_AP,bao_fsigma8

    contains
    procedure :: LogLike => BAO_LnLike
    procedure :: ReadIni => BAO_ReadIni
    procedure, private :: SDSS_dvtors
    procedure, private :: SDSS_dAtors
    procedure, private :: Acoustic
    procedure, private :: BAO_DR7_loglike
    procedure, private :: BAO_DR11_loglike
    procedure, private :: BAO_RSD_loglike
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
    logical :: use_bao_lss  = .false.

    public BAOLikelihood, BAOLikelihood_Add, use_bao_lss
    contains

    subroutine BAOLikelihood_Add(LikeList, Ini)
    class(TLikelihoodList) :: LikeList
    class(TSettingIni) :: ini
    Type(BAOLikelihood), pointer :: this
    integer numbaosets, i

    if (Ini%Read_Logical('use_BAO',.false.)) then
        numbaosets = Ini%Read_Int('bao_numdatasets',0)
        if (numbaosets<1) call MpiStop('Use_BAO but numbaosets = 0')
        if (Ini%Haskey('BAO_fixed_rs')) then
            BAO_fixed_rs= Ini%Read_Double('BAO_fixed_rs',-1._mcp)
        end if
        do i= 1, numbaosets
            allocate(this)
            call this%ReadDatasetFile(Ini%ReadFileName(numcat('bao_dataset',i)))
            this%LikelihoodType = 'BAO'
            this%num_z = Ini%Read_Int('nz_bao',0)
            if (this%num_z==0) then
                this%needs_background_functions = .true.
            else
                this%needs_powerspectra = .true.
                this%max_z = Ini%Read_Double('max_z_bao',1._mcp)
                use_bao_lss = .true.
            end if
            call LikeList%Add(this)
        end do
        if (Feedback>1) write(*,*) 'read BAO data sets'
    end if

    end subroutine BAOLikelihood_Add

    subroutine BAO_ReadIni(this, Ini)
    class(BAOLikelihood) this
    class(TSettingIni) :: Ini
    character(LEN=:), allocatable :: bao_measurements_file, bao_invcov_file
    integer i,iopb
    Type(TTextFile) :: F

    if (Feedback > 0) write (*,*) 'reading BAO data set: '//trim(this%name)
    this%num_bao = Ini%Read_Int('num_bao',0)
    if (this%num_bao.eq.0) write(*,*) ' ERROR: parameter num_bao not set'
    this%type_bao = Ini%Read_Int('type_bao',1)
    if(this%type_bao /= 3 .and. this%type_bao /=2 .and. this%type_bao /=4 .and. this%type_bao /=5) then
        write(*,*) this%type_bao
        write(*,*)'ERROR: Invalid bao type specified in BAO dataset: '//trim(this%name)
        call MPIStop()
    end if

    allocate(this%bao_z(this%num_bao))
    allocate(this%bao_obs(this%num_bao))
    allocate(this%bao_err(this%num_bao))

    bao_measurements_file = Ini%ReadFileName('bao_measurements_file')
    call F%Open(bao_measurements_file)
    if(this%type_bao == 5) then
       if (this%num_bao /= 1) then
          write(*,*)'ERROR: Need 1 bao data point for: '//trim(this%name)
          call MPIStop()
       end if
       allocate(this%bao_AP(this%num_bao))
       allocate(this%bao_fsigma8(this%num_bao))
       read (F%unit,*, iostat=iopb) this%bao_z(1),this%bao_obs(1),this%bao_AP(1),this%bao_fsigma8(1)
    else
        do i=1,this%num_bao
            read (F%unit,*, iostat=iopb) this%bao_z(i),this%bao_obs(i),this%bao_err(i)
        end do
    end if
    call F%Close()

    if (this%name == 'DR7') then
        !don't used observed value, probabilty distribution instead
        call BAO_DR7_init(Ini%ReadFileName('prob_dist'))
    elseif (this%name == 'DR11CMASS') then
        !don't used observed value, probabilty distribution instead
        call BAO_DR11_init(Ini%ReadFileName('prob_dist'))
    elseif (this%type_bao == 5) then
       allocate(this%bao_invcov(3,3))
       this%bao_invcov=0
       if (Ini%HasKey('bao_invcov_file')) then
          bao_invcov_file  = Ini%ReadFileName('bao_invcov_file')
          call File%ReadTextMatrix(bao_invcov_file, this%bao_invcov)
       else
          write(*,*)'ERROR: No inverse covariance matrix for: '//trim(this%name)
          call MPIStop()
       end if
    else
        allocate(this%bao_invcov(this%num_bao,this%num_bao))
        this%bao_invcov=0

        if (Ini%HasKey('bao_invcov_file')) then
            bao_invcov_file  = Ini%ReadFileName('bao_invcov_file')
            call File%ReadTextMatrix(bao_invcov_file, this%bao_invcov)
        else
            do i=1,this%num_bao
                !diagonal, or actually just 1..
                this%bao_invcov(i,i) = 1/this%bao_err(i)**2
            end do
        end if
    end if

    end subroutine BAO_ReadIni

    function Acoustic(this,CMB,z)
    class(BAOLikelihood) :: this
    class(CMBParams) CMB
    real(mcp) Acoustic
    real(mcp), intent(IN) :: z
    real(mcp) omh2,ckm,omegam,h

    omegam = 1.d0 - CMB%omv - CMB%omk
    h = CMB%h0/100
    ckm = const_c/1e3_mcp !JD c in km/s

    omh2 = omegam*h**2.d0
    Acoustic = 100*this%Calculator%BAO_D_v(z)*sqrt(omh2)/(ckm*z)
    end function Acoustic

    function SDSS_dvtors(this, CMB,z)
    !This uses numerical value of D_v/r_s, but re-scales it to match definition of SDSS
    !paper fitting at the fiducial model. Idea being it is also valid for e.g. varying N_eff
    class(BAOLikelihood) :: this
    class(CMBParams) CMB
    real(mcp) SDSS_dvtors
    real(mcp), intent(IN)::z
    real(mcp) rs
    real(mcp), parameter :: rs_rescale = 153.017d0/148.92 !149.0808

    !    rs = SDSS_CMBToBAOrs(CMB)
    rs = rsdrag_theory*rs_rescale !rescaled to match fitting formula for LCDM
    SDSS_dvtors = this%Calculator%BAO_D_v(z)/rs

    end function SDSS_dvtors

   ! HS modified SDSS_dvtors to calculate D_A/rs 
    function SDSS_dAtors(this, CMB,z)
    !This uses numerical value of D_A/r_s, but re-scales it to match definition of SDSS
    !paper fitting at the fiducial model. Idea being it is also valid for e.g. varying N_eff
    class(BAOLikelihood) :: this
    class(CMBParams) CMB
    real(mcp) SDSS_dAtors
    real(mcp), intent(IN)::z
    real(mcp) rs
    real(mcp), parameter :: rs_rescale = 153.017d0/148.92 !149.0808

    !    rs = SDSS_CMBToBAOrs(CMB)
    rs = rsdrag_theory*rs_rescale !rescaled to match fitting formula for LCDM
    SDSS_dAtors = this%Calculator%AngularDiameterDistance(z)/rs
    end function SDSS_dAtors



    !===================================================================================

    function BAO_LnLike(this, CMB, Theory, DataParams)
    Class(BAOLikelihood) :: this
    Class(CMBParams) CMB
    Class(TCosmoTheoryPredictions), target :: Theory
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
    if (this%name=='DR7') then
        BAO_LnLike = this%BAO_DR7_loglike(CMB,this%bao_z(1))
    elseif (this%name=='DR11CMASS') then
        BAO_LnLike = this%BAO_DR11_loglike(CMB,this%bao_z(1))
    elseif (this%type_bao==5) then
        BAO_LnLike = this%BAO_RSD_loglike(CMB, Theory)
    else
        allocate(BAO_theory(this%num_bao))

        if(this%type_bao ==3)then
            do j=1, this%num_bao
                BAO_theory(j) = this%SDSS_dvtors(CMB,this%bao_z(j))
            end do
        else if(this%type_bao ==2)then
            do j=1, this%num_bao
                BAO_theory(j) = this%Acoustic(CMB,this%bao_z(j))
            end do
        else if(this%type_bao ==4)then
            do j=1, this%num_bao
                BAO_theory(j) = this%Calculator%BAO_D_v(this%bao_z(j))
            end do
        else if(this%type_bao ==5)then
            do j=1, this%num_bao
                BAO_theory(j) = this%SDSS_dAtors(CMB,this%bao_z(j))
            end do

        end if

        do j=1, this%num_bao
            do k=1, this%num_bao
                BAO_LnLike = BAO_LnLike +&
                (BAO_theory(j)-this%bao_obs(j))*this%bao_invcov(j,k)*&
                (BAO_theory(k)-this%bao_obs(k))
            end do
        end do
        BAO_LnLike = BAO_LnLike/2.d0

        deallocate(BAO_theory)
    end if

    if(feedback>1) write(*,*) trim(this%name)//' BAO likelihood = ', BAO_LnLike

    end function BAO_LnLike

    ! RSD like for (D_V/r_s, F_AP, f sigma_8)
    function BAO_RSD_loglike(this,CMB,Theory)
      Class(BAOLikelihood) :: this
      Class(CMBParams) CMB
      Class(TCosmoTheoryPredictions), target :: Theory
      real(mcp) BAO_RSD_loglike
      real(mcp) :: f, z,vec(3) 
      z = this%bao_z(1)
      f = -(1+z)/Theory%sigma_8_z%Value(z)*Theory%sigma_8_z%Derivative(z)
      vec(1) = this%Calculator%BAO_D_v(z)/rsdrag_theory - this%bao_obs(1) ! D_V/R_s
      vec(2) = (1+z)*this%Calculator%AngularDiameterDistance(z)*this%Calculator%Hofz(z) -this%bao_AP(1)! F_AP
      vec(3) = f*Theory%sigma_8_z%Value(z) -this%bao_fsigma8(1) ! f sigma_8
      BAO_RSD_loglike = BAO_RSD_loglike + &
           (vec(1)*(vec(1)*this%bao_invcov(1,1)+vec(2)*this%bao_invcov(1,2)+vec(3)*this%bao_invcov(1,3)) + &
            vec(2)*(vec(1)*this%bao_invcov(2,1)+vec(2)*this%bao_invcov(2,2)+vec(3)*this%bao_invcov(2,3)) + &
            vec(3)*(vec(1)*this%bao_invcov(3,1)+vec(2)*this%bao_invcov(3,2)+vec(3)*this%bao_invcov(3,3)))/2.0
    end function BAO_RSD_loglike

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

    function BAO_DR7_loglike(this,CMB,z)
    Class(BAOLikelihood) :: this
    Class(CMBParams) CMB
    real (mcp) z, BAO_DR7_loglike, alpha_chain, prob
    real,parameter :: rs_wmap7=152.7934d0,dv1_wmap7=1340.177  !r_s and D_V computed for wmap7 cosmology
    integer ii
    alpha_chain = (this%SDSS_dvtors(CMB,z))/(dv1_wmap7/rs_wmap7)
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

    function BAO_DR11_loglike(this,CMB,z)
    Class(BAOLikelihood) :: this
    Class(CMBParams) CMB
    real (mcp) z, BAO_DR11_loglike, alpha_perp, alpha_plel, prob
    real,parameter :: rd_fid=149.28,H_fid=93.558,DA_fid=1359.72 !fiducial parameters
    integer ii,jj
    alpha_perp=(this%Calculator%AngularDiameterDistance(z)/rsdrag_theory)/(DA_fid/rd_fid)
    alpha_plel=(H_fid*rd_fid)/((const_c*this%Calculator%Hofz(z)/1.d3)*rsdrag_theory)
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
