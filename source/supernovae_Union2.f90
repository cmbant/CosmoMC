    ! UNION 2.1 Supernovae Ia dataset
    !
    ! This module uses the SCP (Supernova Cosmology Project) Union 2
    ! compilation. Please cite
    ! "Suzuki et al. (SCP) 2011, arXiv:1105.3470 (2011 ApJ Dec 20 issue)".
    ! and the references of other compiled supernovae data are in there.
    !
    ! Originally by A Slosar, heavily based on the original code by A Lewis, S Bridle
    ! and D Rapetti E-mail: Anze Slosar (anze@berkeley.edu) for questions
    ! about the code and David Rubin (rubind@berkeley.edu) for questions
    ! regarding the dataset itself.   Updated by Nao Suzuki (nsuzuki@lbl.gov)
    !
    ! Marginalizes anayltically over H_0 with flat prior.  (equivalent to
    ! marginalizing over M, absolute magnitude; see appendix F of cosmomc
    ! paper). Resultant log likelihood has arbitary origin and is
    ! numerically equal to -chi^2/2 value at the best-fit value.
    !
    ! Update Note :
    !
    ! Union1   (Kowalski et al 2008)  : 307 SNe with SALT1 fit (Guy et al 2005)
    ! Union2   (Amanullah et al 2010) : 557 SNe with SALT2 fit (Guy et al 2007)
    ! Union2.1 (Suzuki et al 2011)    : 580 SNe with SALT2 fit (Guy et al 2007)
    !
    ! The following parameters are used to calculate distance moduli
    ! (see Suzuki et al. 2011 for complete description)
    !
    ! alpha 0.121851859725    ! Stretch Correction Factor
    ! beta 2.46569277393      ! Color   Correction Factor
    ! delta -0.0363405630486
    ! M(h=0.7, statistical only) -19.3182761161 ! Absolute B Magnitue of SNIa
    ! M(h=0.7, with systematics) -19.3081547178 ! Absolute B Magnitue of SNIa
    !
    !  Tips for running cosmomc with SCP UNION2.1 data
    !
    !  1) Place the following 3 data files in your cosmomc data dir (DataDir)
    !     a) sn_z_mu_dmu_plow_union2.1.txt : SN data
    !        SN name, z, distance moduli mu, mu error, host mass weight probability
    !     b) sn_covmat_sys_union2.1.txt    : Covariance Matrix with    systematic error
    !     c) sn_covmat_nosys_union2.1.txt  : Covariance Matrix without systematic error
    !
    !  2) Make sure DataDir is set in your settings.f90
    !     character(LEN=1024) :: DataDir='yourdirpathto/cosmomc/data/'
    !     The default is 'data/' and if it works for you, just leave it as it is
    !
    !  3) Pick SN data 'with' or 'without' systematic error
    !     (default is 'with' systematic error)
    !     Modify the folowing SN_syscovamat=.True. or .False.
    !
    !  4) To make UNION2 as your default,
    !     either rename supernovae_union2.1.f90 as supernovae.f90 and recompile it
    !     or change targets in your Makefile from supernova to supernovae_union2.1
    !
    !  Note: In your default params.ini, there is a line for 'SN_filename', but this
    !     union2 module does not use it.  You can leave it as it is, and cosmomc
    !     runs without any error but that information is not used.
    !     To avoid confusion, you may want to comment it out.
    !
    !   Update Note by Nao Suzuki (LBNL)

    module Union2
    use CosmologyTypes
    use MatrixUtils
    use likelihood
    use CosmoTheory
    use Calculator_Cosmology
    use Likelihood_Cosmology
    implicit none
    private

    integer, parameter :: SN_num = 580

    type, extends(TCosmoCalcLikelihood) :: Union2Likelihood
        double precision :: SN_z(SN_num), SN_moduli(SN_num), SN_modulierr(SN_num), SN_plow(SN_num)
        double precision :: SN_Ninv(SN_num,SN_Num)
    contains
    procedure :: logLikeTheory => SN_LnLike
    end type Union2Likelihood

    public Union2Likelihood,Union2Likelihood_Add
    contains

    subroutine Union2Likelihood_Add(LikeList, Ini)
    class(TLikelihoodList) :: LikeList
    class(TSettingIni) :: ini
    Type(Union2Likelihood), pointer :: this
    character (LEN=20):: name
    integer i
    ! The following line selects which error estimate to use
    ! default .True. = with systematic errors
    logical :: Union_syscovmat = .False.  !! Use covariance matrix with or without systematics
    Type(TTextFile) :: F

    if (.not. Ini%Read_Logical('use_Union',.false.)) return

    allocate(this)
    this%LikelihoodType = 'SN'
    this%name='Union2.1'
    this%needs_background_functions = .true.
    call LikeList%Add(this)

    Union_syscovmat = Ini%read_Logical('Union_syscovmat',Union_syscovmat)

    if (Feedback > 0) write (*,*) 'Reading: supernovae data'
    call F%Open(trim(DataDir)//'sn_z_mu_dmu_plow_union2.1.txt')
    do i=1,  sn_num
        read(F%unit, *) name, this%SN_z(i),this%SN_moduli(i), &
        this%SN_modulierr(i),this%SN_plow(i)
        !     read(unit, *) name, SN_z(i),SN_moduli(i)
    end do
    call F%Close()

    if (Union_syscovmat) then
        call F%Open(trim(DataDir)//'sn_wmat_sys_union2.1.txt')
    else
        call F%Open(trim(DataDir)//'sn_wmat_nosys_union2.1.txt')
    end if

    do i=1, sn_num
        read (F%unit,*) this%sn_ninv (i,1:sn_num)
    end do

    call F%Close()

    end subroutine Union2Likelihood_Add

    function SN_LnLike(this, CMB)
    !Assume this is called just after CAMB with the correct model  use camb
    Class(CMBParams) CMB
    Class(Union2Likelihood) :: this
    real(mcp) SN_LnLike
    integer i
    double precision z
    real(mcp) diffs(SN_num), chisq

    !! This is actually seems to be faster without OMP

    do i=1, SN_num
        z= this%SN_z(i)
        diffs(i) = 5*log10((1+z)**2*this%Calculator%AngularDiameterDistance(z))+25 -this%sn_moduli(i)
    end do

    chisq = dot_product(diffs,matmul(this%sn_ninv,diffs))

    !! H0 normalisation alla Bridle and co.

    if (Feedback > 1) write (*,*) 'SN chisq: ', chisq

    SN_LnLike = chisq/2


    end function SN_LnLike


    end module Union2
