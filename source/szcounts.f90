    !-----------------------------------------------------------!
    ! Module for SZ cluster counts likelihood
    !
    ! Module prepared for cosmomcplanck by A. Bonaldi 2015
    ! This module have been used to produce the results in
    ! Planck 2015 results XXIV: Cosmology from Sunyaev-Zeldovich cluster counts
    ! for questions and problems contact anna.bonaldi@manchester.ac.uk
    !
    !   Original code to compute cluster counts written from 2001
    !
    !                   v 0.1
    !           Jochen Weller and Richard Battye
    !                   4/2/2007
    ! all the necessary modules and subroutines
    !
    !                   v 1.0
    !               Jochen Weller
    !                 10/1/2008
    !              modified for optical cluster counts
    !
    !                   v 2.0
    !               Jochen Weller
    !                 30/12/2009
    !              modified for SZ cluster counts
    !
    !                   v 3.0
    !                Anna Bonaldi 2011-2013
    !      - Include realistic selection function for Planck
    !      - missing redshifts, errors on redshifts
    !               Matthieu Roman 2012-2013
    !      - QA completeness
    !
    !                   v 4.0
    !                Anna Bonaldi 2014-2015
    !      - updates for 2015 analysis
    !      - dN and dzdq likelihood
    !-----------------------------------------------------------!


    module constants_sz
    use PRECISION
    implicit none
    real(dl), PARAMETER :: Mpc = 3.08568025d22 ! Mpc in metres
    real(dl), PARAMETER :: G = 6.67300d-11 ! Newton's constant in m^3kg^-1s^-2
    real(dl), PARAMETER :: c = 3.0d8 ! speed of light in m/s
    real(dl), PARAMETER :: msun = 1.98892d30 ! mass of sun in kg
    real(dl), PARAMETER :: rhocrit0=2.7751973751261264d11 ! rhocrit in units of h^-1 Msun/ h^-3 Mpc^3
    real(dl),PARAMETER :: PI=3.141592653589793238463D0
    end module constants_sz


    ! ####################################################################


    module cosmology
    ! module computing cosmological functions
    USE precision
    USE CONSTANTS_sz
    use Interpolation
    implicit none
    ! note that we switch to the generalized Linder parametrization now
    public
    TYPE cospar
        REAL(dl) :: H0,w0,w1,omegam,omegav,n,sig8,omegak,omegabh2,gamma,ystar,alpha,sigmaM,bias,biasinv,logystar,beta
        class(TCubicSpline), pointer :: sigmaR
    END TYPE cospar
    Type (cospar), SAVE :: cosmopar

    contains


    !-----------------------------------------------------------------------

    function Eh(z)
    ! E=H(z)/H0
    real(dl) :: Eh
    real(dl), intent(in) :: z
    real(dl) ::  a
    a = 1.0/(1.0+z)
    Eh=sqrt(cosmopar%omegam*a**(-3)+cosmopar%omegav*a**(-3.0&
        &*(1.0+cosmopar%w0+cosmopar%w1))*exp(3.0*cosmopar%w1*(a-1.0))+cosmopar%omegak*a**(-2))
    RETURN
    end function Eh


    !-----------------------------------------------------------------------

    function Omegam(z)

    real(dl), intent(in) :: z
    real(dl) :: Omegam

    Omegam=cosmopar%omegam*(1.0+z)**3/(Eh(z))**2

    RETURN
    end function Omegam



    !-----------------------------------------------------------------------

    function Omegade(z)

    real(dl), intent(in) :: z
    real(dl) :: Omegade

    Omegade=cosmopar%omegav*(1.0+z)**(3.0*(1.0+cosmopar%w0&
        &+cosmopar%w1))*exp(-3.0*cosmopar%w1*z/(1.0+z))/(Eh(z))**2

    RETURN
    end function Omegade


    !-----------------------------------------------------------------------

    function Omegak(z)

    real(dl), intent(in) :: z
    real(dl) :: Omegak

    Omegak=1.0_dl-Omegam(z)-Omegade(z)

    RETURN
    end function Omegak

    !-----------------------------------------------------------------------

    function r(z)

    ! coordinate distance , in units of h^-1 Mpc


    real(dl), INTENT(IN) :: z
    real(dl) :: r
    real(dl), PARAMETER :: tol=1.0d-5
    real(dl) :: integral,rombint
    external rombint
    integral=rombint(rint,0._dl,z,tol)
    ! maybe soften this to 10^-5
    if (cosmopar%omegak == 0._dl) then
        r= c*integral/1.0d5
    elseif (cosmopar%omegak > 0._dl) then
        r=c/sqrt(cosmopar%omegak)*dsinh(sqrt(cosmopar%omegak)*integral)/1.0d5
    else
        r=c/sqrt(-cosmopar%omegak)*dsin(sqrt(-cosmopar%omegak)*integral)/1.0d5
    end if


    return
    end function r

    !-----------------------------------------------------------------------
    function rint(z)

    ! integrand for coordinate distance

    real(dl), INTENT(IN) :: z
    real(dl) :: rint
    rint = 1._dl/Eh(z)
    RETURN
    end function rint


    !-----------------------------------------------------------------------

    function dA(z)
    ! angular diameter distance in units of h^-1 Mpc

    real(dl), INTENT(IN) :: z
    real(dl) :: dA
    dA = r(z)/(1._dl+z)
    RETURN
    end function dA
    !-----------------------------------------------------------------------

    function dVdzdO(z)
    ! volume element in units of h^-3 Mpc^3
    real(dl), INTENT(IN) :: z
    real(dl) :: dVdzdO

    dVdzdO = c/1.0d5*r(z)**2/Eh(z)
    end function dVdzdO
    !-----------------------------------------------------------------------


    end module cosmology

    ! ####################################################################


    module survey
    use PRECISION
    implicit none
    TYPE survpar
        ! Note that in the SZ version mlimin is only used for the number of mass bins
        REAL(dl) :: ylimin,deg2,ymaxin
        INTEGER :: ybin,Nscat
        REAL(dl) :: sfid,ab,nb
        REAL(dl), allocatable :: b(:)
        !     REAL(dl), allocatable :: mbias(:),mscatter(:)
    END TYPE survpar
    TYPE (survpar), SAVE :: surveypar

    end module survey

    ! ####################################################################



    module massobservable
    ! scaling relations
    USE PRECISION
    USE SURVEY
    USE COSMOLOGY
    implicit none
    real(dl),parameter :: thetastar=6.997 ! dic 2012 new scalings
    real(dl),parameter :: alpha_theta=1./3.
    REAL(dl),parameter::q=6. ! signal-to-noise threshold of the catalogue
    contains

    function theta500(m,z)
    !size-mass/redshift relation
    real(dl), INTENT(IN) :: z,m
    real(dl) :: theta500,thetastar2,m2
    m2=m*cosmopar%bias !hydro bias
    thetastar2=thetastar*(cosmopar%H0/70.)**(-2./3.)
    theta500=thetastar2*(m2/3.e14*(100/cosmopar%H0))**alpha_theta*Eh(z)**(-2./3)*(100.0*da(z)/500.0/cosmopar%H0)**(-1.)
    RETURN
    end function theta500



    function y500(m,z)
    !y-mass/redshift relation
    real(dl), INTENT(IN) :: z,m
    real(dl) :: y500,ystar2,m2,alpha
    m2=m*cosmopar%bias !hydro bias
    alpha=cosmopar%alpha
    ystar2=cosmopar%ystar
    ystar2=ystar2*(cosmopar%H0/70.)**(-2.+alpha)
    y500=ystar2*(m2/3.0e14*(100./cosmopar%H0))**alpha*Eh(z)**(cosmopar%beta)*(100.0*da(z)/500.0/cosmopar%H0)**(-2.)
    RETURN
    end function y500

    function erf_compl(y,sn,q)
    !completeness with error function
    REAL(dl)::y,sn,arg,erf_compl,q
    arg=(y-q*sn)/sqrt(2.)/sn
    erf_compl=(erf(arg)+1.)/2.
    end function erf_compl


    end module massobservable

    ! ####################################################################

    module power_sz
    use PRECISION
    use cosmology

    use Calculator_Cosmology
    implicit none
    public
    REAL(dl) :: normsig8
    REAL(dl) :: normgrowth
    contains


    !-----------------------------------------------------------------------
    subroutine INIGROWTH
    ! normalize growth factor to 1 today
    normgrowth = 1.0_dl ! need to define this first because it is used in delta(z)
    normgrowth = 1.0/delta(0.0_dl)

    end subroutine INIGROWTH

    !-----------------------------------------------------------------------

    function delta(z)
    ! growth factor
    real(dl), INTENT(IN) :: z
    real(dl) :: delta
    real(dl), PARAMETER:: zmax=1000.0_dl ! infinity for growth integration
    real(dl), PARAMETER :: tol=1.0d-8
    INTEGER, PARAMETER :: n=2
    INTEGER, PARAMETER :: nw = n
    real(dl) :: c(32),w(nw,9)
    real(dl) :: y(n)
    real(dl) :: a,aend
    real :: dummyr
    real(dl) :: integ,gdum,rombint
    INTEGER :: ind
    external dverk,rombint

    y(1) = 1.0_dl
    a = 1.0_dl/(1.0_dl+zmax)
    y(2) = 0.0_dl
    aend=1.0_dl/(1.0_dl+z)
    ind = 1
    ! dverk is defined in CAMB subroutines.f90
    ! dummyr is not assigned any value

    if (cosmopar%gamma.eq.-1.0) then
        call dverk(dummyr, n, ddelta, a, y, aend, tol,ind,c,nw,w)
        if (ind.ne.3) write(*,*) 'Problem in dverk',ind
        delta = normgrowth*y(2)
        !       write(*,*) a,normgrowth*(y(1)/a-y(2)/a/a)
    else
        integ=rombint(growint,a,aend,tol)
        gdum=aend*exp(integ)
        delta=gdum*normgrowth
    end if
    RETURN
    end function delta


    subroutine ddelta(dummyr,n,a,y,yprime)
    ! note dummyr is completely irrelevant but required for CAMB dverk
    real(dl), INTENT(IN) :: a
    real, INTENT(IN) :: dummyr
    INTEGER, INTENT(IN) :: n
    real(dl) :: y(n),yprime(n)
    real(dl) :: z,w,u,del

    z=1.0_dl/a-1.0_dl
    w=cosmopar%w0+cosmopar%w1*(1.0_dl-a)
    u = y(1)
    del = y(2)

    yprime(1) = -1.5_dl*(1.0_dl+omegak(z)/3.0_dl-w*omegade(z))*u/a+1.5_dl*omegam(z)*del/a/a
    yprime(2) = u
    RETURN
    end subroutine ddelta


    function growint(a)
    real(dl), INTENT(IN) :: a
    real(dl) :: growint,z
    z=1.0_dl/a-1.0_dl
    growint=(omegam(z)**cosmopar%gamma-1.0_dl)/a

    return
    end function growint



    end module power_sz

    ! ####################################################################

    module massfunction
    USE PRECISION
    USE CONSTANTS_sz
    USE POWER_SZ
    USE massobservable

    implicit none
    TYPE MASSPAR
        REAL(dl) :: Amf,Bmf,epsmf,dso,pmf
        INTEGER :: psind
    END TYPE MASSPAR
    TYPE (masspar) :: massfnpar
    contains


    function dndlnM(z,M,g)
    !mass function: Jenkins, Tinker (default), Watson
    ! mass in units of h^-1 M_sun
    real(dl), INTENT(in) :: z,M,g
    real(dl) :: dndlnM
    real(dl) :: R,rhom0,rhom
    real(dl) :: dMdR,sR,fJen
    real(dl) :: fTink
    real(dl) :: fWat,p,q,CDelta,Gamma,ddz,z0
    real(dl) :: alpha,dsoz
    real(dl) :: del(1:9),par_aa(1:9),par_a(1:9),par_b(1:9),par_c(1:9)
    real(dl) :: der_aa(1:9),der_a(1:9),der_b(1:9),der_c(1:9)
    real(dl) :: par1,par2,par3,par4
    integer :: total,i

    rhom0=cosmopar%omegam*rhocrit0
    rhom=rhom0*(1.0+z)**3

    total=9
    del(1)=200
    del(2)=300
    del(3)=400
    del(4)=600
    del(5)=800
    del(6)=1200
    del(7)=1600
    del(8)=2400
    del(9)=3200
    par_aa(1)=0.186
    par_aa(2)=0.200
    par_aa(3)=0.212
    par_aa(4)=0.218
    par_aa(5)=0.248
    par_aa(6)=0.255
    par_aa(7)=0.260
    par_aa(8)=0.260
    par_aa(9)=0.260
    par_a(1)=1.47
    par_a(2)=1.52
    par_a(3)=1.56
    par_a(4)=1.61
    par_a(5)=1.87
    par_a(6)=2.13
    par_a(7)=2.30
    par_a(8)=2.53
    par_a(9)=2.66
    par_b(1)=2.57
    par_b(2)=2.25
    par_b(3)=2.05
    par_b(4)=1.87
    par_b(5)=1.59
    par_b(6)=1.51
    par_b(7)=1.46
    par_b(8)=1.44
    par_b(9)=1.41
    par_c(1)=1.19
    par_c(2)=1.27
    par_c(3)=1.34
    par_c(4)=1.45
    par_c(5)=1.58
    par_c(6)=1.80
    par_c(7)=1.97
    par_c(8)=2.24
    par_c(9)=2.44
    der_aa(1)=0.00
    der_aa(2)=0.50
    der_aa(3)=-1.56
    der_aa(4)=3.05
    der_aa(5)=-2.95
    der_aa(6)=1.07
    der_aa(7)=-0.71
    der_aa(8)=0.21
    der_aa(9)=0.00
    der_a(1)=0.00
    der_a(2)=1.19
    der_a(3)=-6.34
    der_a(4)=21.36
    der_a(5)=-10.95
    der_a(6)=2.59
    der_a(7)=-0.85
    der_a(8)=-2.07
    der_a(9)=0.00
    der_b(1)=0.00
    der_b(2)=-1.08
    der_b(3)=12.61
    der_b(4)=-20.96
    der_b(5)=24.08
    der_b(6)=-6.64
    der_b(7)=3.84
    der_b(8)=-2.09
    der_b(9)=0.00
    der_c(1)=0.00
    der_c(2)=0.94
    der_c(3)=-0.43
    der_c(4)=4.61
    der_c(5)=0.01
    der_c(6)=1.21
    der_c(7)=1.43
    der_c(8)=0.33
    der_c(9)=0.00

    do i=1,9
        del(i)=log10(del(i))
    enddo



    ! radius of shell of mass M with density rhom0
    R = (0.75_dl*M/Pi/rhom0)**(1._dl/3._dl) ! R in units of h^-1 Mpc
    ! check if powerspectrum for sigma(R) is actually calculated to large
    ! enough k
    ! maybe put this later outside dndlnM, for lowest mass limit and test ahead

    ! dM/dR
    dMdR = 3*M/R

    sR =cosmopar%sigmaR%Value(R)
    !    write(*,*) cosmopar%omegam,sR
    !g = delta(z)
    ! Parameters from Jenkins et al. for LCDM (Appendix) and SO(324)
    ! Amf = 0.316 ; Bmf = 0.67 ; epsmf = 3.82
    ! Note SO(324) is mass in mass density units or 0.3*324=97.2

    dsoz=massfnpar%dso/Omegam(z)
    if (massfnpar%psind==1) then  !Jenkins et al.
        fJen = massfnpar%Amf*exp(-abs(-dlog(g*sR)+massfnpar%Bmf)**massfnpar%epsmf)
        dndlnM = -rhom0*fJen*cosmopar%sigmaR%Derivative(R)/dMdR/sR
    elseif (massfnpar%psind==2) then  !Tinker et al.
        call SPLINTNR(del,par_aa,der_aa,total,log10(dsoz),par1)
        call SPLINTNR(del,par_a,der_a,total,log10(dsoz),par2)
        call SPLINTNR(del,par_b,der_b,total,log10(dsoz),par3)
        call SPLINTNR(del,par_c,der_c,total,log10(dsoz),par4)

        alpha=10**(-((0.75_dl/dlog10(dsoz/75.0_dl))**1.2_dl))
        massfnpar%Amf = par1*((1.0_dl+z)**(-0.14_dl))
        massfnpar%epsmf = par2*((1.0_dl+z)**(-0.06_dl))
        massfnpar%Bmf = par3*((1.0_dl+z)**(-alpha))
        massfnpar%pmf = par4

        fTink = massfnpar%Amf*((g*sR/massfnpar%Bmf)**(-massfnpar%epsmf)+1.0_dl)*exp(-massfnpar%pmf/sR/sR/g/g)

        dndlnM = -rhom0*fTink*cosmopar%sigmaR%Derivative(R)/dMdR/sR

    elseif (massfnpar%psind==3) then !Watson et al. 2013

        !FOF

        par1=0.282_dl !A
        par2=2.163_dl !alpha
        par3=1.406_dl !beta
        par4=1.210_dl !gamma


        !CPMSO
        !           par1=0.287_dl !A
        !           par2=2.234_dl !alpha
        !           par3=1.478_dl !beta
        !           par4=1.318_dl !gamma



        !!$           !redshift dependence:
        !par1=Omegam(z)*(0.990_dl*(1.0_dl+z)**(-3.216_dl)+0.074_dl) !A
        !           par1=Omegam(z)*(1.097_dl*(1.0_dl+z)**(-3.216_dl)+0.074_dl) !gaz
        !           par2=Omegam(z)*(5.907_dl*(1.0_dl+z)**(-3.058_dl)+2.349_dl) !bz
        !           par3=Omegam(z)*(3.136_dl*(1.0_dl+z)**(-3.599_dl)+2.344_dl) !az

        massfnpar%Amf = par1
        massfnpar%epsmf = par2
        massfnpar%Bmf = par3
        massfnpar%pmf = par4

        fWat = massfnpar%Amf*((g*sR/massfnpar%Bmf)**(-massfnpar%epsmf)+1.0_dl)*exp(-massfnpar%pmf/sR/sR/g/g)



        !Delta dependence:
        dsoz=massfnpar%dso/Omegam(z)
        ddz=-0.456_dl*Omegam(z)-0.139
        CDelta=exp(0.023_dl*(dsoz/178._dl-1.0_dl))*0.947
        p=0.072_dl
        q=2.130_dl


        Gamma=CDelta*(dsoz/178._dl)**ddz*exp(p*(1.0_dl-dsoz/178._dl)/(sR*g)**q)

        fwat=fwat*Gamma

        dndlnM = -rhom0*fWat*cosmopar%sigmaR%Derivative(R)/dMdR/sR
    else
        write(*,*) 'Invalid mass function indicator: ',massfnpar%psind
        stop
    end if
    RETURN
    end function dndlnM


    SUBROUTINE SPLINTNR(XA,YA,Y2A,N,X,Y)
    INTEGER :: N
    REAL(DL) :: XA(N),YA(N),Y2A(N),X,Y,H,A,B
    INTEGER :: KLO,KHI,K
    KLO=1
    KHI=N
1   IF (KHI-KLO.GT.1) THEN
        K=(KHI+KLO)/2
        IF(XA(K).GT.X)THEN
            KHI=K
        ELSE
            KLO=K
        ENDIF
        GOTO 1
    ENDIF
    H=XA(KHI)-XA(KLO)
    IF (H.EQ.0.) stop 'Bad XA input.'
    A=(XA(KHI)-X)/H
    B=(X-XA(KLO))/H
    Y=A*YA(KLO)+B*YA(KHI)+((A**3-A)*Y2A(KLO)+(B**3-B)*Y2A(KHI))*(H**2)/6.
    RETURN
    END SUBROUTINE SPLINTNR

    end module massfunction

    ! ####################################################################


    module numbercounts
    USE PRECISION
    USE CONSTANTS_sz
    USE POWER_SZ
    USE MASSOBSERVABLE
    USE MASSFUNCTION
    IMPLICIT NONE

    contains

    function next_z(zi,binz)
    real(dl)::zi,binz,next_z,dzi,hr

    !bins in redshifts are defined with higher resolution for z<0.2
    hr=0.2
    if (zi <hr) then
        dzi=1.e-3
    else if ((zi >=hr) .and. (zi <=1.)) then
        dzi=1.e-2
    else
        dzi=binz
    endif
    next_z=zi+dzi
    return
    end function next_z


    SUBROUTINE deltaN_yz(Z,Nz,LOGY,Ny,DN,DN2D,skyfracs,thetas,ylims,switch,qa_ytot,erf_list)
    !subroutine computing counts in y and z
    REAL(dl),INTENT(IN) :: z(:),logy(:),thetas(:),skyfracs(:),ylims(:,:),qa_ytot(:,:),erf_list(:,:)
    INTEGER, INTENT(IN):: Ny,Nz,switch
    REAl(dl),allocatable :: grid(:,:),steps_z(:),steps_z2(:),steps_m(:)
    REAl(dl),allocatable ::ytheta5_data(:,:),ytheta5_model(:),Mz_data(:,:)
    REAl(sp),allocatable ::completeness(:,:),completeness_2d(:,:,:)
    REAl(dl),allocatable ::dif1(:),dif2(:,:),mlims(:,:)
    REAL(dl)::  DN(:),DN2D(:,:)
    REAL(dl) :: rombint,DY,sigmaM,bias,dlogy,dz,z1,z2,z_ii,sum,yi
    REAL(dl) :: zi,dzi,lnm,lnmmin,lnmmax,binlogy,binz,window,zj,dzj,biny
    REAL(dl) :: ylim,lnylim,xi,xi1,sq_sigmaM,lny1,lny2,y1,y2,sum2,frac,col1,col2,zmax
    REAL(dl) :: z_min,z_max,m_zmin,m_zmax,zmin2,zmax2,m_zmin2,m_zmax2,m_max,m_min
    REAL(dl) :: test,obs,dif_t,max_z,min_z
    REAL(dl) :: zp(1)
    INTEGER :: N,iostat,nsum,nrows,reason,Nm
    INTEGER :: I,J,K,ii,jj,nsteps_z,nsteps_m,nsteps_z1
    real(dl), PARAMETER :: dlnm = 0.05_dl
    character (LEN=200)::filename
    integer :: iunit,ylim_choice
    integer ::k_ini,k_fin,nthetas,ntab,nsteps_z2,p1(1),p2(1),ini,npatches,p(1,1),pz,pm
    INTEGER, DIMENSION(8,2) :: values_time
    REAL(SP) :: clock_time


    ntab=size(ylims)
    nthetas=size(thetas)
    npatches=size(skyfracs)

    if (npatches==0) npatches=1 !constant ylim case

    lnmmin=31.
    lnmmax=37.
    binz=z(2)-z(1)
    biny=logy(2)-logy(1)

    !size of grid for mass
    nsteps_m=(lnmmax-lnmmin)/dlnm

    !size of grid for z
    !zmax=Z(Nz)+0.5_dl*binz
    min_z=Z(1)-0.5_dl*binz
    max_z=Z(Nz)+0.5_dl*binz

    zi=min_z
    if (zi <0) zi=0.

    nsteps_z=0
    do while (zi <= max_z)
        zi=next_z(zi,binz)
        nsteps_z=nsteps_z+1
    enddo

    allocate(steps_m(nsteps_m),steps_z(nsteps_z),stat=iostat)
    if (iostat/=0) then
        print*,'allocation error'
    endif


    !grid for mass
    lnm=lnmmin
    do i=1,nsteps_m
        steps_m(i)=lnm+dlnm/2.
        lnm=lnm+dlnm
    enddo

    !grid for z
    zi=min_z
    if (zi <0) zi=0.

    do i=1,nsteps_z
        steps_z(i)=zi
        zi=next_z(zi,binz)
    enddo


    if (steps_z(1) ==0) steps_z(1)=1.e-5


    ! from first bin in z to Nz
    SELECT CASE(SWITCH)
    CASE(1)

        allocate(grid(nsteps_m,nsteps_z),completeness(nsteps_m,nsteps_z),stat=iostat)
        if (iostat/=0) then
            print*,'allocation error'
            stop
        endif
        completeness(:,:)=0.


        call grid_C(npatches,steps_z,steps_m,thetas,ylims,skyfracs,completeness,qa_ytot,erf_list) !tabulate completeness for m and z


        grid(:,:)=0.
        call get_grid(nsteps_z,nsteps_m,steps_z,steps_m,grid) !tabulate number counts for m and z
        DN(:)=0.
        DO I=1,Nz
            !limits of bin in z
            z1=Z(I)-0.5_dl*binz
            z2=Z(I)+0.5_dl*binz

            call integrate_m_z(grid,skyfracs,completeness,steps_z,steps_m,nsteps_z,nsteps_m,z1,z2,binz,dlnm,sum)
            ! integrate on mass accounting for completeness
            DN(I)=sum
        ENDDO

        deallocate(grid,steps_m,steps_z,completeness,stat=iostat)
        if (iostat/=0) then
            print*,'deallocation error'
            stop
        endif
    CASE(2)

        allocate(grid(nsteps_m,nsteps_z),completeness_2d(nsteps_m,nsteps_z,ny+1),stat=iostat)
        if (iostat/=0) then
            print*,'allocation error'
            stop
        endif

        completeness_2d(:,:,:)=0.


        !$OMP PARALLEL
        !$OMP do
        do J=1,Ny+1
            call grid_C_2d(npatches,steps_z,steps_m,thetas,ylims,skyfracs,completeness_2d,logy,qa_ytot,erf_list,J)
            !tabulate completeness for m and z and sky patch
        ENDDO
        !$OMP end PARALLEL

        grid(:,:)=0.
        call get_grid(nsteps_z,nsteps_m,steps_z,steps_m,grid) !tabulate number counts for m and z
        DN(:)=0.


        !$OMP PARALLEL
        !$OMP do
        do J=1,Ny+1
            DO I=1,Nz+1
                !limits of bin in z
                call integrate_m_zq(grid,skyfracs,completeness_2d,steps_z,steps_m,nsteps_z,&
                    nsteps_m,Z(I)-0.5_dl*binz,Z(I)+0.5_dl*binz,binz,J,dlnm,DN2D(I,J))
                ! integrate on mass and z within bins accounting for completeness
            enddo
        ENDDO
        !$OMP end PARALLEL

        deallocate(grid,steps_m,steps_z,completeness_2d,stat=iostat)
        if (iostat/=0) then
            print*,'deallocation error'
            stop
        endif


    END SELECT


    RETURN
    end SUBROUTINE DeltaN_yz


    SUBROUTINE integrate_m_z(grid,skyfracs,compl,steps_z,steps_m,nsteps_z,nsteps_m,z1,z2,binz,dlnm,sum)
    ! integration for 1D likelihood
    REAl(dl),intent(in) ::grid(:,:),skyfracs(:)
    REAL (sp),intent(in)::compl(:,:)
    REAl(dl),intent(in) ::steps_z(:),steps_m(:),z1,z2,binz,dlnm
    REAl(dl):: sum,xi,xi1,zi,sigmaM,sq_sigmaM,lnm,window,bias,sum2,dz,sum3
    REAl(dl)::dnnew,dnold,lny1_new,lnylim,angsize,counts
    REAL (dl)::dif,dif_t,frac,ylim,frac_same,limit_m,test_new,h,C
    INTEGER::nsteps_z,nsteps_m,nsum,nthetas,col,k_ini,k_fin,t,k_ininew,col_old,npatches
    REAl(dl):: test(nsteps_z),c1,c2,f1,f2,x1,x2
    INTEGER::k,kk,JJ,II,j1,j2,i,index
    INTEGER::p(1)


    sum=0._dl
    npatches=size(skyfracs)


    ! find range i1->i2 for integration in z
    !allocate(test(nsteps_z))
    test=abs(steps_z-z1)
    p=minloc(test)
    j1=p(1)
    test=abs(steps_z-z2)
    p=minloc(test)
    j2=p(1)



    sum2=0.
    do jj=j1,j2-1
        DO ii=1,nsteps_m
            x1=steps_z(jj)
            x2=steps_z(jj+1)
            f1=grid(ii,jj)
            f2=grid(ii,jj+1)
            c1=compl(ii,jj)
            c2=compl(ii,jj+1)
            sum2=sum2+0.5*(f1*c1+f2*c2)*(x2-x1)*dlnm
        enddo
    enddo
    sum=sum+sum2



    END SUBROUTINE integrate_m_z

    SUBROUTINE integrate_m_zq(grid,skyfracs,compl,steps_z,steps_m,nsteps_z,nsteps_m,z1,z2,binz,biny,dlnm,sum)
    ! integration for 2D likelihood
    REAl(dl),intent(in) ::grid(:,:),skyfracs(:)
    REAL (sp),intent(in)::compl(:,:,:)
    REAl(dl),intent(in) ::steps_z(:),steps_m(:),z1,z2,binz,dlnm
    REAl(dl):: sum,xi,xi1,zi,sigmaM,sq_sigmaM,lnm,window,bias,sum2,dz,sum3
    REAl(dl)::dnnew,dnold,lny1_new,lnylim,angsize,counts
    REAL (dl)::dif,dif_t,frac,ylim,frac_same,limit_m,test_new,h,C
    INTEGER::nsteps_z,nsteps_m,nsum,nthetas,col,k_ini,k_fin,t,k_ininew,col_old,npatches,biny
    REAl(dl):: test(nsteps_z),c1,c2,f1,f2,x1,x2
    INTEGER::k,kk,JJ,II,j1,j2,i,index
    INTEGER::p(1)

    !compl(nteps_m,nsteps_z,npatches)
    sum=0._dl
    npatches=size(skyfracs)


    ! find range i1->i2 for integration in z
    !allocate(test(nsteps_z))
    test=abs(steps_z-z1)
    p=minloc(test)
    j1=p(1)
    test=abs(steps_z-z2)
    p=minloc(test)
    j2=p(1)

    sum2=0.
    do jj=j1,j2-1
        DO ii=1,nsteps_m
            x1=steps_z(jj)
            x2=steps_z(jj+1)
            f1=grid(ii,jj)
            f2=grid(ii,jj+1)
            c1=compl(ii,jj,biny)
            c2=compl(ii,jj+1,biny)
            sum2=sum2+0.5*(f1*c1+f2*c2)*(x2-x1)*dlnm
        enddo
    enddo
    sum=sum+sum2


    END SUBROUTINE integrate_m_zq


    SUBROUTINE grid_C_2d(npatches,steps_z,steps_m,thetas,ylims,skyfracs,completeness,logy,qa_ytot,erf_list,iy)
    !tabulate completeness for 2D likelihood
    REAl(dl),intent(in) :: steps_z(:),steps_m(:),thetas(:),ylims(:,:),skyfracs(:),logy(:),qa_ytot(:,:),erf_list(:,:)
    INTEGER :: nsteps_z,nsteps_m,nthetas,ntab,npatches,iostat,switch_comp,nsteps_y,nsteps_t,nerf,iy
    real(sp):: completeness(:,:,:) !different completeness for different bin in q
    integer ::i,j,ii,jj,index1,index2,count,P(1),k1,k2,l1,l2,k,N,nthetas2,nrows,indt(1),indy(1),it,i1y,i2y,i3y,kk
    integer ::indminy(2),indmaxy(2),iminy,imaxy,nq
    real (dl):: dif_old,dif,max,min,dlm,binz,m_min,m_max,mp,yp,zp,thp,xk1,xk2,xk3,yk1,yk2,yk3,fact,qmin,qmax,dlogy
    real(dl),allocatable:: dif_y(:),dif_theta(:),difft(:),diffy(:)
    real(dl):: min_thetas,max_thetas,min_y,max_y
    real(dl):: c1,c2,th1,th2,c,y12,y1,y2,y,col1,col0,csum,ytiny
    real(dl),allocatable::thetas2(:),ylims2(:,:),erfs(:,:,:)
    real(dl)::win0,win,arg,arg0,y0,py,lnymax,lnymin,lny,dy,fac,mu,int,dlny,fsky
    real(dl):: sigmaM

    sigmaM=cosmopar%sigmaM
    switch_comp = 1 !ERF  (2D not implemented with QA selection function)
    ntab=size(ylims)
    nsteps_z=size(steps_z)
    nsteps_m=size(steps_m)
    nthetas=size(thetas)
    dlogy=logy(2)-logy(1)
    nq=size(logy)
    allocate(dif_y(ntab),dif_theta(nthetas))
    completeness(:,:,iy)=0.

    min_thetas=minval(thetas)
    max_thetas=maxval(thetas)

    min_y=minval(qa_ytot)
    max_y=maxval(qa_ytot)

    nerf =size(qa_ytot,1)


    if (sigmaM==0) then

        do jj=1,nsteps_z
            do ii=1,nsteps_m
                mp=exp(steps_m(ii))
                zp=steps_z(jj)
                thp=theta500(mp,zp)
                yp=y500(mp,zp)
                if (thp > max_thetas) then
                    l1=nthetas
                    l2=nthetas-1
                    th1=thetas(l1)
                    th2=thetas(l2)

                else if  (thp < min_thetas) then
                    l1=1
                    l2=2
                    th1=thetas(l1)
                    th2=thetas(l2)
                else
                    dif_theta=abs(thetas-thp)
                    P=minloc(dif_theta)
                    l1=P(1)
                    th1=thetas(l1)
                    l2=l1+1
                    if (th1 > thp) l2=l1-1
                    th2=thetas(l2)

                endif



                do i=1,npatches
                    y1=ylims(i,l1)
                    y2=ylims(i,l2)
                    y=y1+(y2-y1)/(th2-th1)*(thp-th1)!sigma at the relevant scale for the patch

                    k=iy
                    qmin=logy(k)-dlogy/2.
                    qmax=logy(k)+dlogy/2.
                    qmin=10.**qmin
                    qmax=10.**qmax
                    c2=erf_compl(yp,y,q)*erf_compl(yp,y,qmin)*(1.-erf_compl(yp,y,qmax))
                    if (k==1) c2=erf_compl(yp,y,q)*(1.-erf_compl(yp,y,qmax))
                    if (k==nq) c2=erf_compl(yp,y,q)*erf_compl(yp,y,qmin)
                    completeness(ii,jj,k)=completeness(ii,jj,k)+c2*skyfracs(i)
                enddo
            enddo
        enddo

    else
        fac=1./sqrt(2.*pi*sigmaM**2)
        lnymin=-11.5
        lnymax=10.
        dlny=0.05

        N=(lnymax-lnymin)/dlny

        fsky=0
        do i=1,npatches
            fsky=fsky+skyfracs(i)
        enddo

        allocate(erfs(N,nthetas,nq),stat=iostat)!y,integrated completeness
        if (iostat/=0) then
            print*,'allocation error'
        endif

        erfs(:,:,:)=0.
        do j=1,nthetas
            lny=lnymin
            do jj=1,N
                y0=dexp(lny)
                lny=lny+dlny

                do i=1,npatches
                    y1=ylims(i,j)


                    k=iy
                    qmin=logy(k)-dlogy/2.
                    qmax=logy(k)+dlogy/2.
                    qmin=10.**qmin
                    qmax=10.**qmax
                    c2=erf_compl(y0,y1,q)*erf_compl(y0,y1,qmin)*(1.-erf_compl(y0,y1,qmax))
                    if (k==1)  c2=erf_compl(y0,y1,q)*(1.-erf_compl(y0,y1,qmax))
                    if (k==nq) c2=erf_compl(y0,y1,qmin)*erf_compl(y0,y1,q)

                    erfs(jj,j,k)=erfs(jj,j,k)+c2*skyfracs(i)


                enddo
            enddo
        enddo


        do jj=1,nsteps_z
            do ii=1,nsteps_m

                mp=exp(steps_m(ii))
                zp=steps_z(jj)
                thp=theta500(mp,zp)
                yp=y500(mp,zp)
                if (thp > max_thetas) then
                    l1=nthetas
                    l2=nthetas-1
                    th1=thetas(l1)
                    th2=thetas(l2)
                else if  (thp < min_thetas) then
                    l1=1
                    l2=2
                    th1=thetas(l1)
                    th2=thetas(l2)
                else
                    dif_theta=abs(thetas-thp)
                    P=minloc(dif_theta)
                    l1=P(1)
                    th1=thetas(l1)
                    l2=l1+1
                    if (th1 > thp) l2=l1-1
                    th2=thetas(l2)
                endif
                y=y500(mp,zp)
                mu=dlog(y500(mp,zp))


                kk=iy
                int=0.
                lny=lnymin
                do k=1,N-1
                    y0=dexp(lny)
                    y=dexp(lny+dlny)
                    dy=y-y0
                    arg0=((lny-mu)/(sqrt(2.)*sigmaM))
                    win0=erfs(k,l1,kk)+(erfs(k,l2,kk)-erfs(k,l1,kk))/(th2-th1)*(thp-th1)
                    win=erfs(k+1,l1,kk)+(erfs(k+1,l2,kk)-erfs(k+1,l1,kk))/(th2-th1)*(thp-th1)
                    lny=lny+dlny
                    arg=((lny-mu)/(sqrt(2.)*sigmaM))
                    py=(win0*fac/y0*exp(-arg0**2)+win*fac/y*exp(-arg**2))*0.5
                    int=int+py*dy
                enddo
                if (int > fsky) int=fsky
                if (int < 0.) int=0.
                completeness(ii,jj,kk)=int
            enddo
        enddo

        deallocate(erfs)

    endif


    END SUBROUTINE grid_C_2d

    SUBROUTINE grid_C(npatches,steps_z,steps_m,thetas,ylims,skyfracs,completeness,qa_ytot,erf_list)
    !tabulate completeness for 1D likelihood
    REAl(dl),intent(in) :: steps_z(:),steps_m(:),thetas(:),ylims(:,:),skyfracs(:),qa_ytot(:,:),erf_list(:,:)
    INTEGER :: nsteps_z,nsteps_m,nthetas,ntab,npatches,iostat,switch_comp,nsteps_y,nsteps_t,nerf
    real(sp):: completeness(:,:)

    integer ::i,j,ii,jj,index1,index2,count,P(1),k1,k2,l1,l2,k,N,nthetas2,nrows,indt(1),indy(1),it,i1y,i2y,i3y
    integer ::indminy(2),indmaxy(2),iminy,imaxy
    real (dl):: dif_old,dif,max,min,dlm,binz,m_min,m_max,mp,yp,zp,thp,xk1,xk2,xk3,yk1,yk2,yk3,fact
    real(dl),allocatable:: dif_y(:),dif_theta(:),difft(:),diffy(:)
    real(dl):: min_thetas,max_thetas,min_y,max_y
    real(dl):: c1,c2,th1,th2,c,y12,y1,y2,y,col1,col0
    real(dl),allocatable::thetas2(:),ylims2(:,:),erfs(:,:)
    real(dl)::win0,win,arg,arg0,y0,py,lnymax,lnymin,lny,dy,fac,mu,int,dlny,fsky
    real(dl):: sigmaM

    sigmaM=cosmopar%sigmaM
    !choice to use error function completeness or QA
    switch_comp = 1 !ERF
    !switch_comp = 2 !QA
    !end choice to use error function completeness or QA

    ntab=size(ylims)
    nsteps_z=size(steps_z)
    nsteps_m=size(steps_m)
    nthetas=size(thetas)

    allocate(dif_y(ntab),dif_theta(nthetas))

    min_thetas=minval(thetas)
    max_thetas=maxval(thetas)

    min_y=minval(qa_ytot)
    max_y=maxval(qa_ytot)


    nerf = size(qa_ytot,1)

    SELECT CASE(switch_comp)
    CASE(1)
        !error function
        if (sigmaM==0) then
            ! no scatter in y-m relation
            do jj=1,nsteps_z
                do ii=1,nsteps_m
                    mp=exp(steps_m(ii))
                    zp=steps_z(jj)
                    thp=theta500(mp,zp)
                    yp=y500(mp,zp)
                    if (thp > max_thetas) then
                        l1=nthetas
                        l2=nthetas-1
                        th1=thetas(l1)
                        th2=thetas(l2)

                    else if  (thp < min_thetas) then
                        l1=1
                        l2=2
                        th1=thetas(l1)
                        th2=thetas(l2)
                    else
                        dif_theta=abs(thetas-thp)
                        P=minloc(dif_theta)
                        l1=P(1)
                        th1=thetas(l1)
                        l2=l1+1
                        if (th1 > thp) l2=l1-1
                        th2=thetas(l2)

                    endif
                    completeness(ii,jj)=0.

                    do i=1,npatches
                        y1=ylims(i,l1)
                        y2=ylims(i,l2)
                        y=y1+(y2-y1)/(th2-th1)*(thp-th1)!sigma at the relevant scale for the patch
                        c2=erf_compl(yp,y,q)
                        completeness(ii,jj)=completeness(ii,jj)+c2*skyfracs(i)

                    enddo
                enddo
            enddo



        else
            ! scatter in y-m relation
            fac=1./sqrt(2.*pi*sigmaM**2)
            lnymin=-11.5
            lnymax=10.
            dlny=0.05

            N=(lnymax-lnymin)/dlny

            fsky=0
            do i=1,npatches
                fsky=fsky+skyfracs(i)
            enddo

            allocate(erfs(N,nthetas),stat=iostat)!y,integrated completeness
            if (iostat/=0) then
                print*,'allocation error'
            endif


            do j=1,nthetas
                lny=lnymin
                do k=1,N
                    y0=dexp(lny)
                    lny=lny+dlny
                    win0=0.
                    do i=1,npatches
                        y1=ylims(i,j)
                        win0=win0+erf_compl(y0,y1,q)*skyfracs(i)
                    enddo
                    erfs(k,j)=win0
                    !           if (erfs(k,j)>fsky) erfs(k,j)=fsky
                enddo
            enddo



            do jj=1,nsteps_z
                do ii=1,nsteps_m

                    mp=exp(steps_m(ii))
                    zp=steps_z(jj)
                    thp=theta500(mp,zp)
                    yp=y500(mp,zp)

                    if (thp > max_thetas) then
                        l1=nthetas
                        l2=nthetas-1
                        th1=thetas(l1)
                        th2=thetas(l2)
                    else if  (thp < min_thetas) then
                        l1=1
                        l2=2
                        th1=thetas(l1)
                        th2=thetas(l2)
                    else
                        dif_theta=abs(thetas-thp)
                        P=minloc(dif_theta)
                        l1=P(1)
                        th1=thetas(l1)
                        l2=l1+1
                        if (th1 > thp) l2=l1-1
                        th2=thetas(l2)
                    endif
                    y=y500(mp,zp)
                    mu=dlog(y500(mp,zp))

                    int=0.
                    lny=lnymin
                    do k=1,N-1
                        y0=dexp(lny)
                        y=dexp(lny+dlny)
                        dy=y-y0
                        arg0=((lny-mu)/(sqrt(2.)*sigmaM))
                        win0=erfs(k,l1)+(erfs(k,l2)-erfs(k,l1))/(th2-th1)*(thp-th1)
                        win=erfs(k+1,l1)+(erfs(k+1,l2)-erfs(k+1,l1))/(th2-th1)*(thp-th1)
                        lny=lny+dlny
                        arg=((lny-mu)/(sqrt(2.)*sigmaM))
                        py=(win0*fac/y0*exp(-arg0**2)+win*fac/y*exp(-arg**2))*0.5
                        int=int+py*dy
                    enddo
                    if (int > fsky) int=fsky
                    if (int < 0.) int=0.
                    completeness(ii,jj)=int

                enddo
            enddo

            deallocate(erfs)

        endif

    CASE(2)
        !QA COMPLETENESS
        if (sigmaM==0) then
            ! no scatter in y-m relation
            do jj=1,nsteps_z
                do ii=1,nsteps_m

                    mp=exp(steps_m(ii))
                    zp=steps_z(jj)

                    thp=theta500(mp,zp)
                    yp=y500(mp,zp)

                    allocate(difft(nthetas))
                    difft = abs(thetas-thp)
                    indt = minloc(difft)
                    it = indt(1)

                    allocate(diffy(nerf))
                    diffy = abs(qa_ytot(:,it)-yp)

                    indy = minloc(diffy)

                    i1y = indy(1)-1
                    i2y = indy(1)+1

                    deallocate(difft)
                    deallocate(diffy)

                    indminy = minloc(qa_ytot)
                    indmaxy = maxloc(qa_ytot)

                    if (i1y <=0) then
                        i1y = indminy(1)
                    endif

                    if (i2y>nerf) then
                        i2y = indmaxy(1)
                    endif

                    if (yp < min_y) then
                        i1y = indminy(1)
                        i2y = indminy(1)+1
                    endif

                    if (yp > max_y) then
                        i1y = indmaxy(1)-1
                        i2y = indmaxy(1)
                    endif

                    xk1 = qa_ytot(i1y,it)
                    xk2 = qa_ytot(i2y,it)
                    yk1 = erf_list(i1y,it)
                    yk2 = erf_list(i2y,it)

                    completeness(ii,jj) = 0.

                    !INTERPOLATION POLY ORDER 1
                    completeness(ii,jj) = exp(log(yk1)+((log(yp)-log(xk1))/(log(xk2)-log(xk1)))*(log(yk2)-log(yk1)))

                    if (completeness(ii,jj)> 1.) completeness(ii,jj)=1.
                    if (completeness(ii,jj)==0.) completeness(ii,jj)=1e-20

                    fact = sum(skyfracs)
                    completeness(ii,jj) = completeness(ii,jj)*fact

                enddo
            enddo
        else
            ! scatter in y-m relation
            print*, 'QA completeness and scatter in Y-M relation'
            print*, 'CASE NOT IMPLEMENTED - STOPPING'
            stop
        endif

    END SELECT

    END SUBROUTINE grid_C

    SUBROUTINE get_grid(nz,nm,z,lnm,grid)
    !tabulate counts from theory
    REAL(dl):: z(:),lnm(:),grid(:,:)
    INTEGER :: I,J,nz,nm
    REAL(dl):: Mnew,g,ynew,vol,theta



    DO I=1,Nz
        g = delta(z(I))
        vol=dVdzdO(z(I))
        DO J=1,Nm
            Mnew=exp(lnm(J))
            grid(j,i) = dndlnM(z(I),Mnew,g)*surveypar%deg2*vol
        ENDDO
    ENDDO

    END SUBROUTINE get_grid

    end module numbercounts


    ! ####################################################################



    module szcounts
    use settings
    use CosmologyTypes
    use CosmoTheory
    use Calculator_Cosmology
    use Likelihood_Cosmology
    use IniObjects
    use massfunction
    use massobservable

    implicit none
    private
    logical :: do_sz_init = .true.
    logical :: use_sz = .false.
    integer :: SZ_num
    integer :: Ny,Nz
    integer :: SZ_switch  !1 for N(z), 2 for N(y), 3 for N(z,y)
    integer :: nmiss_switch  !0 for simple rescaling, 1 for MCMC
    integer :: errorz_switch  !0 for simple rescaling, 1 for MCMC
    integer :: print_counts !to print theory counts (for action=4)
    integer :: nmiss,nred2,ncat,nerf,nrows_qa,nrows_erf
    integer :: ylim_switch  !>1 for constant ylim, <1 for variable ylim
    real(dl), allocatable :: DNcat(:,:),DNzcat(:),DNycat(:),DN(:,:),DNz(:),DNy(:)
    real(dl), allocatable :: Z(:),LOGY(:),ylims(:,:),thetas(:),skyfracs(:),erf_list(:,:),qa_ytot(:,:)
    real(dl), allocatable :: SZcat(:,:)
    real(dl) :: clash,wtg,lens,cccp,pns,oh2,palpha,pystar,psigma,pbeta ! switches for priors =0 no prior, =1 prior
    real :: sz_kmax = 4.0
    character(len=256) :: SZ_filename = ''
    character(LEN=*), parameter :: SZ_version =  'Nov_2014'
    REAL(dl), SAVE :: zmaxs,z0,dz
    INTEGER, SAVE :: Nred, Nscat

    type, extends(TCosmoCalcLikelihood) :: SZLikelihood
    contains
    procedure :: LogLike => SZCC_Cash
    end type SZLikelihood

    PUBLIC :: SZLikelihood_Add, SZcc_Cash

    contains


    subroutine SZLikelihood_Add(LikeList, Ini)
    ! interface of the module with cosmomcplanck code
    ! user-defined settings are interfaced here
    implicit none
    class(TLikelihoodList) :: LikeList
    class(TSettingIni) :: ini
    integer count
    Type(SZLikelihood), pointer :: this

    print_counts=0

    if (Ini%Read_Int('action')==4) print_counts=1
    if (Ini%Read_Logical('use_SZ',.false.)) then
        allocate(this)
        this%needs_background_functions = .true.
        this%needs_powerspectra = .true.
        this%kmax = sz_kmax
        this%needs_sigmaR = .true.
        this%version = SZ_version
        call this%loadParamNames(trim(DataDir)//'SZ.paramnames')
        call LikeList%Add(this)
        this%LikelihoodType = 'SZ'
        this%name='SZ'


        clash=0.  !swithes for the priors (get multiplied to priors in SZCC_Cash)
        wtg=0.    !either 0 or 1
        lens=0.
        cccp=0.
        oh2=0.
        pns=0.
        palpha=0.
        pystar=0.
        psigma=0.
        pbeta=0.

        !
        sz_switch=0 !initialization

        if (Ini%Read_Logical('1D',.false.)) then
            sz_switch=1
            print*,'1D SZ likelihood dN/dz'
        endif
        if ((Ini%Read_Logical('2D',.false.)) .and. sz_switch==1 ) then
            print*,'Error 1D and 2D likelihood both selected'
            stop
        endif

        if (Ini%Read_Logical('2D',.false.)) then
            sz_switch=2
            print*,'2D SZ likelihood dN/dzdq'
        endif

        if (sz_switch==0) then
            sz_switch=2 !default
            print*,'2D SZ likelihood dN/dzdq'
        endif

        if (Ini%Read_Logical('prior_alpha_SZ',.false.)) then
            palpha=1.
            print*,'prior on alpha_SZ'
        endif

        if (Ini%Read_Logical('prior_ystar_SZ',.false.)) then
            pystar=1.
            print*,'prior on ystar_SZ'
        endif

        if (Ini%Read_Logical('prior_scatter_SZ',.false.)) then
            psigma=1.
            print*,'prior on scatter_SZ'
        endif

        if (Ini%Read_Logical('prior_beta_SZ',.false.)) then
            pbeta=1.
            print*,'prior on beta_SZ'
        endif

        !!$        if (Ini%Read_Logical('prior_clash',.false.)) then
        !!$            clash=1.
        !!$            print*,'CLASH prior on SZ mass bias'
        !!$        endif

        if (Ini%Read_Logical('prior_wtg',.false.)) then
            wtg=1.
            print*,'WtG prior on SZ mass bias'
        endif

        if (Ini%Read_Logical('prior_lens',.false.)) then
            print*,'Planck lensing prior on SZ mass bias'
            lens=1.
        endif

        if (Ini%Read_Logical('prior_cccp',.false.)) then
            print*,'CCCP prior on SZ mass bias'
            cccp=1.
        endif

        if (Ini%Read_Logical('prior_ns',.false.))then
            print*,'Prior on ns'
            pns=1.
        endif
        if (Ini%Read_Logical('prior_omegabh2',.false.)) then
            oh2=1.
            print*,'Prior on omegabh2'
        endif
        massfnpar%psind=2 !mass function =tinker mass function
        if (Ini%Read_Logical('use_watson',.false.)) then
            massfnpar%psind=3
        endif

        CALL SZ_init

    endif
    end subroutine SZLikelihood_Add

    subroutine SZ_init
    !initialization of SZ module
    implicit none
    character (LEN=20):: name
    character (LEN=200)::cat_filename,skyfracs_filename,ylims_filename,thetas_filename,filename
    real :: dummy
    real :: dzorg
    integer :: i,j,jj,iostat,nrows,reason,ii,iunit,nrows_old,nthetas,npatches
    integer ::Nfact,error,nq
    real(dl) :: ymin,ymax,dlogy,logymax,logymin,yi,col1,col2,sum,col3,col4
    real(dl) :: y_min,y_max,z_min,z_max,factorial,qmin,qmax,dq,qbin,norm


    error=0
    !file names
    cat_filename='data/SZ_cat.txt'
    !cat_filename='data/SZ_cat_intersection.txt' !intersection catalogue
    thetas_filename='data/SZ_thetas.txt'
    skyfracs_filename='data/SZ_skyfracs.txt'
    ylims_filename='data/SZ_ylims.txt'

    massfnpar%dso=500.
    surveypar%ab=0.
    surveypar%nb=0.
    surveypar%sfid=0.
    Nscat=0
    surveypar%deg2=41253.0  !full sky (sky fraction handled in skyfracs file)

    z0=0.d0
    zmaxs=1.
    dz=0.1d0

    surveypar%ylimin=1.e-3
    surveypar%ymaxin=1.e-1
    logymin=0.7 !s/n=6
    logymax=1.5 !s/n=32. (higher s/n in the next to last bin)
    dlogy=0.25

    if (massfnpar%psind ==2) print*,'Using Tinker mass function'
    if (massfnpar%psind ==3) print*,'Using Watson mass function'


    ylim_switch=-1 !>0=constant ylim, <0=variable ylim
    nmiss_switch=0
    errorz_switch=0

    surveypar%deg2 = 3.046174198d-4*surveypar%deg2
    Nz = DINT((zmaxs-z0)/dz)+1
    Ny = DINT((logymax-logymin)/dlogy)+1
    print*,'Ny=',Ny

    !ylims file
    if (ylim_switch < 0) then

        ymin=-1.*ymin

        ! read catalogue
        filename=cat_filename

        nrows=0

        print*,'Reading catalogue'

        open (newunit=iunit,file=filename,status='old',form="formatted",iostat=reason)
        IF (Reason > 0)  THEN
            print*,'Error opening file',filename
        endif

        print*,filename
        nrows=1
        ncat=1 !number of clusters above s/m threshold

        DO
            READ(iunit,*,IOSTAT=Reason)  col1, col2, col3
            !print*,Reason
            IF (Reason > 0)  THEN
                print*,'Error in reading file'
                print*,filename
                stop
            ELSE IF (Reason < 0) THEN
                exit
            ELSE
                nrows=nrows+1
                if (col3 >=q) ncat=ncat+1
            END IF
        END DO

        CLOSE (unit=iunit)
        print*,'done'

        Print*,'Catalogue Number of clusters=',nrows
        Print*,'Catalogue Number of clusters above the S/N threshold of',q,'=',ncat
        ALLOCATE(SZcat(ncat,3),stat=iostat)! z,y
        if (iostat/=0) then
            print*,'allocation error'
        endif

        open (newunit=iunit,file=filename,status='old',form="formatted")
        ii=1
        DO i=1,nrows
            READ(iunit,*,IOSTAT=Reason)  col1, col2, col3
            if (col3 >=q) then
                SZcat(ii,1)=col1  !redshift
                SZcat(ii,2)=col2  !error on redshift
                SZcat(ii,3)=col3  !detection S/N
                ii=ii+1
            endif
        ENDDO
        CLOSE (unit=iunit)

        ! load files describing selection function

        !file with theta
        ! theta, nbins, first_index, last_index,first_index2, last_index2
        call File%LoadTxt(thetas_filename, thetas, n=nthetas)
        print*,'Number of size thetas=',nthetas

        !file with skyfracs
        ! theta, nbins, first_index, last_index,first_index2, last_index2
        call File%LoadTxt(skyfracs_filename, skyfracs, n=npatches)
        print*,'Number of patches=',npatches

        filename=ylims_filename
        nrows=File%TxtFileLines(filename)
        print*,'Number of size y =',nrows
        if (nrows /= npatches*nthetas) then
            print*,'Format error for ylims.txt:'
            print*,'Expected rows:',npatches*nthetas
            print*,'Actual rows:',nrows
            stop
        endif
        allocate(ylims(npatches,nthetas),stat=iostat)
        if (iostat/=0) then
            print*,'allocation error'
        endif

        open (newunit=iunit,file=filename,status='unknown',form="formatted")
        i=1
        j=1
        DO ii=1,nrows
            READ(iunit,*,IOSTAT=Reason)  col1
            ylims(i,j)=col1
            i=i+1
            if (i > npatches) then
                i=1
                j=j+1
            endif
        ENDDO
        CLOSE (unit=iunit)
    endif


    ALLOCATE(Z(Nz),LOGY(Ny+1),stat=iostat)
    if (iostat /= 0) then
        print *, "Cannot allocate work arrays"
        stop
    endif

    ! logy vector
    yi=logymin+dlogy/2.
    DO I=1,Ny+1
        logy(I)=yi
        yi=yi+dlogy
    END DO

    print*,'q=',10.d0**logy

    ! z vector
    DO I=0,Nz-1
        Z(I+1)=z0+I*dz+0.5_dl*dz
    END DO
    if (z0==0._dl) Z(1)=Z(1)+1.e-8 ! for numerical problem when starting from 0.

    ALLOCATE(DNcat(Nz,Ny+1),DNzcat(Nz),DNycat(Ny),DN(Nz+1,Ny+1),DNz(Nz),DNy(Ny),stat=iostat)
    if (iostat /= 0) then
        print *, "Cannot allocate work arrays"
        stop
    endif

    DNcat(:,:)=0.
    DNzcat(:)=0.

    SELECT CASE(SZ_SWITCH)
    CASE(1)      ! N(z)
        nrows=size(SZcat(:,1))

        ! exclude clusters below s/n threshold q

        nmiss=0
        DO ii=1,nrows
            if (SZcat(ii,1) <0.) nmiss=nmiss+1.
        enddo

        z_min=z0
        z_max=z_min+dz
        DO I=1,Nz
            DO ii=1,nrows
                if ((SZcat(ii,1) >= z_min) .and. (SZcat(ii,1) < z_max)) then
                    DNzcat(I)=DNzcat(I)+1.
                endif
            ENDDO
            z_min=z_min+dz
            z_max=z_max+dz
        END DO

        sum=0.
        DO I=1,Nz
            sum=sum+DNzcat(I)
        ENDDO

        nred2=nrows-nmiss
        print*,nrows,nmiss

        ncat=nrows

        print*,'Number of clusters:',ncat
        print*,'Number of clusters with redshift:',nred2
        print*,'Number of clusters with no redshift:',nmiss
        print*,'Counts:',DNzcat

        if (nmiss==0) nmiss_switch=0


        SELECT CASE(nmiss_switch)
        CASE(0)
            print*,'Rescaling for missing redshifts',dble(nrows)/dble(nred2)
        CASE(1)
            print*,'Randomizing for missing redshifts',dble(nrows)/dble(nred2)
        END SELECT
    CASE(2)
        !N(z,q)

        ! compute catalogue counts in z and q
        ! compute P(q|qm) once and for all

        nrows=size(SZcat(:,1))

        nmiss=0
        DO ii=1,nrows
            if (SZcat(ii,1) <0.d0) nmiss=nmiss+1.
        enddo

        z_min=z0
        z_max=z_min+dz
        DO I=1,Nz
            DO J=1,Ny
                y_min=logY(J)-dlogy/2.
                y_max=logY(J)+dlogy/2.
                y_min=10.d0**y_min
                y_max=10.d0**y_max
                DO ii=1,nrows
                    if ((SZcat(ii,1) >= z_min) .and. (SZcat(ii,1) < z_max) .and.&
                        (SZcat(ii,3) < y_max) .and. (SZcat(ii,3) >= y_min)) then
                    DNcat(I,J)=DNcat(I,J)+1.
                    endif
                ENDDO
            ENDDO
            J=Ny+1 ! the last bin contains all S/N greater than what in the previous bin
            y_min=y_max
            DO ii=1,nrows
                if ((SZcat(ii,1) >= z_min) .and. (SZcat(ii,1) < z_max) .and. (SZcat(ii,3) >= y_min)) then
                    DNcat(I,J)=DNcat(I,J)+1.
                endif
            ENDDO

            z_min=z_min+dz
            z_max=z_max+dz
        END DO

        !missing redshifts
        DO J=1,Ny
            y_min=logY(J)-dlogy/2.
            y_max=logY(J)+dlogy/2.
            y_min=10.**y_min
            y_max=10.**y_max
            DO ii=1,nrows
                if ((SZcat(ii,1) == -1) .and. (SZcat(ii,3) < y_max) .and. (SZcat(ii,3) >= y_min)) then

                    norm=0.
                    do jj=1,Nz
                        norm=norm+DNcat(jj,J)
                    enddo
                    DNcat(:,J)=DNcat(:,J)*(norm+1.d0)/norm
                endif
            ENDDO
        ENDDO
        J=Ny+1 ! the last bin contains all S/N greater than what in the previous bin
        y_min=y_max
        DO ii=1,nrows
            if ((SZcat(ii,1) == -1) .and. (SZcat(ii,3) >= y_min)) then
                norm=0.
                do jj=1,Nz
                    norm=norm+DNcat(jj,J)
                enddo
                DNcat(:,J)=DNcat(:,J)*(norm+1.d0)/norm
            endif
        ENDDO
        !end missing redshifts

        sum=0.
        DO I=1,Nz
            DO J=1,Ny+1
                sum=sum+DNcat(I,J)
            ENDDO
        END DO
        print*,'total cat',sum

        nred2=nrows-nmiss
        ncat=nrows

        print*,'Number of clusters:',ncat
        print*,'Number of clusters with redshift:',nred2
        print*,'Number of clusters with no redshift:',nmiss

    END SELECT

    do_sz_init = .false.

    print*,'End SZ initialization'

    end subroutine SZ_init



    function SZCC_Cash(this,CMB,Theory,DataParams)
    ! likelihood computation
    ! SZ nuisance in dataparams
    use cosmology
    use numbercounts
    use, intrinsic :: ieee_arithmetic
    Class(SZLikelihood) :: this
    Class (CMBParams):: CMB
    Class(TCosmoTheoryPredictions), target :: Theory
    real(mcp) DataParams(:)
    real(mcp)  SZCC_Cash
    INTEGER :: N,i,j,Nf,ii
    REAL(DL) :: sum,factorial,ln_factorial,fact!,SZCC_Cash_exp


    SZCC_Cash=logzero


    !Mapping of cosmo parameters
    cosmopar%H0=CMB%H0
    cosmopar%w0=CMB%w
    cosmopar%w1=0.0 ! Note, we do not support evolution of w right now
    cosmopar%omegam=CMB%omc+CMB%omb+CMB%omnu
    cosmopar%omegav=CMB%omv
    cosmopar%omegak=CMB%omk
    cosmopar%n=CMB%InitPower(ns_index)
    cosmopar%sig8=Theory%sigma_8!
    cosmopar%omegabh2=CMB%ombh2
    cosmopar%gamma=-1

    !Mapping of nuisance parameters
    cosmopar%alpha=DataParams(1)
    cosmopar%ystar=10.**(DataParams(2))/(2.**cosmopar%alpha)*0.00472724
    cosmopar%logystar=(DataParams(2))!to set the prior on logy
    cosmopar%bias=DataParams(3)
    cosmopar%biasinv=1./DataParams(3)
    cosmopar%sigmaM=DataParams(4)
    cosmopar%beta=DataParams(5)
    cosmopar%sigmaR=>Theory%sigma_R
    call INIGROWTH

    DN(:,:)=0.
    DNz(:)=0.
    DNy(:)=0.



    SELECT CASE(SZ_SWITCH)

    CASE(1) ! n(Z)

        call deltaN_yz(Z,Nz,LOGY,Ny,DNZ,DN,skyfracs,thetas,ylims,sz_switch,qa_ytot,erf_list)
        ! N(z)

        sum=0.
        DO I=1,Nz
            sum=sum+DNz(I)
        ENDDO
        if (print_counts==1) then
            print*,'predicted counts'
            print*,DNz
            print*,'total counts',sum
        endif


        if (ieee_is_nan(DNZ(1))) then
            print*,'NaN found in theory counts!'
            stop
        endif


        DNzcat=DNzcat/dble(nred2)*dble(ncat)
        sum=0.
        do i=1,Nz
            if (DNz(i) /= 0.) then
                ln_factorial=0.
                if (DNzcat(i) /= 0.) ln_factorial=0.918939+(DNzcat(i)+0.5)*dlog(DNzcat(i))-DNzcat(i) !Stirling
                sum=sum-1.*(DNzcat(i)*dlog(DNz(i))-DNz(i)-ln_factorial)
            end if
        end do
        SZCC_Cash=sum
        print*,'SZCC_Cash',sum


    CASE(2) ! n(z,q)

        call deltaN_yz(Z,Nz,LOGY,Ny,DNZ,DN,skyfracs,thetas,ylims,sz_switch,qa_ytot,erf_list)

        sum=0.
        if (print_counts==1) print*,'predicted counts'
        DO I=1,Nz
            if (print_counts==1) print*,I,DN(I,:)
            do J=1,Ny+1
                sum=sum+DN(I,J)
            enddo
        ENDDO
        if (print_counts==1) print*,'total counts',sum


        sum=0.
        do i=1,Nz
            do j=1,Ny+1
                if (DN(i,j) /= 0.) then
                    ln_factorial=0.
                    if (DNcat(i,j)/=0.) then
                        if (DNcat(i,j) >10.) then
                            ln_factorial=0.918939+(DNcat(i,j)+0.5)*dlog(DNcat(i,j))-DNcat(i,j)
                            !Stirling approximation only for more than 10 elements
                        else
                            !direct computation of factorial
                            fact=1.
                            do ii=1,int(DNcat(i,j))
                                fact=ii*fact
                            enddo
                            ln_factorial=dlog(fact)
                        endif
                    endif
                    sum=sum-1.*(DNcat(i,j)*dlog(DN(i,j))-DN(i,j)-ln_factorial)
                end if
            enddo
        end do
        SZCC_Cash=sum

    END SELECT

    !   Print*,'SZ lnlike = ',SZCC_Cash

    !priors for SZ nuisance params

    SZCC_Cash = SZCC_Cash + pystar*(cosmopar%logystar-(-0.186))**2./(2.*0.021**2.) + &
        palpha*(cosmopar%alpha-1.789)**2./(2.*0.084**2.) + psigma*(cosmopar%sigmaM-0.075)**2./(2.*0.01**2.) !cosmomc_sz
    !print*,'prior on ystar, alpha, scatter'

    SZCC_Cash = SZCC_Cash +pbeta*(cosmopar%beta-0.6666666)**2/(2.*0.5**2.) !prior beta evolution

    ! PRIORS ON NS AND OMEGABH2:
    SZCC_Cash = SZCC_Cash + pns*(CMB%InitPower(ns_index)-0.9624)**2./(2.*0.014 **2.) + &
        oh2*(CMB%ombh2 - 0.022)**2/(2*0.002**2.)  !BBN

    !PRIORS ON THE BIAS:
    SZCC_Cash = SZCC_Cash + lens*(cosmopar%biasinv-0.99)**2./(2.*0.19**2.) !final
    !(1-b)^(-1)=0.99 pm 0.19. !Planck CMB lensing

    !SZCC_Cash = SZCC_Cash + clash*(cosmopar%bias-0.8)**2/(2.*0.11**2)!not final
    !1-b=0.81+-0.08. !clash
    SZCC_Cash = SZCC_Cash + wtg*(cosmopar%bias-0.688)**2./(2.*0.072**2.)!final

    SZCC_Cash = SZCC_Cash + cccp*(cosmopar%bias-0.780)**2./(2.*0.092**2.)!final


    end function SZCC_Cash


    end module szcounts

