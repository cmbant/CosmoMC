!-----------------------------------------------------------!
!Module for SZ cluster counts likelihood
!
! Module prepared for cosmomcplanck by A. Bonaldi 2014       
! 
! Original code to compute cluster counts from theory written from 2001  
!                                                         !
!                    v 0.1                                  !
!           Jochen Weller and Richard Battye                !
!                     4/2/2007                              !
! all the necessary modules and subroutines                 !
!                                                           !
!                   v 1.0                                   !
!               Jochen Weller                               !
!                 10/1/2008                                 !
!              modified for optical cluster counts          ! 
!                                                           !
!                   v 2.0                                   !
!               Jochen Weller                               !
!                 30/12/2009                                !
!              modified for SZ cluster counts               ! 
!
!                   v 3.0 
!                Anna Bonaldi 2011-2013
!      - Include realistic selection function for Planck
!      - missing redsifts, errors on redshifts
!               Matthieu Roman 2012-2013
!      - s/n mass proxy and QA completeness
!-----------------------------------------------------------!


module constants_sz
  use PRECISION
  implicit none
  real(dl), PARAMETER :: Mpc = 3.08568025d22 ! Mpc in metres
  real(dl), PARAMETER :: G = 6.67300d-11 ! Newton's constant in m^3kg^-1s^-2
  real(dl), PARAMETER :: c = 3.0d8 ! speed of light in m/s
  real(dl), PARAMETER :: msun = 1.98892d30 ! mass of sun in kg
  real(dl), PARAMETER :: rhocrit0=2.7751973751261264d11 ! rhocrit in units of h^-1 Msun/ h^-3 Mpc^3

end module constants_sz


! ####################################################################


module cosmology
  USE precision
  USE CONSTANTS_sz
  implicit none
! note that we switch to the generalized Linder parametrization now
  public
  TYPE cospar
     REAL(dl) :: H0,w0,w1,omegam,omegav,n,sig8,omegak,omegabh2,gamma,ystar,alpha,sigmaM,bias,logystar,sigmaR(1000,2)
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
  USE PRECISION
  USE SURVEY
  USE COSMOLOGY
  implicit none
real(dl),parameter :: thetastar=6.997 ! dic 2012 new scallings
real(dl),parameter :: alpha_theta=1./3.

contains

  function theta500(m,z)
    real(dl), INTENT(IN) :: z,m
    real(dl) :: theta500,thetastar2,m2
    m2=m*cosmopar%bias !hydro bias
    thetastar2=thetastar*(cosmopar%H0/70.)**(-2./3.)
    theta500=thetastar2*(m2/3.e14*(100/cosmopar%H0))**alpha_theta*Eh(z)**(-2./3.)*(100.0*da(z)/500.0/cosmopar%H0)**(-1.)
    RETURN
  end function theta500
  


  function y500(m,z)
    real(dl), INTENT(IN) :: z,m
    real(dl) :: y500,ystar2,m2,alpha
    m2=m*cosmopar%bias !hydro bias
    alpha=cosmopar%alpha
    ystar2=cosmopar%ystar
    ystar2=ystar2*(cosmopar%H0/70.)**(-2.+alpha)
    y500=ystar2*(m2/3.0e14*(100./cosmopar%H0))**alpha*Eh(z)**(2./3.)*(100.0*da(z)/500.0/cosmopar%H0)**(-2.)
    RETURN
  end function y500

  function erf_compl(y,sn,q)
    REAL(dl)::y,sn,arg,erf_compl,q   
    arg=(y-q*sn)/sqrt(2.)/sn
    erf_compl=(erf(arg)+1.)/2.
  end function erf_compl


end module massobservable

! ####################################################################    
! This module links to libcamb, initializes it and provides a 
! subroutine to obtain sigma8

module power
  use PRECISION
  use CAMB
  use TRANSFER
  use cosmology
!use cmbmain
!use CosmologyTypes
!use CosmoTheory
!use Calculator_Cosmology
!use Likelihood_Cosmology

use Calculator_Cosmology
  implicit none
  public
  REAL(dl) :: normsig8
  REAL(dl) :: normgrowth
  type(CAMBparams) :: P2 
  contains




  function sigma(R)
! radius in units of h^-1 Mpc
    real(dl), INTENT(IN) :: R
    real(dl) :: sigma,max_r,min_r,s1,s2,r1,r2
    INTEGER::p(1),p2
  !  Type(MatterTransferData) :: MT
 !   Class(TCosmoTheoryPredictions), target :: Theory

   !if (global_error_flag/=0) return


!!$    call Transfer_Get_sigma8(MT, R*1.0_dl)
!!$    sigma = MT%sigma_8(1,1)
!!$
!!$    if (ISNAN(sigma)==.true.) then
!!$       print*,'Sigma is NaN!!!'
!!$
!!$    endif
!!$print*,'old',sigma
    p=minloc(abs(cosmopar%sigmaR(:,1)-R))
    max_R=maxval(cosmopar%sigmaR(:,1))
    min_R=minval(cosmopar%sigmaR(:,1))

    if (R > max_R) then
       sigma=cosmopar%sigmaR(1000,2)
    else if  (R < min_R) then
       sigma=cosmopar%sigmaR(1,2)
    else
       p=minloc(abs(cosmopar%sigmaR(:,1)-R))
       s1=cosmopar%sigmaR(p(1),2)
       r1=cosmopar%sigmaR(p(1),1)
       p2=p(1)+1
       if (R<r1) p2=p(1)-1 
       s2=cosmopar%sigmaR(p2,2)
       r2=cosmopar%sigmaR(p2,2)
       sigma=s1+(s2-s1)/(r2-r1)*(R-r1)
    endif


!!$print*,'new',sigma,R,s1,s2
!!$    stop
    RETURN
  end function sigma
  
  function dsigdR(R)
! Radius in units of h^-1 Mpc
    real(dl), INTENT(IN) :: R
    real(dl) :: dsigdR
    real(dl), PARAMETER :: deps = 0.01
    real(dl) :: s1,s2
    s1 = sigma(R*(1-deps))
    s2 = sigma(R*(1+deps)) 
    dsigdR=(s2-s1)/deps/R/2.0
    RETURN
  end function dsigdR

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



end module power

! ####################################################################

module massfunction
  USE PRECISION
  USE CONSTANTS_sz
  USE POWER
  USE massobservable

  implicit none
  TYPE MASSPAR
     REAL(dl) :: Amf,Bmf,epsmf,dso,pmf
     INTEGER :: psind
  END TYPE MASSPAR
  TYPE (masspar) :: massfnpar
contains


  function dndlnM_new(z,M,g)
! mass in units of h^-1 M_sun
    real(dl), INTENT(in) :: z,M,g
    real(dl) :: dndlnM_new
    real(dl) :: R,rhom0,rhom
    real(dl) :: dMdR,sR,fJen
    real(dl) :: fTink
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
! of run
!!$    if (1._dl/R > P%Transfer%kmax) then
!!$       write(*,*) 'Cluster radius too small, increase P%Transfer%kmas from',P%Transfer%kmax,' to ',1._dl/R,' or possibly larger !'
!!$       stop
!!$    end if
! dM/dR
    dMdR = 3*M/R

    sR =sigma(R)

!    write(*,*) cosmopar%omegam,sR
    !g = delta(z)
! Parameters from Jenkins et al. for LCDM (Appendix) and SO(324)
! Amf = 0.316 ; Bmf = 0.67 ; epsmf = 3.82    
! Note SO(324) is mass in mass density units or 0.3*324=97.2

    dsoz=massfnpar%dso/Omegam(z)
    if (massfnpar%psind==1) then
       fJen = massfnpar%Amf*exp(-abs(-dlog(g*sR)+massfnpar%Bmf)**massfnpar%epsmf)
       dndlnM_new = -rhom0*fJen*dsigdR(R)/dMdR/sR
    elseif (massfnpar%psind==2) then
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
       dndlnM_new = -rhom0*fTink*dsigdR(R)/dMdR/sR
    else
       write(*,*) 'Invalid mass function indicator: ',massfnpar%psind
       stop
    end if
    RETURN
  end function dndlnM_new


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
    IF (H.EQ.0.) PAUSE 'Bad XA input.'
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
  USE POWER
  USE MASSOBSERVABLE
  USE MASSFUNCTION
  IMPLICIT NONE

contains

  function next_z(zi,binz)
    real(dl)::zi,binz,next_z,dzi,hr

    
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



  SUBROUTINE deltaN_yz(Z,Nz,LOGY,Ny,DN,skyfracs,thetas,ylims,switch)
       REAL(dl),INTENT(IN) :: z(:),logy(:),thetas(:),skyfracs(:),ylims(:,:)
       INTEGER, INTENT(IN):: Ny,Nz,switch
       REAl(dl),allocatable :: grid(:,:),steps_z(:),steps_z2(:),steps_m(:)
       REAl(dl),allocatable ::ytheta5_data(:,:),ytheta5_model(:),Mz_data(:,:)
       REAl(sp),allocatable ::completeness(:,:)
       REAl(dl),allocatable ::dif1(:),dif2(:,:),mlims(:,:)
       REAL(dl)::  DN(:)
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
      

       lnmmin=29.9  !1.e13 Msun/h-1
       lnmmax=38.5  !5.e16 Msun/h-1 
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



       allocate(grid(nsteps_m,nsteps_z),completeness(nsteps_m,nsteps_z),stat=iostat)
       if (iostat/=0) then
          print*,'allocation error'
          stop
       endif

       completeness(:,:)=0.       
       call grid_C(npatches,steps_z,steps_m,thetas,ylims,skyfracs,completeness)


       grid(:,:)=0.
       call get_grid(nsteps_z,nsteps_m,steps_z,steps_m,grid)

       
       DN(:)=0.
       ! from first bin in z to Nz                                                              
       SELECT CASE(SWITCH)
       CASE(1)
          DO I=1,Nz
             !limits of bin in z                                                                               
             z1=Z(I)-0.5_dl*binz
             z2=Z(I)+0.5_dl*binz

             call integrate_m_z(grid,skyfracs,completeness,steps_z,steps_m,nsteps_z,nsteps_m,z1,z2,binz,dlnm,sum)
             DN(I)=sum
          ENDDO
          
       CASE(2)
          DO I=1,Ny
             !limits of bin in z                                              
             
             y1=LOGY(I)-0.5_dl*biny
             y2=LOGY(I)+0.5_dl*biny
             y1=10.**y1
             y2=10.**y2
             call integrate_m_y(grid,skyfracs,completeness,steps_z,steps_m,nsteps_z,nsteps_m,y1,y2,biny,dlnm,sum)
             DN(I)=sum

          ENDDO
       END SELECT

       deallocate(grid,steps_m,steps_z,completeness,stat=iostat)
       if (iostat/=0) then
          print*,'deallocation error'
          stop      
       endif
       RETURN
     end SUBROUTINE DeltaN_yz

  SUBROUTINE integrate_m_z(grid,skyfracs,compl,steps_z,steps_m,nsteps_z,nsteps_m,z1,z2,binz,dlnm,sum)

    REAl(dl),intent(in) ::grid(:,:),skyfracs(:)
    REAL (sp),intent(in)::compl(:,:)
    REAl(dl),intent(in) ::steps_z(:),steps_m(:),z1,z2,binz,dlnm
    REAl(dl):: sum,xi,xi1,zi,sigmaM,sq_sigmaM,lnm,window,bias,sum2,dz,sum3
    REAl(dl)::dnnew,dnold,lny1_new,lnylim,angsize,counts
    REAL (dl)::dif,dif_t,frac,ylim,frac_same,limit_m,test_new,h,C
    REAl(dl):: test(nsteps_z),c1,c2,f1,f2,x1,x2
    INTEGER::nsteps_z,nsteps_m,nsum,nthetas,col,k_ini,k_fin,t,k_ininew,col_old,npatches
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


!       do k=1,npatches
          sum2=0. 
!          frac=skyfracs(k) 
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
          sum=sum+sum2!*frac
!       enddo


  END SUBROUTINE integrate_m_z

  SUBROUTINE integrate_m_y(grid,skyfracs,compl,steps_z,steps_m,nsteps_z,nsteps_m,y1,y2,biny,dlnm,sum)

    REAl(dl),intent(in) ::grid(:,:),skyfracs(:)
    REAL (sp),intent(in)::compl(:,:)
    REAl(dl),intent(in) ::steps_z(:),steps_m(:),y1,y2,biny,dlnm
    REAl(dl):: sum,xi,xi1,zi,sigmaM,sq_sigmaM,lnm,window,bias,sum2,dz,sum3,mp,zp,yp
    REAl(dl)::dnnew,dnold,lny1_new,lnylim,angsize,counts,dy
    REAL (dl)::dif,dif_t,frac,ylim,frac_same,limit_m,test_new,h,C
    REAl(dl):: test(nsteps_z),c1,c2,f1,f2,x1,x2
    INTEGER::nsteps_z,nsteps_m,nsum,nthetas,col,k_ini,k_fin,t,k_ininew,col_old,npatches,n
    INTEGER::k,kk,JJ,II,j1,j2,i,index
    INTEGER::p(1)


    sum=0._dl
    npatches=size(skyfracs)


    sum=0
    n=0


    do jj=1,nsteps_z-1                   
       zp=steps_z(jj)   
       DO ii=1,nsteps_m
          mp=exp(steps_m(ii))
          n=0
          yp=y500(mp,zp)
          if ((yp >= y1) .and. (yp < y2)) then
             n=n+1
             x1=steps_z(jj)
             x2=steps_z(jj+1)
             f1=grid(ii,jj)
             f2=grid(ii,jj+1)
             sum2=0. 
             c1=compl(ii,jj)
             c2=compl(ii,jj+1)
             sum2=sum2+0.5*(f1*c1+f2*c2)*(x2-x1)*dlnm
             sum=sum+sum2
          endif

       enddo

    enddo

  END SUBROUTINE integrate_m_y



  SUBROUTINE grid_C(npatches,steps_z,steps_m,thetas,ylims,skyfracs,completeness)
    REAl(dl),intent(in) :: steps_z(:),steps_m(:),thetas(:),ylims(:,:),skyfracs(:)
    INTEGER :: nsteps_z,nsteps_m,nthetas,ntab,npatches,iostat
    real(sp):: completeness(:,:)
    integer ::i,j,ii,jj,index1,index2,count,P(1),k1,k2,l1,l2,k,N,nthetas2
    real (dl):: dif_old,dif,max,min,dlm,binz,m_min,m_max,mp,yp,zp,thp
    real(dl),allocatable:: dif_y(:),dif_theta(:) 
    real(dl):: min_thetas,max_thetas,min_y,max_y
    real(dl):: c1,c2,th1,th2,c,y12,y1,y2,y
    REAL(dl),parameter::q=7.
    real(dl),allocatable::erfs(:,:)!thetas2(:),ylims2(:,:)
    real(dl)::win0,win,arg,arg0,y0,py,lnymax,lnymin,lny,dy,fac,mu,int,dlny,fsky
    real(dl):: sigmaM



    sigmaM=cosmopar%sigmaM

    nsteps_z=size(steps_z)
    nsteps_m=size(steps_m)
    nthetas=size(thetas)
    ntab=size(ylims)

    allocate(dif_y(ntab),dif_theta(nthetas),stat=iostat)
    if (iostat/=0) then
       print*,'allocation error'
    endif

    min_thetas=minval(thetas)
    max_thetas=maxval(thetas)


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
             completeness(ii,jj)=0.
             do i=1,npatches
                y1=ylims(i,l1)
                y2=ylims(i,l2)
                y=y1+(y2-y1)/(th2-th1)*(thp-th1)
                c2=erf_compl(yp,y,q) 
                completeness(ii,jj)=completeness(ii,jj)+c2*skyfracs(i)
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
             completeness(ii,jj)=int
          enddo
       enddo

       deallocate(erfs,stat=iostat)
       if (iostat/=0) then
          print*,'deallocation error'
          stop   
       endif
    endif
    deallocate(dif_y,dif_theta,stat=iostat)
    if (iostat/=0) then
       print*,'deallocation error'
       stop  
    endif

  END SUBROUTINE grid_C


  SUBROUTINE get_grid(nz,nm,z,lnm,grid)
    REAL(dl):: z(:),lnm(:),grid(:,:)
    INTEGER :: I,J,nz,nm
    REAL(dl):: Mnew,g,ynew,vol,theta
    


    DO I=1,Nz
       g = delta(z(I))
       vol=dVdzdO(z(I))
       DO J=1,Nm
          Mnew=exp(lnm(J))
          grid(j,i) = dndlnM_new(z(I),Mnew,g)*surveypar%deg2*vol
       ENDDO
    ENDDO

  END SUBROUTINE get_grid

end module numbercounts





! ####################################################################


module input
  USE PRECISION
  use cosmology
  use survey
  use massfunction
  use massobservable
  implicit none
  REAL(dl), SAVE :: zmaxs,z0,dz
  INTEGER, SAVE :: Nred, Nscat
  CHARACTER*80, SAVE :: outroot
contains


  function rows_number(filename,iunit)
    integer:: reason,iunit,nrows,rows_number
    character (LEN=200)::filename
    real(DL):: col
    open (UNIT=iunit,file=filename,status='unknown',form="formatted")
    print*,filename
    nrows=1
    DO
       READ(iunit,*,IOSTAT=Reason)  col
       IF (Reason > 0)  THEN 
          print*,'Error in reading file'
          print*,filename
          stop
       ELSE IF (Reason < 0) THEN
          exit
       ELSE
          nrows=nrows+1
       END IF
    END DO
   rows_number=nrows-1
    CLOSE (unit=iunit)
    return
  end function rows_number
  
 subroutine read_values(filename,iunit,nrows,vector)
    integer:: reason,iunit,nrows,i
    character (LEN=200)::filename
    real(DL)::col,vector(:)
    open (UNIT=iunit,file=filename,status='unknown',form="formatted")
    DO i=1,nrows
       READ(iunit,*,IOSTAT=Reason)  col
       vector(i)=col
    END DO
    CLOSE (unit=iunit)
    return
  end subroutine read_values

            
end module input
        
! ####################################################################



module szcounts
use MatrixUtils
use input
use settings
use CosmologyTypes
use CosmoTheory
use Calculator_Cosmology
use Likelihood_Cosmology
use IniObjects


implicit none
 private
 logical :: do_sz_init = .true.
 logical :: use_sz = .false.
 integer :: SZ_num 
 integer :: Ny,Nz
 integer :: SZ_switch  !1 fot N(z), 2 for N(y), 3 for N(z,y) 
 integer :: nmiss_switch  !0 for simple rescaling, 1 for MCMC 
 integer :: errorz_switch  !0 for simple rescaling, 1 for MCMC 
 integer :: nmiss,nred2,ncat
 integer :: ylim_switch  !>1 for constant ylim, <1 for variable ylim
 real(dl), allocatable :: DNcat(:,:),DNzcat(:),DNycat(:),DN(:,:),DNz(:),DNy(:)
 real(dl), allocatable :: Z(:),LOGY(:),ylims(:,:),thetas(:),skyfracs(:)
 real(dl), allocatable :: SZcat(:,:)
 real :: sz_kmax = 4.0
 integer :: sz_k_per_logint = 5
 character(len=256) :: SZ_filename = ''
 character(LEN=*), parameter :: SZ_version =  'June_2014'

 type, extends(TCosmoCalcLikelihood) :: SZLikelihood
  contains
 procedure :: LogLike => SZCC_Cash
 end type SZLikelihood

 !PRIVATE:: ran_mwc,get_free_lun,randgauss_boxmuller
 PUBLIC :: SZLikelihood_Add, SZcc_Cash

contains

 
  subroutine SZLikelihood_Add(LikeList, Ini)
    implicit none
    class(TLikelihoodList) :: LikeList
    class(TSettingIni) :: ini
    integer count
    Type(SZLikelihood), pointer :: this
   
    if (.not. Ini%Read_Logical('use_SZ',.false.)) return

    if (Ini%Read_Logical('use_SZ',.false.)) then
       allocate(this)
       this%needs_background_functions = .true.
       this%version = SZ_version
       call this%loadParamNames(trim(DataDir)//'SZ.paramnames')
       call LikeList%Add(this)
       this%LikelihoodType = 'SZ'
       this%name='SZ'
       !this%power_kmax=sz_kmax
       CALL SZ_init
    endif
     end subroutine SZLikelihood_Add

     subroutine SZ_init
!!$       use settings
!!$       use input
       ! reading input file containing data and selection function
       ! computing the data counts
       ! basic settings
       implicit none
       character (LEN=20):: name
       character (LEN=200)::cat_filename,skyfracs_filename,ylims_filename,thetas_filename,filename
       real :: dummy
       real :: dzorg
       integer :: i,j,iostat,nrows,reason,ii,iunit,nrows_old,nthetas,npatches
       integer ::Nfact,error
       real(dl) :: ymin,ymax,dlogy,logymax,logymin,yi,col1,col2,sum,col3,col4
       real(dl) :: y_min,y_max,z_min,z_max,factorial

       !real(dl),allocatable :: SZcat(:,:)

       print*,'Initialising SZ likelihood'
       error=0
       !file names 
       cat_filename='data/SZ_cat.txt'
       thetas_filename='data/SZ_thetas.txt'
       skyfracs_filename='data/SZ_skyfracs.txt'
       ylims_filename='data/SZ_ylims.txt'

       massfnpar%dso=500.
       massfnpar%psind=2 ! Fix to Tinker et al. mass function

       !dz=0.05
       surveypar%ab=0.
       surveypar%nb=0.
       surveypar%sfid=0.
       Nscat=0
       surveypar%deg2=41253.0  !full sky (sky fraction handled in skyfracs file)

       z0=0.
       zmaxs=1.
       dz=0.1

       surveypar%ylimin=1.e-3
       surveypar%ymaxin=1.e-1
       dlogy=0.2


       ymin=surveypar%ylimin!1.e-3
       ymax=surveypar%ymaxin!1.e-1
       logymin=dlog10(ymin)
       logymax=dlog10(ymax)
       sz_switch=1 !1=n(z); 2=N(y); 3=N(z,y)

       ylim_switch=-1 !>0=constant ylim, <0=variable ylim
       nmiss_switch=0
       !nmiss_switch=1
       errorz_switch=0
       !    errorz_switch=1

       surveypar%deg2 = 3.046174198d-4*surveypar%deg2 ! convert to ... units
       Nz = DINT((zmaxs-z0)/dz)+1 ! for some reason here it approximates by defect, at difference with scalar version
       Ny = DINT((logymax-logymin)/dlogy)+1
       !    print*,'Ny=',Ny

       !ylims file
       if (ylim_switch < 0) then

          ymin=-1.*ymin  

          ! read mock catalogue
          filename=cat_filename

          nrows=0

          print*,'Reading catalogue'

          CALL get_free_lun(iunit)
          open (UNIT=iunit,file=filename,status='old',form="formatted",iostat=reason)
          IF (Reason > 0)  THEN 
             print*,'Error opening file',filename
          endif

          print*,filename
          nrows=1
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
             END IF
          END DO
          nrows=nrows-1
          CLOSE (unit=iunit)
          print*,'done'

          Print*,'Catalogue Number of clusters=',nrows
          ALLOCATE(SZcat(nrows,3),stat=iostat)! z,y
          if (iostat/=0) then
             print*,'allocation error'
          endif
          CLOSE (unit=iunit)

          CALL get_free_lun( iunit )
          open (UNIT=iunit,file=filename,status='old',form="formatted")
          DO i=1,nrows
             READ(iunit,*,IOSTAT=Reason)  col1, col2, col3
             !SZcat(i,1)=col1  !redshift
             SZcat(i,1)=col2
             SZcat(i,2)=col2  !error on redshift
             ! SZcat(i,3)=col3  !y flux?
             !     !print*,col1
          ENDDO
          CLOSE (unit=iunit)



          !file with theta
          filename=thetas_filename
          nthetas=rows_number(filename,1)
          print*,'Number of size thetas=',nthetas
          allocate(thetas(nthetas),stat=iostat) 
          if (iostat/=0) then
             print*,'allocation error'
          endif

          ! theta, nbins, first_index, last_index,first_index2, last_index2
          call read_values(filename,1,nthetas,thetas)

          !file with skyfracs
          filename=skyfracs_filename
          npatches=rows_number(filename,1)
          print*,'Number of patches=',npatches
          allocate(skyfracs(npatches),stat=iostat)
          if (iostat/=0) then
             print*,'allocation error'
          endif

          ! theta, nbins, first_index, last_index,first_index2, last_index2
          call read_values(filename,1,npatches,skyfracs)

          filename=ylims_filename
          nrows=rows_number(filename,1)
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

          CLOSE (unit=iunit)
          CALL get_free_lun( iunit )

          open (UNIT=iunit,file=filename,status='unknown',form="formatted")
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

       ALLOCATE(Z(Nz),LOGY(Ny),stat=iostat)
       if (iostat /= 0) then
          print *, "Cannot allocate work arrays"
          stop
       endif

       ! logy vector
       yi=logymin+dlogy/2.
       DO I=1,Ny
          logy(I)=yi
          yi=yi+dlogy
       END DO

       ! z vector
       DO I=0,Nz-1
          Z(I+1)=z0+I*dz+0.5_dl*dz
       END DO
       if (z0==0._dl) Z(1)=Z(1)+1.e-8 ! for numerical problem when starting from 0. 

       ALLOCATE(DNcat(Nz+1,Ny+1),DNzcat(Nz),DNycat(Ny),DN(Nz+1,Ny+1),DNz(Nz),DNy(Ny),stat=iostat)
       if (iostat /= 0) then
          print *, "Cannot allocate work arrays"
          stop
       endif

       DNcat(:,:)=0.
       DNzcat(:)=0.
       DNycat(:)=0.

       SELECT CASE(SZ_SWITCH)
       CASE(1)      ! N(z)
          nrows=size(SZcat(:,1))

          nmiss=0
          DO ii=1,nrows 
             if (SZcat(ii,1) <0.) nmiss=nmiss+1.
          enddo


          DO I=1,Nz
             z_min=Z(I)-dz/2.
             z_max=Z(I)+dz/2.
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

          if (dble(nred2) /= sum) then  
             print*,'Error number of clusters!'
             stop
          endif
          ncat=nrows
!!$         if (nred2+nmiss /=nrows) then
!!$            print*,'Error number of clusters!'
!!$            stop
!!$         endif

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
          DNycat(:)=0.
          DNycat(1) =     8.
          DNycat(2) =     49.
          DNycat(3) =     73.
          DNycat(4) =     63.
          DNycat(5) =     37.
          DNycat(6) =     16.
          DNycat(7) =     9.
          DNycat(8) =     13.
          DNycat(9) =     1.
          DNycat(10) =    1.
          DNycat(11) =    0.

          sum=0.
          DO I=1,Ny+1
             sum=sum+DNycat(I)
          ENDDO

          print*,'DN(y) number of clusters',sum

       END SELECT

       do_sz_init = .false.

       print*,'End SZ initialization'
      
     end subroutine SZ_init


 function SZCC_Cash(this,CMB,Theory,DataParams)
   ! SZ nuisance in dataparams 
   use camb
   use cosmology
   use numbercounts
   Class(SZLikelihood) :: this
   Class (CMBParams):: CMB
   Class(TCosmoTheoryPredictions), target :: Theory
   real(mcp) DataParams(:)
   real(mcp)  SZCC_Cash
   real(mcp) :: SZCC_Cash_ini
   INTEGER :: N,i,j,Nf,ii
   REAL(DL) :: sum,factorial,ln_factorial,SZCC_Cash_exp
   character*80 :: outfile
   INTEGER :: iseed,p(1),jj,nit,it,iostat
   REAL(DL) ::test,dif,difold,z_min,z_max
   REAL(DL),allocatable :: DNzcum(:),DNz_old(:),RANDCAT(:)

   save iseed
!!$   if (do_sz_init) then 
!!$      call SZ_init
!!$   end if

   SZCC_Cash=logzero

   nit=1000
   cosmopar%H0=CMB%H0
   cosmopar%w0=CMB%w
   cosmopar%w1=0.0 ! Note, we do not support evolution of w right now
   cosmopar%omegam=CMB%omc+CMB%omb+CMB%omnu !cosmopar%omegam=0.3
   cosmopar%omegav=CMB%omv!1.-cosmopar%omegam!  cosmopar%omegav=0.7
   cosmopar%omegak=CMB%omk
   cosmopar%n=CMB%InitPower(1)
   cosmopar%sig8=Theory%sigma_8!     
   cosmopar%omegabh2=CMB%ombh2
   cosmopar%gamma=-1

   !new mapping of nuisance parameters
   cosmopar%alpha=DataParams(1)
   cosmopar%ystar=10.**(DataParams(2))/(2.**cosmopar%alpha)*0.00472724
   cosmopar%logystar=(DataParams(2))!to set the prior on logy
   cosmopar%bias=DataParams(3)
   cosmopar%sigmaM=DataParams(4)
   cosmopar%sigmaR=Theory%sigma_R
!!$   print*,'ystar=',cosmopar%ystar
!!$   print*,'n=',cosmopar%n
!!$   print*,'alpha=',cosmopar%alpha
!!$   print*,'scatter=',cosmopar%sigmaM
!!$   print*,'bias=',cosmopar%bias
!!$   print*,'H0=',cosmopar%H0
!!$   print*,'omegam=',cosmopar%omegam
!!$   print*,'omegav=',cosmopar%omegav
!!$   print*,'omegak=',cosmopar%omegak
!!$   print*,'sigma8=',cosmopar%sig8

   call INIGROWTH
 
   DN(:,:)=0.
   DNz(:)=0.
   DNy(:)=0.

   SELECT CASE(SZ_SWITCH)

   CASE(1) ! n(Z)
      !theo counts

      call deltaN_yz(Z,Nz,LOGY,Ny,DNZ,skyfracs,thetas,ylims,sz_switch)
      ! N(z)

      sum=0.
      DO I=1,Nz
         sum=sum+DNz(I)
      ENDDO

      !print*,'theory counts',DNZ
      if (ISNAN(DNZ(1))==.true.) then
         print*,'NaN found in theory counts!'
         stop
      endif
      !cumulative for missing redshifts
      allocate(DNzcum(Nz),DNz_old(Nz),randcat(ncat),stat=iostat)
      if (iostat/=0) then
         print*,'allocation error'
      endif
      
      do jj=1,Nz
         sum=0
         do ii=1,jj
            sum=sum+DNz(ii)
         end do
         DNzcum(jj)=sum
      enddo

      !check likelihood 
      DNz_old=DNzcat
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

      DNzcat=DNz_old

      if (nmiss_switch==0 .and. errorz_switch==0) goto 2

      if (dexp(-sum) /= 0._dl) then
         DNzcat=DNz_old
         !print*,'Initial Loglike=',SZCC_Cash
         SZCC_Cash_ini=SZCC_Cash

        ! print*,'ITERATIONS:'!, dexp(-sum)    
         SZCC_Cash_exp=0._dl

         do it=1,nit
            DNzcat=DNz_old
            if (errorz_switch==1) then
            randcat(:)=-1.
            ! realization of catalogue with redshift
            do i=1,ncat
               if (SZcat(i,1) > 0.) then
                  randcat(i)=randgauss_boxmuller(iseed)*SZcat(i,2)+SZcat(i,1)
                  if (randcat(i) <0.) randcat(i)=1.e-3
               endif
            enddo
            DNzcat(:)=0.
            DO I=1,Nz
               z_min=Z(I)-dz/2.
               z_max=Z(I)+dz/2.
               DO ii=1,ncat 
                  if ((randcat(ii) >0).and.(randcat(ii) >= z_min) .and. (randcat(ii) < z_max)) then
                     DNzcat(I)=DNzcat(I)+1.
                  endif
               ENDDO
               
            END DO

            sum=0.
            DO I=1,Nz
               sum=sum+DNzcat(I)
            ENDDO
            !if (int(sum) /= nred2) then 
            if (abs(sum-dble(nred2)>1.e-5)) then 
               print*,'error number of clusters with redhift'

               print*,nred2,sum
               stop
            endif
         endif
            SELECT CASE(nmiss_switch)
            case(0)
               
               DNzcat=DNZcat/dble(nred2)*dble(ncat)
            case(1)
               
            do i=1,Nmiss
               test=ran_mwc(iseed)*dble(DNzcum(Nz))
               difold=dble(DNzcum(Nz))
               do jj=1,Nz
                  dif=Dnzcum(jj)-test
                  if ((dif>= 0.) .and. (dif < difold)) then
                     j=jj 
                     difold=dif
                  endif
               enddo

               DNzcat(J)=DNZcat(J)+1
            enddo

         END SELECT

         !print*,'with missing z:',DNZcat

            sum=0.
            DO I=1,Nz
               sum=sum+DNzcat(I)
            ENDDO
            if (abs(sum-dble(ncat)>1.e-5)) then 
               print*,'error total number of clusters'
               print*,sum,ncat
               stop
            endif

            ! likelihood computation
            sum=0.
            do i=1,Nz
               if (DNz(i) /= 0.) then
               ln_factorial=0.
                  if (DNzcat(i) /= 0.) ln_factorial=0.918939+(DNzcat(i)+0.5)*dlog(DNzcat(i))-DNzcat(i) !Stirling
      
                  sum=sum-1.*(DNzcat(i)*dlog(DNz(i))-DNz(i)-ln_factorial)
               end if
            end do
            SZCC_Cash_exp = SZCC_Cash_exp+exp(-sum)
         ENDDO
         SZCC_Cash_exp=SZCC_Cash_exp/dble(nit)
         if (SZCC_Cash_exp >0._dl) then  
            !print*,'AVERAGE'
            SZCC_Cash=-1.*dlog(SZCC_Cash_exp)
           print*,SZCC_Cash,SZCC_Cash_ini
         else
            !print*,'INITIAL'
            SZCC_Cash=SZCC_Cash_ini
         endif
         print*,SZCC_Cash,SZCC_Cash_ini
      endif
2   continue
      !print*,'SZ Cash likelihood=',SZCC_Cash
      deallocate(DNzcum,DNZ_old,randcat,stat=iostat)
      if (iostat/=0) then
        print*,'deallocation error'
        stop  
   endif
 
!stop
   CASE(2) ! n(y)

      call deltaN_yz(Z,Nz,LOGY,Ny,DNy,skyfracs,thetas,ylims,sz_switch)

      ! N(y)
      sum=0.
      do i=1,Ny
         if (DNycat(i) /= 0.) then
            ln_factorial=0.918939+(DNy(i)+0.5)*dlog(DNy(i))-DNy(i) !Stirling
            sum=sum-2.*(DNy(i)*dlog(DNycat(i))-DNycat(i)-ln_factorial)
         end if
      end do
      SZCC_Cash = sum

!!$   CASE(3)
!!$
!!$      ! N(z,y)
!!$      sum=0.
!!$      do j=1,Ny
!!$         do i=1,Nz
!!$            if (DNcat(i,j) /= 0.) then
!!$               sum=sum+(DNcat(i,j)-DN(i,j))**2/DNcat(i,j)
!!$            end if
!!$         end do
!!$      enddo
!!$      SZCC_LnLike = sum/2.0
!!$

   END SELECT
   !print*,'like=',SZCC_Cash

   !priors for SZ nuisance params
   SZCC_Cash = SZCC_Cash + (cosmopar%logystar-(-0.186))**2/(2.*0.021**2) + (cosmopar%alpha-1.789)**2/(2.*0.084**2) + (cosmopar%sigmaM-0.075)**2/(2.*0.01**2) !cosmomc_sz
   !print*,'prior on ystar, alpha, scatter'

   ! PRIORS for PLANCK:
   !SZCC_Cash = SZCC_Cash + (CMB%InitPower(1)-0.9624)**2/(2.*0.014 **2) + (CMB%ombh2 - 0.022)**2/(2*0.002**2)  !BBN
   !print*,'BBN prior on Ns and Omegabh2'
   
   !PRIORS for WMAP:
   !GetLogLikePost = GetLogLikePost + ns,omegab
    !print*,'like=',SZCC_Cash
!end priors
!stop
!   if (Feedback>2) Print*,'SZ lnlike = ',SZCC_Cash

   Print*,'SZ lnlike = ',SZCC_Cash

 end function SZCC_Cash


  !=======================================================================
  function ran_mwc(iseed)
    !=======================================================================
    !     This routine implements the Marsaglia "multiply with carry method"
    !     and adds a Marsaglia shift sequence to it.
    !     (cf. references at http://www.cs.hku.hk/)
    !     You are welcome to replace this with your preferred generator
    !     of uniform variates in the interval ]0,1[ (i.e. excluding 0 and 1)

    !     (Re-)initialise with setting iseed to a negative number.
    !     Note that iseed=0 gives the same sequence as iseed=-1
    !     After initialisation iseed becomes abs(iseed) (or 1 if it was 0).

    !     B. D. Wandelt, May 1999
    !=======================================================================
    implicit none
    integer, intent(inout):: iseed
    real(SP) :: ran_mwc

    integer:: i,iseedl,iseedu,mwc,combined
    integer,save :: upper,lower,shifter
    integer,parameter :: mask16=65535,mask30=2147483647
    real(SP),save :: small
    logical, save :: first=.true.
    
    if (first.or.iseed<=0) then
       if(iseed==0) iseed=-1
       iseed=abs(iseed)
       small=nearest(1.0_sp,-1.0_sp)/mask30

       ! Users often enter small seeds - I spread them out using the
       ! Marsaglia shifter a few times.
       shifter=iseed
       do i=1,9
          shifter=ieor(shifter,ishft(shifter,13))
          shifter=ieor(shifter,ishft(shifter,-17))
          shifter=ieor(shifter,ishft(shifter,5))
       enddo

       iseedu=ishft(shifter,-16)
       upper=ishft(iseedu+8765,16)+iseedu !This avoids the fixed points.
       iseedl=iand(shifter,mask16)
       lower=ishft(iseedl+4321,16)+iseedl !This avoids the fixed points.

       first=.false.
    endif

100 continue

    shifter=ieor(shifter,ishft(shifter,13))
    shifter=ieor(shifter,ishft(shifter,-17))
    shifter=ieor(shifter,ishft(shifter,5))
    
    upper=36969*iand(upper,mask16)+ishft(upper,-16)
    lower=18000*iand(lower,mask16)+ishft(lower,-16)

    mwc=ishft(upper,16)+iand(lower,mask16)

    combined=iand(mwc,mask30)+iand(shifter,mask30)

    ran_mwc=small*iand(combined,mask30)
    if(ran_mwc==0._sp) goto 100

    return
  end function ran_mwc


  !=======================================================================
  function randgauss_boxmuller(iseed)
    !=======================================================================
    !     Box-Muller method for converting uniform into Gaussian deviates 
    !=======================================================================
    integer, intent(inout) :: iseed
    real(SP) :: randgauss_boxmuller
    logical, save :: empty=.true.
    real(SP) :: fac,rsq,v1,v2
    real(SP),save :: gset

    if (empty .or. iseed < 0) then ! bug correction, EH, March 13, 2003
1      v1=2.*ran_mwc(iseed)-1.
       v2=2.*ran_mwc(iseed)-1.
       rsq=v1**2+v2**2
       if(rsq.ge.1.or.rsq.eq.0.) goto 1
       fac=sqrt(-2.*log(rsq)/rsq)
       gset=v1*fac
       randgauss_boxmuller=v2*fac
       empty=.false.
    else
       randgauss_boxmuller=gset
       empty=.true.
    endif
    return
  end function randgauss_boxmuller
 

  SUBROUTINE get_free_lun( lun )
    IMPLICIT NONE
    INTEGER, INTENT(out) :: lun

    INTEGER, SAVE :: last_lun = 19
    LOGICAL :: used

    lun = last_lun
    DO
        INQUIRE( unit=lun, opened=used )
        IF ( .NOT. used ) EXIT
        lun = lun + 1
    END DO
    last_lun = lun

    END SUBROUTINE get_free_lun

end module szcounts

