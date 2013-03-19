! Module for using SDSS McDonald et al. (2004) Ly-alpha data

! NOTE:  Not valid in growth function calculation for dark matter equations
! of state away from -1;
! also, this does not use alpha (running) information from McDonald chi^2 table

! If using SDSS Ly-alpha data:  cite
!  P.~McDonald {\it et al.},
!  %``The Linear Theory Power Spectrum from the Lyman-alpha Forest in the Sloan
!  %Digital Sky Survey,''
!  arXiv:astro-ph/0407377.

! implemented for CosmoMC by Kevork Abazajian

! AL: note requires *linear* power spectrum input
! AL: April 2006: added assumption violation traps
! AL: modified for /data/ path consistency

module Lya
use settings
use cmbtypes
use Precision
implicit none

 !Note all units are in k/h here
 logical :: Use_Lya = .false.


 integer, parameter :: c_lya_points = 12

 logical, parameter :: use_sdsslyaP  = .true.

 integer, parameter :: n_SDSSLya = 41

 real(dl) :: SDSSLya_delta(n_SDSSLya)
 real(dl) :: SDSSLya_neff(n_SDSSLya)
 real(dl) :: SDSSLya_chi2(n_SDSSLya,n_SDSSLya)
 real(dl) :: SDSSLya_chi2a(n_SDSSLya,n_SDSSLya)

 real ::  kh_pts0(c_lya_points),pk_pts0(c_lya_points),err_pts0(c_lya_points),z_pts

 logical :: do_lya_init  = .true.
 real(dl) :: ommpass,omlpass
 real :: lya_kmax = 1.5
 real(dl) :: z_SDSSP,kh_SDSSP,minSDSSPchi
contains

  function dgrowth(zz)
    real(dl) dgrowth
    real(dl) zz

    dgrowth = (1.d0+zz)/((ommpass*(1.d0+zz)**3+omlpass)**1.5d0)

  end function dgrowth


 subroutine lya_init
   integer i,j
   real temp
   if (Feedback > 0) write(*,*) 'reading: SDSS Ly-alpha data'

    if(use_SDSSLyaP)then
      minSDSSPchi = 1.d30
      z_SDSSP = 3.d0
      kh_SDSSP = 0.009d0
      call OpenTxtFile(trim(DataDir)//'SDSSPlyachi2.txt', tmp_file_unit)
      do i=1,n_SDSSLya
         do j=1,n_SDSSLya
            read(tmp_file_unit,*)SDSSLya_delta(i),SDSSLya_neff(j),temp,SDSSLya_chi2(i,j)
            if(SDSSLya_chi2(i,j)>0.)minSDSSPchi = min(SDSSLya_chi2(i,j),minSDSSPchi)
         end do
      end do

      if (Feedback > 0) write(*,*) 'done reading: SDSS Ly-alpha data'

      call splie2(SDSSLya_delta,SDSSLya_neff,SDSSLya_chi2,n_SDSSLya,n_SDSSLya,SDSSLya_chi2a)
   end if
   if (Feedback > 0) write(*,*) 'done reading: SDSS Ly-alpha data'

   do_lya_init = .false.

 end subroutine lya_init


 function LSS_Lyalike(CMB, Theory)
   Type (CMBParams) CMB
   Type (TheoryPredictions) Theory
   real LSS_Lyalike

   real hubz,h,omm,oml
   real logdelta,normf
   real(dl) chi2
   real(dl) delta,neff

   real(mcp) D0,D2,khlya,pk0,lnkhlya1,lnkhlya2,lnpk1,lnpk2

   real(dl) rombint
   external rombint

   if (do_lya_init) then
      call lya_init
   endif

    if (CMB%W /= -1. .or. CMB%omnu/=0. .or. CMB%InitPower(3)/= 0.) then
      write (*,*) 'This SDSS Lya module does not support extended models'
      write (*,*) 'for extensions see http://www.slosar.com/aslosar/lya.html'
      stop
    end if


   omm = CMB%omc+CMB%omb+CMB%omnu
   oml = CMB%omv
   h   = CMB%H0/100.

   ommpass = omm
   omlpass = oml
   D0 = rombint(dgrowth,0.d0,10000.d0,1.d-6) !5 omm / 2 not needed since only want ratio
   LSS_Lyalike = 0.

   if(use_sdsslyaP)then  ! use pt. slope

      normf = 1./(19.739209)  ! 1/(2*pi**2)

      hubz= CMB%H0*sqrt(omm*(1.+z_SDSSP)**3+oml)

      D2 =  (hubz/CMB%H0*rombint(dgrowth,z_SDSSP,10000.d0,1.d-6)/D0)**2 ! Note: not valid for w.ne.-1

      khlya = kh_SDSSP * hubz/h/(1.+z_SDSSP)
      pk0 = MatterPowerAt(Theory,khlya)*D2
      lnkhlya1 = log(khlya)*0.99
      lnkhlya2 = log(khlya)*1.01
      lnpk1 = log(MatterPowerAt(Theory,exp(lnkhlya1))*D2)
      lnpk2 = log(MatterPowerAt(Theory,exp(lnkhlya2))*D2)
      neff  = (lnpk2-lnpk1)/(lnkhlya2-lnkhlya1)

      delta = (pk0*khlya**3*normf)

      neff  = (lnpk2-lnpk1)/(lnkhlya2-lnkhlya1)

      call splin2(SDSSLya_delta,SDSSLya_neff,SDSSLya_chi2,SDSSLya_chi2a,n_SDSSLya,n_SDSSLya &
           ,delta,neff,chi2)

      if(chi2<minSDSSPchi)chi2 = 2.d30

      LSS_Lyalike = LSS_Lyalike + chi2/2.
   end if



 end function LSS_lyalike


 SUBROUTINE splie2(x1a,x2a,ya,m,n,y2a)
   INTEGER m,n,NN
   real(dl) x1a(m),x2a(n),y2a(m,n),ya(m,n)
   PARAMETER (NN=100)
   INTEGER j,k
   real(dl) y2tmp(NN),ytmp(NN)
   do j=1,m
      do k=1,n
         ytmp(k)=ya(j,k)
      end do
      call spline(x2a,ytmp,n,1.d30,1.d30,y2tmp)
      do k=1,n
         y2a(j,k)=y2tmp(k)
      end do
   end do
 end SUBROUTINE splie2


 SUBROUTINE splin2(x1a,x2a,ya,y2a,m,n,x1,x2,y)
   INTEGER m,n,NN
   real(dl) x1,x2,y,x1a(m),x2a(n),y2a(m,n),ya(m,n)
   PARAMETER (NN=100)
   INTEGER j,k
   real(dl) y2tmp(NN),ytmp(NN),yytmp(NN)
   do j=1,m
      do k=1,n
         ytmp(k)=ya(j,k)
         y2tmp(k)=y2a(j,k)
      end do
      yytmp(j) = splintt(x2a,ytmp,y2tmp,n,x2)
   end do
   call spline(x1a,yytmp,m,1.d30,1.d30,y2tmp)
   y = splintt(x1a,yytmp,y2tmp,m,x1)
 END SUBROUTINE splin2


 FUNCTION splintt(xa,ya,y2a,n,x)
   !  USE nrtype; USE nrutil, ONLY : assert_eq,nrerror
   !  USE nr, ONLY: locate
   integer, INTENT(IN)  :: n
   real(dl), INTENT(IN) :: xa(n),ya(n),y2a(n)
   real(dl), INTENT(IN) :: x
   real(dl) :: splintt
   integer :: khi,klo
   real(dl) :: a,b,h,y

   klo=max(min(locate(xa,x),n-1),1)
   khi=klo+1
   h=xa(khi)-xa(klo)
   if (h == 0.0)then
      write(*,*) 'bad xa input in splintt'
      splintt = -99.
   else
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      !            /** linear interpolation outside boundary **/
      if((klo==1).and.(b<0.0))then
         y=a*ya(klo)+b*ya(khi)-b*h*h*(2.0*y2a(klo)+y2a(khi))/6.0
      else if((khi==n).and.(a<0.0))then
         y=a*ya(klo)+b*ya(khi)-a*h*h*(y2a(klo)+2.0*y2a(khi))/6.0
      else
         y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.0
      endif
      splintt = y
   endif
 END FUNCTION splintt

FUNCTION locate(xx,x)
  REAL(dl), DIMENSION(:), INTENT(IN) :: xx
  REAL(dl), INTENT(IN) :: x
  INTEGER :: locate
  INTEGER :: n,jl,jm,ju
  LOGICAL :: ascnd
  n=size(xx)
  ascnd = (xx(n) >= xx(1))
  jl=0
  ju=n+1
  do
     if (ju-jl <= 1) exit
     jm=(ju+jl)/2
     if (ascnd .eqv. (x >= xx(jm))) then
        jl=jm
     else
        ju=jm
     end if
  end do
  if (x == xx(1)) then
     locate=1
  else if (x == xx(n)) then
     locate=n-1
  else
     locate=jl
  end if
END FUNCTION locate

end module Lya






























