 !Module for using Viel et al. Lyman-alpha data
!J.Lesgourgues,29/10/04
!AL modification to integrate over calibration numerically,
!   reject models with weird high om_m that probe v small scales
!April 2006: Fixed major bug (introduced by AL)
!            removed unsafe OPENMP; CROFT on by default

module lya
use settings
use cmbtypes
implicit none

 integer, parameter :: lya_points = 9
 real :: lya_k(lya_points)

 ! for LUQAS:
 real :: lya_P1(lya_points), lya_dP1(lya_points)
 logical :: use_LUQAS = .true.

 ! for CROFT:
 real :: lya_P2(lya_points), lya_dP2(lya_points)
 logical :: use_CROFT = .true.

 real :: lya_kmax = 6.
 logical :: do_lya_init  = .true.
 logical :: Use_lya = .false.

contains

 subroutine lya_init
   integer i

    if (Feedback > 0) write(*,*) 'reading: Lyman-alpha data'
    call OpenTxtFile('data/lyman_alpha.dat', tmp_file_unit)

    do i =1, lya_points
      read (tmp_file_unit,*, end = 200, err=200) &
           lya_k(i), lya_P1(i), lya_dP1(i), lya_P2(i), lya_dP2(i)
    end do
    close(tmp_file_unit)

    do_lya_init = .false.

    goto 300

200 stop 'Error reading Lyman-alpha file'

300 return

 end subroutine lya_init

 function LSS_lyalike(CMB, Theory)
    Type (CMBParams) CMB
   Type (CosmoTheory) Theory
   real LSS_lyalike
   real  omegam, g2, omegam_z, omegav_z, th, chisq, minchisq
   real z, z1, z2, coef, D_ratio, acal
   integer, parameter :: ncal=12
   real dif2(-ncal:ncal)
   real, parameter :: dcal= 0.25
   real calweights(-ncal:ncal)
   real, parameter :: cal_sigma = 0.29
   integer i, ical

    if (do_lya_init) call lya_init

    if (CMB%W /= -1. .or. CMB%omnu/=0.) then
      write (*,*) 'Lya.f90 not tested for extended models'
      stop
    end if


    z1=2.125
    z2=2.72
    omegam=1.-CMB%omv-CMB%omk

    dif2 = 0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! compute chi2 for LUQAS

    if (use_LUQAS) then

       z=z1

    ! coefficient from h/Mpc to s/km
       g2=omegam*(1.+z)**3+CMB%omk*(1.+z)**2+CMB%omv
       coef = 100. * sqrt(g2) / (1.+z)
    ! growth factor
       omegam_z= omegam*(1.+z)**3/g2
       omegav_z= CMB%omv/g2
       D_ratio = omegam_z/(1.+z)/(exp(4./7.*log(omegam_z))-omegav_z &
         +(1.+omegam_z/2.)*(1.+omegav_z/70.)) &
         /omegam*(exp(4./7.*log(omegam))-CMB%omv &
         +(1.+omegam/2.)*(1.+CMB%omv/70.))


      if (lya_k(lya_points)*coef > &
        exp(log(matter_power_minkh) + (num_matter_power-1)*matter_power_dlnkh)) then
         !Just thow out if way-off model
          LSS_lyalike = LogZero
          return
       end if

       do i=1, lya_points

         th = MatterPowerAt(Theory,lya_k(i)*coef) *coef**3 &
               *(D_ratio)**2 &
               *(1.+1.4*exp(0.6*log(omegam_z)))/2.4

         do ical=-ncal,ncal
           acal = 1+ical*cal_sigma*dcal

           dif2(ical) = dif2(ical) +  ((th-lya_P1(i)*0.93*acal)**2/ &
                      (lya_dP1(i)*0.93*acal)**2)
         end do

       end do

    end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! compute chi2 for z=z1

    if (use_CROFT) then

       z=z2

    ! coefficient from h/Mpc to s/km
       g2=omegam*(1.+z)**3+CMB%omk*(1.+z)**2+CMB%omv
       coef = 100. * sqrt(g2) / (1.+z)
    ! growth factor
       omegam_z= omegam*(1.+z)**3/g2
       omegav_z= CMB%omv/g2
       D_ratio = omegam_z/(1.+z)/(exp(4./7.*log(omegam_z))-omegav_z &
         +(1.+omegam_z/2.)*(1.+omegav_z/70.)) &
         /omegam*(exp(4./7.*log(omegam))-CMB%omv &
         +(1.+omegam/2.)*(1.+CMB%omv/70.))


       do i=1, lya_points

         th = MatterPowerAt(Theory,lya_k(i)*coef) *coef**3 &
               *(D_ratio)**2 &
               *(1.+1.4*exp(0.6*log(omegam_z)))/2.4

         do ical=-ncal,ncal
           acal = 1+ical*cal_sigma*dcal

           dif2(ical) = dif2(ical) +  ((th-lya_P2(i)*0.93*acal)**2/ &
                      (lya_dP2(i)*0.93*acal)**2)
         end do

       end do

    end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    do i=-ncal,ncal
      acal = 1+i*cal_sigma*dcal
      calweights(i) = exp(-(1-acal)**2/cal_sigma**2/2)
    end do

    minchisq = minval(dif2)
    chisq = sum(exp(-(dif2-minchisq)/2)*calweights)/sum(calweights)

    if (chisq == 0) then
     chisq = 2*LogZero
    else
     chisq =  -2*log(chisq) + minchisq
    end if

    if (Feedback>1) write(*,*) 'Lyman-alpha chi-sq: ', chisq

    LSS_lyalike = chisq/2.

 end function LSS_lyalike

end module lya
