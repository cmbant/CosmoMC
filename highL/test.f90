!===========================================================
program test_likelihood

! test high ell chi2

use highell_options
use highell_likelihood

implicit none

real(8), dimension(:), allocatable :: cl_tt
character(LEN=128) :: filename
real(8)            :: like_tot
real(8)            :: amp_tsz,amp_ksz,xi,aps148,aps217,aps95,aps150,aps220,acib150,acib220,rps0,rps1,rps,rcib
real(8)            :: cas1,cas2,cae1,cae2,cal_1,cal_2,cal_3
integer            :: il, dummy 
!---------------------------------------------------

print *,""
print *,"High ell likelihood chi2 test"
print *,"==================================="

!---------------------------------------------------
! read in test Cls
!---------------------------------------------------
filename = 'data/v2a1s_best_lcdm_6000.txt'
write(*,*)"Reading in Cls from: ",trim(filename)
open(unit=557,file=filename,action='read',status='old')

allocate(cl_tt(2:tt_lmax))
cl_tt(2:tt_lmax)=0.d0

do il=2,tt_lmax_mc
   read(557,*)dummy,cl_tt(il)
enddo
close(557)

call highell_likelihood_init

amp_tsz = 1.0 
amp_ksz = 1.2
xi      = 0.22
aps148  = 9.3
aps217  = 75.0
aps95   = 7.54
aps150  = 8.7
aps220  = 76.3
acib150 = 5.7
acib220 = 61.3
rps0    = 0.90
rps1    = 0.71
rps     = 1.00
rcib    = 0.92
cas1    = 1.00
cas2    = 1.01
cae1    = 1.00
cae2    = 0.99
cal_1   = 1.00
cal_2   = 1.01
cal_3   = 1.01

call highell_likelihood_compute(cl_tt,amp_tsz,amp_ksz,xi,aps148,aps217,aps95,aps150,aps220,acib150,acib220,rps0,rps1,rps,rcib,cas1,cas2,cae1,cae2,cal_1,cal_2,cal_3,like_tot)
print *, "----------------------------------------" 
print *, 'Total High ell chi2  =', 2*like_tot, -like_tot
print *, "----------------------------------------"

end program test_likelihood
