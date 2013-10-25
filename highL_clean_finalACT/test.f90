!===========================================================
program test_likelihood

! test high ell chi2

use highell_options
use highell_likelihood

implicit none

real(8), dimension(:), allocatable :: cl_tt
character(LEN=128) :: filename
real(8)            :: like_tot
real(8)            :: amp_tsz,amp_ksz,xi,aps148,aps217,aps95,aps150,aps220,acib150,acib220,ncib,rps0,rps1,rps,rcib,ags,age
real(8)            :: cas1,cas2,cae1,cae2,cal_1,cal_2,cal_3
integer            :: il, dummy 
!---------------------------------------------------

print *,""
print *,"High ell likelihood chi2 test"
print *,"==================================="

!---------------------------------------------------
! read in test Cls
!---------------------------------------------------
filename = trim(data_dir)//'wmap7_act_lcdm_bestfit_lensedCls_6000.dat'
write(*,*)"Reading in Cls from: ",trim(filename)
open(unit=557,file=filename,action='read',status='old')

allocate(cl_tt(2:tt_lmax))
cl_tt(2:tt_lmax)=0.d0

do il=2,tt_lmax_mc
   read(557,*)dummy,cl_tt(il)
enddo
close(557)

call highell_likelihood_init

amp_tsz = 4.80d0 !tSZ 
amp_ksz = 2.10d0 !kSZ
xi      = 0.20d0 !tSZxCIB
aps148  = 10.0d0 !ACT PS 148 GHz
aps217  = 75.0d0 !ACT PS 218 GHz
aps95   = 8.00d0 !SPT PS 95 GHz
aps150  = 10.0d0 !SPT PS 150 GHz
aps220  = 75.0d0 !SPT PS 220 GHz
acib150 = 5.70d0 !CIB 150 GHz
acib220 = 61.3d0 !CIB 220 GHz
ncib    = 0.80d0 !CIB index 
rps0    = 0.90d0 !rPS 95/150
rps1    = 0.71d0 !rPS 95/220
rps     = 1.00d0 !rPS 150/220
rcib    = 0.92d0 !rCIB 150/220
ags     = 0.40d0 !ACT-S dust
age     = 0.80d0 !ACT-E dust
cas1    = 1.00d0 !ACT-S 148 GHz cal
cas2    = 1.01d0 !ACT-S 218 GHz cal
cae1    = 1.00d0 !ACT-E 148 GHz cal
cae2    = 0.99d0 !ACT-E 218 GHz cal
cal_1   = 1.00d0 !SPT 95 GHz cal
cal_2   = 1.01d0 !SPT 150 GHz cal
cal_3   = 1.01d0 !SPT 220 GHz cal

call highell_likelihood_compute(cl_tt,amp_tsz,amp_ksz,xi,aps148,aps217,aps95,aps150,aps220,acib150,acib220,ncib,rps0,rps1,rps,rcib,ags,age,cas1,cas2,cae1,cae2,cal_1,cal_2,cal_3,like_tot)
print *, "----------------------------------------" 
print *, 'Total High ell chi2  =', 2.d0*like_tot
print *, "----------------------------------------"
print *, 'Expected -2lnL = 776.357229802025' 
print *, "----------------------------------------"

end program test_likelihood
