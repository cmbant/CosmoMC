!
! UNION Supernovae Ia dataset
!
!
! This module uses the Union compilation. Please cite "Kowalski et
! al. (The Supernova Cosmology Project), Ap.J., 2008.".  If the Union
! compilation is included in any other data sets that are distributed,
! please include this citation request there too.
!
! By A Slosar, heavily based on the original code by A Lewis, S Bridle
! and D Rapetti. E-mail: Anze Slosar (anze@berkeley.edu) for questions
! about the code and David Rubin (rubind@berkeley.edu) for questions
! regarding the dataset itself.
!
! Marginalizes anayltically over H_0 with flat prior.  (equivalent to
! marginalizing over M, absolute magnitude; see appendix F of cosmomc
! paper). Resultant log likelihood has arbitary origin and is
! numerically equal to -chi^2/2 value at the best-fit value.


module snovae
use cmbtypes
use MatrixUtils
implicit none

 integer, parameter :: SN_num = 307
 double precision, parameter :: Pi_num = 3.14159265359D0 
 double precision :: SN_z(SN_num), SN_moduli(SN_num)
 double precision :: SN_Ninv(SN_num,SN_Num)
 double precision :: SN_sumninv

 logical, parameter :: SN_marg = .True.
 logical, parameter :: SN_syscovmat = .True.  !! Use covariance matrix with or without systematics
 

contains


 subroutine SN_init
   use settings
   character (LEN=20):: name
   integer i
   real :: tmp_mat(sn_num, sn_num)

   if (Feedback > 0) write (*,*) 'Reading: supernovae data'
   call OpenTxtFile(trim(DataDir)//'sn_z_mu_dmu.txt',tmp_file_unit)
   do i=1,  sn_num
      read(tmp_file_unit, *) name, SN_z(i), SN_moduli(i)
   end do
   close(tmp_file_unit)

   if (SN_syscovmat) then
      call OpenTxtFile(trim(DataDir)//'sn_covmat_sys.txt',tmp_file_unit)
   else
      call OpenTxtFile(trim(DataDir)//'sn_covmat_nosys.txt',tmp_file_unit)
   end if

   do i=1, sn_num
      read (tmp_file_unit,*) tmp_mat (i,1:sn_num)
   end do
   
   close (tmp_file_unit)

   call Matrix_Inverse(tmp_mat)
   sn_ninv = DBLE (tmp_mat)

   SN_sumninv = SUM(sn_ninv)
   
  end subroutine SN_init

 function SN_LnLike(CMB)
   use camb
  !Assume this is called just after CAMB with the correct model  use camb
  implicit none
  type(CMBParams) CMB
  logical, save :: do_SN_init = .true.

  real SN_LnLike
  integer i
  double precision z, AT, BT
  real diffs(SN_num), chisq

  if (do_SN_init) then 
     call SN_init
     do_SN_init = .false.
  end if


!! This is actually seems to be faster without OMP
     do i=1, SN_num
        z= SN_z(i)
        diffs(i) = 5*log10((1+z)**2*AngularDiameterDistance(z))+25 -sn_moduli(i)
     end do
     
     AT = dot_product(diffs,matmul(sn_ninv,diffs))
     BT = SUM(matmul(sn_ninv,diffs))
     
     !! H0 normalisation alla Bridle and co. 
     chisq = AT-BT**2/sn_sumninv

     if (Feedback > 1) write (*,*) 'SN chisq: ', chisq



     SN_LnLike = chisq/2


 end function SN_LnLike


end module snovae
