!
! UNION2 Supernovae Ia dataset
!
! This module uses the SCP (Supernova Cosmology Project) Union 2 
! compilation. Please cite 
! "Amanullah et al. (SCP) 2010, arXiv:1004.1711 (ApJ accepted)".  
! and the references of other compiled supernovae data are in there.
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
!
! Update Note :
!
! Union1 (Kowalski et al 2008)  : 307 SNe with SALT1 fit (Guy et al 2005) 
! Union2 (Amanullah et al 2010) : 557 SNe with SALT2 fit (Guy et al 2007)
!
! The following parameters are used to calculate distance moduli 
! (see Amanullah et al. 2010 for complete description)
!
!  alpha 0.120887675685  ! Stretch Correction Factor
!  beta  2.51356117225   ! Color   Correction Factor
!  M(h=0.7, statistical only) -19.3111817501  ! Absolute B Magnitue of SNIa
!  M(h=0.7, with systematics) -19.3146267582  
!
!  Tips for running cosmomc (ver Jan, 2010) with SCP UNION2 data 
!
!  1) Place the following 3 data files in your cosmomc data dir (DataDir)
!     sn_z_mu_dmu_union2.txt     : SN data : SN name, z, distance moduli mu, mu error             
!     sn_covmat_sys_union2.txt   : Covariance Matrix with    systematic error
!     sn_covmat_nosys_union2.txt : Covariance Matrix without systematic error
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
!     either rename supernovae_union2.f90 as supernovae.f90 and recompile it 
!     or change targets in your Makefile from supernova to supernovae_union2
!
!  Note: In your default params.ini, there is a line for 'SN_filename', but this
!     union2 module does not use it.  You can leave it as it is, and cosmomc
!     runs without any error but that information is not used.
!     To avoid confusion, you may want to comment it out.
!
!   Update Note by Nao Suzuki (LBNL)

module snovae
use cmbtypes
use MatrixUtils
implicit none

 integer, parameter :: SN_num = 557
 double precision, parameter :: Pi_num = 3.14159265359D0 
 double precision :: SN_z(SN_num), SN_moduli(SN_num)
 double precision :: SN_Ninv(SN_num,SN_Num)
 double precision :: SN_sumninv

 logical, parameter :: SN_marg = .True.

! The following line selects which error estimate to use
! default .True. = with systematic errors
 logical, parameter :: SN_syscovmat = .True.  !! Use covariance matrix with or without systematics

! The following line is not used by this module but it is needed for SDSSII supernova module
 character(len=256) :: SN_filename = ''

contains


 subroutine SN_init
   use settings
   character (LEN=20):: name
   integer i
   real :: tmp_mat(sn_num, sn_num)

   if (Feedback > 0) write (*,*) 'Reading: supernovae data'
   call OpenTxtFile(trim(DataDir)//'sn_z_mu_dmu_union2.txt',tmp_file_unit)
   do i=1,  sn_num
      read(tmp_file_unit, *) name, SN_z(i), SN_moduli(i)
   end do
   close(tmp_file_unit)

   if (SN_syscovmat) then
      call OpenTxtFile(trim(DataDir)//'sn_covmat_sys_union2.txt',tmp_file_unit)
   else
      call OpenTxtFile(trim(DataDir)//'sn_covmat_nosys_union2.txt',tmp_file_unit)
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
