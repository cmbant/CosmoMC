!
! SDSS MLSC / SALT-II Supernovae Ia datasets, 
! Including HST, SNLS, ESSENCE and Lowz
!
!
! This module uses the data as presented in arXiv:0908.4274v1 
!
! By Wessel Valkenburg, LAPTH, 2009, containing traces of the original 
! code by A Lewis, S Bridle and D Rapetti. 
! E-mail: wessel.valkenburg@lapp.in2p3.fr

include 'sn_sdss_parser.f90'

module snovae
  use cmbtypes
  use MatrixUtils
  implicit none
  
  private

   integer :: SN_num
   double precision, parameter :: Pi_num = 3.14159265359D0 
   double precision, allocatable :: SN_z(:), SN_moduli(:)
   double precision, allocatable :: SN_Ninv(:,:)
   double precision, allocatable :: allDA(:)
   double precision :: SN_sumninv

   character(len=256) :: SN_filename = ''
   
   logical, parameter :: SN_marg = .True.

   public SN_LnLike, SN_filename

 contains

   subroutine SN_init
     use sn_sdss_parser
     implicit none

     type(sdssdata) :: thisdata

     character (LEN=20):: name
     integer i
     real :: tmp_mat(sn_num, sn_num)
     double precision :: thissigz

     if (Feedback > 0) write (*,*) 'Reading: supernovae data (SDSS).'

     ! Get the data - hardcode the right path here:
     if (SN_filename =='') SN_filename = trim(DataDir) // 'supernovae.dataset'
     call GetSDSSSN(SN_filename, thisdata, feedback)


     ! Set the sigmas etc
     SN_Num = size(thisdata%z)

     allocate(SN_z(SN_Num))
     allocate(allDA(SN_Num))
     allocate(SN_moduli(SN_Num))
     allocate(SN_Ninv(SN_Num,SN_Num))

     SN_z = thisdata%z
     SN_moduli = thisdata%mu

     SN_Ninv = 0.d0
     do i=1,SN_Num
        thissigz = thisdata%sigz(i) * (5.d0/dLog(10.d0)) * (1+thisdata%z(i))/(1+thisdata%z(i)/2.d0)/thisdata%z(i)
        SN_Ninv(i,i) = (thisdata%muerr(i)**2+thisdata%sigint(i)**2+thissigz**2)**(-1)
     end do


     SN_sumninv = SUM(sn_ninv)


   end subroutine SN_init

   function SN_LnLike(CMB)
     use camb
     !Assume this is called just after CAMB with the correct model  use camb
     implicit none
     type(CMBParams) CMB
     logical, save :: do_SN_init = .true.

     real SN_LnLike
     integer i, sni
     double precision z, AT, BT
     real, allocatable :: diffs(:)
     real chisq

     if (do_SN_init) then 
        call SN_init
        do_SN_init = .false.
     end if

     allocate(diffs(SN_Num))

     call DzArray(SN_z,allDA)


     diffs(:) = 5*log10((1+SN_z(:))**2*allDA(:))+25 - SN_moduli(:)


     AT = dot_product(diffs,matmul(sn_ninv,diffs))
     if(SN_marg)then
        BT = SUM(matmul(sn_ninv,diffs))
     else
        BT=0.d0
     end if
  
     deallocate(diffs)
     
     !! H0 normalisation alla Bridle and co. 
     chisq = AT-BT**2/sn_sumninv

     if (Feedback > 1) write (*,*) 'SN chisq: ', chisq

     SN_LnLike = chisq/2.

   end function SN_LnLike

  subroutine DzArray(SNz,ThisDz)
    use camb
    implicit none
    !real(dl) :: z, dum
    real(dl), intent(in) :: SNz(:)
    real(dl), intent(out) :: ThisDz(:)
    integer :: i

    
    do i = 1,size(SNz)
       ThisDz(i) = AngularDiameterDistance(SNz(i))
    end do
    
    ! Or write code here for other distances, in a void etc. etc.
    
  end subroutine DzArray
  
end module snovae
