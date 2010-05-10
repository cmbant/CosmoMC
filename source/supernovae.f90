!Edit SN_usexxx parameters below; SNLS by default
!SN_useRiess: Use the gold sample SNe Ia of Riess et al.(2004)
!   See astro-ph/0402512 
!SN_useSNLS: use SNLS astro-ph/0510447

!May 2004 (modified by David Rapetti)
!May 2006: includes SNLS option (thanks to Anze Slosar)

!Marginalize anayltically over H_0 with flat prior 
!(equivalent to mariganlize over M, absolute magnitude; see appendix F of cosmomc paper) 
!resultant log likelihood has arbitary origin, but is returned equal
!to the best-fit value.
module snovae
use cmbtypes
implicit none

 integer, parameter :: riessN=157, snlsN=115
 integer, parameter :: SN_num = riessN+snlsN
 double precision, parameter :: Pi_num = 3.14159265359D0 
 double precision :: SN_z(SN_num), SN_moduli(SN_num), SN_diagerr(SN_num)
 double precision :: SN_Ninv(SN_num,SN_Num), SN_Ninvmarge(SN_num,SN_Num)
 double precision :: aprima, bprima, cprima, a1
 double precision :: SN_trNinv
 logical :: do_SN_init = .true.

 logical, parameter :: SN_marg = .true.
 logical, parameter :: SN_useRiess=.false.
 logical, parameter :: SN_useSNLS=.true.

contains


 subroutine SN_init
   character (LEN=1200) :: InLine
   character (LEN=20) :: names(186)
   double precision :: input (186,3)
   integer gold(riessN)
   integer i

   if (Feedback > 0) write (*,*) 'reading: supernovae data'
   call OpenTxtFile('data/sn_data_riess.dat',tmp_file_unit)
   read(tmp_file_unit,'(a)') InLine
   do i=1, 186
      read(tmp_file_unit, *) names(i), input(i,:)
   end do
   close(tmp_file_unit)

!Use only those supernovae used for `gold sample'
   gold = (/(I, I=1, 29), 32, 33, 34, 35, 37, (I, I=42,45)&
        ,(I,I=48,57),60,62,63,64,65,68,69,(I,I=71,90),(I,I=92,105)&
        ,107,(I,I=109,121),(I,I=123,131),133,134,135,136,137,138&
        ,(I,I=140,144),(I,I=146,161),163,164,169,170,171,173&
        ,174,(I,I=176,186)/)

   SN_z(1:riessN) = input(gold,1)
   SN_moduli(1:riessN) = input(gold,2)
   SN_diagerr(1:riessN) = 1/input(gold,3)**2 !included the extra velocities (see Riess et al.(2004))

!!! Now add snls
   call OpenTxtFile('data/snls.dat',tmp_file_unit)
   do i=riessN+1, sn_Num
      read(tmp_file_unit, *) SN_z(i), sn_moduli(i), sn_diagerr(i)
      sn_diagerr(i)=1d0/sn_diagerr(i)**2
   end do
   close(tmp_file_unit)
   

   SN_trNinv = sum(SN_diagerr)
   do_SN_init = .false.

 end subroutine SN_init

 function SN_LnLike(CMB)
  !Assume this is called just after CAMB with the correct model
  use camb
  implicit none
  type(CMBParams) CMB
  real SN_LnLike
  integer i,kk,aa,bb
  double precision z
  real diffs(SN_num), chisq

  if (do_SN_init) call SN_init

     !$OMP PARALLEL DO DEFAULT(SHARED),SCHEDULE(STATIC), PRIVATE(z,i)
     do i=1, SN_num
        !Obviously this is not v efficient...

        if ((i.le.riessN).and.(.not.SN_useRiess)) cycle
        if ((i.gt.riessn).and.(.not.SN_useSNLS)) cycle
        z= SN_z(i)
        diffs(i) = 5*log10((1+z)**2*AngularDiameterDistance(z))+25 - SN_moduli(i)
!        print *, 5*log10((1+z)**2*AngularDiameterDistance(z)*1e5/CMB%h), SN_moduli(i), cmb%h
     end do
     !$OMP END PARALLEL DO 


!!! WE do analytical marginalisation separatelly for riess and SNLS

     chisq=0

     do kk=1,2

        if ((kk.eq.1).and.(.not.SN_useRiess)) cycle
        if ((kk.eq.2).and.(.not.SN_useSNLS)) cycle


        if (kk.eq.1) then 
           aa=1
           bb=riessN
        else
           aa=riessN+1
           bb=SN_Num
        end if

        !analytical marginalization over H_0 (equivaltent for M, absolute magnitude) 
        
        
        aprima = sum(diffs(aa:bb)**2*SN_diagerr(aa:bb))
        bprima = sum(diffs(aa:bb)*SN_diagerr(aa:bb))
        cprima = sum(SN_diagerr(aa:bb))
!        print *,aprima, cprima
     
        if (sn_marg) then 
           ! to calculate the abslolute chisquare 
           
           a1=log(cprima/(2*Pi_num))
           chisq=chisq+a1+aprima-((bprima**2)/(cprima))
        else
           chisq=chisq+aprima
        end if

     end do


!     if (Feedback > 1) write (*,*) 'SN chisq: ', chisq
     SN_LnLike = chisq/2

!     stop
 end function SN_LnLike


end module snovae
