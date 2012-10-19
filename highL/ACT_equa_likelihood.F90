!============================================================================
MODULE act_equa_likelihood 
! Parameters are defined in Highell_options module
! ===========================================================================

  use highell_options
  use highell_subroutines
  use foregrounds_loading

  logical :: initialise_act=.true.
  REAL(8), dimension(:,:) :: btt_dat1(nsp11_e,0:nbin11-1),btt_dat12(nsp12_e,0:nbin12-1),btt_dat2(nsp22_e,0:nbin22-1)
  REAL(8), dimension(:,:) :: bval1(nsp11_e,0:nbin11-1),bval12(nsp12_e,0:nbin12-1),bval2(nsp22_e,0:nbin22-1)
  REAL(8) ::  inverse(1:datap_e,1:datap_e)
  REAL(8), dimension (:), allocatable :: cl_src
  REAL(8) :: win_func(nspec,0:tbin-1,1:10000)
  
  PRIVATE
  public :: act_equa_likelihood_init
  public :: act_equa_likelihood_compute
  public :: get_inverse

contains
  
  ! ============================================================================
  SUBROUTINE act_equa_likelihood_init
  ! ============================================================================
    
    IMPLICIT NONE
    
    INTEGER  :: i,j,lun,il
    REAL(8)  :: inverse_e_new(1:datap_e,1:datap_e)
    REAL(8)  :: dummy, ii
    CHARACTER(LEN=240) :: ttfilename(nspec_e), winfilename(nspec)
    LOGICAL  :: good

    allocate(cl_src(2:tt_lmax))

#ifdef TIMING
    call act_equa_timing_start('act_likelihood_init')
#endif
    
    !------------------------------------------------
    ! set file names
    !------------------------------------------------
    
    ttfilename(1) = trim(ACT_data_dir)//'equa/spectrum_148x148_season3exseason3e.dat'
    ttfilename(2) = trim(ACT_data_dir)//'equa/spectrum_148x148_season3exseason4e.dat'
    ttfilename(3) = trim(ACT_data_dir)//'equa/spectrum_148x148_season4exseason4e.dat'
    ttfilename(4) = trim(ACT_data_dir)//'equa/spectrum_148x220_season3exseason3e.dat'
    ttfilename(5) = trim(ACT_data_dir)//'equa/spectrum_148x220_season3exseason4e.dat'
    ttfilename(6) = trim(ACT_data_dir)//'equa/spectrum_148x220_season4exseason3e.dat'
    ttfilename(7) = trim(ACT_data_dir)//'equa/spectrum_148x220_season4exseason4e.dat'
    ttfilename(8) = trim(ACT_data_dir)//'equa/spectrum_220x220_season3exseason3e.dat'
    ttfilename(9) = trim(ACT_data_dir)//'equa/spectrum_220x220_season3exseason4e.dat'
    ttfilename(10) = trim(ACT_data_dir)//'equa/spectrum_220x220_season4exseason4e.dat'

    winfilename(1) = trim(ACT_data_dir)//'equa/BblMean_148x148_season4exseason4e.dat'
    winfilename(2) = trim(ACT_data_dir)//'equa/BblMean_148x220_season4exseason4e.dat'
    winfilename(3) = trim(ACT_data_dir)//'equa/BblMean_220x220_season4exseason4e.dat'

    !----------------------------------------------
    ! load TT data 
    !----------------------------------------------
    do j=1,nspec_e

       inquire(file=ttfilename(j),exist = good)
       if(.not.good)then
          write(*,*) 'cant find', trim(ttfilename(j)), trim(ACT_data_dir)
          stop
       endif
       call get_free_lun( lun )
       open(unit=lun,file=ttfilename(j),form='formatted',status='unknown',action='read')
       if (j .le. nsp11_e) then 
          do i=0,nbin11-1
             read(lun,*) bval1(j,i),btt_dat1(j,i)
          enddo
       else if (j .ge. (nsp11_e+1) .and. j .le. (nsp11_e+nsp12_e)) then  
          do i=0,nbin12-1
             read(lun,*) bval12(j-nsp11_e,i),btt_dat12(j-nsp11_e,i)
          enddo
       else
          do i=0,nbin22-1
             read(lun,*) bval2(j-(nsp11_e+nsp12_e),i),btt_dat2(j-(nsp11_e+nsp12_e),i)
          enddo
       end if
       close(lun)
    enddo

    !----------------------------------------------
    ! read windows for theory 
    !----------------------------------------------
    do j=1,nspec
       inquire (file=winfilename(j),exist = good)
       if(.not.good)then
          write(*,*) 'cant find', trim(winfilename(j)), trim(ACT_data_dir)
          stop
       endif
       call get_free_lun( lun )
       open(unit=lun,file=winfilename(j),form='formatted',status='unknown',action='read')
       win_func(j,0:tbin-1,1:10000)=0.d0
       do il = 2, tt_lmax
          read(lun,*) ii, (win_func(j,i,il), i=0,tbin-1)       
       enddo
       close(lun) 
    enddo
 

   !------------------------------------------------- 
   !Read inverse covariance matrix 
   !-------------------------------------------------

    call get_inverse(inverse_e_new)
    do i=1,datap_e
       do j = 1,datap_e
          inverse(i,j) = inverse_e_new(i,j)
       enddo
    enddo 

    initialise_act = .false.

#ifdef TIMING
    call act_equa_timing_end()
#endif

  END SUBROUTINE act_equa_likelihood_init

 ! ===================================================================================================================================
  SUBROUTINE act_equa_likelihood_compute(cltt,amp_tsz,amp_ksz,xi,aps148,aps217,acib150,acib220,rps,rcib,cae1,cae2,like_acte)
 ! ===================================================================================================================================

    IMPLICIT NONE
    REAL(8), intent(in) :: cltt(2:*), amp_tsz,amp_ksz,xi,aps148,aps217,acib150,acib220,rps,rcib,cae1,cae2
    REAL(8), intent(out) :: like_acte
    INTEGER :: lun,il,i,j,k
    REAL(8) :: cltt_temp(2:tt_lmax)
    REAL(8) :: btt_th(nspec,0:tbin-1)
    REAL(8) :: diffs(datap_e,1),tmp(datap_e,1),diffs2(1,datap_e),chi2(1,1)
    REAL(8) :: f0,f1,f2,fcal_j,beta_c
    REAL(8) :: f1_sz,f1_synch,f1_dust,f2_sz,f2_synch,f2_dust,fp2,fp3
    REAL(8) :: sz_corr, planckfunctionratio_corr,flux2tempratio_corr
    REAL(8) :: planckratiod1,planckratiod2,fluxtempd1,fluxtempd2
 
    fp2  = 143.d0
    fp3  = 217.d0
    
    f1_sz     =146.9d0
    f1_synch  =147.6d0
    f1_dust   =149.7d0
    f2_sz     =220.2d0
    f2_synch  =217.6d0
    f2_dust   =219.6d0

    call sz_func(f1_sz,sz_corr)
    f1 = sz_corr
    call sz_func(f2_sz,sz_corr)
    f2 = sz_corr
    call sz_func(fp2,sz_corr)
    f0 = sz_corr
    call planckfunctionratio(f1_dust,fp2,planckfunctionratio_corr)
    planckratiod1 = planckfunctionratio_corr
    call planckfunctionratio(f2_dust,fp3,planckfunctionratio_corr)
    planckratiod2 = planckfunctionratio_corr
    call flux2tempratio(f1_dust,fp2,flux2tempratio_corr)
    fluxtempd1 = flux2tempratio_corr
    call flux2tempratio(f2_dust,fp3,flux2tempratio_corr)
    fluxtempd2 = flux2tempratio_corr

    beta_c = 2.2d0

    !----------------------------------------------------------------
    ! Calculate theory
    !----------------------------------------------------------------
    
    do j=1,nspec
       cltt_temp(2:tt_lmax)=0.d0
       do il=2,tt_lmax
          if(j==1) then
             cl_src(il) = aps148*cl_p(il)+acib150*cl_c(il)*(f1_dust/fp2)**(2.0*beta_c)*(planckratiod1*fluxtempd1)**2.0 &
                          -2.0*sqrt(acib150*amp_tsz*4.796*f1*f1/f0/f0)*xi*cl_szcib(il)*(f1_dust/fp2)**beta_c*(planckratiod1*fluxtempd1)
             cltt_temp(il) =cltt(il)+cl_src(il)+f1*f1/f0/f0*amp_tsz*cl_tsz(il)+amp_ksz*cl_ksz(il)
          else if (j==2) then
             cl_src(il) = rps*sqrt(aps148*aps217)*cl_p(il)+rcib*sqrt(acib150*acib220)*cl_c(il)*(f1_dust/fp2)**beta_c*(planckratiod1*fluxtempd1)*(f2_dust/fp3)**beta_c*(planckratiod2*fluxtempd2)
             cltt_temp(il) =cltt(il)+cl_src(il)+amp_ksz*cl_ksz(il)-sqrt(acib220*amp_tsz*4.796*f1*f1/f0/f0)*xi*cl_szcib(il)*(f2_dust/fp3)**beta_c*(planckratiod2*fluxtempd2)
          else if(j ==3) then
             cl_src(il) = aps217*cl_p(il)+acib220*cl_c(il)*(f2_dust/fp3)**(2.0*beta_c)*(planckratiod2*fluxtempd2)**2.0
             cltt_temp(il) =cltt(il)+cl_src(il)+amp_ksz*cl_ksz(il)
          endif
          cltt_temp(il) = cltt_temp(il)/((dble(il)*(dble(il)+1.0))/(2*PI))
      enddo  
      btt_th(j,0:tbin-1)=MATMUL(win_func(j,0:tbin-1,2:tt_lmax),cltt_temp(2:tt_lmax))
    enddo 
 
    !--------------------------------------------------------------
    ! Calibrate theory
    !--------------------------------------------------------------

    do j=1,nspec
       if(j ==1 ) fcal_j = cae1*cae1
       if(j ==2 ) fcal_j = cae1*cae2
       if(j ==3 ) fcal_j = cae2*cae2
       btt_th(j,0:tbin-1) = btt_th(j,0:tbin-1)/fcal_j
    enddo

    !--------------------------------------------------------------
    ! chi2 calculation
    !--------------------------------------------------------------

    like_acte = 0.d0

    diffs(datap_e,1) = 0.d0
    diffs2(1,datap_e) = 0.d0

    do i =0,nbin11-1
       do j=0,nsp11_e-1
          diffs(i+1+j*nbin11,1) = btt_dat1(j+1,i) - btt_th(1,i+4)
       enddo
    enddo
    do i =0,nbin12-1
       do j=0,nsp12_e-1
          diffs(i+1+nsp11_e*nbin11+j*nbin12,1) = btt_dat12(j+1,i) - btt_th(2,i+14)
       enddo
    enddo
    do i =0,nbin22-1
       do j =0,nsp22_e-1
          diffs(i+1+nsp11_e*nbin11+nsp12_e*nbin12+j*nbin22,1) = btt_dat2(j+1,i) - btt_th(3,i+14)
       enddo
    enddo

    do i =1,datap_e
       diffs2(1,i) = diffs(i,1)
    enddo

    tmp(:,:) = matmul(inverse(:,:),diffs(:,:))
    chi2(:,:) = matmul(diffs2(:,:),tmp(:,:))

    like_acte = like_acte+chi2(1,1)/2.0

10  continue
    
#ifdef TIMING
    call act_equa_timing_end()
#endif

  end SUBROUTINE act_equa_likelihood_compute
    
  !================================================================================


  ! ============================================================================
  SUBROUTINE get_inverse(inverse_e_new)
  ! ============================================================================

    IMPLICIT NONE

    INTEGER  :: i,j,k,l,lun,n,stat
    REAL(8)  :: dummy
    CHARACTER(LEN=240) :: binfilename,invcovfilenamee
    LOGICAL  :: good
    REAL(8)  :: covmat_e(1:datap_e, 1:datap_e),inverse_e(1:datap_e,1:datap_e),dum_mat_e(datap_e,datap_e)
    REAL(8), INTENT(OUT) :: inverse_e_new(1:datap_e,1:datap_e)
    INTEGER  :: lmin(nspec),lmax(nspec)
    REAL(8)  :: blmean(1:nspec,0:tbin-1),blmin(1:nspec,0:tbin-1),blmax(1:nspec,0:tbin-1),bmin(nspec),bmax(nspec)

    lmin(1)=lmin11
    lmax(1)=lmax11

    lmin(2)=lmin12
    lmax(2)=lmax12

    lmin(3)=lmin22
    lmax(3)=lmax22


   !------------------------------------------------- 
   !Read inverse covariance matrixes 
   !-------------------------------------------------

    invcovfilenamee= trim(ACT_data_dir)//'equa/Inverse_Realistic_Cov_Mat_Equa.dat'
    binfilename= trim(ACT_data_dir)//'equa/binningFile.dat'

    call get_free_lun( lun )
    open(unit=lun,file=invcovfilenamee,form='formatted',status='unknown',action='read')
    do i=1,datap_e
       read(lun,*) inverse_e(i,1:datap_e)
    enddo
    close(lun)

   !------------------------------------------------- 
   !Invert to get covariance matrixes 
   !-------------------------------------------------

   dum_mat_e(1:datap_e,1:datap_e) = inverse_e(1:datap_e, 1:datap_e)
   n = size(inverse_e,1)
   call dpotrf( 'L', n, dum_mat_e, n, stat )
   call dpotri( 'L', n, dum_mat_e, n, stat )
   do i=1,datap_e
      dum_mat_e(i,1:datap_e)=dum_mat_e(1:datap_e,i)
   enddo
   covmat_e(1:datap_e,1:datap_e)=dum_mat_e(1:datap_e,1:datap_e)


   !------------------------------------------------- 
   !Get ell ranges
   !-------------------------------------------------

   do j=1,nspec
      call get_free_lun( lun )
      open(unit=lun,file=binfilename,form='formatted',status='unknown',action='read')
      do i = 0,tbin-1
         read(lun,*) blmin(j,i), blmax(j,i), blmean(j,i)
      end do
      close(lun)

      do i=0,tbin-1
         if(blmean(j,i) <= lmin(j)) bmin(j)=i+1
         if(blmean(j,i) <= lmax(j)) bmax(j)=i
      enddo
   enddo

   bmin(1) = bmin(1)-4
   bmin(2) = bmin(2)-14
   bmin(3) = bmin(3)-14

   bmax(1) = bmax(1)-4
   bmax(2) = bmax(2)-14
   bmax(3) = bmax(3)-14

   !------------------------------------------------- 
   !Change covmat to select the ell range
   !-------------------------------------------------

   do i=1,nsp11_e
      do j=1,nbin11
         k= nbin11*(i-1)+j
         if( k .ge. nbin11*(i-1) .and. k .lt. nbin11*(i-1)+bmin(1)+1) then
             covmat_e(k,1:datap_e) = 0
             covmat_e(1:datap_e,k) = 0
             covmat_e(k,k) = 1E+10
         end if
         if( k .gt. nbin11*(i-1)+bmax(1)+1 .and. k .le. nbin11*(i+1)) then
             covmat_e(k,1:datap_e) = 0
             covmat_e(1:datap_e,k) = 0
             covmat_e(k,k) = 1E+10
         end if
      enddo
   enddo
   do i=1,nsp12_e
      do j=1,nbin12
         k= nbin11*nsp11_e+nbin12*(i-1)+j
         if( k .ge. nbin11*nsp11_e+nbin12*(i-1) .and. k .lt. nbin11*nsp11_e+nbin12*(i-1)+bmin(2)+1) then
             covmat_e(k,1:datap_e) = 0
             covmat_e(1:datap_e,k) = 0
             covmat_e(k,k) = 1E+10
         end if
         if( k .gt. nbin11*nsp11_e+nbin12*(i-1)+bmax(2)+1 .and. k .le. nbin11*nsp11_e+nbin12*(i+1)) then
             covmat_e(k,1:datap_e) = 0
             covmat_e(1:datap_e,k) = 0
             covmat_e(k,k) = 1E+10
         end if
      enddo
   enddo

   do i=1,nsp22_e
      do j=1,nbin22
         k= nbin11*nsp11_e+nbin12*nsp12_e+nbin22*(i-1)+j
         if( k .ge. nbin11*nsp11_e+nbin12*nsp12_e+nbin22*(i-1) .and. k .lt. nbin11*nsp11_e+nbin12*nsp12_e+nbin22*(i-1)+bmin(3)+1) then
             covmat_e(k,1:datap_e) = 0
             covmat_e(1:datap_e,k) = 0
             covmat_e(k,k) = 1E+10
         end if
         if( k .gt. nbin11*nsp11_e+nbin12*nsp12_e+nbin22*(i-1)+bmax(3)+1 .and. k .le. nbin11*nsp11_e+nbin12*nsp12_e+nbin22*(i+1)) then
             covmat_e(k,1:datap_e) = 0
             covmat_e(1:datap_e,k) = 0
             covmat_e(k,k) = 1E+10
         end if
      enddo
   enddo


   !------------------------------------------------- 
   !Get new inverse covmats 
   !-------------------------------------------------

   dum_mat_e(1:datap_e,1:datap_e) = covmat_e(1:datap_e,1:datap_e)
   n = size(covmat_e,1)
   call dpotrf( 'L', n, dum_mat_e, n, stat )
   call dpotri( 'L', n, dum_mat_e, n, stat )
   do i=1,datap_e
      dum_mat_e(i,1:datap_e)=dum_mat_e(1:datap_e,i)
   enddo
   inverse_e_new(1:datap_e,1:datap_e)=dum_mat_e(1:datap_e,1:datap_e)

   END SUBROUTINE get_inverse
  !================================================================================


END MODULE act_equa_likelihood
