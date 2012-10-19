!===========================================================================
MODULE spt_reichardt_likelihood ! Parameters are defined in options module

! ===========================================================================

  use highell_options
  use highell_subroutines
  use foregrounds_loading

  logical :: initialise_spt=.true.
  REAL(8), dimension(:,:) :: btt_dat(nspec_r,0:bmax0_r-1)
  REAL(8), dimension(:,:,:) :: btt_var(nspec_r,nspec_r,0:bmax0_r-1)
  REAL(8) ::  bval(nspec_r,0:bmax0_r-1), inverse(1:datap_r,1:datap_r)
  REAL(8), dimension (:), allocatable :: cl_src
  REAL(8) :: win_func(nspec_r,0:bmax0_r-1,1:10000)
  
  PRIVATE
  public :: spt_reichardt_likelihood_init
  public :: spt_reichardt_likelihood_compute

  contains
  
  ! ===========================================================================
  SUBROUTINE spt_reichardt_likelihood_init
  ! ===========================================================================
    
    IMPLICIT NONE
    
    INTEGER  :: i,j,lun,il
    REAL(8)  :: dummy, ii
    CHARACTER(LEN=240) :: ttfilename(nspec_r),winfilename(nspec_r),invcovfilename
    LOGICAL  :: good

    allocate(cl_src(2:tt_lmax))

#ifdef TIMING
       call spt_timing_start( 'spt_reichardt_likelihood_init' )
#endif
    
    !-----------------------------------------------
    ! set file names
    !-----------------------------------------------
    
    ttfilename(1) = trim(SPT_data_dir)//'spt_high/spt_bandpower_95x95.dat'
    ttfilename(2) = trim(SPT_data_dir)//'spt_high/spt_bandpower_95x150.dat'
    ttfilename(3) = trim(SPT_data_dir)//'spt_high/spt_bandpower_95x220.dat'
    ttfilename(4) = trim(SPT_data_dir)//'spt_high/spt_bandpower_150x150.dat'
    ttfilename(5) = trim(SPT_data_dir)//'spt_high/spt_bandpower_150x220.dat'
    ttfilename(6) = trim(SPT_data_dir)//'spt_high/spt_bandpower_220x220.dat'
    invcovfilename   = trim(SPT_data_dir)//'spt_high/inverse_full.txt'
   
    winfilename(1) = trim(SPT_data_dir)//'spt_high/BblMean_95x95_new.dat'
    winfilename(2) = trim(SPT_data_dir)//'spt_high/BblMean_95x150_new.dat'
    winfilename(3) = trim(SPT_data_dir)//'spt_high/BblMean_95x220_new.dat'  
    winfilename(4) = trim(SPT_data_dir)//'spt_high/BblMean_150x150_new.dat'
    winfilename(5) = trim(SPT_data_dir)//'spt_high/BblMean_150x220_new.dat'
    winfilename(6) = trim(SPT_data_dir)//'spt_high/BblMean_220x220_new.dat'

    !-----------------------------------------------
    ! load TT data 
    !----------------------------------------------
    do j=1,nspec_r

       inquire(file=ttfilename(j),exist = good)
       if(.not.good)then
          write(*,*) 'cant find', trim(ttfilename(j)), trim(SPT_data_dir)
          stop
       endif
       call get_free_lun( lun )
       open(unit=lun,file=ttfilename(j),form='formatted',status='unknown',action='read')
       do i=0,bmax0_r-1
          read(lun,*) bval(j,i),btt_dat(j,i), btt_var(j,j,i)
       enddo
       close(lun)
 
       inquire (file=winfilename(j),exist = good)
       if(.not.good)then
          write(*,*) 'cant find', trim(winfilename(j)), trim(SPT_data_dir)
          stop
       endif
       call get_free_lun( lun )
       open(unit=lun,file=winfilename(j),form='formatted',status='unknown',action='read')
       win_func(j,0:bmax0_r-1,1:10000)=0.d0 
       do il = 2, tt_lmax
          read(lun,*) ii, (win_func(j,i,il), i=0,bmax0_r-1)       
       enddo
       close(lun) 
   enddo
   
   !------------------------------------------------- 
   !Read spt inverse covariance matrix 
   !-------------------------------------------------

    call get_free_lun( lun )
    open(unit=lun,file=invcovfilename,form='formatted',status='unknown',action='read')
    do i=1,datap_r
       read(lun,*) inverse(i,1:datap_r)
    enddo
    close(lun)

    initialise_spt = .false.

#ifdef TIMING
      call spt_timing_end()
#endif

  END SUBROUTINE spt_reichardt_likelihood_init
 
 ! ==========================================================================================================================================
  SUBROUTINE spt_reichardt_likelihood_compute(cltt,amp_tsz,amp_ksz,xi,aps95,aps150,aps220,acib150,acib220,rps0,rps1,rps,rcib,cal_1,cal_2,cal_3,like_sptr)
 ! ==========================================================================================================================================

    IMPLICIT NONE
    REAL(8), intent(in) :: cltt(2:*), amp_tsz,amp_ksz,xi,aps95,aps150,aps220,acib150,acib220,rps0,rps1,rps,rcib,cal_1,cal_2,cal_3
    REAL(8), intent(out) :: like_sptr
    INTEGER :: lun,il,i,j,k
    REAL(8) :: cltt_temp(2:tt_lmax)
    REAL(8) :: btt_th(nspec_r,0:bmax0_r-1)
    REAL(8) :: diffs(datap_r,1),tmp(datap_r,1),diffs2(1,datap_r),chi2(1,1)
    REAL(8) :: gc1, gc2, gc3
    REAL(8) :: fp1,fp2,fp3,f0,f1,f2,f3,fcal_j,beta_c
    REAL(8) :: sz_corr, planckfunctionratio_corr,flux2tempratio_corr
    REAL(8) :: f1_sz,f1_synch,f1_dust,f2_sz,f2_synch,f2_dust,f3_sz,f3_synch,f3_dust
    REAL(8) :: planckratiod1,planckratiod2,planckratiod3,fluxtempd1,fluxtempd2,fluxtempd3


    fp1 = 100.0d0
    fp2 = 143.0d0
    fp3 = 217.0d0
 
    f1_sz     =97.6d0
    f1_synch  =95.3d0
    f1_dust   =97.9d0
    f2_sz     =152.9d0
    f2_synch  =150.2d0
    f2_dust   =153.8d0
    f3_sz     =218.1d0
    f3_synch  =214.1d0
    f3_dust   =219.6d0
 

    call sz_func(f1_sz,sz_corr)
    f1 = sz_corr
    call sz_func(f2_sz,sz_corr)
    f2 = sz_corr
    call sz_func(f3_sz,sz_corr)
    f3 = sz_corr
    call sz_func(fp2,sz_corr)
    f0 = sz_corr
    call planckfunctionratio(f1_dust,fp1,planckfunctionratio_corr)
    planckratiod1 = planckfunctionratio_corr
    call planckfunctionratio(f2_dust,fp2,planckfunctionratio_corr)
    planckratiod2 = planckfunctionratio_corr
    call planckfunctionratio(f3_dust,fp3,planckfunctionratio_corr)    
    planckratiod3 = planckfunctionratio_corr
    call flux2tempratio(f1_dust,fp1,flux2tempratio_corr)
    fluxtempd1 = flux2tempratio_corr
    call flux2tempratio(f2_dust,fp2,flux2tempratio_corr)
    fluxtempd2 = flux2tempratio_corr
    call flux2tempratio(f3_dust,fp3,flux2tempratio_corr)
    fluxtempd3 = flux2tempratio_corr

    gc1 = 0.16d0
    gc2 = 0.21d0
    gc3 = 2.19d0
    beta_c = 2.20d0

    !----------------------------------------------------------------
    ! Calculate theory
    !----------------------------------------------------------------
    
    do j=1,nspec_r
       cltt_temp(2:tt_lmax)=0.d0
       do il=2,tt_lmax
          if(j==1) then 
             cl_src(il) = aps95*cl_p(il)
             cltt_temp(il) = cltt(il)+cl_src(il)+f1*f1/f0/f0*amp_tsz*cl_tsz(il)+amp_ksz*cl_ksz(il)
           else if (j==2) then
             cl_src(il) = rps0*sqrt(aps95*aps150)*cl_p(il)
             cltt_temp(il) = cltt(il)+cl_src(il)+f1/f0*f2/f0*amp_tsz*cl_tsz(il)+amp_ksz*cl_ksz(il)-sqrt(acib150*amp_tsz*4.796*f1*f1/f0/f0)*xi*cl_szcib(il)*(f2_dust/fp2)**beta_c*(planckratiod2*fluxtempd2)
           else if (j==3) then
             cl_src(il) = rps1*sqrt(aps95*aps220)*cl_p(il)
             cltt_temp(il) = cltt(il)+cl_src(il)+amp_ksz*cl_ksz(il)-sqrt(acib220*amp_tsz*4.796*f1*f1/f0/f0)*xi*cl_szcib(il)*(f3_dust/fp3)**beta_c*(planckratiod3*fluxtempd3)
           else if (j==4) then
             cl_src(il) = aps150*cl_p(il)+acib150*cl_c(il)*(f2_dust/fp2)**(2.0*beta_c)*(planckratiod2*fluxtempd2)**2.0 &
                          -2.0*sqrt(acib150*amp_tsz*4.796*f2*f2/f0/f0)*xi*cl_szcib(il)*(f2_dust/fp2)**beta_c*(planckratiod2*fluxtempd2)
             cltt_temp(il) =cltt(il)+cl_src(il)+f2*f2/f0/f0*amp_tsz*cl_tsz(il)+amp_ksz*cl_ksz(il)
          else if (j==5) then
             cl_src(il) = rps*sqrt(aps150*aps220)*cl_p(il)+rcib*sqrt(acib150*acib220)*cl_c(il)*(f2_dust/fp2)**beta_c*(planckratiod2*fluxtempd2)*(f3_dust/fp3)**beta_c*(planckratiod3*fluxtempd3)
             cltt_temp(il) =cltt(il)+cl_src(il)+amp_ksz*cl_ksz(il)-sqrt(acib220*amp_tsz*4.796*f2*f2/f0/f0)*xi*cl_szcib(il)*(f3_dust/fp3)**beta_c*(planckratiod3*fluxtempd3)
          else if(j ==6) then
             cl_src(il) = aps220*cl_p(il)+acib220*cl_c(il)*(f3_dust/fp3)**(2.0*beta_c)*(planckratiod3*fluxtempd3)**2.0
             cltt_temp(il) =cltt(il)+cl_src(il)+amp_ksz*cl_ksz(il)
          endif
       enddo
     btt_th(j,0:bmax0_r-1)=MATMUL(win_func(j,0:bmax0_r-1,2:tt_lmax),cltt_temp(2:tt_lmax))
    enddo
   
    !--------------------------------------------------------------
    ! Calibrate theory
    !--------------------------------------------------------------

    do j=1,nspec_r
       if(j ==1 ) fcal_j = cal_1*cal_1
       if(j ==2 ) fcal_j = cal_1*cal_2
       if(j ==3 ) fcal_j = cal_1*cal_3
       if(j ==4 ) fcal_j = cal_2*cal_2
       if(j ==5 ) fcal_j = cal_2*cal_3
       if(j ==6 ) fcal_j = cal_3*cal_3
       btt_th(j,0:bmax0_r-1) = btt_th(j,0:bmax0_r-1)/fcal_j
    enddo

    !--------------------------------------------------------------
    ! chi2 calculation
    !--------------------------------------------------------------

    like_sptr = 0.d0

    do i =0,bmax0_r-1
       do j =0,nspec_r-1
          diffs(i+1+j*bmax0_r,1) = btt_dat(j+1,i) - btt_th(j+1,i)
       enddo
    enddo

    do i =1,datap_r
       diffs2(1,i) = diffs(i,1)
    enddo


    tmp(:,:) = matmul(inverse(:,:),diffs(:,:))
    chi2(:,:) = matmul(diffs2(:,:),tmp(:,:)) 

    like_sptr = like_sptr+chi2(1,1)/2.0

   10  continue
    
#ifdef TIMING
       call spt_timing_end()
#endif

  end SUBROUTINE spt_reichardt_likelihood_compute
  !================================================================================

END MODULE spt_reichardt_likelihood

