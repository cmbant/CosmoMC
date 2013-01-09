! ===========================================================================
MODULE spt_keisler_likelihood
! Erminia : likelihood AR1 SPT Keisler
! Parameters are defined in SPT_keisler_options module
! ===========================================================================

  use highell_options
  use highell_subroutines
  use foregrounds_loading

  REAL(8) ::  btt_dat(0:bmax0_k-1),btt_var(0:bmax0_k-1)
  REAL(8) ::  inverse(1:bmax0_k,1:bmax0_k),bval(0:bmax0_k-1)
  REAL(8), dimension (:), allocatable :: cl_src, cl_c
  REAL(8), dimension(:,:), allocatable :: win_func
  
  PRIVATE
  public :: spt_keisler_likelihood_init
  public :: spt_keisler_likelihood_compute
  contains
  
  ! ===========================================================================
  SUBROUTINE spt_keisler_likelihood_init
  ! ===========================================================================
    
    IMPLICIT NONE
    
    INTEGER  :: i,j,lun,il
    REAL(8)  :: dummy,ii
    CHARACTER(LEN=240) :: ttfilename, bblfilename, invcovfilename
    LOGICAL  :: good
    
    allocate(cl_c(2:tt_lmax),cl_src(2:tt_lmax_k))

    !-----------------------------------------------
    ! set file names
    !-----------------------------------------------
    
    ttfilename  = trim(SPT_data_dir)//'spt_lowell/spt_kspectrum_150x150.dat'
    invcovfilename = trim(SPT_data_dir)//'spt_lowell/inverse_short.dat'
    bblfilename = trim(SPT_data_dir)//'spt_lowell/BblMean_150x150.dat'

    !-----------------------------------------------------------
    ! load spectrum, covariance and band power window functions 
    !-----------------------------------------------------------
    
    inquire(file=ttfilename,exist = good)
    if(.not.good)then
       write(*,*) 'cant find', trim(ttfilename), trim(SPT_data_dir)
       stop
    endif
    call get_free_lun( lun )
    open(unit=lun,file=ttfilename,form='formatted',status='unknown',action='read')    
    do i=0,bmax0_k-1
       read(lun,*) bval(i),btt_dat(i), btt_var(i)
    enddo
    close(lun)

    call get_free_lun( lun )
    open(unit=lun,file=invcovfilename,form='formatted',status='unknown',action='read')
    do i=1,bmax0_k
       read(lun,*) inverse(i,1:bmax0_k)
    enddo
    close(lun)

    inquire (file=bblfilename,exist = good)
    if(.not.good)then
       write(*,*) 'cant find', trim(bblfilename), trim(SPT_data_dir)
       stop
    endif
    call get_free_lun( lun )
    open(unit=lun,file=bblfilename,form='formatted',status='unknown',action='read')
    allocate(win_func(0:bmax0_k-1,1:tt_lmax_k))
    do il = 2, tt_lmax_k
        read(lun,*) ii, (win_func(i,il), i=0,bmax0_k-1) 
        enddo
    close(lun)

  END SUBROUTINE spt_keisler_likelihood_init
  
  ! ===========================================================================================================================
  SUBROUTINE spt_keisler_likelihood_compute(cltt,amp_tsz,amp_ksz,xi,aps150,acib150,ncib,cal_2,like_sptk)
  ! ===========================================================================================================================
    
    IMPLICIT NONE
    REAL(8), intent(in) :: cltt(2:*), amp_tsz,amp_ksz,xi,aps150,acib150,ncib,cal_2
    REAL(8), intent(out) :: like_sptk 
    INTEGER :: lun,il,i   
    REAL(8) :: cltt_temp(2:tt_lmax_k)
    REAL(8) :: btt_th(0:bmax0_k-1)
    REAL(8) :: diffs(bmax0_k,1),tmp(bmax0_k,1),diffs2(1,bmax0_k),chi2(1,1)
    REAL(8) :: fp2,f2_sz,f2_synch,f2_dust,f2,f0,beta_c
    REAL(8) :: planckratiod2,fluxtempd2
    REAL(8) :: sz_corr, planckfunctionratio_corr, flux2tempratio_corr

    fp2  = 143.d0

    f2_sz     =152.9d0
    f2_synch  =150.2d0
    f2_dust   =153.8d0

    beta_c = 2.20d0

    call sz_func(f2_sz,sz_corr)
    f2 = sz_corr
    call sz_func(fp2,sz_corr)
    f0 = sz_corr
    call planckfunctionratio(f2_dust,fp2,planckfunctionratio_corr)
    planckratiod2 = planckfunctionratio_corr
    call flux2tempratio(f2_dust,fp2,flux2tempratio_corr)
    fluxtempd2 = flux2tempratio_corr

    !----------------------------------------------------------------
    ! Define CIB term
    !----------------------------------------------------------------
    cl_c(2:tt_lmax) = 0.d0
    do il=2,tt_lmax
       cl_c(il)=(il/3000.d0)**ncib
    enddo


    !----------------------------------------------------------------
    ! Calculate theory
    !----------------------------------------------------------------

    do il=2,tt_lmax_k
       cl_src(il) = (aps150+9.2)*cl_p(il)+acib150*cl_c(il)*(f2_dust/fp2)**(2.0*beta_c)*(planckratiod2*fluxtempd2)**2.0 &
                    -2.0*sqrt(acib150*amp_tsz*f2*f2/f0/f0)*xi*cl_szcib(il)*(f2_dust/fp2)**beta_c*(planckratiod2*fluxtempd2)
       cltt_temp(il) =cltt(il)+cl_src(il)+f2*f2/(f0*f0)*amp_tsz*cl_tsz(il)+amp_ksz*cl_ksz(il)

       ! Calibrate theory
       cltt_temp(il) =cltt_temp(il)/(cal_2**2.0)
    enddo
    btt_th(0:bmax0_k-1)=MATMUL(win_func(0:bmax0_k-1,2:tt_lmax_k),cltt_temp(2:tt_lmax_k))
    

    !--------------------------------------------------------------
    ! chi2 calculation
    !--------------------------------------------------------------

    like_sptk = 0.d0

    do i = 1, bmax0_k
       diffs(i,1) = btt_dat(i-1) - btt_th(i-1)
       diffs2(1,i) = diffs(i,1)
    enddo

    tmp(:,:) = matmul(inverse(:,:),diffs(:,:))
    chi2(:,:) = matmul(diffs2(:,:),tmp(:,:))

    like_sptk = like_sptk+chi2(1,1)/2.0

   10  continue
    
  end SUBROUTINE spt_keisler_likelihood_compute  
  !===============================================================================

END MODULE spt_keisler_likelihood

