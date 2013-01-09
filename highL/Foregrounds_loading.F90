! ============================================================================
MODULE foregrounds_loading
! ============================================================================
  use highell_options
  use highell_subroutines

  logical :: initialise_foregrounds=.true.
  REAL(8), dimension (:), allocatable :: cl_tsz,cl_ksz,cl_szcib,cl_cir,cl_cirspt,cl_p
  public :: foregrounds_init

   contains

   ! ============================================================================
   SUBROUTINE foregrounds_init
   ! ============================================================================

    IMPLICIT NONE
    INTEGER  :: lun,il
    REAL(8)  :: dummy

    allocate(cl_tsz(2:tt_lmax),cl_ksz(2:tt_lmax),cl_szcib(2:tt_lmax))
    allocate(cl_p(2:tt_lmax),cl_cir(2:tt_lmax),cl_cirspt(2:tt_lmax))

    
    call get_free_lun(lun)
    open(unit=lun,file=trim(data_dir)//'Fg/tsz_143_eps0.50.dat',form='formatted',status='unknown')
    do il=2,tt_lmax
       read(lun,*) dummy,cl_tsz(il)
       cl_tsz(il) = cl_tsz(il)/4.796 !normalization at ell = 3000
    enddo
    close(lun)

    call get_free_lun(lun)
    open(unit=lun,file=trim(data_dir)//'Fg/cl_ksz_148_trac.dat',form='formatted',status='unknown')
    cl_ksz(2:tt_lmax) = 0.d0
    do il=2,tt_lmax-1
       read(lun,*) dummy,cl_ksz(il)
       cl_ksz(il) = cl_ksz(il)/2.05697 !normalization at ell = 3000
    enddo
    close(lun)
    
    do il=2,tt_lmax
       cl_p(il)=(il/3000.d0)**2.0
    enddo

    call get_free_lun(lun)
    open(unit=lun,file=trim(data_dir)//'Fg/sz_x_cib_template.dat',form='formatted',status='unknown')
    cl_szcib(2:tt_lmax) = 0.d0
    do il=2,tt_lmax-1
       read(lun,*) dummy,cl_szcib(il)
    enddo
    close(lun)

    do il=2,tt_lmax
       cl_cir(il)=(il/3000.d0)**(-0.7)
    enddo

    do il=2,tt_lmax
       cl_cirspt(il)=(il/3000.d0)**(-1.2)
    enddo

    initialise_foregrounds=.false.

    END SUBROUTINE foregrounds_init 
    !================================================================================

END MODULE foregrounds_loading
!====================================================================================
