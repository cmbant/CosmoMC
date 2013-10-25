! ===========================================================================
MODULE HIGHELL_OPTIONS
! This module contains the options and settings used in the likelihood code
!
! ===========================================================================

!Options
!---------------------------------------------------
! location of input data
! ---------------------------------------------------
  character(len=5000) :: data_dir = '/usersVol2/erminia/planck_groupshare/erminia/highL_clean_finalACT/data/'!Define this to point to your data dir
  character(len=5000) :: ACT_data_dir = '/usersVol2/erminia/planck_groupshare/erminia/highL_clean_finalACT/data/data_act/'
  character(len=5000) :: SPT_data_dir = '/usersVol2/erminia/planck_groupshare/erminia/highL_clean_finalACT/data/data_spt/'

!---------------------------------------------------
! change these to include/exclude experiments
!---------------------------------------------------
  logical :: use_act_south    = .true. ! include ACT south data
  logical :: use_act_equa     = .true. ! include ACT equa data
  logical :: use_spt_lowell   = .false. ! include SPT keisler data
  logical :: use_spt_highell  = .true. ! include SPT reichardt data

!---------------------------------------------------
!These can be changed to choose your ell range for ACT
!---------------------------------------------------
  integer :: lmin11  = 1000  !148x148 cannot below 500
  integer :: lmax11  = 10000 !148x148
  integer :: lmin12  = 1500  !148x220 cannot below 1500
  integer :: lmax12  = 10000 !148x220
  integer :: lmin22  = 1500  !220x220 cannot below 1500
  integer :: lmax22  = 10000 !220x220


! Settings (should not be altered)
!---------------------------------------------------
! general settings
!---------------------------------------------------
  integer :: tt_lmax_mc       = 5000
  integer :: tt_lmax          = 10000
  real(8), parameter :: PI    = 3.14159265358979323846264d0
!---------------------------------------------------

!---------------------------------------------------
! likelihood terms from ACT data
!---------------------------------------------------
  integer, parameter :: nbin11  = 30 !data from 500 to 10000
  integer, parameter :: nbin12  = 20 !data from 1000 to 10000
  integer, parameter :: nbin22  = 20 !data from 1000 to 10000
  integer, parameter :: nspec   = 3  !number of theory spectra
  integer, parameter :: tbin    = 33 !theory from 2 to 10000

!south
!---------------------------------------------------
  integer, parameter :: nsp11_s  = 6
  integer, parameter :: nsp12_s  = 9
  integer, parameter :: nsp22_s  = 6
  integer, parameter :: nspec_s  = 21  !number of spectra
  integer, parameter :: datap_s  = 480 !30x6+20x9+20x6 number of bins x season pairs

! equator
!---------------------------------------------------
  integer, parameter :: nsp11_e  = 3
  integer, parameter :: nsp12_e  = 4
  integer, parameter :: nsp22_e  = 3
  integer, parameter :: nspec_e  = 10  !number of spectra
  integer, parameter :: datap_e  = 230 !30x3+20x4+20x3 number of bins x season pairs
!---------------------------------------------------

!---------------------------------------------------
! likelihood terms from SPT Reichardt data
!---------------------------------------------------
  integer, parameter :: nspec_r  = 6 
  integer, parameter :: bmax0_r  = 15 !max nbins in SPT data	
  integer, parameter :: datap_r  = 90 !15 bins x 6 spectra 
!---------------------------------------------------

!---------------------------------------------------
! likelihood terms from SPT Keisler data
!---------------------------------------------------
  integer :: tt_lmax_k = 3300
  integer, parameter :: bmax0_k   = 47 !max nbins in SPT data
!650-2000 selected in the inverse covmat calculation - 27 effective bins
!---------------------------------------------------
 
END MODULE HIGHELL_OPTIONS
