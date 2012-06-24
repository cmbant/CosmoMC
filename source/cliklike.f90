module cliklike
  use clik
  use cmbtypes
  use settings

  implicit none

  logical :: use_clik = .false.  
  logical, save :: initialised = .false.
  logical :: using_CAMspec = .false.
  logical :: using_PLik = .false.

  integer, parameter :: dp = kind(1.d0)
  integer(kind=4),dimension(6) :: clik_has_cl,clik_lmax
  integer :: clik_n,clik_ncl,clik_nnuis
  integer :: i,j,l
  integer, dimension(4) :: mapped_index
  real(dp), dimension(:), allocatable :: clik_cl_and_pars

  character (len=256) :: clik_filename
  character (len=2),dimension(6) :: clnames
  character (len=256), dimension(:), pointer :: names

  type(clik_object) :: clikid
  

  private
  public :: clik_lnlike, use_clik, clik_filename, using_CAMspec, using_PLik
  
contains

  function clik_lnlike(cl,clik_nuis)
    real(dp) :: clik_lnlike
    real(dp), intent(in) :: cl(lmax,num_cls_tot)
    real(dp), intent(in), optional :: clik_nuis(1:num_freq_params-1)   
    integer :: offset

    if (.not. initialised) then
       Print*,'Initialising clik...'
       call clik_init(clikid,clik_filename)
       call clik_get_has_cl(clikid,clik_has_cl)
       call clik_get_lmax(clikid,clik_lmax)

!Safeguard
       if (lmax .lt. maxval(clik_lmax)+500) then
          print*,'lmax too low: it should at least be set to',(maxval(clik_lmax)+500)
          call MpiStop
       end if



!output Cls used
       clnames(1)='TT'
       clnames(2)='EE'
       clnames(3)='BB'
       clnames(4)='TE'
       clnames(5)='TB'
       clnames(6)='EB'
       print*,'Likelihood uses the following Cls:'
       do i=1,6
          if (clik_has_cl(i) .eq. 1) then
             print*,'  ',trim(clnames(i)),' from l=0 to l=',clik_lmax(i)
          end if
       end do
       
       clik_ncl = sum(clik_lmax) + 6 
       
       clik_nnuis = clik_get_extra_parameter_names(clikid,names)
 

!Another safeguard
       if (clik_nnuis .gt. (num_freq_params-1)) then
          Print*,'You are not supplying all required nuisance parameters. Stopping.'
          stop
       end if
        
!WARNING: The following part is specific for use with CAMspec v4.1 and PLik v2.  
!If future versions of these likelihoods use different nuisance parameters,
!it will need to be updated.

       if (clik_nnuis .eq. 11) then
          Print*,'You appear to be using a CAMspec v4.X likelihood file.'
          using_CAMspec = .true.
       else if (clik_nnuis .eq. 2) then
          Print*,'You appear to be using a PLik likelihood file.'
          using_PLik = .true.
       else if (clik_nnuis .eq. 9) then
          Print*,'You appear to be using a CAMspec v3.X or earlier likelihood file &
               & which is no longer supported by cliklike. You will have to modify &
               & cliklike.f90 to ensure that the nuisance parameters are correctly &
               & passed to the likelihood function in calclike.f90. &
               & The following nuisance parameters are required:'
          do i=1,clik_nnuis
             Print*,trim(names(i))
          end do
          stop
       else if (clik_nnuis .eq. 3) then
          Print*,'You appear to be using a PLik v1 or earlier likelihood file &
               & which is no longer supported by cliklike. You will have to modify &
               & cliklike.f90 to ensure that the nuisance parameters are correctly &
               & passed to the likelihood function in calclike.f90. &
               & The following nuisance parameters are required:'
          do i=1,clik_nnuis
             Print*,trim(names(i))
          end do
          stop
       else if (clik_nnuis .ne. 0) then
          Print*,'Unknown likelihood format.  Make sure that nuisance parameters are &
               & correctly passed to the likelihood function in calclike.f90 before &
               & proceeding. The following nuisance parameters are required:'
          do i=1,clik_nnuis
             Print*,trim(names(i))
          end do
          stop
       end if


       if (clik_nnuis .ne. 0) then           
          Print*,'Clik will run with the following nuisance parameters:'
          do i=1,clik_nnuis
             Print*,trim(names(i))
          end do
       end if



!tidying up
       if (clik_nnuis .gt. 0) deallocate(names)
       
       clik_n = clik_ncl + clik_nnuis
       allocate(clik_cl_and_pars(clik_n))
       
!Mapping CosmoMC's power spectrum indices to clic's
       mapped_index(1) = 1
       mapped_index(2) = 3
       mapped_index(3) = 4
       mapped_index(4) = 2
       
       initialised = .true.

    end if

 
!set C_l and parameter vector to zero initially
    clik_cl_and_pars = 0.d0

    j = 1

!TB and EB assumed to be zero
!If your model predicts otherwise, this function will need to be updated
    do i=1,4
       do l=0,clik_lmax(i)
          !skip C_0 and C_1
          if (l .ge. 2) then
             clik_cl_and_pars(j) = cl(l,mapped_index(i))
          end if
          j = j+1
       end do
    end do
 
!Appending nuisance parameters     
    if (clik_nnuis .ne. 0) then 
       if (using_CAMspec) offset = 0
       if (using_PLik) offset = 11
       do i=1,clik_nnuis
          clik_cl_and_pars(j) = clik_nuis(i+offset)
          j = j+1
        end do
    end if   


!Get - ln like needed by CosmoMC
    clik_lnlike = -1.d0*clik_compute(clikid,clik_cl_and_pars)

    Print*,'clik lnlike = ',clik_lnlike


  end function clik_lnlike




end module cliklike
