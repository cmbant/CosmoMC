module cliklike
  use clik
  use cmbtypes
  use settings
  use CMB_Cls, only: CMB_lensing  

  implicit none

  integer, parameter :: dp = kind(1.d0)

  type clik_input_cls
    real(dp), dimension(:), allocatable :: clik_cl_and_pars
  end type clik_input_cls


  integer :: numcliksets
  integer,parameter :: maxcliksets = 4

  logical :: use_clik = .false.  
  logical, save :: initialised = .false.
  logical, dimension(maxcliksets) :: using_CAMspec = .false.
  logical, dimension(maxcliksets) :: using_ACTSPT = .false.
  logical, dimension(maxcliksets) :: is_lensing = .false.


  integer(kind=4), dimension(6,maxcliksets) :: clik_has_cl,clik_lmax
  integer(kind=4), dimension(maxcliksets) :: clik_lensing_lmax
  integer, dimension(maxcliksets) :: clik_n,clik_ncl,clik_nnuis
  integer :: i,j,l
  integer, dimension(4) :: mapped_index



!Number of nuisance parameters expected for different data sets
  integer :: num_CAMspec_ACTSPT_common = 7

  type(clik_input_cls), dimension(maxcliksets) :: clik_input

  character(len=256), dimension(maxcliksets) :: clik_filename
  character(len=2),dimension(6) :: clnames
  character(len=256), dimension(:), pointer :: names

  type(clik_object), dimension(maxcliksets) :: clikid


  private
  public :: clik_readParams, clik_lnlike, use_clik, clik_filename, numcliksets
  
contains

  subroutine clik_readParams(Ini,num)
    Type(TIniFile) Ini
    integer :: num,n

       if (num .gt. maxcliksets) then
          Print*,'too many clik data sets'
          call MPIstop
       end if
       if (feedback .gt. 1) print*,'Using ',num,' clik likelihood files.'
       if (num .eq. 0) then
          Print*,'No clik data sets selected. Stopping.'
          call MPIstop
       end if

       do n=1,num
          clik_filename(n) = Ini_Read_String_File(Ini,numcat('clik_likefile',n), .false.)
          if (feedback .gt. 1) print*,'Using clik with likelihood file ',trim(clik_filename(n))
          call clik_try_lensing(is_lensing(n),clik_filename(n))
       end do

  end subroutine clik_readParams


  subroutine clik_initialise
     implicit none
     integer :: n

    Print*,'Initialising clik...'
    do n=1,numcliksets

       if (is_lensing(n)) then
          call clik_lensing_initialise(n)
          if (.not. CMB_lensing) then
             Print*,'CMB_lensing must be true if you want to use lensing data.'
             call MPIstop
          end if
       else

          call clik_init(clikid(n),clik_filename(n))
          call clik_get_has_cl(clikid(n),clik_has_cl(:,n))
          call clik_get_lmax(clikid(n),clik_lmax(:,n))

!Safeguard
          if ((lmax .lt. maxval(clik_lmax(:,n))+500) .and. (lmax .lt. 4500)) then
             print*,'lmax too low: it should at least be set to',min(4500,(maxval(clik_lmax(:,n))+500))
             call MPIstop
          end if

!output Cls used
          clnames(1)='TT'
          clnames(2)='EE'
          clnames(3)='BB'
          clnames(4)='TE'
          clnames(5)='TB'
          clnames(6)='EB'
          print*,'Likelihood uses the following Cls:'
          do j=1,6
             if (clik_has_cl(j,n) .eq. 1) then
                print*,'  ',trim(clnames(j)),' from l=0 to l=',clik_lmax(j,n)
             end if
          end do

          clik_ncl(n) = sum(clik_lmax(:,n)) + 6

          clik_nnuis(n) = clik_get_extra_parameter_names(clikid(n),names)


!Another safeguard
          if (clik_nnuis(n) .gt. (num_freq_params-1)) then
             Print*,'You are not supplying all required nuisance parameters. Stopping.'
             call MPIstop
          end if


          if (clik_nnuis(n) .eq. num_CAMspec) then
             Print*,'You appear to be using a CAMspec v6+ likelihood file.'
             using_CAMspec(n) = .true.
          else if (clik_nnuis(n) .eq. num_ACTSPT) then
             Print*,'You appear to be using an ACT/SPT likelihood file.'
             using_ACTSPT(n) = .true.
          else if (clik_nnuis(n) .ne. 0) then
             Print*,'Unknown likelihood format.  Make sure that nuisance parameters are &
                  & correctly passed to the likelihood function in calclike.f90 before &
                  & proceeding. The following nuisance parameters are required:'
             do j=1,clik_nnuis(n)
                Print*,trim(names(j))
             end do
             call MPIstop
          end if

          if (clik_nnuis(n) .ne. 0) then
             Print*,'Clik will run with the following nuisance parameters:'
             do j=1,clik_nnuis(n)
                Print*,j,trim(names(j))
             end do
          end if
         
!tidying up
          if (clik_nnuis(n) .gt. 0) deallocate(names)

          clik_n(n) = clik_ncl(n) + clik_nnuis(n)
          allocate(clik_input(n)%clik_cl_and_pars(clik_n(n)))

!Mapping CosmoMC's power spectrum indices to clik's
          mapped_index(1) = 1
          mapped_index(2) = 3
          mapped_index(3) = 4
          mapped_index(4) = 2

       end if

    end do

  end subroutine clik_initialise


  subroutine clik_lensing_initialise(n)
     implicit none
     integer, intent(in) :: n

     call clik_lensing_init(clikid(n),clik_filename(n))
     call clik_lensing_get_lmax(clikid(n),clik_lensing_lmax(n))

     clik_nnuis(n) = clik_lensing_get_extra_parameter_names(clikid(n),names)

     if (clik_nnuis(n) .ne. 0) then
        Print*,'Clik will run with the following nuisance parameters:'
        do j=1,clik_nnuis(n)
           Print*,j,trim(names(j))
        end do
     end if

!tidying up
     if (clik_nnuis(n) .gt. 0) deallocate(names)

     clik_n(n) = 2*(clik_lensing_lmax(n)+1) + clik_nnuis(n)
     allocate(clik_input(n)%clik_cl_and_pars(clik_n(n)))

  end subroutine clik_lensing_initialise


  function clik_lnlike(cl,clik_nuis)
    real(dp) :: clik_lnlike
    real(dp) :: cl(lmax,num_cls_tot)
    real(dp), intent(in), optional :: clik_nuis(1:num_freq_params-1)   
    integer :: offset,n,m

    clik_lnlike = 0.d0

    if (.not. initialised) then
       call clik_initialise
       initialised = .true.
    end if


!loop over all data sets
    do n=1,numcliksets 

       if (is_lensing(n)) then
          clik_lnlike = clik_lnlike + clik_lensing_lnlike(cl,n,clik_nuis)
       else


!set C_l and parameter vector to zero initially
          clik_input(n)%clik_cl_and_pars = 0.d0

          j = 1

!TB and EB assumed to be zero
!If your model predicts otherwise, this function will need to be updated
          do i=1,4
             do l=0,clik_lmax(i,n)
                !skip C_0 and C_1
                if (l .ge. 2) then
                   clik_input(n)%clik_cl_and_pars(j) = cl(l,mapped_index(i))
                end if
                j = j+1
             end do
          end do

 
!Appending nuisance parameters
!Not pretty. Oh well.     
          if (clik_nnuis(n) .ne. 0) then

             if (using_CAMspec(n)) then
                do i=1,clik_nnuis(n)
                   clik_input(n)%clik_cl_and_pars(j) = clik_nuis(i)
                   j = j+1
                end do
             end if


             if (using_actspt(n)) then
               offset = num_CAMspec
               m = 1
               do i=1,clik_nnuis(n)
       	          ! Fill in common nuisance parameters
                  if (i .eq. 1) then
                     clik_input(n)%clik_cl_and_pars(j) = clik_nuis(6)   ! A_sz
                  else if (i .eq. 2) then 
                     clik_input(n)%clik_cl_and_pars(j) = clik_nuis(14)  ! A_ksz
                  else if (i .eq. 3) then
                     clik_input(n)%clik_cl_and_pars(j) = clik_nuis(13)  ! xi_sz_cib
                  else if (i .eq. 9) then
                     clik_input(n)%clik_cl_and_pars(j) = clik_nuis(4)   ! A_cib_143
       	          else if (i .eq. 10) then
                     clik_input(n)%clik_cl_and_pars(j) = clik_nuis(5)   ! A_cib_217
                  else if (i .eq. 11) then
                     clik_input(n)%clik_cl_and_pars(j) = clik_nuis(9)   ! n_Dl_cib                 
                  else if (i .eq. 15) then
                     clik_input(n)%clik_cl_and_pars(j) = clik_nuis(8)   ! r_cib
                  else
                     clik_input(n)%clik_cl_and_pars(j) = clik_nuis(m+offset)
                     m = m+1
       	          end if
! Gaussian prior on a_gs and a_ge
                  if (i .eq. 16) then ! a_gs = 0.4 \pm 0.2
                     clik_lnlike = clik_lnlike + &
                        & ((clik_input(n)%clik_cl_and_pars(j) - 0.4d0)/0.2d0)**2
                  end if       	       	
                  if (i .eq. 17) then ! a_ge = 0.8 \pm 0.2
                     clik_lnlike = clik_lnlike + &
                        & ((clik_input(n)%clik_cl_and_pars(j) - 0.8d0)/0.2d0)**2
       	          end if
                  j = j+1
               end do
             end if

          end if   
   
!   do i=1,clik_nnuis(n)
!     Print*,i,clik_input(n)%clik_cl_and_pars(size(clik_input(n)%clik_cl_and_pars)-clik_nnuis(n)+i)
!   end do

!Get - ln like needed by CosmoMC
          clik_lnlike = clik_lnlike - 1.d0*clik_compute(clikid(n),clik_input(n)%clik_cl_and_pars)

       end if

    end do

    Print*,'clik lnlike = ',clik_lnlike

  end function clik_lnlike



  function clik_lensing_lnlike(cl,n,clik_nuis)
    real(dp) :: clik_lensing_lnlike
    real(dp) :: cl(lmax,num_cls_tot)
    real(dp), intent(in), optional :: clik_nuis(1:num_freq_params-1)
    integer :: n

    !set C_l and parameter vector to zero initially
    clik_input(n)%clik_cl_and_pars = 0.d0

    j = 1


    do l=0,clik_lensing_lmax(n)
       !skip C_0 and C_1
       if (l .ge. 2) then
          clik_input(n)%clik_cl_and_pars(j) = cl(l,num_cls+1)/(dble(l*(l+1)))**2*twopi
       end if
       j = j+1
    end do
     

    do l=0,clik_lensing_lmax(n)
       !skip C_0 and C_1
       if (l .ge. 2) then
          clik_input(n)%clik_cl_and_pars(j) = cl(l,mapped_index(1))
       end if
       j = j+1
    end do

    if (clik_nnuis(n) .ne. 0) then
       Print*,'Lensing nuisance parameters not implemented yet.'
       call MPIstop
    end if


!Get - ln like needed by CosmoMC
     clik_lensing_lnlike = - 1.d0*clik_lensing_compute(clikid(n),clik_input(n)%clik_cl_and_pars)

    Print*,'clik lensing lnlike = ',clik_lensing_lnlike

  end function clik_lensing_lnlike

end module cliklike
