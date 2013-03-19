!Module from Jan Hamann, 4/2010
!Modified by AL, 4/2010
module bbn

  use settings

  implicit none
  private

  integer, parameter :: dp = KIND(1.d0)

  type bbnstuff
     real(dp), dimension(:), pointer :: ombh2,deltan
     real(dp), dimension(:,:), pointer :: yp,ddyp
     integer :: n_ombh2,n_deltan
  end type

  type(bbnstuff) bbn_data

  public bbnini,yp_bbn

contains

  subroutine bbnini
    real(dp), dimension(:,:), allocatable :: bbn_tmp
! Number of grid points in \omega_b h^2 and \Delta N_\nu
    integer, parameter :: num_ombh2 = 26, num_deltan = 11
    integer :: i,j, file_id
    character :: dummy


    if (feedback .ge. 1) print*,'Initialising BBN Helium data...'

    allocate(bbn_tmp(num_ombh2*num_deltan,3))
    allocate(bbn_data%ombh2(num_ombh2))
    allocate(bbn_data%deltan(num_deltan))
    allocate(bbn_data%yp(num_ombh2,num_deltan))
    allocate(bbn_data%ddyp(num_ombh2,num_deltan))

    bbn_data%n_ombh2 = num_ombh2
    bbn_data%n_deltan = num_deltan

   file_id = new_file_unit()
   call OpenTxtFile(trim(DataDir)//'helium.dat', file_id)

!skip data file header
    do i=1,7
       read(file_id,*) dummy
    end do

!read in data
    do i = 1,num_ombh2*num_deltan
       read(file_id,*) bbn_tmp(i,1),bbn_tmp(i,2),bbn_tmp(i,3)
    end do

    call CloseFile(file_id)

    do i = 1,num_ombh2
       bbn_data%ombh2(i) = bbn_tmp(i,1)

       do j = 1,num_deltan
          bbn_data%yp(i,j) = bbn_tmp(num_ombh2*(j-1)+i,3)
       end do
    end do

    do i = 1,num_deltan
       bbn_data%deltan(i) = bbn_tmp((i-1)*num_ombh2+1,2)
    end do

    deallocate(bbn_tmp)

!prepare array of second derivatives needed for spline interpolation
    call bbn_splie2(bbn_data%ombh2,bbn_data%deltan,bbn_data%yp,num_ombh2,num_deltan,bbn_data%ddyp)

  end subroutine bbnini


!interpolate Yp from grid
  real(mcp) function yp_bbn(ombh2_in,deltan_in)
    real(mcp) :: ombh2_in,deltan_in
    logical, save:: initialized = .false.
    real(dp) res

    if (.not. initialized) then
     call bbnini
     initialized = .true.
    end if


!Make sure ombh2 and deltan are within ranges covered by the data file
    if ((ombh2_in .gt. bbn_data%ombh2(size(bbn_data%ombh2))) .or. (ombh2_in .lt. bbn_data%ombh2(1))) then
       Print ('(A,F6.3,A,F6.3,A)'),' Baryon density must be between',bbn_data%ombh2(1) &
            & ,' and',bbn_data%ombh2(size(bbn_data%ombh2)) &
            & ,' if you want to use the BBN consistency relation. Stopping.'
       call MpiStop
    end if
    if ((deltan_in .gt. bbn_data%deltan(size(bbn_data%deltan))) .or. (deltan_in .lt. bbn_data%deltan(1))) then
       Print ('(A,F6.3,A,F7.3,A)'),' Effective number of neutrinos must be between'&
            ,(bbn_data%deltan(1)+3.046),' and',(bbn_data%deltan(size(bbn_data%deltan))+3.046) &
            & ,' if you want to use the BBN consistency relation. Stopping.'
       call MpiStop
    end if


    call bbn_splin2(bbn_data%ombh2,bbn_data%deltan,bbn_data%yp,bbn_data%ddyp,bbn_data%n_ombh2, &
                    bbn_data%n_deltan,real(ombh2_in,dp),real(deltan_in,dp),res)
    yp_bbn = res

  end function yp_bbn



!spline and splint routines, adapted from Numerical Recipes
  subroutine bbn_spline(x,y,n,yp1,ypn,y2)
    integer :: n,i,k
    real(dp) :: yp1,ypn,x(n),y(n),y2(n)
    integer, parameter :: NMAX=500
    real(dp) :: p,qn,sig,un,u(NMAX)

    if (yp1.gt..99e30) then
       y2(1)=0.
       u(1)=0.
    else
       y2(1)=-0.5
       u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
    endif
    do i=2,n-1
       sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
       p=sig*y2(i-1)+2.
       y2(i)=(sig-1.)/p
       u(i)=(6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
     end do
       if (ypn.gt..99e30) then
          qn=0.
          un=0.
       else
          qn=0.5
          un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
       endif
       y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)
       do k=n-1,1,-1
          y2(k)=y2(k)*y2(k+1)+u(k)
       end do
  end subroutine bbn_spline


  subroutine bbn_splie2(x1a,x2a,ya,m,n,y2a)
    integer :: m,n,j,k
    real(dp) :: x1a(m),x2a(n),y2a(m,n),ya(m,n)
    integer, parameter :: NN=100
    real(dp) :: y2tmp(NN),ytmp(NN)

    do j=1,m
       do k=1,n
          ytmp(k)=ya(j,k)
       end do
       call bbn_spline(x2a,ytmp,n,1.d30,1.d30,y2tmp)
       do k=1,n
             y2a(j,k)=y2tmp(k)
       end do
    end do
  end subroutine bbn_splie2


  subroutine bbn_splint(xa,ya,y2a,n,x,y)
    integer :: n,k,khi,klo
    real(dp) :: x,y,xa(n),y2a(n),ya(n),a,b,h

!safeguard against extrapolation added
    if (.not. (((x.ge.xa(n)) .and. (x.le.xa(1))) .or. ((x.ge.xa(1)) .and. (x.lt.xa(n))))) then
       Print*,'Input of bbn_splint out of interpolation range.'
       Print*,xa(n),x,xa(1)
       stop
    end if
    klo=1
    khi=n
    do while (khi-klo.gt.1)
       k=(khi+klo)/2
       if(xa(k).gt.x)then
          khi=k
       else
          klo=k
       endif
    end do
    h=xa(khi)-xa(klo)
    if (h.eq.0.) call MpiStop('bad xa input in bbn_splint')
    a=(xa(khi)-x)/h
    b=(x-xa(klo))/h
    y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.
  end subroutine bbn_splint


  subroutine bbn_splin2(x1a,x2a,ya,y2a,m,n,x1,x2,y)
    integer :: m,n,j,k
    real(dp) :: x1,x2,y,x1a(m),x2a(n),y2a(m,n),ya(m,n)
    integer, parameter :: NN=100
    real(dp) :: y2tmp(NN),ytmp(NN),yytmp(NN)

    do j=1,m
       do k=1,n
          ytmp(k)=ya(j,k)
          y2tmp(k)=y2a(j,k)
       end do
      call bbn_splint(x2a,ytmp,y2tmp,n,x2,yytmp(j))
    end do
    call bbn_spline(x1a,yytmp,m,1.d30,1.d30,y2tmp)
    call bbn_splint(x1a,yytmp,y2tmp,m,x1,y)
  end subroutine bbn_splin2


end module bbn

