    module pol_like_camspec
    !    use settings, only: MPIrank
    !use AMLutils

    implicit none

    private
    integer, parameter :: campc = KIND(1.d0)

    real(campc), dimension(:), allocatable :: X_data, X_theory, Y
    real(campc),  dimension(:,:), allocatable :: c_inv
    integer :: Nspec,nX,num_ells,nXfromdiff, CAMspec_lmax
    integer, dimension(:), allocatable  :: lminX, lmaxX, np, npt

    public pol_like_init,pol_calc_like

    contains

    subroutine pol_like_init(like_file)
    use MatrixUtils
    integer :: i, j,k,l
    logical :: pre_marged
    character(LEN=1024), intent(in) :: like_file
    logical, save :: needinit=.true.
    
    if(.not. needinit) return

    open(48, file=like_file, form='unformatted', status='unknown')

    read(48) Nspec,nX
    allocate(lminX(Nspec))
    allocate(lmaxX(Nspec))
    allocate(np(Nspec))
    allocate(npt(Nspec))

    read(48) (lminX(i), lmaxX(i), np(i), npt(i), i = 1, Nspec)
        allocate(X_data(nX))
        allocate(X_theory(nX))
        allocate(Y(nX))
        allocate(c_inv(nX,nX))
        read(48) (X_data(i), i=1, nX)
        read(48) !skip  covariance
        read(48) ((c_inv(i, j), j = 1, nX), i = 1,  nX) !inv covariance
        close(48)
        print *,' nX ',nX
        do i=1,Nspec
            print *, 'spec ',i, 'lmin,lmax, start ix = ', lminX(i),lmaxX(i), npt(i)
        end do

    needinit=.false.

    end subroutine pol_like_init


    subroutine pol_calc_like(zlike,  cell_cmbx,cell_cmbe)
    real(campc), dimension(0:) :: cell_cmbx, cell_cmbe
    integer ::  i, j, l, ii,jj,kk
    real(campc) zlike

    if (.not. allocated(lminX)) then
        print*, 'like_init should have been called before attempting to call calc_like.'
        stop
    end if
    if(Nspec.ne.2) then
        print*, 'Nspec inconsistent in pol calc_like.'
        stop
    end if

    X_theory=0.d0

      do kk = 1, Nspec
      do l = lminX(kk), lmaxX(kk)
       if(kk.eq.1) then
       X_theory(l - lminX(kk) + npt(kk)) = cell_cmbx(l)
       end if
       if(kk.eq.2) then
       X_theory(l - lminX(kk) + npt(kk)) = cell_cmbe(l)
       end if

      end do
      end do


    Y = X_data - X_theory

      zlike = 0.d+00
       do  i = 1, nX
         do  j = 1, nX
           zlike = zlike + Y(i)*Y(j)*c_inv(j, i)
         end do
       end do

    end subroutine pol_calc_like


    end module pol_like_camspec
