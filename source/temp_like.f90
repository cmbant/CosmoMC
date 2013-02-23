    module temp_like
    !    use settings, only: MPIrank
    !use AMLutils

    implicit none

    integer, parameter :: campc = KIND(1.d0)

    real(campc), dimension(:), allocatable :: X_data
    real(campc),  dimension(:,:), allocatable :: c_inv
    integer :: Nspec,nX,num_ells,nXfromdiff, CAMspec_lmax
    integer, dimension(:), allocatable  :: lminX, lmaxX, np, npt
    integer, parameter :: lmax_sz = 5000
    real(campc) :: sz_143_temp(lmax_sz)
    real(campc) :: ksz_temp(lmax_sz), tszxcib_temp(lmax_sz)

    real(campc), dimension(:,:), allocatable :: beam_cov,beam_cov_inv
    real(campc), dimension(:,:,:), allocatable :: beam_modes ! mode#, l, spec#
    integer :: num_modes_per_beam,beam_lmax,beam_Nspec,cov_dim

    logical :: writebest=.false.
    real(campc) :: bestlike = 1e30_campc
    integer :: tempinstance
    logical :: beam_cov_marge = .true.

    integer, allocatable :: indices(:),indices_reverse(:),keep_indices(:),keep_indices_reverse(:)
    real(campc), allocatable :: beam_conditional_mean(:,:)
    logical, allocatable :: want_marge(:)
    integer marge_num, keep_num

    logical :: storeall=.false.
    integer :: countnum
    !    character*100 storeroot,storename,storenumstring
    !   character*100 :: bestroot,bestnum,bestname
    character(LEN=*), parameter :: CAMSpec_like_version = '6.1/6.2-11012013'

    contains

    !!Does not seem to be actually faster
    !function CAMSpec_Quad(Mat,vec)
    !real(campc), intent(in) :: Mat(:,:)
    !real(campc) vec(:), CAMSpec_Quad
    !real(campc) Outv(nx)
    !real(campc)  mult, beta
    !real(campc) ddot
    !external ddot
    !
    !
    !mult = 1
    !beta = 0
    !call DSYMV('U',nX,mult,Mat,nX,vec, 1,beta, Outv,1)
    !CAMSpec_Quad = ddot(nX, vec, 1, outv, 1)
    !
    !end function CAMSpec_Quad

    subroutine CAMspec_ReadNormSZ(fname, templt)
    character(LEN=*), intent(in) :: fname
    real(campc) :: templt(lmax_sz)
    integer i, dummy
    real(campc) :: renorm

    open(48, file=fname, form='formatted', status='unknown')
    do i=2,lmax_sz
        read(48,*) dummy,templt(i)
        if (dummy/=i) stop 'CAMspec_ReadNormSZ: inconsistency in file read'
    enddo
    close(48)

    renorm=1.d0/templt(3000)
    templt=templt*renorm
    end subroutine CAMspec_ReadNormSZ

    subroutine like_init(like_file, sz143_file, tszxcib_file, ksz_file, beam_file)
    use MatrixUtils
    integer :: i, j,k,l
    character(LEN=1024), intent(in) :: like_file, sz143_file, ksz_file, tszxcib_file, beam_file
    logical, save :: needinit=.true.
    real(campc) , allocatable:: fid_cl(:), beam_cov(:,:), beam_cov_full(:,:)
    integer if1,if2, ie1, ie2, ii, jj, L2

    ! cl_ksz_148_tbo.dat file is in D_l, format l D_l, from l=2 to 10000
    ! tsz_x_cib_template.txt is is (l D_l), from l=2 to 9999, normalized to unity
    !    at l=3000

    if(.not. needinit) return

    open(48, file=like_file, form='unformatted', status='unknown')

    read(48) Nspec,nX
    allocate(lminX(Nspec))
    allocate(lmaxX(Nspec))
    allocate(np(Nspec))
    allocate(npt(Nspec))
    allocate(X_data(nX))
    allocate(c_inv(nX,nX))

    read(48) (lminX(i), lmaxX(i), np(i), npt(i), i = 1, Nspec)
    read(48) (X_data(i), i=1, nX)
    read(48) ((c_inv(i, j), j = 1, nX), i = 1,  nX) !covarianbce
    read(48) !((c_inv(i, j), j = 1, nX), i = 1,  nX) !inver covariuance
    close(48)

    CAMspec_lmax = maxval(lmaxX)
    allocate(fid_cl(CAMspec_lmax))

    open(48, file='./data/base_planck_CAMspec_lowl_lowLike_highL.bestfit_cl', form='formatted', status='unknown')
    do i=2,CAMspec_lmax
        read(48,*) j,fid_cl(i)
        if (j/=i) stop 'error reading fiducial C_l for beams'
        fid_cl(i) = fid_cl(i)/(i*(i+1))!want C_l/2Pi
    enddo
    close(48)

    call CAMspec_ReadNormSZ(sz143_file, sz_143_temp)
    call CAMspec_ReadNormSZ(ksz_file, ksz_temp)
    call CAMspec_ReadNormSZ(tszxcib_file, tszxcib_temp)

    open(48, file=beam_file, form='unformatted', status='unknown')
    read(48) beam_Nspec,num_modes_per_beam,beam_lmax
    if(beam_Nspec.ne.Nspec) stop 'Problem: beam_Nspec != Nspec'
    allocate(beam_modes(num_modes_per_beam,0:beam_lmax,beam_Nspec))
    cov_dim=beam_Nspec*num_modes_per_beam
    allocate(beam_cov_inv(cov_dim,cov_dim))
    allocate(beam_cov_full(cov_dim,cov_dim))
    read(48) (((beam_modes(i,l,j),j=1,Nspec),l=0,beam_lmax),i=1,num_modes_per_beam)
    read(48) ((beam_cov_full(i,j),j=1,cov_dim),i=1,cov_dim)  ! beam_cov
    read(48) ((beam_cov_inv(i,j),j=1,cov_dim),i=1,cov_dim) ! beam_cov_inv
    close(48)

    allocate(want_marge(cov_dim))
    want_marge=.true.
    want_marge(1)=.false.
    marge_num=count(want_marge)
    keep_num=cov_dim-marge_num
    allocate(indices(marge_num))
    allocate(indices_reverse(cov_dim))
    allocate(keep_indices(keep_num))
    allocate(keep_indices_reverse(cov_dim))
    print *,'beam marginalizing:',marge_num,'keeping',keep_num

    j=0
    k=0
    do i=1,cov_dim
        if (want_marge(i)) then
            j=j+1
            indices(j) = i
        else
            k=k+1
            keep_indices(k)=i
        end if
        indices_reverse(i)=j
        keep_indices_reverse(i)=k
    end do
    allocate(beam_cov(marge_num, marge_num))
    beam_cov = beam_cov_inv(indices,indices)
    call Matrix_Inverse(beam_cov)

    if (beam_cov_marge) then 
        do if2=1,beam_Nspec
            do if1=1,beam_Nspec
                do ie2=1,num_modes_per_beam
                    do ie1=1,num_modes_per_beam
                        ii=ie1+num_modes_per_beam*(if1-1)
                        jj=ie2+num_modes_per_beam*(if2-1)
                        if (want_marge(ii) .and. want_marge(jj)) then
                            do L2 = lminX(if2),lmaxX(if2)
                                c_inv(npt(if1):npt(if1)+lmaxX(if1)-lminX(if1),npt(if2)+L2 -lminX(if2) ) = &
                                c_inv(npt(if1):npt(if1)+lmaxX(if1)-lminX(if1),npt(if2)+L2 -lminX(if2) ) + &
                                beam_modes(ie1,lminX(if1):lmaxX(if1),if1)* &
                                beam_cov(indices_reverse(ii),indices_reverse(jj)) * beam_modes(ie2,L2,if2) &
                                *fid_cl(L2)*fid_cl(lminX(if1):lmaxX(if1))
                            end do
                        end if
                    enddo
                enddo
            enddo
        enddo
    end if
    ! print *,'after', c_inv(500,500), c_inv(npt(3)-500,npt(3)-502),  c_inv(npt(4)-1002,npt(4)-1000)

    allocate(beam_conditional_mean(marge_num, keep_num))
    beam_conditional_mean=-matmul(beam_cov, beam_cov_inv(indices,keep_indices))
    deallocate(beam_cov_inv)
    allocate(beam_cov_inv(keep_num,keep_num))
    beam_cov_inv = beam_cov_full(keep_indices,keep_indices)
    call Matrix_inverse(beam_cov_inv)

    call Matrix_inverse(c_inv)

    !tempinstance=MPIrank+1
    !
    !write(bestnum,*) tempinstance
    !bestroot='/home/stg20/src/cosmomc/chains/best'
    !bestname=trim(bestroot)//'_'//trim(adjustl(bestnum))//'.txt'
    !storeroot='/home/stg20/src/cosmomc/chains/store'
    countnum=0

    needinit=.false.

    end subroutine like_init


    subroutine calc_like(zlike,  cell_cmb, freq_params)
    real(campc), intent(in)  :: freq_params(:)
    real(campc), dimension(0:) :: cell_cmb
    integer ::  j, l, ii,jj
    real(campc) , allocatable, save ::  X_beam_corr_model(:), Y(:),  C_foregrounds(:,:)
    real(campc) zlike
    real(campc) A_ps_100, A_ps_143, A_ps_217, A_cib_143, A_cib_217, A_sz, r_ps, r_cib, ncib, cal0, cal1, cal2, xi, A_ksz
    real(campc) zCIB
    real(campc) ztemp
    real(campc) beam_params(cov_dim),beam_coeffs(beam_Nspec,num_modes_per_beam)
    integer :: ie1,ie2,if1,if2, ix
    integer num_non_beam
    real(campc) cl_cib(CAMspec_lmax) !CIB
    real(campc), parameter :: sz_bandpass100_nom143 = 2.022d0
    real(campc), parameter :: cib_bandpass143_nom143 = 1.134d0
    real(campc), parameter :: sz_bandpass143_nom143 = 0.95d0
    real(campc), parameter :: cib_bandpass217_nom217 = 1.33d0
    real(campc), parameter :: ps_scale  = 1.d-6/9.d0
    real(campc) :: A_cib_217_bandpass, A_sz_143_bandpass, A_cib_143_bandpass

    !    real(campc) atime

    if (.not. allocated(lminX)) then
        print*, 'like_init should have been called before attempting to call calc_like.'
        stop
    end if
    if(Nspec.ne.4) then
        print*, 'Nspec inconsistent with foreground corrections in calc_like.'
        stop
    end if
    if (.not. allocated(Y)) then
        allocate(X_beam_corr_model(1:nX))
        allocate(Y(1:nX))
        allocate(C_foregrounds(CAMspec_lmax,Nspec))
        C_foregrounds=0
    end if

    ! atime = MPI_Wtime()

    num_non_beam = 14
    if (size(freq_params) < num_non_beam +  beam_Nspec*num_modes_per_beam) stop 'CAMspec: not enough parameters'
    A_ps_100=freq_params(1)
    A_ps_143 = freq_params(2)
    A_ps_217 = freq_params(3)
    A_cib_143 =freq_params(4)
    A_cib_217 =freq_params(5) 
    A_sz = freq_params(6)  !143
    r_ps = freq_params(7)
    r_cib = freq_params(8)
    ncib = freq_params(9)
    cal0 = freq_params(10)
    cal1 = freq_params(11) 
    cal2 = freq_params(12)
    xi = freq_params(13)
    A_ksz = freq_params(14)


    beam_params = freq_params(num_non_beam+1:num_non_beam+cov_dim)
    !set marged beam parameters to their mean subjected to fixed non-marged modes
    beam_params(indices) = matmul(beam_conditional_mean, freq_params(num_non_beam+keep_indices))

    do ii=1,beam_Nspec
        do jj=1,num_modes_per_beam
            ix = jj+num_modes_per_beam*(ii-1)
            beam_coeffs(ii,jj)=beam_params(ix)
        enddo
    enddo

    if (beam_cov_marge) beam_coeffs= 0

    do l=1, CAMspec_lmax
        cl_cib(l) = (real(l,campc)/3000)**(ncib)
    end do

    !   100 foreground
    !
    do l = lminX(1), lmaxX(1)
        C_foregrounds(l,1)= A_ps_100*ps_scale+  &
        ( A_ksz*ksz_temp(l) + A_sz*sz_bandpass100_nom143*sz_143_temp(l) )/(l*(l+1))
        X_beam_corr_model(l-lminX(1)+1) = ( cell_cmb(l) + C_foregrounds(l,1) )* corrected_beam(1,l)/cal0
    end do

    !   143 foreground
    !
    A_sz_143_bandpass = A_sz * sz_bandpass143_nom143
    A_cib_143_bandpass = A_cib_143 * cib_bandpass143_nom143
    do l = lminX(2), lmaxX(2)
        zCIB = A_cib_143_bandpass*cl_cib(l)
        C_foregrounds(l,2)= A_ps_143*ps_scale + &
        (zCIB +  A_ksz*ksz_temp(l) + A_sz_143_bandpass*sz_143_temp(l) &
        -2.0*sqrt(A_cib_143_bandpass * A_sz_143_bandpass)*xi*tszxcib_temp(l) )/(l*(l+1))
        X_beam_corr_model(l-lminX(2)+npt(2)) =  (cell_cmb(l)+ C_foregrounds(l,2))*corrected_beam(2,l)/cal1
    end do

    !   217 foreground
    !
    A_cib_217_bandpass = A_cib_217 * cib_bandpass217_nom217
    do l = lminX(3), lmaxX(3)
        zCIB = A_cib_217_bandpass*cl_cib(l)
        C_foregrounds(l,3) = A_ps_217*ps_scale + (zCIB + A_ksz*ksz_temp(l) )/(l*(l+1))
        X_beam_corr_model(l-lminX(3)+npt(3)) = (cell_cmb(l)+ C_foregrounds(l,3))* corrected_beam(3,l)/cal2
    end do

    !   143x217 foreground
    !
    do l = lminX(4), lmaxX(4)
        zCIB = sqrt(A_cib_143_bandpass*A_cib_217_bandpass)*cl_cib(l)
        C_foregrounds(l,4) = r_ps*sqrt(A_ps_143*A_ps_217)*ps_scale + &
        ( r_cib*zCIB + A_ksz*ksz_temp(l) -sqrt(A_cib_217_bandpass * A_sz_143_bandpass)*xi*tszxcib_temp(l) )/(l*(l+1))
        X_beam_corr_model(l-lminX(4)+npt(4)) =  ( cell_cmb(l) + C_foregrounds(l,4))*corrected_beam(4,l)/sqrt(cal1*cal2)
    end do

    Y = X_data - X_beam_corr_model

    zlike = 0
    !$OMP parallel do private(j,ztemp) reduction(+:zlike) schedule(static,16)
    do  j = 1, nX
        ztemp= dot_product(Y(j+1:nX), c_inv(j+1:nX, j))
        zlike=zlike+ (ztemp*2 +c_inv(j, j)*Y(j))*Y(j)
    end do
    !    zlike = 0
    !    do  j = 1, nX
    !       ztemp= 0
    !       do  i = 1, nX
    !          ztemp = ztemp + Y(i)*c_inv(i, j)
    !       end do
    !       zlike=zlike+ztemp*Y(j)
    !    end do
    !   zlike = CAMSpec_Quad(c_inv, Y)

    do if2=1,beam_Nspec
        do if1=1,beam_Nspec
            do ie2=1,num_modes_per_beam
                jj=ie2+num_modes_per_beam*(if2-1)
                if (.not. want_marge(jj)) then
                    do ie1=1,num_modes_per_beam
                        ii=ie1+num_modes_per_beam*(if1-1)
                        if (.not. want_marge(ii)) then
                            zlike=zlike+beam_coeffs(if1,ie1)*&
                            beam_cov_inv(keep_indices_reverse(ii),keep_indices_reverse(jj))*beam_coeffs(if2,ie2)
                        end if
                    enddo
                end if
            enddo
        enddo
    enddo

    zlike=zlike+((cal2/cal1-0.9966d0)/0.0015d0)**2  +((cal0/cal1-1.0006d0)/0.0004d0)**2

    ! print *,'CAMspec time:',  MPI_Wtime() - atime

    ! WARNING; weakening prior...
    !    zlike=zlike+((cal2/cal1-1.0056d0)/0.0063d0)**2 &
    !               +((cal0/cal1-1.0127d0)/0.003d0)**2


    ! correction to improve old way of handling cals...
    !zlike=zlike &
    !-2.d0*(lmaxX(1)-lminX(1)+1)*log(cal0) &
    !-2.d0*(lmaxX(2)-lminX(2)+1)*log(cal1) &
    !-2.d0*(lmaxX(3)-lminX(3)+1)*log(cal2) &
    !-2.d0*(lmaxX(4)-lminX(4)+1)*.5d0*(log(cal1)+log(cal2))


    !    if(writebest) then
    !    if(zlike.lt.bestlike) then
    !       open(49,file=bestname,form='formatted',status='unknown')
    !       !do if1=1,beam_Nspec
    !       !write(49,101) (beam_coeffs(if1,ie1),ie1=1,num_modes_per_beam)
    !       !enddo
    !101    format(5(e12.4))
    !       do l = 0,2500
    !          write(49,*) l, l*(l+1)*cell_cmb(l)
    !       enddo
    !       close(49)
    !       bestlike=zlike
    !    endif
    !    endif

    !if(storeall) then
    !countnum=countnum+1
    !write(storenumstring,*) countnum
    !storename=trim(storeroot)//'_'//trim(adjustl(bestnum))//'_' &
    !     //trim(adjustl(storenumstring))//'.txt'
    !   open(49,file=storename,form='formatted',status='unknown')
    !   do l=0,2500
    !   write(49,*) l,l*(l+1)*cell_cmb(l)
    !   enddo
    !   write(49,*)
    !   write(49,*) 'A_ps_100=',A_ps_100
    !   write(49,*) 'A_ps_143=',A_ps_143
    !   write(49,*) 'A_ps_217=',A_ps_217
    !   write(49,*) 'A_cib_143=',A_cib_143
    !   write(49,*) 'A_cib_217=',A_cib_217
    !   write(49,*) 'A_sz=',A_sz
    !   write(49,*) 'r_ps=',r_ps
    !   write(49,*) 'r_cib=',r_cib
    !   write(49,*) 'xi=',xi
    !   write(49,*) 'A_ksz=',A_ksz
    !   write(49,*) 'cal0=',cal0
    !   write(49,*) 'cal1=',cal1
    !   write(49,*) 'cal2=',cal2
    !   do i=1,4
    !   do j=1,num_modes_per_beam
    !      write(49,*) 'beam_coeffs(',i,',',j,')=',beam_coeffs(i,j)
    !      enddo
    !      enddo
    !      write(49,*)
    !      write(49,*) 'zlike(=chi^2)=',zlike
    !
    !
    !endif
    contains

    real(campc) function corrected_beam(spec_num,l)
    integer, intent(in) :: spec_num,l
    integer :: i

    corrected_beam=1.d0
    do i=1,num_modes_per_beam
            corrected_beam=corrected_beam+beam_coeffs(spec_num,i)*beam_modes(i,l,spec_num)
    enddo
    end function corrected_beam

    end subroutine calc_like


    end module temp_like
