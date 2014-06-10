    module temp_like_camspec
    implicit none

    private
    integer, parameter :: campc = KIND(1.d0)

    real(campc), dimension(:), allocatable :: X_data
    real(campc),  dimension(:,:), allocatable :: c_inv
    integer :: Nspec,nX,num_ells,nXfromdiff, CAMspec_lmax
    integer, dimension(:), allocatable  :: lminX, lmaxX, npt
    integer, parameter :: CAMSpec_lmax_foreground = 5000
    integer, parameter :: lmax_sz = CAMSpec_lmax_foreground
    real(campc) :: sz_143_temp(lmax_sz)
    real(campc) :: ksz_temp(lmax_sz), tszxcib_temp(lmax_sz)

    real(campc), dimension(:,:), allocatable :: beam_cov,beam_cov_inv
    real(campc), dimension(:,:,:), allocatable :: beam_modes ! mode#, l, spec#
    integer :: num_modes_per_beam,beam_lmax,beam_Nspec,cov_dim

    integer, allocatable :: marge_indices(:),marge_indices_reverse(:)
    integer, allocatable :: keep_indices(:),keep_indices_reverse(:)
    real(campc), allocatable :: beam_conditional_mean(:,:)
    logical, allocatable :: want_marge(:)
    integer marge_num, keep_num

    logical :: make_cov_marged = .false.
    real(campc) :: beam_factor = 2.7_campc
    integer, parameter :: CamSpec_cib_pivot = 3000
    integer, parameter :: CamSpec_sz_pivot = 3000

    logical :: want_spec(6) = .true.
    logical :: camspec_has_TT = .true.
    integer :: camspec_lmins(6) =0
    integer :: camspec_lmaxs(6) =0
    real(campc)  :: llp1(CAMSpec_lmax_foreground)

    !paramters for systematics testing
    integer :: camspec_nwiggles = 0

    integer :: camspec_beam_mcmc_num = 1

    character(LEN=1024) camspec_fiducial_foregrounds, camspec_fiducial_cl
    character(LEN=80) :: marge_file_variant = ''

    character(LEN=*), parameter :: CAMSpec_like_version = 'CamSpec_v2_wig'
    public like_init,calc_like,CAMSpec_like_version, camspec_beam_mcmc_num, &
        want_spec,camspec_lmins,camspec_lmaxs, make_cov_marged, marge_file_variant,&
        compute_fg, Nspec,CAMSpec_lmax_foreground,camspec_fiducial_foregrounds,camspec_fiducial_cl
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

    open(48, file=fname, form='formatted', status='old')
    do i=2,lmax_sz
        read(48,*) dummy,templt(i)
        if (dummy/=i) stop 'CAMspec_ReadNormSZ: inconsistency in file read'
    enddo
    close(48)

    renorm=1.d0/templt(CamSpec_sz_pivot)
    templt=templt*renorm
    end subroutine CAMspec_ReadNormSZ

    subroutine ReadFiducialCl(fid_cl)
    integer i,j
    real(campc) :: fid_theory_tt,fid_theory_te,fid_theory_ee
    real(campc) :: fid_cl(:,:)

    open(48, file=CAMspec_fiducial_foregrounds, form='formatted', status='old')
    do i=2,CAMspec_lmax
        read(48,*) j, fid_cl(i,1:4)
        if(Nspec.eq.6) then
            fid_cl(i,5)=0.0
            fid_cl(i,6)=0.0
        endif
        if (j/=i) stop 'error reading fiducial foreground C_l for beams'
    enddo
    close(48)
    open(48, file=CAMspec_fiducial_cl, form='formatted', status='old')
    do i=2,CAMspec_lmax
        read(48,*, end=100) j,fid_theory_tt,fid_theory_te,fid_theory_ee
100     if (j/=i) stop 'error reading fiducial C_l for beams'
        if(Nspec.eq.4) then
            fid_cl(i,:) = (fid_cl(i,:) + fid_theory_tt)*llp1(i) !want C_l/2Pi foregrounds+theory
        else if(Nspec.eq.6) then
            fid_cl(i,1:4) = (fid_cl(i,1:4) + fid_theory_tt)*llp1(i) !want C_l/2Pi foregrounds+theory
            fid_cl(i,5) = (fid_cl(i,5) + fid_theory_te)*llp1(i) !want C_l/2Pi foregrounds+theory
            fid_cl(i,6) = (fid_cl(i,6) + fid_theory_ee)*llp1(i) !want C_l/2Pi foregrounds+theory
        endif
    enddo
    close(48)

    end subroutine ReadFiducialCl

    subroutine like_init(pre_marged,like_file, sz143_file, tszxcib_file, ksz_file, beam_file,data_vector)
    use MatrixUtils
    integer :: i, j,k,l
    logical :: pre_marged
    character(LEN=*), intent(in) :: like_file, sz143_file, ksz_file, tszxcib_file, beam_file,data_vector
    logical, save :: needinit=.true.
    real(campc) , allocatable:: fid_cl(:,:), beam_cov(:,:), beam_cov_full(:,:)
    integer if1,if2, ie1, ie2, ii, jj, L2
    real(campc), dimension(:), allocatable :: X_data_in
    real(campc),  dimension(:,:), allocatable :: cov
    integer, allocatable :: indices(:),np(:)
    integer ix, status
    real(campc) dummy, in_data(4)
    real(campc), allocatable :: CL_in(:,:)

    if(.not. needinit) return

    do i=1,CAMSpec_lmax_foreground
        llp1(i) = 1/real(i*(i+1),campc)
    end do

    camspec_has_TT = any(want_spec(1:4))

    open(48, file=like_file, form='unformatted', status='old')

    read(48) Nspec,nX
    allocate(lminX(Nspec))
    allocate(lmaxX(Nspec))
    allocate(np(Nspec))
    allocate(npt(Nspec))

    read(48) (lminX(i), lmaxX(i), np(i), npt(i), i = 1, Nspec)
    if (pre_marged) then
        print *,'Using pre-beam-marged covariance; L ranges etc ignored'
        allocate(X_data(nX))
        allocate(c_inv(nX,nX))
        read(48) (X_data(i), i=1, nX)
        read(48) !skip  covariance
        read(48) ((c_inv(i, j), j = 1, nX), i = 1,  nX) !inv covariance
        close(48)
        do i=1,Nspec
            print *, 'spec ',i, 'lmin,lmax, start ix = ', lminX(i),lmaxX(i), npt(i)
        end do
    else
        !Chop L ranges/frequencies if requested
        allocate(indices(nX))
        allocate(cov(nX,nX))
        allocate(X_data_in(nX))
        read(48) X_data_in !(X_data(i), i=1, nX)
        read(48) cov !((c_inv(i, j), j = 1, nX), i = 1,  nX) !covariance
        read(48) !((c_inv(i, j), j = 1, nX), i = 1,  nX) !inver covariance
        close(48)

        !Cut on L ranges
        print *,'Determining L ranges'
        j=0
        ix=0
        npt(1)=1
        if(.not. want_spec(1) .and. marge_num>0) stop 'One beam mode may not be right here'
        do i=1,Nspec
            do l = lminX(i), lmaxX(i)
                j =j+1
                if (want_spec(i) .and. (camspec_lmins(i)==0 .or. l>=camspec_lmins(i)) &
                    .and. (camspec_lmaxs(i)==0 .or. l<=camspec_lmaxs(i)) ) then
                ix =ix+1
                indices(ix) = j
                end if
            end do
            if (want_spec(i)) then
                if (camspec_lmins(i)/=0) lminX(i) = max(camspec_lmins(i), lminX(i))
                if (camspec_lmaxs(i)/=0) lmaxX(i) = min(camspec_lmaxs(i), lmaxX(i))
            else
                lmaxX(i) =0
            end if
            if (i<NSpec) npt(i+1) = ix+1
            print *, 'spec ',i, 'lmin,lmax, start ix = ', lminX(i),lmaxX(i), npt(i)
        end do
        if (j/=nx) stop 'Error cutting camspec cov matrix'
        nX = ix
        allocate(X_data(nX))
        allocate(c_inv(nX,nX))
        do i=1, nX
            c_inv(:,i) = cov(indices(1:nX),indices(i)) !c_inv is covariance here
            X_data(i) = X_data_in(indices(i))
        end do
        !    call Matrix_inverse(c_inv)
        deallocate(cov)
    end if
    CAMspec_lmax = maxval(lmaxX)

    if (data_vector/='') then
        if (Nspec.eq.6) stop 'not yet modified for te+ee'
        !!override data vector in the main camspec file
        !The input file is in comprehensible L, L(L+1)C_L^i/2pi format
        open(48, file=data_vector, form='formatted', status='old', iostat=status)
        if (status/=0) stop 'Error opening camspec data_vector override'
        print *,'Using CamSpec data vector : '//data_vector
        allocate(CL_in(CAMspec_lmax, 4))
        L=0
        do
            read(48,*,iostat = status) L, in_data
            if (status/=0) exit
            CL_in(L,:) = in_data
        end do
        if (L< CAMspec_lmax) stop 'Error reading camspec data_vector override'
        close(48)
        ix = 0
        do i=1,NSpec
            do L=lminX(i), lmaxX(i)
                ix=ix+1
                X_data(ix) = CL_in(L,i)/(L*(L+1))
            end do
        end do
    end if

    if (.not. pre_marged) then
        allocate(fid_cl(CAMspec_lmax,Nspec))
        call ReadFiducialCl(fid_cl)
    end if

    call CAMspec_ReadNormSZ(sz143_file, sz_143_temp)
    call CAMspec_ReadNormSZ(ksz_file, ksz_temp)
    call CAMspec_ReadNormSZ(tszxcib_file, tszxcib_temp)

    open(48, file=beam_file, form='unformatted', status='unknown')
    read(48) beam_Nspec,num_modes_per_beam,beam_lmax
    ! removed until we decide about how to handle beam errors for TE & EE
    !    if(beam_Nspec.ne.Nspec) stop 'Problem: beam_Nspec != Nspec'
    allocate(beam_modes(num_modes_per_beam,0:beam_lmax,beam_Nspec))
    cov_dim=beam_Nspec*num_modes_per_beam
    allocate(beam_cov_full(cov_dim,cov_dim))
    read(48) (((beam_modes(i,l,j),j=1,beam_Nspec),l=0,beam_lmax),i=1,num_modes_per_beam)
    allocate(beam_cov_inv(cov_dim,cov_dim))
    read(48) ((beam_cov_full(i,j),j=1,cov_dim),i=1,cov_dim)  ! beam_cov
    read(48) ((beam_cov_inv(i,j),j=1,cov_dim),i=1,cov_dim) ! beam_cov_inv
    close(48)

    allocate(want_marge(cov_dim))
    want_marge=.true.
    want_marge(1:camspec_beam_mcmc_num)=.false.
    marge_num=count(want_marge)
    keep_num=cov_dim-marge_num
    allocate(marge_indices(marge_num))
    allocate(marge_indices_reverse(cov_dim))
    allocate(keep_indices(keep_num))
    allocate(keep_indices_reverse(cov_dim))

    j=0
    k=0
    do i=1,cov_dim
        if (want_marge(i)) then
            j=j+1
            marge_indices(j) = i
        else
            k=k+1
            keep_indices(k)=i
        end if
        marge_indices_reverse(i)=j
        keep_indices_reverse(i)=k
    end do

    allocate(beam_conditional_mean(marge_num, keep_num))
    if (marge_num>0) then
        allocate(beam_cov(marge_num, marge_num))
        beam_cov = beam_cov_inv(marge_indices,marge_indices)
        call Matrix_Inverse(beam_cov)
        beam_conditional_mean=-matmul(beam_cov, beam_cov_inv(marge_indices,keep_indices))
    end if
    deallocate(beam_cov_inv)

    if (pre_marged) then
        print *,'using marginalizing modes:',marge_num,'keeping',keep_num
    else
        print *,'beam marginalizing:',marge_num,'keeping',keep_num
        if (marge_num>0) then
            do if2=1,beam_Nspec
                if (want_spec(if2)) then
                    do if1=1,beam_Nspec
                        if (want_spec(if1)) then
                            do ie2=1,num_modes_per_beam
                                do ie1=1,num_modes_per_beam
                                    ii=ie1+num_modes_per_beam*(if1-1)
                                    jj=ie2+num_modes_per_beam*(if2-1)
                                    if (want_marge(ii) .and. want_marge(jj)) then
                                        do L2 = lminX(if2),lmaxX(if2)
                                            !c_inv is covariance here
                                            c_inv(npt(if1):npt(if1)+lmaxX(if1)-lminX(if1),npt(if2)+L2 -lminX(if2) ) = &
                                                c_inv(npt(if1):npt(if1)+lmaxX(if1)-lminX(if1),npt(if2)+L2 -lminX(if2) ) + &
                                                beam_factor**2*beam_modes(ie1,lminX(if1):lmaxX(if1),if1)* &
                                                beam_cov(marge_indices_reverse(ii),marge_indices_reverse(jj)) * beam_modes(ie2,L2,if2) &
                                                *fid_cl(L2,if2)*fid_cl(lminX(if1):lmaxX(if1),if1)
                                        end do
                                    end if
                                enddo
                            enddo
                        end if
                    enddo
                end if
            enddo
            ! print *,'after', c_inv(500,500), c_inv(npt(3)-500,npt(3)-502),  c_inv(npt(4)-1002,npt(4)-1000)

            call Matrix_inverse(c_inv) !now c_inv is the inverse covariance
            if (make_cov_marged .and. marge_num>0) then
                open(48, file=trim(like_file)//'_beam_marged'//trim(marge_file_variant), form='unformatted', status='unknown')
                write(48) Nspec,nX
                write(48) (lminX(i), lmaxX(i), np(i), npt(i), i = 1, Nspec)
                write(48) (X_data(i), i=1, nX)
                dummy=-1
                write(48) dummy !inver covariance, assume not used
                write(48) ((c_inv(i, j), j = 1, nX), i = 1,  nX) !inver covariuance
                close(48)
                write(*,*) 'Marged file done: '//trim(like_file)//'_beam_marged'//trim(marge_file_variant)
                stop
            end if
        else
            call Matrix_inverse(c_inv)
        end if
    end if !pre-marged

    if (keep_num>0) then
        allocate(beam_cov_inv(keep_num,keep_num))
        beam_cov_inv = beam_cov_full(keep_indices,keep_indices)
        call Matrix_inverse(beam_cov_inv)
    end if

    needinit=.false.

    end subroutine like_init


    subroutine compute_fg(C_foregrounds,freq_params, lmax_compute)
    real(8), intent(in)  :: freq_params(:)
    integer, intent(in) :: lmax_compute
    real(8), intent(inout)  ::  C_foregrounds(:,:)
    real(campc) A_ps_100, A_ps_143, A_ps_217, A_cib_143, A_cib_217, asz143, r_ps, r_cib, xi, A_ksz
    real(campc) ncib217, ncib143
    real(campc) zCIB
    integer:: l, i
    real(campc) cl_cib_143(CAMSpec_lmax_foreground), cl_cib_217(CAMSpec_lmax_foreground) !CIB
    real(campc), parameter :: sz_bandpass100_nom143 = 2.022d0
    real(campc), parameter :: cib_bandpass143_nom143 = 1.134d0
    real(campc), parameter :: sz_bandpass143_nom143 = 0.95d0
    real(campc), parameter :: cib_bandpass217_nom217 = 1.33d0
    real(campc), parameter :: ps_scale  = 1.d-6/9.d0
    real(campc) :: A_cib_217_bandpass, A_sz_143_bandpass, A_cib_143_bandpass
    integer lmin(Nspec), lmax(Nspec)
    real(campc) lnrat, nrun_cib
    real(campc) wigamp(2), wiggle_corr, wiggle_center, wiggle_width, wiggle_cl(CAMspec_lmax)
    integer f_ix

    if (lmax_compute/=0) then
        lmin = 2
        lmax = min(CAMSpec_lmax_foreground,lmax_compute)
    else
        lmin = lminX
        lmax = lmaxX
    end if

    A_ps_100=freq_params(1)
    A_ps_143 = freq_params(2)
    A_ps_217 = freq_params(3)
    A_cib_143 =freq_params(4)
    A_cib_217 =freq_params(5)
    asz143 = freq_params(6)
    r_ps = freq_params(7)
    r_cib = freq_params(8)
    ncib143 = freq_params(9)
    ncib217 = freq_params(10)
    nrun_cib = freq_params(11)
    xi = freq_params(12)
    A_ksz = freq_params(13)
    f_ix = 14

    do l=1, maxval(lmax)
        lnrat = log(real(l,campc)/CamSpec_cib_pivot)
        cl_cib_217(l) = exp(ncib217*lnrat + nrun_cib/2*lnrat**2)
        if (ncib143<-9) then
            cl_cib_143(l) = cl_cib_217(l)
        else
            cl_cib_143(l)=exp(ncib143*lnrat + nrun_cib/2*lnrat**2)
        end if
    end do

    !   100 foreground
    !
    do l = lmin(1), lmax(1)
        C_foregrounds(l,1)= A_ps_100*ps_scale+  &
            ( A_ksz*ksz_temp(l) + asz143*sz_bandpass100_nom143*sz_143_temp(l) )*llp1(l)
    end do

    !   143 foreground
    !
    A_sz_143_bandpass = asz143 * sz_bandpass143_nom143
    A_cib_143_bandpass = A_cib_143 * cib_bandpass143_nom143
    do l = lmin(2), lmax(2)
        zCIB = A_cib_143_bandpass*cl_cib_143(l)
        C_foregrounds(l,2)= A_ps_143*ps_scale + &
            (zCIB +  A_ksz*ksz_temp(l) + A_sz_143_bandpass*sz_143_temp(l) &
            -2.0*sqrt(A_cib_143_bandpass * A_sz_143_bandpass)*xi*tszxcib_temp(l) )*llp1(l)
    end do

    !   217 foreground
    !
    A_cib_217_bandpass = A_cib_217 * cib_bandpass217_nom217
    do l = lmin(3), lmax(3)
        zCIB = A_cib_217_bandpass*cl_cib_217(l)
        C_foregrounds(l,3) = A_ps_217*ps_scale + (zCIB + A_ksz*ksz_temp(l) )*llp1(l)
    end do

    !   143x217 foreground
    !
    do l = lmin(4), lmax(4)
        zCIB = sqrt(A_cib_143_bandpass*A_cib_217_bandpass*cl_cib_143(l)*cl_cib_217(l))
        C_foregrounds(l,4) = r_ps*sqrt(A_ps_143*A_ps_217)*ps_scale + &
            ( r_cib*zCIB + A_ksz*ksz_temp(l) -sqrt(A_cib_217_bandpass * A_sz_143_bandpass)*xi*tszxcib_temp(l) )*llp1(l)
    end do

    do i=1,2
        wigamp = freq_params(f_ix:f_ix+1)
        if (any(wigamp/=0)) then
            wiggle_corr = freq_params(f_ix+2)
            wiggle_center = freq_params(f_ix+3)
            wiggle_width = freq_params(f_ix+4)
            wiggle_cl=0
            do L= nint(wiggle_center) - nint(3*wiggle_width),nint(wiggle_center) + nint(3*wiggle_width)
                wiggle_cl(L) = exp(-(L-wiggle_center)**2/real(2*wiggle_width**2,campc))/(L*(L+1))
            end do

            C_foregrounds(lmin(2):lmax(2),2)= C_foregrounds(lmin(2):lmax(2),2)+ wigamp(1)*wiggle_cl(lmin(2):lmax(2))
            C_foregrounds(lmin(3):lmax(3),3)= C_foregrounds(lmin(3):lmax(3),3)+ wigamp(2)*wiggle_cl(lmin(3):lmax(3))
            C_foregrounds(lmin(4):lmax(4),4)= C_foregrounds(lmin(4):lmax(4),4)+ &
                & wiggle_corr*sqrt(wigamp(2)*wigamp(1))*wiggle_cl(lmin(4):lmax(4))
        end if
        f_ix = f_ix + 5
    end do

    end subroutine compute_fg

    subroutine calc_like(zlike,  cell_cmb,cell_cmbx,cell_cmbe, freq_params)
    real(campc), intent(in)  :: freq_params(:)
    real(campc), dimension(0:) :: cell_cmb, cell_cmbx,cell_cmbe
    integer ::  j, l, ii,jj
    real(campc) , allocatable, save ::  X_beam_corr_model(:), Y(:),  C_foregrounds(:,:)
    real(campc) cal0, cal1, cal2
    real(campc) zlike, ztemp
    real(campc) beam_params(cov_dim),beam_coeffs(num_modes_per_beam,beam_Nspec)
    integer :: ie1,ie2,if1,if2, ix
    integer num_non_beam

    if (.not. allocated(lminX)) then
        print*, 'like_init should have been called before attempting to call calc_like.'
        stop
    end if
    if(.not.(Nspec.eq.4.or.Nspec.eq.6)) then
        print*, 'Nspec inconsistent with foreground corrections in calc_like.'
        stop
    end if
    if (.not. allocated(Y)) then
        allocate(X_beam_corr_model(1:nX))
        allocate(Y(1:nX))
        allocate(C_foregrounds(CAMspec_lmax,Nspec))
        C_foregrounds=0
    end if
    if (camspec_has_TT) then
        call compute_fg(C_foregrounds,freq_params, 0)

        cal0 = freq_params(24)
        cal1 = freq_params(25)
        cal2 = freq_params(26)

        num_non_beam = 26
        if (size(freq_params) < num_non_beam +  beam_Nspec*num_modes_per_beam) stop 'CAMspec: not enough parameters'

        if (keep_num>0) then
            beam_params = freq_params(num_non_beam+1:num_non_beam+cov_dim)
            !set marged beam parameters to their mean subject to fixed non-marged modes
            if (marge_num>0) beam_params(marge_indices) = matmul(beam_conditional_mean, beam_params(keep_indices))

            do ii=1,beam_Nspec
                do jj=1,num_modes_per_beam
                    ix = jj+num_modes_per_beam*(ii-1)
                    beam_coeffs(jj,ii)=beam_params(ix)
                enddo
            enddo
        else
            beam_coeffs=0
        end if

        do l = lminX(1), lmaxX(1)
            X_beam_corr_model(l-lminX(1)+1) = ( cell_cmb(l) + C_foregrounds(l,1) )* corrected_beam(1,l)/cal0
        end do

        do l = lminX(2), lmaxX(2)
            X_beam_corr_model(l-lminX(2)+npt(2)) =  (cell_cmb(l)+ C_foregrounds(l,2))*corrected_beam(2,l)/cal1
        end do

        do l = lminX(3), lmaxX(3)
            X_beam_corr_model(l-lminX(3)+npt(3)) = (cell_cmb(l)+ C_foregrounds(l,3))* corrected_beam(3,l)/cal2
        end do

        do l = lminX(4), lmaxX(4)
            X_beam_corr_model(l-lminX(4)+npt(4)) =  ( cell_cmb(l) + C_foregrounds(l,4))*corrected_beam(4,l)/sqrt(cal1*cal2)
        end do
    end if

    if(Nspec.eq.6) then
        ! TE...
        do l = lminX(5), lmaxX(5)
            X_beam_corr_model(l-lminX(5)+npt(5)) = cell_cmbx(l)
        end do

        ! EE...
        do l = lminX(6), lmaxX(6)
            X_beam_corr_model(l-lminX(6)+npt(6)) = cell_cmbe(l)
        end do
    endif

    Y = X_data - X_beam_corr_model

    zlike = 0
    !$OMP parallel do private(j,ztemp) reduction(+:zlike) schedule(static,16)
    do  j = 1, nX
        ztemp= dot_product(Y(j+1:nX), c_inv(j+1:nX, j))
        zlike=zlike+ (ztemp*2 +c_inv(j, j)*Y(j))*Y(j)
    end do

    if (keep_num>0 .and. camspec_has_TT) then
        do if2=1,beam_Nspec
            do if1=1,beam_Nspec
                do ie2=1,num_modes_per_beam
                    jj=ie2+num_modes_per_beam*(if2-1)
                    if (.not. want_marge(jj)) then
                        do ie1=1,num_modes_per_beam
                            ii=ie1+num_modes_per_beam*(if1-1)
                            if (.not. want_marge(ii)) then
                                zlike=zlike+beam_coeffs(ie1,if1)*&
                                    beam_cov_inv(keep_indices_reverse(ii),keep_indices_reverse(jj))*beam_coeffs(ie2,if2)
                            end if
                        enddo
                    end if
                enddo
            enddo
        enddo
    end if

    !Moved these hard coded priors in prior[cal0] and prior[cal2] input parameters
    !    if (want_spec(1)) zlike=zlike+ ((cal0/cal1-1.0006d0)/0.0004d0)**2
    !    if (any(want_spec(3:4))) zlike=zlike+ ((cal2/cal1-0.9966d0)/0.0015d0)**2

    contains

    real(campc) function corrected_beam(spec_num,l)
    integer, intent(in) :: spec_num,l

    corrected_beam=1.d0 + dot_product(beam_coeffs(:,spec_num),beam_modes(:,l,spec_num))*beam_factor

    end function corrected_beam

    end subroutine calc_like


    end module temp_like_camspec
