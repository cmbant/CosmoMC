    module temp_like_camspec4
    !Now supports cleaned likelihood simple power law foreground model
    !and arbitary L cuts (at all frequencies)
    implicit none
    private
    integer, parameter :: campc = KIND(1.d0)

    real(campc), dimension(:), allocatable :: X_data
    real(campc),  dimension(:,:), allocatable :: c_inv
    integer :: Nspec,nX,nX_cut, num_ells,nXfromdiff, CAMspec_lmax
    integer, dimension(:), allocatable  :: lminX, lmaxX, npt
    integer, parameter :: CAMSpec_lmax_foreground = 5000
    integer, parameter :: lmax_sz = CAMSpec_lmax_foreground
    real(campc) :: sz_143_temp(lmax_sz)
    real(campc) :: ksz_temp(lmax_sz), tszxcib_temp(lmax_sz)

    real(campc) :: cib_217_temp(lmax_sz)
    real(campc) :: dust_temp(lmax_sz,4)
    logical :: applydust=.false.
    logical :: cibfromfile=.false.

    logical :: dump_spectra = .false.

    real(campc), dimension(:,:), allocatable :: beam_cov,beam_cov_inv
    real(campc), dimension(:,:,:), allocatable :: beam_modes ! mode#, l, spec#
    integer :: num_modes_per_beam,beam_lmax,beam_Nspec,cov_dim

    integer, allocatable :: marge_indices(:),marge_indices_reverse(:)
    integer, allocatable :: keep_indices(:),keep_indices_reverse(:)
    real(campc), allocatable :: beam_conditional_mean(:,:)
    logical, allocatable :: want_marge(:)
    integer marge_num, keep_num

    logical :: make_cov_marged = .false.
    real(campc) :: beam_factor = 1.0_campc
    integer, parameter :: CamSpec_cib_pivot = 3000
    integer, parameter :: CamSpec_sz_pivot = 3000
    integer, parameter :: CamSpec_powerlaw_pivot = 1500

    logical :: want_spec(6) = .true.
    logical :: camspec_has_TT = .true.
    integer :: camspec_lmins(6) =0
    integer :: camspec_lmaxs(6) =0
    integer, allocatable :: L_cut_indices(:)
    character(LEN=128) :: camspec_ranges = ''
    real(campc)  :: llp1(CAMSpec_lmax_foreground)

    !paramters for systematics testing
    integer :: camspec_nwiggles = 0

    integer :: camspec_beam_mcmc_num = 1

    character(LEN=1024) camspec_fiducial_foregrounds, camspec_fiducial_cl
    character(LEN=80) :: marge_file_variant = ''

    character(LEN=*), parameter :: CAMSpec_like_version = 'CamSpec_v5_cleaned'

    logical :: apply_tight_sz_prior=.false.

    public like_init,calc_like,CAMSpec_like_version, camspec_beam_mcmc_num, dump_spectra, camspec_ranges, &
        want_spec,camspec_lmins,camspec_lmaxs, make_cov_marged, marge_file_variant,&
        compute_fg, Nspec,CAMSpec_lmax_foreground,camspec_fiducial_foregrounds, &
        camspec_fiducial_cl, apply_tight_sz_prior


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

    call CAMspec_ReadNorm(fname, templt, lmax_sz)
    end subroutine CAMspec_ReadNormSZ

    subroutine CAMspec_ReadNorm(fname, templt,thelmax)
    character(LEN=*), intent(in) :: fname
    real(campc) :: templt(lmax_sz)
    integer i, dummy, thelmax
    real(campc) :: renorm
    integer unit

    print *,'reading ', fname, ' in readnorm'

    if(thelmax < CamSpec_sz_pivot) stop 'CAMspec_ReadNorm: lmax too low'

    open(newunit=unit, file=fname, form='formatted', status='old')
    do i=2,thelmax
        read(unit,*) dummy,templt(i)
        if (dummy/=i) stop 'CAMspec_ReadNorm: inconsistency in file read'
    enddo
    close(unit)

    renorm=1.d0/templt(CamSpec_sz_pivot)
    templt=templt*renorm
    end subroutine CAMspec_ReadNorm

    subroutine CAMspec_Read(fname, templt,thelength)
    character(LEN=*), intent(in) :: fname
    real(campc) :: templt(lmax_sz)
    integer i, dummy, readlmax
    integer, optional :: thelength
    integer unit

    print *,'reading ', fname, ' in read'

    if(present(thelength)) then
        readlmax=min(thelength,lmax_sz)
    else
        readlmax=lmax_sz
    endif

    open(newunit=unit, file=fname, form='formatted', status='old')
    do i=2,readlmax
        read(unit,*) dummy,templt(i)
        if (dummy/=i) stop 'CAMspec_Read: inconsistency in file read'
    enddo
    close(unit)
    !    do i=100,readlmax,100
    !        print *, i, templt(i)
    !    enddo
    end subroutine CAMspec_Read

    subroutine ReadFiducialCl(fid_cl)
    integer i,j
    real(campc) :: fid_theory_tt,fid_theory_te,fid_theory_ee
    real(campc) :: fid_cl(:,:)
    integer unit

    open(newunit=unit, file=CAMspec_fiducial_foregrounds, form='formatted', status='old')
    do i=2,CAMspec_lmax
        read(unit,*) j, fid_cl(i,1:4)
        if(Nspec.eq.6) then
            fid_cl(i,5)=0.0
            fid_cl(i,6)=0.0
        endif
        if (j/=i) stop 'error reading fiducial foreground C_l for beams'
    enddo
    close(unit)
    open(newunit=unit, file=CAMspec_fiducial_cl, form='formatted', status='old')
    if (index(CAMspec_fiducial_cl,'minimum.theory_cl')>0) read(unit,*)
    do i=2,CAMspec_lmax
        read(unit,*, end=100) j,fid_theory_tt,fid_theory_te,fid_theory_ee
100     if (j/=i) stop 'error reading fiducial C_l for beams'
        if(Nspec.eq.4) then
            fid_cl(i,:) = (fid_cl(i,:) + fid_theory_tt)*llp1(i) !want C_l/2Pi foregrounds+theory
        else if(Nspec.eq.6) then
            fid_cl(i,1:4) = (fid_cl(i,1:4) + fid_theory_tt)*llp1(i) !want C_l/2Pi foregrounds+theory
            fid_cl(i,5) = (fid_cl(i,5) + fid_theory_te)*llp1(i) !want C_l/2Pi foregrounds+theory
            fid_cl(i,6) = (fid_cl(i,6) + fid_theory_ee)*llp1(i) !want C_l/2Pi foregrounds+theory
        endif
    enddo
    close(unit)

    end subroutine ReadFiducialCl

    subroutine read_range(S, output)
    use ObjectLists
    use StringUtils
    character(LEN=*), intent(in) :: S
    integer, intent(out), allocatable :: output(:)
    Type(TStringList) :: L
    integer i,j,k,ix, out_ix
    character(LEN=:), allocatable :: item
    integer outitems(5000)

    allocate(item, source=S)
    j=0
    out_ix=0
    do i=1, len_trim(S)
        if (verify(S(i:i),'0123456789') == 0) then
            j=j+1
            item(j:j) = S(i:i)
        else
            if (j>0) call L%Add(item(1:j))
            if (S(i:i) == '-') then
                call L%Add('-')
            end if
            j=0
        end if
    end do
    if (j>0) call L%Add(item(1:j))
    ix =1
    do while (ix <= L%Count)
        i =StrToInt(L%Item(ix))
        if (ix<L%Count .and. L%Item(ix+1)=='-') then
            j = StrToInt(L%Item(ix+2))
            ix = ix+2
        else
            j=i
        end if
        do j=i,j
            out_ix=out_ix+1
            if (out_ix>5000) stop 'too many L'
            outitems(out_ix) = j
        end do
        ix =ix+1
    end do
    allocate(output, source = outitems(1:out_ix))
    end subroutine read_range

    subroutine like_init(pre_marged,like_file, sz143_file, tszxcib_file, ksz_file, beam_file,data_vector,&
        cib217_file,dust100_file,dust143_file,dust217_file,dust143x217_file)
    use MatrixUtils
    use FileUtils
    use ArrayUtils
    integer :: i, j,k,l
    logical :: pre_marged
    character(LEN=*), intent(in) :: like_file, sz143_file, ksz_file, tszxcib_file, beam_file,data_vector
    character(LEN=*), intent(in) ::  cib217_file,dust100_file,dust143_file,dust217_file,dust143x217_file
    logical, save :: needinit=.true.
    real(campc) , allocatable:: fid_cl(:,:), beam_cov(:,:), beam_cov_full(:,:)
    integer if1,if2, ie1, ie2, ii, jj, L2
    real(campc), dimension(:), allocatable :: X_data_in
    real(campc),  dimension(:,:), allocatable :: cov
    integer, allocatable :: indices(:),np(:)
    integer ix, status
    real(campc) dummy
    real(campc), allocatable :: CL_in(:,:), in_data(:)
    real(campc), allocatable :: CL_out(:,:), CL_scale(:)
    real, allocatable :: real_cov(:,:)
    character(LEN=7), parameter :: spectrum_names(6) = &
        [character(7)::'100x100','143x143','217x217','143x217','TE','EE']
    Type(TBinaryFile) :: binary_file
    integer, allocatable :: L_keep(:)
    integer unit


    if(cib217_file /='') print*,'cib217 file present in like_init:',cib217_file
    if(dust100_file/='') print*,'dust100 file present in like_init:',dust100_file
    if(dust143_file/='') print*,'dust143 file present in like_init:',dust143_file
    if(dust217_file/='') print*,'dust217 file present in like_init:',dust217_file
    if(dust143x217_file/='') print*,'dust143x217 file present in like_init:',dust143x217_file

    if(.not. needinit) return

    do i=1,CAMSpec_lmax_foreground
        llp1(i) = 1/real(i*(i+1),campc)
    end do

    camspec_has_TT = any(want_spec(1:4))

    open(newunit=unit, file=like_file, form='unformatted', status='old')

    read(unit) Nspec,nX
    allocate(lminX(Nspec))
    allocate(lmaxX(Nspec))
    allocate(np(Nspec))
    allocate(npt(Nspec))

    read(unit) (lminX(i), lmaxX(i), np(i), npt(i), i = 1, Nspec)
    if (pre_marged) then
        print *,'Using pre-beam-marged covariance; L ranges etc ignored'
        allocate(X_data(nX))
        allocate(c_inv(nX,nX))
        nX_cut = nX
        read(unit) (X_data(i), i=1, nX)
        read(unit) !skip  covariance
        read(unit) c_inv !((c_inv(i, j), j = 1, nX), i = 1,  nX) !inv covariance
        close(unit)
        do i=1,Nspec
            print *, spectrum_names(i), ' lmin,lmax, start ix = ', lminX(i),lmaxX(i), npt(i)
        end do
    else
        !Chop L ranges/frequencies if requested
        allocate(indices(nX))
        allocate(cov(nX,nX))
        allocate(X_data_in(nX))
        read(unit) X_data_in !(X_data(i), i=1, nX)
        read(unit) cov !((c_inv(i, j), j = 1, nX), i = 1,  nX) !covariance
        read(unit) !((c_inv(i, j), j = 1, nX), i = 1,  nX) !inver covariance
        close(unit)

        if (dump_spectra) then
            allocate(CL_out(0:maxval(lmaxX(1:Nspec)),Nspec))
            allocate(CL_scale(nX))
            CL_out=0
            ix=0
            do i=1,NSpec
                do L=lminX(i), lmaxX(i)
                    ix=ix+1
                    CL_scale(ix) = L*(L+1)
                    CL_out(L,i)=L*(L+1)*X_data_in(ix)
                end do
            end do
            call File%saveTxt(File%ExtractName(like_file) // '_spectra.txt', CL_out)
            deallocate(CL_out)
            allocate(real_cov(nX,nX))
            do ix=1, nX
                do i = 1, nX
                    real_cov(ix,i) = cov(ix,i)*CL_scale(ix)*CL_scale(i)
                end do
            end do
            call binary_file%CreateFile(File%ExtractName(like_file) // '_cov.bin')
            call binary_file%Write(real_cov)
            call binary_file%Close()
            open(newunit=unit,file=File%ExtractName(like_file) // '_data_ranges.txt', form='formatted', status='replace')
            do i=1, nSpec
                write(unit,*) spectrum_names(i), lminX(i), lmaxX(i)
            end do
            close(unit)
            call MpiStop()
        end if

        !Cut on L ranges
        print *,'Determining L ranges'
        if (camspec_ranges/='') then
            call read_range(camspec_ranges, L_keep)
            print *, 'Only keeping L=', trim(camspec_ranges)
        end if
        j=0
        ix=0
        npt(1)=1
        if(.not. want_spec(1) .and. marge_num>0) stop 'One beam mode may not be right here'
        do i=1,Nspec
            do l = lminX(i), lmaxX(i)
                j =j+1
                if (want_spec(i) .and. (camspec_lmins(i)==0 .or. l>=camspec_lmins(i)) &
                    .and. (camspec_lmaxs(i)==0 .or. l<=camspec_lmaxs(i)) .and. &
                    (.not. allocated(L_keep) .or. indexOf(l,L_keep,size(L_keep))/=0)) then
                ix =ix+1
                indices(ix) = j
                end if
            end do
            if (.not. allocated(L_keep)) then
                if (want_spec(i)) then
                    if (camspec_lmins(i)/=0) lminX(i) = max(camspec_lmins(i), lminX(i))
                    if (camspec_lmaxs(i)/=0) lmaxX(i) = min(camspec_lmaxs(i), lmaxX(i))
                else
                    lmaxX(i) =0
                end if
                if (i<NSpec) npt(i+1) = ix+1
                print *, spectrum_names(i), ' lmin,lmax, start ix = ', lminX(i),lmaxX(i), npt(i)
            elseif (want_spec(i)) then
                print *, spectrum_names(i), 'pre-cut lmin,lmax, start ix = ', lminX(i),lmaxX(i), npt(i)
            end if
        end do
        if (j/=nx) stop 'Error cutting camspec cov matrix'
        if (allocated(L_keep)) then
            allocate(L_cut_indices, source = indices(1:ix))
        else
            nX = ix
        end if
        nX_cut = ix
        allocate(X_data(nX_cut))
        allocate(c_inv(nX_cut,nX_cut))
        do i=1, nX_cut
            c_inv(:,i) = cov(indices(1:nX_cut),indices(i)) !c_inv is covariance here
            X_data(i) = X_data_in(indices(i))
        end do
        !    call Matrix_inverse(c_inv)
        deallocate(cov)
    end if
    CAMspec_lmax = maxval(lmaxX(1:Nspec))

    if (data_vector/='') then
        !!override data vector in the main camspec file
        !The input file is in comprehensible L, L(L+1)C_L^i/2pi format
        open(newunit=unit, file=data_vector, form='formatted', status='old', iostat=status)
        if (status/=0) stop 'Error opening camspec data_vector override'
        print *,'Using CamSpec data vector : '//data_vector
        allocate(CL_in(CAMspec_lmax, Nspec))
        allocate(in_data(Nspec))
        L=0
        do
            read(unit,*,iostat = status) L, in_data
            if (status/=0 .or. L > CAMspec_lmax) exit
            CL_in(L,:) = in_data
        end do
        if (L< CAMspec_lmax) stop 'Error reading camspec data_vector override'
        close(unit)
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


    dust_temp=0
    applydust=.false.

    if(dust100_file/='') then
        call CAMspec_Read(dust100_file, dust_temp(:,1),3000)
        applydust=.true.
    endif

    if(dust143_file/='') then
        call CAMspec_Read(dust143_file, dust_temp(:,2),3000)
        applydust=.true.
    endif

    if(dust217_file/='') then
        call CAMspec_Read(dust217_file, dust_temp(:,3),3000)
        applydust=.true.
    endif

    if(dust143x217_file/='') then
        call CAMspec_Read(dust143x217_file, dust_temp(:,4),3000)
        applydust=.true.
    endif

    if(cib217_file/='') then
        call CAMspec_ReadNormSZ(cib217_file, cib_217_temp)
        cibfromfile=.true.
    endif


    open(newunit=unit, file=beam_file, form='unformatted', status='unknown')
    read(unit) beam_Nspec,num_modes_per_beam,beam_lmax
    ! removed until we decide about how to handle beam errors for TE & EE
    !    if(beam_Nspec.ne.Nspec) stop 'Problem: beam_Nspec != Nspec'
    allocate(beam_modes(num_modes_per_beam,0:beam_lmax,beam_Nspec))
    cov_dim=beam_Nspec*num_modes_per_beam
    allocate(beam_cov_full(cov_dim,cov_dim))
    read(unit) (((beam_modes(i,l,j),j=1,beam_Nspec),l=0,beam_lmax),i=1,num_modes_per_beam)
    allocate(beam_cov_inv(cov_dim,cov_dim))
    read(unit) ((beam_cov_full(i,j),j=1,cov_dim),i=1,cov_dim)  ! beam_cov
    read(unit) ((beam_cov_inv(i,j),j=1,cov_dim),i=1,cov_dim) ! beam_cov_inv
    close(unit)

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
                                                beam_cov(marge_indices_reverse(ii),marge_indices_reverse(jj)) * &
                                                beam_modes(ie2,L2,if2) *fid_cl(L2,if2)*fid_cl(lminX(if1):lmaxX(if1),if1)
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
                open(newunit=unit, file=trim(like_file)//'_beam_marged'//trim(marge_file_variant), &
                    form='unformatted', status='unknown')
                write(unit) Nspec,nX
                write(unit) (lminX(i), lmaxX(i), np(i), npt(i), i = 1, Nspec)
                write(unit) (X_data(i), i=1, nX)
                dummy=-1
                write(unit) dummy !inver covariance, assume not used
                write(unit) ((c_inv(i, j), j = 1, nX), i = 1,  nX) !inver covariuance
                close(unit)
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
    real(campc) A_cib_143, A_cib_217, asz143, r_ps, r_cib, xi, A_ksz
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
    real(campc) dust(4), pointsource(4)
    real(campc) power_law_amp(4), power_law_tilt(4)
    integer f_ix

    if (lmax_compute/=0) then
        lmin = 2
        lmax = min(CAMSpec_lmax_foreground,lmax_compute)
    else
        lmin = lminX
        lmax = lmaxX
    end if

    pointsource(1:3) = freq_params(1:3)
    A_cib_143 =freq_params(4)
    A_cib_217 =freq_params(5)
    asz143 = freq_params(6)
    r_ps = freq_params(7)
    pointsource(4) = r_ps*sqrt(pointsource(2)*pointsource(3)) !143x217 correlation
    r_cib = freq_params(8)
    ncib143 = freq_params(9)
    ncib217 = freq_params(10)
    nrun_cib = freq_params(11)
    xi = freq_params(12)
    A_ksz = freq_params(13)
    dust = freq_params(14:17)
    power_law_amp = freq_params(18:21)
    power_law_tilt = freq_params(22:25)

    f_ix = 26

    do l=1, maxval(lmax)
        lnrat = log(real(l,campc)/CamSpec_cib_pivot)
        cl_cib_217(l) = exp(ncib217*lnrat + nrun_cib/2*lnrat**2)
        if (ncib143<-9) then
            cl_cib_143(l) = cl_cib_217(l)
        else
            cl_cib_143(l)=exp(ncib143*lnrat + nrun_cib/2*lnrat**2)
        end if
    end do

    if(cibfromfile) then
        !In this case take CIB shape from file; unchanged if ncib=nrun_cib=0, otherwise tilted accordingly
        !(ncib etc allowed so we can do consistency checks on the model)
        do l=1,maxval(lmax)
            cl_cib_143(l)=cib_217_temp(l)*cl_cib_143(l)
            cl_cib_217(l)=cib_217_temp(l)*cl_cib_217(l)
        end do
        if (A_cib_143<-0.99) then
            A_cib_143=.094*A_cib_217/cib_bandpass143_nom143*cib_bandpass217_nom217
            !The above came from ratioing Paolo's templates, which were already colour-corrected,
            !and assumed perfect correlation
        end if
        !r_cib=1.0 !assume set in .ini, allow to vary to be able to check consistent with 1
    elseif (A_cib_143 < 0) then
        stop 'A_cib_143 < 0 but no cib from file'
    endif
    !   100 foreground
    !
    do l = lmin(1), lmax(1)
        C_foregrounds(l,1)= A_ksz*ksz_temp(l) + asz143*sz_bandpass100_nom143*sz_143_temp(l)
    end do

    !   143 foreground
    !
    A_sz_143_bandpass = asz143 * sz_bandpass143_nom143
    A_cib_143_bandpass = A_cib_143 * cib_bandpass143_nom143
    do l = lmin(2), lmax(2)
        zCIB = A_cib_143_bandpass*cl_cib_143(l)
        C_foregrounds(l,2)= zCIB +  A_ksz*ksz_temp(l) + A_sz_143_bandpass*sz_143_temp(l) &
            -2.0*sqrt(A_cib_143_bandpass * A_sz_143_bandpass)*xi*tszxcib_temp(l)
    end do

    !   217 foreground
    !
    A_cib_217_bandpass = A_cib_217 * cib_bandpass217_nom217
    do l = lmin(3), lmax(3)
        zCIB = A_cib_217_bandpass*cl_cib_217(l)
        C_foregrounds(l,3) = zCIB + A_ksz*ksz_temp(l)
    end do

    !   143x217 foreground
    !
    do l = lmin(4), lmax(4)
        zCIB = sqrt(A_cib_143_bandpass*A_cib_217_bandpass*cl_cib_143(l)*cl_cib_217(l))
        C_foregrounds(l,4) =  r_cib*zCIB + A_ksz*ksz_temp(l) &
            - sqrt(A_cib_217_bandpass * A_sz_143_bandpass)*xi*tszxcib_temp(l)
    end do


    if(applydust) then
        do i=1,4
            do l=lmin(i),lmax(i)
                C_foregrounds(l,i) = C_foregrounds(l,i)+dust(i)*dust_temp(l,i)
            enddo
        end do
    endif

    do i=1,4
        !Additional power laws, e.g. for use with foreground-cleaned likelihoods
        if (power_law_amp(i) /= 0.0_campc) then
            do l=lmin(i),lmax(i)
                lnrat = log(real(l,campc)/CamSpec_powerlaw_pivot)
                C_foregrounds(l,i)= C_foregrounds(l,i) + power_law_amp(i) * exp(power_law_tilt(i) *lnrat)
            enddo
        end if
        !Divide by L(L+1) factors and add white point sources
        do l=lmin(i),lmax(i)
            C_foregrounds(l,i)=C_foregrounds(l,i)*llp1(l) + ps_scale*pointsource(i)
        enddo
    end do

    do i=1,2
        ! Not used, just testing putting in wiggles.
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
    real(campc) calTE, calEE, calPlanck
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
        allocate(Y(1:nX_cut))
        allocate(C_foregrounds(CAMspec_lmax,Nspec))
        C_foregrounds=0
    end if

    if (camspec_has_TT) call compute_fg(C_foregrounds,freq_params, 0)

    calPlanck = freq_params(36)**2 !Total calibration that scales everything
    cal0 = freq_params(37)*calPlanck
    cal1 = freq_params(38)*calPlanck
    cal2 = freq_params(39)*calPlanck
    calTE = freq_params(40)*calPlanck
    calEE = freq_params(41)*calPlanck
    num_non_beam = 41

    if (camspec_has_TT) then
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
            X_beam_corr_model(l-lminX(5)+npt(5)) = cell_cmbx(l)/calTE
        end do

        ! EE...
        do l = lminX(6), lmaxX(6)
            X_beam_corr_model(l-lminX(6)+npt(6)) = cell_cmbe(l)/calEE
        end do
    endif

    if (allocated(L_cut_indices)) then
        Y = X_data  - X_beam_corr_model(L_cut_indices)
    else
        Y = X_data - X_beam_corr_model
    end if

    zlike = 0
    !$OMP parallel do private(j,ztemp) reduction(+:zlike) schedule(static,16)
    do  j = 1, nX_cut
        ztemp= dot_product(Y(j+1:nX_cut), c_inv(j+1:nX_cut, j))
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

    if(apply_tight_sz_prior) then
        ! not used, prior is put in when running not as part of likelihood
        !       asz143 = freq_params(6)
        !       A_ksz = freq_params(13)
        zlike=zlike+((freq_params(13)+1.6*freq_params(6)-9.5)/3.0)**2
    endif

    contains

    real(campc) function corrected_beam(spec_num,l)
    integer, intent(in) :: spec_num,l

    corrected_beam=1.d0 + dot_product(beam_coeffs(:,spec_num),beam_modes(:,l,spec_num))*beam_factor

    end function corrected_beam

    end subroutine calc_like


    !Just for checking eigenvalues etc. slower than cholesky
    subroutine Matrix_Inverse2(M)
    use MpiUtils
    use MatrixUtils
    !Inverse of symmetric positive definite matrix
    real(campc), intent(inout):: M(:,:)
    integer i, n
    real(campc) w(Size(M,DIM=1))
    real(campc), dimension(:,:), allocatable :: tmp, tmp2
    real(campc), dimension(:), allocatable :: norm

    n = size(M,DIM=1)
    do i=1, size(M,DIM=1)
        if (abs(M(i,i)) < 1d-30) call MpiStop('Matrix_Inverse: very small diagonal'  )
    end do
    allocate(tmp(Size(M,DIM=1),Size(M,DIM=1)))

    n=Size(M,DIM=1)
    if (n<=1) return
    if (Size(M,DIM=2)/=n) call MpiStop('Matrix_Inverse: non-square matrix')

    allocate(norm(n))
    do i=1, n
        norm(i) = sqrt(abs(M(i,i)))
        if (norm(i) < 1d-30) &
            call MpiStop('Matrix_Inverse: very small diagonal'  )
        M(i,:) = M(i,:)/norm(i)
        M(:,i) = M(:,i)/norm(i)
    end do

    call Matrix_Write('smallmat',M,.true.)

    call Matrix_Diagonalize(M,w,n)
    call Matrix_Write('eigenvecs',M,.true.)
    !print *,'evalues:', w

    write (*,*) 'min/max eigenvalues = ', minval(w), maxval(w)
    if (any(w<=0)) then
        write (*,*) 'min/max eigenvalues = ', minval(w), maxval(w)
        call MpiStop('Matrix_Inverse: negative or zero eigenvalues')
    end if
    do i=1, n
        tmp(i,:) = M(:,i)/w(i)
    end do
    allocate(tmp2(Size(M,DIM=1),Size(M,DIM=1)))
    call Matrix_Mult(M,tmp,tmp2)
    M = tmp2
    do i=1, n
        M(i,:) = M(i,:)/norm(i)
        M(:,i) = M(:,i)/norm(i)
    end do
    deallocate(tmp, tmp2)
    deallocate(norm)
    call Matrix_End('Inverse')

    end subroutine Matrix_Inverse2
    end module temp_like_camspec4
