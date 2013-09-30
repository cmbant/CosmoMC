    !Module storing observed matter power spectrum datasets, their points and window functions
    !and routines for computing the likelihood

    !This code is based on that in cmbdata.f90
    !and on Sam Leach's incorporation of Max Tegmark's SDSS code
    !
    !Originally SLB Sept 2004
    !AL April 2006: added covariance matrix support (following 2df 2005)
    !LV_06 : incorporation of LRG DR4 from Tegmark et al . astroph/0608632
    !AL: modified LV SDSS to do Q and b^2 or b^2*Q marge internally as for 2df
    !BR09: added model LRG power spectrum.
    !AL Oct 20: switch to Ini_Read_xxx_File; fortran compatibility changes

    !JD 29/07/2013 removed LRG stuff
    !JD 03/08/2013 fixed compute_scaling_factor and associated functions
    !to work with w_a/=0

    !JD 09/13: Replaced compute_scaling_factor routines with routines that use CAMB's
    !          built in D_V function.


    module mpk
    use precision
    use settings
    use cmbtypes
    use likelihood
    use wigglezinfo
    use wigglez, only : WiggleZLikelihood_Add
    implicit none

    type, extends(CosmologyLikelihood) :: MPKLikelihood
        logical :: use_set
        integer :: num_mpk_points_use ! total number of points used (ie. max-min+1)
        integer :: num_mpk_kbands_use ! total number of kbands used (ie. max-min+1)
        real(mcp), pointer, dimension(:,:) :: N_inv
        real(mcp), pointer, dimension(:,:) :: mpk_W, mpk_invcov
        real(mcp), pointer, dimension(:) :: mpk_P, mpk_sdev, mpk_k
        logical :: use_scaling !as SDSS_lrgDR3
        !JD 09/13  New variables so we can use camb routines to calculate a_scl
        real(mcp) DV_fid   !Fiducial D_V
        real(mcp) redshift !effective redshift of the P(k) for computing a_scl

        !for Q and A see e.g. astro-ph/0501174, astro-ph/0604335
        logical :: Q_marge, Q_flat
        real(mcp) :: Q_mid, Q_sigma, Ag
    contains
    procedure :: LogLike => MPK_Lnlike
    end type MPKLikelihood

    logical :: use_mpk = .false.

    contains

    subroutine MPKLikelihood_Add(LikeList, Ini)
    use IniFile
    use settings
    class(LikelihoodList) :: LikeList
    Type(TIniFile) :: ini
    Type(MPKLikelihood), pointer :: like
    integer nummpksets, i


    use_mpk = (Ini_Read_Logical_File(Ini, 'use_mpk',.false.))

    if(.not. use_mpk) return

    call WiggleZLikelihood_Add(LikeList, Ini)

    nummpksets = Ini_Read_Int('mpk_numdatasets',0)
    do i= 1, nummpksets
        allocate(like)
        call ReadMpkDataset(like, ReadIniFileName(Ini,numcat('mpk_dataset',i)) )
        like%LikelihoodType = 'MPK'
        like%needs_powerspectra = .true.
        call LikeList%Add(like)
    end do
    if (Feedback>1) write(*,*) 'read mpk datasets'

    end subroutine MPKLikelihood_Add

    subroutine mpk_SetTransferRedshifts(redshifts)
    use wigglezinfo
    implicit none
    real(mcp), intent(inout) :: redshifts(matter_power_lnzsteps)
    !input is default log z spacing; can change here; check for consistency with other (e.g. lya)

    !Note internal ordering in CAMB is the opposite to that used in cosmomc transfer arrays (as here)
    !first index here must be redshift zero
    if(use_wigglez_mpk .and. matter_power_lnzsteps < 5) &
    call MpiStop('For WiggleZ MPK, matter_power_lnzsteps should be set to at least 5 (hardcoded in cmbtypes)')
    !wigglez redshifts: z0 = 0.d0, za = 0.22d0, zb = 0.41d0, zc = 0.6d0, zd = 0.78d0

    if(matter_power_lnzsteps==1 .or. .not. use_wigglez_mpk) return

    if(feedback>1) write(*,*) 'mpk_SetTransferRedshifts:  Reordering redshifts'

    izwigglez(1) = 2
    izwigglez(2) = 3
    izwigglez(3) = 4
    izwigglez(4) = 5
    redshifts(izwigglez(1)) = za
    redshifts(izwigglez(2)) = zb
    redshifts(izwigglez(3)) = zc
    redshifts(izwigglez(4)) = zd

    if(redshifts(1) > 0.001) call MpiStop('redshifts(1) should be at z=0!')
    return
    end subroutine mpk_SetTransferRedshifts

    subroutine ReadmpkDataset(mset,gname)
    use MatrixUtils
    type(MPKLikelihood) mset
    character(LEN=*), intent(IN) :: gname
    character(LEN=Ini_max_string_len) :: kbands_file, measurements_file, windows_file, cov_file

    integer i,iopb
    real(mcp) keff,klo,khi,beff
    integer :: num_mpk_points_full ! actual number of bandpowers in the infile
    integer :: num_mpk_kbands_full ! actual number of k positions " in the infile
    integer :: max_mpk_points_use ! in case you don't want the smallest scale modes (eg. sdss)
    integer :: min_mpk_points_use ! in case you don't want the largest scale modes
    integer :: max_mpk_kbands_use ! in case you don't want to calc P(k) on the smallest scales (will truncate P(k) to zero here!)
    integer :: min_mpk_kbands_use ! in case you don't want to calc P(k) on the largest scales (will truncate P(k) to zero here!)
    real(mcp), dimension(:,:), allocatable :: mpk_Wfull, mpk_covfull
    real(mcp), dimension(:), allocatable :: mpk_kfull, mpk_fiducial

    character(80) :: dummychar
    logical bad
    Type(TIniFile) :: Ini
    integer file_unit


    file_unit = new_file_unit()
    call Ini_Open_File(Ini, gname, file_unit, bad, .false.)
    if (bad) then
        write (*,*)  'Error opening dataset file '//trim(gname)
        stop
    end if

    mset%name = Ini_Read_String_File(Ini,'name')
    if (mset%name == 'twodf') stop 'twodf no longer supported - use data/2df_2005.dataset'
    if (mset%name == 'lrg_2009') stop 'lrg_2009 no longer supported'
    Ini_fail_on_not_found = .false.
    mset%use_set =.true.
    if (Feedback > 0) write (*,*) 'reading: '//trim(mset%name)
    num_mpk_points_full = Ini_Read_Int_File(Ini,'num_mpk_points_full',0)
    if (num_mpk_points_full.eq.0) write(*,*) ' ERROR: parameter num_mpk_points_full not set'
    num_mpk_kbands_full = Ini_Read_Int_File(Ini,'num_mpk_kbands_full',0)
    if (num_mpk_kbands_full.eq.0) write(*,*) ' ERROR: parameter num_mpk_kbands_full not set'
    min_mpk_points_use = Ini_Read_Int_File(Ini,'min_mpk_points_use',1)
    min_mpk_kbands_use = Ini_Read_Int_File(Ini,'min_mpk_kbands_use',1)
    max_mpk_points_use = Ini_Read_Int_File(Ini,'max_mpk_points_use',num_mpk_points_full)
    max_mpk_kbands_use = Ini_Read_Int_File(Ini,'max_mpk_kbands_use',num_mpk_kbands_full)
    mset%num_mpk_points_use = max_mpk_points_use - min_mpk_points_use +1
    mset%num_mpk_kbands_use = max_mpk_kbands_use - min_mpk_kbands_use +1

    allocate(mpk_Wfull(num_mpk_points_full,num_mpk_kbands_full))
    allocate(mpk_kfull(num_mpk_kbands_full))
    allocate(mset%mpk_P(mset%num_mpk_points_use))
    allocate(mset%mpk_sdev(mset%num_mpk_points_use))  ! will need to replace with the covmat
    allocate(mset%mpk_k(mset%num_mpk_kbands_use))
    allocate(mset%mpk_W(mset%num_mpk_points_use,mset%num_mpk_kbands_use))
    allocate(mpk_fiducial(mset%num_mpk_points_use))

    kbands_file  = ReadIniFileName(Ini,'kbands_file')
    call ReadVector(kbands_file,mpk_kfull,num_mpk_kbands_full)
    mset%mpk_k(1:mset%num_mpk_kbands_use)=mpk_kfull(min_mpk_kbands_use:max_mpk_kbands_use)
    if (Feedback > 1) then
        write(*,*) 'reading: '//trim(mset%name)//' data'
        write(*,*) 'Using kbands windows between',real(mset%mpk_k(1)),' < k/h < ',real(mset%mpk_k(mset%num_mpk_kbands_use))
    endif
    if  (mset%mpk_k(1) < matter_power_minkh) then
        write (*,*) 'WARNING: k_min in '//trim(mset%name)//'less than setting in cmbtypes.f90'
        write (*,*) 'all k<matter_power_minkh will be set to matter_power_minkh'
    end if

    measurements_file  = ReadIniFileName(Ini,'measurements_file')
    call OpenTxtFile(measurements_file, tmp_file_unit)
    mset%mpk_P=0.
    read (tmp_file_unit,*) dummychar
    read (tmp_file_unit,*) dummychar
    do i= 1, (min_mpk_points_use-1)
        read (tmp_file_unit,*, iostat=iopb) keff,klo,khi,beff,beff,beff
    end do
    if (Feedback > 1 .and. min_mpk_points_use>1) write(*,*) 'Not using bands with keff=  ',real(keff),' or below'
    do i =1, mset%num_mpk_points_use
        read (tmp_file_unit,*, iostat=iopb) keff,klo,khi,mset%mpk_P(i),mset%mpk_sdev(i),mpk_fiducial(i)
    end do
    close(tmp_file_unit)
    if (Feedback > 1) write(*,*) 'bands truncated at keff=  ',real(keff)

    windows_file  = ReadIniFileName(Ini,'windows_file')
    if (windows_file.eq.'') write(*,*) 'ERROR: mpk windows_file not specified'
    call ReadMatrix(windows_file,mpk_Wfull,num_mpk_points_full,num_mpk_kbands_full)
    mset%mpk_W(1:mset%num_mpk_points_use,1:mset%num_mpk_kbands_use)= &
    mpk_Wfull(min_mpk_points_use:max_mpk_points_use,min_mpk_kbands_use:max_mpk_kbands_use)


    cov_file  = ReadIniFileName(Ini,'cov_file')
    if (cov_file /= '') then
        allocate(mpk_covfull(num_mpk_points_full,num_mpk_points_full))
        call ReadMatrix(cov_file,mpk_covfull,num_mpk_points_full,num_mpk_points_full)
        allocate(mset%mpk_invcov(mset%num_mpk_points_use,mset%num_mpk_points_use))
        mset%mpk_invcov=  mpk_covfull(min_mpk_points_use:max_mpk_points_use,min_mpk_points_use:max_mpk_points_use)
        call Matrix_Inverse(mset%mpk_invcov)
        deallocate(mpk_covfull)
    else
        nullify(mset%mpk_invcov)
    end if

    mset%use_scaling = Ini_Read_Logical_File(Ini,'use_scaling',.false.)

    !JD 09/13 Read in fiducial D_V and redshift for use when calculating a_scl
    if(mset%use_scaling) then
        mset%redshift = Ini_Read_Double_File(Ini,'redshift',0.35d0)
        !DV_fid should be in units CMB%H0*BAO_D_v(z)
        mset%DV_fid = Ini_Read_Double_File(Ini,'DV_fid',-1.d0)
        if(mset%DV_fid == -1.d0) then
            write(*,*)'ERROR: use_scaling = T and no DV_fid given '
            write(*,*)'       for dataset '//trim(mset%name)//'.'
            write(*,*)'       Please check your .dataset files.'
            call MPIstop()
        end if
    end if

    mset%Q_marge = Ini_Read_Logical_File(Ini,'Q_marge',.false.)
    if (mset%Q_marge) then
        mset%Q_flat = Ini_Read_Logical_File(Ini,'Q_flat',.false.)
        if (.not. mset%Q_flat) then
            !gaussian prior on Q
            mset%Q_mid = Ini_Read_Real_File(Ini,'Q_mid')
            mset%Q_sigma = Ini_Read_Real_File(Ini,'Q_sigma')
        end if
        mset%Ag = Ini_Read_Real_File(Ini,'Ag', 1.4)
    end if
    if (iopb.ne.0) then
        stop 'Error reading mpk file'
    endif

    call Ini_Close_File(Ini)
    call ClearFileUnit(file_unit)

    deallocate(mpk_Wfull, mpk_kfull,mpk_fiducial)

    end subroutine ReadmpkDataset


    function MPK_LnLike(like,CMB,Theory,DataParams) ! LV_06 added CMB here
    Class(CMBParams) CMB
    Class(MPKLikelihood) :: like
    Class(TheoryPredictions) Theory
    real(mcp) :: DataParams(:)
    real(mcp) MPK_LnLike,LnLike
    real(mcp), dimension(:), allocatable :: mpk_Pth, mpk_k2,mpk_lin,k_scaled !LV_06 added for LRGDR4
    real(mcp), dimension(:), allocatable :: w
    real(mcp), dimension(:), allocatable :: mpk_WPth, mpk_WPth_k2
    real(mcp) :: covdat(like%num_mpk_points_use), covth(like%num_mpk_points_use),  covth_k2(like%num_mpk_points_use)
    real(mcp) :: normV, Q, minchisq
    real(mcp) :: a_scl  !LV_06 added for LRGDR4
    integer :: i, iQ
    logical :: do_marge
    integer, parameter :: nQ=6
    real(mcp) :: tmp, dQ = 0.4
    real(mcp) chisq(-nQ:nQ)
    real(mcp) calweights(-nQ:nQ)
    real(mcp) vec2(2),Mat(2,2)

    allocate(mpk_lin(like%num_mpk_kbands_use) ,mpk_Pth(like%num_mpk_kbands_use))
    allocate(mpk_WPth(like%num_mpk_points_use))
    allocate(k_scaled(like%num_mpk_kbands_use))!LV_06 added for LRGDR4
    allocate(w(like%num_mpk_points_use))

    chisq = 0

    if (.not. like%use_set) then
        LnLike = 0
        return
    end if

    !JD 09/13 new compute_scaling_factor functions
    if(like%use_scaling) then
        call compute_scaling_factor(like%redshift,CMB,like%DV_fid,a_scl)
    else
        a_scl = 1
    end if

    do i=1, like%num_mpk_kbands_use
        !Errors from using matter_power_minkh at lower end should be negligible
        k_scaled(i)=max(matter_power_minkh,a_scl*like%mpk_k(i))
        mpk_lin(i)=MatterPowerAt(Theory,k_scaled(i))/a_scl**3
    end do


    do_marge = like%Q_Marge
    if (do_marge .and. like%Q_flat) then
        !Marginalize analytically with flat prior on b^2 and b^2*Q
        !as recommended by Max Tegmark for SDSS
        allocate(mpk_k2(like%num_mpk_kbands_use))
        allocate(mpk_WPth_k2(like%num_mpk_points_use))

        mpk_Pth=mpk_lin/(1+like%Ag*k_scaled)
        mpk_k2=mpk_Pth*k_scaled**2
        mpk_WPth = matmul(like%mpk_W,mpk_Pth)
        mpk_WPth_k2 = matmul(like%mpk_W,mpk_k2)

        if (associated(like%mpk_invcov)) then
            covdat = matmul(like%mpk_invcov,like%mpk_P)
            covth = matmul(like%mpk_invcov,mpk_WPth)
            covth_k2 = matmul(like%mpk_invcov,mpk_WPth_k2)
        else
            w=1/(like%mpk_sdev**2)
            covdat = like%mpk_P*w
            covth = mpk_WPth*w
            covth_k2 = mpk_WPth_k2*w
        end if

        Mat(1,1) = sum(covth*mpk_WPth)
        Mat(2,2) = sum(covth_k2*mpk_WPth_k2)
        Mat(1,2) = sum(covth*mpk_WPth_k2)
        Mat(2,1) = Mat(1,2)
        LnLike = log( Mat(1,1)*Mat(2,2)-Mat(1,2)**2)
        call inv_mat22(Mat)
        vec2(1) = sum(covdat*mpk_WPth)
        vec2(2) = sum(covdat*mpk_WPth_k2)
        LnLike = (sum(like%mpk_P*covdat) - sum(vec2*matmul(Mat,vec2)) + LnLike ) /2

        deallocate(mpk_k2,mpk_WPth_k2)
    else
        if (like%Q_sigma==0) do_marge = .false.

        do iQ=-nQ,nQ
            Q = like%Q_mid +iQ*like%Q_sigma*dQ

            if (like%Q_marge) then
                mpk_Pth=mpk_lin*(1+Q*k_scaled**2)/(1+like%Ag*k_scaled)
            else
                mpk_Pth = mpk_lin
            end if

            mpk_WPth = matmul(like%mpk_W,mpk_Pth)

            !with analytic marginalization over normalization nuisance (flat prior on b^2)
            !See appendix F of cosmomc paper

            if (associated(like%mpk_invcov)) then
                covdat = matmul(like%mpk_invcov,like%mpk_P)
                covth = matmul(like%mpk_invcov,mpk_WPth)
                normV = sum(mpk_WPth*covth)
                chisq(iQ) = sum(like%mpk_P*covdat)  - sum(mpk_WPth*covdat)**2/normV  + log(normV)
            else
                !with analytic marginalization over normalization nuisance (flat prior on b^2)
                w=1/(like%mpk_sdev**2)
                normV = sum(mpk_WPth*mpk_WPth*w)
                tmp=sum(mpk_WPth*like%mpk_P*w)/normV ! avoid subtracting one large number from another
                chisq(iQ) = sum(like%mpk_P*(like%mpk_P - mpk_WPth*tmp)*w)  + log(normV)
            end if

            if (do_marge) then
                calweights(iQ) = exp(-(iQ*dQ)**2/2)
            else
                LnLike = chisq(iQ)/2
                exit
            end if

        end do

        !without analytic marginalization
        !! chisq = sum((like%mpk_P(:) - mpk_WPth(:))**2*w) ! uncommented for debugging purposes

        if (do_marge) then
            minchisq=minval(chisq)
            LnLike = sum(exp(-(chisq-minchisq)/2)*calweights)/sum(calweights)
            if (LnLike == 0) then
                LnLike = LogZero
            else
                LnLike =  -log(LnLike) + minchisq/2
            end if
        end if
    end if !not analytic over Q

    if (Feedback>1) write(*,*) 'mpk chi-sq:', LnLike*2

    if (LnLike > 1e8) then
        write(*,*) 'Chisq is huge, maybe there is a problem? chisq=',chisq
    end if

    deallocate(mpk_Pth,mpk_lin)
    deallocate(mpk_WPth,k_scaled,w)

    MPK_LnLike=LnLike

    end function MPK_LnLike

    subroutine inv_mat22(M)
    real(mcp) M(2,2), Minv(2,2), det

    det = M(1,1)*M(2,2)-M(1,2)*M(2,1)
    Minv(1,1)=M(2,2)
    Minv(2,2) = M(1,1)
    Minv(1,2) = - M(2,1)
    Minv(2,1) = - M(1,2)
    M = Minv/det

    end subroutine inv_mat22

    !-----------------------------------------------------------------------------
    ! JD 09/13: Replaced compute_scaling_factor routines so we use
    !           D_V calculations from CAMB.  New routines below

    subroutine compute_scaling_factor(z,CMB,DV_fid,a_scl)
    implicit none
    type(CMBParams) CMB
    real(mcp), intent(in) :: z, DV_fid
    real(mcp), intent(out) :: a_scl

    a_scl = DV_x_H0(z,CMB)/DV_fid
    !Like in original code, we need to apply a_scl in the correct direction
    a_scl = 1.0_mcp/a_scl
    end subroutine compute_scaling_factor

    function DV_x_H0(z,CMB)  !Want D_V*H_0
    use CAMB, only : BAO_D_v
    implicit none
    type(CMBParams) CMB
    real(mcp), intent(in) :: z
    real(mcp):: DV_x_H0

    !We calculate H_0*D_V because we dont care about scaling of h since
    !k is in units of h/Mpc
    DV_x_H0 = CMB%H0*BAO_D_v(z)

    end function DV_x_H0

    end module mpk
