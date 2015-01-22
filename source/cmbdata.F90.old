    !Module storing observed datasets, their points and window functions
    !and routines for computing the likelihood
    !Cls passed to these routines are in in MicroK^2, no 2pi or l(l+1)
    !WMAP data is handled separately as a special case

    !If windows_are_bandpowers=T Windows functions read in from files are assumed to be W_l/l.
    !(unless windows_are_bare=T, in which case we assume they are W_l)
    !The band powers are obtained from
    !W_l by \sum_l (l+1/2) W_l C_l/2pi where W_l is normalized so \sum_l (l+1/2)/l(l+1) W_l = 1.
    !If windows_are_normalized=T  the windows are assumed to be already normalized;
    !this must be the case for polarization in which the window files contain
    ! l TT TE EE BB
    !contributions to each band power. Usually all but two columns will be zeros.

    !If windows_are_bandpowers=F then window functions are assumed raw, and directly related
    !to the bandpowers by B_i=sum_l W_{il} C_l

    !This code is very indebted to Sarah's cosmarg
    !Analytic marginalization over calibration and beam from astro-ph/0112114

    !x factors are described in astro-ph/9808264. Essentially they are a first correction
    !for the bandpower errors being non-Gaussian. If no x-factors we assume they are gaussian.
    !Using x-factors we transform to a variable Z = log (bandpower + x), and assume Z has Gaussian errors.
    !When x-factors are used the inverse covariance matrix is assumed to be that for Z.
    !See the RADPACK page for data http://bubba.ucdavis.edu/~knox/radpack.html

    !Mar 04: added first_band parameter to .dataset files, added format for doing exact likelihoods
    !Jul 05: Readdataset_bcp changes for BOOM/CBI data, allowing for band cuts
    !Mar 06: changed to WMAP 3-year likelihood
    !Aug 06: added cal**2 scaling of x-factors
    !Oct 06: edited ReadAllExact to auto-account for number of cls (e.g. missing BB)
    !Oct 07: added Planck-like CMBLikes format
    !Sept 09: modified ReadDataset_bcp for QUaD, allowing beam errors to be explicitly provided for each band
    !         CMBLnLike passes array of parameters for frequency-dependent part of signal
    !Oct 27 Oct 09: fixed bugs using .newdat files
    !Jan 10: switched to support WMAP7
    !May 10: initialized CMBlike=.false.
    !Mar 11: all_l_exact likelihood uses fsky not fsky^2 everywhere; also consistent for datasetdir placeholders
    module cmbdata
    use settings
    use cmbtypes
    use CosmoTheory
    use MatrixUtils
    use CMBLikes
    use Likelihood_Cosmology
    implicit none

    logical :: Use_clik= .false.

    !if CMBdataset%has_xfactors then obs, var and N_inv are for the offset-lognormal variables Z

    Type CMBdatapoint
        integer win_min, win_max
        !Ranges in which window is non-zero
        real(mcp), pointer, dimension(:,:) :: window
        !Normalized window function in l
        real(mcp) obs, err_minus, err_plus, sigma, var
        !Observed value of l(l+1) C_l/2pi bandpowers in MicroK^2, with errors
        real(mcp) beam_err !fractional beam error (file is value in MicroK^2)
        logical inc_pol

    end Type CMBdatapoint

    Type CMBdataset
        logical :: use_set
        logical :: has_pol, windows_are_bandpowers,windows_are_normalized
        logical :: has_sz_template
        real(mcp) :: calib_uncertainty
        logical :: beam_uncertain, has_corr_errors, has_xfactors
        integer :: num_points, file_points
        character(LEN=:), allocatable :: dataset_filename
        Type(CMBdatapoint), pointer, dimension(:) :: points
        real(mcp), pointer, dimension(:,:) :: N_inv
        real(mcp), pointer, dimension(:) :: xfactors
        logical, pointer, dimension(:) :: has_xfactor !whether each bin has one
        logical :: all_l_exact
        logical :: CMBLike !New format
        integer :: all_l_lmax
        real(mcp), pointer, dimension(:,:) :: all_l_obs, all_l_noise
        real(mcp), pointer, dimension(:) :: all_l_fsky
        real(mcp), pointer, dimension(:) :: sz_template
        Type (TCMBLikes), pointer :: CMBLikes
    end Type CMBdataset

    type, extends(TCosmologyLikelihood) :: CMBDataLikelihood
        Type(CMBdataset) dataset
    contains
    procedure :: LogLike => CMBLnLike
    end type CMBDataLikelihood


    logical :: init_MAP = .true.

    integer :: cl_bin_width =1

    integer, private :: num_cls = 4


    integer, parameter :: halfsteps = 5 !do 2*halfsteps+1 steps in numerical marginalization
    real(mcp) margeweights(-halfsteps:halfsteps)
    real(mcp) :: margenorm = 0
    private halfsteps, margeweights, margenorm
    contains


    subroutine ReadWindow(AP, aname, are_bare, aset)
    Type(CMBdatapoint) AP
    character(LEN=*), intent(IN) :: aname
    logical, intent(IN) :: are_bare
    Type (CMBdataset) :: aset
    integer l, ncls
    real(mcp) wpol(1:3),ll, IW,w
    character(LEN=200) tmp
    Type(TTExtFile) :: F

    if (Feedback > 1) write (*,*) 'reading window: '//trim(aname)

    if (aset%has_pol) then
        ncls = 4
    else
        ncls = 1
    endif
    allocate(AP%window(ncls,lmax))

    AP%window = 0

    call F%Open(aname)

    do
        if (aset%has_pol) then
            read(F%unit,'(a)',end=1) tmp
            read(tmp,*, end=1) ll, w, wpol
        else
            read(F%unit,*, end=1) ll, w
        end if

        l=nint(ll)
        if (abs(l-ll) > 1e-4) stop 'non-integer l in window file'
        if (l>=2 .and. l<=lmax) then
            AP%window(1,l) = w
            if(aset%has_pol) then
                AP%window(2:4,l) = wpol
            end if
            if (.not. are_bare) AP%window(:,l) = AP%window(:,l)*l
        else
            if (l>lmax .and. w /= 0) then
                write (*,*) 'Warning: Window function non-zero outside l_max: ',trim(aname)
                write (*,*) 'assuming contribution is negligible'
                exit
            end if
        end if
    end do
1   call F%Close()

    do l=2, lmax
        if (any(AP%window(:,l)/=0)) then
            AP%win_min = l
            exit
        end if
    end do

    do l=lmax,2,-1
        if (any(AP%window(:,l)/=0)) then
            AP%win_max = l
            exit
        end if
    end do

    AP%inc_pol = .false.

    if (aset%windows_are_bandpowers) then
        IW = 0
        do l = AP%win_min, AP%win_max
            IW = IW + ((l+0.5)*AP%window(1,l))/(l*(l+1))
            AP%window(:,l) = AP%window(:,l)*(l+0.5)/(2*pi)
        end do
        if (.not. aset%windows_are_normalized) then
            if (aset%has_pol) stop 'Polarization windows must be normalized'
            AP%window(1,AP%win_min:AP%win_max) = &
            AP%window(1,AP%win_min:AP%win_max)/IW
        end if
    end if
    !If has_pol we are assuming windows are normalized

    if (aset%has_pol) AP%inc_pol = any(AP%window(2:num_cls,:)/=0)

    end subroutine ReadWindow

    subroutine ReadAllExact(Ini,aset)
    class(TSettingIni) :: Ini
    Type (CMBdataset) :: aset
    character(LEN=:), allocatable :: fname
    integer l, idum, ncls, ncol
    real(mcp) inobs(4)
    Type(TTextFile) :: F

    !In this case we have data for every l, and use exact full-sky likelihood expression
    !with some fudge factor fsky_eff^2 to reduce the degrees of freedom: fsky^eff*(2l+1)

    if (Feedback > 0) &
    write(*,*) 'all_l_exact note: all fsky_eff^2 changed to fsky_eff in this version'

    aset%num_points = 0
    aset%all_l_lmax = Ini%Read_Int('all_l_lmax')

    if (aset%all_l_lmax > lmax) stop 'cmbdata.f90::ReadAllExact: all_l_lmax > lmax'
    if (aset%has_pol) then
        allocate(aset%all_l_obs(2:aset%all_l_lmax,num_cls))
        allocate(aset%all_l_noise(2:aset%all_l_lmax,2))
    else
        allocate(aset%all_l_obs(2:aset%all_l_lmax,1))
        allocate(aset%all_l_noise(2:aset%all_l_lmax,1))
    end if
    allocate(aset%all_l_fsky(2:aset%all_l_lmax))
    fname = Ini%ReadFilename('all_l_file',File%ExtractPath(aset%dataset_filename))
    ncol = File%TxtColumns(fname)
    if (ncol==7) then
        ncls = 3
    elseif (ncol==8) then
        ncls = 4
    elseif (ncol==4) then
        ncls=1
        if (aset%has_pol) stop 'cmbdata.f90::ReadAllExact: has_pol wrong'
    else
        stop 'cmbdata.f90::ReadAllExact: wrong number of columns'
    end if

    call F%Open(fname)
    !File format:
    !l C_TT (C_TE C_EE [C_BB]) N_T (N_P) fsky_eff
    !Must have num_cls set to correct value for file
    do l = 2, aset%all_l_lmax
        read (F%unit, *, end=100, err=100) idum,inobs(1:ncls), aset%all_l_noise(l,:), aset%all_l_fsky(l)
        if (idum /= l) stop 'Error reading all_l_file'
        !set BB to pure noise if not in file
        if (aset%has_pol .and. ncls < num_cls) inobs(num_cls) = aset%all_l_noise(l,2)
        aset%all_l_obs(l,:) = inobs(1:num_cls)
    end do
    call F%Close()

    return
100 stop 'Error reading all_l_file file'


    end subroutine ReadAllExact


    function ChiSqExact(cl, aset)
    !Compute -ln(Likelihood)
    real(mcp) cl(CosmoSettings%lmax,num_cls)
    Type(CMBdataset) :: aset
    integer l
    real(mcp) ChiSqExact, chisq, term, CT, CE, CB, CC
    real(mcp) CThat, CChat, CEhat, CBhat
    real(mcp) dof
    integer i

    chisq=0

    do l=2, 30
        dof = aset%all_l_fsky(l)*(2*l+1)
        CT = cl(l,1) + aset%all_l_noise(l,1)
        if (aset%has_pol) then
            CE = cl(l,3) + aset%all_l_noise(l,2)
            term = CT*CE- cl(l,2)**2
            chisq = chisq + dof*( &
            (CT*aset%all_l_obs(l,3) + CE*aset%all_l_obs(l,1) - 2 *cl(l,2)*aset%all_l_obs(l,2))/term &
            + log( term/ (aset%all_l_obs(l,1)*aset%all_l_obs(l,3) - aset%all_l_obs(l,2)**2)) -2)
            if (num_cls>3) then
                !add in BB
                CB = cl(l,num_cls) + aset%all_l_noise(l,2)
                chisq = chisq + dof * (aset%all_l_obs(l,num_cls)/CB &
                +log(CB/aset%all_l_obs(l,num_cls)) - 1)
            end if
        else
            chisq = chisq + dof * (aset%all_l_obs(l,1)/CT &
            +log(CT/aset%all_l_obs(l,1)) - 1)
        end if
    end do


    do l=31 , aset%all_l_lmax, cl_bin_width
        dof = 0
        CT=0
        CE=0
        CC=0
        CB=0
        CThat=0
        CEhat=0
        CChat=0
        CBhat=0
        do i=l,l+cl_bin_width-1
            dof = dof + aset%all_l_fsky(i)*(2*i+1)
            !switched to just fsky here Apr 2011
            CT = CT + (cl(i,1) + aset%all_l_noise(i,1))*(2*i+1)
            CThat = CThat + aset%all_l_obs(i,1)*(2*i+1)
            if (aset%has_pol) then
                CE = CE + (cl(i,3) + aset%all_l_noise(i,2) ) *(2*i+1)
                CC = CC +  cl(i,2)*(2*i+1)
                CEhat = CEhat + aset%all_l_obs(i,3)*(2*i+1)
                CChat = CChat + aset%all_l_obs(i,2)*(2*i+1)
                if (num_cls>3) then
                    !add in BB
                    CB = CB + (cl(i,num_cls) + aset%all_l_noise(i,2))*(2*i+1)
                    CBhat = CBhat + aset%all_l_obs(i,num_cls)*(2*i+1)
                end if
            end if
        end do


        if (aset%has_pol) then
            term = CT*CE - CC**2
            chisq = chisq + dof*( &
            (CT*CEHat + CE*CThat - 2 *CC*CCHat)/term &
            + log( term/ (CTHat*CEHat - CCHat**2)) -2)
            if (num_cls>3) then
                !add in BB
                chisq = chisq + dof * (CBhat/CB +log(CB/CBhat) - 1)
            end if
        else
            chisq = chisq + dof * (CTHat/CT +log(CT/CTHat) - 1)
        end if
    end do

    ChiSqExact = ChiSq

    end function ChiSqExact

    subroutine ReadDataset(like,aname)
    use CMBLikes
    Type(CMBDataLikelihood), target :: like
    character(LEN=*), intent(IN) :: aname
    character(LEN=:), allocatable :: InLine, window_dir, Ninv_file, xfact_file, band_file
    logical bad, windows_are_bare
    Type (CMBdataset), pointer :: aset
    integer i, first_band, use_i
    real(mcp), pointer, dimension(:,:) :: tmp_mat
    real(mcp), pointer, dimension(:) :: tmp_arr
    character(LEN=:), allocatable :: data_format
    Type(TSettingIni) :: Ini
    Type(TTextFile) :: F

    aset=> like%dataset

    aset%has_sz_template = .false.
    aset%CMBlike = .false.
    aset%dataset_filename=aname

    !Special cases
    if (aname == 'MAP' .or. aname == 'WMAP') then
        like%name = 'WMAP'
        like%speed = -1
        return
    elseif( aname(LEN_TRIM(aname)-5:LEN_TRIM(aname)) == 'newdat') then
        !Carlo format for polarized Boomerang et al.
        if (Feedback > 0) write(*,*) 'Reading BCP data set: ' // TRIM(aname)
        call ReadDataset_bcp(like, aname)
        return
    end if

    call Ini%Open(aname,  bad, .false.)
    if (bad) then
        write (*,*)  'Error opening dataset file '//trim(aname)
        stop
    end if

    like%name = Ini%Read_String('name', .true.)
    aset%use_set =.true.
    aset%num_points = 0

    if (Feedback > 0) write (*,*) 'reading: '//trim(like%name)

    data_format  = Ini%Read_String('dataset_format')

    aset%CMBlike = data_format == 'CMBLike'

    aset%all_l_exact = (data_format =='all_l_exact') &
    .or. Ini%Read_Logical('all_l_exact',.false.)
    if (aset%CMBLike) then
        allocate(aset%CMBLikes)
        call aset%CMBLikes%ReadData(Ini, File%ExtractPath(aname))
    else if (aset%all_l_exact) then
        aset%has_pol = Ini%Read_Logical('has_pol',.false.)
        call ReadAllExact(Ini,aset)
    else if (data_format/='') then
        write(*,*) 'Error in '//trim(aname)
        write(*,*) 'Unknown data_format: '//trim(data_format)
        stop
    else
        !Otherwise do usual guassian/offset lognormal stuff

        aset%has_pol = Ini%Read_Logical('has_pol',.false.)

        aset%num_points = Ini%Read_Int('num_points')

        aset%calib_uncertainty = Ini%Read_Double('calib_uncertainty')
        aset%beam_uncertain = Ini%Read_logical('beam_uncertainty')

        window_dir  = Ini%ReadFilename('window_dir')

        windows_are_bare = Ini%Read_Logical('windows_are_bare',.false.)
        aset%windows_are_bandpowers = Ini%Read_Logical('windows_are_bandpowers',.true.)
        aset%windows_are_normalized = Ini%Read_Logical('windows_are_normalized',.false.)

        aset%file_points = Ini%read_Int('file_points',aset%num_points)
        first_band = Ini%read_Int('first_band',1)
        if (first_band + aset%num_points > aset%file_points+1) then
            write (*,*)  'Error with dataset file '//trim(aname)
            write (*,*) 'first_band + num_points > file_points'
            stop
        end if
        !Read in the observed values, errors and beam uncertainties
        allocate(aset%points(aset%num_points))
        band_file = Ini%Read_String('bandpowers')
        if (band_file /= '') call F%Open(band_file)
        do i=1, aset%num_points + first_band -1
            if (band_file /= '') then
                if (.not. F%ReadLine(InLine)) stop 'error reading bandpowers'
            else
                InLine = Ini%Read_String(numcat('data',i), .true.)
            end if
            if (i < first_band) cycle
            use_i = i - first_band + 1
            if (aset%beam_uncertain) then
                read(InLine,*) aset%points(use_i)%obs, aset%points(use_i)%err_minus, aset%points(use_i)%err_plus, &
                aset%points(use_i)%beam_err
                aset%points(use_i)%beam_err = aset%points(use_i)%beam_err/aset%points(use_i)%obs
            else
                read(InLine,*) aset%points(use_i)%obs, aset%points(use_i)%err_minus, aset%points(use_i)%err_plus
                aset%points(use_i)%beam_err = 0
            end if
            aset%points(use_i)%sigma = (aset%points(use_i)%err_minus + aset%points(use_i)%err_plus)/2
            aset%points(use_i)%var = aset%points(use_i)%sigma**2
            call ReadWindow(aset%points(use_i),trim(window_dir)//'/'//trim(numcat(like%name,i)),windows_are_bare,aset)
        end do
        if (band_file /= '') call F%Close()


        !See if the inverse covariance matrix is given (otherwise assume diagonal)

        Ninv_file = Ini%Read_String('N_inv')
        aset%has_corr_errors = Ninv_file /= ''
        if (aset%has_corr_errors) then
            allocate(tmp_mat(aset%file_points,aset%file_points))
            allocate(aset%N_inv(aset%num_points,aset%num_points))
            call File%ReadTextMatrix(Ninv_file, tmp_mat,aset%file_points,aset%file_points)
            if (aset%num_points /= aset%file_points) then
                !!Changed to truncation of covariance matrix, AL: May 2003
                call Matrix_Inverse(tmp_mat)
                aset%N_inv = tmp_mat(first_band:first_band+aset%num_points-1,&
                first_band:first_band + aset%num_points -1)
                call Matrix_Inverse(aset%N_inv)
            else
                aset%N_inv = tmp_mat(1:aset%num_points,1:aset%num_points)
            end if
            deallocate(tmp_mat)
        end if

        if (Ini%Read_Logical('use_hyperparameter',.false.)) stop 'Hyperparameters deprecated'

        !See if xfactors are given

        xfact_file = Ini%Read_String('xfactors')
        aset%has_xfactors = xfact_file /= ''

        if (aset%has_xfactors) then
            allocate(tmp_arr(aset%num_points + first_band -1))
            call File%ReadTextVector(xfact_file, tmp_arr,aset%num_points+ first_band -1)
            allocate(aset%xfactors(aset%num_points))
            allocate(aset%has_xfactor(aset%num_points))
            aset%has_xfactor = .true.
            aset%xfactors = tmp_arr(first_band:first_band+aset%num_points-1)
            deallocate(tmp_arr)
            aset%points(:)%var = aset%points(:)%var/(aset%points(:)%obs +aset%xfactors)**2
            aset%points(:)%obs = log(aset%points(:)%obs +aset%xfactors)
        end if

    end if !not all_l_exact or cut sky unbinned

    call Ini%Close()

    end subroutine ReadDataset

    subroutine ReadSZTemplate(aset, aname, ascale)
    Type (CMBdataset) :: aset
    real(mcp), intent(in) :: ascale
    character(LEN=*), intent(IN) :: aname
    integer l
    real(mcp) sz
    Type(TTextFile) :: F

    allocate(aset%sz_template(2:lmax))
    aset%sz_template = 0
    aset%has_sz_template = .true.
    call F%Open(aname)
    do
        read(F%unit,*,end=2) l, sz
        if (l>=2 .and. l<=lmax) aset%sz_template(l) = ascale * sz/(l*(l+1)/twopi)
    end do

2   call F%Close()
    end subroutine ReadSZTemplate

    function GetWinBandPower(AP, cl)
    real(mcp)  GetWinBandPower
    real(mcp) cl(CosmoSettings%lmax,num_cls)
    Type(CMBdatapoint) AP
    integer l
    real(mcp) bandpower

    bandpower = sum(cl(AP%win_min:AP%win_max,1)*AP%window(1,AP%win_min:AP%win_max))

    if (AP%inc_pol) then
        do l= AP%win_min, AP%win_max
            bandpower = bandpower + sum(cl(l,2:num_cls)*AP%window(2:num_cls,l))
        end do
    endif
    GetWinBandPower  = bandpower

    end function GetWinBandPower

    subroutine InitNumericalMarge
    integer i

    do i= -halfsteps, halfsteps
        margeweights(i) = exp(-(i*3/real(halfsteps))**2/2)
    end do
    margenorm = sum(margeweights)

    end subroutine InitNumericalMarge

    function GetCalibMargexChisq(bandpowers, aset)
    !Numerically integrate over the calibration uncertainty
    !Assume Gaussian prior, as for analytic calculation without x-factors
    !Could also Monte-Carlo
    real(mcp) GetCalibMargexChisq
    Type(CMBdataset) aset
    real(mcp) bandpowers(aset%num_points),beambandpowers(aset%num_points),diffs(aset%num_points)
    real(mcp) calib, chisq(-halfsteps:halfsteps),chisqcalib(-halfsteps:halfsteps)
    real(mcp) minchisq
    integer i,j, ibeam

    if (margenorm == 0) call InitNumericalMarge

    do ibeam= -halfsteps, halfsteps
        if (aset%beam_uncertain) then
            beambandpowers = bandpowers*(1 + aset%points(:)%beam_err*ibeam*3/real(halfsteps))!Go out to 3 sigma
        else
            beambandpowers = bandpowers
        end if

        do i=-halfsteps,halfsteps
            calib = 1 + aset%calib_uncertainty*i*3./halfsteps !Go out to 3 sigma

            if (aset%has_xfactors) then
                do j=1, aset%num_points
                    if (aset%has_xfactor(j)) then
                        diffs(j) = aset%points(j)%obs- log(calib*beambandpowers(j) + aset%xfactors(j))
                    else
                        diffs(j) = aset%points(j)%obs - calib*beambandpowers(j)
                    endif
                end do
            else
                diffs = aset%points(:)%obs - calib*beambandpowers
            end if

            if (aset%has_corr_errors) then
                chisq(i) = SUM(diffs*MATMUL(aset%N_inv,diffs))
            else
                chisq(i) = SUM(diffs**2/aset%points(:)%var)
            end if
        end do

        minchisq = minval(chisq)

        chisqcalib(ibeam) = -2*log(sum(margeweights*exp(max(-30._mcp,-(chisq-minchisq)/2)))/margenorm) + minchisq

        if (.not. aset%beam_uncertain) then
            GetCalibMargexChisq = chisqcalib(ibeam)
            return
        end if
    end do

    minchisq = minval(chisqcalib)
    GetCalibMargexChisq = -2*log(sum(margeweights*exp(max(-30._mcp,-(chisqcalib-minchisq)/2)))/margenorm) + minchisq

    end function GetCalibMargexChisq


    !Routine by Carlo Contaldi to read .newdat file format (Boomerang et al)
    !Modified to account for offset lognormal toggle per band
    !AL July 2005: modified to allow band selection
    !MLB May 09: modified to allow provision of per-band beam errors (Quad)
    SUBROUTINE ReadDataset_bcp(like,aname)
    !File Format:
    !name
    !n_bands_TT n_EE, n_BB, n_EB, n_TE, n_TB
    !has_calib_uncertain calib(amplitude) calib_err(power)
    !has_beam_uncertain beam beam_err
    !ilike (0: Gaussian, 1: all x-factor, 2: specified have x-factor)
    !loop over {
    ! band-types
    ! band info: num obs + - x_factor l_min l_max use_x
    ! correlation-matrix (ignored)
    ! }
    ! covariance matrix
    Type(CMBDataLikelihood), target :: like
    CHARACTER(LEN=*), INTENT(IN) :: aname
    TYPE (CMBdataset), pointer :: aset

    CHARACTER(LEN=100) :: instr
    CHARACTER(LEN=3), DIMENSION(1:6) :: ch_types

    LOGICAL windows_are_bare
    INTEGER i, j, k,use_i, n_types, xin
    real(mcp), POINTER, DIMENSION(:) :: tmp_x
    real(mcp), POINTER, DIMENSION(:,:) :: lb
    real(mcp), allocatable, dimension(:,:) :: tmp_mat
    integer, allocatable, dimension(:) :: used_bands

    INTEGER :: npol(6), minmax(2,6)
    INTEGER :: file_i,ijunk, ilike, ibeam
    real(mcp) :: cal, beam_width, beam_sigma, l_mid
    integer, parameter :: like_xfactorall=1, like_xfactorsome = 2
    !to be compatible with some older CITA output files
    LOGICAL :: FISHER_T_CMB
    Type(TTextFile) :: F

    aset=>like%dataset
    call F%Open(aname)

    READ(F%unit,'(a)') instr
    FISHER_T_CMB = .FALSE.
    IF(instr == 'FISHER_T_CMB') THEN
        FISHER_T_CMB = .TRUE.
        READ(F%unit,'(a)') instr
        WRITE(*,'(a)') 'FISHER_T_CMB is set for :'//TRIM(ADJUSTL(instr))
    ENDIF
    like%name = TRIM(ADJUSTL(instr))
    WRITE(*,*) 'Reading: '//TRIM(ADJUSTL(like%name))
    aset%use_set =.TRUE.

    READ(F%unit,*) npol(1:6)

    aset%has_pol = any(npol(2:6) /=0)

    aset%all_l_exact = .FALSE.
    aset%file_points = SUM(npol)
    aset%num_points = SUM(npol)
    n_types = count(npol /= 0)
    READ(F%unit,'(a)') instr
    IF(instr == 'BAND_SELECTION') THEN
        !list of 'first_band last_band' for each pol type
        !if first_band=0 then ignore that pol type
        aset%num_points = 0
        aset%has_pol = .false.
        if(feedback>0) WRITE(*,*) 'Using selected band ranges'
        do i=1,6
            READ(F%unit,*) minmax(1:2,i)
            if (minmax(1,i)/=0) then
                aset%num_points = aset%num_points + minmax(2,i) - minmax(1,i) + 1
                if (i>1) aset%has_pol = .true.
            else
                minmax(2,i) = 0
            end if
        end do
        READ(F%unit,'(a)') instr
    ELSE
        !use all bands in file
        do i=1,6
            minmax(1,i)=1
            minmax(2,i)=npol(i)
        end do
    ENDIF


    READ(instr,*) ijunk, cal, aset%calib_uncertainty
    IF(ijunk == 0) aset%calib_uncertainty = 0.e0

    READ(F%unit,*) ibeam, beam_width, beam_sigma
    aset%beam_uncertain = ibeam /= 0

    !this agrees with latest windows coming out of MPIlikely
    windows_are_bare = .FALSE.
    aset%windows_are_bandpowers = .TRUE.
    aset%windows_are_normalized = .TRUE.

    ALLOCATE(aset%points(aset%num_points))
    allocate(used_bands(aset%num_points))

    READ(F%unit,*) ilike
    aset%has_xfactors = ilike ==like_xfactorsome .or. ilike==like_xfactorall
    !1 : all bands are offset lognormal
    !2 : only bands specified have offset lognormal
    IF(aset%has_xfactors) then
        ALLOCATE(aset%has_xfactor(1:aset%num_points))
        aset%has_xfactor = .true.
    end if

    ALLOCATE(tmp_x(1:aset%num_points))
    ALLOCATE(lb(1:aset%num_points,1:2))
    k = 0
    use_i = 0
    file_i = 0
    DO j=1,n_types
        READ(F%unit,'(a2)') ch_types(j)
        k = k + 1
        DO i=1,20
            IF(npol(k) == 0 ) THEN
                k = k + 1
            ELSE
                EXIT
            ENDIF
        ENDDO

        if(feedback>1) write(*,*) j, ch_types(j), minmax(1,k), minmax(2,k)

        DO i=1,npol(k)
            file_i = file_i+1
            if (i>=minmax(1,k) .and. i<=minmax(2,k)) then
                use_i = use_i + 1
                used_bands(use_i) = file_i
                IF(ibeam .le. 1) THEN
                    IF(ilike /= like_xfactorsome) THEN
                        READ(F%unit,'(a)') instr
                        READ(instr,*) ijunk, aset%points(use_i)%obs,aset%points(use_i)%err_minus, &
                        aset%points(use_i)%err_plus,tmp_x(use_i),lb(use_i,1),lb(use_i,2)
                    ELSE
                        !like_xfactorsome
                        !read also offset switch per band
                        READ(F%unit,*) ijunk, aset%points(use_i)%obs,aset%points(use_i)%err_minus, &
                        aset%points(use_i)%err_plus,tmp_x(use_i),lb(use_i,1),lb(use_i,2),xin
                        aset%has_xfactor(use_i) = xin/=0
                    ENDIF
                    !beam error in bandpower
                    l_mid = (lb(use_i,2)-lb(use_i,1))/2.d0 + lb(use_i,1)
                    aset%points(use_i)%beam_err = exp(-l_mid*(l_mid+1.d0)*1.526e-8*2.d0*beam_sigma*beam_width)-1.d0
                    aset%points(use_i)%beam_err = abs(aset%points(use_i)%beam_err)
                ELSE !Bandpowers from file, a la Quad
                    IF(ilike /= like_xfactorsome) THEN
                        READ(F%unit,'(a)') instr
                        READ(instr,*) ijunk, aset%points(use_i)%obs,aset%points(use_i)%err_minus, &
                        aset%points(use_i)%err_plus,tmp_x(use_i),lb(use_i,1),lb(use_i,2),aset%points(use_i)%beam_err
                    ELSE
                        !like_xfactorsome
                        !read also offset switch per band
                        READ(F%unit,*) ijunk, aset%points(use_i)%obs,aset%points(use_i)%err_minus, &
                        aset%points(use_i)%err_plus,tmp_x(use_i),lb(use_i,1),lb(use_i,2),xin,aset%points(use_i)%beam_err
                        aset%has_xfactor(use_i) = xin/=0
                    ENDIF
                    l_mid = (lb(use_i,2)-lb(use_i,1))/2.d0 + lb(use_i,1)
                ENDIF
                aset%points(use_i)%sigma = (aset%points(use_i)%err_minus + aset%points(use_i)%err_plus)/2

                if (Feedback>1 ) print*, aset%beam_uncertain, l_mid, aset%points(use_i)%beam_err

                !recalibrate
                aset%points(use_i)%obs = cal**2 * aset%points(use_i)%obs
                aset%points(use_i)%sigma = cal**2 * aset%points(use_i)%sigma

                aset%points(use_i)%var = aset%points(use_i)%sigma**2
                !AL: Oct 08, changed to set path from the dataset path
                CALL ReadWindow(aset%points(use_i), trim(concat(File%ExtractPath(aname),'windows/')) // &
                TRIM(numcat(like%name,file_i)),windows_are_bare,aset)
            else
                !discard band
                READ(F%unit,'(a)') instr
            end if
        ENDDO
        !discard correlation submatrix
        READ(F%unit,'(a)') (instr,i=1,npol(k))
    ENDDO

    !assume always have the matrix
    aset%has_corr_errors = .TRUE.
    allocate(tmp_mat(aset%file_points,aset%file_points))
    ALLOCATE(aset%N_inv(aset%num_points,aset%num_points))
    READ(F%unit,*) (tmp_mat(1:aset%file_points,i),i=1,aset%file_points)
    aset%N_inv = tmp_mat(used_bands,used_bands)
    deallocate(tmp_mat)
    deallocate(used_bands)

    !READ(file_unit,*,err=101,end=101) instr
    !stop 'ReadDataset_bcp:Should be at end of file after reading matrix'
101 call F%Close()

    !recalibrate and change units as required
    !some older output had final fisher matrix in
    !in units of T_CMB whitle bandpowers are in units
    !of \microK^2
    IF(FISHER_T_CMB) THEN
        !      aset%N_inv = cal**4 * aset%N_inv *COBE_CMBTemp**4 * 1.e24
        call MpiStop('FISHER_T_CMB deprecated')
    ELSE
        aset%N_inv = cal**4 * aset%N_inv
    ENDIF

    IF(aset%has_pol) WRITE(*,*) 'has pols: ', ch_types(1:n_types) !removed ADJUSTR, maybe avoid PG fortran

    !transform into Z_B = ln(C_B+x_B) for bandpowers with xfactors
    IF(aset%has_xfactors) THEN
        ALLOCATE(aset%xfactors(1:aset%num_points))
        aset%xfactors(1:aset%num_points) = cal**2*tmp_x(1:aset%num_points)
        DO i=1,aset%num_points
            DO j=1,aset%num_points
                if(aset%has_xfactor(i)) aset%N_inv(i,j) =  aset%N_inv(i,j)/(aset%points(i)%obs + aset%xfactors(i))
                if(aset%has_xfactor(j)) aset%N_inv(i,j) =  aset%N_inv(i,j)/(aset%points(j)%obs + aset%xfactors(j))
            ENDDO
        ENDDO
        DO i=1,aset%num_points
            IF(aset%has_xfactor(i)) THEN
                aset%points(i)%var = aset%points(i)%var/(aset%points(i)%obs +aset%xfactors(i))**2
                aset%points(i)%obs  = LOG(aset%points(i)%obs + aset%xfactors(i))
            ENDIF
        ENDDO
        CALL Matrix_Inverse(aset%N_inv)
    ELSE
        CALL Matrix_Inverse(aset%N_inv)
    ENDIF

    DEALLOCATE(tmp_x)
    deallocate(lb)
    END SUBROUTINE ReadDataset_bcp


    function CalcLnLike(clall, aset)
    !Compute -ln(Likelihood)
    real(mcp) clall(lmax,num_cls_tot), CalcLnLike
    Type(CMBdataset) aset
    integer i
    real(mcp) cl(lmax, num_cls)
    real(mcp) chisq
    real(mcp) chi2op, chi2pp, wpp, wdd
    real(mcp) chi2dd,chi2pd,chi2od
    real(mcp) bandpowers(aset%num_points), diffs(aset%num_points), tmp(aset%num_points), beam(aset%num_points)
    real(mcp) denom

    if (.not. aset%use_set) then
        CalcLnLike = 0
        return
    end if
    cl = clall(:,1:num_cls) !without lensing power spectrum or other extra CL

    if (aset%CMBLike) then
        chisq =  aset%CMBLikes%CMBLike(clall)
    else if (aset%all_l_exact) then
        chisq = ChiSqExact(cl,aset)
    else
        denom = 1 !Assume Prob \propto exp(-chisq/2)/sqrt(denom)

        do i=1, aset%num_points
            bandpowers(i) = GetWinBandPower(aset%points(i), cl)
        end do

        if (aset%has_xfactors .and. (aset%calib_uncertainty > 1e-4 .or. aset%beam_uncertain))  then
            chisq = GetCalibMargexChisq(bandpowers,aset)
        else
            if (aset%has_xfactors) then
                do i=1, aset%num_points
                    if (aset%has_xfactor(i)) then
                        !obs in this case is Z = log(observed + x)
                        diffs(i) = aset%points(i)%obs- log(bandpowers(i) + aset%xfactors(i))
                    else
                        diffs(i) = aset%points(i)%obs - bandpowers(i)
                    end if
                end do
            else
                diffs = aset%points(:)%obs - bandpowers
            end if

            if (aset%has_corr_errors) then
                chisq = SUM(diffs*MATMUL(aset%N_inv,diffs))
            else
                chisq = SUM(diffs**2/aset%points(:)%var)
            end if

            if (aset%calib_uncertainty > 1e-4 .or. aset%beam_uncertain) then
                if (aset%has_corr_errors) then
                    tmp = MATMUL(aset%N_inv,bandpowers)
                else
                    tmp = bandpowers/aset%points(:)%var
                end if
                chi2op = SUM(diffs*tmp)
                chi2pp = SUM(bandpowers*tmp)
                if (aset%beam_uncertain) then
                    beam = aset%points(:)%beam_err*bandpowers
                    if (aset%has_corr_errors) then
                        tmp = MATMUL(aset%N_inv,beam)
                    else
                        tmp = beam/aset%points(:)%var
                    end if
                    chi2dd = SUM(beam*tmp)
                    chi2pd = SUM(bandpowers*tmp)
                    chi2od = SUM(diffs*tmp)
                end if

                if (aset%calib_uncertainty > 1e-4) then
                    !analytic marginalization over calibration uncertainty
                    wpp = 1/(chi2pp+1/aset%calib_uncertainty**2)
                    chisq = chisq - wpp*chi2op**2
                    denom = denom/wpp*aset%calib_uncertainty**2
                else
                    wpp = 0
                end if

                if (aset%beam_uncertain) then
                    !analytic marginalization over beam uncertainty
                    wdd=1/(chi2dd-wpp*chi2pd**2+1)
                    chisq = chisq - wdd*(chi2od-wpp*chi2op*chi2pd)**2
                    denom = denom/wdd
                end if
            end if
        end if

        if (denom /= 1) chisq = chisq + log(denom)
    end if
    CalcLnLike = chisq/2

    end function CalcLnLike


    function CMBLnLike(like, CMB, Theory, DataParams)
    Class(CMBDataLikelihood) :: like
    Class (CMBParams) CMB
    Class(TCosmoTheoryPredictions), target :: Theory
    real(mcp) :: DataParams(:)
    real(mcp) cl(CosmoSettings%lmax,num_cls_tot)
    real(mcp) CMBLnLike
    real(mcp) sznorm, szcl(lmax,num_cls_tot)

    call Theory%ClArray(cl(:,1),1,1)
    call Theory%ClArray(cl(:,2),2,1)
    call Theory%ClArray(cl(:,3),2,2)
    call Theory%ClArray(cl(:,4),3,3)

    call Theory%ClsFromTheoryData(cl)

    szcl= cl
    if (like%dataset%has_sz_template) then
        sznorm = DataParams(1)
        szcl(2:lmax,1) = szcl(2:lmax,1) + sznorm*like%dataset%sz_template(2:lmax)
    end if
    if (like%name == 'WMAP') then
        CMBLnLike = MAPLnLike(szcl)
    else
        CMBLnLike = CalcLnLike(szcl,like%dataset)
    end if

    end function


    end module cmbdata
