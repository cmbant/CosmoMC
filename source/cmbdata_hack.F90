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

module cmbdata
use settings
use cmbtypes
use MatrixUtils
use CMBLikes
use constants
implicit none

 
!if CMBdataset%has_xfactors then obs, var and N_inv are for the offset-lognormal variables Z

  Type CMBdatapoint
    integer win_min, win_max
     !Ranges in which window is non-zero
    real, pointer, dimension(:,:) :: window
     !Normalized window function in l
    real obs, err_minus, err_plus, sigma, var
     !Observed value of l(l+1) C_l/2pi bandpowers in MicroK^2, with errors
    real beam_err !fractional beam error (file is value in MicroK^2)
    logical inc_pol

  end Type CMBdatapoint

  Type CMBdataset
    logical :: use_set
    logical :: has_pol, windows_are_bandpowers,windows_are_normalized
    logical :: has_sz_template
    real :: calib_uncertainty
    logical :: beam_uncertain, has_corr_errors, has_xfactors
    integer :: num_points, file_points
    character(LEN=80) :: name
    Type(CMBdatapoint), pointer, dimension(:) :: points
    real, pointer, dimension(:,:) :: N_inv
    real, pointer, dimension(:) :: xfactors
    logical, pointer, dimension(:) :: has_xfactor !whether each bin has one
    logical :: all_l_exact
    logical :: CMBLike !New format
    integer :: all_l_lmax
    integer :: nuisance_parameters
    real, pointer, dimension(:,:) :: all_l_obs, all_l_noise
    real, pointer, dimension(:) :: all_l_fsky
    real, pointer, dimension(:) :: sz_template
    Type (TCMBLikes), pointer :: CMBLikes
  end Type CMBdataset

  integer :: num_datasets = 0
  Type(CMBdataset) datasets(10)

  logical :: init_MAP = .true.

  integer :: cl_bin_width =1

 
  integer, parameter :: halfsteps = 5 !do 2*halfsteps+1 steps in numerical marginalization
  real margeweights(-halfsteps:halfsteps)
  real :: margenorm = 0 
  private halfsteps, margeweights, margenorm
contains


 subroutine ReadWindow(AP, aname, are_bare, aset)
   Type(CMBdatapoint) AP
   character(LEN=*), intent(IN) :: aname
   logical, intent(IN) :: are_bare
   Type (CMBdataset) :: aset
   integer l, ncls
   real wpol(1:num_cls-1),ll, IW,w
   character(LEN=200) tmp

   if (Feedback > 1) write (*,*) 'reading window: '//trim(aname)

   if (aset%has_pol) then
    ncls = num_cls
   else
    ncls = 1
   endif
   allocate(AP%window(ncls,lmax))

   AP%window = 0

   call OpenTxtFile(aname, tmp_file_unit)

     do
     if (aset%has_pol) then
       read(tmp_file_unit,'(a)',end=1) tmp
       read(tmp,*, end=1) ll, w, wpol
     else
       read(tmp_file_unit,*, end=1) ll, w
     end if

     l=nint(ll)
     if (abs(l-ll) > 1e-4) stop 'non-integer l in window file' 
     if (l>=2 .and. l<=lmax) then
        AP%window(1,l) = w
        if(aset%has_pol) then
         AP%window(2:num_cls,l) = wpol
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
1  close(tmp_file_unit)

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
      Type(TIniFile) :: Ini
      Type (CMBdataset) :: aset
      character(LEN=Ini_max_string_len) :: fname
      integer l, idum, ncls, ncol
      real inobs(4)
      integer file_unit

    !In this case we have data for every l, and use exact full-sky likelihood expression
    !with some fudge factor fsky_eff^2 to reduce the degrees of freedom: fsky^eff*(2l+1)
        
       if (Feedback > 0) &
        write(*,*) 'all_l_exact note: you might want to change fsky_eff^2 factor to fsky_eff' 

       aset%num_points = 0
       aset%all_l_lmax = Ini_Read_Int_File(Ini,'all_l_lmax')

       if (aset%all_l_lmax > lmax) stop 'cmbdata.f90::ReadAllExact: all_l_lmax > lmax'
       if (aset%has_pol) then
        allocate(aset%all_l_obs(2:aset%all_l_lmax,num_cls))
        allocate(aset%all_l_noise(2:aset%all_l_lmax,2))
       else
        allocate(aset%all_l_obs(2:aset%all_l_lmax,1))
        allocate(aset%all_l_noise(2:aset%all_l_lmax,1))
       end if
       allocate(aset%all_l_fsky(2:aset%all_l_lmax))
       fname = trim(Ini_Read_String_File(Ini,'all_l_file'))
       ncol = TxtFileColumns(fname)
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
       
       file_unit =  new_file_unit()

       call OpenTxtFile(fname,file_unit)
!File format:
!l C_TT (C_TE C_EE [C_BB]) N_T (N_P) fsky_eff
!Must have num_cls set to correct value for file
       do l = 2, aset%all_l_lmax
          read (file_unit, *, end=100, err=100) idum,inobs(1:ncls), aset%all_l_noise(l,:), aset%all_l_fsky(l)
          if (idum /= l) stop 'Error reading all_l_file'
            !set BB to pure noise if not in file
          if (aset%has_pol .and. ncls < num_cls) inobs(num_cls) = aset%all_l_noise(l,2)
          aset%all_l_obs(l,:) = inobs(1:num_cls)
        end do
       call CloseFile(file_unit)

       return
100    stop 'Error reading all_l_file file'
     

 end subroutine ReadAllExact


 function ChiSqExact(cl, aset) 
  !Compute -ln(Likelihood) 
   real cl(lmax,num_cls)
   Type(CMBdataset) :: aset
   integer l
   real ChiSqExact, chisq, term, CT, CE, CB, CC
   real CThat, CChat, CEhat, CBhat
   real dof
   integer i

   chisq=0
   
   do l=2, 30 
     dof = aset%all_l_fsky(l)*(2*l+1)
       !Ignoring l correlations but using f_sky^2_eff fudge factor may be a good approx
       !for nearly full sky observations
       !switched to just fsky**1 default Nov 09 since usually more useful
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
      dof = dof + aset%all_l_fsky(i)**2*(2*i+1)
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

 subroutine ReadDataset(aname)
   use CMBLikes
   character(LEN=*), intent(IN) :: aname
   character(LEN=Ini_max_string_len) :: InLine, window_dir, Ninv_file, xfact_file, band_file
   logical bad, windows_are_bare
   Type (CMBdataset) :: aset
   integer i, first_band, use_i
   real, pointer, dimension(:,:) :: tmp_mat
   real, pointer, dimension(:) :: tmp_arr
   character(LEN=Ini_max_string_len) :: data_format
   integer file_unit
   Type(TIniFile) :: Ini
   
   num_datasets = num_datasets + 1
   
   if (num_datasets > 10) stop 'too many datasets'

   aset%has_sz_template = .false.
   aset%nuisance_parameters = 0
   aset%CMBlike = .false.
   
!Special cases
   if (aname == 'MAP' .or. aname == 'WMAP') then 
     aset%name = 'WMAP'
     datasets(num_datasets) = aset
     return   
   elseif( aname(LEN_TRIM(aname)-5:LEN_TRIM(aname)) == 'newdat') then
    !Carlo format for polarized Boomerang et al.
     if (Feedback > 0) write(*,*) 'Reading BCP data set: ' // TRIM(aname)
     call ReadDataset_bcp(aset,aname)
     datasets(num_datasets) = aset
     return 
   end if

   file_unit = new_file_unit()
   
   call Ini_Open_File(Ini,aname, file_unit, bad, .false.)
   if (bad) then
     write (*,*)  'Error opening dataset file '//trim(aname)
     stop
   end if
   Ini_fail_on_not_found = .true.
      
   aset%name = Ini_Read_String_File(Ini,'name')
   aset%use_set =.true.
   aset%num_points = 0
    
   if (Feedback > 0) write (*,*) 'reading: '//trim(aset%name)

   Ini_fail_on_not_found = .false.

   data_format  = Ini_Read_String_File(Ini,'dataset_format')

   aset%CMBlike = data_format == 'CMBLike'

   aset%all_l_exact = (data_format =='all_l_exact') &
            .or. Ini_Read_Logical_File(Ini,'all_l_exact',.false.)
   if (aset%CMBLike) then
      allocate(aset%CMBLikes)
      call CMBLikes_ReadData(aset%CMBLikes, Ini, ExtractFilePath(aname))
      aset%nuisance_parameters = aset%CMBLikes%num_nuisance_parameters 
    
   else if (aset%all_l_exact) then
       aset%has_pol = Ini_Read_Logical_File(Ini,'has_pol',.false.)
       call ReadAllExact(Ini,aset)
   else if (data_format/='') then
        write(*,*) 'Error in '//trim(aname)
        write(*,*) 'Unknown data_format: '//trim(data_format)
        stop
   else
    !Otherwise do usual guassian/offset lognormal stuff

       aset%has_pol = Ini_Read_Logical_File(Ini,'has_pol',.false.)

       aset%num_points = Ini_Read_Int_File(Ini,'num_points')

       aset%calib_uncertainty = Ini_Read_Double_File(Ini,'calib_uncertainty')
       aset%beam_uncertain = Ini_Read_logical_File(Ini,'beam_uncertainty')
        
       window_dir  = ReadIniFilename(Ini,'window_dir')

       windows_are_bare = Ini_Read_Logical_File(Ini,'windows_are_bare',.false.)
       aset%windows_are_bandpowers = Ini_Read_Logical_File(Ini,'windows_are_bandpowers',.true.)
       aset%windows_are_normalized = Ini_Read_Logical_File(Ini,'windows_are_normalized',.false.)
  
       aset%file_points = Ini_read_Int_File(Ini,'file_points',aset%num_points)
       first_band = Ini_read_Int_File(Ini,'first_band',1)
       if (first_band + aset%num_points > aset%file_points+1) then
            write (*,*)  'Error with dataset file '//trim(aname)
            write (*,*) 'first_band + num_points > file_points'
            stop
       end if
    !Read in the observed values, errors and beam uncertainties
       allocate(aset%points(aset%num_points))
       band_file = Ini_Read_String_File(Ini,'bandpowers')
       if (band_file /= '') call OpenTxtFile(band_file, 51)
       Ini_fail_on_not_found = .true.
       do i=1, aset%num_points + first_band -1
          if (band_file /= '') then
            read(51,'(a)') InLine
          else
            InLine = Ini_Read_String_File(Ini,numcat('data',i))
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
          call ReadWindow(aset%points(use_i),trim(window_dir)//'/'//trim(numcat(aset%name,i)),windows_are_bare,aset)
       end do
       if (band_file /= '') Close(51)
  

    !See if the inverse covariance matrix is given (otherwise assume diagonal)
       Ini_fail_on_not_found = .false.
    
       Ninv_file = Ini_Read_String_File(Ini,'N_inv')
       aset%has_corr_errors = Ninv_file /= ''
       if (aset%has_corr_errors) then
          allocate(tmp_mat(aset%file_points,aset%file_points))
          allocate(aset%N_inv(aset%num_points,aset%num_points))
          call ReadMatrix(Ninv_file, tmp_mat,aset%file_points,aset%file_points)
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

       if (Ini_Read_Logical_File(Ini,'use_hyperparameter',.false.)) stop 'Hyperparameters deprecated'

    !See if xfactors are given

       xfact_file = Ini_Read_String_File(Ini,'xfactors')
       aset%has_xfactors = xfact_file /= ''

       if (aset%has_xfactors) then
          allocate(tmp_arr(aset%num_points + first_band -1))
          call ReadVector(xfact_file, tmp_arr,aset%num_points+ first_band -1)
          allocate(aset%xfactors(aset%num_points))
          allocate(aset%has_xfactor(aset%num_points))
          aset%has_xfactor = .true.
          aset%xfactors = tmp_arr(first_band:first_band+aset%num_points-1)
          deallocate(tmp_arr)
          aset%points(:)%var = aset%points(:)%var/(aset%points(:)%obs +aset%xfactors)**2
          aset%points(:)%obs = log(aset%points(:)%obs +aset%xfactors)
        
       end if

   end if !not all_l_exact or cut sky unbinned
   
   call Ini_Close_File(Ini)
   call ClearFileUnit(file_unit)

   datasets(num_datasets) = aset

 end subroutine ReadDataset

 subroutine ReadSZTemplate(aset, aname, ascale)
    Type (CMBdataset) :: aset
    real, intent(in) :: ascale
    character(LEN=*), intent(IN) :: aname
    integer l, unit
    real sz
     allocate(aset%sz_template(2:lmax))
     aset%sz_template = 0
     aset%has_sz_template = .true.
     call OpenTxtFile(aname, tmp_file_unit)
     do
       read(tmp_file_unit,*,end=2) l, sz
       if (l>=2 .and. l<=lmax) aset%sz_template(l) = ascale * sz/(l*(l+1)/twopi)
     end do
     
2    Close(tmp_file_unit)
 end subroutine ReadSZTemplate

 function GetWinBandPower(AP, cl)
   real  GetWinBandPower
   real cl(lmax,num_cls)
   Type(CMBdatapoint) AP
   integer l
   real bandpower

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
   real GetCalibMargexChisq
   Type(CMBdataset) aset
   real bandpowers(aset%num_points),beambandpowers(aset%num_points),diffs(aset%num_points)
   real calib, chisq(-halfsteps:halfsteps),chisqcalib(-halfsteps:halfsteps)
   real minchisq
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

      chisqcalib(ibeam) = -2*log(sum(margeweights*exp(max(-30.,-(chisq-minchisq)/2)))/margenorm) + minchisq

      if (.not. aset%beam_uncertain) then
         GetCalibMargexChisq = chisqcalib(ibeam)
         return
      end if

   end do

   minchisq = minval(chisqcalib)
   GetCalibMargexChisq = -2*log(sum(margeweights*exp(max(-30.,-(chisqcalib-minchisq)/2)))/margenorm) + minchisq
  
 end function GetCalibMargexChisq


  !Routine by Carlo Contaldi to read .newdat file format (Boomerang et al)
  !Modified to account for offset lognormal toggle per band
  !AL July 2005: modified to allow band selection
  !MLB May 09: modified to allow provision of per-band beam errors (Quad)
 SUBROUTINE ReadDataset_bcp(aset,aname)
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
   use constants
 
   CHARACTER(LEN=*), INTENT(IN) :: aname
   TYPE (CMBdataset) :: aset

   CHARACTER(LEN=100) :: instr
   CHARACTER(LEN=3), DIMENSION(1:6) :: ch_types

   LOGICAL windows_are_bare
   INTEGER i, j, k,use_i, n_types, xin
   REAL, POINTER, DIMENSION(:) :: tmp_x
   REAL, POINTER, DIMENSION(:,:) :: lb
   real, allocatable, dimension(:,:) :: tmp_mat
   integer, allocatable, dimension(:) :: used_bands

   INTEGER :: npol(6), minmax(2,6)
   INTEGER :: file_i,ijunk, ilike, ibeam
   REAL :: cal, beam_width, beam_sigma, l_mid
   integer, parameter :: like_xfactorall=1, like_xfactorsome = 2
   !to be compatible with some older CITA output files
   LOGICAL :: FISHER_T_CMB
   integer file_unit
   
   
   file_unit = new_file_unit()
   CALL OpenTxtFile(aname, file_unit)

   READ(file_unit,'(a)') instr
   FISHER_T_CMB = .FALSE.
   IF(instr == 'FISHER_T_CMB') THEN
      FISHER_T_CMB = .TRUE.
      READ(file_unit,'(a)') instr
      WRITE(*,'(a)') 'FISHER_T_CMB is set for :'//TRIM(ADJUSTL(instr))
   ENDIF
   aset%name = TRIM(ADJUSTL(instr))
   WRITE(*,*) 'Reading: '//TRIM(ADJUSTL(aset%name))
   aset%use_set =.TRUE.

   READ(file_unit,*) npol(1:6)

   aset%has_pol = any(npol(2:6) /=0)
 
   aset%all_l_exact = .FALSE.
   aset%file_points = SUM(npol)
   aset%num_points = SUM(npol)
   n_types = count(npol /= 0)
   READ(file_unit,'(a)') instr
   IF(instr == 'BAND_SELECTION') THEN
      !list of 'first_band last_band' for each pol type
      !if first_band=0 then ignore that pol type
      aset%num_points = 0
      aset%has_pol = .false.
      if(feedback>0) WRITE(*,*) 'Using selected band ranges' 
      do i=1,6
       READ(file_unit,*) minmax(1:2,i)
       if (minmax(1,i)/=0) then
        aset%num_points = aset%num_points + minmax(2,i) - minmax(1,i) + 1
        if (i>1) aset%has_pol = .true.
       else
        minmax(2,i) = 0
       end if
      end do
      READ(file_unit,'(a)') instr
   ELSE
     !use all bands in file
      do i=1,6
        minmax(1,i)=1
        minmax(2,i)=npol(i)
      end do  
   ENDIF


   READ(instr,*) ijunk, cal, aset%calib_uncertainty
   IF(ijunk == 0) aset%calib_uncertainty = 0.e0

   READ(file_unit,*) ibeam, beam_width, beam_sigma
   aset%beam_uncertain = ibeam /= 0

   !this agrees with latest windows coming out of MPIlikely
   windows_are_bare = .FALSE.
   aset%windows_are_bandpowers = .TRUE. 
   aset%windows_are_normalized = .TRUE.

   ALLOCATE(aset%points(aset%num_points))
   allocate(used_bands(aset%num_points))

   READ(file_unit,*) ilike
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
      READ(file_unit,'(a2)') ch_types(j)
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
                  READ(file_unit,'(a)') instr
                  READ(instr,*) ijunk, aset%points(use_i)%obs,aset%points(use_i)%err_minus, &
                       aset%points(use_i)%err_plus,tmp_x(use_i),lb(use_i,1),lb(use_i,2)
               ELSE
                  !like_xfactorsome
                  !read also offset switch per band
                  READ(file_unit,*) ijunk, aset%points(use_i)%obs,aset%points(use_i)%err_minus, &
                       aset%points(use_i)%err_plus,tmp_x(use_i),lb(use_i,1),lb(use_i,2),xin
                  aset%has_xfactor(use_i) = xin/=0
               ENDIF
               !beam error in bandpower
               l_mid = (lb(use_i,2)-lb(use_i,1))/2.d0 + lb(use_i,1)
               aset%points(use_i)%beam_err = exp(-l_mid*(l_mid+1.d0)*1.526e-8*2.d0*beam_sigma*beam_width)-1.d0 
               aset%points(use_i)%beam_err = abs(aset%points(use_i)%beam_err)
            ELSE !Bandpowers from file, a la Quad
               IF(ilike /= like_xfactorsome) THEN
                  READ(file_unit,'(a)') instr
                  READ(instr,*) ijunk, aset%points(use_i)%obs,aset%points(use_i)%err_minus, &
                       aset%points(use_i)%err_plus,tmp_x(use_i),lb(use_i,1),lb(use_i,2),aset%points(use_i)%beam_err
               ELSE 
                  !like_xfactorsome
                  !read also offset switch per band
                  READ(file_unit,*) ijunk, aset%points(use_i)%obs,aset%points(use_i)%err_minus, &
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
             CALL ReadWindow(aset%points(use_i), trim(concat(ExtractFilePath(aname),'windows/')) // &
                  TRIM(numcat(aset%name,file_i)),windows_are_bare,aset)

         else
          !discard band
          READ(file_unit,'(a)') instr
         end if
      ENDDO
      !discard correlation submatrix 
      READ(file_unit,'(a)') (instr,i=1,npol(k)) 
   ENDDO

   !assume always have the matrix
   aset%has_corr_errors = .TRUE.
   allocate(tmp_mat(aset%file_points,aset%file_points))
   ALLOCATE(aset%N_inv(aset%num_points,aset%num_points))
   READ(file_unit,*) (tmp_mat(1:aset%file_points,i),i=1,aset%file_points)
   aset%N_inv = tmp_mat(used_bands,used_bands)
   deallocate(tmp_mat)
   deallocate(used_bands)
 
   !READ(file_unit,*,err=101,end=101) instr
   !stop 'ReadDataset_bcp:Should be at end of file after reading matrix'
101 call CloseFile(file_unit)

   !recalibrate and change units as required
   !some older output had final fisher matrix in 
   !in units of T_CMB whitle bandpowers are in units
   !of \microK^2
   IF(FISHER_T_CMB) THEN
      aset%N_inv = cal**4 * aset%N_inv *COBE_CMBTemp**4 * 1.e24
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


 function CalcLnLike(clall, aset,nuisance_params) 
  !Compute -ln(Likelihood) 
   real clall(lmax,num_cls_tot), CalcLnLike
   real, intent(in) :: nuisance_params(:)
   Type(CMBdataset) aset
   integer i
   real cl(lmax, num_cls)
   real chisq
   real chi2op, chi2pp, wpp, wdd
   real chi2dd,chi2pd,chi2od
   real bandpowers(aset%num_points), diffs(aset%num_points), tmp(aset%num_points), beam(aset%num_points)
   real denom

   if (.not. aset%use_set) then
      CalcLnLike = 0
      return
   end if
   cl = clall(:,1:num_cls) !without lensing power spectrum or other extra CL

   if (aset%CMBLike) then

     chisq =  CMBLikes_CMBLike(aset%CMBLikes, clall, nuisance_params) 

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

 function CMBLnLike(cl, freq_params, nuisance_params)
  real, intent(in) ::  cl(lmax,num_cls_tot)
  real CMBLnLike
  real,intent(in) :: freq_params(num_freq_params),nuisance_params(:)
  real sznorm, szcl(lmax,num_cls_tot)
  integer i
  integer nuisance
  real tot(num_datasets)
  
  sznorm = freq_params(1)  
  nuisance =1
  do i=1, num_datasets
     szcl= cl
     if (datasets(i)%has_sz_template) then
      szcl(2:lmax,1) = szcl(2:lmax,1) + sznorm*datasets(i)%sz_template(2:lmax)  
     end if
     if (datasets(i)%name == 'WMAP') then
      tot(i) = MAPLnLike(szcl)
     else
       if (ubound(nuisance_params,1) < 1) then 
        tot(i) = CalcLnLike(szcl,datasets(i), nuisance_params)
       else
        tot(i) = CalcLnLike(szcl,datasets(i), nuisance_params(nuisance:))
       end if
      if (datasets(i)%CMBLike) nuisance = nuisance + datasets(i)%CMBLikes%num_nuisance_parameters
     end if
  end do
  CMBLnLike = SUM(tot) 
 end function


 function MAPLnLike(cl)
#ifndef NOWMAP
  use wmap_likelihood_7yr
  use WMAP_OPTIONS
  use WMAP_UTIL
#endif
  real cl(lmax,num_cls_tot), MAPLnLike
#ifndef NOWMAP

  real(8), dimension(2:ttmax) :: cl_tt,cl_te,cl_ee,cl_bb
  real(8)                     :: like(num_WMAP),like_tot
  integer  l

  if (Init_MAP) then
#ifdef WMAPNOHIGHLTT
   use_TT = .false.
   use_TT_beam_ptsrc = .false.   
#endif
   if (lmax<ttmax) stop 'lmax not large enough for WMAP'
   if (Feedback>0) write(*,*) 'reading WMAP7 data'
   Init_MAP = .false.
  end if

  
  do l = 2, ttmax
     cl_tt(l) = cl(l,1)*l*(l+1)/twopi
     cl_te(l) = cl(l,2)*l*(l+1)/twopi
     cl_ee(l) = cl(l,3)*l*(l+1)/twopi
     if(num_cls == 4) then
        cl_bb(l) = cl(l,num_cls)*l*(l+1)/twopi
     else
        cl_bb(l) = 0.0d0
     end if
  end do

  like=0.0d0
  call wmap_likelihood_compute(cl_tt,cl_te,cl_ee,cl_bb,like)
  !call wmap_likelihood_error_report
  
  if (wmap_likelihood_ok) then
     MAPLnLike = sum(like)
  else
     MAPLnLike = LogZero
  endif
#else
   MAPLnLike=cl(2,1) !just stop unused symbol warnings
   stop 'Compiled without WMAP'
#endif
 end function


!WMAP1
! function MAPLnLike(cl)
!  use WMAP
!  real cl(lmax,num_cls), MAPLnLike
! integer  l
!  real(WMAP_precision) clTT(WMAP_lmax_TT), clTE(WMAP_lmax_TT), ClEE(WMAP_lmax_TT)
!  integer stat 
!  character(LEN=20) :: ttFile = 'WMAP/tt_diag.dat'
!  character(LEN=20) :: ttOffDiag ='WMAP/tt_offdiag.dat'
! character(LEN=20) :: teFile   = 'WMAP/te_diag.dat'
!  character(LEN=20) :: teOffDiag ='WMAP/te_offdiag.dat'
! 
!  if (Init_MAP) then
!   if (lmax<WMAP_lmax_TT) stop 'lmax not large enough for WMAP'
!   if (Feedback>0) write(*,*) 'reading WMAP data'
!   Call WMAP_init(ttFile, ttOffDiag, teFile, teOffDiag, stat)!
!
!   if (stat /=0) stop 'Error reading WMAP files'
!   if (Feedback>0) write(*,*) 'WMAP read'
!   Init_MAP = .false.
!  end if

!  do l = 2, WMAP_lmax_TT
!      clTT(l) = cl(l,1)*l*(l+1)/twopi
!      clTE(l) = cl(l,2)*l*(l+1)/twopi
!      clEE(l) = cl(l,3)*l*(l+1)/twopi
!  end do
!  MAPLnLike = -(WMAP_LnLike_TT(clTT) + WMAP_LnLike_TE(clTT, clTE, clEE))

! end function


end module cmbdata


