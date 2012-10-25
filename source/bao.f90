!!! Generalized BAO module added by J. Dossett
! Copied structure from mpk.f90 and Reid BAO code
!
! When using WiggleZ data set cite Blake et al. arXiv:1108.2635 
!!!!!!!!
!for SDSS data set: http://www.sdss3.org/science/boss_dr9_final_results/README_sdssbaodr9iso
!! for explanation of the changes to the rs expression, see Hamann et al, 
!! http://arxiv.org/abs/1003.3999
!
!AL/JH Oct 2012: encorporate DR9 data into something close to new cosmomc format


module bao
use cmbtypes
use CAMB, only : AngularDiameterDistance  !!angular diam distance also in Mpc no h units
use constants
use Precision
 implicit none
 
  Type baodataset
    logical :: use_set
    integer :: num_bao ! total number of points used
    integer :: type_bao  
     !what type of bao data is used 
     !1: old sdss, no longer
     !2: A(z) =2 ie WiggleZ
     !3 D_V/rs in fitting forumla appox (DR9)
     !4: D_v - 6DF
    character(LEN=20) :: name
    real(dl), pointer, dimension(:) :: bao_z, bao_obs, bao_err
    real(dl), pointer, dimension(:,:) :: bao_invcov
   
  end Type baodataset
  integer :: num_bao_datasets = 0
  Type(baodataset) baodatasets(10)

  real(dl), dimension(10000) :: DR7_alpha_file, DR7_prob_file
  real(dl) DR7_dalpha
  integer DR7_alpha_npoints

contains 

!JD copied structure from mpk.f90
subroutine ReadBaoDataset(gname)   
    use MatrixUtils
    character(LEN=*), intent(IN) :: gname
    character(LEN=Ini_max_string_len) :: bao_measurements_file, bao_invcov_file
    type (baodataset) bset
    integer i,iopb
    logical bad
    Type(TIniFile) :: Ini
    integer file_unit
    
 
    num_bao_datasets = num_bao_datasets + 1
    if (num_bao_datasets > 10) stop 'too many datasets'
    file_unit = new_file_unit()
    call Ini_Open_File(Ini, gname, file_unit, bad, .false.)
    if (bad) then
      write (*,*)  'Error opening data set file '//trim(gname)
      stop
    end if
    
    bset%name = Ini_Read_String_File(Ini,'name')
    
    Ini_fail_on_not_found = .false.
    bset%use_set =.true.
    if (Feedback > 0) write (*,*) 'reading BAO data set: '//trim(bset%name)
    bset%num_bao = Ini_Read_Int_File(Ini,'num_bao',0)
    if (bset%num_bao.eq.0) write(*,*) ' ERROR: parameter num_bao not set'
    bset%type_bao = Ini_Read_Int_File(Ini,'type_bao',1)
    if(bset%type_bao /= 3 .and. bset%type_bao /=2 .and. bset%type_bao /=4) then
        write(*,*) bset%type_bao
        write(*,*)'ERROR: Invalid bao type specified in BAO dataset: '//trim(bset%name)
        call MPIStop()
    end if

   allocate(bset%bao_z(bset%num_bao))
   allocate(bset%bao_obs(bset%num_bao))
   allocate(bset%bao_err(bset%num_bao))

    bao_measurements_file = ReadIniFileName(Ini,'bao_measurements_file')
    call OpenTxtFile(bao_measurements_file, tmp_file_unit)
    do i=1,bset%num_bao
        read (tmp_file_unit,*, iostat=iopb) bset%bao_z(i),bset%bao_obs(i),bset%bao_err(i)
    end do
    close(tmp_file_unit) 

    if (bset%name == 'DR7') then
     !don't used observed value, probabilty distribution instead
      call BAO_DR7_init(ReadIniFileName(Ini,'prob_dist'))    
    else

       allocate(bset%bao_invcov(bset%num_bao,bset%num_bao))
       bset%bao_invcov=0
        
       if (Ini_HasKey_File(Ini,bao_invcov_file)) then

        bao_invcov_file  = ReadIniFileName(Ini,'bao_invcov_file')
        call OpenTxtFile(bao_invcov_file, tmp_file_unit)
        do i=1,bset%num_bao
            read (tmp_file_unit,*, iostat=iopb) bset%bao_invcov(i,:)
        end do
        close(tmp_file_unit)

        if (iopb.ne.0) then
           call MpiStop('Error reading bao file '//trim(bao_invcov_file))
        endif
        
       else
        do i=1,bset%num_bao
           !diagonal, or actually just 1..
            bset%bao_invcov(i,i) = 1/bset%bao_err(i)**2
        end do
       end if
        
    end if
    
   call Ini_Close_File(Ini)
   call ClearFileUnit(file_unit)

   baodatasets(num_bao_datasets) = bset

end subroutine ReadBaoDataset


!JD copied from Reid BAO code
function CMBToBAOrs()
   use settings
   use cmbtypes
   use ModelParams
   use Precision
   use ThermoData, only : z_drag
   implicit none
   real(dl) ::  adrag, atol, rsdrag
   real(dl), external :: dsoundda, rombint
   real(dl) :: CMBToBAOrs
   
   adrag = 1.0d0/(1.0d0+z_drag)
   atol = 1e-6
   rsdrag = rombint(dsound_da,1d-8,adrag,atol)
   CMBToBAOrs = rsdrag

end function CMBToBAOrs

function D_v(z)
    real(dl), intent(IN) :: z
    real(dl)  D_v, Hz, ADD
    real(dl) dtauda 
    external dtauda !CAMB's d eta/da

    ADD = AngularDiameterDistance(z)*(1.d0+z)
    Hz = 1/dtauda(1/(1._dl+z))*(1._dl+z)**2
    D_v = ((ADD)**2.d0*z/Hz)**(1.d0/3.d0)
end function D_v

function Acoustic(CMB,z)
    Type(CMBParams) CMB
    real(dl) Acoustic
    real(dl), intent(IN) :: z
    real(dl) omh2,ckm,omegam,h
    omegam = 1.d0 - CMB%omv - CMB%omk
    h = CMB%h0/100
    ckm = c/1e3_dl !JD c in km/s
    
    omh2 = omegam*h**2.d0
    Acoustic = 100*D_v(z)*sqrt(omh2)/(ckm*z)
end function Acoustic

function SDSS_dvtors(CMB,z)
!This uses numerical value of D_v/r_s, but re-scales it to match definition of SDSS
!paper fitting at the fiducial model. Idea being it is also valid for e.g. varying N_eff
    Type(CMBParams) CMB
    real(dl) SDSS_dvtors
    real(dl), intent(IN)::z
    real(dl) rs
    real(dl), parameter :: rs_rescale = 153.017d0/149.0808
    
!    rs = SDSS_CMBToBAOrs(CMB)
     rs = CMBToBAOrs()*rs_rescale !rescaled to match fitting formula for LCDM
     SDSS_dvtors = D_v(z)/rs

end function SDSS_dvtors


!===================================================================================

function BAO_LnLike(CMB)
   implicit none
   Type(CMBParams) CMB
   integer i,j,k
   real(dl)  BAO_LnLike
   real(dl) tot(num_bao_datasets)
   real(dl), allocatable :: BAO_theory(:)
   
   tot=0
   do i=1,num_bao_datasets
        if (baodatasets(i)%name=='DR7') then
            tot(i) = BAO_DR7_loglike(CMB,baodatasets(i)%bao_z(1))
        else

        allocate(BAO_theory(baodatasets(i)%num_bao))
    
        if(baodatasets(i)%type_bao ==3)then
            do j=1, baodatasets(i)%num_bao
                BAO_theory(j) = SDSS_dvtors(CMB,baodatasets(i)%bao_z(j))
            end do
        else if(baodatasets(i)%type_bao ==2)then
            do j=1, baodatasets(i)%num_bao
                BAO_theory(j) = Acoustic(CMB,baodatasets(i)%bao_z(j))
            end do
        else if(baodatasets(i)%type_bao ==4)then
            do j=1, baodatasets(i)%num_bao
                BAO_theory(j) = D_v(baodatasets(i)%bao_z(j))
            end do
        end if

        do j=1, baodatasets(i)%num_bao
            do k=1, baodatasets(i)%num_bao
                tot(i) = tot(i) +&
                        (BAO_theory(j)-baodatasets(i)%bao_obs(j))* &
                        baodatasets(i)%bao_invcov(j,k)*&
                        (BAO_theory(k)-baodatasets(i)%bao_obs(k))
            end do
        end do
        tot(i) = tot(i)/2.d0
        deallocate(BAO_theory)
        end if
       
       if(feedback>1)write(*,*)'Bao dataset: '//trim(baodatasets(i)%name)//' LnLike = ',tot(i)
   end do
   if (any(tot==logZero)) then
       BAO_LnLike = logZero
   else
       BAO_LnLike = sum(tot)
   end if
   if(feedback>1)write(*,*)'Bao_LnLike = ', Bao_LnLike

end function BAO_LnLike


 subroutine BAO_DR7_init(fname)
 character(LEN=*), intent(in) :: fname
 real(dl) :: tmp0,tmp1,DR7_alpha
 integer ios,ii

  open(unit=7,file=fname,status='old')
  !Read data file
  ios = 0
  ii  = 0
  do while (ios.eq.0)
    read (7,*,iostat=ios) tmp0,tmp1
    if (ios .ne. 0) cycle
    if((ii.gt.1).and.(abs(DR7_dalpha-(tmp0-DR7_alpha)).gt.1e-6)) then
       stop 'binning should be uniform in sdss_baoDR7.txt'
    endif
    ii = ii+1
    DR7_alpha_file(ii) = tmp0
    DR7_prob_file (ii) = tmp1
    DR7_dalpha = tmp0-DR7_alpha
    DR7_alpha  = tmp0
  enddo
  DR7_alpha_npoints = ii
  if (ii.eq.0) call MpiStop('ERROR : reading file')
  close(7)
  !Normalize distribution (so that the peak value is 1.0)
  tmp0=0.0
  do ii=1,DR7_alpha_npoints
    if(DR7_prob_file(ii).gt.tmp0) then
       tmp0=DR7_prob_file(ii)
    endif
  enddo
  DR7_prob_file=DR7_prob_file/tmp0

end subroutine BAO_DR7_init

 function BAO_DR7_loglike(CMB,z)
 Type(CMBParams) CMB
 real (dl) z, BAO_DR7_loglike, alpha_chain, prob
 real,parameter :: rs_wmap7=152.7934d0,dv1_wmap7=1340.177  !r_s and D_V computed for wmap7 cosmology
 integer ii
             alpha_chain = (SDSS_dvtors(CMB,z))/(dv1_wmap7/rs_wmap7)
             if ((alpha_chain.gt.DR7_alpha_file(DR7_alpha_npoints-1)).or.(alpha_chain.lt.DR7_alpha_file(1))) then
                BAO_DR7_loglike = logZero
             else
                ii=1+floor((alpha_chain-DR7_alpha_file(1))/DR7_dalpha)
                prob=DR7_prob_file(ii)+(DR7_prob_file(ii+1)-DR7_prob_file(ii))/ &
                         (DR7_alpha_file(ii+1)-DR7_alpha_file(ii))*(alpha_chain-DR7_alpha_file(ii))
                BAO_DR7_loglike = -log( prob )
             endif
 
 end function BAO_DR7_loglike

 function SDSS_CMBToBAOrs(CMB)
   use settings
   use cmbtypes
   use ModelParams
   use Precision
   implicit none
   Type(CMBParams) CMB
   real(dl) ::  rsdrag
   real(dl) :: SDSS_CMBToBAOrs
   real(dl) :: zeq,zdrag,omh2,obh2,b1,b2
   real(dl) :: rd,req,wkeq

   obh2=CMB%ombh2
   omh2=CMB%ombh2+CMB%omdmh2

   b1     = 0.313*omh2**(-0.419)*(1+0.607*omh2**0.674)
   b2     = 0.238*omh2**0.223
   zdrag  = 1291.*omh2**0.251*(1.+b1*obh2**b2)/(1.+0.659*omh2**0.828)
   zeq    = 2.50e4*omh2*(2.726/2.7)**(-4.)
   wkeq   = 7.46e-2*omh2*(2.726/2.7)**(-2)
   req    = 31.5*obh2*(2.726/2.7)**(-4)*(1e3/zeq)
   rd     = 31.5*obh2*(2.726/2.7)**(-4)*(1e3/zdrag)
   rsdrag = 2./(3.*wkeq)*sqrt(6./req)*log((sqrt(1.+rd)+sqrt(rd+req))/(1.+sqrt(req)))

   SDSS_CMBToBAOrs = rsdrag

 end function SDSS_CMBToBAOrs
 

end module bao