!!! Generalized BAO module added by J. Dossett
! Copied structure from mpk.f90 and Reid BAO code
!
! When using WiggleZ data set cite Blake et al. arXiv:1108.2635 
!!!!!!!!
!for SDSS data set:
!default values from http://arxiv.org/abs/0907.1660
!! for explanation of the changes to the rs expression, see Hamann et al, 
!! http://xxx.lanl.gov/abs/1003.3999



module bao
use cmbtypes
use CAMB, only : AngularDiameterDistance  !!angular diam distance also in Mpc no h units
use constants
use Precision
 implicit none
 
  Type baodataset
    logical :: use_set
    integer :: num_bao ! total number of points used
    integer :: type_bao  !what type of bao data is used rs/D_v =1 ie SDSS; A(z) =2 ie WiggleZ
    character(LEN=20) :: name
    real(dl), pointer, dimension(:) :: bao_z, bao_obs
    real(dl), pointer, dimension(:,:) :: bao_invcov
   
   end Type baodataset
  integer :: num_bao_datasets = 0
  Type(baodataset) baodatasets(10)
 

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
   	if(bset%type_bao /= 1 .and. bset%type_bao /=2 ) then
   		write(*,*) bset%type_bao
   		write(*,*)'ERROR: Invalid bao type specified in BAO dataset: '//trim(bset%name)
   		call MPIStop()
   	end if
   	

    allocate(bset%bao_invcov(bset%num_bao,bset%num_bao))
    allocate(bset%bao_z(bset%num_bao))
    allocate(bset%bao_obs(bset%num_bao))
    

    bao_invcov_file  = ReadIniFileName(Ini,'bao_invcov_file')
    call OpenTxtFile(bao_invcov_file, tmp_file_unit)
    do i=1,bset%num_bao
		read (tmp_file_unit,*, iostat=iopb) bset%bao_invcov(i,:)
	end do
	close(tmp_file_unit)
    bao_measurements_file = ReadIniFileName(Ini,'bao_measurements_file')
    call OpenTxtFile(bao_measurements_file, tmp_file_unit)
	do i=1,bset%num_bao
		read (tmp_file_unit,*, iostat=iopb) bset%bao_z(i),bset%bao_obs(i)
	end do
    close(tmp_file_unit) 

    if (iopb.ne.0) then
       stop 'Error reading mpk file'
    endif
 
   call Ini_Close_File(Ini)
   call ClearFileUnit(file_unit)

   baodatasets(num_bao_datasets) = bset

end subroutine ReadBaoDataset


!JD copied from Reid BAO code
function CMBToBAOrs(CMB)
   use settings
   use cmbtypes
   use ModelParams
   use Precision
   use ThermoData, only : z_drag
   implicit none
   Type(CMBParams) CMB
   real(dl) ::  adrag, atol, rsdrag
   real(dl), external :: dsoundda, rombint
   real(dl) :: CMBToBAOrs
   integer error
   
   adrag = 1.0d0/(1.0d0+z_drag)
   atol = 1e-6
   rsdrag = rombint(dsoundda,1d-8,adrag,atol)
   CMBToBAOrs = rsdrag

end function CMBToBAOrs

function D_v(CMB,z)
    Type(CMBParams) CMB
	real(dl), intent(IN) :: z
    real(dl)  D_v, Hz, ADD, hzoh0,omegam
    
    omegam = 1.d0 - CMB%omv - CMB%omk 
    hzoh0 = sqrt(CMB%omk*(1.d0+z)**2.d0+omegam*(1.d0+z)**3.d0 &
			+ CMB%omv*(1.d0+z)**(3.d0*(1.d0+CMB%w)))
    
    ADD = AngularDiameterDistance(z)*(1.d0+z)
	Hz = CMB%h0*1000.d0*hzoh0
    D_v = ((ADD)**2.d0*c*z/Hz)**(1.d0/3.d0)
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
	Acoustic = 100*D_v(CMB,z)*sqrt(omh2)/(ckm*z)
end function Acoustic

function rstodv(CMB,z)
	Type(CMBparams) CMB
	real(dl) rstodv
	real(dl), intent(IN)::z
	real(dl) rs
	real(dl), parameter :: rs_rescale = 154.6588d0/150.8192d0
	
	rs = CMBToBAOrs(CMB)*rs_rescale
	
	rstodv = rs/D_v(CMB,z)
end function rstodv


!===================================================================================

function BAO_LnLike(CMB)
   implicit none
   Type(CMBParams) CMB
   integer i,j,k
   real(dl)  BAO_LnLike
   real(dl) tot(num_bao_datasets)
   real(dl), allocatable :: BAO_theory(:)
   
   do i=1,num_bao_datasets
   		allocate(BAO_theory(baodatasets(i)%num_bao))
	
		if(baodatasets(i)%type_bao ==1)then
			do j=1, baodatasets(i)%num_bao
				BAO_theory(j) = rstodv(CMB,baodatasets(i)%bao_z(j))
			end do
		else if(baodatasets(i)%type_bao ==2)then
			do j=1, baodatasets(i)%num_bao
				BAO_theory(j) = Acoustic(CMB,baodatasets(i)%bao_z(j))
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
	
		if(feedback>1)write(*,*)'Bao dataset: '//trim(baodatasets(i)%name)//' LnLike = ',tot(i)
   		deallocate(BAO_theory)
   end do
    
   BAO_LnLike = sum(tot)
   if(feedback>1)write(*,*)'Bao_LnLike = ', Bao_LnLike

end function BAO_LnLike


end module bao