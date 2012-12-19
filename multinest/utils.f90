 module RandomNS
 integer :: rand_instNS = 0
 double precision, dimension(:), allocatable :: C, CD, CM, GSET
 double precision, dimension(:,:), allocatable :: U
 integer, dimension(:), allocatable :: I97, J97, ISET
 integer numNodes

 contains

  subroutine initRandomNS(n,i)
  implicit none
  integer n !no. of nodes
  integer, optional, intent(IN) :: i
  integer kl,ij,k
  character(len=10) :: fred
  real :: klr

  !sanity  check
  if(n<=0) then
  	write(*,*)'you have asked for ',n,'nodes'
      stop
  endif

  numNodes=n

  !memory allocation
  allocate(C(n),CD(n),CM(n),U(n,97),I97(n),J97(n),ISET(n),GSET(n))

  ISET=0

  do k=1,n
  	if(present(i)) then
		kl = 9373
		ij = (i+(k-1))*45
  	else
		call system_clock(count=ij)
      	ij = mod(ij + rand_instNS*100, 31328)+(k-1)*45
            call date_and_time(time=fred)
       	read (fred,'(e10.3)') klr
       	kl = mod(int(klr*1000), 30081)
      end if

      !write(*,'(" randomNS seeds:",1I6,",",1I6," rand_instNS:",1I4)")') ij,kl,rand_instNS
      call RMARINNS(ij,kl,k)
  enddo
  end subroutine initRandomNS


  subroutine killRandomNS
      deallocate(C,CD,CM,U,I97,J97,ISET,GSET)
  end subroutine killRandomNS


  double precision function GAUSSIAN1NS(idg)
    implicit none
    integer idg,id !node no.
    double precision urv
    double precision R, V1, V2, FAC

    id=idg+1
    if(ISET(id)==0) then
	R=2
	do while (R >=1.d0)
		urv=ranmarns(id-1)
		V1=2.d0*urv-1.d0
		urv=ranmarns(id-1)
		V2=2.d0*urv-1.d0
		R=V1**2+V2**2
	end do
	FAC=sqrt(-2.d0*log(R)/R)
	GSET(id)=V1*FAC
	gaussian1ns=V2*FAC
	ISET(id)=1
    else
	gaussian1ns=GSET(id)
	ISET(id)=0
    endif

  end function GAUSSIAN1NS



  subroutine RMARINNS(IJ,KL,id)
! This is the initialization routine for the randomNS number generator ranmarns()
! NOTE: The seed variables can have values between:    0 <= IJ <= 31328
!                                                      0 <= KL <= 30081
!The randomNS number sequences created by these two seeds are of sufficient
! length to complete an entire calculation with. For example, if sveral
! different groups are working on different parts of the same calculation,
! each group could be assigned its own IJ seed. This would leave each group
! with 30000 choices for the second seed. That is to say, this randomNS
! number generator can create 900 million different subsequences -- with
! each subsequence having a length of approximately 10^30.
!
! Use IJ = 1802 & KL = 9373 to test the randomNS number generator. The
! subroutine ranmarns should be used to generate 20000 randomNS numbers.
! Then display the next six randomNS numbers generated multiplied by 4096*4096
! If the randomNS number generator is working properly, the randomNS numbers
!    should be:
!           6533892.0  14220222.0  7275067.0
!           6172232.0  8354498.0   10633180.0

	if(IJ<0) IJ=31328+IJ
	if(IJ>31328) IJ=mod(IJ,31328)
	if(KL<0) KL=30081+KL
	if(KL>30081) KL=mod(KL,30081)

  	if(IJ<0  .or.  IJ>31328  .or. KL<0  .or.  KL>30081 ) then
		print '(A)', ' The first randomNS number seed must have a value  between 0 and 31328'
		print '(A)',' The second seed must have a value between 0 and   30081'
		stop
  	endif
  	I = mod(IJ/177, 177) + 2
  	J = mod(IJ    , 177) + 2
  	K = mod(KL/169, 178) + 1
  	L = mod(KL,     169)
  	do II = 1, 97
		S = 0.0
		T = 0.5
		do JJ = 1, 24
			M = mod(mod(I*J, 179)*K, 179)
            	I = J
            	J = K
            	K = M
            	L = mod(53*L+1, 169)
            	if (mod(L*M, 64)>=32) S = S + T
            	T = 0.5 * T
		enddo
         	U(id,II) = S
      enddo
      C(id) = 362436.0 / 16777216.0
      CD(id) = 7654321.0 / 16777216.0
      CM(id) = 16777213.0 /16777216.0
      I97(id) = 97
      J97(id) = 33

  end subroutine RMARINNS



  double precision function ranmarns(idg)
! This is the random number generator proposed by George Marsaglia in
! Florida State University Report: FSU-SCRI-87-50
! It was slightly modified by F. James to produce an array of pseudorandom
! numbers.

      id=idg+1

	UNI = U(id,I97(id)) - U(id,J97(id))
  	if(UNI<0.) UNI = UNI + 1.
  	U(id,I97(id)) = UNI
  	I97(id) = I97(id) - 1
  	if(I97(id)==0) I97(id) = 97
  	J97(id) = J97(id) - 1
  	if(J97(id)==0) J97(id) = 97
  	C(id) = C(id) - CD(id)
  	if(C(id)<0.) C(id) = C(id) + CM(id)
  	UNI = UNI - C(id)
  	if(UNI<0.) UNI = UNI + 1. ! bug?
  	ranmarns = UNI

  end function ranmarns


 end module RandomNS
