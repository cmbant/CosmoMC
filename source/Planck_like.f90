!Pseudo-Cl (or other C_l esimator) based likelihood approximation for cut sky with polarization
!Simple harmonic low-l likelihood
!Obviously this is not a realistic Planck likelihood code
!AL Feb 2008

module CMBLikes
 use settings
 use cmbtypes
 use IniFile
 use AMLutils
 use MatrixUtils !Compiled with MATRIX_SINGLE defined
 implicit none

  integer :: cl_E = 3, cl_B=4 ! stop compile time errors with num_cls=3

  Type TSqMatrix
    real, dimension(:,:), pointer :: M
  end Type TSqMatrix

  Type TLowlLike
    integer nmodes
    integer tmodes, almmodes
    integer polalmmodes, EBmodes
    integer lexact, lmax
    double precision, dimension(:), pointer :: ModeDataVector 
    double precision, dimension(:,:), pointer :: TheoryProj, ReProj, ImProj
    double precision highlScaleT, highlScaleE, highlScaleC, highlScaleB 
    double precision, dimension(:,:) , pointer :: NoiseCov, HighlCov
  end Type TLowlLike

  Type TCMBLikes
     logical highl_cl, lowl_exact !Which bits of likelihood to include
    
     integer nfields !number of fields
     logical use_field(3)
     integer field_index(3) !mapping from 1,2,3 to index actually used
     integer fields(3) !indices (1=T,2=E,3=B) of fields to use
     character field_order(3)
     integer ncl !calculated from above = nfields*(nfields+1)/2
     integer ncl_used !Number of C_l actually used in covariance matrix (others assumed zero)
     integer cl_use_index(6)
     integer cl_lmin, cl_lmax !The range of l to use psuedo-cl-based likelihood
     integer vecsize
     integer like_approx
     real, dimension(:,:), pointer :: ClHat, ClFiducial, ClNoise
     real, dimension(:,:), pointer :: inv_covariance
     Type(TSqMatrix) ,dimension(:), pointer :: sqrt_fiducial
     Type(TLowlLike) :: LowL
  end Type TCMBLikes
  
  character(LEN=3), parameter :: field_names = 'TEB'
  integer, parameter :: like_approx_diag=1, like_approx_fid_gaussian=2, like_approx_fullsky_exact=3
   ! like_approx_diag: new approximation from Hammimeche & Lewis
   ! like_approx_fid_gaussian: fixed covariance matrix, (X-Xhat)^TC^{-1}(X-Xhat)
   ! like_approx_exact: ignore all correlations, use exact full sky likelihood function
  
contains


  subroutine CMBLikes_ReadLowlFile(D,aname)
    Type(TCMBLikes) :: D
    character(LEN=*), intent(in) :: aname
    integer  filemodes, nmodes, file_unit,i,j
    double precision, allocatable :: coupling_row(:)   
    double precision, dimension(:), allocatable :: TModeData, EModeData, BModeData

  !Note this doesn't currently support chopping to requested fields - always uses TEB
  !Also uses binary fortran files rather than e.g. FITS or other endian-independent standard
        file_unit = new_file_unit()

        call OpenFile(aname,file_unit,'unformatted')
   
        read (file_unit) filemodes, D%LowL%tmodes
        if (filemodes /= (D%Lowl%lmax+1)**2) call MpiStop('lowl likelihood lmax mismatch')
        D%LowL%almmodes = (D%Lowl%lexact+1)**2
        allocate(TModeData(D%LowL%tmodes))
        read(file_unit) TModeData
        allocate(coupling_row(filemodes))
        allocate(D%Lowl%TheoryProj(D%Lowl%almmodes , D%LowL%tmodes))
        do i=1, D%LowL%tmodes
         read(file_unit) coupling_row
         D%Lowl%TheoryProj(1:D%Lowl%almmodes,i) = coupling_row(1:D%Lowl%almmodes)
        end do
        deallocate(coupling_row)

        Read(file_unit) D%LowL%highlScaleT
        nmodes = D%LowL%tmodes
  
       
        read(file_unit) D%Lowl%highlScaleE, D%Lowl%highlScaleC, D%Lowl%highlScaleB
        
        read(file_unit) filemodes, D%LowL%EBmodes
        D%Lowl%polalmmodes  =  (D%Lowl%lexact+1)**2-4

        allocate(EModeData(D%LowL%EBmodes))
        allocate(BModeData(D%LowL%EBmodes))
        read(file_unit) EModeData, BModeData
 
        allocate(D%Lowl%ReProj(D%Lowl%polalmmodes, D%LowL%EBmodes))
        allocate(D%Lowl%ImProj(D%Lowl%polalmmodes, D%LowL%EBmodes))
        allocate(coupling_row(filemodes))
        do i=1, D%LowL%EBmodes
          read(file_unit) coupling_row
          D%Lowl%ReProj(1:D%Lowl%polalmmodes,i) = coupling_row(1:D%Lowl%polalmmodes)
          read(file_unit) coupling_row
          D%Lowl%ImProj(1:D%Lowl%polalmmodes,i) = coupling_row(1:D%Lowl%polalmmodes)
        end do
        deallocate(coupling_row)
        nmodes= nmodes + D%LowL%EBmodes*2

       
        allocate(D%LowL%NoiseCov(nmodes,nmodes))
        allocate(D%LowL%HighlCov(nmodes,nmodes))
        do i=1,nmodes
          read(file_unit) D%LowL%NoiseCov(1:i,i)
          read(file_unit) D%LowL%HighlCov(1:i,i)          
        end do
        do i=1,nmodes
         do j=i+1, nmodes
           D%LowL%NoiseCov(j,i) =  D%LowL%NoiseCov(i,j) 
           D%LowL%HighlCov(j,i) =  D%LowL%HighlCov(i,j) 
         end do
        end do 
          
        read(file_unit) i
        if (i/=252353) call MpiStop('Bad low l likelihood data file') 

        call CloseFile(file_unit)

        allocate( D%LowL%ModeDataVector(nmodes))
        D%LowL%ModeDataVector(1:D%LowL%tmodes) = TModeData
        D%LowL%ModeDataVector(D%LowL%tmodes+1:D%LowL%tmodes+D%LowL%EBmodes) = EModeData
        D%LowL%ModeDataVector(D%LowL%tmodes+D%LowL%EBmodes+1:D%LowL%tmodes+2*D%LowL%EBmodes) = BModeData
        deallocate(TModeData,BModeData,EModeData)

        D%LowL%nmodes = nmodes

  end subroutine CMBLikes_ReadLowlFile

 function idx_T(l,m)
  !Conversion from l,m into vectors
  integer, intent(in) :: l,m
  integer idx_T

  idx_T = l*(l+1) + m +1 

 end function

 function idx_P(l,m)
  !Conversion from l,m into vectors
  integer, intent(in) :: l,m
  integer idx_P

  idx_P = l*(l+1) + m +1 -4 

 end function


  subroutine CMBLikes_lowl_GetFullCovariance(D, Cov, cl, lmin, lsum)
   Type(TCMBLikes) :: D
   double precision Cov(:,:)
   real cl(lmax,num_cls)
   integer, intent (in) :: lmin, lsum
   integer i,j
   double precision :: sum1,sum2,tmp,tmpEE,tmpBB, tmpEB, tmpBE
   integer l, mix, mixT

   Cov=0
   
   do i=1, D%Lowl%tmodes
    do j=1, i 
       tmp = 0
       do l=lmin,lsum
        mix = idx_T(l,-l)
        tmp =tmp + dot_product(D%Lowl%TheoryProj(mix:mix+2*l,i),D%Lowl%TheoryProj(mix:mix+2*l,j))*cl(l,1)
       end do      
       Cov(i,j) = tmp
       if (i/=j) Cov(j,i) =  tmp
    end do
   end do 
  
    do i=1, D%Lowl%tmodes
        do j=1, D%Lowl%EBmodes

        !TE
        tmp = 0
        do l=lmin,lsum
            mix = idx_P(l,-l)
            mixT = idx_T(l,-l)
            tmp =tmp + dot_product(D%Lowl%TheoryProj(mixT:mixT+2*l,i),&
                D%Lowl%ReProj(mix:mix+2*l,j))*cl(l,2)
        end do      
        Cov(i,D%Lowl%tmodes+j) = tmp
        Cov(D%Lowl%tmodes+j,i) =  tmp
        
        !TB
        tmp = 0
        do l=lmin,lsum
            mix = idx_P(l,-l)
            mixT = idx_T(l,-l)
            tmp =tmp - dot_product(D%Lowl%TheoryProj(mixT:mixT+2*l,i),&
                D%Lowl%ImProj(mix:mix+2*l,j))*cl(l,2) 
   
        end do      
        Cov(i,D%Lowl%tmodes+D%Lowl%EBmodes+j) = tmp
        Cov(D%Lowl%tmodes+D%Lowl%EBmodes+j,i) =  tmp
        
        
       end do
      end do 


     do i=1, D%Lowl%EBmodes
        do j=1, i 
        
        !EE   
        tmpEE = 0
        tmpBB = 0
        do l=lmin,lsum
            mix = idx_P(l,-l)
            sum1=dot_product(D%Lowl%ReProj(mix:mix+2*l,i),D%Lowl%ReProj(mix:mix+2*l,j))
            sum2=dot_product(D%Lowl%ImProj(mix:mix+2*l,i),D%Lowl%ImProj(mix:mix+2*l,j))
            tmpEE =tmpEE + sum1*cl(l,cl_E)  + sum2*cl(l,cl_B)
            tmpBB =tmpBB + sum1*cl(l,cl_B)  + sum2*cl(l,cl_E)
        end do      
        Cov(D%Lowl%tmodes+i,D%Lowl%tmodes+j) = tmpEE
        if (i/=j) Cov(D%Lowl%tmodes+j,D%Lowl%tmodes+i) =  tmpEE
        Cov(D%Lowl%tmodes+D%Lowl%EBmodes+i,D%Lowl%tmodes+D%Lowl%EBmodes+j) = tmpBB
        if (i/=j) Cov(D%Lowl%tmodes+D%Lowl%EBmodes+j,D%Lowl%tmodes+D%Lowl%EBmodes+i) =  tmpBB
   
        !EB/BE
        tmpEB = 0
        tmpBE=0
        do l=lmin,lsum
            mix = idx_P(l,-l)
            sum1=dot_product(D%Lowl%ReProj(mix:mix+2*l,i),D%Lowl%ImProj(mix:mix+2*l,j))
            sum2=dot_product(D%Lowl%ImProj(mix:mix+2*l,i),D%Lowl%ReProj(mix:mix+2*l,j))
            tmpEB =tmpEB - sum1*cl(l,cl_E) + sum2*cl(l,cl_B)
            tmpBE =tmpBE + sum1*cl(l,cl_B) - sum2*cl(l,cl_E) 
        end do      
        Cov(D%Lowl%tmodes+i,D%Lowl%tmodes+D%Lowl%EBmodes+j) = tmpEB
        Cov(D%Lowl%tmodes+D%Lowl%EBmodes+j,D%Lowl%tmodes+i) =  tmpEB
        Cov(D%Lowl%tmodes+D%Lowl%EBmodes+i,D%Lowl%tmodes+j) = tmpBE
        Cov(D%Lowl%tmodes+j,D%Lowl%tmodes+D%Lowl%EBmodes+i) =  tmpBE
         
        end do
    end do 

    
  end subroutine CMBLikes_lowl_GetFullCovariance


  function CMBLikes_lowl_CMBLike(D, cl) result (chisq)
    Type(TCMBLikes) :: D
    real cl(lmax,num_cls)
    real chisq
    double precision, allocatable :: Cov(:,:)
    integer j
   

    print *,'getting low l'
   
    allocate(Cov(D%Lowl%nmodes,D%Lowl%nmodes))
 
    call CMBLikes_lowl_GetFullCovariance(D, Cov, cl, 2, D%Lowl%lexact)

        do j=1,D%Lowl%nmodes
         Cov(:,j) = Cov(:,j)+D%Lowl%NoiseCov(:,j) 
        end do

   !Scale high l
       Cov(1:D%Lowl%tmodes,1:D%Lowl%tmodes) = &
        Cov(1:D%Lowl%tmodes,1:D%Lowl%tmodes)  &
          + cl(D%Lowl%lexact+1,1)/D%Lowl%highlScaleT*D%Lowl%HighlCov(1:D%Lowl%tmodes,1:D%Lowl%tmodes)

     
        !TE
        Cov(1:D%Lowl%tmodes,D%Lowl%tmodes+1:D%Lowl%tmodes+D%Lowl%EBmodes) = & 
          Cov(1:D%Lowl%tmodes,D%Lowl%tmodes+1:D%Lowl%tmodes+D%Lowl%EBmodes) + & 
           sqrt(cl(D%Lowl%lexact+1,1)/D%Lowl%highlScaleT*cl(D%Lowl%lexact+1,3)/D%Lowl%highlScaleE) * &
          D%Lowl%HighlCov(1:D%Lowl%tmodes,D%Lowl%tmodes+1:D%Lowl%tmodes+D%Lowl%EBmodes)
       
        Cov(D%Lowl%tmodes+1:D%Lowl%tmodes+D%Lowl%EBmodes,1:D%Lowl%tmodes) = &
         transpose(Cov(1:D%Lowl%tmodes,D%Lowl%tmodes+1:D%Lowl%tmodes+D%Lowl%EBmodes)) 
 
       
       !EE
        Cov(D%Lowl%tmodes+1:D%Lowl%tmodes+D%Lowl%EBmodes,D%Lowl%tmodes+1:D%Lowl%tmodes+D%Lowl%EBmodes) = & 
         Cov(D%Lowl%tmodes+1:D%Lowl%tmodes+D%Lowl%EBmodes,D%Lowl%tmodes+1:D%Lowl%tmodes+D%Lowl%EBmodes) + & 
           cl(D%Lowl%lexact+1,3)/D%Lowl%highlScaleE* &
          D%Lowl%HighlCov(D%Lowl%tmodes+1:D%Lowl%tmodes+D%Lowl%EBmodes, &
           D%Lowl%tmodes+1:D%Lowl%tmodes+D%Lowl%EBmodes)
  

       if (num_cls > 3) then
       !BB
        Cov(D%Lowl%tmodes+D%Lowl%EBmodes+1:D%Lowl%tmodes+2*D%Lowl%EBmodes,&
                       D%Lowl%tmodes+D%Lowl%EBmodes+1:D%Lowl%tmodes+2*D%Lowl%EBmodes) = &
              Cov(D%Lowl%tmodes+D%Lowl%EBmodes+1:D%Lowl%tmodes+2*D%Lowl%EBmodes, &
                       D%Lowl%tmodes+D%Lowl%EBmodes+1:D%Lowl%tmodes+2*D%Lowl%EBmodes) + &
              cl(D%Lowl%lexact+1,cl_B)/D%Lowl%highlScaleB* &
               D%Lowl%HighlCov(D%Lowl%tmodes+D%Lowl%EBmodes+1:D%Lowl%tmodes+2*D%Lowl%EBmodes, &
              D%Lowl%tmodes+D%Lowl%EBmodes+1:D%Lowl%tmodes+2*D%Lowl%EBmodes)              

       
        Cov(D%Lowl%tmodes+D%Lowl%EBmodes+1:D%Lowl%tmodes+2*D%Lowl%EBmodes,&
              D%Lowl%tmodes+1:D%Lowl%tmodes+D%Lowl%EBmodes) = &
         Cov(D%Lowl%tmodes+D%Lowl%EBmodes+1:D%Lowl%tmodes+2*D%Lowl%EBmodes,&
              D%Lowl%tmodes+1:D%Lowl%tmodes+D%Lowl%EBmodes) + & 
              cl(D%Lowl%lexact+1,cl_E)/D%Lowl%highlScaleE* &
               D%Lowl%HighlCov(D%Lowl%tmodes+D%Lowl%EBmodes+1:D%Lowl%tmodes+2*D%Lowl%EBmodes,&
              D%Lowl%tmodes+1:D%Lowl%tmodes+D%Lowl%EBmodes)
       !EB    
        Cov(D%Lowl%tmodes+1:D%Lowl%tmodes+D%Lowl%EBmodes, &
                 D%Lowl%tmodes+D%Lowl%EBmodes+1:D%Lowl%tmodes+2*D%Lowl%EBmodes) = &
           transpose(Cov(D%Lowl%tmodes+D%Lowl%EBmodes+1:D%Lowl%tmodes+2*D%Lowl%EBmodes,&
              D%Lowl%tmodes+1:D%Lowl%tmodes+D%Lowl%EBmodes))

        end if

    chisq = 2*Matrix_GaussianLogLikeDouble(Cov,D%LowL%ModeDataVector)
   
    if (Feedback > 1) print *,'lowl chisq = ', chisq
    deallocate(Cov)

  end function CMBLikes_lowl_CMBLike

 function TypeIndex(C)
  character(LEN=*), intent(in) :: C
  integer TypeIndex
  !Get order T, E, B -> 1,2,3
  if (C=='T') then
   TypeIndex=1
  elseif (C=='E') then
   TypeIndex=2
  elseif (C=='B') then
   TypeIndex=3
  else
   call mpiStop('Invalid C_l part, must be T E or B')
  end if

 end function TypeIndex

 subroutine CMBLikes_ReadClArr(D, aname, order, Cl, lmin)
  Type(TCMBLikes) :: D
  character(LEN=*), intent(in) :: aname, order
  integer, intent(in) :: lmin
  real :: Cl(:,lmin:)
  character(LEN=1024) :: tmp
  integer ix, i, j,i1,l,ll
  integer cols(6)
  real norm,tmp_ar(6)
  Type (TStringList) :: Li
  integer file_unit
      
     call TStringList_Init(Li)
     call TStringList_SetFromString(Li,order,field_names)
     ix=0
     cols=0
     do i=1,D%nfields
      do j=1,i
        ix = ix +1
        i1 = TStringList_IndexOf(Li, D%field_order(i)//D%field_order(j))
        if (i1==-1) i1 = TStringList_IndexOf(Li, D%field_order(j)//D%field_order(i)) 
        if (i1/=-1) then
           cols(ix) = i1
        end if  
      end do       
     end do

     file_unit = new_file_unit()
     call OpenTxtFile(aname, file_unit)
     Cl=0
     do
      read(file_unit,'(a)',end=1) tmp
      read(tmp,*, end=1) l, tmp_ar(1:Li%Count)
      ll=l
      if (l>=D%cl_lmin .and. l <=D%cl_lmax) then
        norm = l*(l+1)/twopi
        do ix=1,D%ncl
          if (cols(ix)/=0) Cl(ix,l) = tmp_ar(cols(ix))/norm
        end do
      end if
     end do 
     if (ll<D%cl_lmax) then
       write(*,*) 'CMBLikes_ReadClArr: C_l file does not go up to lmax:', D%cl_lmax
       write (*,*) trim(aname)
       call MpiStop()
     end if
1    call CloseFile(file_unit)

     call TStringList_Clear(Li)

 end subroutine CMBLikes_ReadClArr


  subroutine UseString_to_colIx(D, S, C, totnum)
    Type(TCMBLikes) :: D
   character(LEN=*), intent(in) :: S
   integer, intent(inout), optional :: totnum
   integer :: C(:,:)
   Type(TStringList) :: L
   integer i,i1,i2

    call TStringList_Init(L)
    call TStringList_SetFromString(L, S, 'TEB')
    C=0

     do i=1, L%Count
       if (size(L%Items(i)%P)/=2) call mpiStop('Invalid C_l order')
       i1= D%field_index(TypeIndex(L%Items(i)%P(1)))
       i2= D%field_index(TypeIndex(L%Items(i)%P(2)))
       if (i1/=0 .and. i2/=0) then
        C(i1,i2) = i
        C(i2,i1)= i
       end if
     end do       

    if (present(totnum)) totnum = L%Count
   call TStringList_Clear(L)

  end subroutine UseString_to_colIx


 subroutine CMBLikes_ReadDataFile(D, aname)
  Type(TCMBLikes) :: D
   character(LEN=*) :: aname
   
   Type(TIniFile) :: Ini 
   logical bad
 
   call Ini_Open_File(Ini,aname, 1, bad, .false.)
   if (bad) then
     call MpiStop('Error opening dataset file '//trim(aname))     
   end if

   call CMBLikes_ReadData(D, Ini, ExtractFilePath(aname))

   call Ini_Close_File(Ini)
 
 
 end  subroutine CMBLikes_ReadDataFile 

 subroutine CMBLikes_ReadData(D, Ini,dataset_dir)
  Type(TCMBLikes) :: D
  Type(TIniFile) :: Ini 
  real, dimension(:,:), allocatable :: Cov
  character(LEN=*), intent(in) :: dataset_dir
  integer ix, i
  character(LEN=Ini_max_string_len) :: S, S_order
  integer l,j
  integer, dimension(:,:), allocatable :: indices
  integer lmin_covmat,lmax_covmat, vecsize_in
  integer cov_num_cls
!  character(LEN=Ini_max_string_len) cache_name
   Ini_fail_on_not_found = .true.

   S = Ini_Read_String_File(Ini,'fields_use')
   D%use_field = .false.
   do i=1, len_trim(S)
     if (trim(S(i:i))/='') D%use_field(TypeIndex(S(i:i))) = .true. 
   end do

   D%highl_cl = Ini_Read_Logical_File(Ini,'highl_cl')

   D%nfields=0
   D%vecsize =0
   D%ncl=0
   D%ncl_used=0 
  
   if (D%highl_cl) then

   D%like_approx = Ini_read_Int_File(Ini,'like_approx')

   D%nfields = count(D%use_field)
   ix=0
   D%field_index=0
   do i=1,3
    if (D%use_field(i)) then
        ix=ix+1
        D%field_index(i)=ix
        D%fields(ix) =i
        D%field_order(ix) = field_names(i:i)
    end if
   end do 
   D%ncl = (D%nfields*(D%nfields+1))/2


   D%cl_lmin = Ini_Read_Int_File(Ini,'cl_lmin')
   D%cl_lmax = Ini_Read_Int_File(Ini,'cl_lmax')
   D%vecsize = (D%cl_lmax-D%cl_lmin+1)
   D%ncl_used = 0
   
   allocate(D%ClHat(D%ncl,D%cl_lmin:D%cl_lmax))
   allocate(D%ClFiducial(D%ncl,D%cl_lmin:D%cl_lmax))
   allocate(D%ClNoise(D%ncl,D%cl_lmin:D%cl_lmax))
   
   S = Ini_Read_String_File(Ini,'cl_hat_file')
   call StringReplace('%DATASETDIR%',dataset_dir,S)
   S_order = Ini_read_String_File(Ini,'cl_hat_order')
   call CMBLikes_ReadClArr(D, S,S_order,D%ClHat,D%cl_lmin)

   S = Ini_Read_String_File(Ini,'cl_fiducial_file')
   call StringReplace('%DATASETDIR%',dataset_dir,S)
   S_order = Ini_read_String_File(Ini,'cl_fiducial_order')
   call CMBLikes_ReadClArr(D, S,S_order,D%ClFiducial,D%cl_lmin)

   S = Ini_Read_String_File(Ini,'cl_noise_file')
   call StringReplace('%DATASETDIR%',dataset_dir,S)
   S_order = Ini_read_String_File(Ini,'cl_noise_order')
   call CMBLikes_ReadClArr(D, S,S_order,D%ClNoise,D%cl_lmin)

   if (.not. Ini_Read_Logical_File(Ini,'cl_hat_includes_noise')) then
     D%ClHat =  D%ClHat + D%ClNoise
   end if

   if (D%like_approx /= like_approx_fullsky_exact) then
        
    if (D%like_approx==like_approx_diag) then
        allocate(D%sqrt_fiducial(D%cl_lmin:D%cl_lmax))
        do l=D%cl_lmin,D%cl_lmax
         allocate(D%sqrt_fiducial(l)%M(D%nfields,D%nfields))
         call ElementsToMatrix(D, D%ClFiducial(:,l)+D%ClNoise(:,l), D%sqrt_fiducial(l)%M)
         call Matrix_Root(D%sqrt_fiducial(l)%M, D%nfields, 0.5) 
        end do
    end if

    lmax_covmat = Ini_Read_Int_File(Ini,'covmat_lmax')
    lmin_covmat = Ini_Read_Int_File(Ini,'covmat_lmin')
    if (lmin_covmat > D%cl_lmin) call MpiStop('lmin_covmat must be  <= cl_lmin')
    if (lmax_covmat < D%cl_lmax) call MpiStop('lmax_covmat must be  >= cl_lmax')
        
    S = Ini_Read_String_File(Ini,'covmat_cl')

    allocate(indices(D%nfields,D%nfields))
    call UseString_to_colIx(D, S, indices, cov_num_cls)
    call MatrixToElementsInt(D,indices,D%cl_use_index)
    deallocate(indices)
        
    D%ncl_used = count(D%cl_use_index(1:D%ncl) /=0)
            
    vecsize_in =  (lmax_covmat-lmin_covmat+1)

        
    allocate(D%inv_covariance(D%vecsize*D%ncl_used,  D%vecsize*D%ncl_used))
        
    S = Ini_Read_String_File(Ini,'covmat_fiducial')
    call StringReplace('%DATASETDIR%',dataset_dir,S)
!    cache_name=trim(S)//'.inverse_'//trim(IntToStr(D%cl_lmax))//'-'//trim(IntToStr(D%cl_lmin))
!    j=0
!    do i=1, D%ncl
!        if (D%cl_use_index(i)) j=j+2**i
!    end do
!    cache_name=trim(cache_name)//'_'//trim(IntToStr(j))

!    if (FileExists(cache_name)) then
!        call MatrixSym_Read_Binary(cache_name, D%inv_covariance)
!    else
        allocate(Cov(vecsize_in*cov_num_cls,vecsize_in*cov_num_cls)) 
        call MatrixSym_Read_Binary(S, Cov)

        do i=1, D%ncl_used
            do j=1,D%ncl_used
            
            D%inv_covariance((i-1)*D%vecsize+1:i*D%vecsize,(j-1)*D%vecsize+1:j*D%vecsize) &
            = Cov((i-1)*vecsize_in+(D%cl_lmin-lmin_covmat+1):(i-1)*vecsize_in +(D%cl_lmax-lmin_covmat+1), &
                    (j-1)*vecsize_in+(D%cl_lmin-lmin_covmat+1):(j-1)*vecsize_in +(D%cl_lmax-lmin_covmat+1))
            end do
        end do
        deallocate(Cov)
        call Matrix_inverse(D%inv_covariance)
!        if (IsMainMPI()) then
!         call MatrixSym_Write_Binary(cache_name, D%inv_covariance)
!        end if
!    end if

    end if 

    end if !want high l like

    D%lowl_exact = Ini_Read_Logical_File(Ini,'lowl_exact')
    if (D%lowl_exact) then
     if (num_cls==3) call MpiStop('CMBLikes current untested for only 3 C_l')
     S = Ini_Read_String_File(Ini,'lowl_datafile')
     call StringReplace('%DATASETDIR%',dataset_dir,S)
     D%Lowl%lexact = Ini_Read_Int_File(Ini,'lowl_lexact')
     D%Lowl%lmax = Ini_Read_Int_File(Ini,'lowl_lmax')
     call CMBLikes_ReadLowlFile(D,S)
    end if
  
  
 end subroutine CMBLikes_ReadData



 subroutine CMBLikes_Transform(D, C, Chat, CsHalf)
  !Get  C = C_s^{1/2}  U f(D) U^T C_s^{1/2} where C^{-1/2} CHat C^{-1/2} = U D U^T
  Type(TCMBLikes) :: D
  real C(D%nfields,D%nfields), CHat(D%nfields,D%nfields), CsHalf(D%nfields,D%nfields)
  real :: U(D%nfields,D%nfields), Rot(D%nfields,D%nfields) 
  real :: roots(D%nfields)
  real :: diag(D%nfields)
  integer i

        U = C
        call Matrix_Diagonalize(U,Diag,D%nfields)
        
        Rot= matmul(matmul(transpose(U),CHat),U)
     
        roots = sqrt(Diag)
        
        do i=1, D%nfields
         Rot(i,:)=Rot(i,:)/roots(i)
         Rot(:,i)=Rot(:,i)/roots(i)
        end do
 
        Rot = matmul(U,matmul(Rot,transpose(U))) 
        call Matrix_Diagonalize(Rot,Diag,D%nfields)

        Diag = sign(sqrt(2*max(0.,Diag-log(Diag)-1)),Diag-1)
        !want f(D)-1 to save calculating X-X_s
       
        U = matmul(CsHalf,Rot)
        C = U
        do i=1, D%nfields
          C(:,i) = C(:,i)*Diag(i)
        end do
        C = MatMul(C,transpose(U))


 end subroutine CMBLikes_Transform


 subroutine MatrixToElements(D, M, X)
  Type(TCMBLikes) :: D
  real :: M(D%nfields,D%nfields) 
  real :: X(D%ncl)
  integer ix,i,j

  ix=0
  do i=1, D%nfields
   do j=1,i
     ix = ix+1
     X(ix) = M(i,j)
   end do
  end do     

 end subroutine MatrixToElements

 subroutine MatrixToElementsInt(D, M, X)
  Type(TCMBLikes) :: D
  integer, intent(in) :: M(D%nfields,D%nfields) 
  integer,intent(out) :: X(D%ncl)
  integer ix,i,j

  ix=0
  do i=1, D%nfields
   do j=1,i
     ix = ix+1
     X(ix) = M(i,j)
   end do
  end do     

 end subroutine MatrixToElementsInt


 subroutine ElementsToMatrix(D, X, M)
  Type(TCMBLikes) :: D
  real, intent(out) :: M(D%nfields,D%nfields) 
  real, intent(in) :: X(D%ncl)
  integer ix,i,j

  ix=0
  do i=1, D%nfields
   do j=1,i
     ix = ix+1
     M(i,j) = X(ix)
     M(j,i) = M(i,j)
   end do
  end do     

 end subroutine ElementsToMatrix

 function ExactChiSq(D, C,Chat,l)
  Type(TCMBLikes) :: D
  real, intent(in) :: C(D%nfields,D%nfields), Chat(D%nfields,D%nfields)
  integer, intent(in) :: l   
  real ExactChiSq
  real M(D%nfields,D%nfields)
  
  M = C
  call Matrix_root(M,D%nfields,-0.5)
  M = matmul(M,matmul(Chat,M))
  ExactChiSq = (2*l+1)*(Matrix_Trace(M) - D%nfields - MatrixSym_LogDet(M) )
  
 end function ExactChiSq



 function CMBLikes_CMBLike(D, cl) result (chisq)
  Type(TCMBLikes) :: D
  real cl(lmax,num_cls)
  real chisq
  real C(D%nfields,D%nfields), Chat(D%nfields,D%nfields)
  real vecp(D%ncl)
  real bigX(D%vecsize*D%ncl_used)
  integer l,  i, Ti,Ei,Bi
  logical :: quadratic

  chisq =0

  if (D%highl_cl) then

  Ti = D%field_index(1)
  Ei = D%field_index(2)
  Bi = D%field_index(3)
  
  if (Bi/=0 .and. num_cls<3) call MpiStop('CMBLikes_CMBLike: Need num_cls =4 to use B modes')
   
     
  do l = D%cl_lmin, D%cl_lmax 
   
       call ElementsToMatrix(D,D%ClNoise(:,l), C) 
       if (Ti/=0) C(Ti,Ti) = C(Ti,Ti)+cl(l,1)
       if (Ei/=0) then
         C(Ei,Ei) = C(Ei,Ei)+cl(l,3)
         if (Ti/=0) then
            C(Ei,Ti) = C(Ei,Ti) + cl(l,2)
            C(Ti,Ei) = C(Ei,Ti)
         end if
       end if
       if (Bi/=0) C(Bi,Bi) = C(Bi,Bi) + cl(l,num_cls)
 
             
       call ElementsToMatrix(D,D%ClHat(:,l), Chat) !could cache this
        
       if (D%like_approx == like_approx_diag) then 
        
        call CMBLikes_Transform(D, C, Chat, D%sqrt_fiducial(l)%M)
        call MatrixToElements(D, C, vecp)
        quadratic = .true.
       else if (D%like_approx == like_approx_fid_gaussian) then 
         C = C - Chat
         call MatrixToElements(D, C, vecp)
         quadratic = .true.
       else if (D%like_approx == like_approx_fullsky_exact) then
          quadratic = .false.
          chisq = chisq  + ExactChisq(D, C,Chat,l)
       else
        call MpiStop('Unknown like_approx') 
       end if  

        if (quadratic) then
        do i=1,D%ncl
         if (D%cl_use_index(i)/=0) then
            bigX( (D%cl_use_index(i)-1)*D%vecsize + l-D%cl_lmin+1) = vecp(i)
         end if
        end do
        end if
        

  end do

   if (quadratic) chisq = Matrix_QuadForm(D%inv_covariance,BigX)

   end if

   if (D%lowl_exact) then
       chisq = chisq + CMBLikes_lowl_CMBLike(D, cl) 
   end if
   
 end function CMBLikes_CMBLike


end module CMBLikes
