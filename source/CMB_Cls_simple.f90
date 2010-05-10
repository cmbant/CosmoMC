!Use CAMB
module CMB_Cls
  use cmbtypes
  use CAMB, only : CAMB_GetResults, CAMB_GetAge, CAMBParams, CAMB_SetDefParams,Transfer_GetMatterPower, &
          AccuracyBoost,  Cl_scalar, Cl_tensor, Cl_lensed, outNone, w_lam, &
          CAMBParams_Set, MT, CAMBdata, NonLinear_Pk, Reionization_GetOptDepth, CAMB_GetZreFromTau, &
          CAMB_GetTransfers,CAMB_FreeCAMBdata,CAMB_InitCAMBdata, CAMB_TransfersToPowers, &
          initial_adiabatic,initial_vector,initial_iso_baryon,initial_iso_neutrino, initial_iso_neutrino_vel
          
  use settings
  use snovae
  implicit none
  logical :: Use_SN =.false. !Compute Supernovae likelihoods only when background changes
  logical :: compute_tensors = .false.
  logical :: CMB_lensing = .false.
  logical :: use_nonlinear = .false.

  Type ParamSetInfo  
    
    Type (CosmoTheory) :: Theory
    Type (CAMBdata)    :: Transfers
    Type (CMBParams)   :: LastParams
  end Type ParamSetInfo

  integer, parameter :: ScalClOrder(3) = (/1,3,2/), TensClOrder(4) = (/1,4,2,3/)
      !Mapping of CAMB CL array ordering to TT , TE, EE, BB  
  integer :: ncalls = 0
  type(CAMBParams)  CAMBP 
  logical :: w_is_w  = .true.

contains
  subroutine CMBToCAMB(CMB,P)
    use LambdaGeneral
    type(CMBParams) CMB
    type(CAMBParams)  P
    P = CAMBP
    P%omegab = CMB%omb
    P%omegan = CMB%omnu
    P%omegac = CMB%omc
    P%omegav = CMB%omv
    P%H0 = CMB%H0
    P%Reion%redshift= CMB%zre    
    if (w_is_w) then
     w_lam = CMB%w
    else
     P%InitialConditionVector(initial_iso_baryon) = CMB%w
     w_lam = -1
    end if
  end subroutine CMBToCAMB

 function RecomputeTransfers (A, B)
   logical RecomputeTransfers
   type(CMBParams) A, B

   RecomputeTransfers =  .not. (A%omb == B%omb .and. A%omc == B%omc .and. A%omv == B%omv .and. &
             A%omnu == B%omnu .and. A%zre == B%zre .and. A%omk == B%omk .and. A%w == B%w)
              
 end function RecomputeTransfers


 subroutine GetCls(CMB,Info, Cls, error)
   use ModelParams, only : ThreadNum
   type(CMBParams) CMB
   integer error
   Type(ParamSetInfo) Info
   real Cls(lmax,1:num_Cls)
   type(CAMBParams)  P
   logical NewTransfers
   integer zix
    error = 0
    Newtransfers = .false.
  
    if (RecomputeTransfers(CMB, Info%LastParams))  then
     !Slow parameters have changed
         call CAMB_InitCAMBdata(Info%Transfers)
         call CMBToCAMB(CMB, P) 
      
         if (Feedback > 1) write (*,*) 'Calling CAMB'
         Threadnum =num_threads
  
         call CAMB_GetTransfers(P, Info%Transfers, error)
         NewTransfers = .true.
         Info%LastParams = CMB
         if (Use_SN) then
            Info%Theory%SN_Loglike = SN_LnLike(CMB)
         else
            Info%Theory%SN_Loglike = 0     
         end if  
         ncalls=ncalls+1
         if (mod(ncalls,100)==0 .and. logfile_unit/=0) write (logfile_unit,*) 'CAMB called ',ncalls, ' times'
    
         if (Feedback > 1) write (*,*) 'CAMB done'

    end if 
    
!    if ((error==0) .and. (Newtransfers .or. any(CMB%InitPower /= Info%LastParams%InitPower))) then
       !Use the initial power spectra to get the Cls and matter power spectrum
     if (error == 0) then
      !Always get everything again. Slight waste of time in general, but allows complete mixing of fast
      !parameters, and works with lensing

         call SetCAMBInitPower(Info%Transfers%Params,CMB,1)       

         call CAMB_TransfersToPowers(Info%Transfers)
            !this sets slow CAMB params correctly from value stored in Transfers
           
         call SetTheoryFromCAMB(Info%Theory)
    
         if (any(Info%Theory%cl(:,1) < 0 )) then
            error = 1
            !Kill initial power spectra that go negative
            return
         end if
   
         if (Use_LSS) then
            Info%Theory%sigma_8 = Info%Transfers%MTrans%sigma_8(matter_power_lnzsteps,1)
            if (num_matter_power /= 0) then
                do zix = 1,matter_power_lnzsteps
                 call Transfer_GetMatterPower(Info%Transfers%MTrans,& 
                   Info%Theory%matter_power(:,zix),matter_power_lnzsteps-zix+1,&
                    1,matter_power_minkh, matter_power_dlnkh,num_matter_power) 
                 end do             
             end if
         else
            Info%Theory%sigma_8 = 0
         end if

     end if
     if (error /= 0) return

     call ClsFromTheoryData(Info%Theory, CMB, Cls)
     
 end subroutine GetCls

 subroutine SetTheoryFromCAMB(Theory)
   Type(CosmoTheory) Theory
   real, parameter :: cons =  2.726e6**2*2*pi
   real nm
   integer l

   !The reason we store tensors separately is that can then importance sample re-computing scalars only,
   !using the stored tensor C_l

    do l = 2, lmax

       nm = cons/(l*(l+1))
       if (CMB_Lensing) then
            Theory%cl(l,1:num_clsS) =  nm*Cl_lensed(l,1, TensClOrder(1:num_clsS))
         else 
            Theory%cl(l,1:num_clsS) =  nm*Cl_scalar(l,1, scalClOrder(1:num_clsS))
       end if         
   
       if (compute_tensors .and. l<=lmax_tensor) then
            Theory%cl_tensor(l,1:num_cls) =  nm*Cl_tensor(l,1, TensClOrder(1:num_cls))
       end if
    end do 
    
    if (num_cls>num_clsS) Theory%cl(:,num_clsS+1:num_cls) = 0 


 end subroutine SetTheoryFromCAMB

 subroutine GetClsInfo(CMB, Theory, error, DoCls, DoPk)
   use ModelParams, only : ThreadNum
   type(CMBParams) CMB
   Type(CosmoTheory) Theory
   integer error
   logical, intent(in) :: DoCls, DoPk
   type(CAMBParams)  P
   logical MatterOnly
   integer zix
   error = 0
   Threadnum =num_threads
   call CMBToCAMB(CMB, P)
   P%OnlyTransfers = .false.
   call SetCAMBInitPower(P,CMB,1)   
    
   MatterOnly = .false.
   if (DoPk) then
      P%WantTransfer = .true.
      if (.not. DoCls) then
         MatterOnly = .true.         
         P%WantScalars = .false.
         P%WantTensors = .false.
      end if
   end if
   if (DoCls) then
      !Assume we just want Cls to higher l
      P%WantScalars = .true.
      !P%WantTensors = .false.
      !compute_tensors = .false.
      P%WantTensors = compute_tensors 
       
       if (.not. DoPk) then
        P%WantTransfer = .false.
       end if
   end if
   
   call CAMB_GetResults(P)
   error = 0 !using error optional parameter gives seg faults on SGI
   if (error==0) then
       
      if (DoCls) then
  
       Theory%cl_tensor(2:lmax_tensor,1:num_cls) = 0
       call SetTheoryFromCAMB(Theory)
      end if
  
      if (DoPk) then 
         Theory%sigma_8 = MT%sigma_8(matter_power_lnzsteps,1)
         if (num_matter_power /= 0) then
           do zix = 1,matter_power_lnzsteps
                 call Transfer_GetMatterPower(MT,& 
             Theory%matter_power(:,zix),matter_power_lnzsteps-zix+1,1,matter_power_minkh,&
                        matter_power_dlnkh,num_matter_power) 
           end do
         
         end if    
      end if
      Theory%Age = CAMB_GetAge(P)

   end if
 end subroutine GetClsInfo


 subroutine InitCAMB(CMB,error, DoReion)
   type(CMBParams), intent(in) :: CMB
   logical, optional, intent(in) :: DoReion
   logical WantReion
   type(CAMBParams)  P
   integer error
   
   if (present(DoReion)) then
    WantReion = DoReion
   else
    WantReion = .true.
   end if

   call CMBToCAMB(CMB, P)
   call CAMBParams_Set(P,error,WantReion)

 end subroutine InitCAMB

 function GetOpticalDepth(CMB)
   type(CMBParams) CMB
   real GetOpticalDepth
   type(CAMBParams)  P
   integer error
   
   call CMBToCAMB(CMB, P)
   call CAMBParams_Set(P,error)

   if (error/= 0) then
      GetOpticalDepth = -1 
   else
      GetOpticalDepth = Reionization_GetOptDepth(P%Reion, P%ReionHist) 
   end if
 end  function GetOpticalDepth

 function GetZreFromTau(CMB, tau)
   type(CMBParams) CMB
   real, intent(in) :: tau
   real GetZreFromTau
   type(CAMBParams)  P
   
   call CMBToCAMB(CMB, P)
   GetZreFromTau = CAMB_GetZreFromTau(P,dble(tau))
 
 end  function GetZreFromTau 
 
 function GetAge(CMB, Info)
   !Return <0 if error
   real GetAge
   type(CMBParams) CMB
   Type(ParamSetInfo) Info
   type(CAMBParams)  P
   call CMBToCAMB(CMB, P)
 
   Info%Theory%Age = CAMB_GetAge(P)
 
   GetAge = Info%Theory%Age
 end function GetAge

 subroutine InitCAMBParams(P)
   use lensing
   use ModelParams
   use Lya
   type(CAMBParams)  P 
   integer zix

        Threadnum =num_threads
        w_lam = -1
        call CAMB_SetDefParams(P)
 
        P%OutputNormalization = outNone
      
        P%WantScalars = .true.
        P%WantTensors = compute_tensors
        P%WantTransfer = Use_LSS

        P%Max_l=lmax
        P%Max_eta_k=lmax*2
        P%Max_l_tensor=lmax_tensor
        P%Max_eta_k_tensor=lmax_tensor*5./2
 
        P%Transfer%k_per_logint=0
 
        if (use_nonlinear) then
         P%NonLinear = NonLinear_Pk
         P%Transfer%kmax = 1.2
        else
         P%Transfer%kmax = 0.8
        end if
        if (Use_Lya) P%Transfer%kmax = lya_kmax
        P%Transfer%num_redshifts = matter_power_lnzsteps
        
        if (AccuracyLevel > 1) then
          if (USE_LSS) then
            P%Transfer%high_precision=.true.
            P%Transfer%kmax=P%Transfer%kmax + 0.2
           end if
          AccuracyBoost = AccuracyLevel
          lAccuracyBoost = AccuracyLevel
          lSampleBoost = AccuracyLevel
          P%AccurateReionization = .true.
        end if
                
        if (max_transfer_redshifts < matter_power_lnzsteps) then
          stop 'Need to manually set max_transfer_redshifts larger in CAMB''s modules.f90'
        end if
        if (matter_power_lnzsteps==1) then
          P%Transfer%redshifts(1) = 0
        else
         do zix=1, matter_power_lnzsteps
           !Linear spacing in log(z+1)
           
          P%Transfer%redshifts(zix) = exp( log(matter_power_maxz+1) * &
                real(matter_power_lnzsteps-zix)/(max(2,matter_power_lnzsteps)-1) )-1
               !put in max(2,) to stop compilers complaining of div by zero
         end do
        end if   
        
        P%Num_Nu_Massive = 3
        P%Num_Nu_Massless = 0.04
        P%InitPower%nn = 1
        P%AccuratePolarization = num_cls/=1 
        P%Reion%use_optical_depth = .false.
        P%OnlyTransfers = .true.

        if (CMB_Lensing) then
            P%DoLensing = .true.
            P%Max_l = lmax +250
            P%Max_eta_k = P%Max_l*2 
        end if
        
        lensing_includes_tensors = .false.

        P%Scalar_initial_condition = initial_vector
        P%InitialConditionVector = 0
        P%InitialConditionVector(initial_adiabatic) = -1


 end subroutine InitCAMBParams

 subroutine CMB_Initialize(Info)
   Type(ParamSetInfo) Info
   type(CAMBParams)  P 
        compute_tensors = Ini_Read_Logical('compute_tensors',.false.)
        if (num_cls==3 .and. compute_tensors) write (*,*) 'WARNING: computing tensors with num_cls=3 (BB=0)'
        CMB_lensing = Ini_Read_Logical('CMB_lensing',.false.)

        if (Feedback > 0 ) then
          write (*,*) 'Computing tensors:', compute_tensors
          write (*,*) 'Doing CMB lensing:',CMB_lensing
          write(*,'(" lmax           = ",1I4)') lmax       
          if (compute_tensors) write(*,'(" lmax_tensor    = ",1I4)') lmax_tensor
          write(*,'(" Number of C_ls = ",1I4)') num_cls
        end if

        if (CMB_lensing) num_clsS = num_cls   !Also scalar B in this case
     
        call InitCAMBParams(P)
  
        call CAMB_InitCAMBdata(Info%Transfers)
    
        P%WantTensors = compute_tensors
        Info%LastParams%omb = -1 !Make sure we calculate the CMB first time called
        CAMBP = P
        
 end subroutine CMB_Initialize


 subroutine AcceptReject(accpt, CurParams, Trial)
   logical, intent(in) :: accpt
   Type(ParamSetInfo) CurParams, Trial

   if (.not. associated(CurParams%Transfers%ClTransScal%Delta_p_l_k,&
                        Trial%Transfers%ClTransScal%Delta_p_l_k)) then
    !If they point to same memory don't need to free anything
    if (accpt) then
       call CAMB_FreeCAMBdata(CurParams%Transfers)
    else
       call CAMB_FreeCAMBdata(Trial%Transfers)
    end if

   end if

 end subroutine AcceptReject

end module CMB_Cls

