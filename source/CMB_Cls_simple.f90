!Use CAMB
module CMB_Cls
  use cmbtypes
  use CAMB, only : CAMB_GetResults, CAMB_GetAge, CAMBParams, CAMB_SetDefParams,Transfer_GetMatterPower, &
          AccuracyBoost,  Cl_scalar, Cl_tensor, Cl_lensed, outNone, w_lam, wa_ppf,&
          CAMBParams_Set, MT, CAMBdata, NonLinear_Pk, Reionization_GetOptDepth, CAMB_GetZreFromTau, &
          CAMB_GetTransfers,CAMB_FreeCAMBdata,CAMB_InitCAMBdata, CAMB_TransfersToPowers, &
          initial_adiabatic,initial_vector,initial_iso_baryon,initial_iso_CDM, initial_iso_neutrino, initial_iso_neutrino_vel, &
          HighAccuracyDefault, highL_unlensed_cl_template, ThermoDerivedParams, nthermo_derived
  use Errors !CAMB        
  use settings
  use snovae
  use bao
  use HST
  use IO
  implicit none
  logical :: Use_SN =.false. !Compute Supernovae likelihoods only when background changes
  logical :: Use_HST =.false. !Compute HST likelihoods only when background changes
  logical :: Use_BAO = .false.

  logical :: compute_tensors = .false.
  logical :: CMB_lensing = .false.
  logical :: use_nonlinear = .false.
  
  Type ParamSetInfo  
    
    Type (CosmoTheory) :: Theory
    Type (CAMBdata)    :: Transfers
    Type (CMBParams)   :: LastParams
  end Type ParamSetInfo

  integer :: lmax_computed_cl = lmax !value used for CAMB
  
  integer, parameter :: ScalClOrder(5) = (/1,3,2,4,5/), TensClOrder(4) = (/1,4,2,3/)
      !Mapping of CAMB CL array ordering to TT , TE, EE, BB, phi, phiT  
  integer :: ncalls = 0
  integer :: nerrors = 0
  type(CAMBParams)  CAMBP 
  
  real, allocatable :: highL_lensedCL_template(:,:)
  
contains
  subroutine CMBToCAMB(CMB,P)
    use LambdaGeneral
    use CAMBmain, only : ALens 
    type(CMBParams) CMB
    type(CAMBParams)  P
    P = CAMBP
    P%omegab = CMB%omb
    P%omegan = CMB%omnu
    P%omegac = CMB%omc
    P%omegav = CMB%omv
    P%H0 = CMB%H0
    P%Reion%redshift= CMB%zre    
    P%Reion%delta_redshift = CMB%zre_delta
    w_lam = CMB%w
    wa_ppf = CMB%wa
    ALens = CMB%ALens
    P%InitialConditionVector(initial_iso_CDM) = sign(sqrt(abs(CMB%iso_cdm_correlated) &
           /(1-abs(CMB%iso_cdm_correlated))),CMB%iso_cdm_correlated)
    
!    if (CMB%nnu < 3.04) call MpiStop('CMBToCAMB: nnu < 3.04, would give negative masless neutrinos')
    !Not clear this recipe is the best thing to do in general, but should work for massless case with unusual nnu
    if (num_massive_neutrinos == -1) then
     P%Num_Nu_Massive = int(CMB%nnu)
    else
     P%Num_Nu_Massive = num_massive_neutrinos 
    end if

    P%Num_Nu_Massless = CMB%nnu - P%Num_Nu_Massive !AL Sept 11 for CAMB's new treatment; previously 3.046; we assume three massive
    P%YHe = CMB%YHe
#ifdef COSMOREC    
    if (P%Recomb%fdm/=0.) P%Recomb%runmode = 3
    P%Recomb%fdm = CMB%fdm * 1e-23
#else
    if (CMB%fdm/=0.) call MpiStop('Compile with CosmoRec to use fdm') 
#endif       
  end subroutine CMBToCAMB

 function RecomputeTransfers (A, B)
   logical RecomputeTransfers
   type(CMBParams) A, B

   RecomputeTransfers =  .not. (A%omb == B%omb .and. A%omc == B%omc .and. A%omv == B%omv .and. &
             A%omnu == B%omnu .and. A%zre == B%zre .and. A%omk == B%omk .and. A%w == B%w .and. &
               A%nnu == B%nnu .and. A%YHe == B%YHe.and. A%wa == B%wa .and. &
             A%iso_cdm_correlated == B%iso_cdm_correlated .and. A%zre_delta==B%zre_delta .and. A%ALens == B%ALens)
              
 end function RecomputeTransfers


 subroutine GetCls(CMB,Info, Cls, error)
   use ModelParams, only : ThreadNum
   use InitialPower
#ifdef DR71RG 
   use lrggettheory
   use InitialPower
   real(dl) :: getabstransferscale 
   !! BR09: this variable is for renormalizing the power spectra to the z=0 value; 
   !this is the assumption of the LRG model.
#endif     
   type(CMBParams) CMB
   integer error
   Type(ParamSetInfo) Info
   real Cls(lmax,1:num_cls_tot)
   type(CAMBParams)  P
   logical NewTransfers
   integer zix
   character(LEN=128) :: LogLine
  
  
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
         if (error==0) then
          if (Use_SN) then
            Info%Theory%SN_Loglike = SN_LnLike(CMB)
          else
            Info%Theory%SN_Loglike = 0     
          end if  
          if (Use_BAO) then
            Info%Theory%BAO_loglike = BAO_LnLike(CMB)
          else
            Info%Theory%BAO_loglike = 0
          end if 
          if (Use_HST) then
            Info%Theory%HST_Loglike = HST_LnLike(CMB)
          else
            Info%Theory%HST_Loglike = 0     
          end if
          Info%Theory%numderived = nthermo_derived
          if (nthermo_derived > max_derived_parameters) &
             call MpiStop('nthermo_derived > max_derived_parameters: increase in cmbtypes.f90')
          Info%Theory%derived_parameters(1:nthermo_derived) = ThermoDerivedParams(1:nthermo_derived)
         else
          if (stop_on_error) call MpiStop('CAMB error '//trim(global_error_message))
          if (Feedback > 0) write(*,*) 'CAMB returned error '//trim(global_error_message)
          nerrors=nerrors+1
         end if
         ncalls=ncalls+1
         if (mod(ncalls,100)==0 .and. logfile_unit/=0) then
          write (logLine,*) 'CAMB called ',ncalls, ' times; ', nerrors,' errors'
          call IO_WriteLog(logfile_unit,logLine)
         end if 
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
         if (global_error_flag/=0) then
          error=global_error_flag
          return
         end if 
           
         call SetTheoryFromCAMB(Info%Theory)
         
         if (any(Info%Theory%cl(:,1) < 0 )) then
            error = 1
            !Kill initial power spectra that go negative
            return
         end if
         
         if (compute_tensors) then
             Info%Theory%tensor_ratio_02 = TensorPower(0.002d0,1)/ScalarPower(0.002d0,1)
         else
             Info%Theory%tensor_ratio_02 = 0 
         end if
   
         if (Use_LSS) then
            Info%Theory%sigma_8 = Info%Transfers%MTrans%sigma_8(matter_power_lnzsteps,1)
#ifdef DR71RG 
            !! BR09 get lrgtheory info
            if (num_matter_power /= 0 .and. use_dr7lrg) then
                do zix = 1,matter_power_lnzsteps
                 if(zix .eq. iz0lrg .or. zix .eq. izNEARlrg .or. zix .eq. izMIDlrg .or. zix .eq. izFARlrg) then
                   call Transfer_GetMatterPowerAndNW(Info%Transfers%MTrans,&
                     Info%Theory%matter_power(:,zix),matter_power_lnzsteps-zix+1,&
                      1,matter_power_minkh, matter_power_dlnkh,num_matter_power,&
                       kmindata,getabstransferscale, &
                       Info%Theory%mpk_nw(:,zix),Info%Theory%mpkrat_nw_nl(:,zix))
                     if(zix == iz0lrg) powerscaletoz0(1) = getabstransferscale**2.0d0
                     if(zix == izNEARlrg)   powerscaletoz0(2) = powerscaletoz0(1)/getabstransferscale**2.0d0
                     if(zix == izMIDlrg)   powerscaletoz0(3) = powerscaletoz0(1)/getabstransferscale**2.0d0
                     if(zix == izFARlrg)   powerscaletoz0(4) = powerscaletoz0(1)/getabstransferscale**2.0d0
                  else  !! not an LRG redshift, so call regular function.
                   call Transfer_GetMatterPower(Info%Transfers%MTrans,&
                     Info%Theory%matter_power(:,zix),matter_power_lnzsteps-zix+1,&
                     1,matter_power_minkh, matter_power_dlnkh,num_matter_power)
                  end if
                 end do
                 if(zix == iz0lrg) powerscaletoz0(1) = 1.0d0
             else if (num_matter_power /= 0) then
            !! end BR09 get lrgtheory info
#else
            if (num_matter_power /= 0) then
#endif
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
   real, parameter :: cons =  (COBE_CMBTemp*1e6)**2*2*pi
   real nm
   integer l
   real highL_norm
   
   !The reason we store tensors separately is that can then importance sample re-computing scalars only,
   !using the stored tensor C_l
    Theory%cl=0
    do l = 2, lmax_computed_cl

       nm = cons/(l*(l+1))
       if (CMB_Lensing) then
            Theory%cl(l,1:num_clsS) =  nm*Cl_lensed(l,1, TensClOrder(1:num_clsS))
         else 
            Theory%cl(l,1:num_clsS) =  nm*Cl_scalar(l,1, scalClOrder(1:num_clsS))
       end if         

       if (num_cls>num_clsS) Theory%cl(l,num_clsS+1:num_cls) = 0 
   
       if (compute_tensors .and. l<=lmax_tensor) then
            Theory%cl_tensor(l,1:num_cls) =  nm*Cl_tensor(l,1, TensClOrder(1:num_cls))
       end if
 
       if (num_cls_ext > 0) then
           !CMB lensing potential
           !in camb Cphi is l^4 C_l, we want [l(l+1)]^2Cphi/2pi
           if (.not. CMB_lensing) call MpiStop('Must have lensing on to use lensing potential')
           Theory%cl(l,num_clsS+1) =  Cl_scalar(l,1, scalClOrder(4))*(real(l+1)**2/l**2)/twopi
           if (num_cls_ext>1) then
            !lensing-temp
            if (num_cls_ext>1) call MpiStop('SetTheoryFromCAMB: check defs for num_cls_ext>1')    
            Theory%cl(l,num_clsS+2) =   Cl_scalar(l,1, scalClOrder(5))/real(l)**3      
           end if
       end if
 
    end do 
    
    if (lmax_computed_cl/=lmax) then
        !use template for very high L theory tail, scaled not to be discontinuous
       highL_norm = Theory%cl(lmax_computed_cl,1)/highL_lensedCL_template(lmax_computed_cl,1)
       do l = lmax_computed_cl+1, lmax
          Theory%cl(l,1:num_ClsS) =  highL_norm*highL_lensedCL_template(l,1:num_clsS)
       end do 
     end if

 end subroutine SetTheoryFromCAMB

 subroutine GetClsInfo(CMB, Theory, error, DoCls, DoPk)
   use ModelParams, only : ThreadNum
#ifdef DR71RG
   use lrggettheory
   real(dl) :: getabstransferscale
   !! BR09: this variable is for renormalizing the power spectra to the z=0 value;
   !this is the assumption of the LRG model.
#endif
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
   error = global_error_flag !using error optional parameter gives seg faults on SGI
   if (error==0) then
       
      if (DoCls) then
  
       Theory%cl_tensor(2:lmax_tensor,1:num_cls) = 0
       call SetTheoryFromCAMB(Theory)
      end if

!!BR09 new addition, putting LRGs back here as well, same structure as above.  
      if (DoPk) then 
         Theory%sigma_8 = MT%sigma_8(matter_power_lnzsteps,1)

#ifdef DR71RG
         !! BR09 get lrgtheory info
         if (num_matter_power /= 0 .and. use_dr7lrg) then
             do zix = 1,matter_power_lnzsteps
              if(zix .eq. iz0lrg .or. zix .eq. izNEARlrg .or. zix .eq. izMIDlrg .or. zix .eq. izFARlrg) then
                call Transfer_GetMatterPowerAndNW(MT,&
                  Theory%matter_power(:,zix),matter_power_lnzsteps-zix+1,&
                   1,matter_power_minkh, matter_power_dlnkh,num_matter_power,&
                    kmindata,getabstransferscale, &
                    Theory%mpk_nw(:,zix),Theory%mpkrat_nw_nl(:,zix))
                     if(zix == iz0lrg) powerscaletoz0(1) = getabstransferscale**2.0d0
                     if(zix == izNEARlrg)   powerscaletoz0(2) = powerscaletoz0(1)/getabstransferscale**2.0d0
                     if(zix == izMIDlrg)   powerscaletoz0(3) = powerscaletoz0(1)/getabstransferscale**2.0d0
                     if(zix == izFARlrg)   powerscaletoz0(4) = powerscaletoz0(1)/getabstransferscale**2.0d0
               else  !! not an LRG redshift, so call regular function.
                call Transfer_GetMatterPower(MT,&
                  Theory%matter_power(:,zix),matter_power_lnzsteps-zix+1,&
                  1,matter_power_minkh, matter_power_dlnkh,num_matter_power)
               end if
             end do
             if(zix == iz0lrg) powerscaletoz0(1) = 1.0d0
              else if (num_matter_power /= 0) then
         !! end BR09 get lrgtheory info
#else
            if (num_matter_power /= 0) then
#endif
                do zix = 1,matter_power_lnzsteps
                 call Transfer_GetMatterPower(MT,&
                   Theory%matter_power(:,zix),matter_power_lnzsteps-zix+1,&
                    1,matter_power_minkh, matter_power_dlnkh,num_matter_power)
                 end do
            end if
      end if
      Theory%Age = CAMB_GetAge(P)
      Theory%numderived = nthermo_derived
      if (nthermo_derived > max_derived_parameters) &
        call MpiStop('nthermo_derived > max_derived_parameters: increase in cmbtypes.f90')
      Theory%derived_parameters(1:nthermo_derived) = ThermoDerivedParams(1:nthermo_derived)

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
   use mpk
   type(CAMBParams)  P 
   integer zix
   real redshifts(matter_power_lnzsteps)

        Threadnum =num_threads
        w_lam = -1
        wa_ppf = 0._dl
        call CAMB_SetDefParams(P)
 
        P%OutputNormalization = outNone
      
        P%WantScalars = .true.
        P%WantTensors = compute_tensors
        P%WantTransfer = Use_LSS

        P%Max_l=lmax_computed_cl
        P%Max_eta_k=lmax_computed_cl*2
      
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
        
        if (AccuracyLevel > 1 .or. HighAccuracyDefault) then
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
        if (use_LSS) then
           do zix=1, matter_power_lnzsteps
            if (zix==1) then
             redshifts(1) = 0
            else
            !Default Linear spacing in log(z+1) if matter_power_lnzsteps > 1           
             redshifts(zix) = exp( log(matter_power_maxz+1) * &
                real(zix-1)/(max(2,matter_power_lnzsteps)-1) )-1
               !put in max(2,) to stop compilers complaining of div by zero
            end if
           end do
           
           if (use_mpk) call mpk_SetTransferRedshifts(redshifts) !can modify to use specific redshifts
           if (redshifts(1) > 0.0001) call MpiStop('mpk redshifts: lowest redshift must be zero')
           do zix=1, matter_power_lnzsteps 
            !CAMB's ordering is from highest to lowest
            P%Transfer%redshifts(zix) = redshifts(matter_power_lnzsteps-zix+1)
           end do
        else 
          P%Transfer%num_redshifts = 1
          P%Transfer%redshifts(1) = 0
         end if   
        
        P%Num_Nu_Massive = 3
        P%Num_Nu_Massless = 0.046
        P%InitPower%nn = 1
        P%AccuratePolarization = num_cls/=1 
        P%Reion%use_optical_depth = .false.
        P%OnlyTransfers = .true.

        if (use_BAO) P%want_zdrag = .true. !JH
        P%want_zstar = .false. !set to true if you want CAMB to calculate exact z_star

        if (CMB_Lensing) then
            P%DoLensing = .true.
            P%Max_l = lmax_computed_cl +100 + 50 !+50 in case accuracyBoost>1 and so odd l spacing
            P%Max_eta_k = P%Max_l*2 
        end if
        
        if (HighAccuracyDefault) then
         P%Max_eta_k=max(min(P%max_l,3000)*2.5_dl,P%Max_eta_k)
        end if
        
        lensing_includes_tensors = .false.

        P%Scalar_initial_condition = initial_vector
        P%InitialConditionVector = 0
        P%InitialConditionVector(initial_adiabatic) = -1


 end subroutine InitCAMBParams
 
 
 subroutine LoadFiducialHighLTemplate
 !This should be a lensed scalar CMB power spectrum, e.g. for including at very high L where foregrounds etc. dominate anyway
   integer L
   real array(4), nm
   character(LEN=Ini_max_string_len) :: fname
  
        fname = ReadIniFilename(DefIni,'highL_theory_cl_template',DataDir,.true.)
        allocate(highL_lensedCL_template(2:lmax, num_clsS))
        call OpenTxtFile(fname,tmp_file_unit)
        do
         read(tmp_file_unit,*, end=500) L , array
         if (L>lmax) exit
         nm = 2*pi/(l*(l+1))
         if (L>=2) highL_lensedCL_template(L,1:num_clsS) = nm*array(TensClOrder(1:num_clsS))
        end do
500     close(tmp_file_unit)

        if (highL_lensedCL_template(2,1) < 100) &
           call MpiStop('highL_theory_cl_template must be in muK^2')

        if (L<lmax) call MpiStop('highL_theory_cl_template does not go to lmax')
        if (num_cls_ext>0) write(*,*) 'WARNING: zero padding ext cls in LoadFiducialHighLTemplate'

 end subroutine LoadFiducialHighLTemplate

 subroutine CMB_Initialize(Info)
   Type(ParamSetInfo) Info
   type(CAMBParams)  P 
        compute_tensors = Ini_Read_Logical('compute_tensors',.false.)
        if (num_cls==3 .and. compute_tensors) write (*,*) 'WARNING: computing tensors with num_cls=3 (BB=0)'
        CMB_lensing = Ini_Read_Logical('CMB_lensing',.false.)
        if (CMB_lensing) num_clsS = num_cls   !Also scalar B in this case
        lmax_computed_cl = Ini_Read_Int('lmax_computed_cl',lmax)
        if (lmax_computed_cl /= lmax) then
          if (lmax_tensor > lmax_computed_cl) call MpiStop('lmax_tensor > lmax_computed_cl')
          call LoadFiducialHighLTemplate
        end if
        
        if (Feedback > 0 ) then
          write (*,*) 'Computing tensors:', compute_tensors
          write (*,*) 'Doing CMB lensing:',CMB_lensing
          write(*,'(" lmax              = ",1I4)') lmax       
          write(*,'(" lmax_computed_cl  = ",1I4)') lmax_computed_cl
          if (compute_tensors) write(*,'(" lmax_tensor    = ",1I4)') lmax_tensor
          write(*,'(" Number of C_ls = ",1I4)') num_cls
        end if

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

