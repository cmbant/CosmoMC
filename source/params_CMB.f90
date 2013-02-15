!Parameterization using theta = r_s/D_a instead of H_0, and tau instead of z_re
!and log(A_s) instead of A_s
!Less general, but should give better performance
!
!The well-determined parameter A_s exp(-2tau) should be found by the covariance matrix
!parameter 3 is 100*theta, parameter 4 is tau, others same as params_H except A->log(A)
!Theta is much better constrained than H_0
!
!AL Jul 2005 - fixed bug which screwed up tau values for later importance sampling
!AL Feb 2004 - fixed compiler compatibility issue with **-number, typo in zstar
!AL Dec 2003 - renamed from params_Kowosky, changed to tau - log(A)
!AL Sept 2003 - fixed bug in declaration of dtauda routine
!AL June 2003
!Assumes prior 0.4 < h < 1

   function CMBToTheta(CMB)
     use settings
     use cmbtypes
     use ModelParams
     use CMB_Cls
     use Precision
     implicit none
     Type(CMBParams) CMB
     real CMBToTheta
     integer error
  
     call InitCAMB(CMB,error,.false.)
     CMBToTheta = CosmomcTheta()
!       print *,'z* = ',zstar, 'r_s = ',rs, 'DA = ',DA, rs/DA

  end function CMBToTheta




!Mapping between array of power spectrum parameters and CAMB
     subroutine SetCAMBInitPower(P,CMB,in)
       use camb
       use settings
       use cmbtypes
       use CMB_Cls
       implicit none
       type(CAMBParams)  P
       Type(CMBParams) CMB

       integer, intent(in) :: in

!MODIFIED P(K)
!       if (Power_Name == 'power_tilt') then
       if ((Power_Name == 'power_tilt') .or. (Power_Name == 'power_tilt_modpk')) then
!END MODIFIED P(K)

       P%InitPower%k_0_scalar = pivot_k
       P%InitPower%k_0_tensor = pivot_k

       P%InitPower%ScalarPowerAmp(in) = cl_norm*CMB%norm(norm_As)
       P%InitPower%rat(in) = CMB%norm(norm_amp_ratio)
        
       P%InitPower%an(in) = CMB%InitPower(1)
       P%InitPower%ant(in) = CMB%InitPower(2)
       if (P%InitPower%rat(in)>0 .and. .not. compute_tensors) &
        call MpiStop('computing r>0 but compute_tensors=F')
       P%InitPower%n_run(in) = CMB%InitPower(3)
       if (inflation_consistency) then
         P%InitPower%ant(in) = - CMB%norm(norm_amp_ratio)/8.
          !note input n_T is ignored, so should be fixed (to anything)
       end if
       else
         stop 'params_CMB:Wrong initial power spectrum'
       end if

    end subroutine SetCAMBInitPower
 

 subroutine SetForH(Params,CMB,H0, firsttime)
     use settings
     use cmbtypes
     use CMB_Cls
     use bbn
!MODIFIED P(K)
     use modpkparams, only : instreheat
!END MODIFIED P(K)
     implicit none
     real Params(num_Params)
     logical, intent(in) :: firsttime
     Type(CMBParams) CMB
     real h2,H0
     
  
    CMB%H0=H0
    if (firsttime) then
        CMB%reserved = 0
        CMB%ombh2 = Params(1)    
        CMB%tau = params(4) !tau, set zre later
        CMB%Omk = Params(5)
        if (neutrino_param_mnu) then
        !Params(6) is now mnu, params(2) is omch2
         CMB%omnuh2=Params(6)/93.04
         if (CMB%omnuh2 > 0 .and. (Params(9) < 3 .or. Params(9)>3.1)) &
            call MpiStop('params_CMB: change for non-standard nnu with massive nu')
         CMB%omch2 = Params(2)
         CMB%omdmh2 = CMB%omch2+ CMB%omnuh2
         CMB%nufrac=CMB%omnuh2/CMB%omdmh2
        else
         CMB%omdmh2 = Params(2)
         CMB%nufrac=Params(6)
         CMB%omnuh2 = CMB%omdmh2*CMB%nufrac
         CMB%omch2 = CMB%omdmh2 - CMB%omnuh2
        end if 
        CMB%w = Params(7)
        CMB%wa = Params(8)
        CMB%nnu = Params(9) !3.046
        
        if (bbn_consistency) then
         CMB%YHe = yp_bbn(CMB%ombh2,CMB%nnu  - 3.046)
        else
         !e.g. set from free parameter..
         CMB%YHe  =Params(10)
         !call MpiStop('params_CMB: YHe not free parameter in default parameterization')
        end if
        
        CMB%iso_cdm_correlated =  Params(11)
        CMB%zre_delta = Params(12)
        CMB%ALens = Params(13)
        CMB%fdm = Params(14)
        
!MODIFIED P(K)
        if (.not.instreheat) CMB%N_pivot = Params(index_vpar)
        CMB%vparams(1:num_vpar-1) = Params(index_vpar+1:index_vpar+num_vpar-1)
!END MODIFIED P(K)

        CMB%InitPower(1:num_initpower) = Params(index_initpower:index_initpower+num_initPower-1)
        CMB%norm(1) = exp(Params(index_norm))
        CMB%norm(2:num_norm) = Params(index_norm+1:index_norm+num_norm-1)
        CMB%nuisance(1:num_nuisance_params) = Params(index_nuisance:index_nuisance+num_nuisance_params-1)
    end if
    
    CMB%h = CMB%H0/100
    h2 = CMB%h**2
    CMB%omb = CMB%ombh2/h2
    CMB%omc = CMB%omch2/h2
    CMB%omnu = CMB%omnuh2/h2
    CMB%omdm = CMB%omdmh2/h2
    CMB%omv = 1- CMB%omk - CMB%omb - CMB%omdm

 end  subroutine SetForH

 subroutine ParamsToCMBParams(Params, CMB)
     use settings
     use cmbtypes
     use CMB_Cls

     implicit none
     real Params(num_params)
     real, save :: LastParams(num_params) = 0.
     real, save :: LastH0, Lastzre

     Type(CMBParams) CMB
     real DA
     real  D_b,D_t,D_try,try_b,try_t, CMBToTheta, lasttry
     external CMBToTheta

!MODIFIED P(K)
!     if (all(Params(1:num_hard) == Lastparams(1:num_hard))) then
     if (all(Params(1:num_hard) == Lastparams(1:num_hard)) .and. &
          all(Params(index_vpar:index_vpar+num_vpar-1) == LastParams(index_vpar:index_vpar+num_vpar-1))) then
!END MODIFIED P(K)
       call SetForH(Params,CMB,LastH0, .true.)
       CMB%zre = Lastzre
     else

     DA = Params(3)/100
     try_b = 40
     call SetForH(Params,CMB,try_b, .true.)
     D_b = CMBToTheta(CMB)
     try_t = 100
     call SetForH(Params,CMB,try_t, .false.)
     D_t = CMBToTheta(CMB)
     if (DA < D_b .or. DA > D_t) then
      cmb%H0=0 !Reject it
     else
     lasttry = -1
     do
            call SetForH(Params,CMB,(try_b+try_t)/2, .false.)
            D_try = CMBToTheta(CMB)
               if (D_try < DA) then
                  try_b = (try_b+try_t)/2
               else
                  try_t = (try_b+try_t)/2
               end if
               if (abs(D_try - lasttry)< 1e-7) exit
              lasttry = D_try
     end do

    !!call InitCAMB(CMB,error)
    CMB%zre = GetZreFromTau(CMB, CMB%tau)       
  
    LastH0 = CMB%H0
    Lastzre = CMB%zre
    LastParams = Params
    end if

  
     end if
 
   end subroutine ParamsToCMBParams

   subroutine CMBParamsToParams(CMB, Params)
     use settings
     use cmbtypes
     implicit none
     real Params(num_Params)
     Type(CMBParams) CMB
     real CMBToTheta
     external CMBToTheta
 
      Params(1) = CMB%ombh2 
 
      Params(3) = CMBToTheta(CMB)*100
      Params(4) = CMB%tau
      Params(5) = CMB%omk 
      
      if (neutrino_param_mnu) then
          Params(2) = CMB%omch2
          Params(6) = CMB%omnuh2*93.04
      else
          Params(2) = CMB%omdmh2
          Params(6) = CMB%nufrac
      end if
      Params(7) = CMB%w
      Params(8) = CMB%wa
      Params(9) = CMB%nnu
      if (bbn_consistency) then
          Params(10) = 0
      else
          Params(10) = CMB%YHe
      endif
      Params(11) = CMB%iso_cdm_correlated
      Params(12) = CMB%zre_delta  
      Params(13) = CMB%ALens  
      Params(14) = CMB%fdm  
      
      Params(index_initpower:index_initpower+num_initpower-1) =CMB%InitPower(1:num_initpower) 
      Params(index_norm) = log(CMB%norm(1))
      Params(index_norm+1:index_norm+num_norm-1) = CMB%norm(2:num_norm)
      Params(index_nuisance:index_nuisance+num_nuisance_params-1)=CMB%nuisance(1:num_nuisance_params) 

!MODIFIED P(K)
      Params(index_vpar) = CMB%N_pivot
      Params(index_vpar+1:index_vpar+num_vpar-1) = CMB%vparams(1:num_vpar-1)
!END MODIFIED P(K)

   end subroutine CMBParamsToParams

   subroutine SetParamNames(Names)
    use settings
    use ParamNames
    Type(TParamNames) :: Names
 
    if (ParamNamesFile /='') then
      call ParamNames_init(Names, ParamNamesFile)
    else
     if (generic_mcmc) then
      Names%nnames=0
      if (Feedback>0) write (*,*) 'edit SetParamNames in params_CMB.f90 if you want to use named params'
     else
#ifdef CLIK
       call ParamNames_init(Names, trim(LocalDir)//'clik.paramnames')
#else
       call ParamNames_init(Names, trim(LocalDir)//'params_CMB.paramnames')
#endif
     end if
    end if
   end subroutine SetParamNames

 
  function CalcDerivedParams(P, derived) result (num_derived)
     use settings
     use cmbtypes
     use ParamDef
     use Lists
     implicit none
     Type(real_pointer) :: derived
     Type(ParamSet) P
     Type(CMBParams) CMB
     real r10
     integer num_derived 
     
     num_derived = 12 +  P%Info%Theory%numderived
   
     allocate(Derived%P(num_derived))
   
      call ParamsToCMBParams(P%P,CMB)

      if (lmax_tensor /= 0 .and. compute_tensors) then
          r10 = P%Info%Theory%cl_tensor(10,1)/P%Info%Theory%cl(10,1)
      else
        r10 = 0
      end if

      derived%P(1) = CMB%omv
      derived%P(2) = CMB%omdm+CMB%omb
      derived%P(3) = P%Info%Theory%Sigma_8
      derived%P(4) = CMB%zre
      derived%P(5) = r10
      derived%P(6) = CMB%H0
      derived%P(7) = P%Info%Theory%tensor_ratio_02
      derived%P(8) = cl_norm*CMB%norm(norm_As)*1e9
      derived%P(9) = CMB%omdmh2 + CMB%ombh2
      derived%P(10)= (CMB%omdmh2 + CMB%ombh2)*CMB%h
      derived%P(11)= CMB%Yhe !value actually used, may be set from bbn consistency
      derived%P(12)= derived%P(8)*exp(-2*CMB%tau)  !A e^{-2 tau}
      
      derived%P(13:num_derived) = P%Info%Theory%derived_parameters(1: P%Info%Theory%numderived)
      
  end function CalcDerivedParams
  
  subroutine WriteParams_nest(P, nPar, Cube)
     use settings
     use cmbtypes
     use ParamDef
     use IO
     use Lists
!MODIFIED P(K)
     use modpkparams, only : use_modpk, instreheat
!END MODIFIED P(K)
     implicit none
     Type(ParamSet) P
     Type(CMBParams) CMB
     integer :: mult
!MODIFIED P(K)
    integer il, nPar, i, j
  ! set range in l to write for each model if feedback=3
    integer, parameter :: wp_lmin = 2
    integer, parameter :: wp_lmax = 1000
!END MODIFIED P(K)
     real, allocatable :: output_array(:)
     real*8 :: Cube(1:nPar)
     Type(real_pointer) :: derived
     integer numderived 
     integer CalcDerivedParams
     external CalcDerivedParams
  
     mult = 0
     if (outfile_handle ==0) return
  
    if (generic_mcmc) then

      !call IO_OutputChainRow(outfile_handle, mult, like, P%P)
     
    else
    
      numderived = CalcDerivedParams(P, derived)

!MODIFIED P(K)
      if (use_modpk) then
         allocate(output_array(num_real_params + numderived + nuisance_params_used + 7 ))
      else
         allocate(output_array(num_real_params + numderived + nuisance_params_used ))
      end if
!END MODIFIED P(K)
      do i=1, num_params_used
		Cube(i) = P%P(params_used(i))
      end do

      output_array(1:num_real_params) =  P%P(1:num_real_params)
      output_array(num_real_params+1:num_real_params+numderived) =  derived%P
      deallocate(derived%P)

      j=1
      do i=num_params_used+1, num_params_used+numderived
                Cube(i) = output_array(num_real_params+j)
		j=j+1
      end do

!MODIFIED P(K)
      if (use_modpk) then
         output_array(num_real_params+numderived+1) = P%Info%Theory%modpk_Npivot
         output_array(num_real_params+numderived+2) = P%Info%Theory%modpk_ns
         output_array(num_real_params+numderived+3) = P%Info%Theory%modpk_nt
         output_array(num_real_params+numderived+4) = P%Info%Theory%modpk_nrun
         output_array(num_real_params+numderived+5) = log(1.d10*P%Info%Theory%modpk_As)
         output_array(num_real_params+numderived+6) = P%Info%Theory%modpk_r
         output_array(num_real_params+numderived+7) = P%Info%Theory%modpk_w
         j=numderived+1
         do i=num_params_used+numderived+1, num_params_used+numderived+7
                Cube(i) = output_array(num_real_params+j)
		j=j+1
         end do
     end if
!END MODIFIED P(K)

      if (nuisance_params_used>0) then
!MODIFIED P(K)
         if (use_modpk) then
            output_array(num_real_params+numderived+7+1:num_real_params+numderived+7+nuisance_params_used) = &
                 P%P(num_real_params+1:num_real_params+nuisance_params_used) 
         else
            output_array(num_real_params+numderived+1:num_real_params+numderived+nuisance_params_used) = &
                 P%P(num_real_params+1:num_real_params+nuisance_params_used) 
         end if
         Cube(j:j+nuisance_params_used-1) = &
              P%P(num_real_params+1:num_real_params+nuisance_params_used)
!END MODIFIED P(K)
      end if

!MODIFIED P(K)
      if (use_modpk) Cube(j+nuisance_params_used:nPar) = 0.d0 !zero out crud in the rest of Cube
!END MODIFIED P(K)
 
      !call IO_OutputChainRow(outfile_handle, mult, like, output_array)
      deallocate(output_array)           
    end if

!MODIFIED P(K)
! write large scale C_l's for each sample to files
! (useful for testing or for computing distribution of C_l)
     if (Feedback > 2) then
        open(unit=13,file=trim(rootname)//'.clt',form='formatted',status='unknown',position='append')
        if (compute_tensors) then
           write(13,'(10000E15.5)') mult,(il*(il+1)* &
                (P%Info%Theory%cl(il,1)+P%Info%Theory%cl_tensor(il,1))/(2.*pi), il=wp_lmin,wp_lmax)
        else
           write(13,'(10000E15.5)') mult,(il*(il+1)*P%Info%Theory%cl(il,1)/(2.*pi), il=wp_lmin,wp_lmax)
        end if
        close(13)
        open(unit=13,file=trim(rootname)//'.cle',form='formatted',status='unknown',position='append')
        if (compute_tensors) then
           write(13,'(10000E15.5)') mult,(il*(il+1)* &
                (P%Info%Theory%cl(il,3)+P%Info%Theory%cl_tensor(il,3))/(2.*pi), il=wp_lmin,wp_lmax)
        else
           write(13,'(10000E15.5)') mult,(il*(il+1)*P%Info%Theory%cl(il,3)/(2.*pi), il=wp_lmin,wp_lmax)
        end if
        close(13)
        open(unit=13,file=trim(rootname)//'.clc',form='formatted',status='unknown',position='append')
        if (compute_tensors) then
           write(13,'(10000E15.5)') mult,(il*(il+1)* &
                (P%Info%Theory%cl(il,2)+P%Info%Theory%cl_tensor(il,2))/(2.*pi), il=wp_lmin,wp_lmax)
        else
           write(13,'(10000E15.5)') mult,(il*(il+1)*P%Info%Theory%cl(il,2)/(2.*pi), il=wp_lmin,wp_lmax)
        end if
        close(13)
     end if
!END MODIFIED P(K)

  end  subroutine WriteParams_nest


  subroutine WriteParams(P, mult, like)
     use settings
     use cmbtypes
     use ParamDef
     use IO
     use Lists
!MODIFIED P(K)
     use CAMB
     use modpkparams, only : use_modpk, instreheat, &
          N_pivot, modpk_ns, modpk_nt, modpk_nrun, modpk_As, modpk_r, modpk_w
     use camb_interface, only : pk_bad, pk_initialized
     use InitialPower, only : InitializePowers
     use CMB_Cls, only : Use_CMB
     use mpk, only : use_mpk
!END MODIFIED P(K)
     implicit none
     Type(ParamSet) P
     real, intent(in) :: mult, like
!MODIFIED P(K)
     Type(CMBParams) CMB
     Type(CAMBParams) PP
     integer il
  ! set range in l to write for each model if feedback=3
     integer, parameter :: wp_lmin = 2
     integer, parameter :: wp_lmax = 1000
!END MODIFIED P(K)
     real, allocatable :: output_array(:)
     Type(real_pointer) :: derived
     integer numderived 
     integer CalcDerivedParams
     external CalcDerivedParams
  
    if (outfile_handle ==0) return
  
    if (generic_mcmc) then

      call IO_OutputChainRow(outfile_handle, mult, like, P%P)
     
    else
    
      numderived = CalcDerivedParams(P, derived)

!MODIFIED P(K)
      if (use_modpk) then
         allocate(output_array(num_real_params + numderived + nuisance_params_used + 7 ))
      else
         allocate(output_array(num_real_params + numderived + nuisance_params_used ))
      end if
!END MODIFIED P(K)
      output_array(1:num_real_params) =  P%P(1:num_real_params)
      output_array(num_real_params+1:num_real_params+numderived) =  derived%P
      deallocate(derived%P)
!MODIFIED P(K)
      if (use_modpk) then
         if ((.not.Use_CMB) .and. (.not.use_mpk)) then
         ! compute P(k) to get derived parameters
         ! *** This is computed for the CURRENT parameters, not the *****
         ! *** TRIAL parameters as for the usual sampling, so pk_bad ****
         ! *** reports in the log file will list the wrong parameters ***
            pk_initialized = .false.
            call InitializePowers(PP%InitPower,-CMB%Omk/(299792.458_dl/CMB%H0)**2)
            output_array(num_real_params+numderived+1) = N_pivot
            output_array(num_real_params+numderived+2) = modpk_ns
            output_array(num_real_params+numderived+3) = modpk_nt
            output_array(num_real_params+numderived+4) = modpk_nrun
            output_array(num_real_params+numderived+5) = log(1.d10*modpk_As)
            output_array(num_real_params+numderived+6) = modpk_r
            output_array(num_real_params+numderived+7) = modpk_w
         else
            output_array(num_real_params+numderived+1) = P%Info%Theory%modpk_Npivot
            output_array(num_real_params+numderived+2) = P%Info%Theory%modpk_ns
            output_array(num_real_params+numderived+3) = P%Info%Theory%modpk_nt
            output_array(num_real_params+numderived+4) = P%Info%Theory%modpk_nrun
            output_array(num_real_params+numderived+5) = log(1.d10*P%Info%Theory%modpk_As)
            output_array(num_real_params+numderived+6) = P%Info%Theory%modpk_r
            output_array(num_real_params+numderived+7) = P%Info%Theory%modpk_w
         end if
      end if
!END MODIFIED P(K)

      if (nuisance_params_used>0) then
!MODIFIED P(K)
         if (use_modpk) then
            output_array(num_real_params+numderived+7+1:num_real_params+numderived+7+nuisance_params_used) = &
                 P%P(num_real_params+1:num_real_params+nuisance_params_used) 
         else
            output_array(num_real_params+numderived+1:num_real_params+numderived+nuisance_params_used) = &
                 P%P(num_real_params+1:num_real_params+nuisance_params_used) 
         end if
!END MODIFIED P(K)
      end if
 
 !MODIFIED P(K) 
      if (pk_bad==0) call IO_OutputChainRow(outfile_handle, mult, like, output_array)
!      call IO_OutputChainRow(outfile_handle, mult, like, output_array)
!END MODIFIED P(K)
      deallocate(output_array)           
    end if

!MODIFIED P(K)
! write large scale C_l's for each sample to files
! (useful for testing or for computing distribution of C_l)
     if (Feedback > 2) then
        open(unit=13,file=trim(rootname)//'.clt',form='formatted',status='unknown',position='append')
	if (compute_tensors) then
           write(13,'(10000E15.5)') mult,(il*(il+1)* &
                (P%Info%Theory%cl(il,1)+P%Info%Theory%cl_tensor(il,1))/(2.*pi), il=wp_lmin,wp_lmax)
        else
           write(13,'(10000E15.5)') mult,(il*(il+1)*P%Info%Theory%cl(il,1)/(2.*pi), il=wp_lmin,wp_lmax)
        end if
        close(13)
        open(unit=13,file=trim(rootname)//'.cle',form='formatted',status='unknown',position='append')
	if (compute_tensors) then
           write(13,'(10000E15.5)') mult,(il*(il+1)* &
                (P%Info%Theory%cl(il,3)+P%Info%Theory%cl_tensor(il,3))/(2.*pi), il=wp_lmin,wp_lmax)
        else
           write(13,'(10000E15.5)') mult,(il*(il+1)*P%Info%Theory%cl(il,3)/(2.*pi), il=wp_lmin,wp_lmax)
        end if
        close(13)
        open(unit=13,file=trim(rootname)//'.clc',form='formatted',status='unknown',position='append')
	if (compute_tensors) then
           write(13,'(10000E15.5)') mult,(il*(il+1)* &
                (P%Info%Theory%cl(il,2)+P%Info%Theory%cl_tensor(il,2))/(2.*pi), il=wp_lmin,wp_lmax)
        else
           write(13,'(10000E15.5)') mult,(il*(il+1)*P%Info%Theory%cl(il,2)/(2.*pi), il=wp_lmin,wp_lmax)
        end if
        close(13)
     end if
!END MODIFIED P(K)

  end  subroutine WriteParams




  subroutine WriteParamsAndDat(P, mult, like)
     use settings
     use cmbtypes
     use ParamDef
     use IO
     use Lists
     implicit none
     Type(ParamSet) P
     real, intent(in) :: mult, like
     real,allocatable :: output_array(:)
     Type(real_pointer) :: derived
     integer numderived 
     integer CalcDerivedParams
     external CalcDerivedParams
         
    if (outfile_handle ==0) return

      numderived = CalcDerivedParams(P, derived)

      allocate(output_array(num_real_params + numderived + num_matter_power ))
      output_array(1:num_real_params) =  P%P(1:num_real_params)
      output_array(num_real_params+1:num_real_params+numderived) =  derived%P
      deallocate(derived%P)

      output_array(num_real_params+numderived+1:num_real_params+numderived+num_matter_power) = &
        P%Info%Theory%matter_power(:,1) 

      call IO_OutputChainRow(outfile_handle, mult, like, output_array)
      deallocate(output_array)

  end  subroutine WriteParamsAndDat
