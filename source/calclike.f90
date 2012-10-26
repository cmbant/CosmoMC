module CalcLike
 use CMB_Cls
 use cmbtypes
 use cmbdata
 use mpk
 use Random
 use settings
 use ParamDef
 use snovae
 use WeakLen
 use Lya
#ifdef CLIK
 use cliklike
#endif
 implicit none

 logical :: Use_Age_Tophat_Prior = .true.
 logical :: Use_CMB = .true.
 logical :: Use_Clusters = .false.
 
 integer :: H0_min = 40, H0_max = 100
 real :: Omk_min = -0.3, Omk_max = 0.3
 real :: Use_min_zre = 0
 integer :: Age_min = 10, Age_max = 20
 real :: Temperature  = 1

contains

  function GenericLikelihoodFunction(Params) 
    type(ParamSet)  Params 
    real :: GenericLikelihoodFunction
   
  
   !Used when you want to plug in your own CMB-independent likelihood function:
   !set generic_mcmc=.true. in settings.f90, then write function here returning -Ln(Likelihood)
   !Parameter array is Params%P
    !!!!
   ! GenericLikelihoodFunction = greatLike(Params%P)
   ! GenericLikelihoodFunction = LogZero 
    call MpiStop('GenericLikelihoodFunction: need to write this function!')
    GenericLikelihoodFunction=0

  end function GenericLikelihoodFunction

  
  function TestHardPriors(CMB, Info) 
    real TestHardPriors
    real Age
    Type (CMBParams) CMB
    Type(ParamSetInfo) Info

    TestHardPriors = logZero
 
    if (.not. generic_mcmc) then
     if (CMB%H0 < H0_min .or. CMB%H0 > H0_max) return
     if (CMB%omk < Omk_min .or. CMB%omk > Omk_max .or. CMB%Omv < 0) return
     if (CMB%zre < Use_min_zre) return

     Age = GetAge(CMB, Info)
      !This sets up parameters in CAMB, so do not delete unless you are careful!
 
     if (Use_Age_Tophat_Prior .and. (Age < Age_min .or. Age > Age_max) .or. Age < 0) return
    
    end if
    TestHardPriors = 0
 
  end function TestHardPriors

  function GetLogLike(Params) !Get -Ln(Likelihood)
    type(ParamSet)  Params 
    Type (CMBParams) CMB
    real GetLogLike
 
    if (any(Params%P > Scales%PMax) .or. any(Params%P < Scales%PMin)) then
       GetLogLike = logZero
        return
    end if

    if (generic_mcmc) then
        GetLogLike = GenericLikelihoodFunction(Params) 
        if (GetLogLike /= LogZero) GetLogLike = GetLogLike + getLogPriors(Params%P)
        if (GetLogLike /= LogZero) GetLogLike = GetLogLike/Temperature
    else

     call ParamsToCMBParams(Params%P,CMB)
     GetLogLike = GetLogLikePost(CMB, Params%P, Params%Info,.false.)

    end if 

  end function GetLogLike
  
  function getLogPriors(P) result(logLike)
  integer i
  real, intent(in) :: P(num_params)
  real logLike
  
  logLike=0
  do i=1,num_params
        if (Scales%PWidth(i)/=0 .and. GaussPriors%std(i)/=0) then
          logLike = logLike + ((P(i)-GaussPriors%mean(i))/GaussPriors%std(i))**2
        end if
  end do
  logLike=logLike/2

  end function getLogPriors
  
    
  function GetLogLikePost(CMB, P, Info, HasTheory) 
    real GetLogLikePost
    Type (CMBParams) CMB
    real P(num_params)
    Type(ParamSetInfo) Info
    logical, intent(in) :: HasTheory
    real acl(lmax,num_cls_tot)
    integer error
    logical SlowChanged
    logical, save:: first=.true.

    if (generic_mcmc) stop 'GetLogLikePost: not supported for generic'
    
    SlowChanged = first .or. any(Info%lastParamArray(1:num_hard)/=P(1:num_hard))
!    NormChanged = SlowChanged .or. any(Info%lastParamArray(num_hard)/=P(1:num_hard))
    
    first= .false.

    NewTransfers = .false.
    GetLogLikePost  = TestHardPriors(CMB, Info)
    if ( GetLogLikePost /= logZero) then
       if (Use_CMB .or. Use_LSS) then 
          if (HasTheory) then
           error =0
           call SetTheoryBackground(CMB, Info)
          else
           NewTransfers = GetCMBTheory(CMB, Info, error)
          end if
         if (error /= 0) then
          GetLogLikePost = logZero 
         else
          if (Use_CMB) then
            call ClsFromTheoryData(Info%Theory, CMB, acl)
            GetLogLikePost = &
              CMBLnLike(acl, CMB%norm(norm_freq_ix:norm_freq_ix+num_freq_params-1),CMB%nuisance) + GetLogLikePost
#ifdef CLIK
!Assuming CAMspec nuisance parameters are set as freq_params(2:34), PLik nuisance parameters as 
!freq_params(35:44), ACT/SPT as freq_params(45:65)
            if (Use_clik) then
             GetLogLikePost = GetLogLikePost + clik_lnlike(dble(acl),dble(CMB%norm(norm_freq_ix+1:norm_freq_ix+num_freq_params-1)))
            end if
#endif
          end if
          if (Use_mpk) GetLogLikePost = GetLogLikePost + LSSLnLike(CMB, Info%theory)
          if (Use_WeakLen) GetLogLikePost = GetLogLikePost + WeakLenLnLike(CMB, Info%theory)     
          if (Use_Lya) GetLogLikePost = GetLogLikePost +  LSS_Lyalike(CMB, Info%Theory)
          if (GetLogLikePost >= logZero) GetLogLikePost = logZero
         end if
         if (Use_SN .and. GetLogLikePost /= logZero ) then
            if (Info%Theory%SN_loglike /= 0) then
             GetLogLikePost = GetLogLikePost + Info%Theory%SN_loglike
            else
             GetLogLikePost = GetLogLikePost + SN_LnLike(CMB)
            end if
               !Assume computed only every time hard parameters change
         end if
         if (Use_BAO .and. GetLogLikePost /= logZero ) then
            if (Info%Theory%BAO_loglike /= 0) then
             GetLogLikePost = GetLogLikePost + Info%Theory%BAO_loglike
            else
             GetLogLikePost = GetLogLikePost + BAO_LnLike(CMB)
            end if
               !Assume computed only every time hard parameters change
         end if
         if (Use_HST .and. GetLogLikePost /= logZero) then
            if (Info%Theory%HST_loglike /= 0) then
             GetLogLikePost = GetLogLikePost + Info%Theory%HST_loglike
            else
             GetLogLikePost = GetLogLikePost + HST_LnLike(CMB)
            end if
         end if  
     
       else !No transfer functions
           call SetTheoryBackground(CMB, Info)

           if (Use_SN) GetLogLikePost = GetLogLikePost + SN_LnLike(CMB)
           if (Use_BAO)then
                !From Jason Dosset           
                !JD had to change below so new BAO will work without CMB
                !previous way z_drag was not calculated without calling get_cls.
                if (RecomputeTransfers(CMB, Info%LastParams))  then
                    call CMBToCAMB(CMB, P)
                    call CAMBParams_Set(P)
                    call InitVars
                    Info%Theory%BAO_loglike = Bao_lnLike(CMB)
                    Info%LastParams = CMB
                end if
                GetLogLikePost = GetLogLikePost + Bao_lnLike(CMB)
           end if
           if (Use_HST) GetLogLikePost = GetLogLikePost + HST_LnLike(CMB)
         end if
         
     
      if (Use_Clusters .and. GetLogLikePost /= LogZero) then
          GetLogLikePost = GetLogLikePost + &
                 (Info%Theory%Sigma_8-0.9)**2/(2*0.05**2)
          stop 'Write your cluster prior in calclike.f90 first!'
      end if

     if (GetLogLikePost /= LogZero) GetLogLikePost = GetLogLikePost + getLogPriors(Params%P)
     if (GetLogLikePost /= LogZero) GetLogLikePost = GetLogLikePost/Temperature

    end if

  end function GetLogLikePost

end module CalcLike
