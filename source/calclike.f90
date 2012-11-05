module CalcLike
 use CMB_Cls
 use cmbtypes
 use Random
 use settings
 use ParamDef
 use Likelihood
 implicit none

 logical :: Use_CMB = .true.

 integer :: H0_min = 40, H0_max = 100
 real :: Omk_min = -0.3, Omk_max = 0.3
 real :: Use_min_zre = 0
 integer :: Age_min = 10, Age_max = 20
 real :: Temperature  = 1

 logical changeMask(num_params)
 logical SlowChanged, PowerChanged


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

    end if
    TestHardPriors = 0
 
  end function TestHardPriors

  function GetLogLike(Params) !Get -Ln(Likelihood) for chains
    type(ParamSet)  Params 
    Type (CMBParams) CMB
    real GetLogLike
    logical, save:: first=.true.

    if (any(Params%P > Scales%PMax) .or. any(Params%P < Scales%PMin)) then
       GetLogLike = logZero
       return
    end if

    if (generic_mcmc) then
        GetLogLike = GenericLikelihoodFunction(Params) 
        if (GetLogLike /= LogZero) GetLogLike = GetLogLike + getLogPriors(Params%P)
        if (GetLogLike /= LogZero) GetLogLike = GetLogLike/Temperature
    else
      GetLogLikePost  = TestHardPriors(CMB, Params%Info)
      if (GetLogLikePost == logZero) return
      call ParamsToCMBParams(Params%P,CMB)
      if (first) then
           changeMask = .true.
           first = .false.
      else
           changeMask = Params%Info%lastParamArray/=Params%P
      end if

     if (CalculateRequiredTheoryChanges(CMB, Params)) then
       GetLogLike = GetLogLikeWithTheorySet(CMB, Params%P, Params%Info,.false.)
     else
       GetLogLike = logZero
     end if

     if (GetLogLike/=logZero) Params%Info%lastParamArray = Params%P
    end if 

  end function GetLogLike

  function GetLogLikePost(CMB, P, Info)
  !for importance sampling where theory may be pre-stored
    real GetLogLikePost
    Type (CMBParams) CMB
    Type(ParamSetInfo) Info
    real P(num_params)
   
    !need to init background correctly if new BAO etc

  end function GetLogLikePost
  
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
  
  function CalculateRequiredTheoryChanges(CMB, Params)
    type(ParamSet)  Params 
    type (CMBParams) CMB
    integer error

     SlowChanged = any(changeMask(1:num_hard))
     PowerChanged = any(changeMask(index_initpower:index_initpower+num_init_power-1))
     error=0
     if (Use_CMB .or. Use_LSS) then
         if (SlowChanged) then
           call GetNewTransferData(CMB, Info, error)
           if (error/=0) return .false.
         end if
         if (SlowChanged .or. PowerChanged) then
           call GetNewPowerData(CMB, Info, error)
           if (error/=0) return .false.
         end if
     end if

  end function CalculateRequiredTheoryChanges

  
  function GetLogLikeWithTheorySet(CMB, P, Info) result(logLike)
    real logLike
    Type (CMBParams) CMB
    Type(ParamSetInfo) Info
    integer error
    real P(num_params)
    real itemLike
    Type(DataLikelihood) :: like
    integer i
    logical backgroundSet
    
    backgroundSet = slowChanged

    do i= 1, Likihoods%count
     like = likelihoods%Items(i)
     if (any(like%dependent_params .and. changeMask )) then
          if (like%needs_background_functions .and. .not. backgroundSet) then
              call SetTheoryForBackground(CMB)
              backgroundSet = .true.

          end if
          itemLike = like%LogLike(CMB, Info%Theory)
          if (itemLike == logZero) return logZero
          Info%Likelihooods(i) = itemLike
     end if
    end do
    logLike = sum(Info%Likelihoods(1:Likelihoods%Count))
    if (logLike /= LogZero) logLike = logLike + getLogPriors(Params%P)
    if (logLike /= LogZero) logLike = logLike/Temperature
    
  end function GetLogLikeWithTheorySet


  if (Use_mpk) logLike = logLike + LSSLnLike(CMB, Info%theory)
          if (Use_WeakLen) logLike = logLike + WeakLenLnLike(CMB, Info%theory)     
          if (Use_Lya) logLike = logLike +  LSS_Lyalike(CMB, Info%Theory)
          if (logLike >= logZero) logLike = logZero
         end if
         if (Use_SN .and. logLike /= logZero ) then
            if (Info%Theory%SN_loglike /= 0) then
             logLike = logLike + Info%Theory%SN_loglike
            else
             logLike = logLike + SN_LnLike(CMB)
            end if
               !Assume computed only every time hard parameters change
         end if
         if (Use_BAO .and. logLike /= logZero ) then
            if (Info%Theory%BAO_loglike /= 0) then
             logLike = logLike + Info%Theory%BAO_loglike
            else
             logLike = logLike + BAO_LnLike(CMB)
            end if
               !Assume computed only every time hard parameters change
         end if
         if (Use_HST .and. logLike /= logZero) then
            if (Info%Theory%HST_loglike /= 0) then
             logLike = logLike + Info%Theory%HST_loglike
            else
             logLike = logLike + HST_LnLike(CMB)
            end if
         end if  


end module CalcLike
