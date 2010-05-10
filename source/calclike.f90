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
 implicit none

 logical :: Use_HST = .true.
 logical :: Use_Age_Tophat_Prior = .true.
 logical :: Use_CMB = .true.
 logical :: Use_BBN = .false.
 logical :: Use_Clusters = .false.
 
 integer :: H0_min = 40, H0_max = 100
 real :: Omk_min = -0.3, Omk_max = 0.3
 real :: Use_min_zre = 0
 integer :: Age_min = 10, Age_max = 20
 real :: Temperature  = 1

contains
  
  function GetLogPrior(CMB, Info) !Get -Ln(Prior)
    real GetLogPrior
    real Age
    Type (CMBParams) CMB
    Type(ParamSetInfo) Info

    GetLogPrior = logZero
 
    if (CMB%H0 < H0_min .or. CMB%H0 > H0_max) return
    if (CMB%omk < Omk_min .or. CMB%omk > Omk_max .or. CMB%Omv < 0) return
    if (CMB%zre < Use_min_zre) return

    Age = GetAge(CMB, Info)
     !This sets up parameters in CAMB, so do not delete unless you are careful!
 
    if (Use_Age_Tophat_Prior .and. (Age < Age_min .or. Age > Age_max) .or. Age < 0) return
    GetLogPrior = 0
 
  end function GetLogPrior

  function GetLogLike(Params) !Get -Ln(Likelihood)
    type(ParamSet)  Params 
    Type (CMBParams) CMB
    real GetLogLike
    real dum(1,1)
 
    if (any(Params%P > Scales%PMax) .or. any(Params%P < Scales%PMin)) then
       GetLogLike = logZero
        return
    end if

    call ParamsToCMBParams(Params%P,CMB)

    GetLogLike = GetLogLikePost(CMB, Params%Info,dum,.false.)
   
   end function GetLogLike

    
  function GetLogLikePost(CMB, Info, inCls, HasCls) 
    real GetLogLikePost
    Type (CMBParams) CMB
    Type(ParamSetInfo) Info
    real, intent(in):: inCls(:,:)
    logical, intent(in) :: HasCls
    real acl(lmax,num_cls)
    integer error

    GetLogLikePost  = GetLogPrior(CMB, Info)
    if ( GetLogLikePost >= logZero) then
       GetLogLikePost = logZero
       
    else 

       if (Use_HST) GetLogLikePost = GetLogLikePost + (CMB%H0 - 72)**2/(2*8**2)  !HST 
       if (Use_BBN) GetLogLikePost = GetLogLikePost + (CMB%ombh2 - 0.022)**2/(2*0.002**2) 
          !I'm using increased error bars here
   
       if (Use_CMB .or. Use_LSS) then
          if (HasCls) then
           acl = inCls
           error =0
          else
           call GetCls(CMB, Info, acl, error)
          end if
         if (error /= 0) then
          GetLogLikePost = logZero 
         else
          if (Use_CMB) GetLogLikePost = CMBLnLike(acl) + GetLogLikePost
          if (Use_mpk) GetLogLikePost = GetLogLikePost + LSSLnLike(CMB, Info%theory)
          if (Use_WeakLen) GetLogLikePost = GetLogLikePost + WeakLenLnLike(CMB, Info%theory)     
          if (Use_Lya) GetLogLikePost = GetLogLikePost +  LSS_Lyalike(CMB, Info%Theory)
          if ( GetLogLikePost >= logZero) then
            GetLogLikePost = logZero
          end if
         end if
         if (Use_SN .and. GetLogLikePost /= logZero ) then
            if (Info%Theory%SN_loglike /= 0) then
             GetLogLikePost = GetLogLikePost + Info%Theory%SN_loglike
            else
             GetLogLikePost = GetLogLikePost + SN_LnLike(CMB)
            end if
               !Assume computed only every time hard parameters change
  
  
         end if
   
       else
         if (Use_SN) GetLogLikePost = GetLogLikePost + SN_LnLike(CMB)
       end if

      if (Use_Clusters .and. GetLogLikePost /= LogZero) then
          GetLogLikePost = GetLogLikePost + &
                 (Info%Theory%Sigma_8-0.9)**2/(2*0.05**2)
          stop 'Write your cluster prior in calclike.f90 first!'
      end if

     if (GetLogLikePost /= LogZero) GetLogLikePost = GetLogLikePost/Temperature
   

    end if

  end function GetLogLikePost

end module CalcLike
