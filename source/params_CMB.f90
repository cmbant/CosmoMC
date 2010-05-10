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
     implicit none
     Type(CMBParams) CMB
     double precision zstar, astar, atol, rs, dsoundda, DA, rombint
     external dsoundda, rombint
     real CMBToTheta
     integer error

     call InitCAMB(CMB,error)

!!From Hu & Sugiyama
       zstar =  1048*(1+0.00124*CMB%ombh2**(-0.738))*(1+ &
        (0.0783*CMB%ombh2**(-0.238)/(1+39.5*CMB%ombh2**0.763)) * &
           (CMB%omdmh2+CMB%ombh2)**(0.560/(1+21.1*CMB%ombh2**1.81)))
     
       astar = 1/(1+zstar)
       atol = 1e-6
       rs = rombint(dsoundda,1d-8,astar,atol)
       DA = AngularDiameterDistance(zstar)/astar
       CMBToTheta = rs/DA
!       print *,'z* = ',zstar, 'r_s = ',rs, 'DA = ',DA, rs/DA

  end function CMBToTheta


      function TauToZre(tau)
      use ModelParams
      implicit none
      real TauToZre,tau
! This subroutine calculates the redshift of reionization
! from an Reion%optical_depth and sets the reionization fraction CP%Reion%fraction=1.
! assumes sharp reionization and fraction of 1   
      integer na
      real dz, optd
      real(dl) dtauda  !diff of tau w.r.t a and integration
      external dtauda

      real(dl) z
      
         na=1
         optd=0
         dz=1._dl/2000
         z=0
         do while (optd < tau)
            z=na*dz
            optd=optd+dz*akthom*dtauda(1/(1+z))
            na=na+1
         end do
         TauToZre = z
   
      end function TauToZre


!Mapping between array of power spectrum parameters and CAMB
     subroutine SetCAMBInitPower(P,CMB,in)
       use camb
       use settings
       use cmbtypes
       implicit none
       type(CAMBParams)  P
       Type(CMBParams) CMB

       integer, intent(in) :: in


       if (Power_Name == 'power_tilt') then

       P%InitPower%k_0_scalar = pivot_k
       P%InitPower%k_0_tensor = pivot_k

       P%InitPower%ScalarPowerAmp(in) = cl_norm*CMB%norm(norm_As)
       P%InitPower%rat(in) = CMB%norm(norm_amp_ratio)
        
       P%InitPower%an(in) = CMB%InitPower(1)
       P%InitPower%ant(in) = CMB%InitPower(2)
       P%InitPower%n_run(in) = CMB%InitPower(3)
       if (inflation_consistency) then
         P%InitPower%ant(in) = - CMB%norm(norm_amp_ratio)/8.
          !note input n_T is ignored, so should be fixed (to anything)
       end if
       else
         stop 'params_CMB:Wrong initial power spectrum'
       end if

    end subroutine SetCAMBInitPower
 

 subroutine SetForH(Params,CMB,H0)
     use settings
     use cmbtypes
     use CMB_Cls
     implicit none
     real Params(num_Params)

     Type(CMBParams) CMB
     real h2,H0
  
    CMB%H0=H0
    CMB%ombh2 = Params(1)    
    CMB%omdmh2 = Params(2)
    CMB%zre = Params(4) !!Not actually used.. is tau in this parameterization
    CMB%Omk = Params(5)
    CMB%nufrac = Params(6)
    CMB%w = Params(7)

    CMB%InitPower(1:num_initpower) = Params(index_initpower:index_initpower+num_initPower-1)
    CMB%norm(1) = exp(Params(index_norm))
    CMB%norm(2:num_norm) = Params(index_norm+1:index_norm+num_norm-1)

    CMB%h = CMB%H0/100
    h2 = CMB%h**2
    CMB%omnuh2 = CMB%omdmh2*CMB%nufrac
    CMB%omch2 = CMB%omdmh2 - CMB%omnuh2
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
!SZ  use WMAP_OPTIONS

     implicit none
     real Params(num_params)
     real, save :: LastParams(num_params) = 0.
     real, save :: LastH0, Lastzre

     Type(CMBParams) CMB
     real DA, tauToZre
     integer error
     real  D_b,D_t,D_try,try_b,try_t, CMBToTheta, lasttry,tau
     external TauToZre,CMBToTheta

!SZ   SZ_amp= Params(13)

     if (all(Params(1:num_hard) == Lastparams(1:num_hard))) then
       call SetForH(Params,CMB,LastH0)
       CMB%zre = Lastzre
       CMB%reserved(1) = params(4)
     else

     DA = Params(3)/100
     try_b = 40
     call SetForH(Params,CMB,try_b)
     D_b = CMBToTheta(CMB)
     try_t = 100
     call SetForH(Params,CMB,try_t)
     D_t = CMBToTheta(CMB)
     if (DA < D_b .or. DA > D_t) then
      cmb%H0=0 !Reject it
     else
     lasttry = -1
     do
            call SetForH(Params,CMB,(try_b+try_t)/2)
            D_try = CMBToTheta(CMB)
               if (D_try < DA) then
                  try_b = (try_b+try_t)/2
               else
                  try_t = (try_b+try_t)/2
               end if
               if (abs(D_try - lasttry)< 1e-7) exit
              lasttry = D_try
     end do

    call InitCAMB(CMB,error)
    tau = params(4)
    CMB%zre = TauToZre(tau)       

    LastH0 = CMB%H0
    Lastzre = CMB%zre
    LastParams = Params
    end if

    CMB%reserved = 0
    CMB%reserved(1) = params(4) !tau
  
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
 
      Params(1) =CMB%ombh2 
      Params(2) =CMB%omdmh2
 
      Params(3) = CMBToTheta(CMB)*100
      Params(4) = CMB%reserved(1)
      Params(5) =CMB%omk 
      
      Params(6) =CMB%nufrac 
      Params(7) =CMB%w
      Params(index_initpower:index_initpower+num_initpower-1) =CMB%InitPower(1:num_initpower) 
      Params(index_norm) = log(CMB%norm(1))
      Params(index_norm+1:index_norm+num_norm-1) = CMB%norm(2:num_norm)
  
   end subroutine CMBParamsToParams


  subroutine WriteParams(P, mult, like)
     use settings
     use cmbtypes
     use ParamDef
     implicit none
    Type(ParamSet) P
    real, intent(in) :: mult, like
    character(LEN =30) fmt
    Type(CMBParams) C
    real r10
  
    if (outfile_unit ==0) return
      call ParamsToCMBParams(P%P,C)

      if (lmax_tensor /= 0 .and. compute_tensors) then
          r10 = P%Info%Theory%cl_tensor(10,1)/P%Info%Theory%cl(10,1)
      else
        r10 = 0
      end if

      fmt = trim(numcat('(2E16.7,',num_params))//'E16.7,7E16.7)'
      write (outfile_unit,fmt) mult,like, P%P, C%omv,P%Info%Theory%Age, C%omdm+C%omb, &
          P%Info%Theory%Sigma_8, C%zre,r10,C%H0

     if (flush_write) call FlushFile(outfile_unit)

  end  subroutine WriteParams




  subroutine WriteParamsAndDat(P, mult, like)
     use settings
     use cmbtypes
     use ParamDef
     implicit none
    Type(ParamSet) P
    real, intent(in) :: mult, like
    character(LEN =30) fmt
    Type(CMBParams) C
    real r10
  
    if (outfile_unit ==0) return
      call ParamsToCMBParams(P%P,C)

      if (lmax_tensor /= 0 .and. compute_tensors) then
          r10 = P%Info%Theory%cl_tensor(10,1)/P%Info%Theory%cl(10,1)
      else
        r10 = 0
      end if

     fmt = trim(numcat('(2E16.7,',num_params+num_matter_power))//'E16.7,7E16.7)'
      write (outfile_unit,fmt) mult,like, P%P, C%omv,P%Info%Theory%Age, C%omdm+C%omb, &
          P%Info%Theory%Sigma_8, C%zre,r10,C%H0, P%Info%Theory%matter_power(:,1)


     if (flush_write) call FlushFile(outfile_unit)

  end  subroutine WriteParamsAndDat


  function dsoundda(a)
          use Precision
          use ModelParams
     
          implicit none
          real(dl) dsoundda,dtauda,a,R,cs
          external dtauda

           R=3.0d4*a*CP%omegab*(CP%h0/100.0d0)**2
           cs=1.0d0/sqrt(3*(1+R))
           dsoundda=dtauda(a)*cs
        
  end function dsoundda

