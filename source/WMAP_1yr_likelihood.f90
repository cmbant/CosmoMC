! WMAP likelihood code
! Written by Licia Verde and Hiranya Peiris, Princeton University,
!      December 2002.
!
! F90 version by Antony Lewis Feb 03

!****************************************************************************
! If you use this code in any publication, please cite Verde et al. (2003),
! Hinshaw et al. (2003) and Kogut et al. (2003). 
!****************************************************************************


! History
! 28 Feb 03: AL
!   * Fixed reading in of lmax'th element of TE array
!   * Uses packed array storage (1/2 memory use)
!   * Switch to functions, consistent naming scheme
! 27 Feb 03: AL 
!   * Change to add in off-diagonal TT Fisher only if the model is not a very bad fit
!      -- prevents erroneous high likelihoods being returned for very bad fit models 
!   * F90 translation, misc changes, deletion of separate Fisher routines
!   * Re-organized to use symmetries (2*speed up) 
! 21 Feb 03: WMAP
!   * Fixed bug reading in off-diagonmal TE matrix.
!     -- changed likelihoods slightly, no significant effect on parameter estimates
!     (importance weights of order 3 relative to buggy version)

! This file contains a set of Fortran subroutines that compute likelihoods
! for a set of theoretical cl's:
!    WMAP_init_TT      - Loads TT data files
!                               Call this routine once before calling
!                               WMAP_LnLike_TT
!    WMAP_init_TE      - Loads TE data files.
!                               Call this routine once before calling
!                               WMAP_LnLike_TE
!    WMAP_init         - Calls WMAP_init_TT and WMAP_init_TE.
!    WMAP_LnLike_TT    - Computes temperature data likelihoods.
!    WMAP_LnLike_TE    - Computes TE data likelihoods.
!
! The methods used here are described in the following references:
!    Verde, L., et.al. 2003, ApJ, submitted. astro-ph/0302218
!    Bond, J.R., Jaffe, A.H., and Knox, L., 2002, ApJ, 533
!
! ===========================================================================
!   "WMAP_init_TT" fills the common block used by
!   "WMAP_LnLike_TT"---this is the initialization subroutine.
!
!   Inputs:
!      clFile  - The name of the file containing the cl data.
!      offDiag - The name of the file containing the off-diagonal terms.
!
!   Outputs:
!      stat - A status code:  0=success.
!
!   Written by Licia Verde and Hiranya Peiris, Princeton University,
!      December 2002.
! ===========================================================================

    module WMAP
     implicit none
      private

      integer, parameter :: WMAP_lmax_TT = 900, WMAP_lmax_TE = 450, &
             WMAP_lmax_TE_file = 512 
    
      integer, parameter :: WMAP_lmin_TT = 2, WMAP_lmin_TE = 2
     !Can change l_min here to remove e.g. quadrupole 

      integer, parameter :: wp = KIND(1.0)
      integer, parameter :: WMAP_precision =wp
!TT data
      real(wp) cl_data(WMAP_lmax_TT), neff(WMAP_lmax_TT), fskyeff(WMAP_lmax_TT)
      real(wp) r_off_diag(((WMAP_lmax_TT-1)*(WMAP_lmax_TT-2))/2)
      real(wp) off_diag(((WMAP_lmax_TT-1)*(WMAP_lmax_TT-2))/2)

!TE data
      real(wp) te_data(WMAP_lmax_TE),ntt(WMAP_lmax_TE),nee(WMAP_lmax_TE)
      real(wp) te_off_diag(((WMAP_lmax_TE-1)*(WMAP_lmax_TE-2))/2)
      real(wp) te_tt(WMAP_lmax_TE)

!Public constant and subroutines 
      public  WMAP_lmax_TT, WMAP_lmax_TE, WMAP_precision,&
              WMAP_LnLike_TT,WMAP_LnLike_TE, &
              WMAP_init_TT, WMAP_init_TE,  WMAP_init             
    contains

        SUBROUTINE WMAP_init_TT (clFile, offDiag, stat)
!
        IMPLICIT NONE
!                INPUT
!                
        character (*) clFile, offDiag
!
!                OUTPUT
!
        integer stat
!
!                DATA/INTERNAL
!
        integer idum, l, ll,i,j,ix
! ---------------------------------------------------------------------------
!               Read the CL data.
!
        open (11, file=clFile, status='old', IOStat=stat)
        if (stat .NE. 0) RETURN
        rewind(11)
        do l = 2, WMAP_lmax_TT
          read (11, *) idum, cl_data(l), neff(l), fskyeff(l)
          if (idum /= l) stop 'Error reading TT diag'
        end do
        close (11)
!
!               Read in off diag terms:               
!               In the *covariance* matrix there are 2 type of terms 
!               one that does scale with the clth (due to the mask) and 
!               one that does not (due to beam and point sources marginalization).
!               In the curvature matrix (what we use here) all the off diagonal 
!               terms end up scaling with the clth but the 2 contributions scale
!               in different ways see paper for details
!               thus here we read them in separately
!
        open(11,file=offDiag, status='old', IOStat=stat)
        if (stat .NE. 0) RETURN
        rewind(11)
        ix=1 
        do  l = 2,WMAP_lmax_TT  
           do ll= l+1,WMAP_lmax_TT
              read(11,*) i,j,off_diag(ix),r_off_diag(ix)  
              if (l >= WMAP_lmin_TT) ix=ix+1 !Skip all entries up to lmin            
              if (l.ne.i .or. ll .ne. j) stop 'error reading TT off diag' 
           end do
        end do
        close(11)
! ---------------------------------------------------------------------------
        END SUBROUTINE WMAP_init_TT
! ===========================================================================
!   "WMAP_init_TE" fills the common block used by
!   "compute_mapte_likelihood"---this is the initialization subroutine.
!
!   Inputs:
!      clFile  - The name of the file containing the cl data.
!      offDiag - The name of the file containing the off-diagonal terms.
!
!   Outputs:
!      stat - A status code:  0=success.
!
!   Common blocks:
!        te_dat  - This common block is used to pass/store information
!                  that does not change between calls.
!
!   Written by Licia Verde and Hiranya Peiris, Princeton University,
!      December 2002.
! ===========================================================================
        SUBROUTINE WMAP_init_TE (clFile, offDiag, stat)
!                INPUT
        character (*) clFile, offDiag
!
!                OUTPUT
!
        integer stat
!
!                DATA/INTERNAL
!
        integer idum, l, ll, i,j,ix
        real(wp) tmp
! ---------------------------------------------------------------------------
!               Read the CL data.
!
        if (WMAP_lmax_TE_file < WMAP_lmax_TE) stop 'Wrong WMAP_lmax_TE_file'

        open (11, file=clFile, status='old', IOStat=stat)
        if (stat .NE. 0) RETURN
        rewind(11)
        do  l = 2, WMAP_lmax_TE
           read(11,*) idum,te_data(l),te_tt(l),ntt(l),nee(l)
           if (idum.ne.l) stop 'Error reading TE diag file'
        end do
        close(11)
!
!               Read in off diag terms.
!
        te_off_diag = 0
        open(11,file=offDiag, status='old', IOStat=stat)
        if (stat .NE. 0) RETURN
        rewind(11)
        ix=1
        do  l = 2,WMAP_lmax_TE  
           do ll= l+1,WMAP_lmax_TE_file 
              read(11,*) i,j,tmp
              if (l /= i .or. j /= ll) stop 'Error reading TE file'
              if (l>=WMAP_lmin_TE .and. ll<=WMAP_lmax_TE) then
               te_off_diag(ix) = tmp
               ix=ix+1
              end if
           end do
        end do
        close(11)
! ---------------------------------------------------------------------------
        END SUBROUTINE WMAP_init_TE
! ===========================================================================
!   "WMAP_init" encapsulates "WMAP_init_TT" and
!   "WMAP_init_TE" into a single call.
!
!   Inputs:
!      TclFile   - The name of the file containing the temp. cl data.
!      ToffDiag  - The name of the file containing the temp. off-diagonal terms.
!      TEclFile  - The name of the file containing the TE cl data.
!      TEoffDiag - The name of the file containing the TE off-diagonal terms.
!
!   Outputs:
!      stat - A status code:  0=success.
!
!   Written by Licia Verde and Hiranya Peiris, Princeton University,
!      December 2002.
! ===========================================================================
        SUBROUTINE WMAP_init (TclFile, ToffDiag,&
                                          TEclFile, TEoffDiag, stat)
!                INPUT
!                
        character (*) TclFile, ToffDiag, TEclFile, TEoffDiag
!
!                OUTPUT
!
        integer stat
! ---------------------------------------------------------------------------
        Call WMAP_init_TT (TclFile, ToffDiag, stat)
        If (stat .NE. 0) Return
        Call WMAP_init_TE (TEclFile, TEoffDiag, stat)
! ---------------------------------------------------------------------------
        END SUBROUTINE WMAP_init
! ===========================================================================
!   "WMAP_LnLike_TT" computes the likelihood for temperature data
!    using the form for the likelihood as in Verde et al 2003 sec. 2.1 
!
!    There are 2 contributions to the  off diagonal terms of the fisher matrix: 
!    this is because in the covariance matrix the terms due to marginalization over 
!    point sources and beam uncertainties depend on the power spectrum that's out 
!    there in the sky while the coupling due to the mask depend on the Cl_theory 
!    (i.e. changes as the cosmological parameters change in exploring the 
!    likelihood surface) see Verde et al. 2003 , Hinshaw et al 2003.
!
!    WARNING: a 2% bias around the first peak results from incorrect scaling 
!    of the off diagonal terms.    
!
!   Inputs:
!      clth - An array of cl's describing the theory being tested.
!
!   Outputs:
!      LnLike   - The natural log of the likelihood for l >= 2.
!
!   Common blocks:
!        tt_data - This common block is used to pass/store information
!                  that does not change between calls.
!
!   Written by Licia Verde and Hiranya Peiris, Princeton University,
!      December 2002.
! ===========================================================================
        FUNCTION WMAP_LnLike_TT(clth)
!
        IMPLICIT NONE
!
!                INPUT: clt_theory: l(l+1)C_l/2pi in microK^2
!                
        real(wp) clth(*)
!
!                OUTPUT: log(likelihood)
!
        real(wp) WMAP_LnLike_TT
!
!                DATA/INTERNAL
!
        integer l,ll,ix
        real(wp) chisq, offchisq
        real(wp) dchisq, Fisher, Fdiag(WMAP_lmax_TT)
        real(wp) z(WMAP_lmax_TT), zbar(WMAP_lmax_TT)
        real(wp) off_log_curv
        real(wp) Fdiagsqrt(WMAP_lmax_TT)
! ---------------------------------------------------------------------------
        chisq   = 0
                
! prepare to compute the offset lognormal likelihood as in Bond Jaffe Knox.
! with the difference that the transformation we do on the curvature
! matrix is using clth not cltdata
! this is closer to the equal variance approx (see Bond Jaffe Knox again)  
! form more details see Verde et al .2003 
!
! here Fdiag is the diagonal term of the covariance matrix 
! Fisher denotes the curvature matrix.
!
        do l = WMAP_lmin_TT, WMAP_lmax_TT        
           Fdiag(l) = 2*(Clth(l)+neff(l))**2&
                    / ((2*l+1)*fskyeff(l)**2)
           Fdiagsqrt(l) = 1/sqrt(Fdiag(l))
           z(l)=log(cl_data(l)+neff(l))           
           zbar(l)=log(clth(l)+neff(l))
 
! Get the diagonal terms in the likelihood
           Fisher=1/Fdiag(l)
           off_log_curv=(clth(l)+neff(l))**2*Fisher

           dchisq=2._wp/3*(z(l)-zbar(l))**2*off_log_curv &
                      +1._wp/3*(clth(l)-cl_data(l))**2*Fisher 
           chisq=chisq+dchisq
 
        end do
        
        if (chisq < WMAP_lmax_TT*2) then
!Only get off-diagonal terms if not a really bad fit, otherwise they will 
!be wildly wrong
!In principle the likelihood may not be continuous, but likelihood is so low
!model will always be rejected or contribute zero anyway

        offchisq=0
        ix =1
        do l = WMAP_lmin_TT,WMAP_lmax_TT        
           do ll=l+1,WMAP_lmax_TT
!
! here the two contributions to the off diagonal terms to the curvature matrix 
! are treated separately
!
! see Verde et.al. 2003 for details.
! note off_diag = -epsilon in the papers
            Fisher=r_off_diag(ix)*Fdiagsqrt(l)*Fdiagsqrt(ll) &
                             +off_diag(ix)/(Fdiag(l)*Fdiag(ll))
            off_log_curv=(clth(l)+neff(l))*Fisher*(clth(ll)+neff(ll))

! to correct for residual 0.5% bias around the peak
! see Verde et.al. 2003 for more details
!
! this is an interpolation between Bond Knox Jaffe and Gaussian Likelihood.
! works extremely well on sims (again see paper for details)
!
              dchisq=2._wp/3*(z(l)-zbar(l))*off_log_curv*(z(ll)-zbar(ll))+ &
                    1._wp/3*(clth(l)-cl_data(l))*Fisher*(clth(ll)-cl_data(ll))
      
              offchisq=offchisq+dchisq
              ix=ix+1           
           end do
        end do

! add it twice because of the symmetry of the matrix
        chisq=chisq+offchisq*2
        
        end if 

        WMAP_LnLike_TT=-chisq/2.d0
!
!        write(*,*) ' MAP_T log(LnLike)=',LnLike,'chisq=',chisq
!
! ---------------------------------------------------------------------------
        END FUNCTION WMAP_LnLike_TT


! "WMAP_LnLike_TE" computes the likelihood for TE data.
!
! Uses the expression for the likelihood  as in Kogut et al 2003 (ApJ in press) 
! and Verde et al 2003.
!
! The Curvature matrix has been calibrated from Monte Carlo simulations.  
! fsky is set to 0.85, if not using the P2 cut (i.e. if not using the cl 
! we provided you with) you should change fsky accordingly.
!
!   Inputs:
!      cltt - An array of temperature cl's describing the theory
!             being tested.
!      clte - An array of TE cl's describing the theory
!             being tested.
!      clee - An array of EE cl's describing the theory
!             being tested.
!
!   Outputs:
!      LnLike   - The natural log of the likelihood for all l.
!
!   Common blocks:
!        te_dat - This common block is used to pass/store information
!                 that does not change between calls.
!
!   Written by Licia Verde and Hiranya Peiris, Princeton University,
!      December 2002.
! ===========================================================================
        FUNCTION WMAP_LnLike_TE(cltt, clte, clee)
!
        IMPLICIT NONE

!               INPUT: clte_theory
!                
        real(wp) cltt(*),clte(*),clee(*)
!
!               OUTPUT: log(likelihood)
!
        real(wp) WMAP_LnLike_TE

        integer il,ill,ix
        real(wp) fsky,delta_chisq,chisq, offchisq
        real(wp) Fdiag, Fdiagsqrt(WMAP_lmax_TE),Fisher
      !  real(wp) ct, ce
! ---------------------------------------------------------------------------
        fsky=0.85d0
        chisq=0

!Alternative code that includes TT - TE correlations
!        do il = WMAP_lmin_TE,WMAP_lmax_TE        
!            ct = cltt(il) + ntt(il)
!            ce = clee(il) + nee(il)
!            chisq = chisq+ (te_data(il) - clte(il) - clte(il)/ct*(te_tt(il)-ct))**2/ &
!                    (ct*ce-clte(il)**2)  * (2*il+1)*fsky**2/1.14d0      
!
!        end do
!       WMAP_LnLike_TE=-chisq/2.d0
!       return


! Fdiag is the diagonal element of the covariance matrix
! Fisher denotes the curvature matrix.
!
        do il = WMAP_lmin_TE,WMAP_lmax_TE        
           Fdiag=((cltt(il)+ntt(il))*(clee(il)+nee(il))&
                    +clte(il)*clte(il))&
                    /((2*il+1)*fsky**2/1.14d0)
           Fdiagsqrt(il) = 1/sqrt(Fdiag)

! this correction factor (1.14) has been obtained from a calibration of
! the covariance matrix on monte carlo sims (see papers for more details)
!
           delta_chisq = (clte(il)-te_data(il))**2/Fdiag 

           chisq=chisq+delta_chisq
          end do

!
! include mask effect on off diagonal terms in the curvature matrix
   
        offchisq = 0
        ix = 1
        do il = WMAP_lmin_TE,WMAP_lmax_TE
           do ill=il+1,WMAP_lmax_TE
              Fisher=te_off_diag(ix)*Fdiagsqrt(il)*Fdiagsqrt(ill)
 
              delta_chisq = (clte(il)-te_data(il))*Fisher*&
                              (clte(ill)-te_data(ill))
              offchisq=offchisq+delta_chisq
              ix=ix+1 
           end do
        end do

        chisq = chisq+offchisq*2

        WMAP_LnLike_TE=-chisq/2.d0
!
!        write(*,*) ' MAP_TE log(LnLike)=',LnLike,'chisq=',chisq
!
! ---------------------------------------------------------------------------
        END FUNCTION WMAP_LnLike_TE

  end module WMAP
  
