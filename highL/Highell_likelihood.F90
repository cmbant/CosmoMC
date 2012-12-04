!===========================================================
MODULE Highell_likelihood
! call high ells experiments

  use highell_options
  use foregrounds_loading
  use act_south_likelihood
  use act_equa_likelihood
  use spt_reichardt_likelihood
  use spt_keisler_likelihood

  contains

  ! ===========================================================================
  SUBROUTINE highell_likelihood_init

   call act_south_likelihood_init
   call act_equa_likelihood_init
   call spt_reichardt_likelihood_init
   call spt_keisler_likelihood_init
   call foregrounds_init

   print *, 'Initializing High ell likelihood'  

 END SUBROUTINE highell_likelihood_init
  ! ===========================================================================


  ! ====================================================================================================================================
  SUBROUTINE highell_likelihood_compute(cl_tt,amp_tsz,amp_ksz,xi,aps148,aps217,aps95,aps150,aps220,acib150,acib220,rps0,rps1,rps,rcib,ags,age,cas1,cas2,cae1,cae2,cal_1,cal_2,cal_3,like_tot)

    IMPLICIT NONE
    REAL(8), dimension(2:tt_lmax_mc) :: cl_tt
    REAL(8), intent(in) :: amp_tsz,amp_ksz,xi,aps148,aps217,aps95,aps150,aps220,acib150,acib220,rps0,rps1,rps,rcib,ags,age
    REAL(8), intent(in) :: cas1,cas2,cae1,cae2,cal_1,cal_2,cal_3
    REAL(8)  :: like_sptr,like_sptk,like_acts,like_acte
    REAL(8), intent(out) :: like_tot

    like_acts = 0.d0
    like_acte = 0.d0
    like_sptr = 0.d0
    like_sptk = 0.d0
    like_tot  = 0.d0


    if (use_act_south .eqv. .true.) then
       call act_south_likelihood_compute(cl_tt,amp_tsz,amp_ksz,xi,aps148,aps217,acib150,acib220,rps,rcib,ags,cas1,cas2,like_acts)
       print *, "----------------------------------------"
       print *, 'ACT south chi2 =', 2*like_acts
    end if

    if (use_act_equa .eqv. .true.) then
       call act_equa_likelihood_compute(cl_tt,amp_tsz,amp_ksz,xi,aps148,aps217,acib150,acib220,rps,rcib,age,cae1,cae2,like_acte)
       print *, 'ACT equa chi2 =', 2*like_acte
    end if

    if (use_spt_highell .eqv. .true.) then
        call spt_reichardt_likelihood_compute(cl_tt,amp_tsz,amp_ksz,xi,aps95,aps150,aps220,acib150,acib220,rps0,rps1,rps,rcib,cal_1,cal_2,cal_3,like_sptr)
        print *, 'SPT high ell chi2 =', 2*like_sptr 
    end if

    if (use_spt_lowell .eqv. .true.) then
       call spt_keisler_likelihood_compute(cl_tt,amp_tsz,amp_ksz,xi,aps150,acib150,cal_2,like_sptk)
       print *, 'SPT low ell chi2 =', 2*like_sptk
    end if

    like_tot= like_acts + like_acte + like_sptr + like_sptk


  END SUBROUTINE highell_likelihood_compute
  ! ====================================================================================================================================


END MODULE Highell_likelihood

