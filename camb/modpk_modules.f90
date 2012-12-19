MODULE camb_interface
  INTEGER :: pk_bad
  LOGICAL :: pk_initialized
  LOGICAL :: modpkoutput=.false.
END MODULE camb_interface

MODULE ode_path
  INTEGER*4 :: nok,nbad,kount
  LOGICAL, SAVE :: save_steps=.false.
  LOGICAL :: ode_underflow
  LOGICAL :: ode_ps_output
  LOGICAL :: ode_infl_end
  LOGICAL :: ode_trial_phi=.false. ! SR reconstruction
  REAL*8 :: Nefolds_to_zero=0 ! SR reconstruction
  REAL*8 :: Nefolds_to_kmin=0 ! SR reconstruction
  LOGICAL :: infl_ended
  REAL*8 :: dxsav
  REAL*8, DIMENSION(:), POINTER :: xp
  REAL*8, DIMENSION(:,:), POINTER :: yp
END MODULE ode_path

MODULE modpkparams
  IMPLICIT NONE

  LOGICAL :: use_modpk, modpk_physical_priors, vnderivs, instreheat

! increase max_vparams to use more potential parameters
  INTEGER*4, parameter :: max_vparams = 9

  INTEGER :: vparams_num
  INTEGER :: potential_choice
  REAL*8 :: vparams(max_vparams)
  REAL*8 :: modpk_rho_reheat, modpk_w_primordial_lower, modpk_w_primordial_upper

  INTEGER*4 :: nactual
  INTEGER*4, PARAMETER :: nsteps=1000000
  REAL*8, PARAMETER :: M_Pl=1.0d0
  REAL*8, PARAMETER :: Mpc2Mpl=2.6245d-57
  REAL*8 :: k_pivot, N_pivot
  REAL*8 :: a_init
  REAL*8 :: phi_init0,phi_init,h_init,rescale_factor
  REAL*8 :: phidot_sign
  REAL*8 :: Nefold_max=100000.
  REAL*8 :: lna(nsteps),phiarr(nsteps),dphiarr(nsteps)
  REAL*8 :: hubarr(nsteps),aharr(nsteps)
  LOGICAL :: slowroll_infl_end
  LOGICAL :: slowroll_start=.false.
  REAL*8 :: phi_infl_end=0.
  REAL*8 :: findiffdphi

  REAL*8 :: k_min=1.d-5 ! SR reconstruction
  REAL*8 :: k_max=5 ! SR reconstruction

  REAL*8 :: k_star ! SR reconstruction
  REAL*8 :: H_star, H0 ! SR reconstruction
  LOGICAL :: flag_do_reconstruction=.false. ! SR reconstruction
  REAL*8 :: reconstruction_Nefold_limit ! SR reconstruction

  REAL*8 :: modpk_ns, modpk_nt, modpk_nrun, modpk_As, modpk_r, modpk_w

END MODULE modpkparams


MODULE internals
  IMPLICIT NONE
  REAL, PARAMETER :: PI=3.141592653589793238462643383279502884197
  REAL*8 :: h_ik,pow_ik,powt_ik, phi_ik
  REAL*8 :: k, a_ik
END MODULE internals


MODULE powersp
  IMPLICIT NONE
  INTEGER*4 :: ik
  REAL*8 :: eval_ps,k_start
END MODULE powersp
