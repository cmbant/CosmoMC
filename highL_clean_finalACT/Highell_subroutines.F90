module Highell_subroutines
!===========================================================================
  contains
  !================================================================================
  subroutine sz_func(freq, sz_corr)

  IMPLICIT NONE 
  REAL(8), intent(in) :: freq
  REAL(8), intent(out) :: sz_corr
  REAL(8), parameter ::  t_cmb = 2.725d0, k_b = 1.3806503d-23, h_pl = 6.626068d-34
  REAL(8) :: nu, xx

  sz_corr = 0.d0
  nu=freq*1.d9
  xx=h_pl*nu/(k_b*t_cmb)
  sz_corr=(2.d0-(xx/2.d0)/tanh(xx/2.d0))
    
  end subroutine

  !================================================================================
    
  subroutine planckfunctionratio(freq,fe,planckfunctionratio_corr)
  ! All the constants cancel out when rescaling the Planck function 
  ! to the effective nominal frequencies
  ! B(nu)/B(nu0)
  IMPLICIT NONE
  REAL(8), intent(in) :: freq, fe
  REAL(8), intent(out) :: planckfunctionratio_corr
  REAL(8), parameter ::  t_eff = 9.7d0
  REAL(8), parameter ::  k_b = 1.3806503d-23, h_pl = 6.626068d-34
  REAL(8) :: nu,nu0,xx,xx0 

  planckfunctionratio_corr = 0.d0
  nu=freq*1.d9
  nu0=fe*1.d9
  xx=h_pl*nu/(k_b*t_eff)
  xx0=h_pl*nu0/(k_b*t_eff)
  planckfunctionratio_corr = (nu/nu0)**3.d0*(exp(xx0)-1.d0)/(exp(xx)-1.d0)
   
  end subroutine
  !=================================================================================

  subroutine flux2tempratio(freq,fe,flux2tempratio_corr)
  ! rescaled to nominal frequencies
  IMPLICIT NONE
  REAL(8), intent(in) :: freq,fe
  REAL(8), intent(out) :: flux2tempratio_corr
  REAL(8), parameter ::  t_cmb = 2.725d0, k_b = 1.3806503d-23, h_pl = 6.626068d-34
  REAL(8) :: nu,nu0,xx,xx0

  flux2tempratio_corr = 0.d0
  nu=freq*1.d9
  nu0=fe*1.d9
  xx=h_pl*nu/(k_b*t_cmb)
  xx0=h_pl*nu0/(k_b*t_cmb)
  flux2tempratio_corr = (nu0/nu)**4.d0*(exp(xx0)/exp(xx))*((exp(xx)-1.d0)/(exp(xx0)-1.d0))**2.d0

  end subroutine
  !=================================================================================

  subroutine get_free_lun( lun )

  implicit none
  integer, intent(out) :: lun
  integer, save :: last_lun = 19
  logical :: used
  lun = last_lun
  
  do
    inquire( unit=lun, opened=used )
    if ( .not. used ) exit
        lun = lun + 1
  end do
  
  last_lun = lun
  end subroutine
  !================================================================================
 
end module Highell_subroutines 
