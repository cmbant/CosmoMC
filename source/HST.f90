!HST from http://arxiv.org/abs/0905.0695
!Thanks to Beth Reid, minor mods by AL, Oct 09

!Updated by KM to use Riess et al (2011) value of H0 = 73.8 +/- 2.4 km/s/Mpc
!Riess et al: 1103.2976

module HST
use cmbtypes
use CAMB, only : AngularDiameterDistance  !!physical angular diam distance also in Mpc no h units
use constants
use likelihood
implicit none

type, extends(DataLikelihood) :: HSTLikelihood
    contains
    procedure :: LogLike => HST_LnLike
end type HSTLikelihood

! angdistinveffh0 is the inverse of the angular diameter distance at z = 0.04 for H_0 = 74.2
! and a fiducial cosmology (omega_k = 0, omega_lambda = 0.7, w = -1); this is proportional to 
!H_0 but includes the tiny cosmological dependence of the measurement (primarily on w) correctly.
! angdistinveffh0err = 3.6 / DL(0.04)
!real(dl), parameter :: angdistinvzeffh0 = 6.49405e-3, zeffh0 = 0.04, &
!                       angdistinvzeffh0errsqr = 9.93e-8

real(dl), parameter :: angdistinvzeffh0 = 6.45904e-3, zeffh0 = 0.04, &
                       angdistinvzeffh0errsqr = 4.412e-8

contains

   subroutine HST_init(like, ini)
    use IniFile
    class(HSTLikelihood) :: like
    Type(TIniFile) :: ini
    
    Like%LikelihoodName = 'Hubble'
    call Like%datasets%Add(DataItem('HST'))
    
    end subroutine HST_init
    
real function HST_LnLike(like, CMB, Theory)
  Class(HSTLikelihood) :: like
  type(CMBParams) CMB
  Type(CosmoTheory) Theory
  real(dl) :: theoryval

  theoryval = 1.0/AngularDiameterDistance(zeffh0)
  HST_LnLike = (theoryval - angdistinvzeffh0)**2/(2*angdistinvzeffh0errsqr)
  if (Feedback > 1) print *,'HST_LnLike like: ',HST_LnLike
end function  HST_LnLike

end module HST
