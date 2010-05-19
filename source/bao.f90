! Percival et al 2009 BAO results hard-coded here by Beth Reid March 2009
! Copied structure from supernovae.f90
!
!default values from http://arxiv.org/abs/0907.1660
!! for explanation of the changes to the rs expression, see Hamann et al, 
!! http://xxx.lanl.gov/abs/1003.3999

module bao
use cmbtypes
use CAMB, only : AngularDiameterDistance  !!angular diam distance also in Mpc no h units
use constants
implicit none

real(dl), parameter :: rstodvz1 = 0.190533, z1 = 0.2
real(dl), parameter :: rstodvz2 = 0.109715, z2 = 0.35
real(dl), dimension(2,2) :: invcov

contains

 subroutine BAO_init

invcov(1,1) = 30124.1d0
invcov(1,2) = -17226.9d0
invcov(2,1) = invcov(1,2)
invcov(2,2) = 86976.6d0

 end subroutine BAO_init

!JH: new routines; integrate to get sound horizon rather than using EH98 formula.
function CMBToBAOrs(CMB)
   use settings
   use cmbtypes
   use ModelParams
   use Precision
   use ThermoData, only : z_drag
   implicit none
   Type(CMBParams) CMB
   real(dl) ::  adrag, atol, rsdrag
   real(dl), external :: dsoundda, rombint
   real(dl) :: CMBToBAOrs
   integer error

   adrag = 1.0d0/(1.0d0+z_drag)
   atol = 1e-6
   rsdrag = rombint(dsoundda,1d-8,adrag,atol)
   CMBToBAOrs = rsdrag

end function CMBToBAOrs


real(dl) function BAO_LnLike(CMB)
  use Precision
  type(CMBParams) CMB
  real :: rs, dv1theory, dv2theory, hz1, hz2, omegam
  real :: rstodvz1theorydelta, rstodvz2theorydelta
!JH: ratio of fitting formula vs. exact result for fiducial model of Percival et al., arXiv:0907.1660
  real, parameter :: rs_rescale = 154.6588d0/150.8192d0
  logical, save :: do_BAO_init = .true.

  if(do_BAO_init) then
    call BAO_init
    do_BAO_init = .false.
  end if

!JH: Need to rescale rs because Percival et al. data assume inaccurate fitting formula result for z_drag
!    rescaled rs has correct dependence on all cosmological parameters though (e.g., N_nu, Y_He, ...)
  rs = CMBToBAOrs(CMB)*rs_rescale

  !!AngularDiameterDistance and rs returned in Mpc no h units.
  !! at z <~ 0.5, the neutrinos are nonrelativistic, so they contribute to the matter density, unlike at zdrag.
  !! note for really tiny neutrino masses, this breaks down; see Section 3.3 of Komatsu et al 2010, WMAP7 cosmological interpretation paper.  However, completely negigible given current error bars!
  omegam = 1.d0 - CMB%omv - CMB%omk
  hz1 = sqrt(omegam*(1.0d0+z1)**3.0d0+CMB%omk*(1.0d0+z1)**2.0+CMB%omv*(1.0d0+z1)**(3.0d0*(1.0d0+CMB%w)))
  hz2 = sqrt(omegam*(1.0d0+z2)**3.0d0+CMB%omk*(1.0d0+z2)**2.0+CMB%omv*(1.0d0+z2)**(3.0d0*(1.0d0+CMB%w)))
  dv1theory = ((1.0d0+z1)*AngularDiameterDistance(z1))**2.0d0*c*z1/CMB%H0/hz1/1000.0d0
  dv2theory = ((1.0d0+z2)*AngularDiameterDistance(z2))**2.0d0*c*z2/CMB%H0/hz2/1000.0d0
  dv1theory = dv1theory**(1.0d0/3.0d0)
  dv2theory = dv2theory**(1.0d0/3.0d0)

  rstodvz1theorydelta = rs/dv1theory - rstodvz1
  rstodvz2theorydelta = rs/dv2theory - rstodvz2

  BAO_LnLike = 0.5*((rstodvz1theorydelta) * invcov(1,1) * (rstodvz1theorydelta) &
     & + 2.0d0 * (rstodvz1theorydelta) * invcov(1,2) * (rstodvz2theorydelta) &
     & + (rstodvz2theorydelta) * invcov(2,2) * (rstodvz2theorydelta))

   if (Feedback > 1) print *,'BAO_LnLike: ',BAO_LnLike
end function BAO_LnLike

end module bao
