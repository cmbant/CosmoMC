! Percival et al 2009 BAO results hard-coded here by Beth Reid March 2009
! Copied structure from supernovae.f90
!
!default values from http://arxiv.org/abs/0907.1660

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

function CMBToBAOrs(CMB)
     use settings
     use cmbtypes
     use ModelParams
     use Precision
     implicit none
     Type(CMBParams) CMB
     double precision zdrag, adrag, atol, rsdrag, myomh2val, b1, b2
     double precision, external :: dsoundda, rombint
     real CMBToBAOrs
     integer error

!!From Eisenstein and Hu 1998, formula for zdrag.
!!At zdrag neutrinos are still relativistic, so omega_m h^2 is omch2 + ombh2 (does not include neutrinos).
!!Note be careful if you muck around here -- I had to add parenthesis and/or d0 to get these numbers to come out correctly.
       myomh2val = (CMB%omdmh2*(1.0-CMB%nufrac)+CMB%ombh2)
       b1 = 0.313d0*(myomh2val)**(-0.419d0)*(1.0d0+0.607*(myomh2val)**(0.674d0))
       b2 = 0.238d0*(myomh2val)**0.223d0
       zdrag =  1291.0d0*(myomh2val)**0.251d0/(1.0d0+0.659*(myomh2val)**0.828)*(1.0d0+b1*(CMB%ombh2)**b2)
       adrag = 1.0d0/(1.0d0+zdrag)
       atol = 1e-6
       rsdrag = rombint(dsoundda,1d-8,adrag,atol)
       CMBToBAOrs = rsdrag
  end function CMBToBAOrs

real(dl) function BAO_LnLike(CMB)
  use Precision
  type(CMBParams) CMB
  real :: rs, dv1theory, dv2theory, hz1, hz2, omegam
  real :: rstodvz1theorydelta, rstodvz2theorydelta
  logical, save :: do_BAO_init = .true.

  if(do_BAO_init) then
    call BAO_init
    do_BAO_init = .false.
  end if

  rs = CMBToBAOrs(CMB)

  !!AngularDiameterDistance and rs returned in Mpc no h units.
  !! at z <~ 0.5, the neutrinos are nonrelativistic, so they contribute to the matter density, unlike at zdrag.
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
