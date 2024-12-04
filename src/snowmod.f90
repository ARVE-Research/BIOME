module snowmod

implicit none

public :: Tt

contains

! ---------------------------------------

real(sp) function Tt(Tmax,z,phi)

use parametersmod, only : sp,dp

implicit none

! arguments

real(sp), dimension(:), intent(in) :: Tmax  ! vector of smoothed pseudo-daily maximum air temperature (degC)
real(sp),               intent(in) :: z     ! elevation (m)
real(dp),               intent(in) :: phi   ! latitude (degrees)

! parameters                  line 576      | from paper
!                            splash.point.R | (NB the sign of k9 and k10 is inconsistent between paper and code)
real(sp), parameter :: k7  = -0.4710405934  ! -0.5827
real(sp), parameter :: k8  =  1.0473543991  !  1.319
real(sp), parameter :: k9  = -0.0004596581  ! -4.18e-4
real(sp), parameter :: k10 = -0.0110592101  ! -1.140e-2

! 1/(1+exp(-0.4710405934+1.0473543991*as.numeric(tc)-as.numeric(elev)*0.0004596581-abs(as.numeric(lat))*0.0110592101))

! psnow(d) = 1. / (1. + exp(k7 + k8 * Tmax(d) + k9 * z + k10 * abs(phi)))

! local variables

! integer :: d

real(sp), dimension(size(Tmax)) :: psnow  ! probability of snowfall occurrence

! ----
! analytical solution of Sandoval et al. (2024) eqn 19 for Tt at psnow of 0.5

Tt = (-k7 - k9 * z - k10 * abs(phi)) / k8

end function Tt

! ---------------------------------------


! real(sp) function frain(Tair,phi,z)
! 
! use parametersmod, only : sp
! 
! ! arguments
! 
! real(sp), intent(in) :: Tair  ! air temperature (degC)
! real(sp), intent(in) :: z     ! elevation (m)
! real(sp), intent(in) :: phi   ! latitude (degrees)
! 
! 
! real(sp) :: Tt  ! annual maximum air temperature when the probability of snowfall occurrence meets or exceeds 0.5
! 
! 
! 
! real(sp), parameter :: Tr  = 13.
! 
! Trm = Tt + Tt * sin((mon + 2.) / 1.91)

! ---------------------------------------

end module snowmod
