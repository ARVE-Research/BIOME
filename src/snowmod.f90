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

! parameters

real(sp), parameter :: k7  = -0.5827
real(sp), parameter :: k8  =  1.319
real(sp), parameter :: k9  =  4.18e-4
real(sp), parameter :: k10 =  1.140e-2

! local variables

integer :: d

real(sp), dimension(size(Tmax)) :: psnow  ! probability of snowfall occurrence

! ----

do d = 1,size(Tmax)

  psnow(d) = 1. / (1. + exp(k7 + k8 * Tmax(d) + k9 * z + k10 * phi))

end do

Tt = maxloc(Tmax,dim=1,mask=psnow >= 0.5)

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
