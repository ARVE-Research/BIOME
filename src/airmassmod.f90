module airmassmod

use parametersmod, only : sp

implicit none

public  :: initairmass
public  :: Pjj
public  :: elev_corr
public  :: airmass
private :: m
private :: F

! module calculated parameters

real(sp), dimension(3) :: c00  ! air mass coefficients for solar zenith angle <= 80 degrees
real(sp), dimension(3) :: c80  ! air mass coefficients for solar zenith angle  > 80 degrees

contains

! ----------------------------------------------------------------------------------------------------------------

subroutine initairmass()

! calculate parameters used in the airmass calculations

use parametersmod, only : sp,pir => pir_sp

implicit none

! parameters

real(sp), parameter :: m0  =  1.   ! air mass at  0 degree solar zenith angle
real(sp), parameter :: m80 =  5.6  ! air mass at 80 degree solar zenith angle
real(sp), parameter :: m90 = 39.7  ! air mass at 90 degree solar zenith angle

real(sp), parameter :: cos80 = cos(80. * pir)  ! (degrees)

! ---

c00(1) = 0.008307
c00(2) = (m0 - m80) * (c00(1) + 1.) * (c00(1) + cos80) / (cos80 - 1.)
c00(3) = m0 - c00(2) / (c00(1) + 1.)

c80(1) = 0.037160
c80(2) = (m90 - m80) * c80(1) * (c80(1) + cos80) / cos80
c80(3) = m90 - c80(2) / c80(1)

end subroutine initairmass

! ----------------------------------------------------------------------------------------------------------------

real(sp) function Pjj(mtemp,mprec)

! calculation of the "precipitation equitability index"

use parametersmod, only : sp
use utilitymod,    only : imaxloc,iminloc

implicit none

! arguments

real(sp), dimension(:), intent(in) :: mtemp  ! annual time series of monthly temperature
real(sp), dimension(:), intent(in) :: mprec  ! annual time series of monthly precipitation

! local variables

integer :: wm  ! index position of the warmest month
integer :: cm  ! index position of the warmest month

real(sp) :: p_wm  ! precipitation in the warmest month 
real(sp) :: p_cm  ! precipitation in the coldest month

! ------

wm   = imaxloc(mtemp)
cm   = iminloc(mtemp)
p_wm = mprec(wm)
p_cm = mprec(cm)   

if (p_wm + p_cm > 0.) then

  Pjj = 2. * (p_wm - p_cm) / (p_wm + p_cm)
  Pjj = max(Pjj,0.)

else

  Pjj = 0.

end if

end function Pjj

! ----------------------------------------------------------------------------------------------------------------

real(sp) function elev_corr(elevation)

implicit none

real(sp), intent(in)  :: elevation 

real(sp), parameter :: z0 = 1. / 8000.

! ----

elev_corr = exp(-elevation * z0)

end function elev_corr

! ----------------------------------------------------------------------------------------------------------------

real(sp) function m(cosZ,c)

! Instantaneous air mass m, equation 2.1 in Yin, 1997

implicit none

real(sp),               intent(in) :: cosZ
real(sp), dimension(:), intent(in) :: c

m = c(2) / (c(1) + cosZ) + c(3)

end function m

! ----------------------------------------------------------------------------------------------------------------

real(sp) function F(t1,a,b,c)

! integral air mass function F, equation 2.6b in Yin, 1997
! section inside curly braces only - multiply result by 1/t1 to get mbar

use parametersmod, only : sp,pi => pi_sp,pir => pir_sp

implicit none

! arguments

real(sp),               intent(in) :: t1
real(sp),               intent(in) :: a
real(sp),               intent(in) :: b
real(sp), dimension(:), intent(in) :: c

! parameters

real(sp), parameter :: w   = 15.      ! solar angular velocity (degrees hr-1)
real(sp), parameter :: rw  = pir * w  ! solar angular velocity (radians hr-1)

! local variables

real(sp) :: wt1
real(sp) :: wpi

real(sp) :: e1
real(sp) :: e2

! ----

wpi  = 180. / (pi * w)
wt1  = rw * t1

if (a > b) then
  
  F = wpi * c(2) / sqrt(a**2 - b**2) * acos((b + a * cos(wt1)) / (a + b * cos(wt1))) + c(3) * t1

else if (a < b) then
  
  e1 = sqrt((b + a) * (1. + cos(wt1))) + sqrt((b - a) * (1. - cos(wt1)))
  e2 = sqrt((b + a) * (1. + cos(wt1))) - sqrt((b - a) * (1. - cos(wt1)))
  
  F = wpi * c(2) / sqrt(b**2 - a**2) * log(e1 / e2) + c(3) * t1

else

  F = wpi * c(2) / a * tan(wt1 / 2.) + c(3) * t1

end if

end function F

! ----------------------------------------------------------------------------------------------------------------

subroutine airmass(dlat,delta,dayl,Ratm,air)

! Calculation of the optical air mass from:
! Yin, X. (1997). Optical air mass: Daily integration and its applications. 
! Meteorology and Atmospheric Physics, 63(3-4), 227-233. doi:10.1007/Bf01027387

use parametersmod, only : sp,dp,pir => pir_sp
use typesmod,      only : airmasspars

implicit none

! arguments

real(dp),          intent(in)  :: dlat   ! latitude (degrees), double precision
real(sp),          intent(in)  :: delta  ! solar declination (degrees)
real(sp),          intent(in)  :: dayl   ! day length (hours)
real(sp),          intent(in)  :: Ratm   ! relative atmospheric pressure
type(airmasspars), intent(out) :: air    ! airmass parameters

! parameters

real(sp), parameter :: cos80 = cos(80. * pir)  ! (degrees)
real(sp), parameter :: mindayl = 2. * tiny(0._sp)

real(sp), parameter :: m0  =  1.   ! air mass at  0 degree solar zenith angle
real(sp), parameter :: m80 =  5.6  ! air mass at 80 degree solar zenith angle
real(sp), parameter :: m90 = 39.7  ! air mass at 90 degree solar zenith angle

real(sp), parameter :: w   = 15.     ! solar angular velocity (degrees hr-1)

! local variables

real(sp) :: mbar    ! daytime mean optical air mass (unitless, 1 at equatorial noon)
real(sp) :: mo      ! air mass at cosine zenith angle maximum
real(sp) :: mc      ! air mass at cosine zenith angle medium
real(sp) :: ml      ! air mass at cosine zenith angle bottom quarter range point

real(sp) :: lat     ! latitude (degrees)
real(sp) :: rlat    ! latitude (radians)
real(sp) :: rdelta  ! solar declination (radians)

real(sp) :: t1      ! number of hours between sunrise/sunset and solar noon (hr)
real(sp) :: t80     ! solar hour corresponding to the 80 degree zenith angle

! real(sp) :: t       ! solar hour (hr)

real(sp) :: Z       ! solar zenith angle (degrees)
real(sp) :: Zn      ! lesser of solar zenith angle at sunset or at midnight (degrees)
real(sp) :: Z0      ! zenith angle at solar noon (degrees)
real(sp) :: cosZ    ! cosine solar zenith angle (fraction), used in calculation of instantaneous air mass

! real(sp) :: l

! integer :: steps    ! integer number of time steps
! integer :: i        ! counter

real(sp) :: sinlat
real(sp) :: coslat
real(sp) :: sindel
real(sp) :: cosdel

real(sp)               :: a   ! values in equation 2.6b
real(sp)               :: b
real(sp), dimension(3) :: c

real(sp) :: tmp1
real(sp) :: tmp2
real(sp) :: tmp3

real(sp) :: tinv

real(sp) :: rZ0
real(sp) :: rZn

! -------------------------------------

lat = real(dlat,sp)

! calculate daily mean air mass (mbar)

if (dayl == 0.) then

  mbar = m90
  mc   = m90
  ml   = m90
  mo   = m90

else
  
  ! basic setup

  rlat   = pir * lat
  rdelta = pir * delta
  
  sinlat = sin(rlat)
  sindel = sin(rdelta)
  coslat = cos(rlat)
  cosdel = cos(rdelta)

  ! ------
  
  ! calculate the half-daylength
  
  if (dayl > mindayl) then
    t1 = 0.5 * dayl
  else
    t1 = dayl
  end if
  
  ! Eqn. 2.9
  if (abs(lat + delta) >= 90.) then
    Zn = acos(sinlat * sindel - coslat * cosdel) / pir
  else
    Zn = 90.
  end if

  ! Eqn. 2.10
  if (abs(lat - delta) >= 90.) then
    Z0 = 90.
  else
    Z0 = lat - delta
  end if
  
  rZ0 = Z0 * pir  ! convert to radians
  rZn = Zn * pir
  
  ! --------------------------
  ! calculate daytime mean airmass (mbar)

  b = coslat * cosdel
  
  if (t1 == 0.) then

    mbar = m90
    
  else if (abs(Zn) <= 80.) then
  
    tinv = 1. / t1

    c = c00
    a = c(1) + sinlat * sindel
    mbar = tinv * F(t1,a,b,c)

  else if (abs(Z0) >= 80.) then
  
    tinv = 1. / t1

    c = c80
    a = c(1) + sinlat * sindel
    mbar = tinv * F(t1,a,b,c)

  else
    
    t80 = 1. / w * acos((cos80 - sinlat * sindel) / (coslat * cosdel)) / pir  ! Eqn. 2.8

    c = c00
    a = c(1) + sinlat * sindel
        
    tmp1 = F(t80,a,b,c)

    c = c80
    a = c(1) + sinlat * sindel
    tmp2 = F(t1,a,b,c)

    c = c80
    a = c(1) + sinlat * sindel
    tmp3 = F(t80,a,b,c)

    tinv = 1. / t1

    mbar = tinv * (tmp1 + tmp2 - tmp3)
    
  end if

  ! --------------------------
  ! calculate instantaneous air mass at max, mid, and bottom quarter solar zenith angle (m0, mc, ml)
  
  Z = Z0

  cosZ = cos(Z * pir)

  if (Z <= 80.) then
    c = c00
  else
    c = c80
  end if
  
  mo = m(cosZ,c)

  ! --

  Z = (Z0 + Zn) / 2.
  
  cosz = (cos(rZ0) + cos(rZn)) / 2.

  if (Z <= 80.) then
    c = c00
  else
    c = c80
  end if

  mc = m(cosZ,c)

  ! --

  Z = (Z0 + 3. * Zn) / 4.
  
  cosz = (cos(rZ0) + 3. * cos(rZn)) / 4.

  if (Z <= 80.) then
    c = c00
  else
    c = c80
  end if

  ml = m(cosZ,c)

end if

! -------------------------------------
! correct calculated air mass for elevation

air%mbar = Ratm * mbar
air%mo   = Ratm * mo
air%mc   = Ratm * mc
air%ml   = Ratm * ml

end subroutine airmass

! ----------------------------------------------------------------------------------------------------------------

end module airmassmod