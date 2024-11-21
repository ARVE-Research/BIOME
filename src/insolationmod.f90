module insolationmod

! module with routines to calculate top of the atmosphere insolation
! based on subroutines and functions from the Laskar et al. (2004) insolation calculation programs (insolsub.f)

implicit none

public  :: truelon
public  :: insol
private :: hlm
private :: hl

contains

! ---------------------------------------------------------

subroutine insol(truelon,orbit,phi,solar)

! calculates total daily insolation at the top of the atmosphere for a given day and latitude

use parametersmod, only : sp,dp,pi,pir
use typesmod,      only : orbitpars,solarpars

implicit none

! arguments

real(dp),        intent(in)  :: truelon  ! true solar longitude relative to the equinox (radians) 
type(orbitpars), intent(in)  :: orbit    ! structure with orbital parameters
real(dp),        intent(in)  :: phi      ! latitude (radians)
type(solarpars), intent(out) :: solar    ! structure with toa radiation, daylength, and declination

! parameter

real(dp), parameter :: TSI = 1361._dp  ! total solar irradiance

! local variables

real(dp) :: v
real(dp) :: cosv
real(dp) :: aux
real(dp) :: rho
real(dp) :: sinusdelta
real(dp) :: cho
real(dp) :: ho
real(dp) :: a1
real(dp) :: a2
real(dp) :: eps    ! orbital obliquity (radians)
real(dp) :: pibar  ! geocentric longitude of perihelion (radians)
real(dp) :: w      ! mean daily top-of-the-atmosphere insolation (W m-2)
real(dp) :: dayl   ! day length (h)
real(dp) :: delta  ! solar declination (rad)

! ---

eps = pir * orbit%xob

pibar = pi + (pir * orbit%perh)

v = truelon - pibar  ! true anomaly

! Earth-Sun distance (rho)

cosv = cos(v)

aux = 1._dp + orbit%ecc * cosv

rho = (1._dp - orbit%ecc**2) / aux

! solar declination (delta)

sinusdelta = sin(eps) * sin(truelon)
delta = asin(sinusdelta)

! latitude of sunrise and sunset

aux = pi / 2._dp - abs(delta)

a1 = pi / 2._dp - delta
a2 = pi / 2._dp + delta

if (-aux < phi .and. phi < aux) then

  ! normal sunrise and sunset

  ! hour angle of sunrise and sunset (ho)

  cho = -tan(phi) * tan(delta)
  ho = acos(cho)
  
  dayl = 24._dp * ho / pi 

  ! Insolation (w)

  w = (ho * sin(phi) * sin(delta) + cos(phi) * cos(delta) * sin(ho))

  w = w * TSI / (pi * rho**2)
  
else if (phi > a1 .or. phi < -a2) then

  ! polar day

  w = TSI * sin(phi) * sin(delta) / rho**2
  
  dayl = 24._dp
  
else if (phi <= -a1 .or. phi >= a2) then

  ! polar night
  
  w = 0._dp
  
  dayl = 0._dp
  
else

  stop "error calculating aux or a1 a2"

end if

solar%rad0  = real(w)
solar%dayl  = real(dayl)
solar%delta = real(delta)

end subroutine insol

! ---------------------------------------------------------

real(dp) function truelon(orbit,cal,doy)

! calculate the true longitude of a day of the year expressed as mean longitude relative to the vernal equinox

use parametersmod, only : dp,pi,pir
use typesmod,      only : orbitpars,calendartype

implicit none

! arguments

type(orbitpars),    intent(in) :: orbit  ! structure with orbital parameters
type(calendartype), intent(in) :: cal    ! structure with number of days in the year and other info
integer,            intent(in) :: doy    ! day of the year

! local variables

real(dp) :: pibar  ! geocentric longitude of perihelion (radians) 
real(dp) :: hlm0   ! mean longitude of the vernal equinox (radians)
real(dp) :: hlm1   ! mean longitude of the input day of the year (radians)

! ---
! convert heliocentric longitude of perihelion (deg) to geocentric longitude (radians)

pibar = pi + (pir * orbit%perh)

! mean longitude of the vernal equinox (21 March)

hlm0 = hlm(0._dp,orbit%ecc,pibar)

! mean longitude of the selected date relative to the vernal equinox (radians)

hlm1 = hlm0 + 2._dp * pi * (real(doy) - cal%veqday) / real(cal%ndyr)

! true longitude of the selected date

truelon = hl(hlm1,orbit%ecc,pibar)

end function truelon

! ---------------------------------------------------------

real(dp) function hlm(hl,e,pibar)

! calculates mean longitude (hlm) as a function of true longitude (hl), eccentricity (e), 
! and the longitude of perihelion (pibar)

use parametersmod, only : dp

implicit none

real(dp), intent(in) :: hl     ! true longitude (radians)
real(dp), intent(in) :: e      ! eccentricity
real(dp), intent(in) :: pibar  ! longitude of perihelion (radians)

real(dp) :: eV

! ---

eV = hl - pibar

hlm = hl - 2._dp * e * sin(eV) +                                 &
      (3._dp * e**2 / 4._dp + e**4 / 8._dp) * sin(2._dp * eV) -  &
      (e**3 / 3._dp + e**5 / 8._dp) * sin(3._dp * eV ) +         &
       5._dp * e**4 / 32._dp * sin(4._dp * eV) -                 &
       3._dp * e**5 / 40._dp * sin(5._dp * eV)

end function hlm

! ---------------------------------------------------------

real(dp) function hl(hlm,e,pibar)

! calculates true longitude (hl) as a function of mean longitude (hlm), eccentricity (e), 
! and the longitude of perihelion (pibar)

use parametersmod, only : dp

implicit none

real(dp), intent(in) :: hlm    ! mean longitude (radians)
real(dp), intent(in) :: e      ! eccentricity
real(dp), intent(in) :: pibar  ! longitude of perihelion (radians)

real(dp) :: eM

! ---

eM = hlm - pibar

hl = hlm + (2._dp * e - e**3 / 4._dp + 5._dp * e**5 / 96._dp) * sin(eM) +          &
           (5._dp * e**2 /   4._dp - 11._dp * e**4 / 24._dp) * sin(2._dp * eM) +   &
          (13._dp * e**3 /  12._dp - 43._dp * e**5._dp / 64) * sin(3._dp * eM) +   &
          103._dp * e**4 /  96._dp * sin(4._dp * eM) +                             &
         1097._dp * e**5 / 960._dp * sin(5._dp * eM)

end function hl

! ---------------------------------------------------------

end module insolationmod