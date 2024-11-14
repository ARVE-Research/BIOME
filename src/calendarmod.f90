module calendarmod

use parametersmod, only : dp

implicit none

public  :: initcalendar
private :: kepler_time
private :: getmonlen
private :: monlen
private :: adjust_to_ref_length

contains

! ---------------------------------------------------------

subroutine kepler_time(eccen,T,theta_deg,time)

! travel time along orbit relative to periapsis/perihelion (i.e. 0.0 at perihelion)
! input:  true anomaly (theta, degrees); output:  traverse time since perihelion (time, same units as T)
! (Curtis, H.D. (2014) Orbital Mechanics for Engineering Students, Elsevier, Ch. 3)
! adapted from P.J. Bartlein & S.L. Shafer (2019) Paleo calendar-effect adjustments in time-slice and transient
!         climate-model simulations (PaleoCalAdjust v1.0): impact and strategies for data analysis,
!         Geoscientific Model Development, 12:3889-3913, doi: 10.5194/gmd-12-3889-2019
!   - available from GitHub:  https://github.com/pjbartlein/PaleoCalAdjust or Zenodo:  https://doi.org/10.5281/zenodo.1478824

use parametersmod, only : dp,pi,pir

implicit none

! arguments

real(dp), intent(in)  :: eccen      ! orbital eccentricity (unitless)
real(dp), intent(in)  :: T          ! orbital period (yrlen, days)
real(dp), intent(in)  :: theta_deg  ! "true" anomaly (angle from perihelion (degrees)
real(dp), intent(out) :: time       ! traverse time from periapsis (e.g. perihelion, units?)

! local variables

real(dp) :: theta_rad  ! theta_rad (radians)
real(dp) :: E          ! eccentric anomaly
real(dp) :: M          ! mean anomaly

! ----
    
theta_rad = theta_deg * pir
    
E = 2._dp * atan(sqrt((1._dp - eccen) / (1._dp + eccen)) * tan(theta_rad/ 2._dp) )  ! eq 3.13b

M = E - eccen * sin(E)                                                              ! eq 3.14

time = (M / (2._dp * pi)) * T                                                       ! eq 3.15

if (time < 0._dp) time = time + T
    
end subroutine kepler_time    

! ---------------------------------------------------------

subroutine monlen(yrlen,veqday,imonlen,orbit,rmonlen,rmonbeg,rmonmid,rmonend) !,&
!                  tt_VE_to_SS,tt_SS_to_AE,tt_AE_to_WS,tt_WS_to_VE)

! calculate month lengths, and beginning, middle and ending day times for a particular orbital configuration, 
! calendar and vernal equinox day using a "traverse-time/time-of-flight" expression based on Kepler's equation:
! (Curtis, H.D. (2014) Orbital Mechanics for Engineering Students, Elsevier, Ch. 3)
! adapted from P.J. Bartlein & S.L. Shafer (2019) Paleo calendar-effect adjustments in time-slice and transient
!         climate-model simulations (PaleoCalAdjust v1.0): impact and strategies for data analysis,
!         Geoscientific Model Development, 12:3889-3913, doi: 10.5194/gmd-12-3889-2019
!   - available from GitHub:  https://github.com/pjbartlein/PaleoCalAdjust or Zenodo:  https://doi.org/10.5281/zenodo.1478824

use parametersmod, only : dp,pi,nmos
use typesmod,      only : orbitpars

implicit none

! arguments

real(dp),               intent(in)  :: yrlen    ! total length of year (days)
real(dp),               intent(in)  :: veqday   ! day of year of the vernal equinox 
integer, dimension(:),  intent(in)  :: imonlen  ! number of days in each month at present (calendar dependent)
type(orbitpars),        intent(in)  :: orbit    ! orbital parameters
real(dp), dimension(:), intent(out) :: rmonlen  ! month lengths
real(dp), dimension(:), intent(out) :: rmonbeg  ! month beginning day
real(dp), dimension(:), intent(out) :: rmonmid  ! month middle day
real(dp), dimension(:), intent(out) :: rmonend  ! month ending day

! local variables

real(dp) :: eccen
real(dp) :: perih
real(dp) :: tt_VE_to_SS          ! traverse time, VE to SS
real(dp) :: tt_SS_to_AE          ! traverse time, SS to AE    
real(dp) :: tt_AE_to_WS          ! traverse time, AE to WS
real(dp) :: tt_WS_to_VE          ! traverse time, WS to VE
real(dp) :: angle_VE_to_Jan1     ! angle between vernal equinox and Jan 1
real(dp) :: tt_perih_to_VE       ! traverse time, perihelion to vernal equinox
real(dp) :: angle_perih_to_Jan1  ! angle, perihelion to Jan 1
real(dp) :: angle_perih_to_VE    ! angle, perihelion to vernal (March) equinox
real(dp) :: angle_perih_to_SS    ! angle, perihelion to summer (June) solstice
real(dp) :: angle_perih_to_AE    ! angle, perihelion to autumn (September) equinox
real(dp) :: angle_perih_to_WS    ! angle, perihelion to winter (December) solstice
real(dp) :: tt_perih_to_SS       ! traverse time, perihelion to summer solstice
real(dp) :: tt_perih_to_AE       ! traverse time, perihelion to autumn equinox
real(dp) :: tt_perih_to_WS       ! traverse time, perihelion to winter solstice

real(dp), dimension(nmos) :: mon_angle  ! month angle (degrees)
real(dp), dimension(nmos) :: beg_angle
real(dp), dimension(nmos) :: mid_angle 
real(dp), dimension(nmos) :: end_angle  ! beginning, middle and end of month, relative to Jan 1

! angles and traverse times for beginning, middle and ending days of each month 

real(dp), dimension(nmos) :: perih_angle_month_beg
real(dp), dimension(nmos) :: tt_month_beg
real(dp), dimension(nmos) :: t_month_beg
real(dp), dimension(nmos) :: perih_angle_month_mid
real(dp), dimension(nmos) :: tt_month_mid
real(dp), dimension(nmos) :: t_month_mid
real(dp), dimension(nmos) :: perih_angle_month_end
real(dp), dimension(nmos) :: tt_month_end
real(dp), dimension(nmos) :: t_month_end

integer  :: m  ! indices

! --------

eccen = orbit%ecc
perih = orbit%perh

rmonlen = 0._dp
rmonbeg = 0._dp
rmonmid = 0._dp
rmonend = 0._dp

! ----    
! angle and traverse time, perihelion to vernal equinox 
angle_perih_to_VE = 360._dp - perih

! traverse time (along elliptical orbit) perihelion to vernal equinox
call kepler_time(eccen,yrlen,angle_perih_to_VE,tt_perih_to_VE)

! angle, perihelion to Jan1
angle_perih_to_Jan1 = angle_perih_to_VE - veqday * (360._dp / yrlen)

! angle, vernal equinox to Jan1 (Jan 1 begins at 0.0 (or 360.0 degrees))
angle_VE_to_Jan1 = 360._dp - veqday * (360._dp / yrlen)
    
! angle, perihelion to SS
angle_perih_to_SS = angle_perih_to_VE + 90._dp
    
! traverse time (along elliptical orbit) perihelion to SS
call kepler_time(eccen,yrlen,angle_perih_to_SS,tt_perih_to_SS)
    
! traverse time (along elliptical orbit) VE to SS
tt_VE_to_SS = tt_perih_to_SS - tt_perih_to_VE

if (abs(tt_VE_to_SS) > 180._dp) tt_VE_to_SS = yrlen - abs(tt_VE_to_SS)
    
! angle, perihelion to AE
angle_perih_to_AE = angle_perih_to_VE + 180._dp
    
! traverse time (along elliptical orbit) perihelion to AE
call kepler_time(eccen,yrlen,angle_perih_to_AE,tt_perih_to_AE)
    
! traverse time (along elliptical orbit) SS to AE
tt_SS_to_AE = tt_perih_to_AE - tt_perih_to_SS

if (abs(tt_SS_to_AE) > 180._dp) tt_SS_to_AE = yrlen - abs(tt_SS_to_AE)
    
! angle, perihelion to WS
angle_perih_to_WS = angle_perih_to_VE + 270._dp
    
! traverse time (along elliptical orbit) perihelion to WS
call kepler_time(eccen,yrlen,angle_perih_to_WS,tt_perih_to_WS)
    
! traverse time (along elliptical orbit) AE to WS
tt_AE_to_WS = tt_perih_to_WS - tt_perih_to_AE

if (abs(tt_AE_to_WS) > 180._dp) tt_AE_to_WS = yrlen - abs(tt_AE_to_WS)

! traverse time (along elliptical orbit) WS to VE
tt_WS_to_VE = tt_perih_to_VE - tt_perih_to_WS

if (abs(tt_WS_to_VE) > 180._dp) tt_WS_to_VE = yrlen - abs(tt_WS_to_VE)  

! angles, traverse times, month lengths, etc. for individual months

! month angle (and beginning, middle and end of month, relative to Jan1)

mon_angle(1) = real(imonlen(1),dp) * 360._dp / yrlen

beg_angle(1) = 0._dp

mid_angle(1) = beg_angle(1) + mon_angle(1) / 2._dp

end_angle(1) = mon_angle(1)

do m = 2,nmos

  mon_angle(m) = real(imonlen(m),dp) * 360._dp / yrlen

  beg_angle(m) = beg_angle(m-1) + mon_angle(m-1)

  mid_angle(m) = beg_angle(m) + mon_angle(m) / 2._dp

  end_angle(m) = beg_angle(m) + mon_angle(m)

end do

! ----
! angles from perihelion etc. for month beginning, mid, and ending days

do m = 1,nmos
    
  ! angle from perihelion for each month's beginning, middle and ending days (on circular orbit)

  perih_angle_month_beg(m) = angle_perih_to_Jan1 + beg_angle(m) 
  perih_angle_month_mid(m) = angle_perih_to_Jan1 + mid_angle(m) 
  perih_angle_month_end(m) = angle_perih_to_Jan1 + end_angle(m)  
  
  ! traverse times from perihelion for each month's beginning, middle and ending days (on elliptical orbit)

  call kepler_time(eccen,yrlen,perih_angle_month_beg(m),tt_month_beg(m))
  call kepler_time(eccen,yrlen,perih_angle_month_mid(m),tt_month_mid(m))
  call kepler_time(eccen,yrlen,perih_angle_month_end(m),tt_month_end(m))
    
  ! traverse time from vernal equinox for each month's beginning, middle and ending days (on elliptical orbit)

  t_month_beg(m) = tt_month_beg(m) - tt_perih_to_VE
  t_month_mid(m) = tt_month_mid(m) - tt_perih_to_VE
  t_month_end(m) = tt_month_end(m) - tt_perih_to_VE
    
  ! beginning, middle and ending days of each month (relative to Jan1)

  rmonbeg(m) = dmod(t_month_beg(m) + veqday, yrlen)

  if (perih_angle_month_beg(m) > 360._dp) rmonbeg(m) = rmonbeg(m) + yrlen

  rmonmid(m) = dmod(t_month_mid(m) + veqday, yrlen)

  if (perih_angle_month_mid(m) > 360._dp) rmonmid(m) = rmonmid(m) + yrlen

  rmonend(m) = dmod(t_month_end(m) + veqday, yrlen) 

  if (perih_angle_month_end(m) > 360._dp) rmonend(m) = rmonend(m) + yrlen
    
end do
    
! fixup for ending day of last month (when end of last month appers in beginning of year)

if (rmonend(nmos) < 30._dp) rmonend(nmos) = rmonend(nmos) + yrlen  

! month length (month beginning to next month beginning)

do m = 1,nmos-1

  rmonlen(m) = t_month_beg(m+1) - t_month_beg(m)

  if (rmonlen(m) <= 0._dp) rmonlen(m) = rmonlen(m) + yrlen

end do

rmonlen(nmos) = t_month_beg(1) - t_month_beg(nmos)

if (rmonlen(nmos) <= 0._dp) rmonlen(nmos) = rmonlen(nmos) + yrlen

end subroutine monlen

! ---------------------------------------------------------

subroutine adjust_to_ref_length(rmonlen,rmonlenref,rmonlentarg)

use parametersmod, only : dp

implicit none

! arguments

real(dp), dimension(:), intent(inout) :: rmonlen      ! real month lengths
real(dp), dimension(:), intent(in)    :: rmonlenref   ! reference month lengths (usually calculated 0 ka)
real(dp), dimension(:), intent(in)    :: rmonlentarg  ! target month lengths (usually conventional 0 ka)

! ----
! adjust rmonlen to reference year

rmonlen = rmonlen * (rmonlentarg / rmonlenref)
    
end subroutine adjust_to_ref_length

! ---------------------------------------------------------

subroutine getmonlen(orbit,ndyr,ndm)

! get number of days per month

use parametersmod, only : dp,nmos,nd_365,nd_366,veqday_365,present_mon_noleap,veqday_366,present_mon_leap
use typesmod,      only : monthinfotype,orbitpars

implicit none

! arguments

type(orbitpars),        intent(in)  :: orbit
integer,                intent(in)  :: ndyr
real(dp), dimension(:), intent(out) :: ndm

! local variables

logical :: leapyear

real(dp) :: yrlen
real(dp) :: veqday

integer,  dimension(nmos) :: imonlen
real(dp), dimension(nmos) :: rmonbeg 
real(dp), dimension(nmos) :: rmonmid
real(dp), dimension(nmos) :: rmonend

! ----
! calculate real month lengths

if (ndyr == 365) then
  leapyear = .false.
else if (ndyr == 366) then
  leapyear = .true.
else
  write(0,*)'error getmonlen: ndyr not supported'
  stop
end if

if (leapyear) then

  yrlen = real(nd_366,dp)
  
  veqday = veqday_366
  
  imonlen = nint(present_mon_leap)
  
  call monlen(yrlen,veqday,imonlen,orbit,ndm,rmonbeg,rmonmid,rmonend)

else

  yrlen = real(nd_365,dp)

  veqday = veqday_365
  
  imonlen = nint(present_mon_noleap)
  
  call monlen(yrlen,veqday,imonlen,orbit,ndm,rmonbeg,rmonmid,rmonend)

end if

! ----

end subroutine getmonlen

! ---------------------------------------------------------

subroutine initcalendar(yrbp,orbit,noleap,leapyr)

! given a year BP, calculate the orbital parameters and the month lengths of a standard and leap year

use parametersmod
use orbitmod
use typesmod

implicit none

! arguments

integer,             intent(in)  :: yrbp
type(orbitpars),     intent(out) :: orbit
type(calendartype),  intent(out) :: noleap
type(calendartype),  intent(out) :: leapyr

! local variables

integer :: m

real(dp), dimension(nmos) :: ndm0noleap
real(dp), dimension(nmos) :: ndm0leapyr
real(dp), dimension(nmos) :: ndm1noleap
real(dp), dimension(nmos) :: ndm1leapyr

real(dp) :: yrlen0
real(dp) :: yrlen1

! ----
! initialize the reference month lengths for ~0ka for 365- and 366-day years

call getorbitpars(0,orbit)  ! estimate orbital parameters for 0ka

! get number of days per month (real number) for 0ka 365 and 366-day years

call getmonlen(orbit,365,ndm0noleap)
call getmonlen(orbit,366,ndm0leapyr)

! initialize the orbital parameters for the selected run year

call getorbitpars(yrbp,orbit)  ! estimate orbital parameters for 0ka

! calculate real month lengths for the target year

call getmonlen(orbit,365,ndm1noleap)
call getmonlen(orbit,366,ndm1leapyr)

! calculate real month lengths for the target year by normalizing to 0ka

noleap%ndmr = ndm1noleap / ndm0noleap * present_mon_noleap

leapyr%ndmr = ndm1leapyr / ndm0leapyr * present_mon_leap

! --
! correct for total year length: standard year

yrlen0 = 365._dp
yrlen1 = sum(noleap%ndmr)

noleap%ndmr = noleap%ndmr / (yrlen1 / yrlen0)

noleap%ndmi = nint(noleap%ndmr)

! --
! correct for total year length: leap year

yrlen0 = 366._dp
yrlen1 = sum(leapyr%ndmr)

leapyr%ndmr = leapyr%ndmr / (yrlen1 / yrlen0)

leapyr%ndmi = nint(leapyr%ndmr)

! --

do m = 1,nmos
  write(0,*)m,ndm0noleap(m),ndm1noleap(m),noleap%ndmr(m),noleap%ndmi(m)
end do

write(0,*)

do m = 1,nmos
  write(0,*)m,ndm0leapyr(m),ndm1leapyr(m),leapyr%ndmr(m),leapyr%ndmi(m)
end do

write(0,*)sum(noleap%ndmi),sum(leapyr%ndmi)

end subroutine initcalendar

! ---------------------------------------------------------

end module calendarmod
