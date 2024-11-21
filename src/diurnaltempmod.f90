module diurnaltempmod

! Copied from arve-dgvm https://github.com/jedokaplan/arve-dgvm_original.git
! Edited by Leo O Lai, HKU, 2021

implicit none

public :: diurnaltemp

contains

! ---------------------------------------------------------

subroutine diurnaltemp(met)

! Based on: 
! Cesaraccio, C., Spano, D., Duce, P., & Snyder, R. L. (2001). 
! An improved model for determining degree-day values from daily temperature data. 
! Int J Biometeorol, 45(4), 161-169. doi:10.1007/s004840100104

! NB there is an error in this paper, eqn. 7a should have a sin after the alpha.

use parametersmod, only : i4,sp,pi
use typesmod,      only : metvars_out

implicit none

! arguments

type(metvars_out), dimension(:), intent(inout) :: met

! parameter

real(sp), parameter :: mindayl =  1.     ! minimum daylength for temperature calculations (h)
real(sp), parameter :: maxdayl = 23.     ! maximum daylength for temperature calculations (h)
real(sp), parameter :: tfpk    =  1./6.  ! delay between solar noon and peak temperature (fraction)

! local variables

real(sp) :: dayl0   ! length of the current day (h)
real(sp) :: dayl1   ! length of the next day (h)
real(sp) :: tmin0   ! current day minimum temperature (degC)
real(sp) :: tmin1   ! next day minimum temperature (degC)
real(sp) :: tmax    ! current day maximum temperature (degC)
real(sp) :: tday    ! mean daytime temperature (degC)
real(sp) :: tnight  ! mean nighttime temperature (degC)

real(sp) :: dl0   ! length of the current day, adjusted (h)
real(sp) :: dl1   ! length of the next day, adjusted (h)
real(sp) :: hdl0  ! half day length, current day (h)
real(sp) :: hdl1  ! half-day length, next day (h)

real(sp) :: sunrise0  ! time of sunrise, current day (h)
real(sp) :: sunrise1  ! time of sunrise, next day (h)
real(sp) :: sunset    ! time of sunset (h)
real(sp) :: hpeakt    ! time of peak temperature (h)

real(sp) :: t0    ! temperature at sunset (degC)
real(sp) :: r     ! difference between max temp and temp at sunset (degC)
real(sp) :: a     ! amplitude of temp change (max - min) of the current day (degC)
real(sp) :: b     ! b-factor, Cesaraccio et al. (2001), eqn. 10
real(sp) :: hni   ! length of the night (h)
real(sp) :: morn  ! time from sunrise to peakT (h)

! integrated temperature over intervals (degC)
real(sp) :: tam   ! sunrise to peakT
real(sp) :: tpm   ! peakT to sunset
real(sp) :: ti    ! sunset to midnight
real(sp) :: ti1   ! midnight to sunrise


! integer  :: sunrise
! integer  :: sunset
! integer  :: sunrise_n
! integer  :: dayhour
! integer  :: nighthour
! integer  :: dayhour_n

! -------------------------
! to separate daytime and nighttime: (Leo Lai May 2019)
! define daylength of polar night to be 1 hour
! define daylength of polar day to be 23 hours

dayl0 = met(1)%dayl
dayl1 = met(2)%dayl
tmin0 = met(1)%tmin
tmin1 = met(2)%tmin
tmax  = met(1)%tmax

dl0 = min(max(dayl0,mindayl),maxdayl)
dl1 = min(max(dayl1,mindayl),mindayl)

! these are prep for other calculations below
! tmin1 = minimum temperature of the following day

t0 = tmax - 0.39 * (tmax - tmin1)  ! temperature at sunset
a  = tmax - tmin0
r  = tmax - t0

! sunrise and sunset hour of the current day (relative to local solar noon)

hdl0 = 0.5 * dl0

sunrise0 = 12. - hdl0
sunset   = 12. + hdl0

! sunrise hour of the next day (relative to local solar noon)

hdl1 = 0.5 * dl1

sunrise1 = 12. - hdl1

! hour of peak temperature

hpeakt = 12. + 2. * hdl0 * tfpk

! length of the night

hni = 12. - sunset + sunrise1

! integral from sunset to sunrise of the next day (part III)

b = (tmin1 - t0) / sqrt(hni) ! Cesaraccio et al. eqn. 10

ti  = t0 * sunset  ! temperature at sunset

ti1 = t0 * (sunset + hni) + 2./3. * b * (hni**(3./2.))  ! temperature at sunrise (next day)

tnight = (ti1 - ti) / hni

! hours between sunrise and time of peak temperature

morn = sunrise0 - hpeakt  

! integral from sunrise to peak temperature (part I)

ti  = (tmin0 * sunrise0) + (1. / pi) * (2. * a * morn)
ti1 = (tmin0 * hpeakt)  + (1. / pi) * (2. * a * morn * cos(pi/2. * (sunrise0 - hpeakt) / morn))

tam = (ti1 - ti) / (-morn)

! integral from peak temperature to sunset (part II)

ti  = t0 * hpeakt - (1. / pi) * 8. * r * cos(pi / 8. * (-4.))
ti1 = t0 * sunset - (1. / pi) * 8. * r * cos(pi / 8. * (hpeakt - sunset - 4.))

tpm = (ti1 - ti) / (sunset - hpeakt)

! daytime mean temperature

tday = (tam + tpm) / 2.

! ---

met(1)%tday   = tday
met(1)%tnight = tnight

end subroutine diurnaltemp

! -----------------------------------------------------------

subroutine humidity(tmean,tdew,rhum)

! Calculate relative humidity from dew point temperature
! From Lawrence (2005) The relationship between relative humidity and the dewpoint temperature. A. MetSoc (https://doi.org/10.1175/BAMS-86-2-225)
! Equation (11)

use parametersmod, only : sp,dp

implicit none

real(sp), intent(in)  :: tmean
real(sp), intent(in)  :: tdew
real(sp), intent(out) :: rhum

real(sp) :: tmean_K
real(sp) :: tdew_K

real(sp), parameter :: L  = 2257000
real(sp), parameter :: Rw = 461.5

!------

tmean_K = tmean + 273.15
tdew_K  = tdew + 273.15

!------

rhum = 100. * exp(((1. - (tmean_K / tdew_K)) * (L / Rw)) / tmean_K)

if (rhum > 100.) rhum = 100.    ! When temp is lower than dew point


end subroutine humidity


end module diurnaltempmod